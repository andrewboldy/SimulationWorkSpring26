//----------------------------------------------------------------------------------
//
// twoElectronEventCutterHistogrammer.C
// Written by Andrew Boldy University of South Carolina, 2026
// Assisted by Codex
//
// Purpose:
//   1. Loop over an EventNtuple ROOT file or filelist.
//   2. Find events with exactly two reconstructed e-minus tracks whose
//      trkmcsim MC match is a valid rank-0 electron.
//   3. Write those selected event numbers and angular information to a text file.
//   4. Make ROOT histograms and PDF plots for:
//        - trkmcsim momentum and origin positions for the selected two-electron events.
//        - rank-0-only trkmcsim origin-position plots for those selected events.
//        - trkmcsim theta/phi angular distributions for those events.
//        - reconstructed trksegs momentum at TT_Front, TT_Mid, and TT_Back.
//        - a log-y reconstructed trksegs momentum companion plot for tail checks.
//        - comparisons of all events vs selected two-electron events.
//
// Vocabulary used in this macro:
//   "trkmcsim" is the MC-truth match information attached to a reconstructed track.
//   "trksegs" are reconstructed track-segment states at detector surfaces.
//   "rank 0" is EventNtuple's best MC-truth match to a reconstructed track.
//   "selected two-electron event" means exactly two valid rank-0 trkmcsim electrons
//   were found among the reconstructed e-minus tracks in that event.
//
//----------------------------------------------------------------------------------

// Standard C++ includes.
// These give us strings, vectors, input/output files, math helpers, and min/max tools.
#include <algorithm>
#include <cmath>
#include <string>
#include <fstream>
#include <iostream>
#include <limits>
#include <vector>

// CERN ROOT includes.
// TH1F and TH2F make histograms, TCanvas draws/saves them, TFile writes ROOT output,
// and TLegend/TSystem/TStopwatch handle legends, directories, and timing.
#include <TCanvas.h>
#include <TDirectory.h>
#include <TFile.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TLegend.h>
#include <TPad.h>
#include <TROOT.h>
#include <TSystem.h>
#include <TStopwatch.h>

// Mu2e/EventNtuple helper includes.
// RooUtil opens the ntuple and exposes convenient Event/Track access.
// common_cuts.hh gives cuts like is_e_minus, tracker_entrance, tracker_middle, etc.
#include "EventNtuple/rooutil/inc/RooUtil.hh"
#include "EventNtuple/rooutil/inc/common_cuts.hh"

// Namespace shortcuts.
// This keeps the code readable below by avoiding std:: and rooutil:: everywhere.
using namespace std;
using namespace rooutil;

void twoElectronEventCutterHistogrammer(const string& generatorName, const string& fileName)
{
  //----------------------------------------------------------------------------
  // Basic input check
  //
  // ROOT/RooUtil will do the real ntuple reading later.  This small check catches
  // an obviously wrong filename or filelist path early and prints a clear error.
  //----------------------------------------------------------------------------
  ifstream file(fileName);
  if (!file.is_open())
  {
    cerr << "ERROR: could not open input file or filelist: " << fileName << endl;
    return;
  }
  file.close();

  // Open the ntuple through RooUtil and ask it how many events it can see.
  // This number controls the main event loop below.
  RooUtil util(fileName);
  int numEvents = util.GetNEvents();
  cout << "There are " << numEvents << " in the filelist." << endl;

  // Text output file.  This is a human-readable list of the selected events.
  // The generatorName is included so separate samples do not overwrite each other.
  const string outputFileName = "twoElectronEvents_" + generatorName + ".txt";
  ofstream outputFile(outputFileName);
  if (!outputFile.is_open())
  {
    cerr << "ERROR: could not create output text file: " << outputFileName << endl;
    return;
  }
  outputFile << "entry run subrun event num_rank0_electrons "
             << "electron1_theta_deg electron1_phi_deg "
             << "electron2_theta_deg electron2_phi_deg" << endl;

  // Running ROOT in batch mode avoids opening interactive canvases.
  // That matters on the VM, where lingering GUI canvases can make ROOT slow to quit.
  //
  // TH1::AddDirectory(false) keeps new histograms from automatically attaching
  // themselves to whatever ROOT directory happens to be current.  That makes
  // ownership and cleanup more predictable.
  //
  // The old settings are saved so they can be restored at the end of the macro.
  const bool wasBatchMode = gROOT->IsBatch();
  const bool oldAddDirectoryStatus = TH1::AddDirectoryStatus();
  gROOT->SetBatch(true);
  TH1::AddDirectory(false);
  cout << "Running ROOT plotting in batch mode to avoid slow VM canvas cleanup." << endl;

  //----------------------------------------------------------------------------
  // Monte Carlo truth
  //
  // These histograms are filled only for events that pass the two-electron
  // selection below.  "All electrons" means every valid electron SimInfo object
  // in the selected reconstructed e-minus tracks' trkmcsim vectors.  "Rank 0"
  // keeps the subset that EventNtuple classifies as the best MC match to that
  // reconstructed track.
  //----------------------------------------------------------------------------
  // Plot ranges used by several histograms.
  //
  // Momentum is intentionally capped at 65 MeV/c because the features of interest
  // are below that scale for this study.  The Z origin range is restricted to the
  // region of interest so XZ/YZ maps are less visually compressed.
  const double originXYMin = -1000.0;
  const double originXYMax = 1000.0;
  const double originZMin = -6000.0;
  const double originZMax = 2100.0;
  const int originXYBins = 200;
  const int originZBins = 405;
  const double momentumMin = 0.0;
  const double momentumMax = 65.0;
  const double thetaDegMin = 0.0;
  const double thetaDegMax = 180.0;
  const double phiDegMin = -180.0;
  const double phiDegMax = 180.0;
  const int thetaDegBins = 36;
  const int phiDegBins = 36;
  const double cosThetaMin = -1.0;
  const double cosThetaMax = 1.0;
  const double recoTrkSegMomentumLinearXMin = 20.0;
  const string surfaceSuffixes[3] = {"Front", "Middle", "Back"};
  const string surfaceLabels[3] = {"TT_Front", "TT_Mid", "TT_Back"};

  // Main selected-event trkmcsim momentum histograms.
  // These are filled only after an event passes the exactly-two-rank-0-electron cut.
  TH1F* hMCTAllElectronMomentum = new TH1F(
    "hMCTAllElectronMomentum",
    "Monte Carlo truth: all trkmcsim electrons;trkmcsim momentum [MeV/c];Electrons",
    130, momentumMin, momentumMax);
  TH1F* hMCTRank0ElectronMomentum = new TH1F(
    "hMCTRank0ElectronMomentum",
    "Monte Carlo truth: rank-0 trkmcsim electrons;trkmcsim momentum [MeV/c];Electrons",
    130, momentumMin, momentumMax);

  // Event-selection comparison trkmcsim momentum histograms.
  //
  // Each comparison needs two histograms with the same binning:
  //   - "AllEvents" gets filled before the two-electron selection cut.
  //   - "TwoElectronEvents" gets filled only after the event passes the cut.
  //
  // There are two comparison pairs:
  //   - all valid trkmcsim electrons
  //   - valid trkmcsim electrons with rank == 0
  TH1F* hCompareMCTAllElectronMomentumAllEvents = new TH1F(
    "hCompareMCTAllElectronMomentumAllEvents",
    "Monte Carlo truth all electrons: all events vs selected two-electron events;trkmcsim momentum [MeV/c];Electrons",
    130, momentumMin, momentumMax);
  TH1F* hCompareMCTAllElectronMomentumTwoElectronEvents = new TH1F(
    "hCompareMCTAllElectronMomentumTwoElectronEvents",
    "Monte Carlo truth all electrons: selected two-electron events;trkmcsim momentum [MeV/c];Electrons",
    130, momentumMin, momentumMax);
  TH1F* hCompareMCTRank0ElectronMomentumAllEvents = new TH1F(
    "hCompareMCTRank0ElectronMomentumAllEvents",
    "Monte Carlo truth rank-0 electrons: all events vs selected two-electron events;trkmcsim momentum [MeV/c];Electrons",
    130, momentumMin, momentumMax);
  TH1F* hCompareMCTRank0ElectronMomentumTwoElectronEvents = new TH1F(
    "hCompareMCTRank0ElectronMomentumTwoElectronEvents",
    "Monte Carlo truth rank-0 electrons: selected two-electron events;trkmcsim momentum [MeV/c];Electrons",
    130, momentumMin, momentumMax);

  // 1D origin-position histograms for all valid trkmcsim electrons in selected events.
  // These show the separate x, y, and z distributions.
  TH1F* hMCTAllElectronOriginX = new TH1F(
    "hMCTAllElectronOriginX",
    "Monte Carlo truth: all trkmcsim electron origins X;origin x [mm];Electrons",
    originXYBins, originXYMin, originXYMax);
  TH1F* hMCTAllElectronOriginY = new TH1F(
    "hMCTAllElectronOriginY",
    "Monte Carlo truth: all trkmcsim electron origins Y;origin y [mm];Electrons",
    originXYBins, originXYMin, originXYMax);
  TH1F* hMCTAllElectronOriginZ = new TH1F(
    "hMCTAllElectronOriginZ",
    "Monte Carlo truth: all trkmcsim electron origins Z;origin z [mm];Electrons",
    originZBins, originZMin, originZMax);

  // 1D origin-position histograms for only rank-0 trkmcsim electrons in selected events.
  // These are the same coordinates as above, but restricted to the best MC match.
  TH1F* hMCTRank0ElectronOriginX = new TH1F(
    "hMCTRank0ElectronOriginX",
    "Monte Carlo truth: rank-0 trkmcsim electron origins X;origin x [mm];Electrons",
    originXYBins, originXYMin, originXYMax);
  TH1F* hMCTRank0ElectronOriginY = new TH1F(
    "hMCTRank0ElectronOriginY",
    "Monte Carlo truth: rank-0 trkmcsim electron origins Y;origin y [mm];Electrons",
    originXYBins, originXYMin, originXYMax);
  TH1F* hMCTRank0ElectronOriginZ = new TH1F(
    "hMCTRank0ElectronOriginZ",
    "Monte Carlo truth: rank-0 trkmcsim electron origins Z;origin z [mm];Electrons",
    originZBins, originZMin, originZMax);

  // 2D origin maps for all valid trkmcsim electrons in selected events.
  //
  // The XZ and YZ histograms put z on the horizontal axis.  That is deliberate:
  // it makes the long detector direction easier to read and avoids the old issue
  // where the z coordinate looked artificially stuck at +/-1000 mm.
  TH2F* hMCTAllElectronOriginXY = new TH2F(
    "hMCTAllElectronOriginXY",
    "Monte Carlo truth: all trkmcsim electron origins XY;origin x [mm];origin y [mm]",
    originXYBins, originXYMin, originXYMax, originXYBins, originXYMin, originXYMax);
  TH2F* hMCTAllElectronOriginXZ = new TH2F(
    "hMCTAllElectronOriginXZ",
    "Monte Carlo truth: all trkmcsim electron origins XZ;origin z [mm];origin x [mm]",
    originZBins, originZMin, originZMax, originXYBins, originXYMin, originXYMax);
  TH2F* hMCTAllElectronOriginYZ = new TH2F(
    "hMCTAllElectronOriginYZ",
    "Monte Carlo truth: all trkmcsim electron origins YZ;origin z [mm];origin y [mm]",
    originZBins, originZMin, originZMax, originXYBins, originXYMin, originXYMax);

  // 2D origin maps for only the rank-0 trkmcsim electrons in selected events.
  TH2F* hMCTRank0ElectronOriginXY = new TH2F(
    "hMCTRank0ElectronOriginXY",
    "Monte Carlo truth: rank-0 trkmcsim electron origins XY;origin x [mm];origin y [mm]",
    originXYBins, originXYMin, originXYMax, originXYBins, originXYMin, originXYMax);
  TH2F* hMCTRank0ElectronOriginXZ = new TH2F(
    "hMCTRank0ElectronOriginXZ",
    "Monte Carlo truth: rank-0 trkmcsim electron origins XZ;origin z [mm];origin x [mm]",
    originZBins, originZMin, originZMax, originXYBins, originXYMin, originXYMax);
  TH2F* hMCTRank0ElectronOriginYZ = new TH2F(
    "hMCTRank0ElectronOriginYZ",
    "Monte Carlo truth: rank-0 trkmcsim electron origins YZ;origin z [mm];origin y [mm]",
    originZBins, originZMin, originZMax, originXYBins, originXYMin, originXYMax);

  // Angular histograms for valid trkmcsim electrons in selected two-electron events.
  //
  // "All" keeps every valid electron MC match in a selected event.  "Rank 0"
  // keeps only the two best MC matches that define the selected event.  The 1D
  // theta/phi histograms use coarser bins to reduce low-statistics bin chatter.
  TH1F* hMCTAllElectronTheta = new TH1F(
    "hMCTAllElectronTheta",
    "Monte Carlo truth: all trkmcsim electron #theta;#theta [deg];Electrons",
    thetaDegBins, thetaDegMin, thetaDegMax);
  TH1F* hMCTAllElectronPhi = new TH1F(
    "hMCTAllElectronPhi",
    "Monte Carlo truth: all trkmcsim electron #phi;#phi [deg];Electrons",
    phiDegBins, phiDegMin, phiDegMax);
  TH2F* hMCTAllElectronCosThetaVsPhi = new TH2F(
    "hMCTAllElectronCosThetaVsPhi",
    "Monte Carlo truth: all trkmcsim electron cos(#theta) vs #phi;#phi [deg];cos(#theta);Electrons",
    180, phiDegMin, phiDegMax, 100, cosThetaMin, cosThetaMax);

  TH1F* hMCTRank0ElectronTheta = new TH1F(
    "hMCTRank0ElectronTheta",
    "Monte Carlo truth: rank-0 trkmcsim electron #theta;#theta [deg];Electrons",
    thetaDegBins, thetaDegMin, thetaDegMax);
  TH1F* hMCTRank0ElectronPhi = new TH1F(
    "hMCTRank0ElectronPhi",
    "Monte Carlo truth: rank-0 trkmcsim electron #phi;#phi [deg];Electrons",
    phiDegBins, phiDegMin, phiDegMax);
  TH2F* hMCTRank0ElectronCosThetaVsPhi = new TH2F(
    "hMCTRank0ElectronCosThetaVsPhi",
    "Monte Carlo truth: rank-0 trkmcsim electron cos(#theta) vs #phi;#phi [deg];cos(#theta);Electrons",
    180, phiDegMin, phiDegMax, 100, cosThetaMin, cosThetaMax);
  TH2F* hMCTRank0ElectronTheta1VsTheta2 = new TH2F(
    "hMCTRank0ElectronTheta1VsTheta2",
    "Monte Carlo truth: rank-0 electron #theta pairs;electron 1 #theta [deg];electron 2 #theta [deg];Events",
    180, thetaDegMin, thetaDegMax, 180, thetaDegMin, thetaDegMax);

  //----------------------------------------------------------------------------
  // Reconstructed track-segment momentum
  //
  // These histograms use the reconstructed momentum stored in the trksegs
  // branch at the tracker front, middle, and back surfaces for the same
  // selected two-electron events.  The comparison histograms split trksegs
  // into tracks with any valid electron MC match and tracks with a rank-0
  // electron MC match.
  //----------------------------------------------------------------------------
  // Main selected-event trksegs momentum histograms.
  //
  // These three histograms are the original "selected two-electron events only"
  // reconstructed-momentum plots.  They use all reconstructed e-minus trksegs in
  // the selected events, separated by tracker surface.
  TH1F* hRecoTrkSegMomentumFront = new TH1F(
    "hRecoTrkSegMomentumFront",
    "Reconstructed trksegs momentum at tracker surfaces;trksegs momentum [MeV/c];Segments",
    130, momentumMin, momentumMax);
  TH1F* hRecoTrkSegMomentumMiddle = new TH1F(
    "hRecoTrkSegMomentumMiddle",
    "Reconstructed trksegs momentum at TT_Mid;trksegs momentum [MeV/c];Segments",
    130, momentumMin, momentumMax);
  TH1F* hRecoTrkSegMomentumBack = new TH1F(
    "hRecoTrkSegMomentumBack",
    "Reconstructed trksegs momentum at TT_Back;trksegs momentum [MeV/c];Segments",
    130, momentumMin, momentumMax);

  // Event-selection comparison trksegs histograms.
  //
  // Arrays of length 3 are used here because the same comparison is made for
  // three surfaces:
  //   index 0 -> TT_Front
  //   index 1 -> TT_Mid
  //   index 2 -> TT_Back
  //
  // For each surface there are four histograms:
  //   1. electron-matched trksegs from all events
  //   2. electron-matched trksegs from selected two-electron events
  //   3. rank-0-electron-matched trksegs from all events
  //   4. rank-0-electron-matched trksegs from selected two-electron events
  TH1F* hCompareRecoTrkSegAllElectronMomentumAllEvents[3] = {nullptr, nullptr, nullptr};
  TH1F* hCompareRecoTrkSegAllElectronMomentumTwoElectronEvents[3] = {nullptr, nullptr, nullptr};
  TH1F* hCompareRecoTrkSegRank0ElectronMomentumAllEvents[3] = {nullptr, nullptr, nullptr};
  TH1F* hCompareRecoTrkSegRank0ElectronMomentumTwoElectronEvents[3] = {nullptr, nullptr, nullptr};

  // Build those twelve comparison histograms in a loop instead of copy-pasting
  // three nearly identical blocks.  The name string is what ROOT stores internally;
  // the title string is what appears on the plot.
  for (int i_surface = 0; i_surface < 3; ++i_surface)
  {
    hCompareRecoTrkSegAllElectronMomentumAllEvents[i_surface] = new TH1F(
      ("hCompareRecoTrkSegAllElectronMomentum" + surfaceSuffixes[i_surface] + "AllEvents").c_str(),
      ("All electron-matched trksegs at " + surfaceLabels[i_surface] + ": all events vs selected two-electron events;trksegs momentum [MeV/c];Segments").c_str(),
      130, momentumMin, momentumMax);
    hCompareRecoTrkSegAllElectronMomentumTwoElectronEvents[i_surface] = new TH1F(
      ("hCompareRecoTrkSegAllElectronMomentum" + surfaceSuffixes[i_surface] + "TwoElectronEvents").c_str(),
      ("All electron-matched trksegs at " + surfaceLabels[i_surface] + ": selected two-electron events;trksegs momentum [MeV/c];Segments").c_str(),
      130, momentumMin, momentumMax);
    hCompareRecoTrkSegRank0ElectronMomentumAllEvents[i_surface] = new TH1F(
      ("hCompareRecoTrkSegRank0ElectronMomentum" + surfaceSuffixes[i_surface] + "AllEvents").c_str(),
      ("Rank-0 electron-matched trksegs at " + surfaceLabels[i_surface] + ": all events vs selected two-electron events;trksegs momentum [MeV/c];Segments").c_str(),
      130, momentumMin, momentumMax);
    hCompareRecoTrkSegRank0ElectronMomentumTwoElectronEvents[i_surface] = new TH1F(
      ("hCompareRecoTrkSegRank0ElectronMomentum" + surfaceSuffixes[i_surface] + "TwoElectronEvents").c_str(),
      ("Rank-0 electron-matched trksegs at " + surfaceLabels[i_surface] + ": selected two-electron events;trksegs momentum [MeV/c];Segments").c_str(),
      130, momentumMin, momentumMax);
  }

  //----------------------------------------------------------------------------
  // Bookkeeping counters and temporary summary variables
  //
  // These are not histograms.  They are printed to the terminal so you can sanity
  // check what was filled, especially when a plot looks surprising.
  //----------------------------------------------------------------------------

  // Entries that pass the exactly-two-rank-0-electron selection.
  vector<int> dualElectronEntries;

  // trkmcsim electron counters.
  // "allEvent..." means the count is accumulated before applying the event cut.
  // "selected..." means the count is accumulated only after the event passes.
  long long allEventAllElectronFills = 0;
  long long allEventRank0ElectronFills = 0;
  long long selectedAllElectronFills = 0;
  long long selectedRank0ElectronFills = 0;

  // trksegs counters for electron-matched reconstructed tracks.
  long long allEventRecoTrkSegFrontFills = 0;
  long long allEventRecoTrkSegMiddleFills = 0;
  long long allEventRecoTrkSegBackFills = 0;

  // trksegs counters for all reconstructed e-minus tracks in selected events.
  // These correspond to the older standalone reconstructed trksegs momentum plot.
  long long selectedRecoTrkSegFrontFills = 0;
  long long selectedRecoTrkSegMiddleFills = 0;
  long long selectedRecoTrkSegBackFills = 0;

  // trksegs counters for tracks with any valid trkmcsim electron match.
  long long selectedAllElectronRecoTrkSegFrontFills = 0;
  long long selectedAllElectronRecoTrkSegMiddleFills = 0;
  long long selectedAllElectronRecoTrkSegBackFills = 0;

  // trksegs counters for tracks with a valid rank-0 trkmcsim electron match.
  long long allEventRank0RecoTrkSegFrontFills = 0;
  long long allEventRank0RecoTrkSegMiddleFills = 0;
  long long allEventRank0RecoTrkSegBackFills = 0;
  long long selectedRank0RecoTrkSegFrontFills = 0;
  long long selectedRank0RecoTrkSegMiddleFills = 0;
  long long selectedRank0RecoTrkSegBackFills = 0;

  // Track the angular range of the two selected rank-0 electrons.
  const double pi = acos(-1.0);
  double minElectronAngleDeg = numeric_limits<double>::max();
  double maxElectronAngleDeg = -numeric_limits<double>::max();
  bool haveElectronAngleStats = false;

  // Small helper function:
  //   input:  a momentum vector-like object with .R() and .z()
  //   output: cos(theta), where theta is the angle from the z-axis
  //
  // The clamp to [-1, 1] protects later acos calls from tiny floating-point
  // roundoff errors.  The -2 sentinel means the momentum magnitude was invalid.
  auto getCosTheta = [](const auto& momentum) {
    const double momentumMagnitude = momentum.R();
    if (momentumMagnitude <= 0.0)
    {
      return -2.0;
    }

    double cosTheta = momentum.z() / momentumMagnitude;
    cosTheta = max(-1.0, min(1.0, cosTheta));
    return cosTheta;
  };

  // Convert the polar angle to degrees for the text output and raw theta plots.
  auto getAngleFromZDeg = [pi, &getCosTheta](const auto& momentum) {
    const double cosTheta = getCosTheta(momentum);
    if (cosTheta < -1.0)
    {
      return -1.0;
    }
    return acos(cosTheta) * 180.0 / pi;
  };

  // Azimuthal angle around the z-axis in degrees, using the conventional
  // atan2(y, x) range of [-180, 180].
  auto getPhiDeg = [pi](const auto& momentum) {
    return atan2(momentum.y(), momentum.x()) * 180.0 / pi;
  };

  //----------------------------------------------------------------------------
  // Main event loop
  //
  // This is where almost all of the analysis choices happen.
  //
  // Important ordering:
  //   1. Read one event.
  //   2. Loop over reconstructed e-minus tracks in that event.
  //   3. Fill "all events" comparison histograms immediately.
  //   4. Count rank-0 electrons.
  //   5. If the event has exactly two rank-0 electrons, also fill the selected
  //      two-electron histograms and write the event to the text file.
  //
  // That ordering is what lets the comparison plots show:
  //   all events vs selected two-electron events
  //----------------------------------------------------------------------------
  for (int i_event = 0; i_event < numEvents; i_event++)
  {
    // Number of valid rank-0 trkmcsim electrons found in this event.
    // This is reset for every event.
    int numEventElectrons = 0;

    // Temporary containers for this one event only.
    // They gather information while looping over the event's tracks.
    // If the event later passes the exactly-two-electron cut, these vectors are
    // used to fill the selected-event histograms.
    vector<double> eventRank0ElectronThetaDeg;
    vector<double> eventRank0ElectronPhiDeg;
    vector<double> eventRank0ElectronCosTheta;
    vector<const mu2e::SimInfo*> eventElectrons;
    vector<const mu2e::SimInfo*> eventRank0Electrons;

    // Reconstructed trksegs momenta for all reconstructed e-minus tracks.
    // These feed the original selected-event trksegs overlay plot.
    vector<double> eventRecoTrkSegFrontMomenta;
    vector<double> eventRecoTrkSegMiddleMomenta;
    vector<double> eventRecoTrkSegBackMomenta;

    // Reconstructed trksegs momenta for tracks that have any valid electron
    // trkmcsim match.  These feed the "all electrons" trksegs comparison plots.
    vector<double> eventAllElectronRecoTrkSegFrontMomenta;
    vector<double> eventAllElectronRecoTrkSegMiddleMomenta;
    vector<double> eventAllElectronRecoTrkSegBackMomenta;

    // Reconstructed trksegs momenta for tracks whose MC match is specifically
    // a rank-0 electron.  These feed the rank-0 trksegs comparison plots.
    vector<double> eventRank0RecoTrkSegFrontMomenta;
    vector<double> eventRank0RecoTrkSegMiddleMomenta;
    vector<double> eventRank0RecoTrkSegBackMomenta;

    // Pull the current event from RooUtil.  Then keep only reconstructed tracks
    // fit with the e-minus hypothesis, because this study is about electrons.
    auto& event = util.GetEvent(i_event);
    auto e_minus_tracks = event.GetTracks(is_e_minus);

    // Loop over reconstructed e-minus tracks in this event.
    // The track is not const because Track::GetSegments(...) is not a const method.
    for (auto& track : e_minus_tracks)
    {
      // These two flags summarize what was found in this track's trkmcsim vector.
      //
      // trackHasAnyElectron:
      //   true if this track has at least one valid electron MC match.
      //
      // trackHasRank0Electron:
      //   true if this track has at least one valid electron MC match with rank == 0.
      //
      // The trksegs comparison histograms use these flags after the trkmcsim loop.
      bool trackHasAnyElectron = false;
      bool trackHasRank0Electron = false;

      // trkmcsim can be null if this reconstructed track has no attached MC-truth
      // match vector.  Always check the pointer before dereferencing it.
      if (track.trkmcsim != nullptr)
      {
        // Loop over every MC-truth match attached to this reconstructed track.
        for (const auto& mctrack : *(track.trkmcsim))
        {
          // Keep only real, valid electrons.
          //
          // valid must be true so we ignore empty/default SimInfo objects.
          // pdg == 11 means electron.  This does not require rank == 0 yet.
          if (!(mctrack.valid && mctrack.pdg == 11))
          {
            continue;
          }

          // At this point the SimInfo object is a valid electron.
          // Store a pointer to it for possible selected-event filling later.
          eventElectrons.push_back(&mctrack);
          trackHasAnyElectron = true;

          // Fill the all-events trkmcsim comparison immediately.
          // This happens before the exactly-two-electron event cut.
          ++allEventAllElectronFills;
          hCompareMCTAllElectronMomentumAllEvents->Fill(mctrack.mom.R());

          // rank == 0 selects EventNtuple's best MC match to this reconstructed track.
          if (mctrack.rank == 0)
          {
            trackHasRank0Electron = true;

            // This counter is the event-selection variable.
            // After the track loop, the event passes only if this equals 2.
            numEventElectrons++;

            // Save the rank-0 electron for selected-event histograms.
            eventRank0Electrons.push_back(&mctrack);

            // Compute and store this electron's angles for the text file and
            // selected-event angular histograms.
            eventRank0ElectronThetaDeg.push_back(getAngleFromZDeg(mctrack.mom));
            eventRank0ElectronPhiDeg.push_back(getPhiDeg(mctrack.mom));
            eventRank0ElectronCosTheta.push_back(getCosTheta(mctrack.mom));

            // Fill the rank-0 all-events comparison before the event cut.
            ++allEventRank0ElectronFills;
            hCompareMCTRank0ElectronMomentumAllEvents->Fill(mctrack.mom.R());
          } //end work on the electron pdg tracks
        }//end the interface with the TTree at trkmcsim branch level
      }

      // Pull reconstructed track segments at the three tracker surfaces.
      // has_reco_step ensures we are using reconstructed trkseg momentum, not
      // an MC-only SurfaceStep.
      auto frontSegments = track.GetSegments([](TrackSegment& segment) {
        return has_reco_step(segment) && tracker_entrance(segment);
      });
      auto middleSegments = track.GetSegments([](TrackSegment& segment) {
        return has_reco_step(segment) && tracker_middle(segment);
      });
      auto backSegments = track.GetSegments([](TrackSegment& segment) {
        return has_reco_step(segment) && tracker_exit(segment);
      });

      // TT_Front segment momenta for this track.
      for (auto& segment : frontSegments)
      {
        // trkseg->mom.R() is the reconstructed momentum magnitude at this surface.
        const double momentum = segment.trkseg->mom.R();

        // Always keep the momentum for the original selected-event trksegs plot.
        eventRecoTrkSegFrontMomenta.push_back(momentum);

        // Fill the all-electron-matched comparison only if this track had any
        // valid trkmcsim electron match.
        if (trackHasAnyElectron)
        {
          eventAllElectronRecoTrkSegFrontMomenta.push_back(momentum);
          hCompareRecoTrkSegAllElectronMomentumAllEvents[0]->Fill(momentum);
          ++allEventRecoTrkSegFrontFills;
        }
        if (trackHasRank0Electron)
        {
          eventRank0RecoTrkSegFrontMomenta.push_back(momentum);
          hCompareRecoTrkSegRank0ElectronMomentumAllEvents[0]->Fill(momentum);
          ++allEventRank0RecoTrkSegFrontFills;
        }
      }

      // TT_Mid segment momenta for this track.
      for (auto& segment : middleSegments)
      {
        const double momentum = segment.trkseg->mom.R();
        eventRecoTrkSegMiddleMomenta.push_back(momentum);
        if (trackHasAnyElectron)
        {
          eventAllElectronRecoTrkSegMiddleMomenta.push_back(momentum);
          hCompareRecoTrkSegAllElectronMomentumAllEvents[1]->Fill(momentum);
          ++allEventRecoTrkSegMiddleFills;
        }
        if (trackHasRank0Electron)
        {
          eventRank0RecoTrkSegMiddleMomenta.push_back(momentum);
          hCompareRecoTrkSegRank0ElectronMomentumAllEvents[1]->Fill(momentum);
          ++allEventRank0RecoTrkSegMiddleFills;
        }
      }

      // TT_Back segment momenta for this track.
      for (auto& segment : backSegments)
      {
        const double momentum = segment.trkseg->mom.R();
        eventRecoTrkSegBackMomenta.push_back(momentum);
        if (trackHasAnyElectron)
        {
          eventAllElectronRecoTrkSegBackMomenta.push_back(momentum);
          hCompareRecoTrkSegAllElectronMomentumAllEvents[2]->Fill(momentum);
          ++allEventRecoTrkSegBackFills;
        }
        if (trackHasRank0Electron)
        {
          eventRank0RecoTrkSegBackMomenta.push_back(momentum);
          hCompareRecoTrkSegRank0ElectronMomentumAllEvents[2]->Fill(momentum);
          ++allEventRank0RecoTrkSegBackFills;
        }
      }
    }//end the interface with the reconstructed e-minus track level

    // The main event selection:
    // keep only events with exactly two valid rank-0 trkmcsim electrons.
    //
    // If this is not true, the "all events" comparison histograms have already
    // received their fills, but no selected-event histograms are filled.
    if (numEventElectrons != 2)
    {
      numEventElectrons=0;
      continue;
    }//end the check to continue if the number of rank 0 electrons is not exactly 2

    else
    {
      // These default values make the text output robust even if evtinfo is absent.
      int run = -1;
      int subrun = -1;
      int eventNumber = -1;

      // Print the entry number immediately so terminal output can be matched to
      // the selected event list.
      cout << "Entry " << i_event;
      if (event.evtinfo != nullptr)
      {
        // evtinfo gives the real run/subrun/event identifiers.
        run = event.evtinfo->run;
        subrun = event.evtinfo->subrun;
        eventNumber = event.evtinfo->event;
        cout << " Run " << event.evtinfo->run
             << " Subrun " << event.evtinfo->subrun
             << " Event " << event.evtinfo->event;
      }
      cout << " has exactly " << numEventElectrons << " electrons in MonteCarlo truth able to be reconstructed." << endl;
      cout << "  Electron 1 theta from z-axis: " << eventRank0ElectronThetaDeg.at(0)
           << " degrees, phi: " << eventRank0ElectronPhiDeg.at(0) << " degrees" << endl;
      cout << "  Electron 2 theta from z-axis: " << eventRank0ElectronThetaDeg.at(1)
           << " degrees, phi: " << eventRank0ElectronPhiDeg.at(1) << " degrees" << endl;

      // Write one selected event per line to the text file.
      // The two angle pairs are safe to access with .at(0) and .at(1) because
      // this block only runs when numEventElectrons == 2.
      outputFile << i_event << " "
                 << run << " "
                 << subrun << " "
                 << eventNumber << " "
                 << numEventElectrons << " "
                 << eventRank0ElectronThetaDeg.at(0) << " "
                 << eventRank0ElectronPhiDeg.at(0) << " "
                 << eventRank0ElectronThetaDeg.at(1) << " "
                 << eventRank0ElectronPhiDeg.at(1) << endl;

      // Update the minimum and maximum angle summary over all selected events.
      for (const double electronAngleDeg : eventRank0ElectronThetaDeg)
      {
        if (electronAngleDeg >= 0.0)
        {
          haveElectronAngleStats = true;
          minElectronAngleDeg = min(minElectronAngleDeg, electronAngleDeg);
          maxElectronAngleDeg = max(maxElectronAngleDeg, electronAngleDeg);
        }
      }

      // Save this entry number and update simple printed counters.
      dualElectronEntries.push_back(i_event);
      selectedAllElectronFills += eventElectrons.size();
      selectedRank0ElectronFills += eventRank0Electrons.size();
      selectedRecoTrkSegFrontFills += eventRecoTrkSegFrontMomenta.size();
      selectedRecoTrkSegMiddleFills += eventRecoTrkSegMiddleMomenta.size();
      selectedRecoTrkSegBackFills += eventRecoTrkSegBackMomenta.size();
      selectedAllElectronRecoTrkSegFrontFills += eventAllElectronRecoTrkSegFrontMomenta.size();
      selectedAllElectronRecoTrkSegMiddleFills += eventAllElectronRecoTrkSegMiddleMomenta.size();
      selectedAllElectronRecoTrkSegBackFills += eventAllElectronRecoTrkSegBackMomenta.size();
      selectedRank0RecoTrkSegFrontFills += eventRank0RecoTrkSegFrontMomenta.size();
      selectedRank0RecoTrkSegMiddleFills += eventRank0RecoTrkSegMiddleMomenta.size();
      selectedRank0RecoTrkSegBackFills += eventRank0RecoTrkSegBackMomenta.size();

      // One point per selected event, using the two rank-0 electrons in track-loop
      // order.  This shows the theta-theta correlation between the two tracks.
      if (eventRank0ElectronThetaDeg.at(0) >= 0.0 && eventRank0ElectronThetaDeg.at(1) >= 0.0)
      {
        hMCTRank0ElectronTheta1VsTheta2->Fill(eventRank0ElectronThetaDeg.at(0),
                                              eventRank0ElectronThetaDeg.at(1));
      }

      // Fill selected-event trkmcsim histograms for all valid electrons.
      // This includes rank-0 and non-rank-0 valid electron matches in selected events.
      for (const auto* electron : eventElectrons)
      {
        hMCTAllElectronMomentum->Fill(electron->mom.R());

        // Also fill the selected side of the event-selection comparison.
        hCompareMCTAllElectronMomentumTwoElectronEvents->Fill(electron->mom.R());

        // Fill 1D origin coordinates.
        hMCTAllElectronOriginX->Fill(electron->pos.x());
        hMCTAllElectronOriginY->Fill(electron->pos.y());
        hMCTAllElectronOriginZ->Fill(electron->pos.z());

        // Fill 2D origin maps.
        // XY uses x,y.  XZ and YZ use z on the horizontal axis.
        hMCTAllElectronOriginXY->Fill(electron->pos.x(), electron->pos.y());
        hMCTAllElectronOriginXZ->Fill(electron->pos.z(), electron->pos.x());
        hMCTAllElectronOriginYZ->Fill(electron->pos.z(), electron->pos.y());

        // Fill angular plots for every valid electron match in the selected event.
        const double thetaDeg = getAngleFromZDeg(electron->mom);
        const double cosTheta = getCosTheta(electron->mom);
        if (thetaDeg >= 0.0 && cosTheta >= -1.0)
        {
          const double phiDeg = getPhiDeg(electron->mom);
          hMCTAllElectronTheta->Fill(thetaDeg);
          hMCTAllElectronPhi->Fill(phiDeg);
          hMCTAllElectronCosThetaVsPhi->Fill(phiDeg, cosTheta);
        }
      }

      // Fill selected-event trkmcsim histograms for rank-0 electrons only.
      for (unsigned int i_electron = 0; i_electron < eventRank0Electrons.size(); ++i_electron)
      {
        const auto* electron = eventRank0Electrons.at(i_electron);
        hMCTRank0ElectronMomentum->Fill(electron->mom.R());

        // Selected side of the rank-0 trkmcsim comparison.
        hCompareMCTRank0ElectronMomentumTwoElectronEvents->Fill(electron->mom.R());

        // 1D origin coordinates for rank-0 electrons.
        hMCTRank0ElectronOriginX->Fill(electron->pos.x());
        hMCTRank0ElectronOriginY->Fill(electron->pos.y());
        hMCTRank0ElectronOriginZ->Fill(electron->pos.z());

        // 2D origin maps for rank-0 electrons.
        hMCTRank0ElectronOriginXY->Fill(electron->pos.x(), electron->pos.y());
        hMCTRank0ElectronOriginXZ->Fill(electron->pos.z(), electron->pos.x());
        hMCTRank0ElectronOriginYZ->Fill(electron->pos.z(), electron->pos.y());

        // Fill angular plots for the two rank-0 electrons that define the event.
        const double thetaDeg = eventRank0ElectronThetaDeg.at(i_electron);
        const double phiDeg = eventRank0ElectronPhiDeg.at(i_electron);
        const double cosTheta = eventRank0ElectronCosTheta.at(i_electron);
        if (thetaDeg >= 0.0 && cosTheta >= -1.0)
        {
          hMCTRank0ElectronTheta->Fill(thetaDeg);
          hMCTRank0ElectronPhi->Fill(phiDeg);
          hMCTRank0ElectronCosThetaVsPhi->Fill(phiDeg, cosTheta);
        }
      }

      // Fill the original selected-event reconstructed trksegs momentum histograms.
      // These use all reconstructed e-minus trksegs in the selected event.
      for (const double momentum : eventRecoTrkSegFrontMomenta)
      {
        hRecoTrkSegMomentumFront->Fill(momentum);
      }
      for (const double momentum : eventRecoTrkSegMiddleMomenta)
      {
        hRecoTrkSegMomentumMiddle->Fill(momentum);
      }
      for (const double momentum : eventRecoTrkSegBackMomenta)
      {
        hRecoTrkSegMomentumBack->Fill(momentum);
      }

      // Fill the selected side of the all-electron-matched trksegs comparisons.
      for (const double momentum : eventAllElectronRecoTrkSegFrontMomenta)
      {
        hCompareRecoTrkSegAllElectronMomentumTwoElectronEvents[0]->Fill(momentum);
      }
      for (const double momentum : eventAllElectronRecoTrkSegMiddleMomenta)
      {
        hCompareRecoTrkSegAllElectronMomentumTwoElectronEvents[1]->Fill(momentum);
      }
      for (const double momentum : eventAllElectronRecoTrkSegBackMomenta)
      {
        hCompareRecoTrkSegAllElectronMomentumTwoElectronEvents[2]->Fill(momentum);
      }

      // Fill the selected side of the rank-0-electron-matched trksegs comparisons.
      for (const double momentum : eventRank0RecoTrkSegFrontMomenta)
      {
        hCompareRecoTrkSegRank0ElectronMomentumTwoElectronEvents[0]->Fill(momentum);
      }
      for (const double momentum : eventRank0RecoTrkSegMiddleMomenta)
      {
        hCompareRecoTrkSegRank0ElectronMomentumTwoElectronEvents[1]->Fill(momentum);
      }
      for (const double momentum : eventRank0RecoTrkSegBackMomenta)
      {
        hCompareRecoTrkSegRank0ElectronMomentumTwoElectronEvents[2]->Fill(momentum);
      }

      numEventElectrons=0;
    }
  } //end the interface with the TTree at event level

  //----------------------------------------------------------------------------
  // Print summary information after the event loop
  //
  // All histogram filling is complete here.  These terminal messages help check
  // whether the selected-event counts and all-event comparison counts look sane.
  //----------------------------------------------------------------------------
  outputFile.close();

  // Number of events that passed the exactly-two-rank-0-electron selection.
  const string terminalBoldOn = "\033[1m";
  const string terminalTextReset = "\033[0m";
  cout << "Wrote " << terminalBoldOn << dualElectronEntries.size() << terminalTextReset
       << " selected events with two electron reconstructed tracks to " << outputFileName << endl;

  // Print the angular range of the selected rank-0 electrons.
  if (haveElectronAngleStats)
  {
    cout << "Minimum electron angle from z-axis: " << minElectronAngleDeg << " degrees" << endl;
    cout << "Maximum electron angle from z-axis: " << maxElectronAngleDeg << " degrees" << endl;
  }
  else
  {
    cout << "No selected events, so no minimum or maximum electron angle from z-axis was computed." << endl;
  }

  // Print fill counts for trkmcsim momentum comparisons.
  cout << "Monte Carlo truth electron fills from all events: all valid electrons = "
       << allEventAllElectronFills << ", rank-0 electrons = " << allEventRank0ElectronFills << endl;
  cout << "Monte Carlo truth electron fills from selected events: all valid electrons = "
       << selectedAllElectronFills << ", rank-0 electrons = " << selectedRank0ElectronFills << endl;

  // Print fill counts for reconstructed trksegs comparisons.
  // These distinguish electron-matched tracks, all reconstructed e-minus tracks,
  // and rank-0-electron-matched tracks.
  cout << "All-electron-matched trksegs momentum fills from all events: front = "
       << allEventRecoTrkSegFrontFills
       << ", middle = " << allEventRecoTrkSegMiddleFills
       << ", back = " << allEventRecoTrkSegBackFills << endl;
  cout << "All-electron-matched trksegs momentum fills from selected events: front = "
       << selectedAllElectronRecoTrkSegFrontFills
       << ", middle = " << selectedAllElectronRecoTrkSegMiddleFills
       << ", back = " << selectedAllElectronRecoTrkSegBackFills << endl;
  cout << "All reconstructed e-minus trksegs momentum fills from selected events: front = "
       << selectedRecoTrkSegFrontFills
       << ", middle = " << selectedRecoTrkSegMiddleFills
       << ", back = " << selectedRecoTrkSegBackFills << endl;
  cout << "Rank-0 electron-matched trksegs momentum fills from all events: front = "
       << allEventRank0RecoTrkSegFrontFills
       << ", middle = " << allEventRank0RecoTrkSegMiddleFills
       << ", back = " << allEventRank0RecoTrkSegBackFills << endl;
  cout << "Rank-0 electron-matched trksegs momentum fills from selected events: front = "
       << selectedRank0RecoTrkSegFrontFills
       << ", middle = " << selectedRank0RecoTrkSegMiddleFills
       << ", back = " << selectedRank0RecoTrkSegBackFills << endl;

  //----------------------------------------------------------------------------
  // Prepare ROOT output file and directories
  //
  // The ROOT file keeps the histograms and canvases.  The PDFs are easier to
  // view quickly, but the ROOT file is better if you want to inspect bins later.
  //----------------------------------------------------------------------------
  const string histogramFileName = "twoElectronEvents_" + generatorName + "_histograms.root";
  TFile histogramFile(histogramFileName.c_str(), "RECREATE");

  // ROOT directories keep related objects grouped together in the output file.
  auto* mcTruthDirectory = histogramFile.mkdir("Monte Carlo truth");
  auto* mcTruthAngleDirectory = histogramFile.mkdir("Monte Carlo truth angles");
  auto* recoTrkSegMomentumDirectory = histogramFile.mkdir("Reconstructed track-segment momentum");
  auto* eventSelectionComparisonDirectory = histogramFile.mkdir("Event selection comparisons");

  // Start by writing into the Monte Carlo truth directory.
  if (mcTruthDirectory != nullptr)
  {
    mcTruthDirectory->cd();
  }

  // Re-apply the intended z ranges before saving/drawing.
  // This protects against ROOT autoscaling z-axis histograms to a misleading range.
  hMCTAllElectronOriginZ->GetXaxis()->SetRangeUser(originZMin, originZMax);
  hMCTRank0ElectronOriginZ->GetXaxis()->SetRangeUser(originZMin, originZMax);
  hMCTAllElectronOriginXZ->GetXaxis()->SetRangeUser(originZMin, originZMax);
  hMCTAllElectronOriginYZ->GetXaxis()->SetRangeUser(originZMin, originZMax);
  hMCTRank0ElectronOriginXZ->GetXaxis()->SetRangeUser(originZMin, originZMax);
  hMCTRank0ElectronOriginYZ->GetXaxis()->SetRangeUser(originZMin, originZMax);
  cout << "Monte Carlo truth Z plot range: " << originZMin << " to " << originZMax << " mm" << endl;
  cout << "Momentum histogram entries: all = " << hMCTAllElectronMomentum->GetEntries()
       << ", rank-0 = " << hMCTRank0ElectronMomentum->GetEntries() << endl;
  cout << "Momentum histogram visible integrals: all = "
       << hMCTAllElectronMomentum->Integral(1, hMCTAllElectronMomentum->GetNbinsX())
       << ", rank-0 = "
       << hMCTRank0ElectronMomentum->Integral(1, hMCTRank0ElectronMomentum->GetNbinsX()) << endl;

  // Sanity check:
  // Rank-0 electrons are a subset of all valid electrons.  Therefore the rank-0
  // momentum histogram should not be larger than the all-electron histogram in
  // any bin.  If it is, something is inconsistent in the filling logic.
  int rank0GreaterThanAllBins = 0;
  for (int bin = 0; bin <= hMCTAllElectronMomentum->GetNbinsX() + 1; ++bin)
  {
    const double allBinContent = hMCTAllElectronMomentum->GetBinContent(bin);
    const double rank0BinContent = hMCTRank0ElectronMomentum->GetBinContent(bin);
    if (rank0BinContent > allBinContent)
    {
      ++rank0GreaterThanAllBins;
      if (rank0GreaterThanAllBins <= 5)
      {
        cout << "WARNING: rank-0 momentum bin " << bin
             << " has " << rank0BinContent
             << " entries, greater than all-electron bin content "
             << allBinContent << endl;
      }
    }
  }
  if (rank0GreaterThanAllBins == 0)
  {
    cout << "Momentum sanity check passed: rank-0 is never greater than all electrons bin-by-bin." << endl;
  }

  // Write selected-event trkmcsim momentum, origin, and 2D origin-map histograms.
  // These objects are stored under the "Monte Carlo truth" directory.
  hMCTAllElectronMomentum->Write();
  hMCTRank0ElectronMomentum->Write();
  hMCTAllElectronOriginX->Write();
  hMCTAllElectronOriginY->Write();
  hMCTAllElectronOriginZ->Write();
  hMCTRank0ElectronOriginX->Write();
  hMCTRank0ElectronOriginY->Write();
  hMCTRank0ElectronOriginZ->Write();
  hMCTAllElectronOriginXY->Write();
  hMCTAllElectronOriginXZ->Write();
  hMCTAllElectronOriginYZ->Write();
  hMCTRank0ElectronOriginXY->Write();
  hMCTRank0ElectronOriginXZ->Write();
  hMCTRank0ElectronOriginYZ->Write();

  if (mcTruthAngleDirectory != nullptr)
  {
    mcTruthAngleDirectory->cd();
  }

  // Write selected-event angular histograms.
  hMCTAllElectronTheta->Write();
  hMCTAllElectronPhi->Write();
  hMCTAllElectronCosThetaVsPhi->Write();
  hMCTRank0ElectronTheta->Write();
  hMCTRank0ElectronPhi->Write();
  hMCTRank0ElectronCosThetaVsPhi->Write();
  hMCTRank0ElectronTheta1VsTheta2->Write();

  if (recoTrkSegMomentumDirectory != nullptr)
  {
    recoTrkSegMomentumDirectory->cd();
  }

  // Write the original selected-event reconstructed trksegs momentum histograms.
  hRecoTrkSegMomentumFront->Write();
  hRecoTrkSegMomentumMiddle->Write();
  hRecoTrkSegMomentumBack->Write();

  if (eventSelectionComparisonDirectory != nullptr)
  {
    eventSelectionComparisonDirectory->cd();
  }

  // Write all event-selection comparison histograms:
  //   - four trkmcsim momentum histograms
  //   - twelve trksegs momentum histograms, one comparison set for each surface
  hCompareMCTAllElectronMomentumAllEvents->Write();
  hCompareMCTAllElectronMomentumTwoElectronEvents->Write();
  hCompareMCTRank0ElectronMomentumAllEvents->Write();
  hCompareMCTRank0ElectronMomentumTwoElectronEvents->Write();
  for (int i_surface = 0; i_surface < 3; ++i_surface)
  {
    hCompareRecoTrkSegAllElectronMomentumAllEvents[i_surface]->Write();
    hCompareRecoTrkSegAllElectronMomentumTwoElectronEvents[i_surface]->Write();
    hCompareRecoTrkSegRank0ElectronMomentumAllEvents[i_surface]->Write();
    hCompareRecoTrkSegRank0ElectronMomentumTwoElectronEvents[i_surface]->Write();
  }

  //----------------------------------------------------------------------------
  // Draw and save PDF plots
  //
  // All PDFs go into Plots/.  gSystem->mkdir(..., true) creates that directory
  // if it does not already exist and does nothing harmful if it already exists.
  //----------------------------------------------------------------------------
  const string plotsDirectory = "Plots";
  gSystem->mkdir(plotsDirectory.c_str(), true);

  // Timer for only the drawing/writing-to-PDF part.  This helped diagnose ROOT
  // lag on the VM separately from the event loop itself.
  TStopwatch plotTimer;
  plotTimer.Start();

  // Helper used by the new comparison canvases.
  //
  // It draws two histograms on the same pad:
  //   black/gray = all events
  //   red        = selected two-electron events
  //
  // The y-axis maximum is set from the larger of the two so neither histogram is
  // clipped.  The legend is centered and made smaller to avoid blocking the plot.
  auto drawEventSelectionComparison = [](TH1F* allEventsHist,
                                         TH1F* twoElectronEventsHist,
                                         const string& allEventsLabel,
                                         const string& twoElectronEventsLabel) {
    allEventsHist->SetLineColor(kBlack);
    allEventsHist->SetLineWidth(2);
    allEventsHist->SetFillColor(kGray);
    allEventsHist->SetFillStyle(3004);
    twoElectronEventsHist->SetLineColor(kRed);
    twoElectronEventsHist->SetLineWidth(2);
    const double maxBinContent = max(allEventsHist->GetMaximum(), twoElectronEventsHist->GetMaximum());
    if (maxBinContent > 0.0)
    {
      allEventsHist->SetMaximum(1.15 * maxBinContent);
    }
    allEventsHist->Draw("HIST");
    twoElectronEventsHist->Draw("HIST E SAME");
    TLegend* legend = new TLegend(0.40, 0.44, 0.60, 0.56);
    legend->SetTextSize(0.028);
    legend->AddEntry(allEventsHist, allEventsLabel.c_str(), "l");
    legend->AddEntry(twoElectronEventsHist, twoElectronEventsLabel.c_str(), "l");
    legend->Draw();
  };

  // Helper for 2D color plots.  The right margin leaves room for the color scale,
  // and grid lines make the tightened origin z window easier to read.
  auto drawColorMap = [](TH2F* histogram) {
    if (gPad != nullptr)
    {
      gPad->SetRightMargin(0.15);
      gPad->SetGridx();
      gPad->SetGridy();
    }
    histogram->SetStats(false);
    histogram->Draw("COLZ");
  };

  //----------------------------------------------------------------------------
  // PDF 1: selected-event trkmcsim momentum overlay
  //
  // This plot compares all valid trkmcsim electrons vs rank-0 trkmcsim electrons,
  // but only inside selected two-electron events.
  //----------------------------------------------------------------------------
  TCanvas* cMCTMomentum = new TCanvas("cMCTMomentum", "Monte Carlo truth: trkmcsim momentum", 900, 700);

  // Style the "all electrons" histogram as a gray-filled black outline.
  hMCTAllElectronMomentum->SetLineColor(kBlack);
  hMCTAllElectronMomentum->SetLineWidth(2);
  hMCTAllElectronMomentum->SetFillColor(kGray);
  hMCTAllElectronMomentum->SetFillStyle(3004);

  // Style the rank-0 subset as a red line overlaid on top.
  hMCTRank0ElectronMomentum->SetLineColor(kRed);
  hMCTRank0ElectronMomentum->SetLineWidth(2);

  // Make room above the tallest bin for the legend and error bars.
  const double maxMomentumBinContent = max(hMCTAllElectronMomentum->GetMaximum(),
                                           hMCTRank0ElectronMomentum->GetMaximum());
  if (maxMomentumBinContent > 0.0)
  {
    hMCTAllElectronMomentum->SetMaximum(1.15 * maxMomentumBinContent);
  }
  hMCTAllElectronMomentum->Draw("HIST");
  hMCTRank0ElectronMomentum->Draw("HIST E SAME");

  // Centered, smaller legend.
  TLegend* momentumLegend = new TLegend(0.40, 0.44, 0.60, 0.56);
  momentumLegend->SetTextSize(0.028);
  momentumLegend->AddEntry(hMCTAllElectronMomentum, "All electrons", "l");
  momentumLegend->AddEntry(hMCTRank0ElectronMomentum, "Rank 0 electrons", "l");
  momentumLegend->Draw();
  if (mcTruthDirectory != nullptr)
  {
    mcTruthDirectory->cd();
  }

  // Save both into the ROOT file and into the Plots/ directory as a PDF.
  cMCTMomentum->Write();
  const string momentumPdfName = plotsDirectory + "/twoElectronEvents_" + generatorName + "_MonteCarloTruthMomentum.pdf";
  cMCTMomentum->SaveAs(momentumPdfName.c_str());

  //----------------------------------------------------------------------------
  // PDF 2: trkmcsim momentum event-selection comparison
  //
  // Left pad: all valid trkmcsim electrons.
  // Right pad: rank-0 trkmcsim electrons.
  //
  // Each pad overlays all events against selected two-electron events.
  //----------------------------------------------------------------------------
  TCanvas* cMCTMomentumSelectionComparison = new TCanvas(
    "cMCTMomentumSelectionComparison",
    "Monte Carlo truth momentum: all events vs selected two-electron events",
    1200, 650);
  cMCTMomentumSelectionComparison->Divide(2, 1);
  cMCTMomentumSelectionComparison->cd(1);
  drawEventSelectionComparison(
    hCompareMCTAllElectronMomentumAllEvents,
    hCompareMCTAllElectronMomentumTwoElectronEvents,
    "All events",
    "Selected two-electron events");
  cMCTMomentumSelectionComparison->cd(2);
  drawEventSelectionComparison(
    hCompareMCTRank0ElectronMomentumAllEvents,
    hCompareMCTRank0ElectronMomentumTwoElectronEvents,
    "All events",
    "Selected two-electron events");
  if (eventSelectionComparisonDirectory != nullptr)
  {
    eventSelectionComparisonDirectory->cd();
  }
  cMCTMomentumSelectionComparison->Write();
  const string mcTruthMomentumComparisonPdfName =
    plotsDirectory + "/twoElectronEvents_" + generatorName + "_MonteCarloTruthMomentumEventSelectionComparison.pdf";
  cMCTMomentumSelectionComparison->SaveAs(mcTruthMomentumComparisonPdfName.c_str());

  //----------------------------------------------------------------------------
  // PDF 3: selected-event trkmcsim angular histograms and maps
  //
  // Top row: all valid electron matches in selected events.
  // Middle row: rank-0 electron matches in selected events.
  // Bottom-left: theta of rank-0 electron 1 vs theta of rank-0 electron 2.
  //----------------------------------------------------------------------------
  TCanvas* cMCTAngles = new TCanvas("cMCTAngles", "Monte Carlo truth: trkmcsim angles", 1500, 1200);
  cMCTAngles->Divide(3, 3);

  hMCTAllElectronTheta->SetLineColor(kBlack);
  hMCTAllElectronTheta->SetLineWidth(2);
  hMCTAllElectronPhi->SetLineColor(kBlack);
  hMCTAllElectronPhi->SetLineWidth(2);
  hMCTRank0ElectronTheta->SetLineColor(kRed);
  hMCTRank0ElectronTheta->SetLineWidth(2);
  hMCTRank0ElectronPhi->SetLineColor(kRed);
  hMCTRank0ElectronPhi->SetLineWidth(2);

  cMCTAngles->cd(1);
  hMCTAllElectronTheta->Draw("HIST E");
  cMCTAngles->cd(2);
  hMCTAllElectronPhi->Draw("HIST E");
  cMCTAngles->cd(3);
  drawColorMap(hMCTAllElectronCosThetaVsPhi);

  cMCTAngles->cd(4);
  hMCTRank0ElectronTheta->Draw("HIST E");
  cMCTAngles->cd(5);
  hMCTRank0ElectronPhi->Draw("HIST E");
  cMCTAngles->cd(6);
  drawColorMap(hMCTRank0ElectronCosThetaVsPhi);

  cMCTAngles->cd(7);
  drawColorMap(hMCTRank0ElectronTheta1VsTheta2);
  if (mcTruthAngleDirectory != nullptr)
  {
    mcTruthAngleDirectory->cd();
  }
  cMCTAngles->Write();
  const string anglesPdfName = plotsDirectory + "/twoElectronEvents_" + generatorName + "_MonteCarloTruthAngles.pdf";
  cMCTAngles->SaveAs(anglesPdfName.c_str());

  //----------------------------------------------------------------------------
  // PDF 4: selected-event reconstructed trksegs momentum overlay
  //
  // This is the front/middle/back reconstructed momentum comparison for the
  // selected two-electron events only.  The linear-scale view starts at
  // 20 MeV/c to focus on the populated reconstructed-track region.
  //----------------------------------------------------------------------------
  TCanvas* cRecoTrkSegMomentum = new TCanvas(
    "cRecoTrkSegMomentum",
    "Reconstructed trksegs momentum at tracker surfaces",
    900, 700);

  // Give each tracker surface a distinct color.
  hRecoTrkSegMomentumFront->SetLineColor(kBlue + 1);
  hRecoTrkSegMomentumFront->SetLineWidth(2);
  hRecoTrkSegMomentumMiddle->SetLineColor(kGreen + 2);
  hRecoTrkSegMomentumMiddle->SetLineWidth(2);
  hRecoTrkSegMomentumBack->SetLineColor(kOrange + 7);
  hRecoTrkSegMomentumBack->SetLineWidth(2);

  // The front histogram is drawn first, so set its y-axis maximum high enough
  // for all three overlaid histograms.
  const double maxRecoTrkSegMomentumBinContent = max(
    hRecoTrkSegMomentumFront->GetMaximum(),
    max(hRecoTrkSegMomentumMiddle->GetMaximum(), hRecoTrkSegMomentumBack->GetMaximum()));
  if (maxRecoTrkSegMomentumBinContent > 0.0)
  {
    hRecoTrkSegMomentumFront->SetMaximum(1.15 * maxRecoTrkSegMomentumBinContent);
  }
  hRecoTrkSegMomentumFront->GetXaxis()->SetRangeUser(recoTrkSegMomentumLinearXMin, momentumMax);
  hRecoTrkSegMomentumFront->Draw("HIST E");
  hRecoTrkSegMomentumMiddle->Draw("HIST E SAME");
  hRecoTrkSegMomentumBack->Draw("HIST E SAME");

  // Centered, smaller legend for the three tracker surfaces.
  TLegend* recoTrkSegMomentumLegend = new TLegend(0.40, 0.44, 0.60, 0.56);
  recoTrkSegMomentumLegend->SetTextSize(0.028);
  recoTrkSegMomentumLegend->AddEntry(hRecoTrkSegMomentumFront, "TT_Front", "l");
  recoTrkSegMomentumLegend->AddEntry(hRecoTrkSegMomentumMiddle, "TT_Mid", "l");
  recoTrkSegMomentumLegend->AddEntry(hRecoTrkSegMomentumBack, "TT_Back", "l");
  recoTrkSegMomentumLegend->Draw();
  if (recoTrkSegMomentumDirectory != nullptr)
  {
    recoTrkSegMomentumDirectory->cd();
  }
  cRecoTrkSegMomentum->Write();
  const string recoTrkSegMomentumPdfName =
    plotsDirectory + "/twoElectronEvents_" + generatorName + "_RecoTrkSegMomentum.pdf";
  cRecoTrkSegMomentum->SaveAs(recoTrkSegMomentumPdfName.c_str());

  //----------------------------------------------------------------------------
  // PDF 5: selected-event reconstructed trksegs momentum overlay with log y
  //
  // This companion plot keeps the full 0-65 MeV/c momentum range and uses a
  // logarithmic y-axis so low-statistics tails are easier to see.
  //----------------------------------------------------------------------------
  TCanvas* cRecoTrkSegMomentumLogY = new TCanvas(
    "cRecoTrkSegMomentumLogY",
    "Reconstructed trksegs momentum at tracker surfaces, log y",
    900, 700);
  cRecoTrkSegMomentumLogY->SetLogy();

  hRecoTrkSegMomentumFront->GetXaxis()->SetRangeUser(momentumMin, momentumMax);
  if (maxRecoTrkSegMomentumBinContent > 0.0)
  {
    hRecoTrkSegMomentumFront->SetMinimum(0.5);
    hRecoTrkSegMomentumFront->SetMaximum(10.0 * maxRecoTrkSegMomentumBinContent);
  }
  hRecoTrkSegMomentumFront->Draw("HIST E");
  hRecoTrkSegMomentumMiddle->Draw("HIST E SAME");
  hRecoTrkSegMomentumBack->Draw("HIST E SAME");

  TLegend* recoTrkSegMomentumLogYLegend = new TLegend(0.40, 0.44, 0.60, 0.56);
  recoTrkSegMomentumLogYLegend->SetTextSize(0.028);
  recoTrkSegMomentumLogYLegend->AddEntry(hRecoTrkSegMomentumFront, "TT_Front", "l");
  recoTrkSegMomentumLogYLegend->AddEntry(hRecoTrkSegMomentumMiddle, "TT_Mid", "l");
  recoTrkSegMomentumLogYLegend->AddEntry(hRecoTrkSegMomentumBack, "TT_Back", "l");
  recoTrkSegMomentumLogYLegend->Draw();
  if (recoTrkSegMomentumDirectory != nullptr)
  {
    recoTrkSegMomentumDirectory->cd();
  }
  cRecoTrkSegMomentumLogY->Write();
  const string recoTrkSegMomentumLogYPdfName =
    plotsDirectory + "/twoElectronEvents_" + generatorName + "_RecoTrkSegMomentumLogY.pdf";
  cRecoTrkSegMomentumLogY->SaveAs(recoTrkSegMomentumLogYPdfName.c_str());

  //----------------------------------------------------------------------------
  // PDF 6: reconstructed trksegs event-selection comparison
  //
  // This canvas has six pads:
  //   top row    = all electron-matched tracks at TT_Front, TT_Mid, TT_Back
  //   bottom row = rank-0-electron-matched tracks at TT_Front, TT_Mid, TT_Back
  //
  // Each pad overlays all events against selected two-electron events.
  //----------------------------------------------------------------------------
  TCanvas* cRecoTrkSegMomentumSelectionComparison = new TCanvas(
    "cRecoTrkSegMomentumSelectionComparison",
    "Reconstructed trksegs momentum: all events vs selected two-electron events",
    1500, 900);
  cRecoTrkSegMomentumSelectionComparison->Divide(3, 2);
  for (int i_surface = 0; i_surface < 3; ++i_surface)
  {
    // Top row: all electron-matched tracks at this surface.
    cRecoTrkSegMomentumSelectionComparison->cd(i_surface + 1);
    drawEventSelectionComparison(
      hCompareRecoTrkSegAllElectronMomentumAllEvents[i_surface],
      hCompareRecoTrkSegAllElectronMomentumTwoElectronEvents[i_surface],
      "All events",
      "Selected two-electron events");

    // Bottom row: rank-0-electron-matched tracks at this surface.
    cRecoTrkSegMomentumSelectionComparison->cd(i_surface + 4);
    drawEventSelectionComparison(
      hCompareRecoTrkSegRank0ElectronMomentumAllEvents[i_surface],
      hCompareRecoTrkSegRank0ElectronMomentumTwoElectronEvents[i_surface],
      "All events",
      "Selected two-electron events");
  }
  if (eventSelectionComparisonDirectory != nullptr)
  {
    eventSelectionComparisonDirectory->cd();
  }
  cRecoTrkSegMomentumSelectionComparison->Write();
  const string recoTrkSegMomentumComparisonPdfName =
    plotsDirectory + "/twoElectronEvents_" + generatorName + "_RecoTrkSegMomentumEventSelectionComparison.pdf";
  cRecoTrkSegMomentumSelectionComparison->SaveAs(recoTrkSegMomentumComparisonPdfName.c_str());

  //----------------------------------------------------------------------------
  // PDF 7: selected-event trkmcsim origin 1D histograms
  //
  // Six pads:
  //   top row    = all valid electron origins x, y, z
  //   bottom row = rank-0 electron origins x, y, z
  //----------------------------------------------------------------------------
  TCanvas* cMCTOrigins = new TCanvas("cMCTOrigins", "Monte Carlo truth: trkmcsim origins", 1200, 800);
  cMCTOrigins->Divide(3, 2);

  // Top row: all valid electrons.
  cMCTOrigins->cd(1);
  hMCTAllElectronOriginX->Draw("HIST E");
  cMCTOrigins->cd(2);
  hMCTAllElectronOriginY->Draw("HIST E");
  cMCTOrigins->cd(3);
  hMCTAllElectronOriginZ->Draw("HIST E");

  // Bottom row: rank-0 electrons.
  cMCTOrigins->cd(4);
  hMCTRank0ElectronOriginX->Draw("HIST E");
  cMCTOrigins->cd(5);
  hMCTRank0ElectronOriginY->Draw("HIST E");
  cMCTOrigins->cd(6);
  hMCTRank0ElectronOriginZ->Draw("HIST E");
  if (mcTruthDirectory != nullptr)
  {
    mcTruthDirectory->cd();
  }
  cMCTOrigins->Write();
  const string originsPdfName = plotsDirectory + "/twoElectronEvents_" + generatorName + "_MonteCarloTruthOrigins.pdf";
  cMCTOrigins->SaveAs(originsPdfName.c_str());

  //----------------------------------------------------------------------------
  // PDF 8: selected-event trkmcsim origin 2D maps
  //
  // Six pads:
  //   top row    = all valid electron XY, XZ, YZ origin maps
  //   bottom row = rank-0 electron XY, XZ, YZ origin maps
  //
  // COLZ draws a colored density map with a z-color scale.
  //----------------------------------------------------------------------------
  // The z span is still about four times the x/y span, so a wider canvas helps
  // the XZ/YZ pads without making the XY pads unreadably small.
  TCanvas* cMCTOriginMaps = new TCanvas("cMCTOriginMaps", "Monte Carlo truth: trkmcsim origin maps", 1800, 850);
  cMCTOriginMaps->Divide(3, 2);

  // Top row: all valid electrons.
  cMCTOriginMaps->cd(1);
  drawColorMap(hMCTAllElectronOriginXY);
  cMCTOriginMaps->cd(2);
  drawColorMap(hMCTAllElectronOriginXZ);
  cMCTOriginMaps->cd(3);
  drawColorMap(hMCTAllElectronOriginYZ);

  // Bottom row: rank-0 electrons.
  cMCTOriginMaps->cd(4);
  drawColorMap(hMCTRank0ElectronOriginXY);
  cMCTOriginMaps->cd(5);
  drawColorMap(hMCTRank0ElectronOriginXZ);
  cMCTOriginMaps->cd(6);
  drawColorMap(hMCTRank0ElectronOriginYZ);
  if (mcTruthDirectory != nullptr)
  {
    mcTruthDirectory->cd();
  }
  cMCTOriginMaps->Write();
  const string originMapsPdfName = plotsDirectory + "/twoElectronEvents_" + generatorName + "_MonteCarloTruthOriginMaps.pdf";
  cMCTOriginMaps->SaveAs(originMapsPdfName.c_str());

  //----------------------------------------------------------------------------
  // PDF 9: selected-event rank-0 trkmcsim origin 1D histograms
  //
  // This version omits the all-electron row and shows only the rank-0 electrons
  // that define the selected two-electron events.
  //----------------------------------------------------------------------------
  TCanvas* cMCTRank0Origins = new TCanvas(
    "cMCTRank0Origins",
    "Monte Carlo truth: rank-0 trkmcsim origins",
    1200, 450);
  cMCTRank0Origins->Divide(3, 1);
  cMCTRank0Origins->cd(1);
  hMCTRank0ElectronOriginX->Draw("HIST E");
  cMCTRank0Origins->cd(2);
  hMCTRank0ElectronOriginY->Draw("HIST E");
  cMCTRank0Origins->cd(3);
  hMCTRank0ElectronOriginZ->Draw("HIST E");
  if (mcTruthDirectory != nullptr)
  {
    mcTruthDirectory->cd();
  }
  cMCTRank0Origins->Write();
  const string rank0OriginsPdfName =
    plotsDirectory + "/twoElectronEvents_" + generatorName + "_MonteCarloTruthRank0Origins.pdf";
  cMCTRank0Origins->SaveAs(rank0OriginsPdfName.c_str());

  //----------------------------------------------------------------------------
  // PDF 10: selected-event rank-0 trkmcsim origin 2D maps
  //
  // This version omits the all-electron row and shows only the rank-0 XY, XZ,
  // and YZ origin maps.
  //----------------------------------------------------------------------------
  TCanvas* cMCTRank0OriginMaps = new TCanvas(
    "cMCTRank0OriginMaps",
    "Monte Carlo truth: rank-0 trkmcsim origin maps",
    1800, 500);
  cMCTRank0OriginMaps->Divide(3, 1);
  cMCTRank0OriginMaps->cd(1);
  drawColorMap(hMCTRank0ElectronOriginXY);
  cMCTRank0OriginMaps->cd(2);
  drawColorMap(hMCTRank0ElectronOriginXZ);
  cMCTRank0OriginMaps->cd(3);
  drawColorMap(hMCTRank0ElectronOriginYZ);
  if (mcTruthDirectory != nullptr)
  {
    mcTruthDirectory->cd();
  }
  cMCTRank0OriginMaps->Write();
  const string rank0OriginMapsPdfName =
    plotsDirectory + "/twoElectronEvents_" + generatorName + "_MonteCarloTruthRank0OriginMaps.pdf";
  cMCTRank0OriginMaps->SaveAs(rank0OriginMapsPdfName.c_str());

  // Stop the plot timer after the last PDF is saved.
  plotTimer.Stop();

  // Close the ROOT file after all histograms and canvases are written.
  histogramFile.Close();

  //----------------------------------------------------------------------------
  // Cleanup
  //
  // The macro created histograms and canvases with new, so it deletes them before
  // returning.  This keeps an interactive ROOT session from holding onto lots of
  // objects after the macro finishes.
  //----------------------------------------------------------------------------

  // Delete canvases first.  They are the GUI/drawing objects.
  delete cMCTMomentum;
  delete cMCTMomentumSelectionComparison;
  delete cMCTAngles;
  delete cRecoTrkSegMomentum;
  delete cRecoTrkSegMomentumLogY;
  delete cRecoTrkSegMomentumSelectionComparison;
  delete cMCTOrigins;
  delete cMCTOriginMaps;
  delete cMCTRank0Origins;
  delete cMCTRank0OriginMaps;

  // Delete trkmcsim momentum and comparison histograms.
  delete hMCTAllElectronMomentum;
  delete hMCTRank0ElectronMomentum;
  delete hCompareMCTAllElectronMomentumAllEvents;
  delete hCompareMCTAllElectronMomentumTwoElectronEvents;
  delete hCompareMCTRank0ElectronMomentumAllEvents;
  delete hCompareMCTRank0ElectronMomentumTwoElectronEvents;

  // Delete trkmcsim origin histograms and 2D maps.
  delete hMCTAllElectronOriginX;
  delete hMCTAllElectronOriginY;
  delete hMCTAllElectronOriginZ;
  delete hMCTRank0ElectronOriginX;
  delete hMCTRank0ElectronOriginY;
  delete hMCTRank0ElectronOriginZ;
  delete hMCTAllElectronOriginXY;
  delete hMCTAllElectronOriginXZ;
  delete hMCTAllElectronOriginYZ;
  delete hMCTRank0ElectronOriginXY;
  delete hMCTRank0ElectronOriginXZ;
  delete hMCTRank0ElectronOriginYZ;

  // Delete trkmcsim angular histograms and maps.
  delete hMCTAllElectronTheta;
  delete hMCTAllElectronPhi;
  delete hMCTAllElectronCosThetaVsPhi;
  delete hMCTRank0ElectronTheta;
  delete hMCTRank0ElectronPhi;
  delete hMCTRank0ElectronCosThetaVsPhi;
  delete hMCTRank0ElectronTheta1VsTheta2;

  // Delete the original selected-event reconstructed trksegs histograms.
  delete hRecoTrkSegMomentumFront;
  delete hRecoTrkSegMomentumMiddle;
  delete hRecoTrkSegMomentumBack;

  // Delete the comparison trksegs histograms for all three surfaces.
  for (int i_surface = 0; i_surface < 3; ++i_surface)
  {
    delete hCompareRecoTrkSegAllElectronMomentumAllEvents[i_surface];
    delete hCompareRecoTrkSegAllElectronMomentumTwoElectronEvents[i_surface];
    delete hCompareRecoTrkSegRank0ElectronMomentumAllEvents[i_surface];
    delete hCompareRecoTrkSegRank0ElectronMomentumTwoElectronEvents[i_surface];
  }

  // Restore ROOT settings to the values they had before this macro changed them.
  TH1::AddDirectory(oldAddDirectoryStatus);
  gROOT->SetBatch(wasBatchMode);

  // Final terminal output: file names and timing information.
  cout << "Writing and saving plots took "
       << plotTimer.RealTime() << " seconds of wall time." << endl;
  cout << "Wrote Monte Carlo truth histograms to " << histogramFileName << endl;
  cout << "Wrote Monte Carlo truth PDF plots to " << momentumPdfName
       << ", " << anglesPdfName
       << ", " << originsPdfName
       << ", " << originMapsPdfName
       << ", " << rank0OriginsPdfName
       << ", and " << rank0OriginMapsPdfName << endl;
  cout << "Wrote reconstructed trksegs momentum PDF plots to "
       << recoTrkSegMomentumPdfName
       << " and " << recoTrkSegMomentumLogYPdfName << endl;
  cout << "Wrote event-selection comparison PDF plots to "
       << mcTruthMomentumComparisonPdfName
       << " and " << recoTrkSegMomentumComparisonPdfName << endl;

}//end twoElectronEventCutterHistogrammer method

