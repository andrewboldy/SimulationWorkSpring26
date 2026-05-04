
//----------------------------------------------------------------------------------
//
// plotRecoCommonOriginCandidates.C
// Written by Andrew Boldy University of South Carolina, 2026
// Assisted by Codex
//
// Purpose:
// This macro performs a first-pass study of whether two reconstructed tracks
// could have originated from a common point upstream of the tracker.
//
// Important idea:
// We are NOT doing a full two-track KinKal vertex fit here.
// Instead, we use information that is already stored in the EventNtuple:
//
//   1. Reconstructed track segments (`trksegs`)
//   2. Reconstructed momentum vectors at those segments
//   3. Known tracker surfaces such as TT_Front, TT_Mid, and TT_Back
//
// For each selected downstream electron reconstructed track, we build a simple
// local line approximation near the tracker entrance.  Then, for each pair of
// selected tracks in the event, we use the Mu2e helper class `TwoLinePCA_XYZ` to
// compute the distance of closest approach (DCA) between those two extrapolated
// lines.
//
// The midpoint of that closest-approach segment is used here as a simple
// "candidate common origin" estimate.
//
// This is best viewed as:
//   - a learning / exploration tool
//   - a screening tool for possible common-origin pairs
//   - a stepping stone toward a more rigorous vertexing study
//
// It is NOT yet:
//   - a precision external-vertex fitter
//   - a replacement for a full curved-trajectory vertex fit
//
// Example:
// root -l -b -q 'CreatedCode/plotRecoCommonOriginCandidates.C++("B2BCeEndpoint","myfiles.txt")'
//
// Suggested outputs:
//   - pair DCA histogram
//   - candidate vertex Z histogram
//   - candidate vertex radius histogram
//   - 2D map of candidate origin X vs Z
//   - comparison of all selected downstream-e pairs vs pairs whose closest-
//     approach point lies upstream of the tracker entrance for both tracks
//
//----------------------------------------------------------------------------------

//Standard Inclusions
#include <algorithm>
#include <cctype>
#include <cmath>
#include <fstream>
#include <iostream>
#include <limits>
#include <memory>
#include <string>
#include <vector>

//CERN ROOT Inclusions
#include "TCanvas.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TLegend.h"
#include "TLine.h"

//Mu2e Inclusions
#include "Offline/DataProducts/inc/GenVector.hh"
#include "Offline/DataProducts/inc/PDGCode.hh"
#include "Offline/DataProducts/inc/SurfaceId.hh"

#include "EventNtuple/rooutil/inc/RooUtil.hh"
#include "EventNtuple/rooutil/inc/common_cuts.hh"

//Namespace
using namespace std;
using namespace rooutil;

namespace {

  //----------------------------------------------------------------------------
  // Small helper structure used to store the information needed to represent
  // one reconstructed track as a line near the tracker entrance.
  //
  // Why a line?
  // A true Mu2e track in the magnetic field is curved, but over a short region
  // near a chosen surface it is often useful to approximate it with a local
  // tangent line.  That is enough to perform an initial "do these two tracks
  // seem to point back to the same place?" study.
  //----------------------------------------------------------------------------
  struct EntranceLine {
    bool valid = false;                 // Did we successfully build a line for this track?
    XYZVectorF point = XYZVectorF();    // A point on the line, usually the TT_Front segment position
    XYZVectorF dir = XYZVectorF();      // Unit direction vector of the local line approximation
    XYZVectorF refPos = XYZVectorF();   // Same as point, saved separately for readability in later logic
    int pdg = 0;                        // Track fit hypothesis PDG code
    int trackIndex = -1;                // Which track in the event this came from
    bool rank0Electron = false;         // True if the track's MC rank-0 match is an electron
  };

  //----------------------------------------------------------------------------
  // Small helper structure that stores the result of a closest-approach
  // calculation between two 3D lines.
  //
  // This replaces the direct dependency on Mu2eUtilities/TwoLinePCA_XYZ.
  // The reason for doing that is practical: ACLiC can compile the macro itself,
  // but in your environment it did not automatically link the external library
  // object that contains the TwoLinePCA_XYZ implementation.
  //
  // By implementing the math directly here, the macro becomes much easier to
  // run as a standalone ROOT analysis macro.
  //----------------------------------------------------------------------------
  struct ClosestApproachResult {
    bool valid = false;                 // Did the calculation succeed?
    bool nearlyParallel = false;        // Were the lines treated as nearly parallel?
    double dca = -1.0;                  // 3D distance of closest approach
    XYZVectorF point1 = XYZVectorF();   // Closest point on line 1
    XYZVectorF point2 = XYZVectorF();   // Closest point on line 2
  };

  //----------------------------------------------------------------------------
  // Clean strings before using them in output file names.
  //----------------------------------------------------------------------------
  string sanitizeName(const string& input) {
    string output = input;
    for (char& c : output) {
      if (!(isalnum(static_cast<unsigned char>(c)) || c == '_' || c == '-')) {
        c = '_';
      }
    }
    return output;
  }

  //----------------------------------------------------------------------------
  // Decide whether the input string is a ROOT file path or a plain-text file
  // list.
  //
  // This macro was originally written assuming a text file list, but in
  // practice you ran it directly on a single ROOT file.  That is perfectly
  // reasonable, so we support both cases here.
  //----------------------------------------------------------------------------
  bool looksLikeRootFile(const string& input) {
    if (input.size() < 5) {
      return false;
    }
    return input.substr(input.size() - 5) == ".root";
  }

  //----------------------------------------------------------------------------
  // Count how many input files are represented by the user argument.
  //
  // If the user passed a single ROOT file, the answer is just 1.
  // If the user passed a text file list, count non-empty lines.
  //----------------------------------------------------------------------------
  int countInputFiles(const string& inputPath) {
    if (looksLikeRootFile(inputPath)) {
      return 1;
    }

    ifstream listOfFiles(inputPath);
    if (!listOfFiles.is_open()) {
      return 0;
    }

    int fileCount = 0;
    string line;
    while (getline(listOfFiles, line)) {
      if (!line.empty()) {
        ++fileCount;
      }
    }
    return fileCount;
  }

  //----------------------------------------------------------------------------
  // Helper to decide whether a reconstructed track segment is at a specific
  // named surface.
  //
  // We only use segments that really have reconstructed information attached.
  //----------------------------------------------------------------------------
  bool hasRecoSegmentAtSurface(TrackSegment& segment, int sid) {
    return has_reco_step(segment) && segment.trkseg->sid == sid;
  }

  //----------------------------------------------------------------------------
  // Determine whether a reconstructed track is matched to a rank-0 electron.
  //
  // In RooUtil / EventNtuple language:
  //   - `is_track_particle` means "the MC particle with the most hits on the
  //     track" and also requires `rank == 0`
  //   - we then additionally require that this best-matched MC particle is an
  //     electron
  //
  // This is useful because it suppresses many tracker-ionization byproducts and
  // lets us study a cleaner subset of reconstructed electron-like tracks.
  //----------------------------------------------------------------------------
  bool isRank0ElectronTrack(Track& track) {
    auto matchedParticles = track.GetMCParticles(is_track_particle);
    if (matchedParticles.empty() || matchedParticles.at(0).mcsim == nullptr) {
      return false;
    }

    const auto pdg = static_cast<mu2e::PDGCode::type>(matchedParticles.at(0).mcsim->pdg);
    return pdg == mu2e::PDGCode::e_minus || pdg == mu2e::PDGCode::e_plus;
  }

  //----------------------------------------------------------------------------
  // Select the reconstructed fit hypothesis used by this study.
  //
  // We want downstream electron fits only.  In EventNtuple, the electron part is
  // the fit particle hypothesis (`trk->pdg == e_minus`).  The downstream part is
  // inferred from the reconstructed momentum direction at a tracker surface.
  //
  // Prefer TT_Front because this macro anchors its line approximation there.
  // Fall back to TT_Mid and TT_Back so the selector remains usable if a given
  // ntuple does not contain every sampled tracker surface.
  //----------------------------------------------------------------------------
  bool isDownstreamElectronTrack(Track& track) {
    if (!is_e_minus(track)) {
      return false;
    }

    const vector<int> directionSurfaces = {
      mu2e::SurfaceIdDetail::TT_Front,
      mu2e::SurfaceIdDetail::TT_Mid,
      mu2e::SurfaceIdDetail::TT_Back
    };

    for (const int sid : directionSurfaces) {
      auto segments = track.GetSegments([sid](TrackSegment& segment) {
        return hasRecoSegmentAtSurface(segment, sid);
      });

      if (segments.empty()) {
        continue;
      }

      return any_of(segments.begin(), segments.end(), [](TrackSegment& segment) {
        return segment.trkseg != nullptr && segment.trkseg->mom.z() > 0.0;
      });
    }

    return false;
  }

  //----------------------------------------------------------------------------
  // Compute the closest approach between two infinite 3D lines.
  //
  // Inputs:
  //   line 1: p1 + s * t1
  //   line 2: p2 + u * t2
  //
  // where:
  //   p1, p2 are points on the two lines
  //   t1, t2 are direction vectors
  //
  // The method used here is the standard analytic minimum-distance solution.
  // If the lines are nearly parallel, we fall back to a simpler projection
  // method so that the macro remains numerically stable.
  //
  // This is exactly the kind of "geometry helper" that is useful to keep
  // commented and local while you are still learning the analysis chain.
  //----------------------------------------------------------------------------
  ClosestApproachResult computeClosestApproach(const XYZVectorF& p1,
                                               const XYZVectorF& t1Input,
                                               const XYZVectorF& p2,
                                               const XYZVectorF& t2Input,
                                               double parallelCut = 1.0e-8) {
    ClosestApproachResult result;

    const double t1Mag = t1Input.R();
    const double t2Mag = t2Input.R();
    if (t1Mag < 1.0e-12 || t2Mag < 1.0e-12) {
      return result;
    }

    const XYZVectorF t1 = (1.0 / t1Mag) * t1Input;
    const XYZVectorF t2 = (1.0 / t2Mag) * t2Input;
    const XYZVectorF w0 = p1 - p2;

    const double a = t1.Dot(t1);   // = 1 for unit vector, but written explicitly for clarity
    const double b = t1.Dot(t2);
    const double c = t2.Dot(t2);   // = 1 for unit vector
    const double d = t1.Dot(w0);
    const double e = t2.Dot(w0);
    const double denom = a * c - b * b;

    double s = 0.0;
    double u = 0.0;

    if (std::fabs(denom) < parallelCut) {
      // The lines are nearly parallel.
      //
      // In that limit, there is not a unique shortest connecting segment in a
      // numerically stable sense.  A simple and useful fallback is:
      //   - keep point1 at p1
      //   - project the vector from p2 to p1 onto line 2 to get point2
      //
      // That still gives us a meaningful separation estimate.
      result.nearlyParallel = true;
      s = 0.0;
      u = e / c;
    } else {
      // Standard analytic solution for the closest points on two skew lines.
      s = (b * e - c * d) / denom;
      u = (a * e - b * d) / denom;
    }

    result.point1 = p1 + static_cast<float>(s) * t1;
    result.point2 = p2 + static_cast<float>(u) * t2;
    result.dca = (result.point1 - result.point2).R();
    result.valid = true;
    return result;
  }

  //----------------------------------------------------------------------------
  // Return the first reconstructed segment at the requested surface.
  //
  // We keep this helper simple and explicit because it is easier to read than
  // repeating the same lambda-heavy logic throughout the macro.
  //----------------------------------------------------------------------------
  TrackSegment* findFirstRecoSegmentAtSurface(Track& track, int sid) {
    auto segments = track.GetSegments([sid](TrackSegment& segment) {
      return hasRecoSegmentAtSurface(segment, sid);
    });

    if (segments.empty()) {
      return nullptr;
    }

    // RooUtil returns TrackSegment wrapper objects by value.
    // Each wrapper still points at the underlying ntuple object, so it is safe
    // for us to inspect the stored `trkseg` payload through the returned copy.
    //
    // We cannot return `&segments.at(0)` here because that would point to a
    // temporary vector element.  Instead we just use this helper only for
    // extracting data immediately in the calling code.
    return new TrackSegment(segments.at(0));
  }

  //----------------------------------------------------------------------------
  // Build a local line approximation for a reconstructed track.
  //
  // Strategy:
  //   1. Prefer the reconstructed track segment at TT_Front as the point where
  //      we anchor the line.
  //   2. Use a second reconstructed segment downstream (prefer TT_Mid, then
  //      TT_Back) to define the line direction geometrically.
  //   3. If that geometric direction is unavailable or degenerate, fall back to
  //      the momentum vector stored at TT_Front.
  //
  // Why this approach?
  // The segment positions let us estimate the actual reconstructed path
  // direction in a way that is easy to understand.  The momentum fallback keeps
  // the code useful even when only one reconstructed tracker surface is saved.
  //----------------------------------------------------------------------------
  EntranceLine buildEntranceLine(Track& track, int trackIndex) {
    EntranceLine result;
    result.trackIndex = trackIndex;
    result.pdg = track.trk->pdg;
    result.rank0Electron = isRank0ElectronTrack(track);

    unique_ptr<TrackSegment> front(findFirstRecoSegmentAtSurface(track, mu2e::SurfaceIdDetail::TT_Front));
    if (!front || front->trkseg == nullptr) {
      return result;
    }

    result.point = front->trkseg->pos;
    result.refPos = front->trkseg->pos;

    // Try to define the local direction from two reconstructed positions.
    // This is usually the clearest approximation to "which way is the track
    // heading as it passes through the tracker entrance?"
    XYZVectorF direction = XYZVectorF(0.f, 0.f, 0.f);

    unique_ptr<TrackSegment> middle(findFirstRecoSegmentAtSurface(track, mu2e::SurfaceIdDetail::TT_Mid));
    unique_ptr<TrackSegment> back(findFirstRecoSegmentAtSurface(track, mu2e::SurfaceIdDetail::TT_Back));

    if (middle && middle->trkseg != nullptr) {
      direction = middle->trkseg->pos - front->trkseg->pos;
    } else if (back && back->trkseg != nullptr) {
      direction = back->trkseg->pos - front->trkseg->pos;
    }

    // If the direction from positions is too small, use the reconstructed
    // momentum vector at TT_Front instead.
    //
    // This fallback is nice because the segment already stores the fit momentum
    // at that surface.  It is still only a local approximation, but it carries
    // more physical meaning than inventing a direction from nothing.
    if (direction.R() < 1.0e-6) {
      direction = front->trkseg->mom;
    }

    if (direction.R() < 1.0e-6) {
      return result;
    }

    result.dir = (1.0 / direction.R()) * direction;
    result.valid = true;
    return result;
  }

  //----------------------------------------------------------------------------
  // Signed distance along the chosen line direction from the reference point to
  // some test point.
  //
  // Interpretation in this macro:
  // If the returned value is negative, then the test point lies "upstream"
  // relative to the local momentum / track direction that we stored.
  //
  // That makes it a simple diagnostic for whether the candidate common-origin
  // point lies behind the tracker entrance, which is exactly what we want for a
  // backward extrapolation study.
  //----------------------------------------------------------------------------
  double signedDistanceAlongLine(const XYZVectorF& origin,
                                 const XYZVectorF& dir,
                                 const XYZVectorF& testPoint) {
    return (testPoint - origin).Dot(dir);
  }

} // namespace

void plotRecoCommonOriginCandidates(const string& sampleName, const string& filelist) {
  //----------------------------------------------------------------------------
  // Read the file list first.
  //
  // This is not strictly needed for the reconstruction logic, but it is useful
  // bookkeeping and matches the style of your existing ROOT macro.
  //----------------------------------------------------------------------------
  const bool inputIsRootFile = looksLikeRootFile(filelist);

  ifstream listOfFiles;
  if (!inputIsRootFile) {
    listOfFiles.open(filelist);
  }

  if (!inputIsRootFile && !listOfFiles.is_open()) {
    cerr << "Could not open file list: " << filelist << endl;
    return;
  }

  const int fileCount = countInputFiles(filelist);

  //----------------------------------------------------------------------------
  // Set up output names.
  //----------------------------------------------------------------------------
  const string safeSampleName = sanitizeName(sampleName);

  //----------------------------------------------------------------------------
  // Histograms.
  //
  // We keep both "all selected downstream-e track pairs" and a more restrictive
  // "upstream-only" category.  The upstream-only category asks whether the
  // estimated closest-approach point lies behind the tracker entrance for both
  // tracks.
  //
  // That is not yet a full vertexing requirement, but it is a physically useful
  // first filter for the kind of common-origin question you asked about.
  //----------------------------------------------------------------------------
  TH1F* hPairDCAAll = new TH1F(
      "hPairDCAAll",
      (sampleName + " Downstream e- Pair DCA from Backward Extrapolation;Pair DCA [mm];Entries").c_str(),
      100, 0.0, 500.0);

  TH1F* hPairDCAUpstream = new TH1F(
      "hPairDCAUpstream",
      (sampleName + " Downstream e- Pair DCA (Closest Approach Upstream of Tracker);Pair DCA [mm];Entries").c_str(),
      100, 0.0, 500.0);

  TH1F* hVertexZAll = new TH1F(
      "hVertexZAll",
      (sampleName + " Downstream e- Candidate Common-Origin Z;Candidate Z [mm];Entries").c_str(),
      120, -7000.0, 3000.0);

  TH1F* hVertexZUpstream = new TH1F(
      "hVertexZUpstream",
      (sampleName + " Downstream e- Candidate Common-Origin Z (Upstream Only);Candidate Z [mm];Entries").c_str(),
      120, -7000.0, 3000.0);

  TH1F* hVertexRAll = new TH1F(
      "hVertexRAll",
      (sampleName + " Downstream e- Candidate Common-Origin Radius;Candidate Radius [mm];Entries").c_str(),
      120, 0.0, 1000.0);

  TH1F* hVertexRUpstream = new TH1F(
      "hVertexRUpstream",
      (sampleName + " Downstream e- Candidate Common-Origin Radius (Upstream Only);Candidate Radius [mm];Entries").c_str(),
      120, 0.0, 1000.0);

  TH2F* hVertexXZAll = new TH2F(
      "hVertexXZAll",
      (sampleName + " Downstream e- Candidate Common-Origin Map;z [mm];x [mm]").c_str(),
      120, -7000.0, 3000.0,
      120, -1000.0, 1000.0);

  TH2F* hVertexXZUpstream = new TH2F(
      "hVertexXZUpstream",
      (sampleName + " Downstream e- Candidate Common-Origin Map (Upstream Only);z [mm];x [mm]").c_str(),
      120, -7000.0, 3000.0,
      120, -1000.0, 1000.0);

  //----------------------------------------------------------------------------
  // Additional histograms for the cleaner subset in which BOTH reconstructed
  // tracks are matched to rank-0 electrons.
  //
  // This is meant to reduce the flood of tracker-produced secondaries and focus
  // more directly on pairs whose best MC association is electron-like.
  //----------------------------------------------------------------------------
  TH1F* hPairDCARank0Electron = new TH1F(
      "hPairDCARank0Electron",
      (sampleName + " Pair DCA for Rank-0 Electron Track Pairs;Pair DCA [mm];Entries").c_str(),
      100, 0.0, 500.0);

  TH1F* hPairDCARank0ElectronUpstream = new TH1F(
      "hPairDCARank0ElectronUpstream",
      (sampleName + " Pair DCA for Rank-0 Electron Track Pairs (Upstream Only);Pair DCA [mm];Entries").c_str(),
      100, 0.0, 500.0);

  TH1F* hVertexZRank0Electron = new TH1F(
      "hVertexZRank0Electron",
      (sampleName + " Candidate Common-Origin Z for Rank-0 Electron Track Pairs;Candidate Z [mm];Entries").c_str(),
      120, -7000.0, 3000.0);

  TH1F* hVertexZRank0ElectronUpstream = new TH1F(
      "hVertexZRank0ElectronUpstream",
      (sampleName + " Candidate Common-Origin Z for Rank-0 Electron Track Pairs (Upstream Only);Candidate Z [mm];Entries").c_str(),
      120, -7000.0, 3000.0);

  TH1F* hVertexRRank0Electron = new TH1F(
      "hVertexRRank0Electron",
      (sampleName + " Candidate Common-Origin Radius for Rank-0 Electron Track Pairs;Candidate Radius [mm];Entries").c_str(),
      120, 0.0, 1000.0);

  TH1F* hVertexRRank0ElectronUpstream = new TH1F(
      "hVertexRRank0ElectronUpstream",
      (sampleName + " Candidate Common-Origin Radius for Rank-0 Electron Track Pairs (Upstream Only);Candidate Radius [mm];Entries").c_str(),
      120, 0.0, 1000.0);

  TH2F* hVertexXZRank0Electron = new TH2F(
      "hVertexXZRank0Electron",
      (sampleName + " Candidate Common-Origin Map for Rank-0 Electron Track Pairs;z [mm];x [mm]").c_str(),
      120, -7000.0, 3000.0,
      120, -1000.0, 1000.0);

  TH2F* hVertexXZRank0ElectronUpstream = new TH2F(
      "hVertexXZRank0ElectronUpstream",
      (sampleName + " Candidate Common-Origin Map for Rank-0 Electron Track Pairs (Upstream Only);z [mm];x [mm]").c_str(),
      120, -7000.0, 3000.0,
      120, -1000.0, 1000.0);

  hPairDCAAll->SetLineColor(kBlue + 1);
  hPairDCAAll->SetLineWidth(2);
  hPairDCAUpstream->SetLineColor(kRed + 1);
  hPairDCAUpstream->SetLineWidth(2);

  hVertexZAll->SetLineColor(kBlue + 1);
  hVertexZAll->SetLineWidth(2);
  hVertexZUpstream->SetLineColor(kRed + 1);
  hVertexZUpstream->SetLineWidth(2);

  hVertexRAll->SetLineColor(kBlue + 1);
  hVertexRAll->SetLineWidth(2);
  hVertexRUpstream->SetLineColor(kRed + 1);
  hVertexRUpstream->SetLineWidth(2);

  hPairDCARank0Electron->SetLineColor(kMagenta + 2);
  hPairDCARank0Electron->SetLineWidth(2);
  hPairDCARank0ElectronUpstream->SetLineColor(kOrange + 7);
  hPairDCARank0ElectronUpstream->SetLineWidth(2);

  hVertexZRank0Electron->SetLineColor(kMagenta + 2);
  hVertexZRank0Electron->SetLineWidth(2);
  hVertexZRank0ElectronUpstream->SetLineColor(kOrange + 7);
  hVertexZRank0ElectronUpstream->SetLineWidth(2);

  hVertexRRank0Electron->SetLineColor(kMagenta + 2);
  hVertexRRank0Electron->SetLineWidth(2);
  hVertexRRank0ElectronUpstream->SetLineColor(kOrange + 7);
  hVertexRRank0ElectronUpstream->SetLineWidth(2);

  //----------------------------------------------------------------------------
  // Open the EventNtuple through RooUtil.
  //----------------------------------------------------------------------------
  RooUtil util(filelist);
  const int numEvents = util.GetNEvents();

  //----------------------------------------------------------------------------
  // Counters printed at the end.
  //----------------------------------------------------------------------------
  int numTracksSeen = 0;
  int numUsableLines = 0;
  int numPairsTested = 0;
  int numUpstreamPairs = 0;
  int numRank0ElectronTracks = 0;
  int numRank0ElectronPairs = 0;
  int numRank0ElectronUpstreamPairs = 0;

  if (inputIsRootFile) {
    cout << "Reading a single ROOT file: " << filelist << endl;
    cout << "There are " << numEvents << " events in this ROOT file." << endl;
  } else {
    cout << "Read " << fileCount << " files from file list " << filelist << endl;
    cout << "There are " << numEvents << " events across that file list." << endl;
  }

  //----------------------------------------------------------------------------
  // Event loop.
  //----------------------------------------------------------------------------
  for (int i_event = 0; i_event < numEvents; ++i_event) {
    auto& event = util.GetEvent(i_event);

    //------------------------------------------------------------------------
    // Pull only downstream electron reconstructed fit hypotheses.
    //------------------------------------------------------------------------
    auto tracks = event.GetTracks(isDownstreamElectronTrack);

    vector<EntranceLine> lines;
    lines.reserve(tracks.size());

    //------------------------------------------------------------------------
    // Convert each reconstructed track into a local entrance-line
    // approximation.
    //------------------------------------------------------------------------
    for (size_t i_track = 0; i_track < tracks.size(); ++i_track) {
      ++numTracksSeen;
      EntranceLine entranceLine = buildEntranceLine(tracks.at(i_track), static_cast<int>(i_track));
      if (entranceLine.valid) {
        ++numUsableLines;
        if (entranceLine.rank0Electron) {
          ++numRank0ElectronTracks;
        }
        lines.push_back(entranceLine);
      }
    }

    //------------------------------------------------------------------------
    // Pair loop.
    //
    // For each pair of usable downstream electron reconstructed tracks, compute
    // the closest approach between their local backward-extrapolation line
    // approximations.
    //------------------------------------------------------------------------
    for (size_t i = 0; i < lines.size(); ++i) {
      for (size_t j = i + 1; j < lines.size(); ++j) {
        ++numPairsTested;

        const EntranceLine& line1 = lines.at(i);
        const EntranceLine& line2 = lines.at(j);

        // Build the Mu2e two-line closest-approach object.
        const ClosestApproachResult pca = computeClosestApproach(
            line1.point, line1.dir, line2.point, line2.dir);

        if (!pca.valid) {
          continue;
        }

        // Endpoints on the two lines at the point of closest approach.
        const XYZVectorF pcaPoint1 = pca.point1;
        const XYZVectorF pcaPoint2 = pca.point2;

        // Distance of closest approach between the two extrapolated lines.
        const double dca = pca.dca;

        // A common simple choice for a candidate vertex estimate is the
        // midpoint between the two PCA endpoints.
        const XYZVectorF vertexEstimate = 0.5f * (pcaPoint1 + pcaPoint2);

        // Signed distances tell us where the closest-approach point lies
        // relative to the tracker entrance point and the local track direction.
        //
        // Negative means "behind" the entrance point along the momentum
        // direction, which is exactly what we want when looking for an upstream
        // common origin by backward extrapolation.
        const double s1 = signedDistanceAlongLine(line1.refPos, line1.dir, pcaPoint1);
        const double s2 = signedDistanceAlongLine(line2.refPos, line2.dir, pcaPoint2);

        const bool upstreamForTrack1 = (s1 < 0.0);
        const bool upstreamForTrack2 = (s2 < 0.0);
        const bool upstreamForBoth = upstreamForTrack1 && upstreamForTrack2;
        const bool bothRank0Electrons = line1.rank0Electron && line2.rank0Electron;

        // Fill inclusive histograms first.
        hPairDCAAll->Fill(dca);
        hVertexZAll->Fill(vertexEstimate.Z());
        hVertexRAll->Fill(vertexEstimate.Rho());
        hVertexXZAll->Fill(vertexEstimate.Z(), vertexEstimate.X());

        // Fill the more restrictive upstream-only category.
        if (upstreamForBoth) {
          ++numUpstreamPairs;
          hPairDCAUpstream->Fill(dca);
          hVertexZUpstream->Fill(vertexEstimate.Z());
          hVertexRUpstream->Fill(vertexEstimate.Rho());
          hVertexXZUpstream->Fill(vertexEstimate.Z(), vertexEstimate.X());
        }

        // Fill the rank-0-electron-only histograms if BOTH tracks in the pair
        // are matched to rank-0 electrons.
        if (bothRank0Electrons) {
          ++numRank0ElectronPairs;
          hPairDCARank0Electron->Fill(dca);
          hVertexZRank0Electron->Fill(vertexEstimate.Z());
          hVertexRRank0Electron->Fill(vertexEstimate.Rho());
          hVertexXZRank0Electron->Fill(vertexEstimate.Z(), vertexEstimate.X());

          if (upstreamForBoth) {
            ++numRank0ElectronUpstreamPairs;
            hPairDCARank0ElectronUpstream->Fill(dca);
            hVertexZRank0ElectronUpstream->Fill(vertexEstimate.Z());
            hVertexRRank0ElectronUpstream->Fill(vertexEstimate.Rho());
            hVertexXZRank0ElectronUpstream->Fill(vertexEstimate.Z(), vertexEstimate.X());
          }
        }
      }
    }
  }

  //----------------------------------------------------------------------------
  // Print summary information.
  //----------------------------------------------------------------------------
  cout << "Selection: downstream e- reconstructed fit hypotheses only." << endl;
  cout << "Processed " << numTracksSeen << " selected downstream e- track hypotheses." << endl;
  cout << "Built " << numUsableLines << " usable local entrance-line approximations." << endl;
  cout << numRank0ElectronTracks << " of those usable tracks are matched to rank-0 electrons." << endl;
  cout << "Tested " << numPairsTested << " selected downstream e- track pairs." << endl;
  cout << numUpstreamPairs << " pairs had their closest-approach points upstream of TT_Front for both tracks." << endl;
  cout << numRank0ElectronPairs << " tested pairs had both tracks matched to rank-0 electrons." << endl;
  cout << numRank0ElectronUpstreamPairs << " rank-0-electron pairs also had upstream closest approach for both tracks." << endl;

  //----------------------------------------------------------------------------
  // Draw pair DCA comparison.
  //----------------------------------------------------------------------------
  TCanvas* cPairDCA = new TCanvas("cPairDCA", "Pair DCA Comparison", 900, 700);
  double maxDCA = max(hPairDCAAll->GetMaximum(), hPairDCAUpstream->GetMaximum());
  hPairDCAAll->SetMaximum(1.15 * maxDCA);
  hPairDCAAll->Draw("HIST E");
  hPairDCAUpstream->Draw("HIST E SAME");
  TLegend* legDCA = new TLegend(0.56, 0.72, 0.88, 0.88);
  legDCA->AddEntry(hPairDCAAll, "Selected downstream e- track pairs", "l");
  legDCA->AddEntry(hPairDCAUpstream, "Selected pairs with upstream closest approach", "l");
  legDCA->Draw();
  cPairDCA->SaveAs(("RecoCommonOrigin_PairDCA_" + safeSampleName + ".pdf").c_str());

  //----------------------------------------------------------------------------
  // Draw candidate vertex Z comparison.
  //----------------------------------------------------------------------------
  TCanvas* cVertexZ = new TCanvas("cVertexZ", "Candidate Vertex Z", 900, 700);
  double maxVertexZ = max(hVertexZAll->GetMaximum(), hVertexZUpstream->GetMaximum());
  hVertexZAll->SetMaximum(1.15 * maxVertexZ);
  hVertexZAll->Draw("HIST E");
  hVertexZUpstream->Draw("HIST E SAME");
  TLegend* legZ = new TLegend(0.56, 0.72, 0.88, 0.88);
  legZ->AddEntry(hVertexZAll, "Selected downstream e- track pairs", "l");
  legZ->AddEntry(hVertexZUpstream, "Selected pairs with upstream closest approach", "l");
  legZ->Draw();
  cVertexZ->SaveAs(("RecoCommonOrigin_VertexZ_" + safeSampleName + ".pdf").c_str());

  //----------------------------------------------------------------------------
  // Draw candidate vertex radius comparison.
  //----------------------------------------------------------------------------
  TCanvas* cVertexR = new TCanvas("cVertexR", "Candidate Vertex Radius", 900, 700);
  double maxVertexR = max(hVertexRAll->GetMaximum(), hVertexRUpstream->GetMaximum());
  hVertexRAll->SetMaximum(1.15 * maxVertexR);
  hVertexRAll->Draw("HIST E");
  hVertexRUpstream->Draw("HIST E SAME");
  TLegend* legR = new TLegend(0.56, 0.72, 0.88, 0.88);
  legR->AddEntry(hVertexRAll, "Selected downstream e- track pairs", "l");
  legR->AddEntry(hVertexRUpstream, "Selected pairs with upstream closest approach", "l");
  legR->Draw();
  cVertexR->SaveAs(("RecoCommonOrigin_VertexR_" + safeSampleName + ".pdf").c_str());

  //----------------------------------------------------------------------------
  // Draw 2D origin maps.
  //
  // The dashed vertical line at TT_Front z = 0 is NOT guaranteed to be the true
  // tracker-front z value in global Mu2e coordinates.  It is only a visual aid
  // in case your ntuple convention happens to center the tracker coordinates in
  // a convenient way.
  //
  // The more reliable "upstream" classification in this macro is the signed
  // distance test performed along each local line, not the absolute Z value.
  //----------------------------------------------------------------------------
  TCanvas* cMapAll = new TCanvas("cMapAll", "Candidate Common-Origin Map (Selected Downstream e- Pairs)", 900, 700);
  hVertexXZAll->Draw("COLZ");
  cMapAll->SaveAs(("RecoCommonOrigin_MapAll_" + safeSampleName + ".pdf").c_str());

  TCanvas* cMapUpstream = new TCanvas("cMapUpstream", "Candidate Common-Origin Map (Upstream Pairs)", 900, 700);
  hVertexXZUpstream->Draw("COLZ");
  cMapUpstream->SaveAs(("RecoCommonOrigin_MapUpstream_" + safeSampleName + ".pdf").c_str());

  //----------------------------------------------------------------------------
  // Draw pair DCA comparison for rank-0-electron-only pairs.
  //----------------------------------------------------------------------------
  TCanvas* cPairDCARank0 = new TCanvas("cPairDCARank0", "Pair DCA Comparison Rank-0 Electrons", 900, 700);
  double maxDCARank0 = max(hPairDCARank0Electron->GetMaximum(), hPairDCARank0ElectronUpstream->GetMaximum());
  hPairDCARank0Electron->SetMaximum(1.15 * maxDCARank0);
  hPairDCARank0Electron->Draw("HIST E");
  hPairDCARank0ElectronUpstream->Draw("HIST E SAME");
  TLegend* legDCARank0 = new TLegend(0.48, 0.72, 0.88, 0.88);
  legDCARank0->AddEntry(hPairDCARank0Electron, "Rank-0 electron track pairs", "l");
  legDCARank0->AddEntry(hPairDCARank0ElectronUpstream, "Rank-0 electron pairs with upstream closest approach", "l");
  legDCARank0->Draw();
  cPairDCARank0->SaveAs(("RecoCommonOrigin_PairDCA_Rank0Electron_" + safeSampleName + ".pdf").c_str());

  //----------------------------------------------------------------------------
  // Draw candidate vertex Z comparison for rank-0-electron-only pairs.
  //----------------------------------------------------------------------------
  TCanvas* cVertexZRank0 = new TCanvas("cVertexZRank0", "Candidate Vertex Z Rank-0 Electrons", 900, 700);
  double maxVertexZRank0 = max(hVertexZRank0Electron->GetMaximum(), hVertexZRank0ElectronUpstream->GetMaximum());
  hVertexZRank0Electron->SetMaximum(1.15 * maxVertexZRank0);
  hVertexZRank0Electron->Draw("HIST E");
  hVertexZRank0ElectronUpstream->Draw("HIST E SAME");
  TLegend* legZRank0 = new TLegend(0.48, 0.72, 0.88, 0.88);
  legZRank0->AddEntry(hVertexZRank0Electron, "Rank-0 electron track pairs", "l");
  legZRank0->AddEntry(hVertexZRank0ElectronUpstream, "Rank-0 electron pairs with upstream closest approach", "l");
  legZRank0->Draw();
  cVertexZRank0->SaveAs(("RecoCommonOrigin_VertexZ_Rank0Electron_" + safeSampleName + ".pdf").c_str());

  //----------------------------------------------------------------------------
  // Draw candidate vertex radius comparison for rank-0-electron-only pairs.
  //----------------------------------------------------------------------------
  TCanvas* cVertexRRank0 = new TCanvas("cVertexRRank0", "Candidate Vertex Radius Rank-0 Electrons", 900, 700);
  double maxVertexRRank0 = max(hVertexRRank0Electron->GetMaximum(), hVertexRRank0ElectronUpstream->GetMaximum());
  hVertexRRank0Electron->SetMaximum(1.15 * maxVertexRRank0);
  hVertexRRank0Electron->Draw("HIST E");
  hVertexRRank0ElectronUpstream->Draw("HIST E SAME");
  TLegend* legRRank0 = new TLegend(0.48, 0.72, 0.88, 0.88);
  legRRank0->AddEntry(hVertexRRank0Electron, "Rank-0 electron track pairs", "l");
  legRRank0->AddEntry(hVertexRRank0ElectronUpstream, "Rank-0 electron pairs with upstream closest approach", "l");
  legRRank0->Draw();
  cVertexRRank0->SaveAs(("RecoCommonOrigin_VertexR_Rank0Electron_" + safeSampleName + ".pdf").c_str());

  //----------------------------------------------------------------------------
  // Draw 2D origin maps for rank-0-electron-only pairs.
  //----------------------------------------------------------------------------
  TCanvas* cMapRank0 = new TCanvas("cMapRank0", "Candidate Common-Origin Map Rank-0 Electrons", 900, 700);
  hVertexXZRank0Electron->Draw("COLZ");
  cMapRank0->SaveAs(("RecoCommonOrigin_MapAll_Rank0Electron_" + safeSampleName + ".pdf").c_str());

  TCanvas* cMapRank0Upstream = new TCanvas("cMapRank0Upstream", "Candidate Common-Origin Map Rank-0 Electrons Upstream", 900, 700);
  hVertexXZRank0ElectronUpstream->Draw("COLZ");
  cMapRank0Upstream->SaveAs(("RecoCommonOrigin_MapUpstream_Rank0Electron_" + safeSampleName + ".pdf").c_str());

  //----------------------------------------------------------------------------
  // Final note printed to the terminal so it is obvious what this macro is and
  // is not doing.
  //----------------------------------------------------------------------------
  cout << endl;
  cout << "Study complete." << endl;
  cout << "Reminder: this macro uses a local straight-line approximation near the tracker entrance." << endl;
  cout << "It is a useful exploratory common-origin study, but it is not yet a full curved-track vertex fit." << endl;
}

