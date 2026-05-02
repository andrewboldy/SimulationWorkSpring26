//----------------------------------------------------------------------------------
//
// plotRecoTrackSegmentHists.C
//
// Read reconstructed track-segment information from EventNtuple `trksegs`
// and make momentum/energy histograms at multiple tracker surfaces.
//
// This follows the same RooUtil-based analysis style as
// CreatedCode/HistogramMakers/plotElectronMomHistsNEW.C, but it uses
// reconstructed track intersections instead of MC-truth genealogy.
//
// Example:
// root -l -b -q 'CreatedCode/HistogramMakers/plotRecoTrackSegmentHists.C+("CeEndpoint","myfiles.txt")'
//
//----------------------------------------------------------------------------------

//Standard Inclusions
#include <algorithm>
#include <cctype>
#include <cmath>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

//CERN ROOT Inclusions
#include "TCanvas.h"
#include "TH1F.h"
#include "TLegend.h"

//Mu2e Inclusions
#include "Offline/DataProducts/inc/PDGCode.hh"
#include "Offline/DataProducts/inc/SurfaceId.hh"

#include "EventNtuple/rooutil/inc/RooUtil.hh"
#include "EventNtuple/rooutil/inc/common_cuts.hh"

//Namespace
using namespace std;
using namespace rooutil;

namespace {

  //Helper structure that stores the information for each tracker point we want to study
  struct SurfaceConfig {
    string key;
    string title;
    int sid;
    int color;
  };

  //Helper method to determine the particle mass from the fit hypothesis
  //This is used so that the code can compute reconstructed total energy from the reconstructed momentum
  double particleMassMeV(mu2e::PDGCode::type pdg) {
    switch (pdg) {
      case mu2e::PDGCode::e_minus:
      case mu2e::PDGCode::e_plus:
        return 0.51099895;
      case mu2e::PDGCode::mu_minus:
      case mu2e::PDGCode::mu_plus:
        return 105.6583755;
      case mu2e::PDGCode::pi_minus:
      case mu2e::PDGCode::pi_plus:
        return 139.57039;
      case mu2e::PDGCode::proton:
      case mu2e::PDGCode::anti_proton:
        return 938.27208816;
      default:
        return 0.0;
    }
  }

  //Helper method to determine whether a reconstructed track segment is at the requested tracker surface
  bool hasRecoSegmentAtSurface(TrackSegment& segment, int sid) {
    return has_reco_step(segment) && segment.trkseg->sid == sid;
  }

  //Helper method to clean up strings before using them in output file names
  string sanitizeName(const string& input) {
    string output = input;
    for (char& c : output) {
      if (!(isalnum(static_cast<unsigned char>(c)) || c == '_' || c == '-')) {
        c = '_';
      }
    }
    return output;
  }

} // namespace

void plotRecoTrackSegmentHists(const string& generatorName, const string& filelist) {
  //begin the filelist reading
  ifstream listOfFiles(filelist);
  if (!listOfFiles.is_open()) {
    cerr << "Could not open file list: " << filelist << endl;
    return;
  }

  //counter for files in the list
  int fileCount = 0;
  string line;
  while (getline(listOfFiles, line)) {
    if (!line.empty()) {
      ++fileCount;
    }
  }

  //Define the tracker surfaces where we want to evaluate the reconstructed track
  //These are the tracker entrance, the middle of the tracker, and the tracker exit
  const vector<SurfaceConfig> surfaces = {
    {"TrackerEntrance", "Tracker Entrance (TT_Front)", mu2e::SurfaceIdDetail::TT_Front, kBlue + 1},
    {"TrackerMiddle",   "Tracker Middle (TT_Mid)",     mu2e::SurfaceIdDetail::TT_Mid,   kGreen + 2},
    {"TrackerExit",     "Tracker Exit (TT_Back)",      mu2e::SurfaceIdDetail::TT_Back,  kRed + 1}
  };

  //Vectors to store the momentum and energy histograms for the different tracker points
  vector<TH1F*> momentumHists;
  vector<TH1F*> energyHists;
  momentumHists.reserve(surfaces.size());
  energyHists.reserve(surfaces.size());

  //Clean up the generator name so it can be used in the output PDF names
  const string safeGeneratorName = sanitizeName(generatorName);

  //Initialize one momentum histogram and one energy histogram for each requested tracker point
  for (const auto& surface : surfaces) {
    const string momName = "hRecoMom_" + surface.key;
    const string eneName = "hRecoEnergy_" + surface.key;
    const string momTitle = generatorName + " Reco Momentum at " + surface.title + ";p [MeV/c];Entries";
    const string eneTitle = generatorName + " Reco Energy at " + surface.title + ";E [MeV];Entries";

    TH1F* hMom = new TH1F(momName.c_str(), momTitle.c_str(), 80, 90.0, 115.0);
    TH1F* hEne = new TH1F(eneName.c_str(), eneTitle.c_str(), 80, 90.0, 115.0);
    hMom->SetLineColor(surface.color);
    hMom->SetLineWidth(2);
    hEne->SetLineColor(surface.color);
    hEne->SetLineWidth(2);
    momentumHists.push_back(hMom);
    energyHists.push_back(hEne);
  }

  //open up the RooUtil and check number of events contained in the filelist
  RooUtil util(filelist);
  const int numEvents = util.GetNEvents();
  int numRecoTracks = 0;
  vector<int> segmentCounts(surfaces.size(), 0);

  cout << "Read " << fileCount << " files from " << filelist << endl;
  cout << "There are " << numEvents << " events in the file list." << endl;

  //Dig into the actual TTrees
  for (int i_event = 0; i_event < numEvents; ++i_event) {
    //Get the next event from the EventNtuple
    auto& event = util.GetEvent(i_event);

    //Get the reconstructed e_minus track hypotheses from the event
    //This mirrors the style of the previous macro, but now we use reconstructed tracks instead of MC truth
    auto eMinusTracks = event.GetTracks(is_e_minus);

    //Loop through the reconstructed e_minus tracks
    for (auto& track : eMinusTracks) {
      ++numRecoTracks;

      //Get the mass for the particle hypothesis used in the track fit
      //This is needed to compute the reconstructed total energy from the reconstructed momentum
      const double mass = particleMassMeV(static_cast<mu2e::PDGCode::type>(track.trk->pdg));

      //Loop through each tracker point we want to study
      for (size_t i_surface = 0; i_surface < surfaces.size(); ++i_surface) {
        const auto sid = surfaces[i_surface].sid;

        //Select only the reconstructed segments at this specific tracker surface
        //This is the central difference from the MC-based code: we use trksegs through track.GetSegments(...)
        auto recoSegments = track.GetSegments([sid](TrackSegment& segment) {
          return hasRecoSegmentAtSurface(segment, sid);
        });

        //Loop through the reconstructed segments found at this tracker surface
        for (auto& segment : recoSegments) {
          //Read the reconstructed momentum magnitude directly from the reconstructed segment
          const double p = segment.trkseg->mom.R();

          //Compute the reconstructed total energy from the reconstructed momentum and fit mass hypothesis
          const double e = std::sqrt(p * p + mass * mass);

          //Fill the histograms for this tracker point
          momentumHists[i_surface]->Fill(p);
          energyHists[i_surface]->Fill(e);
          ++segmentCounts[i_surface];
        }
      }
    }
  }

  //Print out results for the number of reconstructed tracks and reconstructed segments used
  cout << "Processed " << numRecoTracks << " reconstructed e- track hypotheses." << endl;
  for (size_t i_surface = 0; i_surface < surfaces.size(); ++i_surface) {
    cout << surfaces[i_surface].title << ": filled " << segmentCounts[i_surface]
         << " reconstructed segments." << endl;
  }

  //Draw and save the momentum histograms with one pad per tracker point
  TCanvas* cMom = new TCanvas("cRecoMom", "Reco Track Segment Momentum", 1500, 450);
  cMom->Divide(static_cast<int>(surfaces.size()), 1);
  for (size_t i_surface = 0; i_surface < surfaces.size(); ++i_surface) {
    cMom->cd(static_cast<int>(i_surface) + 1);
    momentumHists[i_surface]->Draw("HIST E");
  }
  cMom->SaveAs(("RecoTrackSegmentMomentum_" + safeGeneratorName + ".pdf").c_str());

  //Draw and save the energy histograms with one pad per tracker point
  TCanvas* cEnergy = new TCanvas("cRecoEnergy", "Reco Track Segment Energy", 1500, 450);
  cEnergy->Divide(static_cast<int>(surfaces.size()), 1);
  for (size_t i_surface = 0; i_surface < surfaces.size(); ++i_surface) {
    cEnergy->cd(static_cast<int>(i_surface) + 1);
    energyHists[i_surface]->Draw("HIST E");
  }
  cEnergy->SaveAs(("RecoTrackSegmentEnergy_" + safeGeneratorName + ".pdf").c_str());

  //Create an overlay plot for the reconstructed momentum at the different tracker points
  TCanvas* cOverlayMom = new TCanvas("cOverlayMom", "Reco Momentum Comparison", 900, 700);
  double maxMom = 0.0;
  for (auto* hist : momentumHists) {
    maxMom = max(maxMom, hist->GetMaximum());
  }
  for (size_t i_surface = 0; i_surface < momentumHists.size(); ++i_surface) {
    momentumHists[i_surface]->SetMaximum(1.15 * maxMom);
    momentumHists[i_surface]->Draw(i_surface == 0 ? "HIST E" : "HIST E SAME");
  }
  TLegend* momLeg = new TLegend(0.58, 0.68, 0.88, 0.88);
  for (size_t i_surface = 0; i_surface < surfaces.size(); ++i_surface) {
    momLeg->AddEntry(momentumHists[i_surface], surfaces[i_surface].title.c_str(), "l");
  }
  momLeg->Draw();
  cOverlayMom->SaveAs(("RecoTrackSegmentMomentumOverlay_" + safeGeneratorName + ".pdf").c_str());

  //Create an overlay plot for the reconstructed energy at the different tracker points
  TCanvas* cOverlayEnergy = new TCanvas("cOverlayEnergy", "Reco Energy Comparison", 900, 700);
  double maxEnergy = 0.0;
  for (auto* hist : energyHists) {
    maxEnergy = max(maxEnergy, hist->GetMaximum());
  }
  for (size_t i_surface = 0; i_surface < energyHists.size(); ++i_surface) {
    energyHists[i_surface]->SetMaximum(1.15 * maxEnergy);
    energyHists[i_surface]->Draw(i_surface == 0 ? "HIST E" : "HIST E SAME");
  }
  TLegend* energyLeg = new TLegend(0.58, 0.68, 0.88, 0.88);
  for (size_t i_surface = 0; i_surface < surfaces.size(); ++i_surface) {
    energyLeg->AddEntry(energyHists[i_surface], surfaces[i_surface].title.c_str(), "l");
  }
  energyLeg->Draw();
  cOverlayEnergy->SaveAs(("RecoTrackSegmentEnergyOverlay_" + safeGeneratorName + ".pdf").c_str());
}

