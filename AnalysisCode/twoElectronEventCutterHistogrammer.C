//----------------------------------------------------------------------------------
//
// twoElectronEventCutterHistogrammer.C
// Written by Andrew Boldy University of South Carolina, 2026
// Assisted by Codex
//
// Purpose: Print, store, and write to txt the number of all events that have two electrons classified as rank 0 according to the trkmcsim branch (MC truth)
// Eventually want to create a histogram that takes these specific events and creates a histogram of the events trkmcsim and trksegs specficially.
//
//----------------------------------------------------------------------------------

//Standard Inclusions
#include <algorithm>
#include <cmath>
#include <string>
#include <fstream>
#include <iostream>
#include <limits>
#include <vector>

//CERN ROOT Inclusions
#include <TCanvas.h>
#include <TDirectory.h>
#include <TFile.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TLegend.h>
#include <TROOT.h>
#include <TSystem.h>
#include <TStopwatch.h>

//Mu2e Inclusions
#include "EventNtuple/rooutil/inc/RooUtil.hh"
#include "EventNtuple/rooutil/inc/common_cuts.hh"

//Namespace
using namespace std;
using namespace rooutil;

void twoElectronEventCutterHistogrammer(const string& generatorName, const string& fileName)
{
  ifstream file(fileName);
  if (!file.is_open())
  {
    cerr << "ERROR: could not open input file or filelist: " << fileName << endl;
    return;
  }
  file.close();

  //open up the RooUtil and check number of events contained in the filelist
  RooUtil util(fileName);
  int numEvents = util.GetNEvents();
  cout << "There are " << numEvents << " in the filelist." << endl;

  const string outputFileName = "twoElectronEvents_" + generatorName + ".txt";
  ofstream outputFile(outputFileName);
  if (!outputFile.is_open())
  {
    cerr << "ERROR: could not create output text file: " << outputFileName << endl;
    return;
  }
  outputFile << "entry run subrun event num_rank0_electrons electron1_angle_deg electron2_angle_deg" << endl;

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
  const double originXYMin = -1000.0;
  const double originXYMax = 1000.0;
  const double originZMin = -10000.0;
  const double originZMax = 16000.0;

  TH1F* hMCTAllElectronMomentum = new TH1F(
    "hMCTAllElectronMomentum",
    "Monte Carlo truth: all trkmcsim electrons;trkmcsim momentum [MeV/c];Electrons",
    200, 0.0, 200.0);
  TH1F* hMCTRank0ElectronMomentum = new TH1F(
    "hMCTRank0ElectronMomentum",
    "Monte Carlo truth: rank-0 trkmcsim electrons;trkmcsim momentum [MeV/c];Electrons",
    200, 0.0, 200.0);

  TH1F* hMCTAllElectronOriginX = new TH1F(
    "hMCTAllElectronOriginX",
    "Monte Carlo truth: all trkmcsim electron origins X;origin x [mm];Electrons",
    200, originXYMin, originXYMax);
  TH1F* hMCTAllElectronOriginY = new TH1F(
    "hMCTAllElectronOriginY",
    "Monte Carlo truth: all trkmcsim electron origins Y;origin y [mm];Electrons",
    200, originXYMin, originXYMax);
  TH1F* hMCTAllElectronOriginZ = new TH1F(
    "hMCTAllElectronOriginZ",
    "Monte Carlo truth: all trkmcsim electron origins Z;origin z [mm];Electrons",
    260, originZMin, originZMax);

  TH1F* hMCTRank0ElectronOriginX = new TH1F(
    "hMCTRank0ElectronOriginX",
    "Monte Carlo truth: rank-0 trkmcsim electron origins X;origin x [mm];Electrons",
    200, originXYMin, originXYMax);
  TH1F* hMCTRank0ElectronOriginY = new TH1F(
    "hMCTRank0ElectronOriginY",
    "Monte Carlo truth: rank-0 trkmcsim electron origins Y;origin y [mm];Electrons",
    200, originXYMin, originXYMax);
  TH1F* hMCTRank0ElectronOriginZ = new TH1F(
    "hMCTRank0ElectronOriginZ",
    "Monte Carlo truth: rank-0 trkmcsim electron origins Z;origin z [mm];Electrons",
    260, originZMin, originZMax);

  TH2F* hMCTAllElectronOriginXY = new TH2F(
    "hMCTAllElectronOriginXY",
    "Monte Carlo truth: all trkmcsim electron origins XY;origin x [mm];origin y [mm]",
    200, originXYMin, originXYMax, 200, originXYMin, originXYMax);
  TH2F* hMCTAllElectronOriginXZ = new TH2F(
    "hMCTAllElectronOriginXZ",
    "Monte Carlo truth: all trkmcsim electron origins XZ;origin z [mm];origin x [mm]",
    260, originZMin, originZMax, 200, originXYMin, originXYMax);
  TH2F* hMCTAllElectronOriginYZ = new TH2F(
    "hMCTAllElectronOriginYZ",
    "Monte Carlo truth: all trkmcsim electron origins YZ;origin z [mm];origin y [mm]",
    260, originZMin, originZMax, 200, originXYMin, originXYMax);

  TH2F* hMCTRank0ElectronOriginXY = new TH2F(
    "hMCTRank0ElectronOriginXY",
    "Monte Carlo truth: rank-0 trkmcsim electron origins XY;origin x [mm];origin y [mm]",
    200, originXYMin, originXYMax, 200, originXYMin, originXYMax);
  TH2F* hMCTRank0ElectronOriginXZ = new TH2F(
    "hMCTRank0ElectronOriginXZ",
    "Monte Carlo truth: rank-0 trkmcsim electron origins XZ;origin z [mm];origin x [mm]",
    260, originZMin, originZMax, 200, originXYMin, originXYMax);
  TH2F* hMCTRank0ElectronOriginYZ = new TH2F(
    "hMCTRank0ElectronOriginYZ",
    "Monte Carlo truth: rank-0 trkmcsim electron origins YZ;origin z [mm];origin y [mm]",
    260, originZMin, originZMax, 200, originXYMin, originXYMax);

  vector<int> dualElectronEntries;
  long long selectedAllElectronFills = 0;
  long long selectedRank0ElectronFills = 0;
  const double pi = acos(-1.0);
  double minElectronAngleDeg = numeric_limits<double>::max();
  double maxElectronAngleDeg = -numeric_limits<double>::max();
  bool haveElectronAngleStats = false;

  auto getAngleFromZDeg = [pi](const auto& momentum) {
    const double momentumMagnitude = momentum.R();
    if (momentumMagnitude <= 0.0)
    {
      return -1.0;
    }

    double cosTheta = momentum.z() / momentumMagnitude;
    cosTheta = max(-1.0, min(1.0, cosTheta));
    return acos(cosTheta) * 180.0 / pi;
  };

  //Dig into the actual TTrees
  for (int i_event = 0; i_event < numEvents; i_event++)
  {
    int numEventElectrons = 0;
    vector<double> eventElectronAnglesDeg;
    vector<const mu2e::SimInfo*> eventElectrons;
    vector<const mu2e::SimInfo*> eventRank0Electrons;
    auto& event = util.GetEvent(i_event);
    const auto& e_minus_tracks = event.GetTracks(is_e_minus);

    for (auto& track : e_minus_tracks)
    {
      if (track.trkmcsim == nullptr)
      {
        continue;
      }

      for (const auto& mctrack : *(track.trkmcsim))
      {
        if (!(mctrack.valid && mctrack.pdg == 11))
        {
          continue;
        }

        eventElectrons.push_back(&mctrack);
        if (mctrack.rank == 0)
        {
          numEventElectrons++;
          eventRank0Electrons.push_back(&mctrack);
          eventElectronAnglesDeg.push_back(getAngleFromZDeg(mctrack.mom));
        } //end work on the electron pdg tracks
      }//end the interface with the TTree at trkmcsim branch level
    }//end the interface with the reconstructed e-minus track level

    if (numEventElectrons != 2)
    {
      numEventElectrons=0;
      continue;
    }//end the check to continue if the number of rank 0 electrons is not exactly 2

    else
    {
      int run = -1;
      int subrun = -1;
      int eventNumber = -1;
      cout << "Entry " << i_event;
      if (event.evtinfo != nullptr)
      {
        run = event.evtinfo->run;
        subrun = event.evtinfo->subrun;
        eventNumber = event.evtinfo->event;
        cout << " Run " << event.evtinfo->run
             << " Subrun " << event.evtinfo->subrun
             << " Event " << event.evtinfo->event;
      }
      cout << " has exactly " << numEventElectrons << " electrons in MonteCarlo truth able to be reconstructed." << endl;
      cout << "  Electron 1 angle from z-axis: " << eventElectronAnglesDeg.at(0) << " degrees" << endl;
      cout << "  Electron 2 angle from z-axis: " << eventElectronAnglesDeg.at(1) << " degrees" << endl;
      outputFile << i_event << " "
                 << run << " "
                 << subrun << " "
                 << eventNumber << " "
                 << numEventElectrons << " "
                 << eventElectronAnglesDeg.at(0) << " "
                 << eventElectronAnglesDeg.at(1) << endl;
      for (const double electronAngleDeg : eventElectronAnglesDeg)
      {
        if (electronAngleDeg >= 0.0)
        {
          haveElectronAngleStats = true;
          minElectronAngleDeg = min(minElectronAngleDeg, electronAngleDeg);
          maxElectronAngleDeg = max(maxElectronAngleDeg, electronAngleDeg);
        }
      }
      dualElectronEntries.push_back(i_event);
      selectedAllElectronFills += eventElectrons.size();
      selectedRank0ElectronFills += eventRank0Electrons.size();

      for (const auto* electron : eventElectrons)
      {
        hMCTAllElectronMomentum->Fill(electron->mom.R());
        hMCTAllElectronOriginX->Fill(electron->pos.x());
        hMCTAllElectronOriginY->Fill(electron->pos.y());
        hMCTAllElectronOriginZ->Fill(electron->pos.z());
        hMCTAllElectronOriginXY->Fill(electron->pos.x(), electron->pos.y());
        hMCTAllElectronOriginXZ->Fill(electron->pos.z(), electron->pos.x());
        hMCTAllElectronOriginYZ->Fill(electron->pos.z(), electron->pos.y());
      }

      for (const auto* electron : eventRank0Electrons)
      {
        hMCTRank0ElectronMomentum->Fill(electron->mom.R());
        hMCTRank0ElectronOriginX->Fill(electron->pos.x());
        hMCTRank0ElectronOriginY->Fill(electron->pos.y());
        hMCTRank0ElectronOriginZ->Fill(electron->pos.z());
        hMCTRank0ElectronOriginXY->Fill(electron->pos.x(), electron->pos.y());
        hMCTRank0ElectronOriginXZ->Fill(electron->pos.z(), electron->pos.x());
        hMCTRank0ElectronOriginYZ->Fill(electron->pos.z(), electron->pos.y());
      }

      numEventElectrons=0;
    }
  } //end the interface with the TTree at event level

  outputFile.close();
  cout << "Wrote " << dualElectronEntries.size() << " selected events to " << outputFileName << endl;
  if (haveElectronAngleStats)
  {
    cout << "Minimum electron angle from z-axis: " << minElectronAngleDeg << " degrees" << endl;
    cout << "Maximum electron angle from z-axis: " << maxElectronAngleDeg << " degrees" << endl;
  }
  else
  {
    cout << "No selected events, so no minimum or maximum electron angle from z-axis was computed." << endl;
  }
  cout << "Monte Carlo truth electron fills from selected events: all valid electrons = "
       << selectedAllElectronFills << ", rank-0 electrons = " << selectedRank0ElectronFills << endl;

  const string histogramFileName = "twoElectronEvents_" + generatorName + "_histograms.root";
  TFile histogramFile(histogramFileName.c_str(), "RECREATE");
  auto* mcTruthDirectory = histogramFile.mkdir("Monte Carlo truth");
  if (mcTruthDirectory != nullptr)
  {
    mcTruthDirectory->cd();
  }

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

  const string plotsDirectory = "Plots";
  gSystem->mkdir(plotsDirectory.c_str(), true);

  TStopwatch plotTimer;
  plotTimer.Start();

  TCanvas* cMCTMomentum = new TCanvas("cMCTMomentum", "Monte Carlo truth: trkmcsim momentum", 900, 700);
  hMCTAllElectronMomentum->SetLineColor(kBlack);
  hMCTAllElectronMomentum->SetLineWidth(2);
  hMCTAllElectronMomentum->SetFillColor(kGray);
  hMCTAllElectronMomentum->SetFillStyle(3004);
  hMCTRank0ElectronMomentum->SetLineColor(kRed);
  hMCTRank0ElectronMomentum->SetLineWidth(2);
  const double maxMomentumBinContent = max(hMCTAllElectronMomentum->GetMaximum(),
                                           hMCTRank0ElectronMomentum->GetMaximum());
  if (maxMomentumBinContent > 0.0)
  {
    hMCTAllElectronMomentum->SetMaximum(1.15 * maxMomentumBinContent);
  }
  hMCTAllElectronMomentum->Draw("HIST");
  hMCTRank0ElectronMomentum->Draw("HIST E SAME");
  TLegend* momentumLegend = new TLegend(0.60, 0.72, 0.88, 0.88);
  momentumLegend->AddEntry(hMCTAllElectronMomentum, "All electrons", "l");
  momentumLegend->AddEntry(hMCTRank0ElectronMomentum, "Rank 0 electrons", "l");
  momentumLegend->Draw();
  cMCTMomentum->Write();
  const string momentumPdfName = plotsDirectory + "/twoElectronEvents_" + generatorName + "_MonteCarloTruthMomentum.pdf";
  cMCTMomentum->SaveAs(momentumPdfName.c_str());

  TCanvas* cMCTOrigins = new TCanvas("cMCTOrigins", "Monte Carlo truth: trkmcsim origins", 1200, 800);
  cMCTOrigins->Divide(3, 2);
  cMCTOrigins->cd(1);
  hMCTAllElectronOriginX->Draw("HIST E");
  cMCTOrigins->cd(2);
  hMCTAllElectronOriginY->Draw("HIST E");
  cMCTOrigins->cd(3);
  hMCTAllElectronOriginZ->Draw("HIST E");
  cMCTOrigins->cd(4);
  hMCTRank0ElectronOriginX->Draw("HIST E");
  cMCTOrigins->cd(5);
  hMCTRank0ElectronOriginY->Draw("HIST E");
  cMCTOrigins->cd(6);
  hMCTRank0ElectronOriginZ->Draw("HIST E");
  cMCTOrigins->Write();
  const string originsPdfName = plotsDirectory + "/twoElectronEvents_" + generatorName + "_MonteCarloTruthOrigins.pdf";
  cMCTOrigins->SaveAs(originsPdfName.c_str());

  TCanvas* cMCTOriginMaps = new TCanvas("cMCTOriginMaps", "Monte Carlo truth: trkmcsim origin maps", 1400, 850);
  cMCTOriginMaps->Divide(3, 2);
  cMCTOriginMaps->cd(1);
  hMCTAllElectronOriginXY->Draw("COLZ");
  cMCTOriginMaps->cd(2);
  hMCTAllElectronOriginXZ->Draw("COLZ");
  cMCTOriginMaps->cd(3);
  hMCTAllElectronOriginYZ->Draw("COLZ");
  cMCTOriginMaps->cd(4);
  hMCTRank0ElectronOriginXY->Draw("COLZ");
  cMCTOriginMaps->cd(5);
  hMCTRank0ElectronOriginXZ->Draw("COLZ");
  cMCTOriginMaps->cd(6);
  hMCTRank0ElectronOriginYZ->Draw("COLZ");
  cMCTOriginMaps->Write();
  const string originMapsPdfName = plotsDirectory + "/twoElectronEvents_" + generatorName + "_MonteCarloTruthOriginMaps.pdf";
  cMCTOriginMaps->SaveAs(originMapsPdfName.c_str());
  plotTimer.Stop();

  histogramFile.Close();

  delete cMCTMomentum;
  delete cMCTOrigins;
  delete cMCTOriginMaps;

  delete hMCTAllElectronMomentum;
  delete hMCTRank0ElectronMomentum;
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

  TH1::AddDirectory(oldAddDirectoryStatus);
  gROOT->SetBatch(wasBatchMode);

  cout << "Writing and saving Monte Carlo truth plots took "
       << plotTimer.RealTime() << " seconds of wall time." << endl;
  cout << "Wrote Monte Carlo truth histograms to " << histogramFileName << endl;
  cout << "Wrote Monte Carlo truth PDF plots to " << momentumPdfName
       << ", " << originsPdfName
       << ", and " << originMapsPdfName << endl;

}//end twoElectronEventCutterHistogrammer method

