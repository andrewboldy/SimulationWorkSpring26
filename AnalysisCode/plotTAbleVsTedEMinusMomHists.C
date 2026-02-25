//----------------------------------------------------------------------------------

//plotElectronMomHist(string filelist)
//Written by Andrew Boldy
//University of South Carolina
//Fall 2025

//----------------------------------------------------------------------------------

//Plots the electron momentum for a given non-mixed filelist in a histogram.
//Must be run using C++ so compile mode in order to rectify issues that are generated when using common_cuts.hh

//----------------------------------------------------------------------------------

//My Inclusions

//Standard Inclusions
#include <string>
#include <fstream>
#include <iostream>
#include <vector>

//CERN ROOT Inclusions
#include <TCanvas.h>
#include <TH1F.h>
#include <TLegend.h>

//Mu2e Inclusions
#include "EventNtuple/rooutil/inc/RooUtil.hh"
#include "EventNtuple/rooutil/inc/common_cuts.hh"

//Personal Inclusions (if any)

//Namespace
using namespace std;

void plotElectronMomHist(const string &triggeredFilelist, const string &triggerableFilelist)
{
  //Initilizing the histograms and the TCanvas
  TCanvas* c1 = new TCanvas("c1","c1");
  TH1F* eMinusAllTriggeredStartMomHist = new TH1F("eMinusAllTriggeredStartMomHist","eMinus Start Momentum (Triggered)", 50, 10, 130);
  eMinusAllTriggeredStartMomHist->GetXaxis()->SetTitle("Start Momentum (MeV/c)");
  eMinusAllTriggeredStartMomHist->GetYaxis()->SetTitle("Counts");
  eMinusAllTriggeredStartMomHist->SetStats(0);
  TH1F* eMinusAllTriggerableStartMomHist = new TH1F("eMinusAllTriggerableStartMomHist", "eMinus Start Momentum (Triggerable)", 50, 10, 130);
  eMinusAllTriggerableStartMomHist->GetXaxis()->SetTitle("Start Momentum (MeV/c)");
  eMinusAllTriggerableStartMomHist->GetYaxis()->SetTitle("Counts");
  eMinusAllTriggerableStartMomHist->SetStats(0);
  TH1F* eMinusRank0TriggeredStartMomHist = new TH1F("eMinusRank0TriggeredStartMomHist", "Rank 0 eMinus Start Momentum (Triggered)", 50, 10, 130);
  eMinusRank0TriggeredStartMomHist->GetXaxis()->SetTitle("Start Momentum (MeV/c)");
  eMinusRank0TriggeredStartMomHist->GetYaxis()->SetTitle("Counts");
  eMinusRank0TriggeredStartMomHist->SetStats(0);
  TH1F* eMinusRank0TriggerableStartMomHist = new TH1F("eMinusRank0TriggerableStartMomHist", "Rank 0 eMinus Start Momentum (Triggerable)", 50, 10, 130);
  eMinusRank0TriggerableStartMomHist->GetXaxis()->SetTitle("Start Momentum (MeV/c)");
  eMinusRank0TriggerableStartMomHist->GetYaxis()->SetTitle("Counts");
  eMinusRank0TriggerableStartMomHist->SetStats(0);
  
  TH1F* eMinusAllTriggeredEndMomHist = new TH1F("eMinusAllTriggeredEndMomHist","eMinus End Momentum (Triggered)", 10, -10, 10);
  eMinusAllTriggeredEndMomHist->GetXaxis()->SetTitle("End Momentum (MeV/c)");
  eMinusAllTriggeredEndMomHist->GetYaxis()->SetTitle("Counts");
  eMinusAllTriggeredEndMomHist->SetStats(0);
  TH1F* eMinusAllTriggerableEndMomHist = new TH1F("eMinusAllTriggerableEndMomHist", "eMinus End Momentum (Triggerable)", 10, -10, 10);
  eMinusAllTriggerableEndMomHist->GetXaxis()->SetTitle("End Momentum (MeV/c)");
  eMinusAllTriggerableEndMomHist->GetYaxis()->SetTitle("Counts");
  eMinusAllTriggerableEndMomHist->SetStats(0);
  TH1F* eMinusRank0TriggeredEndMomHist = new TH1F("eMinusRank0TriggeredEndMomHist", "Rank 0 eMinus End Momentum (Triggered)", 10, -10, 10);
  eMinusRank0TriggeredEndMomHist->GetXaxis()->SetTitle("End Momentum (MeV/c)");
  eMinusRank0TriggeredEndMomHist->GetYaxis()->SetTitle("Counts");
  eMinusRank0TriggeredEndMomHist->SetStats(0);
  TH1F* eMinusRank0TriggerableEndMomHist = new TH1F("eMinusRank0TriggerableEndMomHist", "Rank 0 eMinus End Momentum (Triggerable)", 10, -10, 10);
  eMinusRank0TriggerableEndMomHist->GetXaxis()->SetTitle("End Momentum (MeV/c)");
  eMinusRank0TriggerableEndMomHist->GetYaxis()->SetTitle("Counts");
  eMinusRank0TriggerableEndMomHist->SetStats(0);

  //Work on the triggered file
  cout << "Beginning the analysis of the triggered file." << endl;

  RooUtil triggeredUtil(triggeredFilelist);
  int numTriggeredEvents = triggeredUtil.GetNEvents();
  cout << "This triggered filelist has " << numTriggeredEvents << " events." << endl;

  ifstream listOfFilesA(triggeredFilelist);
  int triggeredFileCount = 0;
  string lineA;
  while (getline(listOfFilesA,lineA))
  {
    if (!lineA.empty())
    {
      triggeredFileCount++;
    }
  }
  cout << "The triggered filelist has: " << triggeredFileCount << " files stored in it." << endl;

  for (int i_event = 0; i_event < numTriggeredEvents; i_event++)
  {
    auto& event = triggeredUtil.GetEvent(i_event);
    auto e_minus_tracks = event.GetTracks(is_e_minus);
    for (auto& track : e_minus_tracks)
    {
        if (track.trkmcsim != nullptr)
        {
          for (const auto& mctrack : *(track.trkmcsim))
          {
            if (mctrack.pdg == 11)
            {
              eMinusAllTriggeredStartMomHist->Fill(mctrack.mom.R());
              eMinusAllTriggeredEndMomHist->Fill(mctrack.endmom.R());
              if (mctrack.rank == 0)
              {
                eMinusRank0TriggeredStartMomHist->Fill(mctrack.mom.R());
                eMinusRank0TriggeredEndMomHist->Fill(mctrack.endmom.R());
              }
            }
          }
        }
     }
  }

  cout << "Beginning the analysis of the triggerable files." << endl;
  ifstream listOfFilesB(triggerableFilelist);
  int triggerableFileCount = 0;
  string lineB;
  while (getline(listOfFilesB,lineB))
  {
    if (!lineB.empty())
    {
      triggerableFileCount++;
    }
  }
  cout << "The triggerable filelist has: " << triggerableFileCount << " files stored in it." << endl;

  RooUtil triggerableUtil(triggerableFilelist);
  int numTriggerableEvents = triggerableUtil.GetNEvents();
  cout << "This triggerable filelist has " << numTriggerableEvents << " events." << endl;
  
  for (int j_event = 0; j_event < numTriggerableEvents; j_event++)
  {
    auto& event = triggerableUtil.GetEvent(j_event);
    auto e_minus_tracks = event.GetTracks(is_e_minus);
    for (auto& track : e_minus_tracks)
    {
        if (track.trkmcsim != nullptr)
        {
          for (const auto& mctrack : *(track.trkmcsim))
          {
            if (mctrack.pdg == 11)
            {
              eMinusAllTriggerableStartMomHist->Fill(mctrack.mom.R());
              eMinusAllTriggerableEndMomHist->Fill(mctrack.endmom.R());
              if (mctrack.rank == 0)
              {
                eMinusRank0TriggerableStartMomHist->Fill(mctrack.mom.R());
                eMinusRank0TriggerableEndMomHist->Fill(mctrack.endmom.R());
              }
            }
          }
        }
     }
  }
  // Print which filelists were passed in
  cout << "Triggered filelist path: " << triggeredFilelist << endl;
  cout << "Triggerable filelist path: " << triggerableFilelist << endl;

  // Print counts returned by RooUtil
  cout << "RooUtil reported triggered events: " << numTriggeredEvents << endl;
  cout << "RooUtil reported triggerable events: " << numTriggerableEvents << endl;

  // Print file counts you computed
  cout << "triggeredFileCount = " << triggeredFileCount << ", triggerableFileCount = " << triggerableFileCount << endl;

  // Print histogram bookkeeping
  cout << "Entries: eMinusAllTriggeredStartMomHist = " << eMinusAllTriggeredStartMomHist->GetEntries()
       << ", eMinusAllTriggerableStartMomHist = " << eMinusAllTriggerableStartMomHist->GetEntries() << endl;

  cout << "Integral: eMinusAllTriggeredStartMomHist = " << eMinusAllTriggeredStartMomHist->Integral()
       << ", eMinusAllTriggerableStartMomHist = " << eMinusAllTriggerableStartMomHist->Integral() << endl;


cout << "Creating and saving 'all' histograms." << endl;
eMinusAllTriggeredStartMomHist->Draw();
c1->SaveAs("eMinusAllTriggeredStartMomHist.pdf");
c1->Clear();
eMinusAllTriggerableStartMomHist->Draw();
c1->SaveAs("eMinusAllTriggerableStartMomHist.pdf");
c1->Clear();
//eMinusAllTriggeredEndMomHist->Draw();
//c1->SaveAs("eMinusAllTriggeredEndMomHist.pdf");
//c1->Clear();
//eMinusAllTriggerableEndMomHist->Draw();
//c1->SaveAs("eMinusAllTriggerableEndMomHist.pdf");
//c1->Clear();

cout << "Creating and saving 'rank 0' histograms." << endl;
eMinusRank0TriggeredStartMomHist->Draw();
c1->SaveAs("eMinusRank0TriggeredStartMomHist.pdf");
c1->Clear();
eMinusRank0TriggerableStartMomHist->Draw();
c1->SaveAs("eMinusRank0TriggerableStartMomHist.pdf");
c1->Clear();
//eMinusRank0TriggeredEndMomHist->Draw();
//c1->SaveAs("eMinusRank0TriggeredEndMomHist.pdf");
//c1->Clear();
//eMinusRank0TriggerableEndMomHist->Draw();
//c1->SaveAs("eMinusRank0TriggerableEndMomHist.pdf");
//c1->Clear();

cout << "Creating and saving combined triggered and triggerable histograms for 'all' eMinus." << endl;
eMinusAllTriggerableStartMomHist->SetTitle("All eMinus Start Momentum Comparison (MeV/c)");
eMinusAllTriggeredStartMomHist->SetLineColor(kBlue);
eMinusAllTriggerableStartMomHist->SetLineColor(kRed);
eMinusAllTriggerableStartMomHist->Draw("HIST");
eMinusAllTriggeredStartMomHist->Draw("HIST SAME");
TLegend* legendA = new TLegend(0.45, 0.7, 0.73, 0.88); // moved left
legendA->SetTextSize(0.03);  // smaller text
legendA->SetBorderSize(0);   // optional: remove border
legendA->SetFillStyle(0);    // optional: transparent background
legendA->AddEntry(eMinusAllTriggeredStartMomHist, "Triggered", "l");
legendA->AddEntry(eMinusAllTriggerableStartMomHist, "Triggerable", "l");
legendA->Draw();
c1->SaveAs("eMinusAllTriggeredAndTriggerableStartMomHist.pdf");
c1->Clear();
//eMinusAllTriggerableEndMomHist->SetTitle("All eMinus End Momentum Comparison (MeV/c)");
//eMinusAllTriggeredEndMomHist->SetLineColor(kBlue);
//eMinusAllTriggerableEndMomHist->SetLineColor(kRed);
//eMinusAllTriggerableEndMomHist->Draw("HIST");
//eMinusAllTriggeredEndMomHist->Draw("HIST SAME");
//TLegend* legendB = new TLegend(0.6, 0.7, 0.88, 0.88); // shifted left
//legendB->SetTextSize(0.03); // smaller text
//legendB->SetBorderSize(0);  // optional: remove border
//legendB->SetFillStyle(0);   // optional: transparent background
//legendB->AddEntry(eMinusAllTriggeredEndMomHist, "Triggered", "l");
//legendB->AddEntry(eMinusAllTriggerableEndMomHist, "Triggerable", "l");
//legendB->Draw();
//c1->SaveAs("eMinusAllTriggeredAndTriggerableEndMomHist.pdf");
//c1->Clear();

cout << "Creating and saving combined triggered and triggerable histograms for 'rank 0' eMinus." << endl;
eMinusRank0TriggerableStartMomHist->SetTitle("Rank 0 eMinus Start Momentum Comparison (MeV/c)");
eMinusRank0TriggeredStartMomHist->SetLineColor(kBlue);
eMinusRank0TriggerableStartMomHist->SetLineColor(kRed);
eMinusRank0TriggerableStartMomHist->Draw("HIST");
eMinusRank0TriggeredStartMomHist->Draw("HIST SAME");
TLegend* legendC = new TLegend(0.45, 0.7, 0.73, 0.88); // shifted left
legendC->SetTextSize(0.03);  // smaller, balanced text
legendC->SetBorderSize(0);   // optional: removes border
legendC->SetFillStyle(0);    // optional: transparent background
legendC->AddEntry(eMinusRank0TriggeredStartMomHist, "Triggered", "l");
legendC->AddEntry(eMinusRank0TriggerableStartMomHist, "Triggerable", "l");
legendC->Draw();
c1->SaveAs("eMinusRank0TriggeredAndTriggerableStartMomHist.pdf");
c1->Clear();
//eMinusRank0TriggerableEndMomHist->SetTitle("Rank 0 eMinus End Momentum Comparison (MeV/c)");
//eMinusRank0TriggeredEndMomHist->SetLineColor(kBlue);
//eMinusRank0TriggerableEndMomHist->SetLineColor(kRed);
//eMinusRank0TriggerableEndMomHist->Draw("HIST");
//eMinusRank0TriggeredEndMomHist->Draw("HIST SAME");
//TLegend* legendD = new TLegend(0.6, 0.7, 0.88, 0.88); // moved left
//legendD->SetTextSize(0.03);  // smaller text
//legendD->SetBorderSize(0);   // remove border (optional)
//legendD->SetFillStyle(0);    // transparent background (optional)
//legendD->AddEntry(eMinusRank0TriggeredEndMomHist, "Triggered", "l");
//legendD->AddEntry(eMinusRank0TriggerableEndMomHist, "Triggerable", "l");
//legendD->Draw();
//c1->SaveAs("eMinusRank0TriggeredAndTriggerableEndMomHist.pdf");
//c1->Clear();

cout << "Creating and saving the combined triggered histograms for 'all' and 'rank 0' histograms for eMinus." << endl;
eMinusAllTriggerableStartMomHist->SetTitle("Triggerable 'all' and 'rank 0' eMinus Start Momentum Comparison (MeV/c)");
eMinusAllTriggerableStartMomHist->SetLineColor(kBlue);
eMinusRank0TriggerableStartMomHist->SetLineColor(kRed);
eMinusAllTriggerableStartMomHist->Draw("HIST");
eMinusRank0TriggerableStartMomHist->Draw("HIST SAME");
TLegend* legendE = new TLegend(0.45, 0.7, 0.73, 0.88); // shifted left
legendE->SetTextSize(0.03);  // smaller, balanced text
legendE->SetBorderSize(0);   // optional: removes border
legendE->SetFillStyle(0);    // optional: transparent background
legendE->AddEntry(eMinusAllTriggerableStartMomHist, "All", "l");
legendE->AddEntry(eMinusRank0TriggerableStartMomHist, "Rank 0", "l");
legendE->Draw();
//gPad->SetLogy();
c1->SaveAs("eMinusTriggerableAllAndRank0StartMomHist.pdf");
c1->Clear();
eMinusAllTriggeredStartMomHist->SetTitle("Triggered 'all' and 'rank 0' eMinus Start Momentum Comparison (MeV/c)");
eMinusAllTriggeredStartMomHist->SetLineColor(kBlue);
eMinusRank0TriggeredStartMomHist->SetLineColor(kRed);
eMinusAllTriggeredStartMomHist->Draw("HIST");
eMinusRank0TriggeredStartMomHist->Draw("HIST SAME");
TLegend* legendF = new TLegend(0.45, 0.7, 0.73, 0.88); // shifted left
legendF->SetTextSize(0.03);  // smaller, balanced text
legendF->SetBorderSize(0);   // optional: removes border
legendF->SetFillStyle(0);    // optional: transparent background
legendF->AddEntry(eMinusAllTriggeredStartMomHist, "All", "l");
legendF->AddEntry(eMinusRank0TriggeredStartMomHist, "Rank 0", "l");
legendF->Draw();
//gPad->SetLogy();
c1->SaveAs("eMinusTriggeredAllAndRank0StartMomHist.pdf");
c1->Clear();
//eMinusAllTriggerableEndMomHist->SetTitle("Triggerable 'all' and 'rank 0' eMinus End Momentum Comparison (MeV/c)");
//eMinusAllTriggerableEndMomHist->SetLineColor(kBlue);
//eMinusRank0TriggerableEndMomHist->SetLineColor(kRed);
//eMinusAllTriggerableEndMomHist->Draw("HIST");
//eMinusRank0TriggerableEndMomHist->Draw("HIST SAME");
//TLegend* legendG = new TLegend(0.6, 0.7, 0.88, 0.88); // shifted left
//legendG->SetTextSize(0.03);  // smaller, balanced text
//legendG->SetBorderSize(0);   // optional: removes border
//legendG->SetFillStyle(0);    // optional: transparent background
//legendG->AddEntry(eMinusAllTriggerableEndMomHist, "All", "l");
//legendG->AddEntry(eMinusRank0TriggerableEndMomHist, "Rank 0", "l");
//legendG->Draw();
//c1->SaveAs("eMinusTriggerableAllAndRank0EndMomHist.pdf");
//c1->Clear();
//eMinusAllTriggeredEndMomHist->SetTitle("Triggered 'all' and 'rank 0' eMinus End Momentum Comparison (MeV/c)");
//eMinusAllTriggeredEndMomHist->SetLineColor(kBlue);
//eMinusRank0TriggeredEndMomHist->SetLineColor(kRed);
//eMinusAllTriggeredEndMomHist->Draw("HIST");
//eMinusRank0TriggeredEndMomHist->Draw("HIST SAME");
//TLegend* legendH = new TLegend(0.6, 0.7, 0.88, 0.88); // shifted left
//legendH->SetTextSize(0.03);  // smaller, balanced text
//legendH->SetBorderSize(0);   // optional: removes border
//legendH->SetFillStyle(0);    // optional: transparent background
//legendH->AddEntry(eMinusAllTriggeredEndMomHist, "All", "l");
//legendH->AddEntry(eMinusRank0TriggeredEndMomHist, "Rank 0", "l");
//legendH->Draw();
//c1->SaveAs("eMinusTriggeredAllAndRank0EndMomHist.pdf");
//c1->Clear();

cout << "Creating and saving the combined histogram for triggered and triggerable for all eMinus in a track and specifically rank 0 eMinus." << endl;
eMinusAllTriggeredStartMomHist->SetTitle("Triggerable and Triggered eMinus Signals for Start Momentum for 'all' and 'rank 0' electrons (MeV/c).");
eMinusAllTriggerableStartMomHist->SetLineColor(kBlue);
eMinusAllTriggeredStartMomHist->SetLineColor(kRed);
eMinusRank0TriggerableStartMomHist->SetLineColor(kBlack);
eMinusRank0TriggeredStartMomHist->SetLineColor(kGreen);
eMinusAllTriggerableStartMomHist->Draw("HIST");
eMinusAllTriggeredStartMomHist->Draw("HIST SAME");
eMinusRank0TriggerableStartMomHist->Draw("HIST SAME");
eMinusRank0TriggeredStartMomHist->Draw("HIST SAME");
TLegend* legendI = new TLegend(0.30, 0.7, 0.58, 0.88);
legendI->SetTextSize(0.03);
legendI->SetBorderSize(0);
legendI->SetFillStyle(0);
legendI->AddEntry(eMinusAllTriggerableStartMomHist,"Triggerable All","l");
legendI->AddEntry(eMinusAllTriggeredStartMomHist,"Triggered All","l");
legendI->AddEntry(eMinusRank0TriggerableStartMomHist,"Triggerable Rank 0","l");
legendI->AddEntry(eMinusRank0TriggeredStartMomHist,"Triggered Rank 0","l");
legendI->Draw();
c1->SaveAs("eMinusFullParameterStartMomHist.pdf");
c1->Clear();
//eMinusAllTriggeredEndMomHist->SetTitle("Triggerable and Triggered eMinus Signals for End Momentum for 'all' and 'rank 0' electrons (MeV/c).");
//eMinusAllTriggerableEndMomHist->SetLineColor(kBlue);
//eMinusAllTriggeredEndMomHist->SetLineColor(kRed);
//eMinusRank0TriggerableEndMomHist->SetLineColor(kBlack);
//eMinusRank0TriggeredEndMomHist->SetLineColor(kGreen);
//eMinusAllTriggerableEndMomHist->Draw("HIST");
//eMinusAllTriggeredEndMomHist->Draw("HIST SAME");
//eMinusRank0TriggerableEndMomHist->Draw("HIST SAME");
//eMinusRank0TriggeredEndMomHist->Draw("HIST SAME");
//TLegend* legendJ = new TLegend(0.6, 0.7, 0.88, 0.88);
//legendJ->SetTextSize(0.03);
//legendJ->SetBorderSize(0);
//legendJ->SetFillStyle(0);
//legendJ->AddEntry(eMinusAllTriggerableEndMomHist,"Triggerable All","l");
//legendJ->AddEntry(eMinusAllTriggeredEndMomHist,"Triggered All","l");
//legendJ->AddEntry(eMinusRank0TriggerableEndMomHist,"Triggerable Rank 0","l");
//legendJ->AddEntry(eMinusRank0TriggeredEndMomHist,"Triggered Rank 0","l");
//legendJ->Draw();
//c1->SaveAs("eMinusFullParameterEndMomHist.pdf");
c1->Clear();
}
