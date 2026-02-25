//----------------------------------------------------------------------------------

//plotElectronMomHist(string filelist)
//Written by Andrew Boldy
//University of South Carolina
//Spring 2025

//----------------------------------------------------------------------------------

//Plots the electron momentum for a given non-mixed generator filelist
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

void plotElectronMomHists(const string &generatorName, const string &filelist)
{
  //Initilizing the histograms and the TCanvas
  TCanvas* c1 = new TCanvas("c1","c1");
  string histTitle = "Momentum (MeV/c) of Electrons Produced by " + generatorName + " Generator";
  TH1F* ElectronMomHist = new TH1F("ElectronMomHist", "Momentum (MeV/c) of Electrons Produced by the Generator", 50, 0, 115);

  //begin the filelist reading
  ifstream listOfFiles(filelist);
  int fileCount = 0;
  string line;

  //counter for files in the list
  while (getline(listOfFiles,line))
  {
    if (!line.empty())
    {
      fileCount++;
    } //end if line empty
  } //end while getLine

  //open up the RooUtil and check number of events contained in the filelist
  RooUtil util(filelist);
  int numEvents = util.GetNEvents();
  int numElectrons = 0; 
  cout << "There are " << numEvents << " in the filelist." << endl;

  //Dig into the actual TTrees 
  for (int i_event = 0; i_event < numEvents; i_event++)
  {
    auto& event = util.GetEvent(i_event);
    const auto& e_minus_tracks = event.GetTracks(is_e_minus);

    for (auto &track : e_minus_tracks)
    {
      for (const auto &mctrack : *(track.trkmcsim))
      {
        if (mctrack.pdg == 11 && mctrack.rank == 0)
        {
          ElectronMomHist->Fill(mctrack.mom.R()); //fill the electron momentum hist with momenta for each track 
          numElectrons++;
        } //end work on the electron pdg tracks 
      }//end the interface with the TTree at trkmcsim branch level
    }//end the interface with the TTree at track level
  } //end the interface with the TTree at event level


//Print out results: 
cout << "There are " << numElectrons << " electrons in " << numEvents << " events for this dataset." << endl;

//Draw and save histograms
ElectronMomHist->Draw();
string pdfName = "ElectronMomentumHistogram" + generatorName + ".pdf";
c1->SaveAs(pdfName.c_str());
} //end plotElectronMomHist method

