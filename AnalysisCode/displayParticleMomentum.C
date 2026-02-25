//----------------------------------------------------------------------------------

//displayParticleMomentum(string filelist)
//Written by Andrew Boldy
//University of South Carolina
//Fall 2025

//----------------------------------------------------------------------------------

//Displays particle momentum for a list of files.
//Must be run using C++ so compile mode in order to rectify issues that are generated when using common_cuts.hh

//----------------------------------------------------------------------------------

//My Inclusions

//Standard Inclusions
#include <string>
#include <fstream>
#include <iostream>
#include <vector>

//CERN ROOT Inclusions

//Mu2e Inclusions
#include "EventNtuple/rooutil/inc/RooUtil.hh"
#include "EventNtuple/rooutil/inc/common_cuts.hh"

//Personal Inclusions (if any)

void displayParticleMomentum(string filelist)
{
  ifstream listOfFiles(filelist);
  int fileCount = 0;
  string line;
  while (getline(listOfFiles,line))
  {
    if (!line.empty())
    {
      fileCount++;
    }
  }

  //Starting values:

  int electronCounter = 0;
  int positronCounter = 0;
  int muMinusCounter = 0;
  int protonCounter = 0;

  float eMinusStartMomTotal = 0;
  float ePlusStartMomTotal = 0;
  float muMinusStartMomTotal = 0;
  float protonStartMomTotal = 0;
  float eMinusEndMomTotal = 0;
  float ePlusEndMomTotal = 0;
  float muMinusEndMomTotal = 0;
  float protonEndMomTotal = 0;

  //Starting Maxes (all zero)
  float protonMaxStartMom = 0.0f;
  float protonMaxEndMom   = 0.0f;
  float eMinusMaxStartMom = 0.0f;
  float eMinusMaxEndMom   = 0.0f;
  float ePlusMaxStartMom  = 0.0f;
  float ePlusMaxEndMom    = 0.0f;
  float muMaxStartMom     = 0.0f;
  float muMaxEndMom       = 0.0f;
  float deutMaxStartMom   = 0.0f;
  float deutMaxEndMom     = 0.0f;

  //Starting Mins (largest possible float, thank you internet for this format)
  const float BIG = numeric_limits<float>::max();

  float protonMinStartMom = BIG;
  float protonMinEndMom   = BIG;
  float eMinusMinStartMom = BIG;
  float eMinusMinEndMom   = BIG;
  float ePlusMinStartMom  = BIG;
  float ePlusMinEndMom    = BIG;
  float muMinStartMom     = BIG;
  float muMinEndMom       = BIG;
  float deutMinStartMom   = BIG;
  float deutMinEndMom     = BIG;

  int numBins = 100;

  RooUtil util(filelist);
  int numEvents = util.GetNEvents();
  cout << "There are " << numEvents << " entries. Printing PDG, process code, number of hits, and rank for each." << endl;
  for (int i_event = 0; i_event < numEvents; i_event++)
  {
    auto& event = util.GetEvent(i_event);
    auto e_minus_tracks = event.GetTracks(is_e_minus);
    for (auto& track : e_minus_tracks)
    {
      if (track.trkmcsim != nullptr)
      {
        for (const auto& mctrack : *(track.trkmcsim))
        {
          if (mctrack.pdg == 11)
          {
            eMinusStartMomTotal = eMinusStartMomTotal + mctrack.mom.R();
            eMinusEndMomTotal = eMinusEndMomTotal + mctrack.endmom.R();
              if (mctrack.mom.R() > eMinusMaxStartMom)
              {
                  eMinusMaxStartMom = mctrack.mom.R();
              }
              if (mctrack.endmom.R() > eMinusMaxEndMom)
              {
                  eMinusMaxEndMom = mctrack.endmom.R();
              }
              if (mctrack.mom.R() < eMinusMinStartMom)
              {
                  eMinusMinStartMom = mctrack.mom.R();
              }
              if (mctrack.endmom.R() < eMinusMinEndMom)
              {
                  eMinusMinEndMom = mctrack.endmom.R();
              }
              electronCounter++;
          }

          if (mctrack.pdg == -11)
          {
            ePlusStartMomTotal = ePlusStartMomTotal + mctrack.mom.R();
            ePlusEndMomTotal = ePlusEndMomTotal + mctrack.endmom.R();
              if (mctrack.mom.R() > ePlusMaxStartMom)
              {
                  ePlusMaxStartMom = mctrack.mom.R();
              }
              if (mctrack.endmom.R() > ePlusMaxEndMom)
              {
                  ePlusMaxEndMom = mctrack.endmom.R();
              }
              if (mctrack.mom.R() < ePlusMinStartMom)
              {
                  ePlusMinStartMom = mctrack.mom.R();
              }
              if (mctrack.endmom.R() < ePlusMinEndMom)
              {
                  ePlusMinEndMom = mctrack.endmom.R();
              }
              positronCounter++;
          }
          if (mctrack.pdg == 13)
          {
            muMinusStartMomTotal = muMinusStartMomTotal + mctrack.mom.R();
            muMinusEndMomTotal = muMinusEndMomTotal + mctrack.endmom.R();
              if (mctrack.mom.R() > muMaxStartMom)
              {
                  muMaxStartMom = mctrack.mom.R();
              }
              if (mctrack.endmom.R() > muMaxEndMom)
              {
                  muMaxEndMom = mctrack.endmom.R();
              }
              if (mctrack.mom.R() < muMinStartMom)
              {
                  muMinStartMom = mctrack.mom.R();
              }
              if (mctrack.endmom.R() < muMinEndMom)
              {
                  muMinEndMom = mctrack.endmom.R();
              }
              muMinusCounter++;
          }
          if (mctrack.pdg == 2212)
          {
            protonStartMomTotal = protonStartMomTotal + mctrack.mom.R();
            protonEndMomTotal = protonEndMomTotal + mctrack.endmom.R();
              if (mctrack.mom.R() > protonMaxStartMom)
              {
                  protonMaxStartMom = mctrack.mom.R();
              }
              if (mctrack.endmom.R() > protonMaxEndMom)
              {
                  protonMaxEndMom = mctrack.endmom.R();
              }
              if (mctrack.mom.R() < protonMinStartMom)
              {
                  protonMinStartMom = mctrack.mom.R();
              }
              if (mctrack.endmom.R() < protonMinEndMom)
              {
                  protonMinEndMom = mctrack.endmom.R();
              }
              protonCounter++;
          }
        }
      }
    }
  }
  float eMinusAverageStartMom = eMinusStartMomTotal / electronCounter;
  float eMinusAverageEndMom = eMinusEndMomTotal / electronCounter;
  float ePlusAverageStartMom = ePlusStartMomTotal / positronCounter;
  float ePlusAverageEndMom = ePlusEndMomTotal / positronCounter;
  float muAverageStartMom = muMinusStartMomTotal / muMinusCounter;
  float muAverageEndMom = muMinusEndMomTotal / muMinusCounter;
  float protonAverageStartMom = protonStartMomTotal / protonCounter;
  float protonAverageEndMom = protonEndMomTotal / protonCounter;
   
  cout << "Calculation complete!" << endl;
  cout << "Printing Results..." << endl;
  cout << "| Min Start | Max Start | Min End | Max End | Momenta by Particle in MeV/c:      --> Average: Start | End " << endl;
  cout << "Electron: | " << eMinusMinStartMom << " | " << eMinusMaxStartMom << " | " << eMinusMinEndMom << " | " << eMinusMaxEndMom << " | --> Average: " << eMinusAverageStartMom << " | " << eMinusAverageEndMom << endl;
  cout << "Positron: | " << ePlusMinStartMom << " | " << ePlusMaxStartMom << " | " << ePlusMinEndMom << " | " << ePlusMaxEndMom << " | --> Average: " << ePlusAverageStartMom << " | " << ePlusAverageEndMom << endl;
  cout << "Muon: | " << muMinStartMom << " | " << muMaxStartMom << " | " << muMinEndMom << " | " << muMaxEndMom << " | --> Average: " << muAverageStartMom << " | " << muAverageEndMom << endl;
  cout << "Proton: | " << protonMinStartMom << " | " << protonMaxStartMom << " | " << protonMinEndMom << " | "  << protonMaxEndMom << " | --> Average: " << protonAverageStartMom << " | " << protonAverageEndMom << endl;
}
