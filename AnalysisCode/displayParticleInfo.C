//----------------------------------------------------------------------------------

//displayParticleInfo(string filelist)
//Written by Andrew Boldy
//University of South Carolina
//Fall 2025

//----------------------------------------------------------------------------------

//Displays particle information for a list of files, specifically PDG number, number of hits on the detector, process code, and rank.
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

void displayParticleInfo(string filelist)
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
  int eMinusCounter = 0;
  int ePlusCounter = 0;
  int protonCounter = 0;
  int muMinusCounter = 0;
  int neutronCounter = 0;
  int piMinusCounter = 0;
  int piPlusCounter = 0;
  int piZeroCounter = 0;
  int lambdaCounter = 0;
  int sigmaPlusCounter = 0;
  int sigmaZeroCounter = 0;
  int sigmaMinusCounter = 0;
  int kShortCounter = 0;
  int kLongCounter = 0;
  int kMinusCounter = 0;
  int kPlusCounter = 0;
  int photonCounter = 0;
  int cascadeZeroCounter = 0;
  int cascadeMinusCounter = 0;
  vector<int> extraParticleIDArray;
  
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
        cout << "New track!" << endl;
        for (const auto& mctrack : *(track.trkmcsim))
        {
            cout << "Entry Number: " << i_event
                 << " PDG Number: " << mctrack.pdg
                 << " Tracker Hits: " << mctrack.nhits
                 << " Process Code: " << mctrack.startCode
                 << " Rank: " << mctrack.rank << endl;
                 
            if (mctrack.pdg == 11)
            {
              eMinusCounter++;
            }
            if (mctrack.pdg == -11)
            {
              ePlusCounter++;
            }
            if (mctrack.pdg == 2212)
            {
              protonCounter++;
            }
            if (mctrack.pdg == 2112)
            {
              neutronCounter++;
            }
            if (mctrack.pdg == 13)
            {
              muMinusCounter++;
            }
            if (mctrack.pdg == -211)
            {
              piMinusCounter++;
            }
            if (mctrack.pdg == 211)
            {
              piPlusCounter++;
            }
            if (mctrack.pdg == 111)
            {
              piZeroCounter++;
            }
            if (mctrack.pdg == 3122)
            {
              lambdaCounter++;
            }
            if (mctrack.pdg == 3222)
            {
              sigmaPlusCounter++;
            }
            if (mctrack.pdg == 3212)
            {
              sigmaZeroCounter++;
            }
            if (mctrack.pdg == 3112)
            {
              sigmaMinusCounter++;
            }
            if (mctrack.pdg == 310)
            {
              kShortCounter++;
            }
            if (mctrack.pdg == 130)
            {
              kLongCounter++;
            }
            if (mctrack.pdg == -321)
            {
              kMinusCounter++;
            }
            if (mctrack.pdg == 321)
            {
              kPlusCounter++;
            }
            if (mctrack.pdg == 22)
            {
              photonCounter++;
            }
            if (mctrack.pdg == 3322)
            {
              cascadeZeroCounter++;
            }
	    if (mctrack.pdg == 3312)
	    {
	      cascadeMinusCounter++;
	    }

            if (mctrack.pdg != 3312 && mctrack.pdg != 3322 && mctrack.pdg != 321 && mctrack.pdg != 22 && mctrack.pdg != -321 && mctrack.pdg != 310 && mctrack.pdg != 130 && mctrack.pdg != 3222 && mctrack.pdg != 3112 && mctrack.pdg != 3212 && mctrack.pdg != 3122 && mctrack.pdg != 13 && mctrack.pdg != 11 && mctrack.pdg != -11 && mctrack.pdg != 2212 && mctrack.pdg != 2112 && mctrack.pdg != -211 && mctrack.pdg != 211 && mctrack.pdg != 111)
            {
              extraParticleIDArray.push_back(mctrack.pdg);
            }
          }
      }
    }
  }
  if (extraParticleIDArray.size() == 0)
  {
    cout << "All particles accounted for!" << endl;
  }
  else
  {
    cout << "Extraneous particles pdgs: " << endl;
    for (size_t particle_i = 0; particle_i < extraParticleIDArray.size(); particle_i++)
    {
      cout << particle_i << ". Particle type id: " << extraParticleIDArray[particle_i] << endl; 
    }
  }
  cout << "Number of electrons in tracks: " << eMinusCounter << endl;
  cout << "Number of protons in tracks: " << protonCounter << endl;
  cout << "Number of neutrons in trakcs: " << neutronCounter << endl;
  cout << "Number of muons in tracks: " << muMinusCounter << endl;
  cout << "Number of Positive | Zero | Negative pi mesons in tracks: " << piPlusCounter << " | " << piZeroCounter << " | " << piMinusCounter << endl;
  cout << "Number of lambda baryons in tracks: " << lambdaCounter << endl;
  cout << "Number of Positive | Neutral | Negative sigma baryons in tracks: " << sigmaPlusCounter << " | " << sigmaZeroCounter << " | " << sigmaMinusCounter << endl;
  cout << "Number of positrons in tracks: " << ePlusCounter << endl;
  cout << "Number of Positive | Negative | Short-lived | Long-lived Kaons in tracks: " << kPlusCounter << " | " << kMinusCounter << " | " << kShortCounter << " | " << kLongCounter << endl;
  cout << "Number of photons in tracks: " << photonCounter << endl;
  cout << "Number of Zero | Minus Cascade Baryons in tracks: " << cascadeZeroCounter << " | " << cascadeMinusCounter << endl;
}
