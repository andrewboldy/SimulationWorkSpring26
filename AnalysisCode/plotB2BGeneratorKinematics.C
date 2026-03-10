//----------------------------------------------------------------------------------

//plotB2BGeneratorKinematics(string filelist)
//Written by Andrew Boldy + Codex
//University of South Carolina
//Spring 2026

//----------------------------------------------------------------------------------

//Plots momentum and position distributions for generator-produced particles only
//from DTS art files. Uses SimParticleCollection branch directly from Events tree.

//Usage: root -l 'plotB2BGeneratorKinematics.C++("filelist")'

//----------------------------------------------------------------------------------

//My Inclusions

//Standard Inclusions
#include <string>
#include <fstream>
#include <iostream>
#include <vector>

//CERN ROOT Inclusions
#include "TCanvas.h"
#include "TChain.h"
#include "TClass.h"
#include "TH1F.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"

//Mu2e Inclusions
#include "Offline/MCDataProducts/inc/SimParticle.hh"

//Personal Inclusions (if any)

//Namespace
using namespace std;

void plotB2BGeneratorKinematics(const string& filelist,
                                const string& simParticleBranch = "mu2e::SimParticlemv_compressDetStepMCs__Primary.obj")
{
  //Load dictionary if needed (when running as ROOT macro)
  if (!TClass::GetClass("mu2e::SimParticle"))
  {
    gSystem->Load("libMCDataProducts");
  }

  //Read input filelist
  ifstream listOfFiles(filelist);
  if (!listOfFiles.good())
  {
    cout << "ERROR: Could not open filelist: " << filelist << endl;
    return;
  }

  vector<string> inputFiles;
  string line;
  while (getline(listOfFiles, line))
  {
    if (!line.empty())
    {
      inputFiles.push_back(line);
    }
  }

  if (inputFiles.size() == 0)
  {
    cout << "ERROR: No files found in filelist." << endl;
    return;
  }

  //Build chain of Events trees
  TChain* eventsChain = new TChain("Events");
  for (size_t i_file = 0; i_file < inputFiles.size(); i_file++)
  {
    eventsChain->Add(inputFiles[i_file].c_str());
  }

  long long numEntries = eventsChain->GetEntries();
  if (numEntries <= 0)
  {
    cout << "ERROR: Events chain is empty." << endl;
    return;
  }

  //Initialize histograms and canvases
  TH1F* hP = new TH1F("hP", "Generator Primaries: Momentum Magnitude;|p| [MeV/c];Counts", 120, 0.0, 120.0);
  TH1F* hPx = new TH1F("hPx", "Generator Primaries: p_{x};p_{x} [MeV/c];Counts", 120, -120.0, 120.0);
  TH1F* hPy = new TH1F("hPy", "Generator Primaries: p_{y};p_{y} [MeV/c];Counts", 120, -120.0, 120.0);
  TH1F* hPz = new TH1F("hPz", "Generator Primaries: p_{z};p_{z} [MeV/c];Counts", 120, -120.0, 120.0);

  TH1F* hX = new TH1F("hX", "Generator Primaries: Start x;x [mm];Counts", 120, -4000.0, 4000.0);
  TH1F* hY = new TH1F("hY", "Generator Primaries: Start y;y [mm];Counts", 120, -4000.0, 4000.0);
  TH1F* hZ = new TH1F("hZ", "Generator Primaries: Start z;z [mm];Counts", 160, -8000.0, 8000.0);

  //Loop through entries and SimParticles
  TTreeReader reader(eventsChain);
  TTreeReaderValue<mu2e::SimParticleCollection> simParticles(reader, simParticleBranch.c_str());

  int primaryCounter = 0;
  int electronCounter = 0;
  int muMinusCounter = 0;
  int otherCounter = 0;

  long long i_event = -1;
  while (reader.Next())
  {
    i_event++;

    for (mu2e::SimParticleCollection::const_iterator i_particle = simParticles->begin();
         i_particle != simParticles->end();
         i_particle++)
    {
      const mu2e::SimParticle& sim = i_particle->second;

      //Generator-produced particles only
      if (!sim.isPrimary())
      {
        continue;
      }

      primaryCounter++;

      int pdg = (int)sim.pdgId();
      if (pdg == 11) electronCounter++;
      else if (pdg == 13) muMinusCounter++;
      else otherCounter++;

      double px = sim.startMomXYZT().x();
      double py = sim.startMomXYZT().y();
      double pz = sim.startMomXYZT().z();
      double p = sim.startMomXYZT().R();

      double x = sim.startPosXYZ().x();
      double y = sim.startPosXYZ().y();
      double z = sim.startPosXYZ().z();

      hP->Fill(p);
      hPx->Fill(px);
      hPy->Fill(py);
      hPz->Fill(pz);
      hX->Fill(x);
      hY->Fill(y);
      hZ->Fill(z);

      //Simple printout for quick validation
      if (primaryCounter <= 25)
      {
        cout << "Entry: " << i_event
             << " PDG: " << pdg
             << " p=( " << px << ", " << py << ", " << pz << " )"
             << " |p|=" << p
             << " start=( " << x << ", " << y << ", " << z << " )" << endl;
      }
    }
  }

  cout << "Finished processing " << numEntries << " entries from " << inputFiles.size() << " files." << endl;
  cout << "Generator-primary SimParticles counted: " << primaryCounter << endl;
  cout << "  electrons (pdg 11): " << electronCounter << endl;
  cout << "  mu- (pdg 13): " << muMinusCounter << endl;
  cout << "  other pdg: " << otherCounter << endl;

  //Draw and save momentum plots
  TCanvas* cMom = new TCanvas("cMom", "cMom", 1400, 900);
  cMom->Divide(2,2);
  cMom->cd(1); hP->Draw();
  cMom->cd(2); hPx->Draw();
  cMom->cd(3); hPy->Draw();
  cMom->cd(4); hPz->Draw();
  cMom->SaveAs("B2BCeEndpointDiffEnergy_GeneratorMomentum.pdf");

  //Draw and save position plots
  TCanvas* cPos = new TCanvas("cPos", "cPos", 1400, 900);
  cPos->Divide(2,2);
  cPos->cd(1); hX->Draw();
  cPos->cd(2); hY->Draw();
  cPos->cd(3); hZ->Draw();
  cPos->SaveAs("B2BCeEndpointDiffEnergy_GeneratorPosition.pdf");
}
