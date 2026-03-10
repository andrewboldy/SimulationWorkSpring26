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
#include <TFile.h>
#include <TTree.h>
#include <TCanvas.h>
#include <TH1F.h>
#include <TLegend.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <TChain.h>

//Mu2e Inclusions
#include "Offline/MCDataProducts/inc/GenParticle.hh"


  //Using statements for readability 
  using std::string;
  using std::ifstream;
  using std::cout;
  using std::endl;
  using std::vector;
    

void dtsDisplayInfo(const string& fileName, const string& genParticleBranch = "mu2e::GenParticles_compressDetStepMCs__Primary.obj")
{

  //---------------------------------------------------------------------------------------------------------------------------
  //Checking Filelist and getting information about the trees before actively digging into each one
  //---------------------------------------------------------------------------------------------------------------------------
  ifstream listOfFiles(fileName);
  if (!listOfFiles.good())
  {
    cout << "ERROR: Could not open filelist: " << fileName << endl;
    return;
  } //end file opening check

  vector<string> inputFiles;
  string line;
  while (getline(listOfFiles, line))
  {
    cout << "Adding " << line << " to the vector of input files." << endl;
    inputFiles.push_back(line);
  } //end iteration over the list of files and addition to the vector of files

  if (inputFiles.size() == 0)
  {
    cout << "ERROR: input filelist appears to be empty." << endl;
    return;
  } //end empty inputFile vector check

  TChain* eventsChain = new TChain("Events"); //this defines a TChain we call eventsChain that hits every Events TTree in each of the files in the list
  for (size_t i_file = 0; i_file < inputFiles.size(); i_file++)
  {
    eventsChain->Add(inputFiles[i_file].c_str()); //adding TTrees to the TChain
  } //end addition of TTrees to the TChain

  long long nEntries = eventsChain->GetEntries();

  if (nEntries <= 0)
  {
    cout << "ERROR: Events TChain is empty." << endl;
    return;
  } //end empty chain check

  //---------------------------------------------------------------------------------------------------------------------------
  // Initialize Variables of Interest for the Actual Plotting
  //---------------------------------------------------------------------------------------------------------------------------

  //Initialize Counters
  int primaryCounter = 0; //initialize a counter that will count Primaries
  int electronCounter = 0; //initialize a counter that will count the number of electrons
  int muMinusCounter = 0; //initialize a counter that will count the number of muons
  int otherParticleCounter = 0; //initialize a counter that will count particles that are not electrons or muMinuses
  long long i_event = -1;

  //Initialize Histograms
  TH1F* hThrownP = new TH1F("hThrownP", "Generator Primaries: Magnitude of the Thrown Momentum;|p| [MeV/c]; Counts", 120, 0.0, 120.0);
  TH1F* hThrownPx = new TH1F("hThrownPx", "Generator Primaries: X Component of the Thrown Momentum p_{x};p_{x}[MeV/c]; Counts", 120, -120.0, 120.0);
  TH1F* hThrownPy = new TH1F("hThrownPy", "Generator Primaries: Y Component of the Thrown Momentum p_{y};p_{y} [MeV/c]; Counts", 120, -120.0, 120.0);
  TH1F* hThrownPz = new TH1F("hThrownPz", "Generator Primaries: Z Component of the Thrown Momentum p_{z};p_{z}[MeV/c]; Counts", 120, -120.0, 120.0);

  TH1F* hStartX = new TH1F("hStartX", "Generator Primaries: X Component of the Starting Position; x [mm] Counts", 120, -4000.0, 4000.0);
  TH1F* hStartY = new TH1F("hStartY", "Generator Primaries: Y Component of the Starting Position; y [mm]; Counts", 120, -4000.0, 4000.0);
  TH1F* hStartZ = new TH1F("hStartZ", "Generator Primaries: Z Component of the Starting Position; z [mm]; Counts", 120, -4000.0, 4000.0);

  //initialize TTreeReader information
  TTreeReader reader(eventsChain); //initializing a TTreeReader called reader that reads through eventsChain
  TTreeReaderValue<mu2e::GenParticleCollection> genParts(reader, genParticleBranch.c_str()); //setting the target of the reader and calling this specific handle genParts

  //Initialize the
  vector<mu2e::GenParticle> particles;
  //---------------------------------------------------------------------------------------------------------------------------
  // Begin Reading
  //---------------------------------------------------------------------------------------------------------------------------

  while (reader.Next())
  {
    i_event++; //increment at the start of the while loop
    cout << "======================================================= \n";
    cout << "Event number: " << i_event << endl;
    for (auto const& sim : *genParts)
    {
      if (sim.pdgId()==11) electronCounter++;
      else if (sim.pdgId()==13) muMinusCounter++;
      else otherParticleCounter++;

      cout << "Particle ID: " << sim.pdgId()
      << " Process ID: " << sim.generatorId()
      << " GenId name: " << sim.generatorId().name()
      << " Momentum: " << sim.momentum()
      << " Position: " << sim.position() << endl;

      particles.push_back(sim);
    } //end the for loop going through each of the particles in the GenParticle colleciton that we have created
  } //end reading through each event in eventsChain
} //end dtsDisplayInfo

