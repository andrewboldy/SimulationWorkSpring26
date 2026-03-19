//----------------------------------------------------------------------------------

//GenBasicTTree(string inputFile, bool skipFirstLine = false)
//Written by Andrew Boldy + Codex
//University of South Carolina
//Spring 2026

//----------------------------------------------------------------------------------

//Reads through a parsed Generator Output and creates a TTree holding all of the information for a generator that produces two particles.
//The input file is expected to have two lines per event: momentum line 1 and momentum line 2.
//The first line can optionally be skipped if it is a header.

//----------------------------------------------------------------------------------

//Parts of the Structure for each event:
//EventNumber (Which Event we are looking at)
//MomNumber1 (First momentum line number: 1)
//CartesianVector1 (First ThreeVector represented in Cartesian Coordinates)
//SphericalVector1 (First ThreeVector represented in Spherical Coordinates)
//FourMag1 (Magnitude of the first 4 vector)
//Energy1 (Energy magnitude of the first electron)
//MomNumber2 (Second momentum line number: 2)
//CartesianVector2 (Second ThreeVector represented in Cartesian Coordinates)
//SphericalVector2 (Second ThreeVector represented in Spherical Coordinates)
//FourMag2 (Magnitude of the second 4 vector)
//Energy2 (Energy magnitude of the second electron)
//b2b (True if vectors are back-to-back within tolerance)

//----------------------------------------------------------------------------------

//Example of what we are reading
//MomNumber: 1 || Cartesian vector: (43.7989, 15.8803, 24.1679) || Spherical Vector: (52.4844, 1.09226, 0.347831) || Magnitude of 4 Vector: 0.510999 || Energy: 52.4869
//MomNumber: 2 || Cartesian vector: (-43.7989, -15.8803, -24.1679) || Spherical Vector: (52.4844, 2.04933, -2.79376) || Magnitude of 4 Vector: 0.510999 || Energy: 52.4869

//----------------------------------------------------------------------------------

//Also want to print out and show the generator outputs and do a true or false value to determine if they are truly back to back then plot specific colored points as a histogram
//lets say red for "forward" and blue for "backwards" directions. May need to consult jackson for calculating which ones would make it to the detector effectively.

//----------------------------------------------------------------------------------

//My Inclusions

//Standard Inclusions
#include <string>
#include <fstream>
#include <iostream>
#include <sstream>
#include <vector>
#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstdlib>

//CERN ROOT Inclusions
#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TCanvas.h>
#include <TLegend.h>

//CLHEP Inclusions (If I'm feeling spicy)

//Mu2e Inclusions

//----------------------------------------------------------------------------------

//Using statements for readability
using std::string;
using std::ifstream;
using std::cout;
using std::cerr;
using std::endl;
using std::vector;

//----------------------------------------------------------------------------------
//Helper to extract all numbers from a line (handles negatives and decimals)
//This function is intentionally dumb and flexible: it scans the line and grabs any
//floating-point number it finds in order of appearance.
//----------------------------------------------------------------------------------
static vector<double> extractNumbers(const string& line)
{
  vector<double> values; //vector to store the numbers that are found in the line
  const char* cstr = line.c_str(); //raw char pointer for scanning the line
  char* endptr = nullptr; //end pointer used by strtod to locate where the number stops

  while (*cstr != '\0')
  {
    //Skip to start of a number
    if ((*cstr >= '0' && *cstr <= '9') || *cstr == '-' || *cstr == '.')
    {
      double val = std::strtod(cstr, &endptr);
      if (endptr != cstr)
      {
        values.push_back(val);
        cstr = endptr;
        continue;
      }
    }
    ++cstr;
  }

  return values;
}

//----------------------------------------------------------------------------------
//Main method
//inputFile: path to the parsed generator output
//skipFirstLine: set to true if the first line of the file is a header
//----------------------------------------------------------------------------------
void GenBasicTTree(const string& inputFile, bool skipFirstLine = false)
{
  //---------------------------------------------------------------------------------------------------------------------------
  //Open input file and verify it is readable
  //---------------------------------------------------------------------------------------------------------------------------
  ifstream inFile(inputFile);
  if (!inFile)
  {
    cerr << "ERROR: Cannot open input file: " << inputFile << endl;
    return;
  } //end opening file check

  //---------------------------------------------------------------------------------------------------------------------------
  //Initialize TTree and output file
  //---------------------------------------------------------------------------------------------------------------------------
  TFile *file = new TFile("GenOutputs.root","RECREATE"); //create (or recreate) the output root file for generators
  TTree *GeneratorBasicTTree = new TTree("GeneratorBasicTTree","Basic Mom TTree"); //create the TTree that will store information for the generators

  //---------------------------------------------------------------------------------------------------------------------------
  //Declare variables that will become branches
  //---------------------------------------------------------------------------------------------------------------------------
  int event;
  int momNumber1;
  int momNumber2;
  vector<double> cartesianVector1;
  vector<double> sphericalVector1;
  double fourMag1;
  double energy1;
  vector<double> cartesianVector2;
  vector<double> sphericalVector2;
  double fourMag2;
  double energy2;
  bool b2b;
  double cosTheta;

  //---------------------------------------------------------------------------------------------------------------------------
  //Create branches in the TTree
  //---------------------------------------------------------------------------------------------------------------------------
  GeneratorBasicTTree->Branch("event", &event, "event/I");
  GeneratorBasicTTree->Branch("momNumber1", &momNumber1, "momNumber1/I");
  GeneratorBasicTTree->Branch("momNumber2", &momNumber2, "momNumber2/I");
  GeneratorBasicTTree->Branch("cartesianVector1", &cartesianVector1);
  GeneratorBasicTTree->Branch("sphericalVector1", &sphericalVector1);
  GeneratorBasicTTree->Branch("fourMag1", &fourMag1, "fourMag1/D");
  GeneratorBasicTTree->Branch("energy1", &energy1, "energy1/D");
  GeneratorBasicTTree->Branch("cartesianVector2", &cartesianVector2);
  GeneratorBasicTTree->Branch("sphericalVector2", &sphericalVector2);
  GeneratorBasicTTree->Branch("fourMag2", &fourMag2, "fourMag2/D");
  GeneratorBasicTTree->Branch("energy2", &energy2, "energy2/D");
  GeneratorBasicTTree->Branch("cosTheta", &cosTheta, "cosTheta/D");
  GeneratorBasicTTree->Branch("b2b", &b2b, "b2b/O");

  //---------------------------------------------------------------------------------------------------------------------------
  //Initialize histograms for quick diagnostic plots
  //---------------------------------------------------------------------------------------------------------------------------
  TH2D* hForward = new TH2D("hForward", "Forward electrons;#theta [rad];#phi [rad]", 100, 0.0, M_PI, 100, -M_PI, M_PI);
  TH2D* hBackward = new TH2D("hBackward", "Backward electrons;#theta [rad];#phi [rad]", 100, 0.0, M_PI, 100, -M_PI, M_PI);
  TH1D* hCosTheta = new TH1D("hCosTheta", "Back-to-back check;cos(#theta);Counts", 200, -1.1, 1.1);
  TH1D* hOpeningAngle = new TH1D("hOpeningAngle", "Opening angle;cos(#theta_{12});Counts", 200, -1.0, 1.0);
  TH2D* hCosThetaPhiAll = new TH2D("hCosThetaPhiAll", "Single-electron angular distribution;cos(#theta);#phi [rad]", 200, -1.0, 1.0, 100, -M_PI, M_PI);

  //---------------------------------------------------------------------------------------------------------------------------
  //Read the input file, two lines per event
  //---------------------------------------------------------------------------------------------------------------------------
  string line1;
  string line2;
  int i_event = 0;
  int b2bCount = 0;

  //Optional skip of the first line (used for header line cases)
  if (skipFirstLine)
  {
    if (!std::getline(inFile, line1))
    {
      cerr << "WARNING: Input file is empty after skipping first line." << endl;
      return;
    }
  }

  //---------------------------------------------------------------------------------------------------------------------------
  //Main loop over the input file
  //---------------------------------------------------------------------------------------------------------------------------
  while (std::getline(inFile, line1))
  {
    //Try to read the second line of the event; if missing, stop.
    if (!std::getline(inFile, line2))
    {
      cerr << "WARNING: Unpaired line at end of file. Stopping." << endl;
      break;
    }

    //Extract all numeric values from both lines
    vector<double> nums1 = extractNumbers(line1);
    vector<double> nums2 = extractNumbers(line2);

    //Sanity check: we expect 9 numbers per line for this format
    if (nums1.size() < 9 || nums2.size() < 9)
    {
      cerr << "WARNING: Could not parse line pair for event " << i_event << endl;
      continue;
    }

    //-------------------------------------------------------------------------------------------------------------------------
    //Map the values into their physical meanings
    //-------------------------------------------------------------------------------------------------------------------------
    event = i_event;
    momNumber1 = static_cast<int>(nums1[0]);
    momNumber2 = static_cast<int>(nums2[0]);

    cartesianVector1 = {nums1[1], nums1[2], nums1[3]};
    sphericalVector1 = {nums1[4], nums1[5], nums1[6]};
    fourMag1 = nums1[7];
    energy1 = nums1[8];

    cartesianVector2 = {nums2[1], nums2[2], nums2[3]};
    sphericalVector2 = {nums2[4], nums2[5], nums2[6]};
    fourMag2 = nums2[7];
    energy2 = nums2[8];

    //-------------------------------------------------------------------------------------------------------------------------
    //Compute the opening angle between the two 3-vectors
    //-------------------------------------------------------------------------------------------------------------------------
    double p1 = std::sqrt(cartesianVector1[0]*cartesianVector1[0] +
                          cartesianVector1[1]*cartesianVector1[1] +
                          cartesianVector1[2]*cartesianVector1[2]);
    double p2 = std::sqrt(cartesianVector2[0]*cartesianVector2[0] +
                          cartesianVector2[1]*cartesianVector2[1] +
                          cartesianVector2[2]*cartesianVector2[2]);
    double dot = cartesianVector1[0]*cartesianVector2[0] +
                 cartesianVector1[1]*cartesianVector2[1] +
                 cartesianVector1[2]*cartesianVector2[2];

    //If both magnitudes are valid, compute cosine of the opening angle
    cosTheta = (p1 > 0.0 && p2 > 0.0) ? (dot / (p1 * p2)) : 1.0;

    //-------------------------------------------------------------------------------------------------------------------------
    //Back-to-back test and histogram fills
    //-------------------------------------------------------------------------------------------------------------------------
    const double b2bTolerance = 1.0e-3; //cosTheta <= -1 + tol
    b2b = (cosTheta <= (-1.0 + b2bTolerance));
    if (b2b)
    {
      b2bCount++;
    }

    hCosTheta->Fill(cosTheta);
    hOpeningAngle->Fill(cosTheta);

    //Fill angular distribution for each electron separately
    hCosThetaPhiAll->Fill(std::cos(sphericalVector1[1]), sphericalVector1[2]);
    hCosThetaPhiAll->Fill(std::cos(sphericalVector2[1]), sphericalVector2[2]);

    //-------------------------------------------------------------------------------------------------------------------------
    //Forward/backward definition based on pz of electron 1 (adjust if needed)
    //-------------------------------------------------------------------------------------------------------------------------
    if (cartesianVector1[2] >= 0.0)
    {
      hForward->Fill(sphericalVector1[1], sphericalVector1[2]);
    }
    else
    {
      hBackward->Fill(sphericalVector1[1], sphericalVector1[2]);
    }

    //-------------------------------------------------------------------------------------------------------------------------
    //Print the first few events for spot checks
    //-------------------------------------------------------------------------------------------------------------------------
    if (i_event < 10)
    {
      cout << "Event " << i_event << endl;
      cout << "  " << line1 << endl;
      cout << "  " << line2 << endl;
      cout << "  cosTheta = " << cosTheta << "  b2b = " << b2b << endl;
    }

    //-------------------------------------------------------------------------------------------------------------------------
    //Fill the tree and move to the next event
    //-------------------------------------------------------------------------------------------------------------------------
    GeneratorBasicTTree->Fill();
    i_event++;
  } //finish reading and filling

  //---------------------------------------------------------------------------------------------------------------------------
  //Write out the TTree and plots
  //---------------------------------------------------------------------------------------------------------------------------
  GeneratorBasicTTree->Write();

  //Draw quick diagnostic plot for forward/backward
  TCanvas* cDir = new TCanvas("cDir", "Forward/Backward", 1200, 600);
  cDir->Divide(2,1);
  cDir->cd(1);
  hForward->SetMarkerColor(kRed + 1);
  hForward->SetMarkerStyle(20);
  hForward->Draw("P");
  cDir->cd(2);
  hBackward->SetMarkerColor(kBlue + 1);
  hBackward->SetMarkerStyle(20);
  hBackward->Draw("P");
  cDir->SaveAs("GenBasicTTree_ForwardBackward.pdf");

  //Back-to-back cosine distribution
  TCanvas* cB2B = new TCanvas("cB2B", "Back-to-Back", 800, 600);
  hCosTheta->Draw();
  cB2B->SaveAs("GenBasicTTree_CosTheta.pdf");

  //Opening angle histogram
  TCanvas* cOpeningAngle = new TCanvas("cOpeningAngle", "Opening Angle", 900, 700);
  hOpeningAngle->Draw();
  cOpeningAngle->SaveAs("GenBasicTTree_OpeningAngle.pdf");

  //Single-electron angular distribution
  TCanvas* cCosThetaPhiAll = new TCanvas("cCosThetaPhiAll", "All Electrons CosTheta vs Phi", 900, 700);
  hCosThetaPhiAll->Draw("COLZ");
  cCosThetaPhiAll->SaveAs("GenBasicTTree_AllElectrons_CosThetaPhi.pdf");

  //Print event counts
  cout << "Total events processed: " << i_event << endl;
  cout << "Back-to-back pairs: " << b2bCount << endl;

  //---------------------------------------------------------------------------------------------------------------------------
  //Close the output file and clean up
  //---------------------------------------------------------------------------------------------------------------------------
  file->Close();

  delete file;
  return;
} //end GenBasicTTree method

