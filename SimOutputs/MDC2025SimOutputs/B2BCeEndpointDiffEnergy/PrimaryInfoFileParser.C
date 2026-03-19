//----------------------------------------------------------------------------------

//PrimaryInfoFileParser(string filelist)
//Written by Andrew Boldy
//University of South Carolina
//Spring 2025

//----------------------------------------------------------------------------------

//Parses through the output of a .txt file to make sure that the output does not have extraneous information.

//----------------------------------------------------------------------------------

//My Inclusions

//Standard Inclusions
#include <string>
#include <fstream>
#include <iostream>
#include <vector>

//CERN ROOT Inclusions


//Mu2e Inclusions



  //Using statements for readability
  using std::string;
  using std::ifstream;
  using std::ofstream;
  using std::cout;
  using std::cerr;
  using std::endl;
  using std::vector;

void PrimaryInfoFileParser(const string& inputFile, const string& outputFile)
{
  ifstream inFile(inputFile);
  ofstream outFile(outputFile);
  
  if (!inFile) {
    cerr << "Error: Could not open input file: " << inputFile << endl;
    return;
  } //end the check if infile was located correctly
  if (!outFile) {
    cerr << "Error: Could not open output file: " << outputFile << endl;
    return; 
  } //end the check if outfile was located correctly
  
  int lineCount=0;
  int cutLineCount=0;
  int remainingLineCount=0;
  string inLine;
  outFile << "Printing the following information for each pair of Lorentz vectors being added to the stack: " 
  << "momNumber (1 or 2) || threeMom vector Cartesian (x,y,z) || threeMomSpherical (r,theta,phi) || magnitude (of 4 vector) || Energy (MeV/c)" << endl;

  while (std::getline(inFile,inLine))
  {
    ++lineCount; //using ++Variable for preIncrement, though should look the same
    if (inLine.rfind("MomNumber: ", 0) == 0)
    {
      outFile << inLine << endl;
      ++remainingLineCount;
    }
    else 
    {
      ++cutLineCount;
      continue;
    }
  }
  cout << "Parsed " << lineCount << " lines: kept " << remainingLineCount
       << ", removed " << cutLineCount << "." << endl;
} //end PrimaryInfoFileParser
