/*-----------------------------------------------------------------------------------------------*\
License
  This file is part of AccelerateCFD_Community_Edition.

  AccelerateCFD_Community_Edition is free software: you can redistribute it 
  and/or modify it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or (at your option) 
  any later version.

  AccelerateCFD_Community_Edition is distributed in the hope that it will be useful, 
  but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
  FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with AccelerateCFD_Community_Edition.  If not, see <http://www.gnu.org/licenses/>.
  
Application
  podROM
  
Description
  This application reads data output of application "podPrecompute" and calculates
  time coefficients for POD reduced order model (POD-ROM). These time varying coefficients
  are then used to reconstruct the velocity by application "podFlowReconstruct" 
  
Author
  Illinois Rocstar LLC
  AccelerateCFD Development Team
  Copyright (C) 2017-2019
 
\*---------------------------------------------------------------------------------------------*/

#include <iostream>
#include <iomanip>
#include <math.h>
#include <string>
#include <sstream>
#include <fstream>
#include <vector>
#include <map>
#include <iterator>

using namespace std;

using vec    = vector<double>;
using matrix = vector<vec>;
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Function that reads CSV files line by line 
matrix readCSV( string filename )
{
  matrix M;

  ifstream in( filename );
  string line;
  while ( getline( in, line ) )
  {
    stringstream ss( line );
    vec row;
    string data;
    while ( getline( ss, data, ',' ) )
    {
      row.push_back( stod( data ) );
    }
    if ( row.size() > 0 ) M.push_back( row );
   }
   return M;
}

void copyrightnotice()
{
  std::cout << 
  "*******************************************************************************" 
  << std::endl;

  std::cout << 
  "*******************************************************************************" 
  << std::endl;

  std::cout << "* AccelerateCFD_Community_Edition" << std::endl;
  std::cout << "* Copyright (C) 2017-2019 Illinois Rocstar LLC" << std::endl;
  std::cout << "* GNU GENERAL PUBLIC LICENSE VERSION 3 (2007)" << std::endl;
  std::cout << "* www.Illinoisrocstar.com" << std::endl;

  std::cout << 
  "*******************************************************************************" 
  << std::endl;

  std::cout << 
  "*******************************************************************************" 
  << std::endl;
}

bool is_numeric(std::string &strng)
{
    int sizeOfString = strng.length();
    for (int i = 0; i < strng.length(); i++)
      if (isdigit(strng[i]) == false)
        return false;

      return true;
}

int main(int argc, char *argv[])
{

  copyrightnotice();

  if (argc > 2) {
    std::cerr << "Only two arguments allowed!" << std::endl;
    std::cout << "For Help --> " << argv[0] << " -h" << std::endl;
    throw;
  }

  std::vector<string> args;
  for (int i=0; i<argc; i++)
    args.push_back(argv[i]);

  int udfDim = 0;

  if ((args.size() > 1) && (args[1] == "-h")) {
    std::cout << "Usage: " << args[0] << std::endl; 
    std::cout << "For Help --> " << args[0] << " -h" << std::endl;
    std::cout << "Providing ROM dimension --> " << args[0] << " <num of modes>" << std::endl;
    return 0; 
  } else if ((args.size() > 1) && (is_numeric(args[1]))) {
    udfDim = std::atoi(argv[1]);
  }

  // To get processor clocktime
  std::clock_t start;
  double duration;
  start = std::clock();

  //Reading user defined values for ROM;
  std::vector<double> A(9,0.0);
  ifstream in("podInfo.csv");
  string line;
  int i = 0;
  while(getline(in,line))
  {
    double num = stod(line);
    A[i] = num;
    i++;
  }

  int nDim = (int) A[0];
  double nu = A[1];
  int writeFreq = (int) A[2];
  double tEnd = A[3];
  double dt = A[4];
  double nRows = A[5];
  double caseRunT = A[6];
  double numDirs = A[7];
  double startTime = A[8];

  if (udfDim != 0) {
    if (udfDim <= nDim)  {
      nDim = udfDim;
    }
    else {
      std::cerr << "User defined ROM dimension must not exceed " << nDim << "!"
                << std::endl;
      throw;
    }
  }

  cout << "Reading output from podPrecompute" << endl;
  matrix con = readCSV("constant.csv");
  matrix lin = readCSV("linear.csv");
  matrix quad = readCSV("quadratic.csv");
  matrix aprev = readCSV("prevVals.csv");

  std::vector<double> constant(con.size());
  std::vector<double> linear(lin.size());
  std::vector<double> quadratic(quad.size());
  std::vector<double> prevAvals(nDim);

  for (int i=0; i<con.size(); i++)
    constant[i] = con[i][1];

  for (int i=0; i<lin.size(); i++)
    linear[lin[i][0]] = lin[i][1];

  for (int i=0; i<quad.size(); i++)
    quadratic[quad[i][0]] = quad[i][1];

  for (int i=0; i<nDim; i++){
    prevAvals[i] = aprev[i][0];
  }

  double nSteps = (tEnd-startTime)/dt;  // total time steps to loop through
  double timeElapsed = startTime;       // time counter

  // 2D vector for storing a values at every time step  
  std::vector<std::vector<double>> avalues(nSteps+1, std::vector<double>(nDim,0.0));

  // 1D vector to store all a's in one time iteration
  std::vector<double> avals(nDim,0.0);

  // Now writing a values form 2D vector to CSV file for use in..
  // "podFlowReconstruct"
  int writeSteps = 0;

  if(writeFreq == 0){
    writeSteps = nSteps/numDirs;}
  else{
    writeSteps = writeFreq;}

  ofstream afiles;
  afiles.open("avals.csv");

  for (int t=0; t<nSteps+1; t++){
    cout << "t = " << timeElapsed << endl; // Case progress info in terminal
    for (int i=0; i<nDim; i++) {
      double da = constant[i];
      for (int j=0; j<nDim; j++) {
        da += linear[i+j*nDim]*prevAvals[j];
        for (int k=0; k<nDim; k++) {
          da += quadratic[i+j*nDim+k*nDim*nDim]*prevAvals[k]*prevAvals[j];
        }
      }
      avals[i] = prevAvals[i] + da*dt;
    }

    for (int i=0; i<nDim; i++){
      prevAvals[i] = avals[i]; // updating previous avals
      avalues[t][i] = avals[i]; // storing a values for in 2D vector for later
    }

    if(t%writeSteps == 0){
      double tcol = startTime + dt*t;
      afiles << tcol << ",";
      for (int j=0; j<nDim; j++){
        afiles << std::fixed << std::setprecision(16) << avalues[t][j] << ",";
      }
      afiles << endl << std::flush;
    }

    timeElapsed += dt; 
  }
  afiles.close();

  // processor clock time info displays when program ends
  duration = (std::clock() - start ) / (double) CLOCKS_PER_SEC;

  cout << "runtime = " << duration << " seconds" << endl;

  return 0;
}


// ************************************************************************* //
