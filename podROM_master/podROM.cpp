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

int main(int argc, char *argv[])
{

  copyrightnotice();

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
  double firstDir = A[8];


  cout << "Reading output from podPrecompute" << endl;
  matrix con = readCSV("constant.csv");
  matrix lin = readCSV("linear.csv");
  matrix quad = readCSV("quadratic.csv");
  matrix aprev = readCSV("prevVals.csv");

  // Converting these matrices into C++ vectors for ease of calculations.

  std::vector<double> prevAvals(nDim);
  std::vector<double> constant(nDim,0.0);
  std::vector<std::vector<double>> linear(nDim, std::vector<double>(nDim,0.0));
  std::vector<std::vector<std::vector<double>>> quadratic(nDim, std::vector<std::vector<double>>(nDim, std::vector<double>(nDim,0.0)));

  for (int i=0; i<nDim; i++){
      constant[i] = con[i][0];
    for (int j=0; j<nDim; j++){
      linear[i][j] = lin[i][j];
    }
  }

  int row = 0;
  for (int i=0; i<nDim; i++){
    for (int j=0; j<nDim; j++){
      for (int k=0; k<nDim; k++){
        quadratic[i][j][k] = quad[row][k];
      }
      row = row + 1;
    }
  }

  for (int i=0; i<nDim; i++){
    prevAvals[i] = aprev[i][0];
  }

  // if you have defined any other value of tEnd than zero in podDict, this will get skipped.
  if(tEnd == 0){
    tEnd = caseRunT;
  }
  double nSteps = (tEnd/dt);  // total time steps to loop through
  double timeElapsed = 0.0;   // time counter

  // 2D vector for storing a values at every time step  
  std::vector<std::vector<double>> avalues(nSteps+1, std::vector<double>(nDim,0.0));

  // 1D vector to store all a's in one time iteration
  std::vector<double> avals(nDim,0.0);

  double da = 0.0; // Initializing da

  for (int t=0; t<nSteps; t++){
  
    for (int k=0; k<nDim; k++){
      da = constant[k];
     for (int m=0; m<nDim; m++){
        da += linear[k][m]*prevAvals[m];
        for (int n=0; n<nDim; n++){
           da += quadratic[k][m][n]*prevAvals[n]*prevAvals[m]; 
        }
     }
     avals[k] = prevAvals[k] + da*dt; // updating current a values
    }

    for (int i=0; i<nDim; i++){
      prevAvals[i] = avals[i]; // updating previous avals
      avalues[t][i] = avals[i]; // storing a values for in 2D vector for later
    }

    timeElapsed += dt; 
    cout << "t = " << timeElapsed << endl; // Case progress info in terminal
  }


  // Now writing a values form 2D vector to CSV file for use in..
  // "podFlowReconstruct"
  int writeSteps = 0;

  if(writeFreq == 0){
    writeSteps = nSteps/numDirs;}
  else{int writeSteps = writeFreq;}

  ofstream afiles;
  afiles.open("avals.csv");

  for (int i=0; i<nSteps; i++){
    if(i%writeSteps == 0){
      //double tcol = ceil(dt*i*100)/100;
      double tcol = dt*i;
      afiles << tcol << ",";
      for (int j=0; j<nDim; j++){
        afiles << avalues[i][j] << ",";
      }

      afiles << endl << std::flush;
    }
  }
  afiles.close();

  // processor clock time info displays when program ends
  duration = (std::clock() - start ) / (double) CLOCKS_PER_SEC;

  cout << "runtime = " << duration << " seconds" << endl;

  return 0;
}


// ************************************************************************* //
