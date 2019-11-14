/*---------------------------------------------------------------------------*\
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
  podFlowReconstruct

Description
  Reconstructs velocity using basis and time coefficients and writes it in every directory.

Author
  Illinois Rocstar
  AccelerateCFD Development Team
  Copyright (C) 2017-2019
    
References
  Zhu Wang et al. (2011)
  Proper Orthogonal Decomposition Closure Models For Turbulent Flows: A Numerical Comparison

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "timeSelector.H"
#include "ReadFields.H"
#include "volFields.H"
#include "IFstream.H"
#include "OFstream.H"
#include "fvc.H"
#include <math.h> 
#include <sstream>
#include <iostream>
#include <vector>
using namespace Foam;
#define PI 3.14159265
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Innnerproduct of 2 vectors
dimensionedScalar innerProductPOD(volVectorField v1, volVectorField v2, volScalarField cellVols) // inner product of 2 vectors
{
  return gSum(volScalarField((v1&v2)*cellVols));
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

  #include "setRootCase.H"       
  #include "createTime.H"
  #include "createNamedMesh.H"

  instantList timeDirs = timeSelector::select0(runTime, args);
  runTime.setTime(0.0,0); 

  // Reading data from podDict
  IOdictionary podDict
  (
    IOobject
    (
      "podDict",
      runTime.constant(),
      mesh,
      IOobject::MUST_READ,
      IOobject::NO_WRITE
    )
  );

  dimensionedScalar nDim_ds
  (
  podDict.lookup("nDim")
  );

  int nDim = (int) nDim_ds.value(); //# of modes used
    
  runTime.setTime(timeDirs.last(),0);

  volVectorField UMean // read mean velocity field from last timestep
  (
    IOobject
    (
      "UMean",
      runTime.timeName(),
      mesh,
      IOobject::MUST_READ,
      IOobject::NO_WRITE
    ),
    mesh
  );

  std::vector<volVectorField> sigs; // vector for storing velocity modes

  for (int iSig=0; iSig<nDim; iSig++) // reading POD basis from last case directory
  {
    std::string sigmaName;
    sigmaName = "sigma_" + std::to_string(iSig);
        
    sigs.push_back
    (
      volVectorField
      (
        IOobject
        (
          sigmaName,
          runTime.timeName(),
          mesh,
          IOobject::MUST_READ,
          IOobject::NO_WRITE
        ),
      mesh
      )
    );
  }

  volVectorField Uinit // reading initial velocity
  (
    IOobject
    (
      "U",
      "0",
      mesh,
      IOobject::MUST_READ,
      IOobject::NO_WRITE
    ),
    mesh
  );
    
  //Reading a values and composing 2D vector
  string line;
  ifstream afile ("avals.csv");
  std::vector<double> time;
  std::vector<std::vector<double>> aVals;

  double dval = 0.0;
    
  while(getline (afile,line)){
    std::stringstream ss(line);
    std::vector<double> vect;
    int lineInd = 0;
    while(ss >> dval){
      if (lineInd == 0){
        time.push_back(dval);
      }else{
        vect.push_back(dval);
      }
      lineInd++;
      if (ss.peek() == ',')
        ss.ignore();
      }
    aVals.push_back(vect);
  }

  int i=0;

  // this loops through all time directories in case and writes reconstructed..
  // velocity (Urom)

  forAll(timeDirs,timei)
  {    
    std::vector<double> aVect = aVals[i];
    Info << "t = " << time[i] << nl; // Case progres info in terminal

    std::string UName;
    UName = "Urom";
    volVectorField Urom
    (
      IOobject
      (
        UName,
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
      ),
      UMean
    );
          
    Urom = dimensionedVector("0", dimLength/dimTime, Zero);
       
    Urom = UMean;

    for (int i=0; i<nDim; i++){
      Urom += aVect[i]*sigs[i];
    }
    
    runTime.setTime(timeDirs[timei],0);

    Urom.write(); // writes Urom in every time directory 
    i++;
  }

  // processor clock time info displays when program ends
  duration = (std::clock() - start ) / (double) CLOCKS_PER_SEC;

  Info << "runtime = " << duration << " seconds" << nl;
  
  return 0;
}


// ************************************************************************* //
