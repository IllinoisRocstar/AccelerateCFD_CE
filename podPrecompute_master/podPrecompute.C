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
  podPrecompute

Description
  Computes and writes all the necessary data for reduced order model.

Author
  Illinois Rocstar LLC
  AccelerateCFD Development Team
  Copyright (C) 2017-2019

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "timeSelector.H"
#include "ReadFields.H"
#include "volFields.H"
#include "IFstream.H"
#include "OFstream.H"
#include "fvc.H"
#include <vector>
#include <iostream>
#include <iomanip>
#include <math.h>
#include <string>
#include <sstream>
#include <fstream>
#include <map>
#include <iterator>

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// inner product of 2 vectors
double innerProductPOD(volVectorField v1, volVectorField v2, volScalarField cellVols)
{
  return gSum(volScalarField((v1&v2)*cellVols));
}

// inner product of 2 tensors (double dot product)
double innerProductPOD2(volTensorField v1, volSymmTensorField v2, volScalarField cellVols)
{
  return gSum(volScalarField((v1&&v2)*cellVols));
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

volVectorField generateMeshField(Foam::Time &runTime, Foam::fvMesh &mesh,
    const std::string &fieldName) {
    
  volVectorField mshField
  (
    IOobject
    (
      fieldName,
      runTime.timeName(),
      mesh,
      IOobject::MUST_READ,
      IOobject::NO_WRITE
    ),
    mesh
  );

  return mshField;
}

IOobject generateCustomField(Foam::Time &runTime, Foam::fvMesh &mesh,
    const std::string &fieldName) {
    
    return IOobject
    (
      fieldName,
      runTime.timeName(),
      mesh,
      IOobject::NO_READ,
      IOobject::NO_WRITE
    );
}

int main(int argc, char *argv[])
{
  copyrightnotice();

  // To get processor clocktime
  std::clock_t start;
  double duration;
  start = std::clock();

  timeSelector::addOptions();

  #include "setRootCase.H"       
  #include "createTime.H"
  #include "createNamedMesh.H"

  instantList timeDirs = timeSelector::select0(runTime, args);

  // read initial velocity
  runTime.setTime(timeDirs[0],0);
  volVectorField U = generateMeshField(runTime,mesh,"U");

  // Reads user-defined data from podDict in constant directory
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

  dimensionedScalar nDim_ds(podDict.lookup("nDim"));
  dimensionedScalar nu_ds(podDict.lookup("nu"));
  dimensionedScalar writeFreq_ds(podDict.lookup("writeFreq"));
  dimensionedScalar tEnd_ds(podDict.lookup("tEnd"));
  dimensionedScalar dt_ds(podDict.lookup("dt"));
  dimensionedScalar closure_nu(podDict.lookup("artificial_nu"));

  int nDim = static_cast<int>(nDim_ds.value());              //# of modes used
  double nu = static_cast<double>(nu_ds.value());            //nu value
  int writeFreq = static_cast<int>(writeFreq_ds.value());    //Closure model update frequency
  double tEnd = static_cast<double>(tEnd_ds.value());        //closure model runtime end
  double dt = static_cast<double>(dt_ds.value());            //FDM time step
  double nu_tilda = static_cast<double>(closure_nu.value()); //closure viscosity

  // Reads important data for reduced order model from controlDict
  IOdictionary controlDict
  (
    IOobject
    (
      "controlDict",
      runTime.system(),
      mesh,
      IOobject::MUST_READ,
      IOobject::NO_WRITE
    )
  );

  double staTime = timeDirs.first().value();
  double endTime = timeDirs.last().value();
  scalar deltaT( readScalar(runTime.controlDict().lookup("deltaT")));
  scalar writeInterval( readScalar(runTime.controlDict().lookup("writeInterval")));
  word writeControl( word(runTime.controlDict().lookup("writeControl")));

  if (tEnd == 0)
    tEnd = endTime;

  double tSteps = (endTime - staTime)/deltaT;
  double runCase = (endTime - staTime);
  double numDirs;
  if (writeControl == "runTime" || writeControl == "adjustableRunTime")
    numDirs = (endTime - staTime)/writeInterval;
  else if (writeControl == "timeStep")
    numDirs = tSteps/writeInterval;

  // Reading cellVolumes from mesh
  runTime.setTime(timeDirs.last(),0);

  volScalarField cellVolume
  (
    IOobject
    (
      "cellVolume",
      runTime.timeName(),
      mesh,
      IOobject::NO_READ,
      IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("zero",dimVolume,0.0)
  );
  cellVolume.ref() = mesh.V();

  int nRows = cellVolume.size(); // Total number of cells

  // Reading mean velocity from last time step
  volVectorField UMean = generateMeshField(runTime,mesh,"UMean");

  std::vector<volVectorField> sigs; // vector for storing POD basis
  std::vector<volTensorField> gradSigs; // tensor for storing gradient of POD basis
  std::vector<volVectorField> laplSigs; // vector for laplacian of POD basis

  // Reads all POD basis from last case directory
  for (int iSig=0; iSig<nDim; iSig++)
  {
    std::string sigmaName;
    sigmaName = "sigma_" + std::to_string(iSig);
    sigs.push_back(generateMeshField(runTime,mesh,sigmaName));

    volVectorField tmp1 = generateMeshField(runTime,mesh,sigmaName);

    volTensorField tmp2(fvc::grad(tmp1));

    std::string laplaName;
    laplaName = "Lapla_" + std::to_string(iSig);
    volVectorField tmp3(generateCustomField(runTime,mesh,laplaName),
                        fvc::laplacian(tmp1));
  
    gradSigs.push_back(tmp2);
    laplSigs.push_back(tmp3);
  } 


  //velocity fluctuations
  volVectorField UPrime(generateCustomField(runTime,mesh,"UPrime"),U-UMean);

  //Calculations of gradients required for reduced order model (ROM)
  volTensorField gradU(fvc::grad(UMean));

  volVectorField UgradU(generateCustomField(runTime,mesh,"UgradU"),UMean&gradU);

  volVectorField laplUMean(generateCustomField(runTime,mesh,"laplUMean"),
                           fvc::laplacian(UMean));

  // Precompute all volVector/Tensor terms in ROM loop
  std::vector<volVectorField> sigGradSigs;
  for (int i=0; i<nDim*nDim; i++) {
    std::string sigGradSigsName;
    sigGradSigsName = "sigGradSigsName_" + std::to_string(i);

    volVectorField tmpName(generateCustomField(runTime,mesh,sigGradSigsName),
                           mesh,
                           dimensionedVector("zero",dimensionSet(0,1,-2,0,0,0,0),Zero));

    sigGradSigs.push_back(tmpName);
  }

  std::vector<volVectorField> uGradSigs;
  std::vector<volVectorField> sigGradUs;
  for (int i=0; i<nDim; i++) {
    std::string uGradSigName;
    uGradSigName = "uGradSigName_" + std::to_string(i);
    volVectorField uGradSig(generateCustomField(runTime,mesh,uGradSigName),
                            UMean&gradSigs[i]);

    uGradSigs.push_back(uGradSig);

    std::string sigGradUName;
    sigGradUName = "sigGradUName_" + std::to_string(i);
    volVectorField sigGradU(generateCustomField(runTime,mesh,sigGradUName),
                            sigs[i]&gradU);
    sigGradUs.push_back(sigGradU);
  }

  for (int i=0; i<nDim; i++) {
    for (int j=0; j<nDim; j++) {
      std::string sigGradSigName;
      sigGradSigName = "sigGradSigName_" + std::to_string(i+nDim*j);
      volVectorField sigGradSig(generateCustomField(runTime,mesh,sigGradSigName),
                                sigs[i]&gradSigs[j]);
      sigGradSigs[i+nDim*j] = sigGradSig;
    }
  }

  std::vector<double> constant(nDim,0.0);
  std::vector<std::vector<double>> linear(nDim, std::vector<double>(nDim,0.0));
  std::vector<std::vector<std::vector<double>>> quadratic(nDim, std::vector<std::vector<double>>(nDim, std::vector<double>(nDim,0.0)));

  // this loop calculates the Galerkin System matrices Q L C for the ROM equation. 
  // constant term, linear term, and quadratic term
  for (int k=0; k<nDim; k++) {
    constant[k] = -1*innerProductPOD(sigs[k],UgradU,cellVolume)+ (nu+nu_tilda)*innerProductPOD(sigs[k],laplUMean,cellVolume);
    for (int m=0; m<nDim; m++) {
      linear[k][m] = -1*innerProductPOD(sigs[k],uGradSigs[m],cellVolume)
      - innerProductPOD(sigs[k],sigGradUs[m],cellVolume) + (nu+nu_tilda)*innerProductPOD(sigs[k],laplSigs[m],cellVolume);
      for (int n=0; n<nDim; n++) {
        quadratic[k][m][n] = -1*innerProductPOD(sigs[k],sigGradSigs[m+nDim*n],cellVolume);
      }
    }
  }

  // Outputs data read from podDict and controlDict for reduced order model (ROM)
  ofstream myfile;
  myfile.open ("podInfo.csv");
  myfile << nDim << nl;
  myfile << nu << nl;
  myfile << writeFreq << nl;
  myfile << tEnd << nl;
  myfile << dt << nl;
  myfile << nRows << nl;
  myfile << runCase << nl;
  myfile << numDirs << nl;
  myfile << staTime << nl;
  myfile.close();

  // calculating initial time coefficients a from initial velocity fluctuation field
  std::vector<double> avalsPrev(nDim,0.0);
  for (int i=0; i<nDim; i++){
    avalsPrev[i] = innerProductPOD(sigs[i],UPrime,cellVolume);
  }

  // All necessary data is calculated. Now writing everything in CSV files so that it can..
  // be read by podROM.cpp program to calculate time varying coefficients of ROM.

  ofstream oldA;
  oldA.open("prevVals.csv");
  ofstream con;
  con.open("constant.csv");
  ofstream lin;
  lin.open("linear.csv");
  ofstream quad;
  quad.open("quadratic.csv");

  // writing avalues
  for (int i=0; i<nDim; i++){
    oldA << std::fixed << std::setprecision(16) << "" << avalsPrev[i] << nl;
  }

  for (int i=0; i<nDim; i++){
    con << std::fixed << std::setprecision(16) << i << "," << constant[i] << nl;
  }

  for (int i=0; i<nDim; i++){
    for (int j=0; j<nDim;j++){
      lin << std::fixed << std::setprecision(16) << i + nDim*j << "," << linear[i][j] << nl;
    }
  }

  for (int i=0; i<nDim; i++){
    for (int j=0; j<nDim; j++){
      for (int k=0; k<nDim;k++){
        quad << std::fixed << std::setprecision(16) << i+j*nDim+k*nDim*nDim << "," << quadratic[i][j][k] << nl;
      }
    }
  }

  oldA.close();
  con.close();
  lin.close();
  quad.close();

  // processor clock time info displays when program ends
  duration = (std::clock() - start ) / (double) CLOCKS_PER_SEC;

  cout << "runtime = " << duration << " seconds" << nl;

  return 0;
}


// ************************************************************************* //
