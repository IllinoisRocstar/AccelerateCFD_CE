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
  podBasisCalc

Description
  Performs proper orthogonal decomposition on the velocity field. Before
  running make sure to leave only non-correlated time step folders.
  Initial version of an application to perform orthogonal decomposition 
  on the velocity field. Must be run on the root directory and will use
  all time steps to perform calculations. At the end the eigenvalues,
  (normalized) eigenvectors and POD modes of the solution will be calculated.
  Requires UMean file in last time step folder.

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
#include <Eigen/Dense>
#include <vector>

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

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
      IOobject::AUTO_WRITE
    );
}

int main(int argc, char *argv[])
{

  copyrightnotice();

  // To get processor clocktime
  std::clock_t start;
  double duration;
  start = std::clock();

  argList::validArgs.append("basisToWrite");
  argList::addOption
  (
    "basisToWrite",
    "amount",
    "Specify number of basis to write"
  );

  timeSelector::addOptions();

  // Add a function to only select user defined times as snapshots

  Foam::argList args(argc, argv);
  if (!args.checkRootCase())
  {
    Foam::FatalError.exit();
  }

  int numBasis = 0;

  if (!args.check())
    numBasis = 0;
  else
    numBasis = std::stoi(args[1]);  
  
  #include "createTime.H" 
  #include "createNamedMesh.H"  
  instantList timeDirs = timeSelector::select0(runTime, args);

  // defining correlation matrix Cmn
  int nDim(timeDirs.size());

  if (numBasis > nDim) {
    std::cerr << "Please select number of basis between 1 and " << nDim << std::endl;
    throw;
  }

  Eigen::MatrixXd Cmn(nDim, nDim);
   
  int m,n;
  m = 0;
    
  runTime.setTime(timeDirs.last(),0); 

  //Let's start collecting some data for calculations
  //Reading cell volumes from mesh
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
    
  // Reading mean velocity for case from last time step. This will be used for calculating..
  // basis as well as for reduced order model
  volVectorField UMean = generateMeshField(runTime,mesh,"UMean");

  // Reading and storing all velocities from every time directories into a vector.
  std::vector<volVectorField> vels;
  Info<< "Reading fields U" << nl;
    
  forAll(timeDirs, timei)
  {
    runTime.setTime(timeDirs[timei], timei);
    volVectorField U = generateMeshField(runTime,mesh,"U");   
    vels.push_back(volVectorField(U-UMean)); // Extracting mean velocity from the flow to..
                             // get velocity fluctuations. POD basis will..
                             // represent these fluctuations in velocities
  }

  // Now that we have collected all the data, let's start calculations.

  // Correlation matrix (Cmn) is used to calculate eigenvalues and eigenvectors associated..
  // ..with fluctuations in velocities. Later these eigenvectors will be used to calculate POD basis 
  Info<< "Assembling matrix Cmn" << nl;
       
  forAll(timeDirs, timei)
  {
    n = 0;
    forAll(timeDirs, timej)
    {
       Cmn(m, n) = 0.0;
      // applying symmetry 
      if (n < m)
      {
        Cmn(m, n) = Cmn(n, m);
        n++;
        continue;
      }

      volScalarField U1dotU2(generateCustomField(runTime,mesh,"U1dotU2"),
                             (vels[timei]&vels[timej])*cellVolume);

      Cmn(m,n) = Cmn(m,n) + gSum(U1dotU2);
      n++;
    }
    m++;
  }
    
  // Normalization of correlation matrix by dividing with total number of velocities used (nDim)
  Cmn = Cmn/nDim;
    
  // Self Adjoint Eigen Solver is used here to solve for eigenvalue problem using Eigen C++ library.
  Info<< "Solving eigenvalue problem" << nl;
     
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(Cmn);

  if (es.info()!=Eigen::Success)
  {
    Info << "Eigen value calculations failed" << endl;
    return(-1);
  }

  // Calculating energy contained in each POD basis and writing to csv file, knowing the energy..
  // contained in basis helps us decide how many basis to use for reduced order model. 
  Eigen::VectorXd eigVal(nDim);
  Eigen::VectorXd indEnergy(nDim);
  Eigen::VectorXd totalEnergy(nDim);
  Eigen::VectorXd podNum(nDim);
  totalEnergy[0] = 0.0;
  eigVal = es.eigenvalues().reverse();
  double sumeig = 0.0;
  double sum = 0.0;

  for (int i=0; i<nDim; i++)
  {
    sumeig = sumeig + eigVal[i];
    podNum[i] = i+1;
  }

  for (int i=0; i<nDim; i++)
  {
    indEnergy[i] = (eigVal[i]/sumeig)*100;
    sum = sum + indEnergy[i];
    totalEnergy[i] = sum;
  }

  // Writing energy contained in basis to CSV file
  ofstream myfile;
  myfile.open ("podEnergy.csv");
  myfile << "Basis#,Individual_Energy_in_Basis(%),Cummulative_Energy_in_Basis_upto_Current_Basis(%), Eigen_Values" << nl;
  for (int i=0; i<nDim; i++)
    myfile << podNum[i] << "," << indEnergy[i] << "," << totalEnergy[i] << ", " << eigVal[i] << nl;

  myfile.close();

  // Storing eigenvectors and eigenvalues
  scalarField eigenVal(nDim);
  Eigen::VectorXd::Map(&eigenVal[0], eigenVal.size()) = es.eigenvalues().reverse();

  List<scalar> eigenVec(nDim);
  List<List<scalar>> eigenVecList(nDim,eigenVec);

  forAll(eigenVecList,i) 
  {
    List<scalar>& vList = eigenVecList[nDim-i-1];
    forAll(vList,j)
    {
      scalar& s = vList[j];
      s = es.eigenvectors()(j,i);
    }
  }

  // normalize eigenvectors
  forAll(timeDirs, i)
  {
    scalar norm = eigenVal[i]*nDim;
    eigenVecList[i] = eigenVecList[i]*std::sqrt(norm);
  }

  // Calculation of POD basis using eigenvectors and velocities
  Info << "Saving pod basis in " << runTime.timeName() << endl;

  if (numBasis == 0)
    numBasis = nDim;

  if (numBasis == nDim)
    numBasis = nDim;
  
  for (int iSig=0; iSig<numBasis; iSig++)
  {
    std::string sigmaName;
    sigmaName = "sigma_" + std::to_string(iSig);

    volVectorField sigma(generateCustomField(runTime,mesh,sigmaName),mesh,
                         dimensionedVector("0",dimLength/dimTime,Zero));
       
    // Constructing POD basis from eigen vectors and velocities
    forAll(timeDirs, timei)
    {          
      sigma = sigma + eigenVecList[iSig][timei]*vels[timei];
    }
       
    sigma = sigma/(nDim*eigenVal[iSig]); //normalizing modes as per as per Grau (2007) eq. 6

    //POD modes written to sigma_0, sigma_1, etc in last time directory of case.
    sigma.write();
  }

  // processor clock time info displays when program ends
  duration = (std::clock() - start ) / (double) CLOCKS_PER_SEC;

  Info << "runtime = " << duration << " seconds" << endl << nl;
    
  return 0;
}


// ************************************************************************* //

