# FindOpenFOAM
# ------------
#
# Find OpenFOAM, a CFD library with several solvers and simulation service modules.
#
#
# To provide the module with a hint about where to find your OpenFOAM
# installation, you can set the environment variable by sourcing the proper
# OpenFOAM environment setup scripts. The find module will then look in a few
# environment variables when searching for OpenFOAM executables, paths, and
# libraries.
#
# TODO:
# In addition to finding the headers and libraries required to compile an
# OpenFOAM-based application, this module also makes an effort to find tools
# that come with the distribution.
#
# This module will define the following variables:
#
# ::
#
#   OPNF_INCLUDE_DIRS - location of the OpenFOAM headers
#   OPNF_LIBRARIES - required libraries for all requested bindings
#   OPNF_FOUND - true if OpenFOAM was found on the system
#   OPNF_VERSION - OpenFOAM version in format Major.Minor.Release
#   OPNF_LIBRARY_DIRS - the full set of library directories
#   OPNF_COMPILE_DEFINITIONS - compile definition needed by OpenFOAM
#
#===============================================================================
# Copyright 2019 Masoud Safdari
#
# This software is distributed WITHOUT ANY WARRANTY; without even the
# implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
# See the License for more information.
#===============================================================================

include(${CMAKE_ROOT}/Modules/SelectLibraryConfigurations.cmake)
include(${CMAKE_ROOT}/Modules/FindPackageHandleStandardArgs.cmake)

# List of minimum required OpenFOAM library components
set(OPNF_COMPONENT_BINDINGS)
set(OPNF_VALID_COMPONENTS
    barotropicCompressibilityModel
    blockMesh
    chemistryModel
    coalCombustion
    combustionModels
    compressibleEulerianInterfacialModels
    compressibleMultiphaseEulerianInterfacialModels
    compressibleTransportModels
    compressibleTurbulenceModels
    compressibleTwoPhaseSystem
    conversion
    decompose
    decompositionMethods
    distributed
    distributionModels
    DPMTurbulenceModels
    driftFluxRelativeVelocityModels
    driftFluxTransportModels
    DSMC
    dynamicFvMesh
    dynamicMesh
    edgeMesh
    engine
    extrude2DMesh
    extrudeModel
    fieldFunctionObjects
    fileFormats
    finiteVolume
    fluidThermophysicalModels
    foamToVTK
    forces
    fvMotionSolvers
    fvOptions
    genericPatchFields
    helpTypes
    immiscibleIncompressibleTwoPhaseMixture
    incompressibleTransportModels
    incompressibleTurbulenceModels
    interfaceProperties
    lagrangianFunctionObjects
    lagrangianIntermediate
    lagrangian
    lagrangianSpray
    lagrangianTurbulence
    laminarFlameSpeedModels
    liquidMixtureProperties
    liquidProperties
    meshTools
    molecularMeasurements
    molecule
    multiphaseInterFoam
    multiphaseMixtureThermo
    multiphaseReactingTurbulenceModels
    multiphaseSystem
    ODE
    OpenFOAM
    pairPatchAgglomeration
    phaseChangeTwoPhaseMixtures
    phaseCompressibleTurbulenceModels
    potential
    Pstream
    pyrolysisModels
    radiationModels
    randomProcesses
    reactingEulerianInterfacialCompositionModels
    reactingEulerianInterfacialModels
    reactingMultiphaseSystem
    reactingPhaseSystem
    reactingTwoPhaseSystem
    reactionThermophysicalModels
    reconstruct
    regionCoupled
    regionCoupling
    regionModels
    renumberMethods
    rhoCentralFoam
    rigidBodyDynamics
    rigidBodyMeshMotion
    sampling
    scotchDecomp
    sixDoFRigidBodyMotion
    SLGThermo
    SloanRenumber
    snappyHexMesh
    solidChemistryModel
    solidMixtureProperties
    solidParticle
    solidProperties
    solidSpecie
    solidThermo
    solverFunctionObjects
    specie
    surfaceFilmDerivedFvPatchFields
    surfaceFilmModels
    surfMesh
    tabulatedWallFunctions
    thermalBaffleModels
    thermophysicalFunctions
    topoChangerFvMesh
    triSurface
    turbulenceModels
    twoPhaseMixture
    twoPhaseMixtureThermo
    twoPhaseProperties
    twoPhaseReactingTurbulenceModels
    userd-foam
    utilityFunctionObjects
    transportModels
)

# set initial parameters from OpenFOAM environment variables
# check existance and read variables, preps for search step
if(EXISTS $ENV{FOAM_SRC})
  # search standard environment variables
  set(OPNF_VERSION $ENV{WM_PROJECT_VERSION})
  set(OPNF_INST_DIR $ENV{WM_PROJECT_DIR})
  set(OPNF_SRC_DIR $ENV{FOAM_SRC})
  set(OPNF_LIB_DIR $ENV{FOAM_LIBBIN})
  set(OPNF_APP_DIR $ENV{FOAM_APPBIN})
  set(OPNF_SITE_APP_DIR $ENV{FOAM_SITE_APPBIN})
  set(OPNF_SITE_LIB_DIR $ENV{FOAM_SITE_LIBBIN})
else()
  # search if root directory is set
  message(FATAL_ERROR "Unable to determine openfoam installation location.")
endif()

# if no component is defined all OpenFOAM components will
# be added to the list otherwise, validate the list of find components
if(NOT OpenFOAM_FIND_COMPONENTS)
  set(OPNF_COMPONENT_BINDINGS ${OPNF_VALID_COMPONENTS})
else()
  # add the extra specified components, ensuring that they are valid.
  foreach(component ${OpenFOAM_FIND_COMPONENTS})
    list(FIND OPNF_VALID_COMPONENTS ${component} component_location)
    if(${component_location} EQUAL -1)
      message(FATAL_ERROR "\"${component}\" is not a valid OpenFOAM component.")
    else()
      find_library(${component}_lib ${component} PATH ${OPNF_LIB_DIR})
      list(APPEND OPNF_COMPONENT_BINDINGS ${${component}_lib})
    endif()
  endforeach()
endif()

# OpenFOAM header files are distributed in several folders, making it hard to
# add their paths into a single variable (will cause long arg list problems
# during compilation step). The following lines only add headers needed for the
# compilation of the current project's OpenFOAM dependent modules. If used for
# outside this project these lines should be changed to include the
# dependenices needed.
set(OPNF_INC_DIR
    ${OPNF_INST_DIR}/src/OSspecific/POSIX/lnInclude
    ${OPNF_INST_DIR}/src/OpenFOAM/lnInclude
    ${OPNF_INST_DIR}/src/finiteVolume/lnInclude
    ${OPNF_INST_DIR}/src/meshTools/lnInclude
    ${OPNF_INST_DIR}/src/triSurface/lnInclude
    ${OPNF_INST_DIR}/src/parallel/decompose/decompositionMethods/lnInclude
    ${OPNF_INST_DIR}/src/mesh/snappyHexMesh/lnInclude
    ${OPNF_INST_DIR}/src/fileFormats/lnInclude
    ${OPNF_INST_DIR}/src/surfMesh/lnInclude
    ${OPNF_INST_DIR}/src/dynamicMesh/lnInclude
    ${OPNF_INST_DIR}/src/lagrangian/basics/lnInclude
    ${OPNF_INST_DIR}/applications/utilities/postProcessing/dataConversion/foamToVTK/foamToVTK/lnInclude
    ${OPNF_INST_DIR}/src/edgeMesh/lnInclude
    ${OPNF_INST_DIR}/src/mesh/blockMesh/lnInclude
    ${OPNF_INST_DIR}/src/transportModels/incompressible/singlePhaseTransportModel
    ${OPNF_INST_DIR}/src/TurbulenceModels/turbulenceModels/lnInclude
    ${OPNF_INST_DIR}/src/TurbulenceModels/incompressible/lnInclude
    ${OPNF_INST_DIR}/src/TurbulenceModels/compressible/lnInclude
    ${OPNF_INST_DIR}/src/transportModels/compressible/lnInclude
    ${OPNF_INST_DIR}/src/thermophysicalModels/basic/lnInclude
    ${OPNF_INST_DIR}/src/finiteVolume/lnInclude
    ${OPNF_INST_DIR}/src/fvOptions/lnInclude
    ${OPNF_INST_DIR}/src/meshTools/lnInclude
    ${OPNF_INST_DIR}/src/sampling/lnInclude
)

# OpenFOAM has custom compiler definitions. All sourced from the OpenFOAM env.
set(OPNF_COMPILE_DEFINITIONS
    $ENV{WM_ARCH}
    WM_ARCH_OPTION=$ENV{WM_ARCH_OPTION}
    WM_$ENV{WM_PRECISION_OPTION}
    WM_LABEL_SIZE=$ENV{WM_LABEL_SIZE}

    NoRepository
)

# setting the output variables
set(OPNF_INCLUDE_DIRS ${OPNF_INC_DIR})
set(OPNF_LIBRARIES ${OPNF_COMPONENT_BINDINGS})
set(OPNF_LIBRARY_DIRS ${OPNF_LIB_DIR})

# We may have picked up some duplicates in various lists during the above
# process for the language bindings (both the C and C++ bindings depend on
# libz for example). Remove the duplicates. It appears that the default
# CMake behavior is to remove duplicates from the end of a list. However,
# for link lines, this is incorrect since unresolved symbols are searched
# for down the link line. Therefore, we reverse the list, remove the
# duplicates, and then reverse it again to get the duplicates removed from
# the beginning.
macro(_remove_duplicates_from_beginning _list_name)
  list(REVERSE ${_list_name})
  list(REMOVE_DUPLICATES ${_list_name})
  list(REVERSE ${_list_name})
endmacro()

if(OPNF_INCLUDE_DIRS)
  _remove_duplicates_from_beginning(OPNF_INCLUDE_DIRS)
endif()
if(OPNF_LIBRARY_DIRS)
  _remove_duplicates_from_beginning(OPNF_LIBRARY_DIRS)
endif()

# final processing
message(STATUS "OpenFOAM library location ${OPNF_LIBRARY_DIRS}")
find_package_handle_standard_args(OpenFOAM
    REQUIRED_VARS OPNF_LIBRARIES OPNF_INCLUDE_DIRS
    VERSION_VAR OPNF_VERSION)