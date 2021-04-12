cmake_minimum_required(VERSION 3.1)

# Set project metadata #########################################################
set(PROJECT_NAME "AccelerateCFD_CE")
project(${PROJECT_NAME} CXX)

SET(DEFAULT_INSTALL_PREFIX "$ENV{FOAM_USER_APPBIN}/../" CACHE PATH "Prefix prepended to
install directories")
SET(CMAKE_INSTALL_PREFIX "${DEFAULT_INSTALL_PREFIX}" CACHE INTERNAL "Prefix
prepended to install directories" FORCE)

set(ACFD_CE_MAJOR_VERSION 0)
set(ACFD_CE_MINOR_VERSION 3)
set(ACFD_CE_PATCH_VERSION 0)
set(ACFD_CE_EXTRA_VERSION "" CACHE STRING "AccelerateCFD_CE extra version string")
set(ACFD_CE_VERSION "${ACFD_CE_MAJOR_VERSION}.${ACFD_CE_MINOR_VERSION}")
set(ACFD_CE_VERSION
    "${ACFD_CE_VERSION}.${ACFD_CE_PATCH_VERSION}${ACFD_CE_EXTRA_VERSION}")
set(ACFD_CE_SHORT_LICENSE "GNU General Public License v3 or later (GPLv3+)")
set(ACFD_CE_PUBLIC_URL "https://github.com/IllinoisRocstar/AccelerateCFD_CE")
set(ACFD_CE_AUTHOR
    "Taiyo Wilson (lead)"
    "Timothy Grunloh"
    "Masoud Safdari"
    "Akash Patel"
    "Luke Calian"
    "Pavle Simonovic")
set(ACFD_CE_AUTHOR_EMAIL "tpgrunloh@illinoisrocstar.com")
set(ACFD_CE_DESCRIPTION "AccelerateCFD uses Proper Orthogonal Decomposition (POD), \
principle to reduce the full order CFD flow into a reduced order model which runs \
several magnitudes faster and reconstructs the flow fields.")
set(PROJECT_NAME "AccelerateCFD_CE ${ACFD_CE_VERSION}")
project(${PROJECT_NAME} CXX C)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_EXTENSIONS OFF)

add_definitions(
-std=c++11
-Dlinux64
-m64
-DWM_ARCH_OPTION=64
-DWM_DP
-DWM_LABEL_SIZE=32
-Wall
-Wextra
-Wold-style-cast
-Wnon-virtual-dtor
-Wno-unused-parameter
-Wno-invalid-offsetof
-O3
-DNoRepository
-ftemplate-depth-100
-fpermissive
-DNoRepository
-fPIC
)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_SOURCE_DIR}/cmake/Modules/)

set(DEFAULT ON CACHE INTERNAL "Default value for enabled-by-default options")
macro(opt OPTION HELP VALUE)
  option(ENABLE_${OPTION} ${HELP} ${VALUE})
  set(OPT_TEXI "${OPT_TEXI}\n@item ENABLE_${OPTION}\n${HELP} (default: ${VALUE})")
endmacro(opt)

# set options
opt(TESTING "Enable testing" ${DEFAULT})

set( CMAKE_CXX_FLAGS "-Wall -Wextra -Wold-style-cast -Wnon-virtual-dtor -Wno-unused-parameter -Wno-invalid-offsetof -ftemplate-depth-100 -DOMPI_SKIP_MPICXX -O3 -fPIC" )
set( CMAKE_MODULE_LINKER_FLAGS "-Xlinker --copy-dt-needed-entries -Xlinker --no-as-needed" )
set( CMAKE_SHARED_LINKER_FLAGS "-Xlinker --copy-dt-needed-entries -Xlinker --no-as-needed" )
set( CMAKE_EXE_LINKER_FLAGS "-Xlinker --copy-dt-needed-entries -Xlinker --no-as-needed" )

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/externalLibraries)
include_directories(
  $ENV{FOAM_SRC}/OSspecific/POSIX/lnInclude
  $ENV{FOAM_SRC}/OpenFOAM/lnInclude
  $ENV{FOAM_SRC}/triSurface/lnInclude
  $ENV{FOAM_SRC}/parallel/decompose/decompositionMethods/lnInclude
  $ENV{FOAM_SRC}/mesh/snappyHexMesh/lnInclude
  $ENV{FOAM_SRC}/fileFormats/lnInclude
  $ENV{FOAM_SRC}/surfMesh/lnInclude
  $ENV{FOAM_SRC}/dynamicMesh/lnInclude
  $ENV{FOAM_SRC}/lagrangian/basics/lnInclude
  $ENV{FOAM_APP}/utilities/postProcessing/dataConversion/foamToVTK/foamToVTK/lnInclude
  $ENV{FOAM_SRC}/edgeMesh/lnInclude
  $ENV{FOAM_SRC}/mesh/blockMesh/lnInclude
  $ENV{FOAM_SRC}/transportModels/incompressible/singlePhaseTransportModel
  $ENV{FOAM_SRC}/TurbulenceModels/turbulenceModels/lnInclude
  $ENV{FOAM_SRC}/TurbulenceModels/incompressible/lnInclude
  $ENV{FOAM_SRC}/TurbulenceModels/compressible/lnInclude
  $ENV{FOAM_SRC}/transportModels/compressible/lnInclude
  $ENV{FOAM_SRC}/thermophysicalModels/basic/lnInclude
  $ENV{FOAM_SRC}/finiteVolume/lnInclude
  $ENV{FOAM_SRC}/fvOptions/lnInclude
  $ENV{FOAM_SRC}/meshTools/lnInclude
  $ENV{FOAM_SRC}/sampling/lnInclude
)

set(ACFD_CE_SRC
  utilities/podBasisCalc.C
  utilities/podPrecompute.C
  utilities/podROM.C
  utilities/podFlowReconstruct.C
  utilities/podPostProcess.C
)

add_library(AccelerateCFD_CE ${ACFD_CE_SRC})

target_include_directories(AccelerateCFD_CE
  PUBLIC ${OPNF_INCLUDE_DIRS})

target_link_libraries(AccelerateCFD_CE
  PUBLIC
  $ENV{FOAM_LIBBIN}/libOpenFOAM.so
  $ENV{FOAM_LIBBIN}/libturbulenceModels.so
  $ENV{FOAM_LIBBIN}/libincompressibleTurbulenceModels.so
  $ENV{FOAM_LIBBIN}/libcompressibleTurbulenceModels.so
  $ENV{FOAM_LIBBIN}/libincompressibleTransportModels.so
  $ENV{FOAM_LIBBIN}/libcompressibleTransportModels.so
  $ENV{FOAM_LIBBIN}/libfluidThermophysicalModels.so
  $ENV{FOAM_LIBBIN}/libspecie.so
  $ENV{FOAM_LIBBIN}/libfiniteVolume.so
  $ENV{FOAM_LIBBIN}/libfvOptions.so
  $ENV{FOAM_LIBBIN}/libgenericPatchFields.so
  $ENV{FOAM_LIBBIN}/libmeshTools.so
  $ENV{FOAM_LIBBIN}/libsampling.so
  $ENV{FOAM_LIBBIN}/libdistributed.so
)

add_executable(podBasisCalc utilities/podBasisCalc.C)
add_executable(podPrecompute utilities/podPrecompute.C)
add_executable(podROM utilities/podROM.C)
add_executable(podFlowReconstruct utilities/podFlowReconstruct.C)
add_executable(podPostProcess utilities/podPostProcess.C)

target_link_libraries(podBasisCalc AccelerateCFD_CE)
target_link_libraries(podPrecompute AccelerateCFD_CE)
target_link_libraries(podROM AccelerateCFD_CE)
target_link_libraries(podFlowReconstruct AccelerateCFD_CE)
target_link_libraries(podPostProcess AccelerateCFD_CE)

install(TARGETS podBasisCalc DESTINATION bin)
install(TARGETS podPrecompute DESTINATION bin)
install(TARGETS podROM DESTINATION bin)
install(TARGETS podFlowReconstruct DESTINATION bin)
install(TARGETS podPostProcess DESTINATION bin)

