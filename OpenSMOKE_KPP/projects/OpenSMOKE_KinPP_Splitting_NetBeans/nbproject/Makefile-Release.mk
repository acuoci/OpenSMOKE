#
# Generated Makefile - do not edit!
#
# Edit the Makefile in the project folder instead (../Makefile). Each target
# has a -pre and a -post target defined where you can add customized code.
#
# This makefile implements configuration specific macros and targets.


# Environment
MKDIR=mkdir
CP=cp
GREP=grep
NM=nm
CCADMIN=CCadmin
RANLIB=ranlib
CC=gcc
CCC=mpic++
CXX=mpic++
FC=gfortran
AS=as

# Macros
CND_PLATFORM=GNU-Linux
CND_DLIB_EXT=so
CND_CONF=Release
CND_DISTDIR=dist
CND_BUILDDIR=build

# Include project Makefile
include Makefile

# Object Directory
OBJECTDIR=${CND_BUILDDIR}/${CND_CONF}/${CND_PLATFORM}

# Object Files
OBJECTFILES= \
	${OBJECTDIR}/_ext/56252444/OpenSMOKE_KPP.o \
	${OBJECTDIR}/_ext/56252444/OpenSMOKE_KPP_Definitions.o \
	${OBJECTDIR}/_ext/bac85f60/OpenSMOKE_KPP_BlockMatrixNetwork.o \
	${OBJECTDIR}/_ext/bac85f60/OpenSMOKE_KPP_Communicator.o \
	${OBJECTDIR}/_ext/bac85f60/OpenSMOKE_KPP_ConvectiveNetworkStatistics.o \
	${OBJECTDIR}/_ext/bac85f60/OpenSMOKE_KPP_DataManager.o \
	${OBJECTDIR}/_ext/bac85f60/OpenSMOKE_KPP_Dictionary.o \
	${OBJECTDIR}/_ext/bac85f60/OpenSMOKE_KPP_NewtonMethod_Manager.o \
	${OBJECTDIR}/_ext/bac85f60/OpenSMOKE_KPP_ODE_Manager.o \
	${OBJECTDIR}/_ext/bac85f60/OpenSMOKE_KPP_ReactorNetwork.o \
	${OBJECTDIR}/_ext/bac85f60/OpenSMOKE_KPP_ReactorNetwork_Residuals.o \
	${OBJECTDIR}/_ext/bac85f60/OpenSMOKE_KPP_SingleReactor.o \
	${OBJECTDIR}/_ext/bac85f60/OpenSMOKE_KPP_SingleReactorStatistics.o \
	${OBJECTDIR}/_ext/bac85f60/OpenSMOKE_KPP_SingleReactor_KineticsManager.o \
	${OBJECTDIR}/_ext/bc30b0d1/OpenSMOKE_DirectLinearSolver_Unsymmetric.o \
	${OBJECTDIR}/_ext/bc30b0d1/OpenSMOKE_LAPACK_Dense.o \
	${OBJECTDIR}/_ext/bc30b0d1/OpenSMOKE_LIS_Unsymmetric.o \
	${OBJECTDIR}/_ext/bc30b0d1/OpenSMOKE_PARDISO_Unsymmetric.o


# C Compiler Flags
CFLAGS=

# CC Compiler Flags
CCFLAGS=-fopenmp -m64 -fexpensive-optimizations -funroll-loops -DUSE_MPI -Wall -Wwrite-strings -Wno-strict-aliasing -Wno-unknown-pragmas -Wfatal-errors -g3
CXXFLAGS=-fopenmp -m64 -fexpensive-optimizations -funroll-loops -DUSE_MPI -Wall -Wwrite-strings -Wno-strict-aliasing -Wno-unknown-pragmas -Wfatal-errors -g3

# Fortran Compiler Flags
FFLAGS=

# Assembler Flags
ASFLAGS=

# Link Libraries and Options
LDLIBSOPTIONS=/home/astagni/NumericalLibraries/petsc/petsc-3.7.3/lib/libpetsc.a -lX11 -Wl,--start-group ../../../../NumericalLibraries/intel/mkl/lib/intel64/libmkl_intel_lp64.a ../../../../NumericalLibraries/intel/mkl/lib/intel64/libmkl_sequential.a ../../../../NumericalLibraries/intel/mkl/lib/intel64/libmkl_core.a -Wl,--end-group -lpthread ../../../../NumericalLibraries/intel/mkl/lib/intel64/libmkl_core.a ../../../../NumericalLibraries/intel/mkl/lib/intel64/libmkl_intel_lp64.a ../../../../NumericalLibraries/intel/mkl/lib/intel64/libmkl_sequential.a ../../../../NumericalLibraries/OpenSMOKE/lib/linux/gnu/libOpenSMOKE_Basic_GNU.a ../../../../NumericalLibraries/OpenSMOKE/lib/linux/gnu/libOpenSMOKE_Engine_GNU.a ../../../../NumericalLibraries/OpenSMOKE/lib/linux/gnu/libOpenSMOKE_AddOns_GNU.a ../../../../NumericalLibraries/OpenSMOKE/lib/linux/gnu/libOpenSMOKE_Distributions_GNU.a ../../../../NumericalLibraries/OpenSMOKE/lib/linux/gnu/libOpenSMOKE_SymbolicKinetics_GNU.a ../../../../NumericalLibraries/BzzMath/BzzMath6-dev/lib/release/libBzzMath6_GNU.a ../../../../NumericalLibraries/NumericalRecipes/lib/release/libNumericalRecipes_GNU.a ../../../../NumericalLibraries/lis-1.2.122/lis-nested-complete/lib/liblis.a -lgfortran -lpthread

# Build Targets
.build-conf: ${BUILD_SUBPROJECTS}
	"${MAKE}"  -f nbproject/Makefile-${CND_CONF}.mk ../../exe/OpenSMOKE_KinPP_Splitting_NetBeans.sh

../../exe/OpenSMOKE_KinPP_Splitting_NetBeans.sh: /home/astagni/NumericalLibraries/petsc/petsc-3.7.3/lib/libpetsc.a

../../exe/OpenSMOKE_KinPP_Splitting_NetBeans.sh: ../../../../NumericalLibraries/intel/mkl/lib/intel64/libmkl_core.a

../../exe/OpenSMOKE_KinPP_Splitting_NetBeans.sh: ../../../../NumericalLibraries/intel/mkl/lib/intel64/libmkl_intel_lp64.a

../../exe/OpenSMOKE_KinPP_Splitting_NetBeans.sh: ../../../../NumericalLibraries/intel/mkl/lib/intel64/libmkl_sequential.a

../../exe/OpenSMOKE_KinPP_Splitting_NetBeans.sh: ../../../../NumericalLibraries/OpenSMOKE/lib/linux/gnu/libOpenSMOKE_Basic_GNU.a

../../exe/OpenSMOKE_KinPP_Splitting_NetBeans.sh: ../../../../NumericalLibraries/OpenSMOKE/lib/linux/gnu/libOpenSMOKE_Engine_GNU.a

../../exe/OpenSMOKE_KinPP_Splitting_NetBeans.sh: ../../../../NumericalLibraries/OpenSMOKE/lib/linux/gnu/libOpenSMOKE_AddOns_GNU.a

../../exe/OpenSMOKE_KinPP_Splitting_NetBeans.sh: ../../../../NumericalLibraries/OpenSMOKE/lib/linux/gnu/libOpenSMOKE_Distributions_GNU.a

../../exe/OpenSMOKE_KinPP_Splitting_NetBeans.sh: ../../../../NumericalLibraries/OpenSMOKE/lib/linux/gnu/libOpenSMOKE_SymbolicKinetics_GNU.a

../../exe/OpenSMOKE_KinPP_Splitting_NetBeans.sh: ../../../../NumericalLibraries/BzzMath/BzzMath6-dev/lib/release/libBzzMath6_GNU.a

../../exe/OpenSMOKE_KinPP_Splitting_NetBeans.sh: ../../../../NumericalLibraries/NumericalRecipes/lib/release/libNumericalRecipes_GNU.a

../../exe/OpenSMOKE_KinPP_Splitting_NetBeans.sh: ../../../../NumericalLibraries/lis-1.2.122/lis-nested-complete/lib/liblis.a

../../exe/OpenSMOKE_KinPP_Splitting_NetBeans.sh: ${OBJECTFILES}
	${MKDIR} -p ../../exe
	mpic++ -o ../../exe/OpenSMOKE_KinPP_Splitting_NetBeans.sh ${OBJECTFILES} ${LDLIBSOPTIONS} -O3 -fopenmp

${OBJECTDIR}/_ext/56252444/OpenSMOKE_KPP.o: ../../src/OpenSMOKE_KPP.cpp 
	${MKDIR} -p ${OBJECTDIR}/_ext/56252444
	${RM} "$@.d"
	$(COMPILE.cc) -O2 -I../../../../NumericalLibraries/BzzMath/BzzMath6-dev/hpp -I../../../../NumericalLibraries/OpenSMOKE/hpp -I../../../../NumericalRecipes/hpp -I../../../../MUMPS_4.9.2/libraries/include -I../../../../NumericalLibraries/intel/mkl/include -I../../src -I../../../../NumericalLibraries/lis-1.2.122/lis-nested-complete/include -I/home/astagni/NumericalLibraries/petsc/petsc-3.7.3/include -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/_ext/56252444/OpenSMOKE_KPP.o ../../src/OpenSMOKE_KPP.cpp

${OBJECTDIR}/_ext/56252444/OpenSMOKE_KPP_Definitions.o: ../../src/OpenSMOKE_KPP_Definitions.cpp 
	${MKDIR} -p ${OBJECTDIR}/_ext/56252444
	${RM} "$@.d"
	$(COMPILE.cc) -O2 -I../../../../NumericalLibraries/BzzMath/BzzMath6-dev/hpp -I../../../../NumericalLibraries/OpenSMOKE/hpp -I../../../../NumericalRecipes/hpp -I../../../../MUMPS_4.9.2/libraries/include -I../../../../NumericalLibraries/intel/mkl/include -I../../src -I../../../../NumericalLibraries/lis-1.2.122/lis-nested-complete/include -I/home/astagni/NumericalLibraries/petsc/petsc-3.7.3/include -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/_ext/56252444/OpenSMOKE_KPP_Definitions.o ../../src/OpenSMOKE_KPP_Definitions.cpp

${OBJECTDIR}/_ext/bac85f60/OpenSMOKE_KPP_BlockMatrixNetwork.o: ../../src/kpp/OpenSMOKE_KPP_BlockMatrixNetwork.cpp 
	${MKDIR} -p ${OBJECTDIR}/_ext/bac85f60
	${RM} "$@.d"
	$(COMPILE.cc) -O2 -I../../../../NumericalLibraries/BzzMath/BzzMath6-dev/hpp -I../../../../NumericalLibraries/OpenSMOKE/hpp -I../../../../NumericalRecipes/hpp -I../../../../MUMPS_4.9.2/libraries/include -I../../../../NumericalLibraries/intel/mkl/include -I../../src -I../../../../NumericalLibraries/lis-1.2.122/lis-nested-complete/include -I/home/astagni/NumericalLibraries/petsc/petsc-3.7.3/include -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/_ext/bac85f60/OpenSMOKE_KPP_BlockMatrixNetwork.o ../../src/kpp/OpenSMOKE_KPP_BlockMatrixNetwork.cpp

${OBJECTDIR}/_ext/bac85f60/OpenSMOKE_KPP_Communicator.o: ../../src/kpp/OpenSMOKE_KPP_Communicator.cpp 
	${MKDIR} -p ${OBJECTDIR}/_ext/bac85f60
	${RM} "$@.d"
	$(COMPILE.cc) -O2 -I../../../../NumericalLibraries/BzzMath/BzzMath6-dev/hpp -I../../../../NumericalLibraries/OpenSMOKE/hpp -I../../../../NumericalRecipes/hpp -I../../../../MUMPS_4.9.2/libraries/include -I../../../../NumericalLibraries/intel/mkl/include -I../../src -I../../../../NumericalLibraries/lis-1.2.122/lis-nested-complete/include -I/home/astagni/NumericalLibraries/petsc/petsc-3.7.3/include -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/_ext/bac85f60/OpenSMOKE_KPP_Communicator.o ../../src/kpp/OpenSMOKE_KPP_Communicator.cpp

${OBJECTDIR}/_ext/bac85f60/OpenSMOKE_KPP_ConvectiveNetworkStatistics.o: ../../src/kpp/OpenSMOKE_KPP_ConvectiveNetworkStatistics.cpp 
	${MKDIR} -p ${OBJECTDIR}/_ext/bac85f60
	${RM} "$@.d"
	$(COMPILE.cc) -O2 -I../../../../NumericalLibraries/BzzMath/BzzMath6-dev/hpp -I../../../../NumericalLibraries/OpenSMOKE/hpp -I../../../../NumericalRecipes/hpp -I../../../../MUMPS_4.9.2/libraries/include -I../../../../NumericalLibraries/intel/mkl/include -I../../src -I../../../../NumericalLibraries/lis-1.2.122/lis-nested-complete/include -I/home/astagni/NumericalLibraries/petsc/petsc-3.7.3/include -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/_ext/bac85f60/OpenSMOKE_KPP_ConvectiveNetworkStatistics.o ../../src/kpp/OpenSMOKE_KPP_ConvectiveNetworkStatistics.cpp

${OBJECTDIR}/_ext/bac85f60/OpenSMOKE_KPP_DataManager.o: ../../src/kpp/OpenSMOKE_KPP_DataManager.cpp 
	${MKDIR} -p ${OBJECTDIR}/_ext/bac85f60
	${RM} "$@.d"
	$(COMPILE.cc) -O2 -I../../../../NumericalLibraries/BzzMath/BzzMath6-dev/hpp -I../../../../NumericalLibraries/OpenSMOKE/hpp -I../../../../NumericalRecipes/hpp -I../../../../MUMPS_4.9.2/libraries/include -I../../../../NumericalLibraries/intel/mkl/include -I../../src -I../../../../NumericalLibraries/lis-1.2.122/lis-nested-complete/include -I/home/astagni/NumericalLibraries/petsc/petsc-3.7.3/include -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/_ext/bac85f60/OpenSMOKE_KPP_DataManager.o ../../src/kpp/OpenSMOKE_KPP_DataManager.cpp

${OBJECTDIR}/_ext/bac85f60/OpenSMOKE_KPP_Dictionary.o: ../../src/kpp/OpenSMOKE_KPP_Dictionary.cpp 
	${MKDIR} -p ${OBJECTDIR}/_ext/bac85f60
	${RM} "$@.d"
	$(COMPILE.cc) -O2 -I../../../../NumericalLibraries/BzzMath/BzzMath6-dev/hpp -I../../../../NumericalLibraries/OpenSMOKE/hpp -I../../../../NumericalRecipes/hpp -I../../../../MUMPS_4.9.2/libraries/include -I../../../../NumericalLibraries/intel/mkl/include -I../../src -I../../../../NumericalLibraries/lis-1.2.122/lis-nested-complete/include -I/home/astagni/NumericalLibraries/petsc/petsc-3.7.3/include -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/_ext/bac85f60/OpenSMOKE_KPP_Dictionary.o ../../src/kpp/OpenSMOKE_KPP_Dictionary.cpp

${OBJECTDIR}/_ext/bac85f60/OpenSMOKE_KPP_NewtonMethod_Manager.o: ../../src/kpp/OpenSMOKE_KPP_NewtonMethod_Manager.cpp 
	${MKDIR} -p ${OBJECTDIR}/_ext/bac85f60
	${RM} "$@.d"
	$(COMPILE.cc) -O2 -I../../../../NumericalLibraries/BzzMath/BzzMath6-dev/hpp -I../../../../NumericalLibraries/OpenSMOKE/hpp -I../../../../NumericalRecipes/hpp -I../../../../MUMPS_4.9.2/libraries/include -I../../../../NumericalLibraries/intel/mkl/include -I../../src -I../../../../NumericalLibraries/lis-1.2.122/lis-nested-complete/include -I/home/astagni/NumericalLibraries/petsc/petsc-3.7.3/include -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/_ext/bac85f60/OpenSMOKE_KPP_NewtonMethod_Manager.o ../../src/kpp/OpenSMOKE_KPP_NewtonMethod_Manager.cpp

${OBJECTDIR}/_ext/bac85f60/OpenSMOKE_KPP_ODE_Manager.o: ../../src/kpp/OpenSMOKE_KPP_ODE_Manager.cpp 
	${MKDIR} -p ${OBJECTDIR}/_ext/bac85f60
	${RM} "$@.d"
	$(COMPILE.cc) -O2 -I../../../../NumericalLibraries/BzzMath/BzzMath6-dev/hpp -I../../../../NumericalLibraries/OpenSMOKE/hpp -I../../../../NumericalRecipes/hpp -I../../../../MUMPS_4.9.2/libraries/include -I../../../../NumericalLibraries/intel/mkl/include -I../../src -I../../../../NumericalLibraries/lis-1.2.122/lis-nested-complete/include -I/home/astagni/NumericalLibraries/petsc/petsc-3.7.3/include -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/_ext/bac85f60/OpenSMOKE_KPP_ODE_Manager.o ../../src/kpp/OpenSMOKE_KPP_ODE_Manager.cpp

${OBJECTDIR}/_ext/bac85f60/OpenSMOKE_KPP_ReactorNetwork.o: ../../src/kpp/OpenSMOKE_KPP_ReactorNetwork.cpp 
	${MKDIR} -p ${OBJECTDIR}/_ext/bac85f60
	${RM} "$@.d"
	$(COMPILE.cc) -O2 -I../../../../NumericalLibraries/BzzMath/BzzMath6-dev/hpp -I../../../../NumericalLibraries/OpenSMOKE/hpp -I../../../../NumericalRecipes/hpp -I../../../../MUMPS_4.9.2/libraries/include -I../../../../NumericalLibraries/intel/mkl/include -I../../src -I../../../../NumericalLibraries/lis-1.2.122/lis-nested-complete/include -I/home/astagni/NumericalLibraries/petsc/petsc-3.7.3/include -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/_ext/bac85f60/OpenSMOKE_KPP_ReactorNetwork.o ../../src/kpp/OpenSMOKE_KPP_ReactorNetwork.cpp

${OBJECTDIR}/_ext/bac85f60/OpenSMOKE_KPP_ReactorNetwork_Residuals.o: ../../src/kpp/OpenSMOKE_KPP_ReactorNetwork_Residuals.cpp 
	${MKDIR} -p ${OBJECTDIR}/_ext/bac85f60
	${RM} "$@.d"
	$(COMPILE.cc) -O2 -I../../../../NumericalLibraries/BzzMath/BzzMath6-dev/hpp -I../../../../NumericalLibraries/OpenSMOKE/hpp -I../../../../NumericalRecipes/hpp -I../../../../MUMPS_4.9.2/libraries/include -I../../../../NumericalLibraries/intel/mkl/include -I../../src -I../../../../NumericalLibraries/lis-1.2.122/lis-nested-complete/include -I/home/astagni/NumericalLibraries/petsc/petsc-3.7.3/include -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/_ext/bac85f60/OpenSMOKE_KPP_ReactorNetwork_Residuals.o ../../src/kpp/OpenSMOKE_KPP_ReactorNetwork_Residuals.cpp

${OBJECTDIR}/_ext/bac85f60/OpenSMOKE_KPP_SingleReactor.o: ../../src/kpp/OpenSMOKE_KPP_SingleReactor.cpp 
	${MKDIR} -p ${OBJECTDIR}/_ext/bac85f60
	${RM} "$@.d"
	$(COMPILE.cc) -O2 -I../../../../NumericalLibraries/BzzMath/BzzMath6-dev/hpp -I../../../../NumericalLibraries/OpenSMOKE/hpp -I../../../../NumericalRecipes/hpp -I../../../../MUMPS_4.9.2/libraries/include -I../../../../NumericalLibraries/intel/mkl/include -I../../src -I../../../../NumericalLibraries/lis-1.2.122/lis-nested-complete/include -I/home/astagni/NumericalLibraries/petsc/petsc-3.7.3/include -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/_ext/bac85f60/OpenSMOKE_KPP_SingleReactor.o ../../src/kpp/OpenSMOKE_KPP_SingleReactor.cpp

${OBJECTDIR}/_ext/bac85f60/OpenSMOKE_KPP_SingleReactorStatistics.o: ../../src/kpp/OpenSMOKE_KPP_SingleReactorStatistics.cpp 
	${MKDIR} -p ${OBJECTDIR}/_ext/bac85f60
	${RM} "$@.d"
	$(COMPILE.cc) -O2 -I../../../../NumericalLibraries/BzzMath/BzzMath6-dev/hpp -I../../../../NumericalLibraries/OpenSMOKE/hpp -I../../../../NumericalRecipes/hpp -I../../../../MUMPS_4.9.2/libraries/include -I../../../../NumericalLibraries/intel/mkl/include -I../../src -I../../../../NumericalLibraries/lis-1.2.122/lis-nested-complete/include -I/home/astagni/NumericalLibraries/petsc/petsc-3.7.3/include -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/_ext/bac85f60/OpenSMOKE_KPP_SingleReactorStatistics.o ../../src/kpp/OpenSMOKE_KPP_SingleReactorStatistics.cpp

${OBJECTDIR}/_ext/bac85f60/OpenSMOKE_KPP_SingleReactor_KineticsManager.o: ../../src/kpp/OpenSMOKE_KPP_SingleReactor_KineticsManager.cpp 
	${MKDIR} -p ${OBJECTDIR}/_ext/bac85f60
	${RM} "$@.d"
	$(COMPILE.cc) -O2 -I../../../../NumericalLibraries/BzzMath/BzzMath6-dev/hpp -I../../../../NumericalLibraries/OpenSMOKE/hpp -I../../../../NumericalRecipes/hpp -I../../../../MUMPS_4.9.2/libraries/include -I../../../../NumericalLibraries/intel/mkl/include -I../../src -I../../../../NumericalLibraries/lis-1.2.122/lis-nested-complete/include -I/home/astagni/NumericalLibraries/petsc/petsc-3.7.3/include -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/_ext/bac85f60/OpenSMOKE_KPP_SingleReactor_KineticsManager.o ../../src/kpp/OpenSMOKE_KPP_SingleReactor_KineticsManager.cpp

${OBJECTDIR}/_ext/bc30b0d1/OpenSMOKE_DirectLinearSolver_Unsymmetric.o: ../../src/linear_solvers/OpenSMOKE_DirectLinearSolver_Unsymmetric.cpp 
	${MKDIR} -p ${OBJECTDIR}/_ext/bc30b0d1
	${RM} "$@.d"
	$(COMPILE.cc) -O2 -I../../../../NumericalLibraries/BzzMath/BzzMath6-dev/hpp -I../../../../NumericalLibraries/OpenSMOKE/hpp -I../../../../NumericalRecipes/hpp -I../../../../MUMPS_4.9.2/libraries/include -I../../../../NumericalLibraries/intel/mkl/include -I../../src -I../../../../NumericalLibraries/lis-1.2.122/lis-nested-complete/include -I/home/astagni/NumericalLibraries/petsc/petsc-3.7.3/include -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/_ext/bc30b0d1/OpenSMOKE_DirectLinearSolver_Unsymmetric.o ../../src/linear_solvers/OpenSMOKE_DirectLinearSolver_Unsymmetric.cpp

${OBJECTDIR}/_ext/bc30b0d1/OpenSMOKE_LAPACK_Dense.o: ../../src/linear_solvers/OpenSMOKE_LAPACK_Dense.cpp 
	${MKDIR} -p ${OBJECTDIR}/_ext/bc30b0d1
	${RM} "$@.d"
	$(COMPILE.cc) -O2 -I../../../../NumericalLibraries/BzzMath/BzzMath6-dev/hpp -I../../../../NumericalLibraries/OpenSMOKE/hpp -I../../../../NumericalRecipes/hpp -I../../../../MUMPS_4.9.2/libraries/include -I../../../../NumericalLibraries/intel/mkl/include -I../../src -I../../../../NumericalLibraries/lis-1.2.122/lis-nested-complete/include -I/home/astagni/NumericalLibraries/petsc/petsc-3.7.3/include -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/_ext/bc30b0d1/OpenSMOKE_LAPACK_Dense.o ../../src/linear_solvers/OpenSMOKE_LAPACK_Dense.cpp

${OBJECTDIR}/_ext/bc30b0d1/OpenSMOKE_LIS_Unsymmetric.o: ../../src/linear_solvers/OpenSMOKE_LIS_Unsymmetric.cpp 
	${MKDIR} -p ${OBJECTDIR}/_ext/bc30b0d1
	${RM} "$@.d"
	$(COMPILE.cc) -O2 -I../../../../NumericalLibraries/BzzMath/BzzMath6-dev/hpp -I../../../../NumericalLibraries/OpenSMOKE/hpp -I../../../../NumericalRecipes/hpp -I../../../../MUMPS_4.9.2/libraries/include -I../../../../NumericalLibraries/intel/mkl/include -I../../src -I../../../../NumericalLibraries/lis-1.2.122/lis-nested-complete/include -I/home/astagni/NumericalLibraries/petsc/petsc-3.7.3/include -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/_ext/bc30b0d1/OpenSMOKE_LIS_Unsymmetric.o ../../src/linear_solvers/OpenSMOKE_LIS_Unsymmetric.cpp

${OBJECTDIR}/_ext/bc30b0d1/OpenSMOKE_PARDISO_Unsymmetric.o: ../../src/linear_solvers/OpenSMOKE_PARDISO_Unsymmetric.cpp 
	${MKDIR} -p ${OBJECTDIR}/_ext/bc30b0d1
	${RM} "$@.d"
	$(COMPILE.cc) -O2 -I../../../../NumericalLibraries/BzzMath/BzzMath6-dev/hpp -I../../../../NumericalLibraries/OpenSMOKE/hpp -I../../../../NumericalRecipes/hpp -I../../../../MUMPS_4.9.2/libraries/include -I../../../../NumericalLibraries/intel/mkl/include -I../../src -I../../../../NumericalLibraries/lis-1.2.122/lis-nested-complete/include -I/home/astagni/NumericalLibraries/petsc/petsc-3.7.3/include -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/_ext/bc30b0d1/OpenSMOKE_PARDISO_Unsymmetric.o ../../src/linear_solvers/OpenSMOKE_PARDISO_Unsymmetric.cpp

# Subprojects
.build-subprojects:

# Clean Targets
.clean-conf: ${CLEAN_SUBPROJECTS}
	${RM} -r ${CND_BUILDDIR}/${CND_CONF}
	${RM} ../../exe/OpenSMOKE_KinPP_Splitting_NetBeans.sh

# Subprojects
.clean-subprojects:

# Enable dependency checking
.dep.inc: .depcheck-impl

include .dep.inc
