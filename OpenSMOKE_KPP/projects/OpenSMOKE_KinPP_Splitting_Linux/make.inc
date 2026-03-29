############################################################################
#
#  Program:         OpenSMOKE
#
#  Module:          make.inc
#
#  Purpose:         Top-level Definitions
#
#  Creation date:   May 28 2012
#
############################################################################

#  Compiler:  intel || gnu

COMPILER_NAME  	       = gnu

#Paths
BZZ_PATH               =  /software/chimica2/libraries/bzzmath/bzzmath-6.0/
OPENSMOKE_PATH         =  ${HOME}/Development/ExternalNumericalLibraries/OpenSMOKE/OpenSMOKE/
NUMERICALRECIPES_PATH  =  /software/chimica2/libraries/nr/nr-gcc-4.8.2/
MKL_PATH               =  /software/chimica2/libraries/intel/oneapi/mkl/latest/
LIS_PATH               =  /software/chimica2/libraries/lis/lis-1.7.29/mpi-long-gcc-6.3.0/
METIS_PATH             =  /software/chimica2/libraries/metis/metis-5.1.0-gcc-6.3.0/
BLAS_PATH              =  /software/chimica2/libraries/intel/oneapi/mkl/latest/
PETSC_PATH             =  /software/chimica2/libraries/petsc/petsc-3.7.7-openmpi-1.1.0-64bit
#MUMPS_PATH            =  ${HOME}/NumericalLibraries/MUMPS/MUMPS_4.9.2_SEQ_MKL_PORD

#Options
OPTIONS_INTEL   =  -O2 
OPTIONS_GNU     =  -O2 -m64 -fexpensive-optimizations -funroll-loops -DUSE_MPI -D_LONGLONG -c -Wall -Wwrite-strings -Wno-strict-aliasing -Wno-unknown-pragmas -g3 -Wfatal-errors

# ------------------------------------------------------------------------------------------------
# LIBRARY NAMES
# ------------------------------------------------------------------------------------------------

#BzzMath Libraries
BZZ_LIBS_INTEL  = -lBzzMath6_INTEL
BZZ_LIBS_GNU    = -lBzzMath6_GNU

#NR libraries
NR_LIBS_INTEL   =  -lnr_INTEL   
NR_LIBS_GNU     =  -lnr_GNU

#MUMPS libraries
#MUMPS_LIBS      = -ldmumps -lmumps_common -lpord $(MUMPS_PATH)/libseq/libmpiseq.a 
LIS_LIBS        = -llis
METIS_LIBS      = -lmetis
BLAS_LIBS       = -Wl,-rpath,$(BLAS_PATH)/lib/intel64 -L$(BLAS_PATH)/lib/intel64 -lmkl_intel_lp64 -lmkl_gnu_thread -lmkl_core
BLAS_LIBS 	= -L/nothing/ 
#PETSC libraries
PETSC_LIBS	= -lpetsc -lpthread

# ------------------------------------------------------------------------------------------------
# COMPILERS
# ------------------------------------------------------------------------------------------------
ifeq ($(COMPILER_NAME),intel)
OBJ_PATH    = intel
STRING_NAME = INTEL
CXX         = icpc
CXX_FLAGS   = $(OPTIONS_INTEL)
CXX_ADD_FLAGS = -openmp
BZZ_LIBS    = $(BZZ_LIBS_INTEL)
NR_LIBS     = $(NR_LIBS_INTEL)
endif	

ifeq ($(COMPILER_NAME),gnu)
OBJ_PATH    = gnu
STRING_NAME = GNU
CXX         = mpic++
CXX_FLAGS   = $(OPTIONS_GNU)
CXX_ADD_FLAGS = -fopenmp
BZZ_LIBS    = $(BZZ_LIBS_GNU) 
NR_LIBS     = $(NR_LIBS_GNU)
endif	

# ------------------------------------------------------------------------------------------------
# INCLUDE
# ------------------------------------------------------------------------------------------------
BZZ_INCLUDE          = -I$(BZZ_PATH)/hpp/release-4.8.5
MKL_INCLUDE          = -I$(MKL_PATH)/include
MKL_ROOT	     = ${MKL_PATH}
NR_INCLUDE           = -I$(NUMERICALRECIPES_PATH)/hpp
OPENSMOKE_INCLUDE    = -I$(OPENSMOKE_PATH)/hpp
#MUMPS_INCLUDE       = -I$(MUMPS_PATH)/include
LIS_INCLUDE          = -I$(LIS_PATH)/include
METIS_INCLUDE        = -I$(METIS_PATH)/include
PETSC_INCLUDE	     = -I$(PETSC_PATH)/include -I$(PETSC_PATH)/linux-release/include

# ------------------------------------------------------------------------------------------------
# LIBRARY PATHS
# ------------------------------------------------------------------------------------------------
BZZ_PATH_LIB           = $(BZZ_PATH)/lib/release-4.8.5
NR_PATH_LIB            = $(NUMERICALRECIPES_PATH)/lib/
OPENSMOKE_PATH_LIB     = $(OPENSMOKE_PATH)/lib/linux/gnu_mkl
MKL_PATH_LIB           = $(MKL_PATH)/lib/intel64
#MUMPS_PATH_LIB         = $(MUMPS_PATH)/lib
LIS_PATH_LIB           = $(LIS_PATH)/lib
METIS_PATH_LIB         = $(METIS_PATH)/lib
BLAS_PATH_LIB          = $(BLAS_PATH)/lib/intel64
PETSC_PATH_LIB	       = $(PETSC_PATH)/lib

#MKL_LIBS        = -Wl,--start-group $(MKL_ROOT)/lib/intel64/libmkl_intel_lp64.a $(MKL_ROOT)/lib/intel64/libmkl_gnu_thread.a $(MKL_ROOT)/lib/intel64/libmkl_core.a -Wl,--end-group -lgomp -lpthread
MKL_LIBS =  -m64  -L$(MKL_ROOT)/lib/intel64  -Wl,--no-as-needed -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lgomp -lpthread -lm -ldl

