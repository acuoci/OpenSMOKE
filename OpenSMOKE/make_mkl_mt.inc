############################################################################
#
#  Program:         OpenSMOKE
#
#  Module:          make.inc
#
#  Purpose:         Top-level Definitions
#
#  Creation date:   January 9 2010
#
############################################################################

#Paths 
MYHOME		       =  $(HOME)
BZZ_PATH 	       =  $(MYHOME)/NumericalLibraries/BzzMath6
OPENSMOKE_PATH         =  $(MYHOME)/NumericalLibraries/OpenSMOKE/OpenSMOKE
NUMERICALRECIPES_PATH  =  $(MYHOME)/NumericalLibraries/NumericalRecipes

#Options
OPTIONS_INTEL    =  -O3 -Wall -m64 -DUSE_INTEL_MKL=1 -openmp -I/usr/include/x86_64-linux-gnu/c++/4.8
OPTIONS_GNU      =  -O3 -Wall -m64 -DUSE_INTEL_MKL=1 -fopenmp

# ------------------------------------------------------------------------------------------------
# LIBRARY NAMES
# ------------------------------------------------------------------------------------------------

#BzzMath Libraries
BZZ_LIBS_INTEL  = -lBzzMath6_INTEL_MKL_MT
BZZ_LIBS_GNU    = -lBzzMath6_GNU_MKL_MT

#NR libraries
NR_LIBS_INTEL   =  -lNumericalRecipes_INTEL   
NR_LIBS_GNU     =  -lNumericalRecipes_GNU


# ------------------------------------------------------------------------------------------------
# COMPILERS
# ------------------------------------------------------------------------------------------------
ifeq ($(COMPILER_NAME),intel)
OBJ_PATH    = intel_mkl_mt
STRING_NAME = INTEL_MKL_MT
CXX         = icpc
CXX_FLAGS   = $(OPTIONS_INTEL)
BZZ_LIBS    = $(BZZ_LIBS_INTEL)
NR_LIBS     = $(NR_LIBS_INTEL)
endif	

ifeq ($(COMPILER_NAME),gnu)
OBJ_PATH    = gnu_mkl_mt
STRING_NAME = GNU_MKL_MT
CXX         = g++
CXX_FLAGS   = $(OPTIONS_GNU)
BZZ_LIBS    = $(BZZ_LIBS_GNU) 
NR_LIBS     = $(NR_LIBS_GNU)
endif	

# ------------------------------------------------------------------------------------------------
# INCLUDE
# ------------------------------------------------------------------------------------------------
BZZ_INCLUDE          = -I$(BZZ_PATH)/hpp
MKL_INCLUDE          = -I$(MKL_PATH)/include
NR_INCLUDE           = -I$(NUMERICALRECIPES_PATH)/hpp
OPENSMOKE_INCLUDE    = -I$(OPENSMOKE_PATH)/hpp

# ------------------------------------------------------------------------------------------------
# LIBRARY PATHS
# ------------------------------------------------------------------------------------------------
BZZ_PATH_LIB           = $(BZZ_PATH)/lib/release
NR_PATH_LIB            = $(NUMERICALRECIPES_PATH)/lib/release
OPENSMOKE_PATH_LIB     = $(OPENSMOKE_PATH)/lib/linux/$(OBJ_PATH)
MKL_PATH_LIB           = $(MKL_PATH)/lib/intel64

# ------------------------------------------------------------------------------------------------
# COMPILER INSTRUCTIONS
# ------------------------------------------------------------------------------------------------
CCP = $(CXX) $(CXX_FLAGS) -c $(BZZ_INCLUDE)  $(OPENSMOKE_INCLUDE) $(MKL_INCLUDE) $(NR_INCLUDE) 



