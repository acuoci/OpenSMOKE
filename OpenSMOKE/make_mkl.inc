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
BZZ_PATH 	       =  /software/chimica2/libraries/bzzmath/bzzmath-6.0
OPENSMOKE_PATH         =  $(HOME)/Development/ExternalNumericalLibraries/OpenSMOKE/OpenSMOKE
NUMERICALRECIPES_PATH  =  /software/chimica2/libraries/nr/nr-gcc-6.3.0
MKL_PATH               =  $(HOME)/intel/mkl

#Options
#OPTIONS_INTEL   =  -O2 -static 
#OPTIONS_GNU     =  -O2 -static -Wall -fexpensive-optimizations -funroll-loops
OPTIONS_INTEL    =  -O3 -Wall -m64 -DUSE_INTEL_MKL=1
OPTIONS_GNU      =  -O3 -Wall -m64 -DUSE_INTEL_MKL=1

# ------------------------------------------------------------------------------------------------
# LIBRARY NAMES
# ------------------------------------------------------------------------------------------------

#BzzMath Libraries
BZZ_LIBS_INTEL  = -lBzzMath6_INTEL
BZZ_LIBS_GNU    = -lBzzMath6_GNU

#NR libraries
NR_LIBS_INTEL   =  -lnr_INTEL   
NR_LIBS_GNU     =  -lnr_GNU


# ------------------------------------------------------------------------------------------------
# COMPILERS
# ------------------------------------------------------------------------------------------------
ifeq ($(COMPILER_NAME),intel)
OBJ_PATH    = intel_mkl
STRING_NAME = INTEL_MKL
CXX         = icpc
CXX_FLAGS   = $(OPTIONS_INTEL)
BZZ_LIBS    = $(BZZ_LIBS_INTEL)
NR_LIBS     = $(NR_LIBS_INTEL)
endif	

ifeq ($(COMPILER_NAME),gnu)
OBJ_PATH    = gnu_mkl
STRING_NAME = GNU_MKL
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
NR_INCLUDE           = -I$(NUMERICALRECIPES_PATH)/include
OPENSMOKE_INCLUDE    = -I$(OPENSMOKE_PATH)/hpp

# ------------------------------------------------------------------------------------------------
# LIBRARY PATHS
# ------------------------------------------------------------------------------------------------
BZZ_PATH_LIB           = $(BZZ_PATH)/lib/release-6.3.0
NR_PATH_LIB            = $(NUMERICALRECIPES_PATH)/lib/release
OPENSMOKE_PATH_LIB     = $(OPENSMOKE_PATH)/lib/linux/$(OBJ_PATH)
MKL_PATH_LIB           = $(MKL_PATH)/lib/intel64

# ------------------------------------------------------------------------------------------------
# COMPILER INSTRUCTIONS
# ------------------------------------------------------------------------------------------------
CCP = $(CXX) $(CXX_FLAGS) -c $(BZZ_INCLUDE)  $(OPENSMOKE_INCLUDE) $(MKL_INCLUDE) $(NR_INCLUDE) 



