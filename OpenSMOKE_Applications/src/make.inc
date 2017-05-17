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

#  Compiler:   intel || gnu

COMPILER_NAME  	      = gnu

MYHOME		       =  $(HOME)
BZZ_PATH 	       =  /software/chimica2/libraries/bzzmath/bzzmath-6.0
OPENSMOKE_PATH         =  $(HOME)/Development/ExternalNumericalLibraries/OpenSMOKE/OpenSMOKE
NUMERICALRECIPES_PATH  =  /software/chimica2/libraries/nr/nr-gcc-6.3.0
MKL_PATH               =  $(HOME)/intel/mkl
#SYSTEM_LIB_PATH       =  /usr/lib/x86_64-redhat-linux5E/lib64/
SYSTEM_LIB_PATH        =  /

#Paths
LIS_PATH               =  /software/chimica2/libraries/lis/lis-1.7.29/serial-gcc-6.3.0
MUMPS_PATH             =  /software/chimica2/libraries/mumps/mumps-5.1.1/serial-pord-gcc-6.3.0
METIS_PATH             =  /software/chimica2/libraries/metis/metis-5.1.0-gcc-6.3.0

#Options
OPTIONS_INTEL   =  -O3 -DUSE_INTEL_MKL=0 
OPTIONS_GNU     =  -O3 -m64 -Wall -fexpensive-optimizations -funroll-loops -DUSE_INTEL_MKL=0	

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
MUMPS_LIBS      = -ldmumps -lmumps_common -lpord -lmpiseq
LIS_LIBS        = -llis
METIS_LIBS      = -lmetis

# ------------------------------------------------------------------------------------------------
# COMPILERS
# ------------------------------------------------------------------------------------------------
ifeq ($(COMPILER_NAME),intel)

OBJ_PATH    = intel
STRING_NAME = INTEL
CXX         = icpc
CXX_FLAGS   = $(OPTIONS_INTEL)
BZZ_LIBS    = $(BZZ_LIBS_INTEL)
NR_LIBS     = $(NR_LIBS_INTEL)
endif	

ifeq ($(COMPILER_NAME),gnu)

OBJ_PATH    = gnu
STRING_NAME = GNU
CXX         = g++
CXX_FLAGS   = $(OPTIONS_GNU) -static
BZZ_LIBS    = $(BZZ_LIBS_GNU) 
NR_LIBS     = $(NR_LIBS_GNU)
endif	

# ------------------------------------------------------------------------------------------------
# INCLUDE
# ------------------------------------------------------------------------------------------------
BZZ_INCLUDE          = -I$(BZZ_PATH)/hpp
MKL_INCLUDE          = -I$(MKL_PATH)/include
MKL_ROOT	     = ${MKL_PATH}
NR_INCLUDE           = -I$(NUMERICALRECIPES_PATH)/include
OPENSMOKE_INCLUDE    = -I$(OPENSMOKE_PATH)/hpp
MUMPS_INCLUDE        = -I$(MUMPS_PATH)/include
LIS_INCLUDE          = -I$(LIS_PATH)/include
METIS_INCLUDE        = -I$(METIS_PATH)/include

# ------------------------------------------------------------------------------------------------
# LIBRARY PATHS
# ------------------------------------------------------------------------------------------------
BZZ_PATH_LIB           = $(BZZ_PATH)/lib/release-6.3.0
NR_PATH_LIB            = $(NUMERICALRECIPES_PATH)/lib
OPENSMOKE_PATH_LIB     = $(OPENSMOKE_PATH)/lib/linux/$(OBJ_PATH)
MKL_PATH_LIB           = $(MKL_PATH)/lib/intel64
MUMPS_PATH_LIB         = $(MUMPS_PATH)/lib
LIS_PATH_LIB           = $(LIS_PATH)/lib
METIS_PATH_LIB         = $(METIS_PATH)/lib

MKL_LIBS               = -L$(SYSTEM_LIB_PATH)

