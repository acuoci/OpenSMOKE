include ../../../make.inc

# ------------------------------------------------------------------------------------------------
# OBJECTS
# ------------------------------------------------------------------------------------------------
OBJS     =  $(OBJ_PATH)/NOxPostProcessor.o \
            $(OBJ_PATH)/MyOdeSystemObjectOneCSTR.o \
            $(OBJ_PATH)/MyOdeSystemObjectAllCSTR.o \
            $(OBJ_PATH)/myBzzCSTRNetwork.o

# ------------------------------------------------------------------------------------------------
# LIBRARIES
# ------------------------------------------------------------------------------------------------
ifeq ($(COMPILER_NAME),intel)
LOC_LIBS 	=  -lOpenSMOKE_Engine_INTEL           \
                   -lOpenSMOKE_Basic_INTEL         \
                   -lOpenSMOKE_AddOns_INTEL                      	     
endif	  

ifeq ($(COMPILER_NAME),gnu)
LOC_LIBS	=  -lOpenSMOKE_Basic_GNU           \
                   -lOpenSMOKE_Engine_GNU          \
                   -lOpenSMOKE_AddOns_GNU                 
endif
   
  
# ------------------------------------------------------------------------------------------------
# INCLUDE FILES
# ------------------------------------------------------------------------------------------------ 
INCDIRS     = -I../../src $(BZZ_INCLUDE) $(MKL_INCLUDE) $(OPENSMOKE_INCLUDE) $(NR_INCLUDE)
  
# ------------------------------------------------------------------------------------------------
# LIBRARIES
# ------------------------------------------------------------------------------------------------
LIBS_PATH 	      = -L$(BZZ_PATH_LIB) -L$(MKL_PATH_LIB) -L$(OPENSMOKE_PATH_LIB) -L$(NR_PATH_LIB)
LIBS 		      = $(LOC_LIBS) $(BZZ_LIBS) $(NR_LIBS) 

# ------------------------------------------------------------------------------------------------
# COMPILER
# ------------------------------------------------------------------------------------------------
CPP       = $(CXX) 	$(CXX_FLAGS) -c 	$(INCDIRS)
LINKER    = $(CXX) 	$(CXX_FLAGS) 		$(INCDIRS) 	$(LIBS_PATH)

# ------------------------------------------------------------------------------------------------
# DETAILS
# ------------------------------------------------------------------------------------------------

default: OpenSMOKE_KinPP_September2007

$(OBJ_PATH)/NOxPostProcessor.o: 	../../src/NOxPostProcessor.cpp 
				$(CPP) ../../src/NOxPostProcessor.cpp  -o $(OBJ_PATH)/NOxPostProcessor.o

$(OBJ_PATH)/MyOdeSystemObjectOneCSTR.o: 	../../src/MyOdeSystemObjectOneCSTR.cpp 
				$(CPP) ../../src/MyOdeSystemObjectOneCSTR.cpp  -o $(OBJ_PATH)/MyOdeSystemObjectOneCSTR.o

$(OBJ_PATH)/MyOdeSystemObjectAllCSTR.o: 	../../src/MyOdeSystemObjectAllCSTR.cpp 
				$(CPP) ../../src/MyOdeSystemObjectAllCSTR.cpp  -o $(OBJ_PATH)/MyOdeSystemObjectAllCSTR.o

$(OBJ_PATH)/myBzzCSTRNetwork.o: 	../../src/myBzzCSTRNetwork.cpp 
				$(CPP) ../../src/myBzzCSTRNetwork.cpp  -o $(OBJ_PATH)/myBzzCSTRNetwork.o
	
OpenSMOKE_KinPP_September2007:  	$(OBJS)  
			$(LINKER) $(OBJS) -o ../../exe/OpenSMOKE_KinPP_September2007_$(STRING_NAME).sh $(LIBS)

clean: ; @rm $(OBJ_PATH)/*.o ../../exe/OpenSMOKE_KinPP_September2007_$(STRING_NAME).sh
