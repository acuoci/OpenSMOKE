# ------------------------------------------------------------------------------------------------
# OBJECTS
# ------------------------------------------------------------------------------------------------
OBJS     =  $(OBJ_PATH)/main.o \
            $(OBJ_PATH)/CRECK_Optimizer_Class.o \
            $(OBJ_PATH)/ExperimentClass.o

# ------------------------------------------------------------------------------------------------
# LIBRARIES
# ------------------------------------------------------------------------------------------------
LOC_LIBS	=  -lOpenSMOKE_PlugFlowReactor_$(STRING_NAME) \
                   -lOpenSMOKE_CSTR_$(STRING_NAME) \
                   -lOpenSMOKE_Flame1D_$(STRING_NAME)         \
		   -lOpenSMOKE_LiquidProperties_$(STRING_NAME)         \
                   -lOpenSMOKE_Basic_$(STRING_NAME)           \
                   -lOpenSMOKE_PreProcessing_$(STRING_NAME)           \
                   -lOpenSMOKE_CHEMKINInterpreter_$(STRING_NAME)           \
                   -lOpenSMOKE_Engine_$(STRING_NAME)          \
                   -lOpenSMOKE_AddOns_$(STRING_NAME)          \
                   -lOpenSMOKE_Flamelet_$(STRING_NAME)        \
                   -lOpenSMOKE_IdealReactors_$(STRING_NAME)   \
                   -lOpenSMOKE_Interfaces_$(STRING_NAME)      \
                   -lOpenSMOKE_QMOM_$(STRING_NAME)                           

# ------------------------------------------------------------------------------------------------
# INCLUDE FILES
# ------------------------------------------------------------------------------------------------ 
INCDIRS     = -I../../src $(BZZ_INCLUDE) $(MKL_INCLUDE) $(OPENSMOKE_INCLUDE) $(NR_INCLUDE)
  
# ------------------------------------------------------------------------------------------------
# LIBRARIES
# ------------------------------------------------------------------------------------------------
LIBS_PATH 	      = -L$(BZZ_PATH_LIB) -L$(MKL_PATH_LIB) -L$(OPENSMOKE_PATH_LIB) -L$(NR_PATH_LIB)
LIBS 		      = $(LOC_LIBS) $(BZZ_LIBS) $(NR_LIBS) $(MKL_LIBS)

# ------------------------------------------------------------------------------------------------
# COMPILER
# ------------------------------------------------------------------------------------------------
CPP       = $(CXX) 	$(CXX_FLAGS) -c 	$(INCDIRS)
LINKER    = $(CXX) 	$(CXX_FLAGS) 		$(INCDIRS) 	$(LIBS_PATH)

# ------------------------------------------------------------------------------------------------
# DETAILS
# ------------------------------------------------------------------------------------------------

default: OpenSMOKE_Optimizer

$(OBJ_PATH)/main.o: 	../../src/main.cpp 
					$(CPP) ../../src/main.cpp  -o $(OBJ_PATH)/main.o

$(OBJ_PATH)/ExperimentClass.o: 	../../src/ExperimentClass.cpp ../../src/ExperimentClass.h
					$(CPP) ../../src/ExperimentClass.cpp  -o $(OBJ_PATH)/ExperimentClass.o


$(OBJ_PATH)/CRECK_Optimizer_Class.o: ../../src/CRECK_Optimizer_Class.cpp ../../src/CRECK_Optimizer_Class.h
					$(CPP) ../../src/CRECK_Optimizer_Class.cpp  -o $(OBJ_PATH)/CRECK_Optimizer_Class.o

	
OpenSMOKE_Optimizer:  	$(OBJS)  
			$(LINKER) $(OBJS) -o ../../exe/OpenSMOKE_Optimizer_$(STRING_NAME).sh $(LIBS)

clean: ; @rm $(OBJ_PATH)/*.o ../../exe/OpenSMOKE_Optimizer_$(STRING_NAME).sh
