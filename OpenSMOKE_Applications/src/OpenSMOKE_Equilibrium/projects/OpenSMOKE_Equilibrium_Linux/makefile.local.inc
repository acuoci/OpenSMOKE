# ------------------------------------------------------------------------------------------------
# OBJECTS
# ------------------------------------------------------------------------------------------------
OBJS         = 	$(OBJ_PATH)/OpenSMOKE_Equilibrium.o $(OBJ_PATH)/OpenSMOKE_EquilibriumModule.o


# ------------------------------------------------------------------------------------------------
# LIBRARIES
# ------------------------------------------------------------------------------------------------
LOC_LIBS	=  -lOpenSMOKE_IdealReactors_$(STRING_NAME)   \
                   -lOpenSMOKE_Flame1D_$(STRING_NAME)         \
                   -lOpenSMOKE_Engine_$(STRING_NAME)          \
                   -lOpenSMOKE_PreProcessing_$(STRING_NAME)   \
			 -lOpenSMOKE_MixManager_$(STRING_NAME)      \
                   -lOpenSMOKE_Flamelet_$(STRING_NAME)        \
			 -lOpenSMOKE_AddOns_$(STRING_NAME)          \
                   -lOpenSMOKE_Basic_$(STRING_NAME)           \
			 -lOpenSMOKE_IdealReactors_$(STRING_NAME)   \
                   -lOpenSMOKE_Interfaces_$(STRING_NAME)      \
		       -lOpenSMOKE_PlugFlowReactor_$(STRING_NAME) \
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

default: OpenSMOKE_Equilibrium

$(OBJ_PATH)/OpenSMOKE_Equilibrium.o: 	../../src/OpenSMOKE_Equilibrium.cpp 
				$(CPP) ../../src/OpenSMOKE_Equilibrium.cpp  -o $(OBJ_PATH)/OpenSMOKE_Equilibrium.o

$(OBJ_PATH)/OpenSMOKE_EquilibriumModule.o: 	../../src/OpenSMOKE_EquilibriumModule.cpp 
				$(CPP) ../../src/OpenSMOKE_EquilibriumModule.cpp  -o $(OBJ_PATH)/OpenSMOKE_EquilibriumModule.o	

OpenSMOKE_Equilibrium:  	$(OBJS)  
		$(LINKER) $(OBJS) -o ../../exe/OpenSMOKE_Equilibrium_$(STRING_NAME).sh $(LIBS)

clean: ; @rm $(OBJ_PATH)/*.o ../../exe/OpenSMOKE_Equilibrium_$(STRING_NAME).sh
