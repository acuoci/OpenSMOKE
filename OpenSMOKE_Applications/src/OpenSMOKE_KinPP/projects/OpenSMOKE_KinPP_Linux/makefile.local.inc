FILE_NAME    =  OpenSMOKE_KinPP

# ------------------------------------------------------------------------------------------------
# OBJECTS
# ------------------------------------------------------------------------------------------------
OBJS         = 	$(OBJ_PATH)/$(FILE_NAME).o

# ------------------------------------------------------------------------------------------------
# LIBRARIES
# ------------------------------------------------------------------------------------------------
LOC_LIBS	=  -lOpenSMOKE_IdealReactors_$(STRING_NAME)   \
                   -lOpenSMOKE_Basic_$(STRING_NAME)           \
                   -lOpenSMOKE_Engine_$(STRING_NAME)          \
                   -lOpenSMOKE_AddOns_$(STRING_NAME)          \
                   -lOpenSMOKE_Interfaces_$(STRING_NAME)      \
		   -lOpenSMOKE_KinPP_$(STRING_NAME)      \
                   -lOpenSMOKE_Distributions_$(STRING_NAME) \
		   -lOpenSMOKE_SymbolicKinetics_$(STRING_NAME) 

  
# ------------------------------------------------------------------------------------------------
# INCLUDE FILES
# ------------------------------------------------------------------------------------------------ 
INCDIRS     = -I../../src $(BZZ_INCLUDE) $(MKL_INCLUDE) $(OPENSMOKE_INCLUDE) $(NR_INCLUDE)
  
# ------------------------------------------------------------------------------------------------
# LIBRARIES
# ------------------------------------------------------------------------------------------------
LIBS_PATH 	      = -L$(BZZ_PATH_LIB) -L$(MKL_PATH_LIB) -L$(OPENSMOKE_PATH_LIB) -L$(NR_PATH_LIB)
LIBS 		      = $(LOC_LIBS) $(BZZ_LIBS) $(NR_LIBS) $(MKL_LIBS) -lpthread -lm -ldl 

# ------------------------------------------------------------------------------------------------
# COMPILER
# ------------------------------------------------------------------------------------------------
CPP       = $(CXX) 	$(CXX_FLAGS) -c 	$(INCDIRS)
LINKER    = $(CXX) 	$(CXX_FLAGS) 		$(INCDIRS) 	$(LIBS_PATH)

# ------------------------------------------------------------------------------------------------
# DETAILS
# ------------------------------------------------------------------------------------------------

default: $(FILE_NAME)

$(OBJ_PATH)/$(FILE_NAME).o: 	../../src/$(FILE_NAME).cpp 
				$(CPP) ../../src/$(FILE_NAME).cpp  -o $(OBJ_PATH)/$(FILE_NAME).o
	
$(FILE_NAME):  	$(OBJS)  
		$(LINKER) $(OBJS) -o ../../exe/$(FILE_NAME)_$(STRING_NAME).sh $(LIBS)

clean: ; @rm $(OBJ_PATH)/*.o ../../exe/$(FILE_NAME)_$(STRING_NAME).sh
