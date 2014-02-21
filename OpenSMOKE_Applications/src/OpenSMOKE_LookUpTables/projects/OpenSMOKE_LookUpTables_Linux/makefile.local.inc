FILE_NAME    =  OpenSMOKE_LookUpTables

# ------------------------------------------------------------------------------------------------
# OBJECTS
# ------------------------------------------------------------------------------------------------
OBJS     =  $(OBJ_PATH)/Main.o \
            $(OBJ_PATH)/Utilities.o \
            $(OBJ_PATH)/OpenSMOKE_LookUp_Table_Flame.o \
            $(OBJ_PATH)/sources_in_flamelet_library.o \
            $(OBJ_PATH)/OpenSMOKE_LookUp_Table_Manager.o \
            $(OBJ_PATH)/OpenSMOKE_LookUp_Table_Executables.o \
            $(OBJ_PATH)/flamelet_group.o \
            $(OBJ_PATH)/flamelet.o


# ------------------------------------------------------------------------------------------------
# LIBRARIES
# ------------------------------------------------------------------------------------------------
LOC_LIBS	=  -lOpenSMOKE_Basic_$(STRING_NAME)          \
                   -lOpenSMOKE_Engine_$(STRING_NAME)         \
                   -lOpenSMOKE_Distributions_$(STRING_NAME)   \
                   -lOpenSMOKE_AddOns_$(STRING_NAME)          \
                   -lOpenSMOKE_Flamelet_$(STRING_NAME)       \
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

default: OpenSMOKE_LookUpTables

$(OBJ_PATH)/Main.o: 	../../src/Main.cpp 
			$(CPP) ../../src/Main.cpp  -o $(OBJ_PATH)/Main.o

$(OBJ_PATH)/Utilities.o: 	../../src/Utilities.cpp ../../src/Utilities.h
				$(CPP) ../../src/Utilities.cpp  -o $(OBJ_PATH)/Utilities.o


$(OBJ_PATH)/OpenSMOKE_LookUp_Table_Flame.o: ../../src/OpenSMOKE_LookUp_Table_Flame.cpp ../../src/OpenSMOKE_LookUp_Table_Flame.h
					$(CPP) ../../src/OpenSMOKE_LookUp_Table_Flame.cpp  -o $(OBJ_PATH)/OpenSMOKE_LookUp_Table_Flame.o

$(OBJ_PATH)/sources_in_flamelet_library.o: ../../src/sources_in_flamelet_library.cpp ../../src/sources_in_flamelet_library.h
					$(CPP) ../../src/sources_in_flamelet_library.cpp  -o $(OBJ_PATH)/sources_in_flamelet_library.o

$(OBJ_PATH)/OpenSMOKE_LookUp_Table_Manager.o: ../../src/OpenSMOKE_LookUp_Table_Manager.cpp ../../src/OpenSMOKE_LookUp_Table_Manager.h
					$(CPP) ../../src/OpenSMOKE_LookUp_Table_Manager.cpp  -o $(OBJ_PATH)/OpenSMOKE_LookUp_Table_Manager.o

$(OBJ_PATH)/OpenSMOKE_LookUp_Table_Executables.o: ../../src/OpenSMOKE_LookUp_Table_Executables.cpp ../../src/OpenSMOKE_LookUp_Table_Executables.h
					$(CPP) ../../src/OpenSMOKE_LookUp_Table_Executables.cpp  -o $(OBJ_PATH)/OpenSMOKE_LookUp_Table_Executables.o

$(OBJ_PATH)/flamelet_group.o: ../../src/flamelet_group.cpp ../../src/flamelet_group.h
					$(CPP) ../../src/flamelet_group.cpp  -o $(OBJ_PATH)/flamelet_group.o
	
$(OBJ_PATH)/flamelet.o: ../../src/flamelet.cpp ../../src/flamelet.h
					$(CPP) ../../src/flamelet.cpp  -o $(OBJ_PATH)/flamelet.o

OpenSMOKE_LookUpTables:  $(OBJS)  
			$(LINKER) $(OBJS) -o ../../exe/OpenSMOKE_LookUpTables_$(STRING_NAME).sh $(LIBS)

clean: ; @rm $(OBJ_PATH)/*.o ../../exe/$(FILE_NAME)_$(STRING_NAME).sh
