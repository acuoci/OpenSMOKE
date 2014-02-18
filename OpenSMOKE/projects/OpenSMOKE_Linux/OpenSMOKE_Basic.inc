# COMPILER
PATH_OBJECTS 	= OpenSMOKE_Basic/$(OBJ_PATH)
PATH_CPP 	= ../../cpp/basic
PATH_HPP 	= ../../hpp/basic

# OBJECT FILES
OBJS =  $(PATH_OBJECTS)/OpenSMOKE_Conversions.o $(PATH_OBJECTS)/OpenSMOKE_Grid1D.o $(PATH_OBJECTS)/OpenSMOKE_Utilities.o $(PATH_OBJECTS)/OpenSMOKE_AdaptiveGrid.o $(PATH_OBJECTS)/OpenSMOKE_Dictionary.o $(PATH_OBJECTS)/OpenSMOKE_LogFile.o  $(PATH_OBJECTS)/OpenSMOKE_WarningFile.o

# MAIN TARGET
../../lib/linux/$(OBJ_PATH)/libOpenSMOKE_Basic_$(STRING_NAME).a : $(OBJS)
	ar rc ../../lib/linux/$(OBJ_PATH)/libOpenSMOKE_Basic_$(STRING_NAME).a $(OBJS)


# CHEMISTRY
$(PATH_OBJECTS)/OpenSMOKE_Conversions.o : $(PATH_CPP)/OpenSMOKE_Conversions.cpp $(PATH_HPP)/OpenSMOKE_Conversions.h
		 	 $(CCP) $(PATH_CPP)/OpenSMOKE_Conversions.cpp -o $(PATH_OBJECTS)/OpenSMOKE_Conversions.o

$(PATH_OBJECTS)/OpenSMOKE_Grid1D.o : $(PATH_CPP)/OpenSMOKE_Grid1D.cpp $(PATH_HPP)/OpenSMOKE_Grid1D.h
		 	 $(CCP) $(PATH_CPP)/OpenSMOKE_Grid1D.cpp -o $(PATH_OBJECTS)/OpenSMOKE_Grid1D.o

$(PATH_OBJECTS)/OpenSMOKE_Utilities.o : $(PATH_CPP)/OpenSMOKE_Utilities.cpp $(PATH_HPP)/OpenSMOKE_Utilities.h
		 	 $(CCP) $(PATH_CPP)/OpenSMOKE_Utilities.cpp -o $(PATH_OBJECTS)/OpenSMOKE_Utilities.o

$(PATH_OBJECTS)/OpenSMOKE_Dictionary.o : $(PATH_CPP)/OpenSMOKE_Dictionary.cpp $(PATH_HPP)/OpenSMOKE_Dictionary.h
		 	 $(CCP) $(PATH_CPP)/OpenSMOKE_Dictionary.cpp -o $(PATH_OBJECTS)/OpenSMOKE_Dictionary.o

$(PATH_OBJECTS)/OpenSMOKE_AdaptiveGrid.o : $(PATH_CPP)/OpenSMOKE_AdaptiveGrid.cpp $(PATH_HPP)/OpenSMOKE_AdaptiveGrid.h
		 	 $(CCP) $(PATH_CPP)/OpenSMOKE_AdaptiveGrid.cpp -o $(PATH_OBJECTS)/OpenSMOKE_AdaptiveGrid.o

$(PATH_OBJECTS)/OpenSMOKE_LogFile.o : $(PATH_CPP)/OpenSMOKE_LogFile.cpp $(PATH_HPP)/OpenSMOKE_LogFile.h
		 	 $(CCP) $(PATH_CPP)/OpenSMOKE_LogFile.cpp -o $(PATH_OBJECTS)/OpenSMOKE_LogFile.o

$(PATH_OBJECTS)/OpenSMOKE_WarningFile.o : $(PATH_CPP)/OpenSMOKE_WarningFile.cpp $(PATH_HPP)/OpenSMOKE_WarningFile.h
		 	 $(CCP) $(PATH_CPP)/OpenSMOKE_WarningFile.cpp -o $(PATH_OBJECTS)/OpenSMOKE_WarningFile.o
			 	 
clean:
	-rm $(PATH_OBJECTS)/*.o ../../lib/linux/$(OBJ_PATH)/libOpenSMOKE_Basic_$(STRING_NAME).a

