#COMPILER
PATH_OBJECTS 	= OpenSMOKE_Interfaces/$(OBJ_PATH)
PATH_CPP 	= ../../cpp/interfaces
PATH_HPP 	= ../../hpp/interfaces

# OBJECT FILES
OBJS =  $(PATH_OBJECTS)/OpenSMOKE_GnuPlotInterface.o $(PATH_OBJECTS)/OpenSMOKE_GnuPlotWrap.o $(PATH_OBJECTS)/OpenSMOKE_LatexInterface.o 

# MAIN TARGET
../../lib/linux/$(OBJ_PATH)/libOpenSMOKE_Interfaces_$(STRING_NAME).a : $(OBJS)
	ar rc ../../lib/linux/$(OBJ_PATH)/libOpenSMOKE_Interfaces_$(STRING_NAME).a $(OBJS)


# MAKE
$(PATH_OBJECTS)/OpenSMOKE_GnuPlotInterface.o : $(PATH_CPP)/OpenSMOKE_GnuPlotInterface.cpp $(PATH_HPP)/OpenSMOKE_GnuPlotInterface.h
		 	 $(CCP) $(PATH_CPP)/OpenSMOKE_GnuPlotInterface.cpp -o $(PATH_OBJECTS)/OpenSMOKE_GnuPlotInterface.o

$(PATH_OBJECTS)/OpenSMOKE_GnuPlotWrap.o : $(PATH_CPP)/OpenSMOKE_GnuPlotWrap.cpp $(PATH_HPP)/OpenSMOKE_GnuPlotWrap.h
		 	 $(CCP) $(PATH_CPP)/OpenSMOKE_GnuPlotWrap.cpp -o $(PATH_OBJECTS)/OpenSMOKE_GnuPlotWrap.o

$(PATH_OBJECTS)/OpenSMOKE_LatexInterface.o : $(PATH_CPP)/OpenSMOKE_LatexInterface.cpp $(PATH_HPP)/OpenSMOKE_LatexInterface.h
		 	 $(CCP) $(PATH_CPP)/OpenSMOKE_LatexInterface.cpp -o $(PATH_OBJECTS)/OpenSMOKE_LatexInterface.o
			 	 
clean:
	-rm $(PATH_OBJECTS)/*.o ../../lib/linux/$(OBJ_PATH)/libOpenSMOKE_Interfaces_$(STRING_NAME).a

