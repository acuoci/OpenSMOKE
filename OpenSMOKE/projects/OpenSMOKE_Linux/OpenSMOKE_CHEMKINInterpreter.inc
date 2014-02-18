#COMPILER
PATH_OBJECTS 	= OpenSMOKE_CHEMKINInterpreter/$(OBJ_PATH)
PATH_CPP 	= ../../cpp/preprocessing_chemkin
PATH_HPP 	= ../../hpp/preprocessing_chemkin

# OBJECT FILES
OBJS =  $(PATH_OBJECTS)/OpenSMOKE_CHEMKINInterpreter.o $(PATH_OBJECTS)/OpenSMOKE_CHEMKINInterpreter_ElementsData.o $(PATH_OBJECTS)/OpenSMOKE_CHEMKINInterpreter_KineticsData.o $(PATH_OBJECTS)/OpenSMOKE_CHEMKINInterpreter_ReactionData.o $(PATH_OBJECTS)/OpenSMOKE_CHEMKINInterpreter_SpeciesData.o $(PATH_OBJECTS)/OpenSMOKE_CHEMKINInterpreter_ThermoData.o $(PATH_OBJECTS)/OpenSMOKE_CHEMKINInterpreter_ThermoSpecies.o $(PATH_OBJECTS)/OpenSMOKE_CHEMKINInterpreter_TransportData.o  $(PATH_OBJECTS)/OpenSMOKE_CHEMKINInterpreter_TransportSpecies.o $(PATH_OBJECTS)/OpenSMOKE_CHEMKINInterpreter_UnitsData.o $(PATH_OBJECTS)/OpenSMOKE_CHEMKINInterpreter_SurfaceReactionData.o $(PATH_OBJECTS)/OpenSMOKE_CHEMKINInterpreter_SurfaceKineticsData.o 


# MAIN TARGET
../../lib/linux/$(OBJ_PATH)/libOpenSMOKE_CHEMKINInterpreter_$(STRING_NAME).a : $(OBJS)
	ar rc ../../lib/linux/$(OBJ_PATH)/libOpenSMOKE_CHEMKINInterpreter_$(STRING_NAME).a $(OBJS)


# MAKE

$(PATH_OBJECTS)/OpenSMOKE_CHEMKINInterpreter.o : $(PATH_CPP)/OpenSMOKE_CHEMKINInterpreter.cpp $(PATH_HPP)/OpenSMOKE_CHEMKINInterpreter.h
		 	 $(CCP) $(PATH_CPP)/OpenSMOKE_CHEMKINInterpreter.cpp -o $(PATH_OBJECTS)/OpenSMOKE_CHEMKINInterpreter.o

$(PATH_OBJECTS)/OpenSMOKE_CHEMKINInterpreter_ElementsData.o : $(PATH_CPP)/OpenSMOKE_CHEMKINInterpreter_ElementsData.cpp $(PATH_HPP)/OpenSMOKE_CHEMKINInterpreter_ElementsData.h
		 	 $(CCP) $(PATH_CPP)/OpenSMOKE_CHEMKINInterpreter_ElementsData.cpp -o $(PATH_OBJECTS)/OpenSMOKE_CHEMKINInterpreter_ElementsData.o

$(PATH_OBJECTS)/OpenSMOKE_CHEMKINInterpreter_KineticsData.o : $(PATH_CPP)/OpenSMOKE_CHEMKINInterpreter_KineticsData.cpp $(PATH_HPP)/OpenSMOKE_CHEMKINInterpreter_KineticsData.h
		 	 $(CCP) $(PATH_CPP)/OpenSMOKE_CHEMKINInterpreter_KineticsData.cpp -o $(PATH_OBJECTS)/OpenSMOKE_CHEMKINInterpreter_KineticsData.o

$(PATH_OBJECTS)/OpenSMOKE_CHEMKINInterpreter_ReactionData.o : $(PATH_CPP)/OpenSMOKE_CHEMKINInterpreter_ReactionData.cpp $(PATH_HPP)/OpenSMOKE_CHEMKINInterpreter_ReactionData.h
		 	 $(CCP) $(PATH_CPP)/OpenSMOKE_CHEMKINInterpreter_ReactionData.cpp -o $(PATH_OBJECTS)/OpenSMOKE_CHEMKINInterpreter_ReactionData.o

$(PATH_OBJECTS)/OpenSMOKE_CHEMKINInterpreter_SpeciesData.o : $(PATH_CPP)/OpenSMOKE_CHEMKINInterpreter_SpeciesData.cpp $(PATH_HPP)/OpenSMOKE_CHEMKINInterpreter_SpeciesData.h
		 	 $(CCP) $(PATH_CPP)/OpenSMOKE_CHEMKINInterpreter_SpeciesData.cpp -o $(PATH_OBJECTS)/OpenSMOKE_CHEMKINInterpreter_SpeciesData.o

$(PATH_OBJECTS)/OpenSMOKE_CHEMKINInterpreter_ThermoData.o : $(PATH_CPP)/OpenSMOKE_CHEMKINInterpreter_ThermoData.cpp $(PATH_HPP)/OpenSMOKE_CHEMKINInterpreter_ThermoData.h
		 	 $(CCP) $(PATH_CPP)/OpenSMOKE_CHEMKINInterpreter_ThermoData.cpp -o $(PATH_OBJECTS)/OpenSMOKE_CHEMKINInterpreter_ThermoData.o

$(PATH_OBJECTS)/OpenSMOKE_CHEMKINInterpreter_ThermoSpecies.o : $(PATH_CPP)/OpenSMOKE_CHEMKINInterpreter_ThermoSpecies.cpp $(PATH_HPP)/OpenSMOKE_CHEMKINInterpreter_ThermoSpecies.h
		 	 $(CCP) $(PATH_CPP)/OpenSMOKE_CHEMKINInterpreter_ThermoSpecies.cpp -o $(PATH_OBJECTS)/OpenSMOKE_CHEMKINInterpreter_ThermoSpecies.o

$(PATH_OBJECTS)/OpenSMOKE_CHEMKINInterpreter_TransportData.o : $(PATH_CPP)/OpenSMOKE_CHEMKINInterpreter_TransportData.cpp $(PATH_HPP)/OpenSMOKE_CHEMKINInterpreter_TransportData.h
		 	 $(CCP) $(PATH_CPP)/OpenSMOKE_CHEMKINInterpreter_TransportData.cpp -o $(PATH_OBJECTS)/OpenSMOKE_CHEMKINInterpreter_TransportData.o

$(PATH_OBJECTS)/OpenSMOKE_CHEMKINInterpreter_TransportSpecies.o : $(PATH_CPP)/OpenSMOKE_CHEMKINInterpreter_TransportSpecies.cpp $(PATH_HPP)/OpenSMOKE_CHEMKINInterpreter_TransportSpecies.h
		 	 $(CCP) $(PATH_CPP)/OpenSMOKE_CHEMKINInterpreter_TransportSpecies.cpp -o $(PATH_OBJECTS)/OpenSMOKE_CHEMKINInterpreter_TransportSpecies.o

$(PATH_OBJECTS)/OpenSMOKE_CHEMKINInterpreter_UnitsData.o : $(PATH_CPP)/OpenSMOKE_CHEMKINInterpreter_UnitsData.cpp $(PATH_HPP)/OpenSMOKE_CHEMKINInterpreter_UnitsData.h
		 	 $(CCP) $(PATH_CPP)/OpenSMOKE_CHEMKINInterpreter_UnitsData.cpp -o $(PATH_OBJECTS)/OpenSMOKE_CHEMKINInterpreter_UnitsData.o

$(PATH_OBJECTS)/OpenSMOKE_CHEMKINInterpreter_SurfaceReactionData.o : $(PATH_CPP)/OpenSMOKE_CHEMKINInterpreter_SurfaceReactionData.cpp $(PATH_HPP)/OpenSMOKE_CHEMKINInterpreter_SurfaceReactionData.h
		 	 $(CCP) $(PATH_CPP)/OpenSMOKE_CHEMKINInterpreter_SurfaceReactionData.cpp -o $(PATH_OBJECTS)/OpenSMOKE_CHEMKINInterpreter_SurfaceReactionData.o

$(PATH_OBJECTS)/OpenSMOKE_CHEMKINInterpreter_SurfaceKineticsData.o : $(PATH_CPP)/OpenSMOKE_CHEMKINInterpreter_SurfaceKineticsData.cpp $(PATH_HPP)/OpenSMOKE_CHEMKINInterpreter_SurfaceKineticsData.h
		 	 $(CCP) $(PATH_CPP)/OpenSMOKE_CHEMKINInterpreter_SurfaceKineticsData.cpp -o $(PATH_OBJECTS)/OpenSMOKE_CHEMKINInterpreter_SurfaceKineticsData.o

			 	 
clean:
	-rm $(PATH_OBJECTS)/*.o ../../lib/linux/$(OBJ_PATH)/libOpenSMOKE_CHEMKINInterpreter_$(STRING_NAME).a

