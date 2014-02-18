#COMPILER
PATH_OBJECTS 	= OpenSMOKE_LiquidProperties/$(OBJ_PATH)
PATH_CPP 	= ../../cpp/liquid
PATH_HPP 	= ../../hpp/liquid


# OBJECT FILES
OBJS =  $(PATH_OBJECTS)/OpenSMOKE_LiquidCriticalConstants_Dictionary.o \
	$(PATH_OBJECTS)/OpenSMOKE_LiquidDensity_Dictionary.o \
	$(PATH_OBJECTS)/OpenSMOKE_LiquidProperties_Database.o \
	$(PATH_OBJECTS)/OpenSMOKE_LiquidSpecies.o \
	$(PATH_OBJECTS)/OpenSMOKE_LiquidSpecificHeat_Dictionary.o \
	$(PATH_OBJECTS)/OpenSMOKE_LiquidThermalConductivity_Dictionary.o \
	$(PATH_OBJECTS)/OpenSMOKE_LiquidVaporizationHeat_Dictionary.o \
	$(PATH_OBJECTS)/OpenSMOKE_LiquidVaporPressure_Dictionary.o 

# MAIN TARGET
../../lib/linux/$(OBJ_PATH)/libOpenSMOKE_LiquidProperties_$(STRING_NAME).a : $(OBJS)
	ar rc ../../lib/linux/$(OBJ_PATH)/libOpenSMOKE_LiquidProperties_$(STRING_NAME).a $(OBJS)


# CHEMISTRY

#1
$(PATH_OBJECTS)/OpenSMOKE_LiquidCriticalConstants_Dictionary.o : $(PATH_CPP)/OpenSMOKE_LiquidCriticalConstants_Dictionary.cpp $(PATH_HPP)/OpenSMOKE_LiquidCriticalConstants_Dictionary.h
		 	 $(CCP) $(PATH_CPP)/OpenSMOKE_LiquidCriticalConstants_Dictionary.cpp -o $(PATH_OBJECTS)/OpenSMOKE_LiquidCriticalConstants_Dictionary.o

#2
$(PATH_OBJECTS)/OpenSMOKE_LiquidDensity_Dictionary.o : $(PATH_CPP)/OpenSMOKE_LiquidDensity_Dictionary.cpp $(PATH_HPP)/OpenSMOKE_LiquidDensity_Dictionary.h
		 	 $(CCP) $(PATH_CPP)/OpenSMOKE_LiquidDensity_Dictionary.cpp -o $(PATH_OBJECTS)/OpenSMOKE_LiquidDensity_Dictionary.o

#3
$(PATH_OBJECTS)/OpenSMOKE_LiquidProperties_Database.o : $(PATH_CPP)/OpenSMOKE_LiquidProperties_Database.cpp $(PATH_HPP)/OpenSMOKE_LiquidProperties_Database.h
		 	 $(CCP) $(PATH_CPP)/OpenSMOKE_LiquidProperties_Database.cpp -o $(PATH_OBJECTS)/OpenSMOKE_LiquidProperties_Database.o

#4
$(PATH_OBJECTS)/OpenSMOKE_LiquidSpecies.o : $(PATH_CPP)/OpenSMOKE_LiquidSpecies.cpp $(PATH_HPP)/OpenSMOKE_LiquidSpecies.h
		 	 $(CCP) $(PATH_CPP)/OpenSMOKE_LiquidSpecies.cpp -o $(PATH_OBJECTS)/OpenSMOKE_LiquidSpecies.o

#5
$(PATH_OBJECTS)/OpenSMOKE_LiquidSpecificHeat_Dictionary.o : $(PATH_CPP)/OpenSMOKE_LiquidSpecificHeat_Dictionary.cpp $(PATH_HPP)/OpenSMOKE_LiquidSpecificHeat_Dictionary.h
		 	 $(CCP) $(PATH_CPP)/OpenSMOKE_LiquidSpecificHeat_Dictionary.cpp -o $(PATH_OBJECTS)/OpenSMOKE_LiquidSpecificHeat_Dictionary.o

#6
$(PATH_OBJECTS)/OpenSMOKE_LiquidThermalConductivity_Dictionary.o : $(PATH_CPP)/OpenSMOKE_LiquidThermalConductivity_Dictionary.cpp $(PATH_HPP)/OpenSMOKE_LiquidThermalConductivity_Dictionary.h
		 	 $(CCP) $(PATH_CPP)/OpenSMOKE_LiquidThermalConductivity_Dictionary.cpp -o $(PATH_OBJECTS)/OpenSMOKE_LiquidThermalConductivity_Dictionary.o

#7
$(PATH_OBJECTS)/OpenSMOKE_LiquidVaporizationHeat_Dictionary.o : $(PATH_CPP)/OpenSMOKE_LiquidVaporizationHeat_Dictionary.cpp $(PATH_HPP)/OpenSMOKE_LiquidVaporizationHeat_Dictionary.h
		 	 $(CCP) $(PATH_CPP)/OpenSMOKE_LiquidVaporizationHeat_Dictionary.cpp -o $(PATH_OBJECTS)/OpenSMOKE_LiquidVaporizationHeat_Dictionary.o

#8
$(PATH_OBJECTS)/OpenSMOKE_LiquidVaporPressure_Dictionary.o : $(PATH_CPP)/OpenSMOKE_LiquidVaporPressure_Dictionary.cpp $(PATH_HPP)/OpenSMOKE_LiquidVaporPressure_Dictionary.h
		 	 $(CCP) $(PATH_CPP)/OpenSMOKE_LiquidVaporPressure_Dictionary.cpp -o $(PATH_OBJECTS)/OpenSMOKE_LiquidVaporPressure_Dictionary.o
		 	 
clean:
	-rm $(PATH_OBJECTS)/*.o ../../lib/linux/$(OBJ_PATH)/libOpenSMOKE_LiquidProperties_$(STRING_NAME).a

