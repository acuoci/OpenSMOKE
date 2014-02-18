#COMPILER 
PATH_OBJECTS 	= OpenSMOKE_Engine/$(OBJ_PATH)
PATH_CPP 	= ../../cpp/engine
PATH_HPP 	= ../../hpp/engine

# OBJECT FILES
OBJS =  $(PATH_OBJECTS)/OpenSMOKE_EquilibriumStanjan.o $(PATH_OBJECTS)/OpenSMOKE_ReactingGas.o $(PATH_OBJECTS)/OpenSMOKE_GlobalKinetics.o $(PATH_OBJECTS)/OpenSMOKE_Kinetics.o $(PATH_OBJECTS)/OpenSMOKE_IdealGas.o 

# MAIN TARGET
../../lib/linux/$(OBJ_PATH)/libOpenSMOKE_Engine_$(STRING_NAME).a : $(OBJS)
	ar rc ../../lib/linux/$(OBJ_PATH)/libOpenSMOKE_Engine_$(STRING_NAME).a $(OBJS)


# MAKE

$(PATH_OBJECTS)/OpenSMOKE_EquilibriumStanjan.o : $(PATH_CPP)/OpenSMOKE_EquilibriumStanjan.cpp $(PATH_HPP)/OpenSMOKE_EquilibriumStanjan.h
		 	 $(CCP) $(PATH_CPP)/OpenSMOKE_EquilibriumStanjan.cpp -o $(PATH_OBJECTS)/OpenSMOKE_EquilibriumStanjan.o

$(PATH_OBJECTS)/OpenSMOKE_ReactingGas.o : $(PATH_CPP)/OpenSMOKE_ReactingGas.cpp $(PATH_HPP)/OpenSMOKE_ReactingGas.h
		 	 $(CCP) $(PATH_CPP)/OpenSMOKE_ReactingGas.cpp -o $(PATH_OBJECTS)/OpenSMOKE_ReactingGas.o

$(PATH_OBJECTS)/OpenSMOKE_GlobalKinetics.o : $(PATH_CPP)/OpenSMOKE_GlobalKinetics.cpp $(PATH_HPP)/OpenSMOKE_GlobalKinetics.h
		 	 $(CCP) $(PATH_CPP)/OpenSMOKE_GlobalKinetics.cpp -o $(PATH_OBJECTS)/OpenSMOKE_GlobalKinetics.o

$(PATH_OBJECTS)/OpenSMOKE_Kinetics.o : $(PATH_CPP)/OpenSMOKE_Kinetics.cpp $(PATH_HPP)/OpenSMOKE_Kinetics.h
		 	 $(CCP) $(PATH_CPP)/OpenSMOKE_Kinetics.cpp -o $(PATH_OBJECTS)/OpenSMOKE_Kinetics.o

$(PATH_OBJECTS)/OpenSMOKE_IdealGas.o : $(PATH_CPP)/OpenSMOKE_IdealGas.cpp $(PATH_HPP)/OpenSMOKE_IdealGas.h
		 	 $(CCP) $(PATH_CPP)/OpenSMOKE_IdealGas.cpp -o $(PATH_OBJECTS)/OpenSMOKE_IdealGas.o
			 	 
clean:
	-rm $(PATH_OBJECTS)/*.o ../../lib/linux/$(OBJ_PATH)/libOpenSMOKE_Engine_$(STRING_NAME).a

