#COMPILER
PATH_OBJECTS 	= OpenSMOKE_PreProcessing/$(OBJ_PATH)
PATH_CPP 	= ../../cpp/preprocessing
PATH_HPP 	= ../../hpp/preprocessing

# OBJECT FILES
OBJS =  $(PATH_OBJECTS)/OpenSMOKE_PreProcessorIdealGas.o $(PATH_OBJECTS)/OpenSMOKE_PreProcessorKinetics.o $(PATH_OBJECTS)/OpenSMOKE_PreProcessorReactingGas.o 

# MAIN TARGET
../../lib/linux/$(OBJ_PATH)/libOpenSMOKE_PreProcessing_$(STRING_NAME).a : $(OBJS)
	ar rc ../../lib/linux/$(OBJ_PATH)/libOpenSMOKE_PreProcessing_$(STRING_NAME).a $(OBJS)


# MAKE
$(PATH_OBJECTS)/OpenSMOKE_PreProcessorIdealGas.o : $(PATH_CPP)/OpenSMOKE_PreProcessorIdealGas.cpp $(PATH_HPP)/OpenSMOKE_PreProcessorIdealGas.h
		 	 $(CCP) $(PATH_CPP)/OpenSMOKE_PreProcessorIdealGas.cpp -o $(PATH_OBJECTS)/OpenSMOKE_PreProcessorIdealGas.o

$(PATH_OBJECTS)/OpenSMOKE_PreProcessorKinetics.o : $(PATH_CPP)/OpenSMOKE_PreProcessorKinetics.cpp $(PATH_HPP)/OpenSMOKE_PreProcessorKinetics.h
		 	 $(CCP) $(PATH_CPP)/OpenSMOKE_PreProcessorKinetics.cpp -o $(PATH_OBJECTS)/OpenSMOKE_PreProcessorKinetics.o

$(PATH_OBJECTS)/OpenSMOKE_PreProcessorReactingGas.o : $(PATH_CPP)/OpenSMOKE_PreProcessorReactingGas.cpp $(PATH_HPP)/OpenSMOKE_PreProcessorReactingGas.h
		 	 $(CCP) $(PATH_CPP)/OpenSMOKE_PreProcessorReactingGas.cpp -o $(PATH_OBJECTS)/OpenSMOKE_PreProcessorReactingGas.o
			 	 
clean:
	-rm $(PATH_OBJECTS)/*.o ../../lib/linux/$(OBJ_PATH)/libOpenSMOKE_PreProcessing_$(STRING_NAME).a

