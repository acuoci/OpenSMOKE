#COMPILER
PATH_OBJECTS 	= OpenSMOKE_ICEM/$(OBJ_PATH)
PATH_CPP 	= ../../cpp/idealreactors/icem
PATH_HPP 	= ../../hpp/idealreactors/icem

# OBJECT FILES
OBJS =  $(PATH_OBJECTS)/OpenSMOKE_ICEM.o $(PATH_OBJECTS)/OpenSMOKE_ICEM_MultiZone.o

# MAIN TARGET
../../lib/linux/$(OBJ_PATH)/libOpenSMOKE_ICEM_$(STRING_NAME).a : $(OBJS)
	ar rc ../../lib/linux/$(OBJ_PATH)/libOpenSMOKE_ICEM_$(STRING_NAME).a $(OBJS)


# MAKE
$(PATH_OBJECTS)/OpenSMOKE_ICEM.o : $(PATH_CPP)/OpenSMOKE_ICEM.cpp $(PATH_HPP)/OpenSMOKE_ICEM.h
		 	 $(CCP) $(PATH_CPP)/OpenSMOKE_ICEM.cpp -o $(PATH_OBJECTS)/OpenSMOKE_ICEM.o

$(PATH_OBJECTS)/OpenSMOKE_ICEM_MultiZone.o : $(PATH_CPP)/OpenSMOKE_ICEM_MultiZone.cpp $(PATH_HPP)/OpenSMOKE_ICEM_MultiZone.h
		 	 $(CCP) $(PATH_CPP)/OpenSMOKE_ICEM_MultiZone.cpp -o $(PATH_OBJECTS)/OpenSMOKE_ICEM_MultiZone.o
			 	 
clean:
	-rm $(PATH_OBJECTS)/*.o ../../lib/linux/$(OBJ_PATH)/libOpenSMOKE_ICEM_$(STRING_NAME).a

