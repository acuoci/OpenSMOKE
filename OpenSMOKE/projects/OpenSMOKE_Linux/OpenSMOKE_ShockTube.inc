#COMPILER
PATH_OBJECTS 	= OpenSMOKE_ShockTube/$(OBJ_PATH)
PATH_CPP 	= ../../cpp/idealreactors/shocktube
PATH_HPP 	= ../../hpp/idealreactors/shocktube

# OBJECT FILES
OBJS =  $(PATH_OBJECTS)/OpenSMOKE_ShockTube.o \
	$(PATH_OBJECTS)/OpenSMOKE_ShockTube_InitialConditions.o

# MAIN TARGET
../../lib/linux/$(OBJ_PATH)/libOpenSMOKE_ShockTube_$(STRING_NAME).a : $(OBJS)
	ar rc ../../lib/linux/$(OBJ_PATH)/libOpenSMOKE_ShockTube_$(STRING_NAME).a $(OBJS)


# MAKE
$(PATH_OBJECTS)/OpenSMOKE_ShockTube.o : $(PATH_CPP)/OpenSMOKE_ShockTube.cpp $(PATH_HPP)/OpenSMOKE_ShockTube.h
		 	 $(CCP) $(PATH_CPP)/OpenSMOKE_ShockTube.cpp -o $(PATH_OBJECTS)/OpenSMOKE_ShockTube.o

$(PATH_OBJECTS)/OpenSMOKE_ShockTube_InitialConditions.o : $(PATH_CPP)/OpenSMOKE_ShockTube_InitialConditions.cpp $(PATH_HPP)/OpenSMOKE_ShockTube_InitialConditions.h
		 	 $(CCP) $(PATH_CPP)/OpenSMOKE_ShockTube_InitialConditions.cpp -o $(PATH_OBJECTS)/OpenSMOKE_ShockTube_InitialConditions.o

			 	 
clean:
	-rm $(PATH_OBJECTS)/*.o ../../lib/linux/$(OBJ_PATH)/libOpenSMOKE_ShockTube_$(STRING_NAME).a
