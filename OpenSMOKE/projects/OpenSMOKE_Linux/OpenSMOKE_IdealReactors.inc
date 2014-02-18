#COMPILER
PATH_OBJECTS 	= OpenSMOKE_IdealReactors/$(OBJ_PATH)
PATH_CPP 	= ../../cpp/idealreactors
PATH_HPP 	= ../../hpp/idealreactors

# OBJECT FILES
OBJS =  $(PATH_OBJECTS)/OpenSMOKE_0DReactor.o $(PATH_OBJECTS)/OpenSMOKE_GasStream.o $(PATH_OBJECTS)/OpenSMOKE_UD_Profile.o 

# MAIN TARGET
../../lib/linux/$(OBJ_PATH)/libOpenSMOKE_IdealReactors_$(STRING_NAME).a : $(OBJS)
	ar rc ../../lib/linux/$(OBJ_PATH)/libOpenSMOKE_IdealReactors_$(STRING_NAME).a $(OBJS)


# MAKE
$(PATH_OBJECTS)/OpenSMOKE_0DReactor.o : $(PATH_CPP)/OpenSMOKE_0DReactor.cpp $(PATH_HPP)/OpenSMOKE_0DReactor.h
		 	 $(CCP) $(PATH_CPP)/OpenSMOKE_0DReactor.cpp -o $(PATH_OBJECTS)/OpenSMOKE_0DReactor.o

$(PATH_OBJECTS)/OpenSMOKE_GasStream.o : $(PATH_CPP)/OpenSMOKE_GasStream.cpp $(PATH_HPP)/OpenSMOKE_GasStream.h
		 	 $(CCP) $(PATH_CPP)/OpenSMOKE_GasStream.cpp -o $(PATH_OBJECTS)/OpenSMOKE_GasStream.o

$(PATH_OBJECTS)/OpenSMOKE_UD_Profile.o : $(PATH_CPP)/OpenSMOKE_UD_Profile.cpp $(PATH_HPP)/OpenSMOKE_UD_Profile.h
		 	 $(CCP) $(PATH_CPP)/OpenSMOKE_UD_Profile.cpp -o $(PATH_OBJECTS)/OpenSMOKE_UD_Profile.o
			 	 
clean:
	-rm $(PATH_OBJECTS)/*.o ../../lib/linux/$(OBJ_PATH)/libOpenSMOKE_IdealReactors_$(STRING_NAME).a

