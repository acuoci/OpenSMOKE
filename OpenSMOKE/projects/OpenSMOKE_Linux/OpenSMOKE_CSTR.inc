#COMPILER
PATH_OBJECTS 	= OpenSMOKE_CSTR/$(OBJ_PATH)
PATH_CPP 	= ../../cpp/idealreactors/cstr
PATH_HPP 	= ../../hpp/idealreactors/cstr

# OBJECT FILES
OBJS =  $(PATH_OBJECTS)/OpenSMOKE_CSTR.o $(PATH_OBJECTS)/OpenSMOKE_CSTR_Geometry.o

# MAIN TARGET
../../lib/linux/$(OBJ_PATH)/libOpenSMOKE_CSTR_$(STRING_NAME).a : $(OBJS)
	ar rc ../../lib/linux/$(OBJ_PATH)/libOpenSMOKE_CSTR_$(STRING_NAME).a $(OBJS)


# MAKE
$(PATH_OBJECTS)/OpenSMOKE_CSTR.o : $(PATH_CPP)/OpenSMOKE_CSTR.cpp $(PATH_HPP)/OpenSMOKE_CSTR.h
		 	 $(CCP) $(PATH_CPP)/OpenSMOKE_CSTR.cpp -o $(PATH_OBJECTS)/OpenSMOKE_CSTR.o

$(PATH_OBJECTS)/OpenSMOKE_CSTR_Geometry.o : $(PATH_CPP)/OpenSMOKE_CSTR_Geometry.cpp $(PATH_HPP)/OpenSMOKE_CSTR_Geometry.h
		 	 $(CCP) $(PATH_CPP)/OpenSMOKE_CSTR_Geometry.cpp -o $(PATH_OBJECTS)/OpenSMOKE_CSTR_Geometry.o

			 	 
clean:
	-rm $(PATH_OBJECTS)/*.o ../../lib/linux/$(OBJ_PATH)/libOpenSMOKE_CSTR_$(STRING_NAME).a

