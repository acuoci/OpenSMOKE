#COMPILER
PATH_OBJECTS 	= OpenSMOKE_PlugFlowReactor/$(OBJ_PATH)
PATH_CPP 	= ../../cpp/idealreactors/pfr
PATH_HPP 	= ../../hpp/idealreactors/pfr

# OBJECT FILES
OBJS =  $(PATH_OBJECTS)/OpenSMOKE_PFR.o \
	$(PATH_OBJECTS)/OpenSMOKE_PFR_Geometry.o

# MAIN TARGET
../../lib/linux/$(OBJ_PATH)/libOpenSMOKE_PlugFlowReactor_$(STRING_NAME).a : $(OBJS)
	ar rc ../../lib/linux/$(OBJ_PATH)/libOpenSMOKE_PlugFlowReactor_$(STRING_NAME).a $(OBJS)


# MAKE
$(PATH_OBJECTS)/OpenSMOKE_PFR.o : $(PATH_CPP)/OpenSMOKE_PFR.cpp $(PATH_HPP)/OpenSMOKE_PFR.h
		 	 $(CCP) $(PATH_CPP)/OpenSMOKE_PFR.cpp -o $(PATH_OBJECTS)/OpenSMOKE_PFR.o

$(PATH_OBJECTS)/OpenSMOKE_PFR_Geometry.o : $(PATH_CPP)/OpenSMOKE_PFR_Geometry.cpp $(PATH_HPP)/OpenSMOKE_PFR_Geometry.h
		 	 $(CCP) $(PATH_CPP)/OpenSMOKE_PFR_Geometry.cpp -o $(PATH_OBJECTS)/OpenSMOKE_PFR_Geometry.o

			 	 
clean:
	-rm $(PATH_OBJECTS)/*.o ../../lib/linux/$(OBJ_PATH)/libOpenSMOKE_PlugFlowReactor_$(STRING_NAME).a

