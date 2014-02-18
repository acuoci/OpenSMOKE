#COMPILER 
PATH_OBJECTS 	= OpenSMOKE_Flamelet/$(OBJ_PATH)
PATH_CPP 	= ../../cpp/idealreactors/flamelet
PATH_HPP 	= ../../hpp/idealreactors/flamelet

# OBJECT FILES
OBJS =  $(PATH_OBJECTS)/OpenSMOKE_Flamelet.o \
	$(PATH_OBJECTS)/OpenSMOKE_Flamelet_DataManager.o \
	$(PATH_OBJECTS)/OpenSMOKE_Flamelet_ODE_Objects.o \
	$(PATH_OBJECTS)/OpenSMOKE_Flamelet_ScheduleClass.o

# MAIN TARGET
../../lib/linux/$(OBJ_PATH)/libOpenSMOKE_Flamelet_$(STRING_NAME).a : $(OBJS)
	ar rc ../../lib/linux/$(OBJ_PATH)/libOpenSMOKE_Flamelet_$(STRING_NAME).a $(OBJS)


# MAKE
$(PATH_OBJECTS)/OpenSMOKE_Flamelet.o : $(PATH_CPP)/OpenSMOKE_Flamelet.cpp $(PATH_HPP)/OpenSMOKE_Flamelet.h
		 	 $(CCP) $(PATH_CPP)/OpenSMOKE_Flamelet.cpp -o $(PATH_OBJECTS)/OpenSMOKE_Flamelet.o

$(PATH_OBJECTS)/OpenSMOKE_Flamelet_DataManager.o : $(PATH_CPP)/OpenSMOKE_Flamelet_DataManager.cpp $(PATH_HPP)/OpenSMOKE_Flamelet_DataManager.h
		 	 $(CCP) $(PATH_CPP)/OpenSMOKE_Flamelet_DataManager.cpp -o $(PATH_OBJECTS)/OpenSMOKE_Flamelet_DataManager.o

$(PATH_OBJECTS)/OpenSMOKE_Flamelet_ODE_Objects.o : $(PATH_CPP)/OpenSMOKE_Flamelet_ODE_Objects.cpp $(PATH_HPP)/OpenSMOKE_Flamelet_ODE_Objects.h
		 	 $(CCP) $(PATH_CPP)/OpenSMOKE_Flamelet_ODE_Objects.cpp -o $(PATH_OBJECTS)/OpenSMOKE_Flamelet_ODE_Objects.o

$(PATH_OBJECTS)/OpenSMOKE_Flamelet_ScheduleClass.o : $(PATH_CPP)/OpenSMOKE_Flamelet_ScheduleClass.cpp $(PATH_HPP)/OpenSMOKE_Flamelet_ScheduleClass.h
		 	 $(CCP) $(PATH_CPP)/OpenSMOKE_Flamelet_ScheduleClass.cpp -o $(PATH_OBJECTS)/OpenSMOKE_Flamelet_ScheduleClass.o
			 	 
clean:
	-rm $(PATH_OBJECTS)/*.o ../../lib/linux/$(OBJ_PATH)/libOpenSMOKE_Flamelet_$(STRING_NAME).a

