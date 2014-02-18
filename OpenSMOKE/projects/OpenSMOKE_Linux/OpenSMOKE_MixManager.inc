#COMPILER
PATH_OBJECTS 	= OpenSMOKE_MixManager/$(OBJ_PATH)
PATH_CPP 	= ../../cpp/mixmanager
PATH_HPP 	= ../../hpp/mixmanager

# OBJECT FILES
OBJS =  $(PATH_OBJECTS)/OpenSMOKE_MixManager.o

# MAIN TARGET
../../lib/linux/$(OBJ_PATH)/libOpenSMOKE_MixManager_$(STRING_NAME).a : $(OBJS)
	ar rc ../../lib/linux/$(OBJ_PATH)/libOpenSMOKE_MixManager_$(STRING_NAME).a $(OBJS)


# CHEMISTRY
$(PATH_OBJECTS)/OpenSMOKE_MixManager.o : $(PATH_CPP)/OpenSMOKE_MixManager.cpp $(PATH_HPP)/OpenSMOKE_MixManager.h
		 	 $(CCP) $(PATH_CPP)/OpenSMOKE_MixManager.cpp -o $(PATH_OBJECTS)/OpenSMOKE_MixManager.o
			 	 
clean:
	-rm $(PATH_OBJECTS)/*.o ../../lib/linux/$(OBJ_PATH)/libOpenSMOKE_MixManager_$(STRING_NAME).a

