#COMPILER
PATH_OBJECTS 	= OpenSMOKE_BatchReactor/$(OBJ_PATH)
PATH_CPP 	= ../../cpp/idealreactors/batch
PATH_HPP 	= ../../hpp/idealreactors/batch

# OBJECT FILES
OBJS =  $(PATH_OBJECTS)/OpenSMOKE_Batch.o

# MAIN TARGET
../../lib/linux/$(OBJ_PATH)/libOpenSMOKE_BatchReactor_$(STRING_NAME).a : $(OBJS)
	ar rc ../../lib/linux/$(OBJ_PATH)/libOpenSMOKE_BatchReactor_$(STRING_NAME).a $(OBJS)


# CHEMISTRY
$(PATH_OBJECTS)/OpenSMOKE_Batch.o : $(PATH_CPP)/OpenSMOKE_Batch.cpp $(PATH_HPP)/OpenSMOKE_Batch.h
		 	 $(CCP) $(PATH_CPP)/OpenSMOKE_Batch.cpp -o $(PATH_OBJECTS)/OpenSMOKE_Batch.o
			 	 
clean:
	-rm $(PATH_OBJECTS)/*.o ../../lib/linux/$(OBJ_PATH)/libOpenSMOKE_BatchReactor_$(STRING_NAME).a

