#COMPILER
PATH_OBJECTS 	= OpenSMOKE_QMOM/$(OBJ_PATH)
PATH_CPP 	= ../../cpp/qmom
PATH_HPP 	= ../../hpp/qmom

# OBJECT FILES
OBJS =  $(PATH_OBJECTS)/OpenSMOKE_Distributions.o $(PATH_OBJECTS)/OpenSMOKE_DQMOM.o $(PATH_OBJECTS)/OpenSMOKE_PDgordon.o $(PATH_OBJECTS)/OpenSMOKE_PhysicalModels.o $(PATH_OBJECTS)/OpenSMOKE_QMOM.o $(PATH_OBJECTS)/OpenSMOKE_QMOM_Module.o $(PATH_OBJECTS)/OpenSMOKE_Soot_Models.o 

# MAIN TARGET
../../lib/linux/$(OBJ_PATH)/libOpenSMOKE_QMOM_$(STRING_NAME).a : $(OBJS)
	ar rc ../../lib/linux/$(OBJ_PATH)/libOpenSMOKE_QMOM_$(STRING_NAME).a $(OBJS)


# MAKE
$(PATH_OBJECTS)/OpenSMOKE_Distributions.o : $(PATH_CPP)/OpenSMOKE_Distributions.cpp $(PATH_HPP)/OpenSMOKE_Distributions.h
		 	 $(CCP) $(PATH_CPP)/OpenSMOKE_Distributions.cpp -o $(PATH_OBJECTS)/OpenSMOKE_Distributions.o

$(PATH_OBJECTS)/OpenSMOKE_DQMOM.o : $(PATH_CPP)/OpenSMOKE_DQMOM.cpp $(PATH_HPP)/OpenSMOKE_DQMOM.h
		 	 $(CCP) $(PATH_CPP)/OpenSMOKE_DQMOM.cpp -o $(PATH_OBJECTS)/OpenSMOKE_DQMOM.o

$(PATH_OBJECTS)/OpenSMOKE_PDgordon.o : $(PATH_CPP)/OpenSMOKE_PDgordon.cpp $(PATH_HPP)/OpenSMOKE_PDgordon.h
		 	 $(CCP) $(PATH_CPP)/OpenSMOKE_PDgordon.cpp -o $(PATH_OBJECTS)/OpenSMOKE_PDgordon.o

$(PATH_OBJECTS)/OpenSMOKE_PhysicalModels.o : $(PATH_CPP)/OpenSMOKE_PhysicalModels.cpp $(PATH_HPP)/OpenSMOKE_PhysicalModels.h
		 	 $(CCP) $(PATH_CPP)/OpenSMOKE_PhysicalModels.cpp -o $(PATH_OBJECTS)/OpenSMOKE_PhysicalModels.o

$(PATH_OBJECTS)/OpenSMOKE_QMOM.o : $(PATH_CPP)/OpenSMOKE_QMOM.cpp $(PATH_HPP)/OpenSMOKE_QMOM.h
		 	 $(CCP) $(PATH_CPP)/OpenSMOKE_QMOM.cpp -o $(PATH_OBJECTS)/OpenSMOKE_QMOM.o

$(PATH_OBJECTS)/OpenSMOKE_QMOM_Module.o : $(PATH_CPP)/OpenSMOKE_QMOM_Module.cpp $(PATH_HPP)/OpenSMOKE_QMOM_Module.h
		 	 $(CCP) $(PATH_CPP)/OpenSMOKE_QMOM_Module.cpp -o $(PATH_OBJECTS)/OpenSMOKE_QMOM_Module.o

$(PATH_OBJECTS)/OpenSMOKE_Soot_Models.o : $(PATH_CPP)/OpenSMOKE_Soot_Models.cpp $(PATH_HPP)/OpenSMOKE_Soot_Models.h
		 	 $(CCP) $(PATH_CPP)/OpenSMOKE_Soot_Models.cpp -o $(PATH_OBJECTS)/OpenSMOKE_Soot_Models.o
			 	 
clean:
	-rm $(PATH_OBJECTS)/*.o ../../lib/linux/$(OBJ_PATH)/libOpenSMOKE_QMOM_$(STRING_NAME).a

