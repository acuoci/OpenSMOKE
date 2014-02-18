#COMPILER
PATH_OBJECTS 	= OpenSMOKE_KinPP/$(OBJ_PATH)
PATH_CPP 	= ../../cpp/kinpp
PATH_HPP 	= ../../hpp/kinpp

# OBJECT FILES
OBJS =  $(PATH_OBJECTS)/OpenSMOKE_CSTRNetwork.o \
	$(PATH_OBJECTS)/OpenSMOKE_CSTRNetwork_MyOdeSystemObjectAllCSTR.o \
	$(PATH_OBJECTS)/OpenSMOKE_CSTRNetwork_MyOdeSystemObjectOneCSTR.o \
	$(PATH_OBJECTS)/OpenSMOKE_DirectLinearSolver_Unsymmetric.o \
	$(PATH_OBJECTS)/OpenSMOKE_PARDISO_Unsymmetric.o 

# MAIN TARGET
../../lib/linux/$(OBJ_PATH)/libOpenSMOKE_KinPP_$(STRING_NAME).a : $(OBJS)
	ar rc ../../lib/linux/$(OBJ_PATH)/libOpenSMOKE_KinPP_$(STRING_NAME).a $(OBJS)


# MAKE
$(PATH_OBJECTS)/OpenSMOKE_CSTRNetwork.o : $(PATH_CPP)/OpenSMOKE_CSTRNetwork.cpp $(PATH_HPP)/OpenSMOKE_CSTRNetwork.h
		 	 $(CCP) $(PATH_CPP)/OpenSMOKE_CSTRNetwork.cpp -o $(PATH_OBJECTS)/OpenSMOKE_CSTRNetwork.o

$(PATH_OBJECTS)/OpenSMOKE_CSTRNetwork_MyOdeSystemObjectAllCSTR.o : $(PATH_CPP)/OpenSMOKE_CSTRNetwork_MyOdeSystemObjectAllCSTR.cpp $(PATH_HPP)/OpenSMOKE_CSTRNetwork_MyOdeSystemObjectAllCSTR.h
		 	 $(CCP) $(PATH_CPP)/OpenSMOKE_CSTRNetwork_MyOdeSystemObjectAllCSTR.cpp -o $(PATH_OBJECTS)/OpenSMOKE_CSTRNetwork_MyOdeSystemObjectAllCSTR.o

$(PATH_OBJECTS)/OpenSMOKE_CSTRNetwork_MyOdeSystemObjectOneCSTR.o : $(PATH_CPP)/OpenSMOKE_CSTRNetwork_MyOdeSystemObjectOneCSTR.cpp $(PATH_HPP)/OpenSMOKE_CSTRNetwork_MyOdeSystemObjectOneCSTR.h
		 	 $(CCP) $(PATH_CPP)/OpenSMOKE_CSTRNetwork_MyOdeSystemObjectOneCSTR.cpp -o $(PATH_OBJECTS)/OpenSMOKE_CSTRNetwork_MyOdeSystemObjectOneCSTR.o

$(PATH_OBJECTS)/OpenSMOKE_DirectLinearSolver_Unsymmetric.o : $(PATH_CPP)/OpenSMOKE_DirectLinearSolver_Unsymmetric.cpp $(PATH_HPP)/OpenSMOKE_DirectLinearSolver_Unsymmetric.h
		 	 $(CCP) $(PATH_CPP)/OpenSMOKE_DirectLinearSolver_Unsymmetric.cpp -o $(PATH_OBJECTS)/OpenSMOKE_DirectLinearSolver_Unsymmetric.o

$(PATH_OBJECTS)/OpenSMOKE_PARDISO_Unsymmetric.o : $(PATH_CPP)/OpenSMOKE_PARDISO_Unsymmetric.cpp $(PATH_HPP)/OpenSMOKE_PARDISO_Unsymmetric.h
		 	 $(CCP) $(PATH_CPP)/OpenSMOKE_PARDISO_Unsymmetric.cpp -o $(PATH_OBJECTS)/OpenSMOKE_PARDISO_Unsymmetric.o
			 	 
clean:
	-rm $(PATH_OBJECTS)/*.o ../../lib/linux/$(OBJ_PATH)/libOpenSMOKE_KinPP_$(STRING_NAME).a

