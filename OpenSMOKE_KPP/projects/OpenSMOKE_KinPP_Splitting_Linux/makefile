include make.inc

# ------------------------------------------------------------------------------------------------
# OBJECTS
# ------------------------------------------------------------------------------------------------
OBJS     =  $(OBJ_PATH)/OpenSMOKE_KPP.o \
	    $(OBJ_PATH)/OpenSMOKE_KPP_Definitions.o \
            $(OBJ_PATH)/linear_solvers/OpenSMOKE_DirectLinearSolver_Unsymmetric.o \
            $(OBJ_PATH)/linear_solvers/OpenSMOKE_LAPACK_Dense.o \
            $(OBJ_PATH)/linear_solvers/OpenSMOKE_LIS_Unsymmetric.o \
            $(OBJ_PATH)/linear_solvers/OOpenSMOKE_PARDISO_Unsymmetric.o \
            $(OBJ_PATH)/kpp/OpenSMOKE_KinPP_SingleReactor.o \
            $(OBJ_PATH)/kpp/OpenSMOKE_KinPP_DataManager.o \
            $(OBJ_PATH)/kpp/OpenSMOKE_Dictionary_KinPP.o \
            $(OBJ_PATH)/kpp/OpenSMOKE_ReactorNetwork.o \
            $(OBJ_PATH)/kpp/OpenSMOKE_KinPP_ConvectiveNetworkStatistics.o \
            $(OBJ_PATH)/kpp/OpenSMOKE_KinPP_SingleReactorStatistics.o \
            $(OBJ_PATH)/kpp/OpenSMOKE_KPP_NewtonMethod_Manager.o \
            $(OBJ_PATH)/kpp/OpenSMOKE_KPP_ODE_Manager.o \
            $(OBJ_PATH)/kpp/OpenSMOKE_KPP_ReactorNetwork_Residuals.o \
            $(OBJ_PATH)/kpp/OpenSMOKE_KPP_SingleReactor_KineticsManager.o \
            $(OBJ_PATH)/kpp/OpenSMOKE_KPP_Communicator.o \

# ------------------------------------------------------------------------------------------------
# LIBRARIES
# ------------------------------------------------------------------------------------------------
ifeq ($(COMPILER_NAME),intel)
LOC_LIBS 	=  -lOpenSMOKE_Basic_INTEL           \
                   -lOpenSMOKE_Engine_INTEL          \
                   -lOpenSMOKE_AddOns_INTEL          \
		   -lOpenSMOKE_Distributions_INTEL   \
		   -lOpenSMOKE_SymbolicKinetics_INTEL      
endif	  

ifeq ($(COMPILER_NAME),gnu)
LOC_LIBS	=  -lOpenSMOKE_Basic_GNU_MKL           \
                   -lOpenSMOKE_Engine_GNU_MKL          \
                   -lOpenSMOKE_AddOns_GNU_MKL          \
		   -lOpenSMOKE_Distributions_GNU_MKL   \
                   -lOpenSMOKE_SymbolicKinetics_GNU_MKL          
     
#CXX_ADD_FLAGS = 
#CXX_ADD_FLAGS =
    
endif
  
# ------------------------------------------------------------------------------------------------
# INCLUDE FILES
# ------------------------------------------------------------------------------------------------ 
INCDIRS     = -I../../src $(BZZ_INCLUDE) $(MKL_INCLUDE) $(OPENSMOKE_INCLUDE) $(NR_INCLUDE) $(LIS_INCLUDE) $(PETSC_INCLUDE)
  
# ------------------------------------------------------------------------------------------------
# LIBRARIES
# ------------------------------------------------------------------------------------------------
LIBS_PATH 	      = -L$(BZZ_PATH_LIB) -L$(MKL_PATH_LIB) -L$(OPENSMOKE_PATH_LIB) -L$(NR_PATH_LIB) -L$(LIS_PATH_LIB) -L$(BLAS_PATH_LIB) -L$(METIS_PATH_LIB) -L$(PETSC_PATH_LIB)
LIBS 		      =   $(LOC_LIBS) $(BZZ_LIBS) $(NR_LIBS) $(LIS_LIBS) $(BLAS_LIBS) $(MKL_LIBS) $(PETSC_LIBS) -lX11 -lgfortran -lpthread  -lssl -lcrypto

# ------------------------------------------------------------------------------------------------
# COMPILER
# ------------------------------------------------------------------------------------------------
CPP       = $(CXX) 	$(CXX_FLAGS) -c $(CXX_ADD_FLAGS)	$(INCDIRS)
LINKER    = $(CXX) 	-O3	 	$(CXX_ADD_FLAGS)	$(INCDIRS) 	$(LIBS_PATH)

# ------------------------------------------------------------------------------------------------
# DETAILS
# ------------------------------------------------------------------------------------------------

OBJS     =  $(OBJ_PATH)/OpenSMOKE_KPP.o \
	    $(OBJ_PATH)/OpenSMOKE_KPP_Definitions.o \
            $(OBJ_PATH)/linear_solvers/OpenSMOKE_DirectLinearSolver_Unsymmetric.o \
            $(OBJ_PATH)/linear_solvers/OpenSMOKE_LAPACK_Dense.o \
            $(OBJ_PATH)/linear_solvers/OpenSMOKE_LIS_Unsymmetric.o \
            $(OBJ_PATH)/linear_solvers/OpenSMOKE_PARDISO_Unsymmetric.o \
            $(OBJ_PATH)/kpp/OpenSMOKE_KPP_SingleReactor.o \
            $(OBJ_PATH)/kpp/OpenSMOKE_KPP_DataManager.o \
            $(OBJ_PATH)/kpp/OpenSMOKE_KPP_Dictionary.o \
            $(OBJ_PATH)/kpp/OpenSMOKE_KPP_ReactorNetwork.o \
            $(OBJ_PATH)/kpp/OpenSMOKE_KPP_ConvectiveNetworkStatistics.o \
            $(OBJ_PATH)/kpp/OpenSMOKE_KPP_SingleReactorStatistics.o \
            $(OBJ_PATH)/kpp/OpenSMOKE_KPP_NewtonMethod_Manager.o \
            $(OBJ_PATH)/kpp/OpenSMOKE_KPP_ODE_Manager.o \
            $(OBJ_PATH)/kpp/OpenSMOKE_KPP_ReactorNetwork_Residuals.o \
            $(OBJ_PATH)/kpp/OpenSMOKE_KPP_SingleReactor_KineticsManager.o \
	    $(OBJ_PATH)/kpp/OpenSMOKE_KPP_BlockMatrixNetwork.o \
            $(OBJ_PATH)/kpp/OpenSMOKE_KPP_Communicator.o \



default: OpenSMOKE_KinPP_Splitting

$(OBJ_PATH)/OpenSMOKE_KPP.o: 	../../src/OpenSMOKE_KPP.cpp 
					$(CPP) ../../src/OpenSMOKE_KPP.cpp  -o $(OBJ_PATH)/OpenSMOKE_KPP.o

$(OBJ_PATH)/OpenSMOKE_KPP_Definitions.o: 	../../src/OpenSMOKE_KPP_Definitions.cpp 
					$(CPP) ../../src/OpenSMOKE_KPP_Definitions.cpp  -o $(OBJ_PATH)/OpenSMOKE_KPP_Definitions.o

$(OBJ_PATH)/linear_solvers/OpenSMOKE_DirectLinearSolver_Unsymmetric.o: ../../src/linear_solvers/OpenSMOKE_DirectLinearSolver_Unsymmetric.cpp ../../src/linear_solvers/OpenSMOKE_DirectLinearSolver_Unsymmetric.h
					$(CPP) ../../src/linear_solvers/OpenSMOKE_DirectLinearSolver_Unsymmetric.cpp  -o $(OBJ_PATH)/linear_solvers/OpenSMOKE_DirectLinearSolver_Unsymmetric.o

$(OBJ_PATH)/linear_solvers/OpenSMOKE_LAPACK_Dense.o: ../../src/linear_solvers/OpenSMOKE_LAPACK_Dense.cpp ../../src/linear_solvers/OpenSMOKE_LAPACK_Dense.h
					$(CPP) ../../src/linear_solvers/OpenSMOKE_LAPACK_Dense.cpp  -o $(OBJ_PATH)/linear_solvers/OpenSMOKE_LAPACK_Dense.o

$(OBJ_PATH)/linear_solvers/OpenSMOKE_PARDISO_Unsymmetric.o: ../../src/linear_solvers/OpenSMOKE_PARDISO_Unsymmetric.cpp ../../src/linear_solvers/OpenSMOKE_PARDISO_Unsymmetric.h
					$(CPP) ../../src/linear_solvers/OpenSMOKE_PARDISO_Unsymmetric.cpp  -o $(OBJ_PATH)/linear_solvers/OpenSMOKE_PARDISO_Unsymmetric.o

$(OBJ_PATH)/linear_solvers/OpenSMOKE_LIS_Unsymmetric.o: ../../src/linear_solvers/OpenSMOKE_LIS_Unsymmetric.cpp ../../src/linear_solvers/OpenSMOKE_LIS_Unsymmetric.h
					$(CPP) ../../src/linear_solvers/OpenSMOKE_LIS_Unsymmetric.cpp  -o $(OBJ_PATH)/linear_solvers/OpenSMOKE_LIS_Unsymmetric.o

$(OBJ_PATH)/kpp/OpenSMOKE_KPP_ReactorNetwork.o: 	../../src/kpp/OpenSMOKE_KPP_ReactorNetwork.cpp ../../src/kpp/OpenSMOKE_KPP_ReactorNetwork.h
					$(CPP) ../../src/kpp/OpenSMOKE_KPP_ReactorNetwork.cpp  -o $(OBJ_PATH)/kpp/OpenSMOKE_KPP_ReactorNetwork.o


$(OBJ_PATH)/kpp/OpenSMOKE_KPP_SingleReactor.o: ../../src/kpp/OpenSMOKE_KPP_SingleReactor.cpp ../../src/kpp/OpenSMOKE_KPP_SingleReactor.h
					$(CPP) ../../src/kpp/OpenSMOKE_KPP_SingleReactor.cpp  -o $(OBJ_PATH)/kpp/OpenSMOKE_KPP_SingleReactor.o

$(OBJ_PATH)/kpp/OpenSMOKE_KPP_DataManager.o: ../../src/kpp/OpenSMOKE_KPP_DataManager.cpp ../../src/kpp/OpenSMOKE_KPP_DataManager.h
					$(CPP) ../../src/kpp/OpenSMOKE_KPP_DataManager.cpp  -o $(OBJ_PATH)/kpp/OpenSMOKE_KPP_DataManager.o

$(OBJ_PATH)/kpp/OpenSMOKE_KPP_Dictionary.o: ../../src/kpp/OpenSMOKE_KPP_Dictionary.cpp ../../src/kpp/OpenSMOKE_KPP_Dictionary.h
					$(CPP) ../../src/kpp/OpenSMOKE_KPP_Dictionary.cpp  -o $(OBJ_PATH)/kpp/OpenSMOKE_KPP_Dictionary.o

$(OBJ_PATH)/kpp/OpenSMOKE_KPP_ConvectiveNetworkStatistics.o: ../../src/kpp/OpenSMOKE_KPP_ConvectiveNetworkStatistics.cpp ../../src/kpp/OpenSMOKE_KPP_ConvectiveNetworkStatistics.h
					$(CPP) ../../src/kpp/OpenSMOKE_KPP_ConvectiveNetworkStatistics.cpp  -o $(OBJ_PATH)/kpp/OpenSMOKE_KPP_ConvectiveNetworkStatistics.o

$(OBJ_PATH)/kpp/OpenSMOKE_KPP_SingleReactorStatistics.o: ../../src/kpp/OpenSMOKE_KPP_SingleReactorStatistics.cpp ../../src/kpp/OpenSMOKE_KPP_SingleReactorStatistics.h
					$(CPP) ../../src/kpp/OpenSMOKE_KPP_SingleReactorStatistics.cpp  -o $(OBJ_PATH)/kpp/OpenSMOKE_KPP_SingleReactorStatistics.o

$(OBJ_PATH)/kpp/OpenSMOKE_KPP_NewtonMethod_Manager.o: ../../src/kpp/OpenSMOKE_KPP_NewtonMethod_Manager.cpp ../../src/kpp/OpenSMOKE_KPP_NewtonMethod_Manager.h
					$(CPP) ../../src/kpp/OpenSMOKE_KPP_NewtonMethod_Manager.cpp  -o $(OBJ_PATH)/kpp/OpenSMOKE_KPP_NewtonMethod_Manager.o

$(OBJ_PATH)/kpp/OpenSMOKE_KPP_ODE_Manager.o: ../../src/kpp/OpenSMOKE_KPP_ODE_Manager.cpp ../../src/kpp/OpenSMOKE_KPP_ODE_Manager.h
					$(CPP) ../../src/kpp/OpenSMOKE_KPP_ODE_Manager.cpp  -o $(OBJ_PATH)/kpp/OpenSMOKE_KPP_ODE_Manager.o

$(OBJ_PATH)/kpp/OpenSMOKE_KPP_ReactorNetwork_Residuals.o: ../../src/kpp/OpenSMOKE_KPP_ReactorNetwork_Residuals.cpp ../../src/kpp/OpenSMOKE_KPP_ReactorNetwork_Residuals.h
					$(CPP) ../../src/kpp/OpenSMOKE_KPP_ReactorNetwork_Residuals.cpp  -o $(OBJ_PATH)/kpp/OpenSMOKE_KPP_ReactorNetwork_Residuals.o

$(OBJ_PATH)/kpp/OpenSMOKE_KPP_SingleReactor_KineticsManager.o: ../../src/kpp/OpenSMOKE_KPP_SingleReactor_KineticsManager.cpp ../../src/kpp/OpenSMOKE_KPP_SingleReactor_KineticsManager.h
					$(CPP) ../../src/kpp/OpenSMOKE_KPP_SingleReactor_KineticsManager.cpp  -o $(OBJ_PATH)/kpp/OpenSMOKE_KPP_SingleReactor_KineticsManager.o

$(OBJ_PATH)/kpp/OpenSMOKE_KPP_BlockMatrixNetwork.o: ../../src/kpp/OpenSMOKE_KPP_BlockMatrixNetwork.cpp ../../src/kpp/OpenSMOKE_KPP_BlockMatrixNetwork.h
					$(CPP) ../../src/kpp/OpenSMOKE_KPP_BlockMatrixNetwork.cpp  -o $(OBJ_PATH)/kpp/OpenSMOKE_KPP_BlockMatrixNetwork.o

$(OBJ_PATH)/kpp/OpenSMOKE_KPP_Communicator.o: ../../src/kpp/OpenSMOKE_KPP_Communicator.cpp ../../src/kpp/OpenSMOKE_KPP_Communicator.h
					$(CPP) ../../src/kpp/OpenSMOKE_KPP_Communicator.cpp  -o $(OBJ_PATH)/kpp/OpenSMOKE_KPP_Communicator.o


OpenSMOKE_KinPP_Splitting:  	$(OBJS)  
			$(LINKER) $(OBJS) -o ../../exe/OpenSMOKE_KinPP_Splitting_$(STRING_NAME).sh $(LIBS)

clean: ; @rm $(OBJ_PATH)/*.o $(OBJ_PATH)/kpp/*.o $(OBJ_PATH)/linear_solvers/*.o ../../exe/OpenSMOKE_KinPP_Splitting_$(STRING_NAME).sh
