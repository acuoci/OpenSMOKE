# ------------------------------------------------------------------------------------------------
# OBJECTS
# ------------------------------------------------------------------------------------------------
OBJS     =  $(OBJ_PATH)/OpenSMOKE_KPP.o \
	    $(OBJ_PATH)/OpenSMOKE_KPP_Definitions.o \
            $(OBJ_PATH)/linear_solvers/OpenSMOKE_DirectLinearSolver_Unsymmetric.o \
            $(OBJ_PATH)/linear_solvers/OpenSMOKE_LAPACK_Dense.o \
            $(OBJ_PATH)/linear_solvers/OpenSMOKE_LIS_Unsymmetric.o \
            $(OBJ_PATH)/linear_solvers/OpenSMOKE_MUMPS_Unsymmetric.o \
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
            $(OBJ_PATH)/kpp/OpenSMOKE_KPP_SingleReactor_KineticsManager.o

# ------------------------------------------------------------------------------------------------
# LIBRARIES
# ------------------------------------------------------------------------------------------------

LOC_LIBS	=  -lOpenSMOKE_Basic_$(STRING_NAME)           \
                   -lOpenSMOKE_Engine_$(STRING_NAME)          \
                   -lOpenSMOKE_AddOns_$(STRING_NAME)          \
		   -lOpenSMOKE_Distributions_$(STRING_NAME)         \
                   -lOpenSMOKE_SymbolicKinetics_$(STRING_NAME)          
    
ifeq ($(COMPILER_NAME),gnu)
CXX_ADD_FLAGS = -fopenmp    
endif

ifeq ($(COMPILER_NAME),intel)
CXX_ADD_FLAGS = -openmp    
endif
  
# ------------------------------------------------------------------------------------------------
# INCLUDE FILES
# ------------------------------------------------------------------------------------------------ 
INCDIRS     = -I../../src $(BZZ_INCLUDE) $(MKL_INCLUDE) $(OPENSMOKE_INCLUDE) $(NR_INCLUDE) $(LIS_INCLUDE) $(MUMPS_INCLUDE)
  
# ------------------------------------------------------------------------------------------------
# LIBRARIES
# ------------------------------------------------------------------------------------------------
LIBS_PATH 	      = -L$(BZZ_PATH_LIB) -L$(MKL_PATH_LIB) -L$(OPENSMOKE_PATH_LIB) -L$(NR_PATH_LIB) -L$(LIS_PATH_LIB) -L$(MUMPS_PATH_LIB) -L$(BLAS_PATH_LIB) -L$(METIS_PATH_LIB) 
LIBS 		      =   $(LOC_LIBS) $(BZZ_LIBS) $(NR_LIBS) $(LIS_LIBS) $(MUMPS_LIBS) $(BLAS_LIBS) $(MKL_LIBS) -lgfortran -lpthread

# ------------------------------------------------------------------------------------------------
# COMPILER
# ------------------------------------------------------------------------------------------------
CPP       = $(CXX) 	$(CXX_FLAGS) -c $(CXX_ADD_FLAGS)	$(INCDIRS)
LINKER    = $(CXX) 	-O3 -static 	$(CXX_ADD_FLAGS)	$(INCDIRS) 	$(LIBS_PATH)

# ------------------------------------------------------------------------------------------------
# DETAILS
# ------------------------------------------------------------------------------------------------

OBJS     =  $(OBJ_PATH)/OpenSMOKE_KPP.o \
	    $(OBJ_PATH)/OpenSMOKE_KPP_Definitions.o \
            $(OBJ_PATH)/linear_solvers/OpenSMOKE_DirectLinearSolver_Unsymmetric.o \
            $(OBJ_PATH)/linear_solvers/OpenSMOKE_LAPACK_Dense.o \
            $(OBJ_PATH)/linear_solvers/OpenSMOKE_LIS_Unsymmetric.o \
            $(OBJ_PATH)/linear_solvers/OpenSMOKE_MUMPS_Unsymmetric.o \
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
	    $(OBJ_PATH)/kpp/OpenSMOKE_KPP_BlockMatrixNetwork.o



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

$(OBJ_PATH)/linear_solvers/OpenSMOKE_MUMPS_Unsymmetric.o: ../../src/linear_solvers/OpenSMOKE_MUMPS_Unsymmetric.cpp ../../src/linear_solvers/OpenSMOKE_MUMPS_Unsymmetric.h
					$(CPP) ../../src/linear_solvers/OpenSMOKE_MUMPS_Unsymmetric.cpp  -o $(OBJ_PATH)/linear_solvers/OpenSMOKE_MUMPS_Unsymmetric.o

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


OpenSMOKE_KinPP_Splitting:  	$(OBJS)  
			$(LINKER) $(OBJS) -o ../../exe/OpenSMOKE_KinPP_Splitting_$(STRING_NAME).sh $(LIBS)

clean: ; @rm $(OBJ_PATH)/*.o $(OBJ_PATH)/kpp/*.o $(OBJ_PATH)/linear_solvers/*.o ../../exe/OpenSMOKE_KinPP_Splitting_$(STRING_NAME).sh
