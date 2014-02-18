#ifndef BZZCSTRNETWORKHPP
#define BZZCSTRNETWORKHPP

#include "engine/OpenSMOKE_ReactingGas.h"

class myBzzCSTRNetwork : public BzzBaseClass
{
	// Declaration of friend classes
	// ----------------------------------------------------------------------------
	friend class BzzNonLinearSystemDoubleSparse;
	friend class MyOdeSystemObjectOneCSTR;
	friend class MyOdeSystemObjectAllCSTR;
	
private:

	// CONSTANTS
	// ----------------------------------------------------------------------------
	static const double			R_CTOT;
	static const double			UR_CTOT;
	static const char *const	BZZ_ERROR;
	
	// Static Variables
	// ----------------------------------------------------------------------------
	static int count;
	static int countInScope;	
	int	whoAmI;
	int printTasks;
	int printSubTasks;


	// Memory Management options
	// ----------------------------------------------------------------------------
	int memoTemperature;	// 1 salva su file 0 memorizza

	// Tolerances
	// ----------------------------------------------------------------------------
	double tolRel;
	double tolAbs;


	// REACTOR NETWORK DIMENSION
	// ----------------------------------------------------------------------------
	int	numComponents;
	int numReactions;
	int numCSTRReactors;
	int numComponentsReduced;
	int originalNumCSTRReactors;


	// REACTOR NETWORK VARIABLES
	// ----------------------------------------------------------------------------
	BzzVector	temperature;		// Reactor temperature [K]
	BzzVector	pressure;			// Reactor Pressure [atm]
	BzzVector	volume;				// Reactor Volume [m3]
	BzzVector	qanm;				// Reactor Squared Normalized temperature variance [-]
	BzzVector	massRate;			// Reactor convective mass flow rate [kg/s]	
	BzzVector	massInput;			// Reactor External mass flow rate Input [kg/s]
	BzzVector	massOutput;			// Reactor External mass flow rate Output [kg/s]
	BzzVector	cTotm;				// Reactor Total concentration [kmol/m3]

	BzzVector	logTm;				// Reactor Auxiliary variable	log(T)
	BzzVector	loguRTm;			// Reactor Auxiliary variable	log(1/RT)
	BzzVector	feedInkRreactor;

	// Matrix Version
	BzzMatrix	massFractionsInReactorsSolution,		// Mass fractions
					massFractionsInReactors,				// Mass Fractions
					volumeMolecularWeightT,					// V * PM
					feed;									// Input mass flow rate
					
	// Vector Version
	BzzVector	massFractionsInReactorsSolution_Vector,	// Mass fractions
					massFractionsInReactors_Vector,			// Mass fractions
					volumeMolecularWeightT_Vector,			// V * PM
					feed_Vector;							// Input mass flow rate


	// COEFFICIENTS FOR NON LINEAR SYSTEM OF EQUATIONS
	// ----------------------------------------------------------------------------
	BzzMatrixSparse	Mg;
	BzzMatrixSparse	Ms;
	BzzVector			sM;
	BzzMatrixSparseLockedByRows	Ls;
	BzzMatrixSparseLockedByRows	Ld;
	BzzFactorizedSparseLockedByRowsGauss Fg;


	// CONVERGENCE VARIABLES
	// ----------------------------------------------------------------------------
	double	F1Stop;							// Stop criterium on the residual						
	int		maxNumberOfIterations_Global;	// Maximum number of iterations in Get First
	

	// DUMMY VARIABLES
	// ----------------------------------------------------------------------------
	BzzMatrix dummyMatrix;			// Auxiliar Matrix
	BzzVector dummyVector;			// Auxiliar Vector


	// SINGLE REACTOR VARIABLES
	// ----------------------------------------------------------------------------
	int		kReactor;	// Reactor Index
	double	T,			// Reactor temperature			(Fixed)		
			P,			// Reactor pressure				(Fixed)
			logT,		// Reactor log(T)				(Fixed)
			loguRT,		// Reactor log(1/RT)			(Fixed)
			cTot,		// Reactor Total Concentration	(Fixed)
			wM;			// Reactor Molecular weight		(Variable)
	

	// SINGLE REACTOR VARIABLES
	// ----------------------------------------------------------------------------
	BzzVector omegaReactor,		// mass fractions
					cReactor,			// concentrations
					xReactor,			// mole fractions
					RReactor;			// reaction rates


	// CLUSTERING ALGORITHM
	// ----------------------------------------------------------------------------
	BzzVectorInt		cluster;
	BzzVectorInt		giveClusterIndexFromCellIndex;
	BzzVectorIntArray	cstrOut;
	BzzVectorIntArray	cstrConnect;
	BzzVectorInt		iSpec;			// Reduced species indices in the complete scheme
	BzzVectorInt		cstrSequence;	// Reactor Sequence


	// GAS PHASE MIXTURE
	// ----------------------------------------------------------------------------
	OpenSMOKE_ReactingGas Reactions;


	// INITIALIZATION ALGORITHM
	// ----------------------------------------------------------------------------
	void read_first_guess(char *first);
	void readKineticSheme();
	void reading_the_network_file(char *fileNetwork, int &countInput);
	void build_cluster(int cicloCluster, int numCSTRReactors, int countInput, int &numCluster, BzzVectorInt &cstrClusterSize);
	void calculate_initial_massfractions_in_clusters(int numCluster, BzzMatrix &omegaReduced, BzzMatrix &aux);
	void calculate_clustering(char *fileNetwork, int cicloDiffusion, int numCluster, BzzVectorInt &cstrClusterSize, BzzMatrix &auxMass, BzzMatrix &aux);
	void read_tolerances(int numCSTRReactors, int cicloCluster);
	void complete_clustering(int countInput, int cicloDiffusion, int relaxation, int iaia, BzzMatrix &auxMass);
	void MemoTemperatureFunctions(void);


	// TOLERANCES FOR CLUSTERING ALGORITHM
	// ----------------------------------------------------------------------------	
	int		numTolT;
	BzzVector tInf; 
	BzzVector tSup;
	BzzVector dt;
	BzzVector dtBase;
	double	tolRelC;
	double	tolAbsC;
	double	tolRelCBase;
	double	tolAbsCBase;


	// SOLVE NETWORK
	// ----------------------------------------------------------------------------	
	void	start_solving_network();
	void	globalAlgorithm(void);
	int		GetFirst(void);
	void	GetSecond(void);
	int		GetThird(void);


	// GET RESIDUALS
	// ----------------------------------------------------------------------------
	void GetResiduals(BzzMatrix &massFractionsInReactors,BzzMatrix &residuals);
	void GetAllResiduals(BzzVector &m,BzzVector &residuals);
	void GetResiduals(BzzVector &omega,BzzVector &residuals);
	void GetJacobian(BzzVector &x,BzzMatrix &JJ);
	void GetJacobianForSingleReactor(int iReactor, BzzVector &cRes, BzzVector &R, BzzMatrix &dRC);


	// GET REACTION RATES
	// ----------------------------------------------------------------------------	
	void	GetReactionsRateInAllReactorsFromMassFractions
				(BzzMatrix &massFractionsInReactors, 
				BzzMatrix &reactionsRateInReactors);
	void	GetReactionsRateInAllReactorsFromMassFractions
				(BzzVector &massFractionsInReactorsV, 
				BzzVector &reactionsRateInReactors);
	void	GetReactionRatesInSingleReactor(int iReactor,
				BzzVector &omega, BzzVector &molefractions,
				BzzVector &concentrations, BzzVector &ReactionRates);


	// OTHER FUNCTIONS
	// ----------------------------------------------------------------------------
	void GetDiagonalFactoredMatricesForLinearizedSistem(BzzMatrix &massFractionsInReactors);
	void GetDiagonalMatricesForLinearizedSistem(char *file,BzzVector &mfV);
	void SwapMassFractionsInReactors(BzzMatrix *massFractionsInReactors);


	void myBzzCSTRNetwork::calculate_massflowrate_in_reactors(BzzVector &inputMassFlowRate,
														  BzzVector &outputMassFlowRate);


public:

	// CONSTRUCTORS
	// ----------------------------------------------------------------------------
	myBzzCSTRNetwork(void); 

	// OPERATORS
	// ----------------------------------------------------------------------------
	void operator()(char *fileNetwork,char *first, int cicloCluster,int cicloDiffusion,int iaia,int relaxation);
	
	// SAVE NETWORK
	// ----------------------------------------------------------------------------
	void Save(char *file);		// formatted
	void Save(char,char *file);	// binary
	
	// UTILITIES
	// ----------------------------------------------------------------------------
	int		WhoAmI(void) const;
	void	SetTolRel(double tolr);
	void	SetTolAbs(double tola);
	void	SetMemoTemperature(void);
	int		GetNumComponents(void);
	int		GetNumReactions(void);
	int		GetNumCSTRReactors(void);
	void	SetTasksPrint(void);
	void	SetSubTasksPrint(void);
	void	SetSolution(BzzMatrix &omega);

	// STATIC FUNCTIONS
	// ----------------------------------------------------------------------------
	static int ObjectCount(void)		{	return count;	}
	static int ObjectCountInScope(void)	{	return countInScope;	}
	
	// PRINT FUNCTIONS
	// ----------------------------------------------------------------------------
	virtual void ObjectBzzPrint(void);
	void OutputPrint(int ciclo);

	// ASSIGNEMENT OPERATORS
	// ----------------------------------------------------------------------------
	myBzzCSTRNetwork &operator = (const myBzzCSTRNetwork &rval);

	// MODIFYING FUNCTIONS
	// ----------------------------------------------------------------------------
	friend void Delete(myBzzCSTRNetwork *result);

	// DESTRUCTORS
	// ----------------------------------------------------------------------------
	~myBzzCSTRNetwork(void);
};


#endif // BZZ_CSTR_NETWORK_HPP

/*
****************************************************
                    Final Cicle
                CSTR Number  13

****************************************************

====================================================
              Exit for all elements
====================================================
                    RATE [kg/s]      MASS FRACT.
    1   H          1.285326e-007    4.423336e-006
    2   O2         5.853392e-003    2.014393e-001
    3   OH         1.622825e-005    5.584807e-004
    4   O          1.352896e-005    4.655871e-004
    5   H2         4.292254e-007    1.477141e-005
    6   H2O        4.086805e-004    1.406438e-002
    7   CO         8.866866e-005    3.051453e-003
    8   CO2        3.692253e-004    1.270656e-002
    9   HO2        5.741170e-008    1.975773e-006
   10   AR         0.000000e+000    0.000000e+000
   11   C2H6       4.324261e-007    1.488156e-005
   12   H2O2       8.786216e-009    3.023698e-007
   13   CH4        2.440110e-005    8.397423e-004
   14   CHO        1.825016e-010    6.280632e-009
   15   CH2O       2.664517e-006    9.169700e-005
   16   CH3        1.421532e-008    4.892078e-007
   17   T-CH2      7.753194e-012    2.668193e-010
   18   S-CH2      1.162198e-012    3.999602e-011
   19   C2H4       1.049443e-006    3.611565e-005
   20   CH3O       2.924155e-011    1.006322e-009
   21   C2H5       2.644289e-012    9.100086e-011
   22   CH3OH      8.142478e-009    2.802162e-007
   23   CH2OH      4.809354e-011    1.655096e-009
   24   CH         1.165233e-013    4.010047e-012
   25   C2H2       7.182578e-008    2.471821e-006
   26   C2H3       1.624088e-011    5.589155e-010
   27   CH2CHO     5.753496e-011    1.980015e-009
   28   C2H4O      8.415729e-009    2.896199e-007
   29   CH2CO      1.098995e-007    3.782092e-006
   30   HCCO       4.104197e-012    1.412423e-010
   31   C2H        2.560755e-013    8.812610e-012
   32   C3H3       9.342417e-011    3.215110e-009
   33   C3H4       5.674296e-007    1.952759e-005
   34   C3H5       1.063565e-008    3.660164e-007
   35   C3H6       3.725083e-007    1.281954e-005
   36   NC3H7      9.533942e-012    3.281022e-010
   37   C3H8       1.112382e-006    3.828164e-005
   38   IC3H7      5.581840e-011    1.920941e-009
   39   N2         2.216656e-002    7.628426e-001
   40   NO         5.272386e-005    1.814445e-003
   41   N          5.172198e-010    1.779966e-008
   42   NH         3.875875e-011    1.333848e-009
   43   N2O        5.910369e-008    2.034001e-006
   44   N2H        9.727103e-012    3.347496e-010
   45   HCN        1.259573e-008    4.334709e-007
   46   NCO        1.267083e-011    4.360553e-010
   47   HNCO       2.799183e-009    9.633140e-008
   48   NH2        3.546712e-012    1.220570e-010
   49   CN         4.160022e-014    1.431635e-012
   50   HNO        3.588918e-008    1.235094e-006
   51   NH3        1.270180e-011    4.371213e-010
   52   NO2        5.727672e-005    1.971128e-003

   Total           2.905785e-002    1.000000e+000
=====================================================
Cpu Seconds for complete solution: 1.965827e+001
User Seconds for complete solution: 1.845654e+001
Kernel Seconds for complete solution: 1.201728e+000
>>>>>>Press any key to continue
Press any key to continue
*/
