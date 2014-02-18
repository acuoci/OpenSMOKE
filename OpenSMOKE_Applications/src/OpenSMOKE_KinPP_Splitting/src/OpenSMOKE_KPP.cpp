/***************************************************************************
 *   Copyright (C) 2011 by Alberto Cuoci								   *
 *   alberto.cuoci@polimi.it                                               *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 *   This program is distributed in the hope that it will be useful,       *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 *   GNU General Public License for more details.                          *
 *                                                                         *
 *   You should have received a copy of the GNU General Public License     *
 *   along with this program; if not, write to the                         *
 *   Free Software Foundation, Inc.,                                       *
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             *
 ***************************************************************************/

#include "lis.h"
#include <iomanip>
#include <omp.h>
#include "kpp/OpenSMOKE_KPP_Dictionary.h"
#include "kpp/OpenSMOKE_KPP_DataManager.h"
#include "kpp/OpenSMOKE_KPP_ReactorNetwork.h"
#include "linear_solvers/OpenSMOKE_LAPACK_Dense.h"
#include "kpp/OpenSMOKE_KPP_BlockMatrixNetwork.h"
#include "kpp/OpenSMOKE_KPP_ReactorNetwork_Residuals.h"



void InputOutputAnalysis(OpenSMOKE_ReactingGas& mix, OpenSMOKE_KPP_ReactorNetwork& network);
void MassUmbalances(OpenSMOKE_KPP_ReactorNetwork& network);
int  SequenceCSTR(OpenSMOKE_KPP_ReactorNetwork& network, OpenSMOKE_KPP_DataManager& data);
void PredictorCorrector(OpenSMOKE_KPP_ReactorNetwork& network, OpenSMOKE_KPP_DataManager& data);
int GlobalODE(OpenSMOKE_KPP_ReactorNetwork& network, OpenSMOKE_KPP_DataManager& data);
int GlobalNewtonMethod(OpenSMOKE_KPP_ReactorNetwork& network, OpenSMOKE_KPP_DataManager& data);
void prova();

ofstream fCPUTime;
ofstream fLog;
ofstream fWarning;
int globalIteration;
double 	globalNLSLinearSystemAbsoluteTolerance;
double 	globalNLSLinearSystemRelativeTolerance;
double 	globalODELinearSystemAbsoluteTolerance;
double 	globalODELinearSystemRelativeTolerance;

// Main
int main(int argc, char* argv[])
{
	// Disabling Bzz OpenMP
	//bzzOpenMP = 0;

	// Dictionary setup
	OpenSMOKE_KPP_Dictionary dictionary;
	dictionary.ParseFile("Input.inp");

	// Read input file
	OpenSMOKE_KPP_DataManager data(dictionary);

	// OpenMP
    //omp_set_dynamic(false);
	//omp_set_num_threads(12);


    cout << "Number of user defined threads: " << data.nThreads()           << endl;
	cout << "Number of current threads:      " << omp_get_num_threads() 	<< endl;
	cout << "Master thread id:               " << omp_get_thread_num() 	<< endl;
	cout << "Number of processes:            " << omp_get_num_procs() 	<< endl;
	cout << "Max number of threads:          " << omp_get_max_threads() 	<< endl;



	// Initialize LIS
	lis_initialize(&argc, &argv);

	// Gas mixture setup
    	cout << "Start" << endl; //getchar();
	OpenSMOKE_ReactingGas* mix;
	if (data.iSaveKineticConstants() == false)
	{
		mix = new OpenSMOKE_ReactingGas[data.nThreads()];
		for (int k=0;k<data.nThreads();k++)
			mix[k].SetupBinary(dictionary.kinetics());
	}
	else
	{
		mix = new OpenSMOKE_ReactingGas[1];
		mix[0].SetupBinary(dictionary.kinetics());
	}	
    	cout << "Done: Mixture" << endl; //getchar();
	
	// Reactor network setup
	OpenSMOKE_KPP_ReactorNetwork network(*mix, data, fLog, fWarning);
        cout << "Done: Network" << endl; //getchar();
	network.ReadFirstGuess();
        cout << "Done: FirstGuess" << endl; //getchar();
	network.ReadTopology();
        cout << "Done: Topology" << endl; //getchar();
	network.BuildNetwork();

	// Initial Analysis
	InputOutputAnalysis(*mix, network);
	MassUmbalances(network);
	
	string log_string = data.nameFolderOutput() + "/Log.out";
	fLog.open( log_string.c_str(), ios::out);
	fLog.setf(ios::scientific);

	string warning_string = data.nameFolderOutput() + "/Warning.out";
	fWarning.open( warning_string.c_str(), ios::out);
	fWarning.setf(ios::scientific);



	globalIteration = 1;

	// Without reactions
	if (data.iReactions() == false)
	{
		network.SolveWithoutReactions();
		network.WriteMassFractionMap();
	}
	
	// With reactions
	if (data.iReactions() == true)
	{
		double timeStartTotal = BzzGetCpuTime();
	
	//	SequenceCSTR(network, data);
	//	PredictorCorrector(network, data);
	//	GlobalODE(network, data);
	//	GlobalNewtonMethod(network, data);

		globalODELinearSystemAbsoluteTolerance = 1.e-11;
	 	globalODELinearSystemRelativeTolerance = 1.e-9;
		globalNLSLinearSystemAbsoluteTolerance = 1.e-15;
	 	globalNLSLinearSystemRelativeTolerance = 1.e-12;

		int iFlagSequence;
		int iFlagODE;
		int iFlagNewtonMethod;

		for(int jj=1;jj<=3;jj++)
		{
			iFlagSequence	  = SequenceCSTR(network, data);
			iFlagODE          = GlobalODE(network, data);

			globalODELinearSystemAbsoluteTolerance /= 10.;
			globalODELinearSystemRelativeTolerance /= 10.;
		}

		iFlagNewtonMethod	= GlobalNewtonMethod(network, data);

		for(int jj=1;jj<=2;jj++)
		{			
			globalODELinearSystemAbsoluteTolerance /= 10.;
			globalODELinearSystemRelativeTolerance /= 10.;
		
			iFlagSequence	  	= SequenceCSTR(network, data);
			iFlagODE          	= GlobalODE(network, data);	
			iFlagNewtonMethod	= GlobalNewtonMethod(network, data);
		}

		double timeEndTotal = BzzGetCpuTime();

		// Final Analysis
		InputOutputAnalysis(*mix, network);
		MassUmbalances(network);

		cout << "Total CPU Time: " << timeEndTotal - timeStartTotal << endl;
	}

	
}

void MassUmbalances(OpenSMOKE_KPP_ReactorNetwork& network)
{
	BzzVector umbalances(network.NumberOfReactors());
	for(int k=1;k<=network.NumberOfReactors();k++)
		umbalances[k] = fabs(network.reactors(k).MassUmbalance());

	BzzVector massFlowIn(network.NumberOfReactors());
	for(int k=1;k<=network.NumberOfReactors();k++)
		massFlowIn[k] = fabs(network.reactors(k).MassFlowIn());

	cout << endl;
	cout << "// *********************************************************** // " << endl;
	cout << "//                     UMBALANCE ANALYSIS                      // " << endl;
	cout << "// *********************************************************** // " << endl;
	cout << setw(20) << left << "Max umbalance:";
	cout << setw(16) << left << umbalances.Max();
	cout << endl;
	cout << setw(20) << left << "Mean umbalance:";
	cout << setw(16) << left << Mean(umbalances);
	cout << endl;
	cout << setw(20) << left << "Max in flow:";
	cout << setw(16) << left << massFlowIn.Max();
	cout << endl;
	cout << setw(20) << left << "Min in flow:";
	cout << setw(16) << left << massFlowIn.Min();
	cout << endl;
	cout << setw(20) << left << "Mean in flow:";
	cout << setw(16) << left << Mean(massFlowIn);
	cout << endl;

}

void InputOutputAnalysis(OpenSMOKE_ReactingGas& mix, OpenSMOKE_KPP_ReactorNetwork& network)
{
	double massFlowIn, massFlowOut;
	
	BzzVector omegaFeeds(mix.NumberOfSpecies());
        BzzVector omegaMin(mix.NumberOfSpecies());
        BzzVector omegaMax(mix.NumberOfSpecies());
	BzzVector omegaOutput(mix.NumberOfSpecies());
	BzzVector omegaElementalFeeds(mix.NumberOfElements());
	BzzVector omegaElementalOutput(mix.NumberOfElements());

	network.ExternalFeeds(massFlowIn, omegaFeeds, omegaElementalFeeds);
	network.ExternalOutput(massFlowOut, omegaOutput, omegaElementalOutput);
        network.MinMaxMassFractions(omegaMin, omegaMax);

	cout << endl;
	cout << "// *********************************************************** // " << endl;
	cout << "//              IN/OUT ANALYSIS: MASS FLOW RATES               // " << endl;
	cout << "// *********************************************************** // " << endl;
	cout << setw(16) << left << "Input:";
	cout << setw(16) << left << massFlowIn;
	cout << endl;
	cout << setw(16) << left << "Output:";
	cout << setw(16) << left << massFlowOut;
	cout << endl;
	cout << setw(16) << left << "Difference(%):";
	cout << setw(16) << left << (massFlowOut-massFlowIn)/massFlowIn*100.;
	cout << endl;

	
	cout << endl;
	cout << "// *********************************************************** // " << endl;
	cout << "//           IN/OUT ANALYSIS: SPECIES MASS FRACTIONS           // " << endl;
	cout << "// *********************************************************** // " << endl;
	for (int j=1;j<=mix.NumberOfSpecies();j++)
		if (omegaFeeds[j]>1.e-12 || omegaOutput[j]>1.e-12)		
		{
			cout << setw(16) << left << mix.names[j];
			cout << setw(16) << left << omegaFeeds[j];
			cout << setw(16) << left << omegaOutput[j];
			cout << endl;
		}

	cout << endl;
	cout << "// *********************************************************** // " << endl;
	cout << "//           IN/OUT ANALYSIS: ELEMENT MASS FRACTIONS           // " << endl;
	cout << "// *********************************************************** // " << endl;
	for (int j=1;j<=mix.NumberOfElements();j++)
		if (omegaElementalFeeds[j]>1.e-12 || omegaElementalOutput[j]>1.e-12)		
		{
			cout << setw(16) << left << mix.list_of_elements[j-1];
			cout << setw(16) << left << omegaElementalFeeds[j];
			cout << setw(16) << left << omegaElementalOutput[j];
			cout << endl;
		}
        
        string fout_name = network.data().nameFolderOutput() + "/Solution/Summary.out";
        ofstream fout;
        fout.open(fout_name.c_str(), ios::out);
        fout.setf(ios::scientific);
        
	fout << endl;
	fout << "// *********************************************************** // " << endl;
	fout << "//           IN/OUT ANALYSIS: SPECIES MASS FRACTIONS           // " << endl;
	fout << "// *********************************************************** // " << endl;
	for (int j=1;j<=mix.NumberOfSpecies();j++)
        {
                fout << setw(16) << left << mix.names[j];
                fout << setw(16) << left << omegaFeeds[j];
                fout << setw(16) << left << omegaOutput[j];
                fout << endl;
        }

	fout << endl;
	fout << "// *********************************************************** // " << endl;
	fout << "//           IN/OUT ANALYSIS: ELEMENT MASS FRACTIONS           // " << endl;
	fout << "// *********************************************************** // " << endl;
	for (int j=1;j<=mix.NumberOfElements();j++)	
        {
                fout << setw(16) << left << mix.list_of_elements[j-1];
                fout << setw(16) << left << omegaElementalFeeds[j];
                fout << setw(16) << left << omegaElementalOutput[j];
                fout << endl;
        }
        
	fout << endl;
	fout << "// *********************************************************** // " << endl;
	fout << "//          MIN/MAX ANALYSIS: SPECIES MASS FRACTIONS           // " << endl;
	fout << "// *********************************************************** // " << endl;
	for (int j=1;j<=mix.NumberOfSpecies();j++)
        {
                fout << setw(16) << left << mix.names[j];
                fout << setw(16) << left << omegaMin[j];
                fout << setw(16) << left << omegaMax[j];
                fout << endl;
        }        
        
        fout.close();
}

int SequenceCSTR(OpenSMOKE_KPP_ReactorNetwork& network, OpenSMOKE_KPP_DataManager& data)
{
	fLog << "****************************************************************************** " << endl; 
	fLog << "*                              SEQUENCE                                      * " << endl;
	fLog << "****************************************************************************** " << endl; 

	int info = network.SequenceCSTR();

	fLog << endl << endl;

	// Write Mass Fraction Maps
	if (info>=0)
	{
		network.WriteBackupFile(data.nameOutputBackupFile());
		network.WriteMassFractionMap();
	}

	return info;
}

void PredictorCorrector(OpenSMOKE_KPP_ReactorNetwork& network, OpenSMOKE_KPP_DataManager& data)
{
	data.SetNetworkStatus(KPP_NETWORK_STATUS_SPLITTING);

	double deltat = data.PredictorCorrector_InitialTimeStep();

	// Convection-Diffusion Matrix
	network.SetTimeStep(deltat);

	// Loop
	for (int i=1;i<=data.PredictorCorrector_MaxIterations();i++)
	{
		double timeStartGlobal = BzzGetCpuTime();
	
		network.PredictorCorrector(deltat);
				
		double timeEndGlobal = BzzGetCpuTime();

		network.ResidualsAnalysis();

		network.TimeStepPolicy(deltat);
	}

	// Write Mass Fraction Maps
	//if (info>=0)
	{
		network.WriteBackupFile(data.nameOutputBackupFile());
		network.WriteMassFractionMap();
	}
}

int GlobalODE(	OpenSMOKE_KPP_ReactorNetwork& network, OpenSMOKE_KPP_DataManager& data)
{
	fLog << "****************************************************************************** " << endl; 
	fLog << "*                            GLOBAL ODE                                      * " << endl;
	fLog << "****************************************************************************** " << endl; 

	// Initialize ODE System
	data.SetNetworkStatus(KPP_NETWORK_STATUS_GLOBALODE);
	network.InitializeGlobal(data.GlobalODE_SparseLinearSolver(), globalODELinearSystemAbsoluteTolerance, globalODELinearSystemRelativeTolerance);

	// Solve ODE System
	int info = network.GlobalODE();

	fLog << endl << endl;

	// Write Mass Fraction Maps
//	if (info>=0)
	{
		network.WriteBackupFile(data.nameOutputBackupFile());
		network.WriteMassFractionMap();
	}

	return info;
}

int GlobalNewtonMethod(OpenSMOKE_KPP_ReactorNetwork& network, OpenSMOKE_KPP_DataManager& data)
{
	fLog << "****************************************************************************** " << endl; 
	fLog << "*                            GLOBAL NLS                                      * " << endl;
	fLog << "****************************************************************************** " << endl; 
	
	// Initialize Non Linear System
	data.SetNetworkStatus(KPP_NETWORK_STATUS_GLOBALNLS);
	network.InitializeGlobal(data.GlobalNLS_SparseLinearSolver(), globalNLSLinearSystemAbsoluteTolerance, globalNLSLinearSystemRelativeTolerance);

	// Solve Non Linear System
	int info = network.GlobalNLS();

	fLog << endl << endl;

	// Write Mass Fraction Maps
	{
		network.WriteBackupFile(data.nameOutputBackupFile());
		network.WriteMassFractionMap();
	}

	return info;
}
/*
//#include "mpi.h"
#include "dmumps_c.h"
#define JOB_INIT -1
#define JOB_END -2
#define USE_COMM_WORLD -987654

int main(int argc, char ** argv)
{
  DMUMPS_STRUC_C id;
  int n = 2;
  int nz = 2;
  int irn[] = {1,2};
  int jcn[] = {1,2};
  double a[2];
  double rhs[2];

  int myid, ierr;

 // ierr = MPI_Init(&argc, &argv);
//  ierr = MPI_Comm_rank(MPI_COMM_WORLD, &myid);
  // Define A and rhs 
  rhs[0]=1.0;rhs[1]=4.0;
  a[0]=1.0;a[1]=2.0;

  // Initialize a MUMPS instance. Use MPI_COMM_WORLD 
  id.job=JOB_INIT; id.par=1; id.sym=0;id.comm_fortran=USE_COMM_WORLD;
  dmumps_c(&id);
  // Define the problem on the host 
  if (myid == 0) {
    id.n = n; id.nz =nz; id.irn=irn; id.jcn=jcn;
    id.a = a; //id.rhs = rhs;
  }
#define ICNTL(I) icntl[(I)-1]
// No outputs 
  id.ICNTL(1)=-1; id.ICNTL(2)=-1; id.ICNTL(3)=-1; id.ICNTL(4)=0;
// Call the MUMPS package.
  id.job=1; 
  dmumps_c(&id);
printf("Ciao"); 
  id.job=JOB_END; dmumps_c(&id);
  if (myid == 0) {
    printf("Solution is : (%8.2f  %8.2f)\n", rhs[0],rhs[1]);
  }
 // ierr = MPI_Finalize();
  return 0;

}
*/

