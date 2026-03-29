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

#include "mpi.h"
#include "lis.h"
#include <iomanip>
#include <omp.h>
#include <petscvec.h>
#include "kpp/OpenSMOKE_KPP_Dictionary.h"
#include "kpp/OpenSMOKE_KPP_DataManager.h"
#include "kpp/OpenSMOKE_KPP_ReactorNetwork.h"
#include "kpp/OpenSMOKE_KPP_Communicator.h"
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

std::ofstream fCPUTime;
std::ofstream fLog;
std::ofstream fWarning;
int globalIteration;
double 	globalNLSLinearSystemAbsoluteTolerance;
double 	globalNLSLinearSystemRelativeTolerance;
double 	globalODELinearSystemAbsoluteTolerance;
double 	globalODELinearSystemRelativeTolerance;
static char help[] = "Implementation of Petsc Vectors\n\n";

// Main
int main(int argc, char* argv[])
{
	// Disabling Bzz OpenMP
	//bzzOpenMP = 0;

	// Initialize LIS
//	lis_initialize(&argc, &argv);
	PetscInitialize(&argc, &argv, (char*)0, help);
        
        LIS_INT argc_lis = argc;
	lis_initialize(&argc_lis, &argv);
        

	// Dictionary setup
	OpenSMOKE_KPP_Dictionary dictionary;
	dictionary.ParseFile("Input.inp");

	// Read input file
	OpenSMOKE_KPP_DataManager data(dictionary);
	
	//Initialize MPI Interface
	int procrank, nprocs, MASTER;

	nprocs = MPI::COMM_WORLD.Get_size();
	procrank = MPI::COMM_WORLD.Get_rank();
	MASTER = 0;

	double timeStartTotal;
	if(procrank == 0)
	    timeStartTotal = omp_get_wtime();

	std::string filename = data.nameFolderOutput() + "/Additional/OverallTime.out";
	std::ofstream timefile(filename.c_str());
	
	// OpenMP
        omp_set_dynamic(false);
	omp_set_num_threads(data.nThreads());

	if(procrank == 0)
	{
	    std::cout << "Number of user defined threads: " << data.nThreads()           << std::endl;
	    std::cout << "Number of current threads:      " << omp_get_num_threads()     << std::endl;
	    std::cout << "Master thread id:               " << omp_get_thread_num() 	    << std::endl;
	    std::cout << "Number of processes:            " << omp_get_num_procs() 	    << std::endl;
	    std::cout << "Max number of threads:          " << omp_get_max_threads()     << std::endl;
	}


	// Gas mixture setup
	if(procrank == 0)
    	    std::cout << "Start" << std::endl; //getchar();

	OpenSMOKE_ReactingGas* mix;

	{
	    mix = new OpenSMOKE_ReactingGas[data.nThreads()];
	    for (int k=0;k<data.nThreads();k++)
		mix[k].SetupBinary(dictionary.kinetics());
	}

	/*if (data.iSaveKineticConstants() == false)
	{
		mix = new OpenSMOKE_ReactingGas[data.nThreads()];
		for (int k=0;k<data.nThreads();k++)
			mix[k].SetupBinary(dictionary.kinetics());
	}
	else
	{
		mix = new OpenSMOKE_ReactingGas[data.nThreads()];
		for (int k = 0; k < data.nThreads(); k++)
		    mix[0].SetupBinary(dictionary.kinetics());
	}*/	

	if(procrank == 0)
    	    std::cout << "Done: Mixture" << std::endl; //getchar();
	
	MPI::COMM_WORLD.Barrier();

	// Reactor network setup
	OpenSMOKE_KPP_ReactorNetwork network(mix, data, fLog, fWarning);
	OpenSMOKE_KPP_Communicator communicator(nprocs, procrank, network);
	network.SetCommunicator(&communicator);

	if(procrank == 0)
            std::cout << "Done: Network" << std::endl;

	network.ReadFirstGuess();
	if(procrank == 0)
            std::cout << "Done: FirstGuess" << std::endl;

	network.ReadTopology();
	if(procrank == 0)
            std::cout << "Done: Topology" << std::endl;

	network.BuildNetwork();
	if(procrank == 0)
            std::cout << "Network built!!!" << std::endl;


	
	// Initial Analysis
        data.SetNetworkStatus(KPP_NETWORK_STATUS_START);
	InputOutputAnalysis(*mix, network);
	MassUmbalances(network);
	
	if(procrank == 0)
	{
	    std::string log_string = data.nameFolderOutput() + "/Log.out";
	    fLog.open( log_string.c_str(), std::ios::out);
	    fLog.setf(std::ios::scientific);

	    std::string warning_string = data.nameFolderOutput() + "/Warning.out";
	    fWarning.open( warning_string.c_str(), std::ios::out);
	    fWarning.setf(std::ios::scientific);
	}

	MPI::COMM_WORLD.Barrier();

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
		globalODELinearSystemAbsoluteTolerance = 1.e-11;
	 	globalODELinearSystemRelativeTolerance = 1.e-9;
		globalNLSLinearSystemAbsoluteTolerance = 1.e-15;
	 	globalNLSLinearSystemRelativeTolerance = 1.e-12;

		int iFlagSequence;
		int iFlagODE;
		int iFlagNewtonMethod;
		double StartTime;

		for(int jj=1;jj<=data.FirstLoopIterations();jj++)
		{
			iFlagSequence	  = SequenceCSTR(network, data);

			iFlagODE          = GlobalODE(network, data);
			if(iFlagODE == 4)	break;

			globalODELinearSystemAbsoluteTolerance /= 10.;
			globalODELinearSystemRelativeTolerance /= 10.;
		}
		
//		iFlagSequence	  	= SequenceCSTR(network, data);
		iFlagNewtonMethod	= GlobalNewtonMethod(network, data);


		for(int jj=1;jj<=data.SecondLoopIterations();jj++)
		{		
			if(iFlagNewtonMethod == 4)				//Means that convergence is reached
			{
			    if(procrank == 0)
				std::cout << std::endl << "Convergence reached!" << std::endl;

			    break;
			}
			
			globalODELinearSystemAbsoluteTolerance /= 10.;
			globalODELinearSystemRelativeTolerance /= 10.;
			
			iFlagSequence	  	= SequenceCSTR(network, data);

			iFlagODE          	= GlobalODE(network, data);

			iFlagNewtonMethod	= GlobalNewtonMethod(network, data);
		}

		double timeEndTotal;

		if(procrank == 0)
		    timeEndTotal = omp_get_wtime();

		// Final Analysis
		InputOutputAnalysis(*mix, network);
		MassUmbalances(network);

		if(procrank == 0)
		{
		    timefile << "Total CPU Time: " << timeEndTotal - timeStartTotal << " seconds" << std::endl;
		    std::cout << "Total CPU Time: " << timeEndTotal - timeStartTotal << " seconds" << std::endl;
		}
	}

//	MPI::Finalize();
//	lis_finalize();
	PetscFinalize();
}

void MassUmbalances(OpenSMOKE_KPP_ReactorNetwork& network)
{
	int procrank_ = MPI::COMM_WORLD.Get_rank();

	if(procrank_ == 0)
	{
	    BzzVector umbalances(network.NumberOfReactors());
	    for(int k=1;k<=network.NumberOfReactors();k++)
		umbalances[k] = fabs(network.Mass_Umbalance()[k]);

	BzzVector massFlowIn(network.NumberOfReactors());
	for(int k=1;k<=network.NumberOfReactors();k++)
		massFlowIn[k] = fabs(network.Mass_FlowIn()[k]);

	std::cout << std::endl;
	std::cout << "// *********************************************************** // " << std::endl;
	std::cout << "//                     UMBALANCE ANALYSIS                      // " << std::endl;
	std::cout << "// *********************************************************** // " << std::endl;
	std::cout << std::setw(20) << std::left << "Max umbalance:";
	std::cout << std::setw(16) << std::left << umbalances.Max();
	std::cout << std::endl;
	std::cout << std::setw(20) << std::left << "Mean umbalance:";
	std::cout << std::setw(16) << std::left << Mean(umbalances);
	std::cout << std::endl;
	std::cout << std::setw(20) << std::left << "Max in flow:";
	std::cout << std::setw(16) << std::left << massFlowIn.Max();
	std::cout << std::endl;
	std::cout << std::setw(20) << std::left << "Min in flow:";
	std::cout << std::setw(16) << std::left << massFlowIn.Min();
	std::cout << std::endl;
	std::cout << std::setw(20) << std::left << "Mean in flow:";
	std::cout << std::setw(16) << std::left << Mean(massFlowIn);
	std::cout << std::endl;
	}

}

void InputOutputAnalysis(OpenSMOKE_ReactingGas& mix, OpenSMOKE_KPP_ReactorNetwork& network)
{
	int procrank_ = MPI::COMM_WORLD.Get_rank();
	double massFlowIn, massFlowOut;
	
	BzzVector omegaFeeds(mix.NumberOfSpecies());
	BzzVector omegaMolarFeeds(mix.NumberOfSpecies());
        BzzVector omegaMin(mix.NumberOfSpecies());
        BzzVector omegaMax(mix.NumberOfSpecies());
	BzzVector omegaOutput(mix.NumberOfSpecies());
	BzzVector omegaMolarOutput(mix.NumberOfSpecies());
	BzzVector omegaElementalFeeds(mix.NumberOfElements());
	BzzVector omegaElementalOutput(mix.NumberOfElements());

	network.ExternalFeeds(massFlowIn, omegaFeeds, omegaElementalFeeds);
	network.ExternalOutput(massFlowOut, omegaOutput, omegaElementalOutput);
        network.MinMaxMassFractions(omegaMin, omegaMax);
	
	if(procrank_ == 0)
	{

	    double MWmixFeeds = mix.GetMWFromMassFractions(omegaFeeds);
	    double MWmixOutput = mix.GetMWFromMassFractions(omegaOutput);

	    mix.GetMoleFractionsFromMassFractionsAndMW(omegaMolarFeeds, omegaFeeds, MWmixFeeds);
	    mix.GetMoleFractionsFromMassFractionsAndMW(omegaMolarOutput, omegaOutput, MWmixOutput);

	    std::cout << std::endl;
	    std::cout << "// *********************************************************** // " << std::endl;
	    std::cout << "//              IN/OUT ANALYSIS: MASS FLOW RATES               // " << std::endl;
	    std::cout << "// *********************************************************** // " << std::endl;
	    std::cout << std::setw(16) << std::left << "Input:";
	    std::cout << std::setw(16) << std::left << massFlowIn;
	    std::cout << std::endl;
	    std::cout << std::setw(16) << std::left << "Output:";
	    std::cout << std::setw(16) << std::left << massFlowOut;
	    std::cout << std::endl;
	    std::cout << std::setw(16) << std::left << "Difference(%):";
	    std::cout << std::setw(16) << std::left << (massFlowOut-massFlowIn)/massFlowIn*100.;
	    std::cout << std::endl;

	
	    std::cout << std::endl;
	    std::cout << "// *********************************************************** // " << std::endl;
	    std::cout << "//           IN/OUT ANALYSIS: SPECIES MASS FRACTIONS           // " << std::endl;
	    std::cout << "// *********************************************************** // " << std::endl;
	    for (int j=1;j<=mix.NumberOfSpecies();j++)
		if (omegaFeeds[j]>1.e-12 || omegaOutput[j]>1.e-12)		
		{
			std::cout << std::setw(16) << std::left << mix.names[j];
			std::cout << std::setw(16) << std::left << omegaFeeds[j];
			std::cout << std::setw(16) << std::left << omegaOutput[j];
			std::cout << std::endl;
		}

	    std::cout << std::endl;
	    std::cout << "// *********************************************************** // " << std::endl;
	    std::cout << "//           IN/OUT ANALYSIS: ELEMENT MASS FRACTIONS           // " << std::endl;
	    std::cout << "// *********************************************************** // " << std::endl;
	    for (int j=1;j<=mix.NumberOfElements();j++)
		if (omegaElementalFeeds[j]>1.e-12 || omegaElementalOutput[j]>1.e-12)		
		{
			std::cout << std::setw(16) << std::left << mix.list_of_elements[j-1];
			std::cout << std::setw(16) << std::left << omegaElementalFeeds[j];
			std::cout << std::setw(16) << std::left << omegaElementalOutput[j];
			std::cout << std::endl;
		}
        
            std::string fout_name = network.data().nameFolderOutput() + "/Solution/Summary.out";
            std::ofstream fout;
            fout.open(fout_name.c_str(), std::ios::out);
            fout.setf(std::ios::scientific);
        
	    fout << std::endl;
	    fout << "// ************************************************* //     ";
	    fout << "// ************************************************* //"  << std::endl;     
	    fout << "//      IN/OUT ANALYSIS: SPECIES MASS FRACTIONS      //     ";
	    fout << "//      IN/OUT ANALYSIS: SPECIES MOLE FRACTIONS      // " << std::endl;
	    fout << "// ************************************************* //     ";
	    fout << "// ************************************************* // " << std::endl;
	    for (int j=1;j<=mix.NumberOfSpecies();j++)
            {
                fout << std::setw(16) << std::left << mix.names[j];
                fout << std::setw(16) << std::left << omegaFeeds[j];
                fout << std::setw(30) << std::left << omegaOutput[j];
                fout << std::setw(16) << std::left << mix.names[j];
		fout << std::setw(16) << std::left << omegaMolarFeeds[j];
		fout << std::setw(16) << std::left << omegaMolarOutput[j];
                fout << std::endl;
            }

	    fout << std::endl;
	    fout << "// *********************************************************** // " << std::endl;
	    fout << "//           IN/OUT ANALYSIS: ELEMENT MASS FRACTIONS           // " << std::endl;
	    fout << "// *********************************************************** // " << std::endl;
	    for (int j=1;j<=mix.NumberOfElements();j++)	
            {
                fout << std::setw(16) << std::left << mix.list_of_elements[j-1];
                fout << std::setw(16) << std::left << omegaElementalFeeds[j];
                fout << std::setw(16) << std::left << omegaElementalOutput[j];
                fout << std::endl;
            }
        
	    fout << std::endl;
	    fout << "// *********************************************************** // " << std::endl;
	    fout << "//          MIN/MAX ANALYSIS: SPECIES MASS FRACTIONS           // " << std::endl;
	    fout << "// *********************************************************** // " << std::endl;
	    for (int j=1;j<=mix.NumberOfSpecies();j++)
            {
                fout << std::setw(16) << std::left << mix.names[j];
                fout << std::setw(16) << std::left << omegaMin[j];
                fout << std::setw(16) << std::left << omegaMax[j];
                fout << std::endl;
            }     
        
            fout.close();
	}
        network.CheckEnergyBalances();
}

int SequenceCSTR(OpenSMOKE_KPP_ReactorNetwork& network, OpenSMOKE_KPP_DataManager& data)
{
	int procrank = MPI::COMM_WORLD.Get_rank();
	if(procrank == 0)
	{
	    fLog << "****************************************************************************** " << std::endl; 
	    fLog << "*                              SEQUENCE                                      * " << std::endl;
	    fLog << "****************************************************************************** " << std::endl; 
	}

	int info = network.SequenceCSTR();

	fLog << std::endl << std::endl;

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
	int procrank = MPI::COMM_WORLD.Get_rank();
	if(procrank == 0)
	{
	    fLog << "****************************************************************************** " << std::endl; 
	    fLog << "*                            GLOBAL ODE                                      * " << std::endl;
	    fLog << "****************************************************************************** " << std::endl; 
	}

	// Initialize ODE System

	data.SetNetworkStatus(KPP_NETWORK_STATUS_GLOBALODE);
	network.InitializeGlobal(data.GlobalODE_SparseLinearSolver(), globalODELinearSystemAbsoluteTolerance, globalODELinearSystemRelativeTolerance);

	// Solve ODE System

	int info = network.GlobalODE();

	if(procrank == 0)
	    fLog << std::endl << std::endl;

	// Write Mass Fraction Maps
	if (info>=0)
	{
		network.WriteBackupFile(data.nameOutputBackupFile());
		network.WriteMassFractionMap();
	}

	return info;
}

int GlobalNewtonMethod(OpenSMOKE_KPP_ReactorNetwork& network, OpenSMOKE_KPP_DataManager& data)
{
	int procrank = MPI::COMM_WORLD.Get_rank();
	bool conv_index = false;

	if(procrank == 0)
	{
	    fLog << "****************************************************************************** " << std::endl; 
	    fLog << "*                            GLOBAL NLS                                      * " << std::endl;
	    fLog << "****************************************************************************** " << std::endl; 
	}
	
	// Initialize Non Linear System
	data.SetNetworkStatus(KPP_NETWORK_STATUS_GLOBALNLS);

	network.InitializeGlobal(data.GlobalNLS_SparseLinearSolver(), globalNLSLinearSystemAbsoluteTolerance, globalNLSLinearSystemRelativeTolerance);

	// Solve Non Linear System
	int info = network.GlobalNLS();

	if(procrank == 0)
	fLog << std::endl << std::endl;

	// Write Mass Fraction Maps
	{
		network.WriteBackupFile(data.nameOutputBackupFile());
		network.WriteMassFractionMap();
	}

	conv_index = network.newtonconvergence();

	if(conv_index == true) info = 4;

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

