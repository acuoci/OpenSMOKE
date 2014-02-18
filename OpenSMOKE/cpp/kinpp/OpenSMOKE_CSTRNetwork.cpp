/***************************************************************************
 *   Copyright (C) 2003-2008 by                                            *
 *   Guido Buzzi-Ferraris, Alessio Frassoldati and Alberto Cuoci		   *
 *                                                                         *
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

#include <vector>
#include "basic/OpenSMOKE_Constants.h"
#include "kinpp/OpenSMOKE_CSTRNetwork.h"
#include "kinpp/OpenSMOKE_CSTRNetwork_MyOdeSystemObjectOneCSTR.h"
#include "kinpp/OpenSMOKE_CSTRNetwork_MyOdeSystemObjectAllCSTR.h"
#include "symbolickinetics/gri12/OpenSMOKE_SymbolicKinetics_GRI12.h"
#include "kinpp/OpenSMOKE_DirectLinearSolver_Unsymmetric.h"
#include "kinpp/OpenSMOKE_PARDISO_Unsymmetric.h"

/*
#include "symbolickinetics/OpenSMOKE_SymbolicKinetics_Polimi_C1C3HTNOX_0810.h"
#include "symbolickinetics/OpenSMOKE_SymbolicKinetics_Polimi_C1C3HTNOX_AVIO.h"
#include "symbolickinetics/OpenSMOKE_SymbolicKinetics_Fluent_Glarborg152.h"
#include "symbolickinetics/OpenSMOKE_SymbolicKinetics_Fluent_DRM22_Polimi.h"
#include "symbolickinetics/OpenSMOKE_SymbolicKinetics_Fluent_DRM22_Polimi_NOX.h"
#include "symbolickinetics/OpenSMOKE_SymbolicKinetics_Fluent_DRM22_Polimi_ThermalNOX.h"
#include "symbolickinetics/OpenSMOKE_SymbolicKinetics_GRI30.h"
#include "symbolickinetics/OpenSMOKE_SymbolicKinetics_SanDiego_AVIO.h"
*/
#include "distributions/OpenSMOKE_BetaDistribution.h"
#include "distributions/OpenSMOKE_ClippedGaussianAccurateDistribution.h"
#include "distributions/OpenSMOKE_DoubleDeltaDiracDistribution.h"
#include "distributions/OpenSMOKE_SinDistribution.h"
#include "distributions/OpenSMOKE_SinIntegralDistribution.h"

const bool   iEnergyAnalysis	= true;
const double MAX_SUM_F			= 1.e-10;

const double OpenSMOKE_CSTRNetwork::R_CTOT = 0.0820578337034;
const double OpenSMOKE_CSTRNetwork::UR_CTOT = 1./R_CTOT;


int OpenSMOKE_CSTRNetwork::count		= 0;
int OpenSMOKE_CSTRNetwork::countInScope = 0;

#if SYMBOLIC_KINETICS==1
	OpenSMOKE_SymbolicKinetics* reactor[5000];
#endif

#define DEBUG_LOADING 0

void OpenSMOKE_CSTRNetwork::ErrorMessage(const string message)
{
    cout << endl;
    cout << "Class:  OpenSMOKE_CSTRNetwork"	<< endl;
    cout << "Object: " << name_object		<< endl;
    cout << "Error:  " << message           << endl;
    cout << "Press a key to continue... "   << endl;
    getchar();
    exit(-1);
}

void OpenSMOKE_CSTRNetwork::WarningMessage(const string message)
{
    cout << endl;
    cout << "Class:  OpenSMOKE_CSTRNetwork"	<< endl;
    cout << "Object: "		<< name_object	<< endl;
    cout << "Warning:  "	<< message      << endl;
    cout << "Press a key to continue... "   << endl;
    getchar();
}

void OpenSMOKE_CSTRNetwork::SetSolution(BzzMatrix &omega)
{
	massFractionsInReactorsSolution = omega;
}

int OpenSMOKE_CSTRNetwork::WhoAmI(void) const
{
	return whoAmI;
}

void OpenSMOKE_CSTRNetwork::SetTolRel(double tolr)
{
	tolRel = tolr;
}

void OpenSMOKE_CSTRNetwork::SetTolAbs(double tola)
{
	tolAbs = tola;
}

void OpenSMOKE_CSTRNetwork::SetMaxTolRel(double maxtolr)
{
	MaxTolRel = maxtolr;
}

void OpenSMOKE_CSTRNetwork::SetMaxTolAbs(double maxtola)
{
	MaxTolAbs = maxtola;
}

void OpenSMOKE_CSTRNetwork::SetMaxCountNewtonIterations(int maxCountNewtonIterations)
{
	MaxCountNewtonIterations = maxCountNewtonIterations;
}

void OpenSMOKE_CSTRNetwork::SetF1Stop(double f1Stop)
{
	F1Stop = f1Stop;
}

void OpenSMOKE_CSTRNetwork::SetOdeOnly(const bool value)
{
	OnlyODE = value;
}

void OpenSMOKE_CSTRNetwork::SetClusteringOnly(const bool value)
{
	OnlyClustering = value;
}

void OpenSMOKE_CSTRNetwork::SetMemoTemperature(void)
{
	memoTemperature = 0;	// 1 salva su file 0 memorizza
}

int OpenSMOKE_CSTRNetwork::GetNumComponents(void)
{
	return Reactions->NumberOfSpecies();
}

int OpenSMOKE_CSTRNetwork::GetNumReactions(void)
{
	return Reactions->NumberOfReactions();
}

int OpenSMOKE_CSTRNetwork::GetNumCSTRReactors(void)
{
	return numCSTRReactors;
}

void OpenSMOKE_CSTRNetwork::SetTasksPrint(void)
{
	printTasks = 1;
}

void OpenSMOKE_CSTRNetwork::SetSubTasksPrint(void)
{
	printSubTasks = 1;
}

OpenSMOKE_CSTRNetwork::~OpenSMOKE_CSTRNetwork(void)
{
	if(Reactions->NumberOfReactions() == 0 || Reactions->NumberOfSpecies() == 0)
		return;
}

void OpenSMOKE_CSTRNetwork::SetMaxCorrectionCoefficient(const double _MaxCoeffCorr)
{
	MaxCoeffCorr = _MaxCoeffCorr;
}

void OpenSMOKE_CSTRNetwork::SetDeltaTFluctuationsMaxDelta(const double _Fluctuations_DeltaMax)
{
		Fluctuations_DeltaMax	= _Fluctuations_DeltaMax;
}

void OpenSMOKE_CSTRNetwork::SetDeltaTFluctuationsMaxUserDefined(const double _Fluctuations_TMax)
{
	Fluctuations_TMax		= _Fluctuations_TMax;
}

void OpenSMOKE_CSTRNetwork::SetDeltaTFluctuationsMaxLocal(const double _Fluctuations_CcMax)
{
	Fluctuations_CcMax		= _Fluctuations_CcMax;
}

void OpenSMOKE_CSTRNetwork::SetFluctuationsList(const string speciesFluctuations)
{
	ifstream fInput;
	openInputFileAndControl(fInput, speciesFluctuations);

	for(;;)
	{
		string dummy;
		fInput >> dummy;
		if (dummy == "#END")
			break;
		list_of_fluctuating_species.push_back(dummy);
	}

	fInput.close();
}

OpenSMOKE_CSTRNetwork::OpenSMOKE_CSTRNetwork(void)
{
	count++;
	countInScope++;

	whoAmI				= count;

	
	numCSTRReactors		= 1;
	sinExpansion		= 16;
	MaxCoeffCorr		= 1.e10;
	
	Fluctuations_DeltaMax	= 0.;
	Fluctuations_TMax		= 0.;
	Fluctuations_CcMax		= 0.;

	name_object = "[Name not assigned]";

	MaxCountNewtonIterations	= 10;
	F1Stop						= 1.e-6;

	tolRel						= 1.e-2;
	tolAbs						= 1.e-8;
	MaxTolRel					= 1.e-5;
	MaxTolAbs					= 1.e-9;

	OnlyODE	= false;
	OnlyClustering	= false;

	cout.setf(ios::scientific);
}

void OpenSMOKE_CSTRNetwork::AssignKineticScheme(OpenSMOKE_ReactingGas &_Reactions)
{
    Reactions		= &_Reactions;
	numReactions	= Reactions->NumberOfReactions();
	numComponents	= Reactions->NumberOfSpecies();
}

void OpenSMOKE_CSTRNetwork::CheckInputfile(ifstream &file, string fileName)
{
	if (file.fail())
		ErrorMessage("It was impossible to correctly read some values in the file " + fileName);
}

void OpenSMOKE_CSTRNetwork::reading_the_network_file(string fileNetwork, int &countInput)
{
	int			i,j,k;
	ifstream	fInput;

	BzzVectorInt	bOut,bki;
	BzzVectorInt	countCSTRConnect;
	BzzVector	massExternalOutputCFDCells;
	BzzVector	massExternalInputCFDCells;
	BzzVector	massInputCFDCells;
	BzzVector	massOutputCFDCells;

	int indexCell;
	int nInput, nOutput;
	int numberOfSpeciesCFD;
	int indexSpecies;
	int destination_cell;
	double omegaSpecies;
	double temperatureCell, pressureCell, volumeCell, TvarianceCell;
	double inputMassFlowRate, outputMassFlowRate;
	double outputDiffusionRate;

	// Opening the network file
	openInputFileAndControl(fInput, fileNetwork);

	// Reading the number of original CFD cells
	fInput >> numCSTRReactors;
	CheckInputfile(fInput, fileNetwork);

	// Memory Allocation loacal variables
	ChangeDimensions(numCSTRReactors, &bki);
	ChangeDimensions(numCSTRReactors, &bOut);
	ChangeDimensions(numCSTRReactors, &countCSTRConnect);
	ChangeDimensions(numCSTRReactors, &massExternalOutputCFDCells);
	ChangeDimensions(numCSTRReactors, &massExternalInputCFDCells);
	ChangeDimensions(numCSTRReactors, &massInputCFDCells);
	ChangeDimensions(numCSTRReactors, &massOutputCFDCells);


	// Memory Allocation Global Variables
	ChangeDimensions(numCSTRReactors, &temperature);
	ChangeDimensions(numCSTRReactors, &qanm);
	ChangeDimensions(numCSTRReactors, &pressure);
	ChangeDimensions(numCSTRReactors, &volume);

	// Connection Vectors - Initialization
	bki = 1;
	bOut = 1;

	// Number of CFD cells which are connected to an external input
	countInput = 0;

	// Loop on each CFD cell
	for(k = 1;k <= numCSTRReactors;k++)
	{
		#if DEBUG_LOADING==1
			cout << "-----------------------------------------------" << endl;
			cout << "           Loading reactor #" << k               << endl;
			cout << "-----------------------------------------------" << endl;
		#endif

		#if DEBUG_LOADING==1
			cout << " * main data (T, P, V, sigma)" << endl;
		#endif
		
		fInput >> indexCell;
		fInput >> temperatureCell;
		fInput >> pressureCell;
		fInput >> volumeCell;
		fInput >> TvarianceCell;
		CheckInputfile(fInput, fileNetwork);

		// The cells must be ordered
		//if(k != indexCell)	
		//	ErrorMessage("Wrong CSTR sequence");

		// Checking if the cell numbers are correct
		if(indexCell < 1 || indexCell > numCSTRReactors)
			ErrorMessage("Reactor out of range!");

		temperature[indexCell] = temperatureCell;
		pressure[indexCell] = pressureCell;
		volume[indexCell] = volumeCell;
		qanm[indexCell] = TvarianceCell;

		// Number of external input in the CFD cell
		// ------------------------------------------------------------------------------
		#if DEBUG_LOADING==1
			cout << " * input stream from the external environment" << endl;
		#endif

		fInput >> nInput;
		CheckInputfile(fInput, fileNetwork);
		if(nInput != 0)
		{
			countInput++;

			for(j=1;j<=nInput;j++)
			{
				#if DEBUG_LOADING==1
					cout << "   - input #" << j << endl;
				#endif

				fInput >> inputMassFlowRate;
				CheckInputfile(fInput, fileNetwork);

				fInput >> numberOfSpeciesCFD;
				CheckInputfile(fInput, fileNetwork);

				for(i = 1;i <= numberOfSpeciesCFD;i++)
				{
					fInput >> indexSpecies;
					fInput >> omegaSpecies;
					CheckInputfile(fInput, fileNetwork);
				}
			}
		}

		// Number of output (internal) in the CFD cell
		// ------------------------------------------------------------------------------
		#if DEBUG_LOADING==1
			cout << " * internal (reactor<->reactor) streams" << endl;
		#endif

		fInput >> nOutput;
		CheckInputfile(fInput, fileNetwork);
		for(i = 1;i <= nOutput;i++)
		{
			#if DEBUG_LOADING==1
				cout << "   - stream #" << i << endl;
			#endif

			fInput >> destination_cell;
			fInput >> outputMassFlowRate;
			fInput >> outputDiffusionRate;
			CheckInputfile(fInput, fileNetwork);

			// Connection number of each CFD cell
			if(destination_cell > 0)
			{
				countCSTRConnect[indexCell]++;
				countCSTRConnect[destination_cell]++;
			}
		}

		#if DEBUG_LOADING==1
			cout << endl;
		#endif
	}
	fInput.close();

	#if DEBUG_LOADING==1
		cout << "CFDNetwork.bzz file correctly loaded!" << endl;
		cout << "Press enter to continue..." << endl;
		getchar();
	#endif

	// Min and Max Temperature
	TminGlobal = temperature.Min();
	TmaxGlobal = temperature.Max();

	// Memory allocation for the sparse matrix of connections
	cstrConnect(countCSTRConnect);

	// Saving information on CFD Network
	BzzSave saveNetworkInput("Temp/CFDNetworkInput.tmp");
	saveNetworkInput << countInput;

	// Opening the network file
	openInputFileAndControl(fInput, fileNetwork);

	// Reading the number of original CFD cells
	fInput >> numCSTRReactors;
	CheckInputfile(fInput, fileNetwork);


	// Cycle on every original CFD cell
	for(k = 1;k <= numCSTRReactors;k++)
	{
		fInput >> indexCell;
		fInput >> temperatureCell;
		fInput >> pressureCell;
		fInput >> volumeCell;
		fInput >> TvarianceCell;
		CheckInputfile(fInput, fileNetwork);

		// Number of external input in the CFD cell
		// ------------------------------------------------------------------------------
		fInput >> nInput;
		CheckInputfile(fInput, fileNetwork);
		if(nInput != 0)
		{
			double totalInputMassFlowRate = 0.;
			BzzVector speciesInletFeed(numComponents);

			saveNetworkInput << indexCell;

			for(j=1; j<=nInput; j++)
			{
				fInput >> inputMassFlowRate;
				CheckInputfile(fInput, fileNetwork);

				fInput >> numberOfSpeciesCFD;
				CheckInputfile(fInput, fileNetwork);

				totalInputMassFlowRate += inputMassFlowRate;
				if(j == nInput)
					saveNetworkInput << totalInputMassFlowRate << numberOfSpeciesCFD;

				double omegaSum = 0.;
				for(i = 1;i <= numberOfSpeciesCFD;i++)
				{
					fInput >> indexSpecies;
					fInput >> omegaSpecies;
					CheckInputfile(fInput, fileNetwork);

					omegaSum += omegaSpecies;
					speciesInletFeed[indexSpecies] += inputMassFlowRate * omegaSpecies;

					if(j == nInput)
						saveNetworkInput << indexSpecies << speciesInletFeed[indexSpecies];
				}

				if(omegaSum < .999999 || omegaSum > 1.000001)
					cout << "WARNING: Wrong fraction in " << indexCell << " reactor: " << omegaSum << endl;
			}

			massExternalInputCFDCells[indexCell] = totalInputMassFlowRate;
		}

		// Number of output (internal) in the CFD cell
		// ------------------------------------------------------------------------------
		fInput >> nOutput;
		CheckInputfile(fInput, fileNetwork);
		for(i = 1;i <= nOutput;i++)
		{
			fInput >> destination_cell;
			fInput >> outputMassFlowRate;
			fInput >> outputDiffusionRate;
			CheckInputfile(fInput, fileNetwork);

			// In this case the cell is internal
			if(destination_cell > 0)
			{
				cstrConnect(indexCell,bki[indexCell]) = destination_cell;
				cstrConnect(destination_cell,bki[destination_cell]) = indexCell;
				bki[indexCell]++;
				bki[destination_cell]++;

				massOutputCFDCells[indexCell]		+=outputMassFlowRate;
				massInputCFDCells[destination_cell]	+=outputMassFlowRate;
			}

			// In this case the cell is external
			else if(destination_cell == 0)
			{
				massExternalOutputCFDCells[indexCell] = outputMassFlowRate;
			}
		}
	}
	saveNetworkInput.End();
	fInput.close();

	ofstream fOutput;
	double meanError=0.;
	double maxError=0.;
	openOutputFileAndControl(fOutput, "Temp/MassFluxesCFDNetwork.tmp");
	for(k = 1;k <= numCSTRReactors;k++)
	{
		double mIN  = massExternalInputCFDCells[k] + massInputCFDCells[k];
		double mOUT = massExternalOutputCFDCells[k] + massOutputCFDCells[k];
		if (mIN+mOUT<1.e-16)
			continue;
		double relative_error = fabs((mOUT-mIN)/(mIN+mOUT));
		meanError+=relative_error;
		if (relative_error>maxError)
			maxError = relative_error;

		fOutput	<< k									<< "\t"
				<< mIN									<< "\t"
				<< mOUT									<< "\t"
				<< relative_error						<< "\t"
				<< massExternalInputCFDCells[k]			<< "\t"
				<< massExternalOutputCFDCells[k]		<< "\t"
				<< massInputCFDCells[k]					<< "\t"
				<< massOutputCFDCells[k]				<< "\t"
				<< endl;
	}

	meanError /= numCSTRReactors;
	cout << "        Mean relative error in mass balances on CFD cells: " << meanError << endl;
	cout << "        Max  relative error in mass balances on CFD cells: " << maxError << endl;

}

void OpenSMOKE_CSTRNetwork::read_tolerances(const string fileTolerances, int numCSTRReactors, double cicloCluster)
{
	int i;

	BzzLoad loadTol(fileTolerances);
	loadTol >> tolRelCBase >> tolAbsCBase >> numTolT;

	ChangeDimensions(numTolT, &tInf);
	ChangeDimensions(numTolT, &tSup);
	ChangeDimensions(numTolT, &dt);
	ChangeDimensions(numTolT, &dtBase);

	for(i = 1;i <= numTolT;i++)
		loadTol >> tInf(i) >> tSup(i) >> dtBase(i);

	loadTol.End();

	tolRelC = cicloCluster * tolRelCBase;
	tolAbsC = cicloCluster * tolAbsCBase;
	for(i = 1;i <= numTolT;i++)
		dt[i] = cicloCluster * dtBase(i);
}

void OpenSMOKE_CSTRNetwork::read_first_guess(string first)
{
	BzzLoad load(first);
	load >> numComponents >> iSpec >> massFractionsInReactors;
	load.End();
	numComponentsReduced = iSpec.Size();
}

void OpenSMOKE_CSTRNetwork::build_cluster(double cicloCluster, int numCSTRReactors, int countInput, int &numCluster, BzzVectorInt &cstrClusterSize)
{
	cout << "I am building the cluster... " << endl;

	int		i,j,k;
	int		ki, kii;
	int		iBoss, iLast, iTry, kiSize;
	int		clusterOK;
	int		totalNumCSTRInCluster;
	double	tBoss, dtBoss;
	int		maxCluster;
	int		divisorForMaxCluster;

	BzzVectorInt	cstrCluster,
					cstrConnectVector,
					auxil;

	BzzVector massFractionsBoss,
					massFractionsTry;


	// Initialize values
	divisorForMaxCluster = 10;
	totalNumCSTRInCluster = 0;
	maxCluster = numCSTRReactors / divisorForMaxCluster;
	numCluster = 0;
	if(cicloCluster == 0)	maxCluster = 1;

	// Initialize values
	BzzSave saveCluster("Temp/NetworkCluster.tmp");
	ChangeDimensions(numCSTRReactors,&giveClusterIndexFromCellIndex);

	// Starting clustering algorithm
	double start = BzzGetCpuTime();

	// Loop on every CFD cell
	for(ki = 1;ki <= numCSTRReactors;ki++)
	{
		kiSize = cstrConnect.GetSize(ki);
		if(kiSize == 0)
			continue;

		numCluster++;
		ChangeDimensions(numCSTRReactors,&cstrCluster);
		cstrConnect.GetVector(ki,&cstrConnectVector);
		iLast = 0;
		iBoss = ki;
		cstrCluster[++iLast] = iBoss;

		tBoss = temperature(iBoss);
		massFractionsInReactors.UseMatrixRowAsVector(iBoss,&massFractionsBoss);
		if(tBoss < tSup(1))
			dtBoss = dt(1);
		else if(tBoss > tInf(numTolT))
			dtBoss = dt(numTolT);
		else
		{
			for(i = 1;i <= numTolT;i++)
			{
				if(tBoss >= tInf(i) && tBoss <= tSup(i))
				{
					dtBoss = dt[i];
					break;
				}
			}
		}

		for(j = 1;j <= kiSize;j++)
		{
			iTry = cstrConnectVector[j];
			if(iTry == iBoss)
				continue;
			if(cstrConnect.GetSize(iTry) == 0)
				continue;
			clusterOK = 0;
			if(fabs(tBoss - temperature[iTry]) < dtBoss)
				clusterOK = 1;
			if(clusterOK == 1) // verifica la composizione
			{
				massFractionsInReactors.UseMatrixRowAsVector(iTry,&massFractionsTry);
				for(i = 1;i <= numComponentsReduced;i++)
				{
					if(fabs(massFractionsTry[i] - massFractionsBoss[i]) >
						tolRelC * massFractionsBoss[i] + tolAbsC)
						{
							clusterOK = 0;
							break;
						}
				}
			}

			if(iLast >= maxCluster)
				break;

			// verifca se metterlo nel cluster
			if(clusterOK == 1) // lo deve mettere
				cstrCluster[++iLast] = iTry;
		}

		cstrConnect.DeleteVector(iBoss);
		for(k = 1;k <= iLast;k++)
		{
			kii = cstrCluster[k];
			kiSize = cstrConnect.GetSize(kii);
			if(kiSize == 0)
				continue;
			cstrConnect.GetVector(kii,&cstrConnectVector);
			for(j = 1;j <= kiSize;j++)
			{
				iTry = cstrConnectVector[j];
				if(iTry == iBoss)
					continue;

				if(cstrConnect.GetSize(iTry) == 0)
					continue;

				clusterOK = 0;

				if(fabs(tBoss - temperature[iTry]) >= dtBoss)
					continue;

				clusterOK = 1;

				for(i = 1;i <= iLast;i++)
				{
					if(iTry == cstrCluster[i])
					{
						clusterOK = 0;
						break;
					}
				}

				if(clusterOK == 1) // verifica la composizione
				{
					massFractionsInReactors.UseMatrixRowAsVector(iTry,&massFractionsTry);
					for(i = 1;i <= numComponentsReduced;i++)
					{
						if(fabs(massFractionsTry[i] - massFractionsBoss[i]) >
							tolRelC * massFractionsBoss[i] + tolAbsC)
							{
								clusterOK = 0;
								break;
							}
					}
				}

				if(iLast >= maxCluster)
					break;

				// verifica se metterlo nel cluster
				if(clusterOK == 1) // lo deve mettere
					cstrCluster[++iLast] = iTry;
			}

			cstrConnect.DeleteVector(kii);
		}

		auxil.GetBzzVector(iLast,1,cstrCluster);
		totalNumCSTRInCluster += iLast;
		cstrClusterSize.Append(auxil.Size());
		for(i = 1;i <= iLast;i++)
			giveClusterIndexFromCellIndex[auxil[i]] = numCluster;
		saveCluster << auxil;
	}

	saveCluster.End();

	cout << "I built the cluster... " << endl;
}

void OpenSMOKE_CSTRNetwork::calculate_initial_massfractions_in_clusters(int numCluster, BzzMatrix &omegaReduced_Cluster, BzzMatrix &omegaReduced_CFD)
{
	int i, k, j;
	int ispecSize;
	int numC;

	BzzVector clusterVolume(numCluster);

	cout << "Calculate initial massfractions in clusters..." << endl;

	BzzLoad load(firstGuessFile);
	load >> numC >> iSpec >> omegaReduced_CFD;
	ispecSize = iSpec.Size();
	load.End();

	ChangeDimensions(numCluster, ispecSize, &omegaReduced_Cluster);

	for(i = 1;i <= numCSTRReactors;i++)
	{
		j = giveClusterIndexFromCellIndex[i];
		clusterVolume[j] += volume[i];
		for(k = 1;k <= ispecSize;k++)
			omegaReduced_Cluster[j][k] += volume[i] * omegaReduced_CFD[i][k];
	}

	for(i = 1;i <= numCluster;i++)
		for(k = 1;k <= ispecSize;k++)
			omegaReduced_Cluster[i][k] /= clusterVolume[i];

	BzzVector sumR;
	omegaReduced_Cluster.GetRowsSum(&sumR);
	for(i = 1;i <= numCluster;i++)
		for(k = 1;k <= ispecSize;k++)
			omegaReduced_Cluster[i][k] /= sumR[i];

	cout << "Calculated initial massfractions in clusters..." << endl;
	cout << "Number of reactors: " << numCluster << endl;
}

void OpenSMOKE_CSTRNetwork::calculate_clustering(string fileNetwork, int cicloDiffusion, int numCluster,
											BzzVectorInt &cstrClusterSize,
											BzzMatrix &omegaReduced_Cluster,
											BzzMatrix &omegaReduced_CFD)
{
	cout << "Calculate clustering..." << endl;

	int i,j,k;
	ifstream fInput;

	int indexCell,destinationCell;
	double temperatureCell, pressureCell, volumeCell, TvarianceCell;
	double inputMassFlowRate, outputMassFlowRate, outputDiffusionRate;
	double omegaSpecies;
	int indexSpecies;
	int numberOfSpeciesCFD; 
	int nInput, nOutput;
	int is,isMax;
	double maxSpec,mspec;

	double diffusionCorrector;
	diffusionCorrector = 1. / BzzPow2(double(cicloDiffusion));

	openInputFileAndControl(fInput, fileNetwork);
	fInput >> numCSTRReactors;
	CheckInputfile(fInput, fileNetwork);

	originalNumCSTRReactors = numCSTRReactors;
	numCSTRReactors = numCluster;

	ChangeDimensions(numCSTRReactors,numComponents, &feed);
	ChangeDimensions(numCSTRReactors, numCSTRReactors,&Ms);
	ChangeDimensions(numCSTRReactors, numCSTRReactors,&Mg);

	ChangeDimensions(numCSTRReactors, &temperature);
	ChangeDimensions(numCSTRReactors, &qanm);
	ChangeDimensions(numCSTRReactors, &pressure);
	ChangeDimensions(numCSTRReactors, &volume);
	ChangeDimensions(numCSTRReactors, &massInput);
	ChangeDimensions(numCSTRReactors, &massOutput);
	ChangeDimensions(numCSTRReactors, &sM);
	ChangeDimensions(numCSTRReactors, &logTm);
	ChangeDimensions(numCSTRReactors, &loguRTm);
	ChangeDimensions(numCSTRReactors, &cTotm);
	ChangeDimensions(numCSTRReactors, &massRate);

	ChangeDimensions(numCSTRReactors, &Tk_20000);
	ChangeDimensions(numCSTRReactors, &Tk_40000);
	ChangeDimensions(numCSTRReactors, &Tk_60000);

	BzzVector volumeCellOriginal(originalNumCSTRReactors);

	Mg.SetDiagonal(0,1.);


	// Vettori di connessione dei flussi massivi
	FromClusterToCluster_Index				= new BzzVectorInt[numCSTRReactors+1];
	FromClusterToCluster_MassFlowRate		= new BzzVector[numCSTRReactors+1];
	FromClusterToCluster_DiffusionFlowRate	= new BzzVector[numCSTRReactors+1];
	ChangeDimensions(numCSTRReactors, &Hinput);

	// Ciclo su tutte le celle della simulazione CFD originale
	for(k = 1;k <= originalNumCSTRReactors;k++)
	{
		fInput >> indexCell;
		fInput >> temperatureCell;
		fInput >> pressureCell;
		fInput >> volumeCell;
		fInput >> TvarianceCell;
		CheckInputfile(fInput, fileNetwork);

		// Recupero indice cluster a cui appartiene la cella
		int indexCluster = giveClusterIndexFromCellIndex[indexCell];

		// Calcolo delle proprieta' medie del cluster
		temperature[indexCluster] += temperatureCell;
		pressure[indexCluster] += pressureCell;
		volume[indexCluster] += volumeCell;
		qanm[indexCluster] += TvarianceCell;

		fInput >> nInput;
		CheckInputfile(fInput, fileNetwork);

		volumeCellOriginal[k] = volumeCell;
		// CUOCI
		// Nel caso in cui ci siano degli input esterni alla cella di calcolo
		if(nInput != 0)
		{
			for(j = 1;j <= nInput;j++)
			{
				BzzVector omega_input(numComponents);

				fInput >> inputMassFlowRate;
				CheckInputfile(fInput, fileNetwork);

				fInput >> numberOfSpeciesCFD;
				CheckInputfile(fInput, fileNetwork);

				massInput[indexCluster] += inputMassFlowRate;

				double sumOmega = 0.;
				for(i = 1;i <= numberOfSpeciesCFD;i++)
				{
					fInput >> indexSpecies;
					fInput >> omegaSpecies;
					CheckInputfile(fInput, fileNetwork);

					omega_input[indexSpecies] = omegaSpecies;

					sumOmega += omegaSpecies;
				}

				if(sumOmega < .999999 || sumOmega > 1.000001)
					cout << "WARNING: Wrong fraction in " << indexCell << " reactor: " << sumOmega << endl;
				
				// Enthalpy
				BzzVector h_input(numComponents);
				Reactions->SpeciesEnthalpy(temperatureCell);								// [-]
				Product((Constants::R_J_kmol*temperatureCell), Reactions->h,  &h_input);	// [J/kmol]
				ElementByElementProduct(h_input,  Reactions->uM, &h_input);				// [J/kg]
				Hinput[indexCluster] += Dot(h_input, omega_input)*inputMassFlowRate;	// [J/s]
			}
		}

		// Numero di celle collegate alla cella in esame
		fInput >> nOutput;
		CheckInputfile(fInput, fileNetwork);

		// Ciclo su tutte le celle vicine
		for(i = 1;i <= nOutput;i++)
		{
			fInput >> destinationCell;
			fInput >> outputMassFlowRate;
			fInput >> outputDiffusionRate;
			CheckInputfile(fInput, fileNetwork);

			// Recupero indice del cluster in cui si trova la cella di destinazione
			int destinationCluster = giveClusterIndexFromCellIndex[destinationCell];

			// Correzione dei coefficienti diffusivi e costruzione delle matrici dei flussi
			// convettivi e diffusivi
			if(destinationCell > 0)
			{
				if(cicloDiffusion == 4)
					outputDiffusionRate = 0.;
				else
				{
					// Viene calcolata la massima differenza tra le frazioni massive delle specie
					// della simulazione CFD originale tra il cluster in cui si trova la cella
					// principale e il cluster in cui si trova quella di destinazione
					maxSpec = 0.;
					for(is = 1;is <= iSpec.Size();is++)
					{
						mspec = fabs(omegaReduced_Cluster[indexCluster][is] - omegaReduced_Cluster[destinationCluster][is]);
						if(mspec > maxSpec)
						{
							maxSpec = mspec;
							isMax = is;
						}
					}

					// Nel caso in cui la differenza sia maggiore di zero (questo non puo' verificarsi)
					// se la cella principale e quella di destinazione sono nello stesso cluster
					// Il coefficiente diffusivo viene corretto
					if(maxSpec > 0.)
					{
						mspec = outputDiffusionRate * fabs(omegaReduced_CFD[destinationCell][isMax] - omegaReduced_CFD[indexCell][isMax]) / maxSpec;
						outputDiffusionRate = mspec * diffusionCorrector;
					}
				}

				// Nella matrice Mg vengono inseriti i flussi massivi uscenti da ciascun
				// cluster: il cluster i riceve la portata massiva outputMassFlowRate dal
				// cluster j - Tutti gli elementi sono positivi perche' si tratta dei flussi uscenti
				// (o entranti a seconda del punto di vista)
				Mg(destinationCluster,indexCluster) += outputMassFlowRate;

				// Sulla diagonale vengono sommati tutti i flussi diffusivi normalizzati
				// che interessano il cluster in esame e che sono pari alla somma di tutti quelli
				// che interessano le celle che li compongono
				Ms(indexCluster,indexCluster)				+= outputDiffusionRate;
				Ms(destinationCluster,destinationCluster)	+= outputDiffusionRate;

				// Al di fuori della diagonale vengono inseriti i flussi diffusivi
				// normalizzati cambiati di segno; se eventualmente i due cluster
				// coincidono vengono ad essere eliminati i termini aggiunti sopra
				Ms(indexCluster,destinationCluster) -= outputDiffusionRate;
				Ms(destinationCluster,indexCluster) -= outputDiffusionRate;

				// TODO
				sM(indexCluster)		+= outputDiffusionRate;
				sM(destinationCluster)	+= outputDiffusionRate;
			}

			// Nel caso in cui la cella di destinazione sia un outflow viene aggiornata
			// la portata massiva uscente totale dal cluster
			else if(destinationCell == 0)
			{
				massOutput[indexCluster] += outputMassFlowRate;
			}

			// Vettori di connessione portate massive convettive
			if(destinationCluster != indexCluster)
			{
				if(destinationCell == 0)	FromClusterToCluster_Index[indexCluster].Append(0);
				else						FromClusterToCluster_Index[indexCluster].Append(destinationCluster);
				FromClusterToCluster_MassFlowRate[indexCluster].Append(outputMassFlowRate);
				FromClusterToCluster_DiffusionFlowRate[indexCluster].Append(outputDiffusionRate);
			}
		}
	}
	fInput.close();

	cout << "Calculated clustering..." << endl;

	// Mean Values for each cluster
	cout << "Mean values..." << endl;
	double uNumberOfCFDCells;
	for(i = 1;i <= numCSTRReactors;i++)
	{
		uNumberOfCFDCells	 = 1. / double(cstrClusterSize[i]);
		temperature[i]		*= uNumberOfCFDCells;
		qanm[i]				*= uNumberOfCFDCells;
		pressure[i]			*= uNumberOfCFDCells;
	}

	ProductT(volume, Reactions->M, &volumeMolecularWeightT);

	cout << "Mean values done..." << endl;

	if (OnlyClustering == true)
	{
		cout << "I am writing ClusteringTopology.out" << endl;
		BzzSave fClusteringTopology("ClusteringTopology.out");
		fClusteringTopology << giveClusterIndexFromCellIndex;
		fClusteringTopology << volumeCellOriginal;
		fClusteringTopology.End();
	}
}

void OpenSMOKE_CSTRNetwork::calculate_massflowrate_in_reactors(BzzVector &inputMassFlowRate,
														  BzzVector &outputMassFlowRate)
{
	int i, j, k;
	ElementBzzMatrixSparse *elem;

	// Calcolo della portata massiva complessivamente entrante nel reattore
	inputMassFlowRate	= massInput;
	for(i = 1;i <= numCSTRReactors;i++)
	{
		elem = Mg.GetStartingElementInRow(i);
		while(elem)
		{
			j = elem->column;
			if(i != j)
				inputMassFlowRate[i] += elem->value;
			elem = elem->next;
		}
	}

	// Calcolo della portata massiva complessivamente uscente dalreattore
	outputMassFlowRate	= massOutput;
	Transpose(&Mg);
	for(i = 1;i <= numCSTRReactors;i++)
	{
		elem = Mg.GetStartingElementInRow(i);
		while(elem)
		{
			j = elem->column;
			if(i != j)
				outputMassFlowRate[i] += elem->value;
			elem = elem->next;
		}
	}

	// Calcolo dell'errore massimo e medio (relativo) nei bilanci massivi sui singoli reattori
	double meanError=0.;
	double maxError=0.;
	for(k = 1;k <= numCSTRReactors;k++)
	{
		double mIN  = inputMassFlowRate[k];
		double mOUT = outputMassFlowRate[k];
		if (mIN+mOUT<1.e-16)
			continue;

		double relative_error = fabs((mOUT-mIN)/(mIN+mOUT));
		meanError+=relative_error;
		if (relative_error>maxError)
			maxError = relative_error;
	}

	meanError /= numCSTRReactors;
	cout << "        Mean relative error in mass balances on clusters: " << meanError << endl;
	cout << "        Max  relative error in mass balances on clusters: " << maxError << endl;
}

void OpenSMOKE_CSTRNetwork::complete_clustering(int countInput, int cicloDiffusion, int relaxation, int iaia, BzzMatrix &omegaReduced_Cluster)
{
	int i, j, k;
	double *ptrVal;
	double val;
	ElementBzzMatrixSparse *elem;
	BzzVector provisionalInputMassRate;
	BzzVector provisionalOutputMassRate;

	// ----------------------------------------------------------------------------------------------
	// KinPP Summary
	//-----------------------------------------------------------------------------------------------
	cout << endl << endl;
	cout << "------------------------------------------------------------------------------"	<< endl;
	cout << "                      KINETIC POST PROCESSOR (KinPP)                          "	<< endl;
	cout << "                         Version: September 2009	                           "	<< endl;
	cout << "                                 MORE	                                       "	<< endl;
	cout << "------------------------------------------------------------------------------"	<< endl;	
	cout << " - Total number of reactors:  " << numCSTRReactors									<< endl;
	cout << " - Original number of cells:  " << originalNumCSTRReactors							<< endl;
	cout << " - Clustering ratio [%]:      " << double(numCSTRReactors)/
												double(originalNumCSTRReactors)*100.			<< endl;
	cout << " - Number of species:         " << numComponents									<< endl;
	cout << " - Number of variables:       " << numComponents*numCSTRReactors					<< endl;
	cout << " - Jacobian size [MB]:        " << BzzPow2(numComponents*numCSTRReactors)*8./1.e6	<< endl;
	cout << endl;

	// Calculate provisional mass flow rate
	cout << " - 01 - Relative Errors before normalization " << endl;
	calculate_massflowrate_in_reactors(provisionalInputMassRate, provisionalOutputMassRate);

	// Vengono normalizzate le portate massive della matrice Mg e ne viene cambiato il segno e viene rifatta la trasposizione;
	// i termini sulla diagonale principale vengono settati pari al valore unitario
	cout << " - 02 - Normalizing mass balances on clusters... " << endl;
	cout << "        Building Mg matrix... " << endl;
	for(i = 1;i <= numCSTRReactors;i++)
	{
		elem = Mg.GetStartingElementInRow(i);
		while(elem)
		{
			j = elem->column;
			if(i != j)
				elem->value = -elem->value / provisionalOutputMassRate[i];
			elem = elem->next;
		}
	}
	Transpose(&Mg);
	Mg.SetDiagonal(0,1.);

	// Viene normalizzato anche il flusso massivo uscente
	cout << "        Output mass flow rate normalization... " << endl;
	for(i = 1;i <= numCSTRReactors;i++)
	{
		if(provisionalOutputMassRate[i] != 0.)
			massOutput[i] /=  provisionalOutputMassRate[i];
	}

	// Vengono calcolati i nuovi flussi convettivi che interessano ciscun cluster
	// vengono inseriti nel vettore massRate
	cout << "        Updating input mass flow rates... " << endl;
	if(relaxation == 0 || relaxation == 2 || (relaxation == 1 && cicloDiffusion == 4))
	{
		bool BzzLinearSolvers = false;
		bool PardisoLinearSolvers = true;

		if (BzzLinearSolvers == true)
		{ 
			cout << "I am using BzzMath..." << endl;
			BzzFactorizedSparseGauss Fg;		
			Fg = Mg;
			massRate = massInput;
			Solve(&Fg,&massRate);
			BzzSave massRateFile('*', "Temp/MassRate.tmp");		
			massRateFile << massRate;
			massRateFile.End();
		}
		else if (PardisoLinearSolvers == true)
		{
			cout << "I am using Pardiso..." << endl;
			OpenSMOKE_PARDISO_Unsymmetric Fg(OPENSMOKE_DIRECTSOLVER_SQUAREMATRIX);
			Fg.SetSparsityPattern(Mg);
			Fg.UpdateMatrix(Mg);
			Fg.NumericalFactorization();
			Fg.Solve(massInput, massRate);
			Fg.Delete();
			BzzSave massRateFile('*', "Temp/MassRate.tmp");		
			massRateFile << massRate;
			massRateFile.End();
		}
	}
	else
	{
		BzzLoad massRateFile('*', "Temp/MassRate.tmp");
		massRateFile >> massRate;
		massRateFile.End();
	}

	// Nella matrice Ms e nel vettore sM vengono aggiunti i flussi convettivi
	// La matrice Ms(i,j) e[ definita nella maniera opposta rispetto alla matrice Mg,
	// cioe' come la portata massiva che dal reattore i va a finire nel reattore j
	// Quindi la somma degli elementi lungo una riga e' la portata massiva uscente dal
	// corrispondente reattore; la somma lungo una colonna e' la portata massiva entrante
	cout << "        Updating Mg and Ms matrices... " << endl;
	Mg.BeginScanning();
	while(ptrVal = Mg.Scanning(&i,&j,&val))
	{
		if(i == j)
		{
			sM[i]	+= massRate[i];
			Ms(i,i)	+= massRate[i];
			Mg(i,i)	 = massRate[i];
		}
		else
		{
			Ms(i,j) += val*massRate[j];
			Mg(i,j)  = val*massRate[j];
		}
	}

	for(i = 1;i <= numCSTRReactors;i++)
		massOutput[i] *= massRate[i];

	// Calculate provisional mass flow rate
	// The transposition is necessary because the Ms matrix is defined in the opposite way
	cout << " - 03 - Relative Errors after normalization " << endl;
	Transpose(&Mg);
	calculate_massflowrate_in_reactors(provisionalInputMassRate, provisionalOutputMassRate);
	Delete(&Mg);


	// Updating External Feeds
	cout << " - 04 - Updating external feeds... " << endl;
	{
		int indexCell;
		double massIn;
		int numOriginalSpecies, iCluster;
		BzzVector massFeed(numCSTRReactors);
		
		ChangeDimensions(numCSTRReactors,numComponents,&feed);
		BzzLoad loadClusterInput("Temp/CFDNetworkInput.tmp");

		loadClusterInput >> countInput;
		for(i = 1;i <= countInput;i++)
		{
			loadClusterInput >> indexCell >> massIn >> numOriginalSpecies;
			iCluster = giveClusterIndexFromCellIndex[indexCell];
			massFeed[iCluster] += massIn;

			for(j = 1;j <= numOriginalSpecies;j++)
			{
				int indexSpecies;
				loadClusterInput >> indexSpecies;
				loadClusterInput >> massIn;
				feed[iCluster][indexSpecies] += massIn;
			}
		}
		loadClusterInput.End();
	}

	if (OnlyClustering == false)
	{
		// Cleaning Ms (sparse)
		// The elements very small are erased to increase the sparsity pattern of this
		// matrix without affecting the accuracy
		Ms.CleanMatrix(1.e-100);

		// Costruzione della matrice Ls contenente tutte le informazioni per la risoluzione
		// del sistema non lineare; questa matrice contiene anche gli elementi sulla diagonale
		// principale
		Ls = Ms;

		// Eliminazione della diagonale principale dalla matrice Ms, in modo da aumentarne
		// ulteriormente la struttura di sparsita'
		for(i = 1;i <= numCSTRReactors;i++)
			Ms.RemoveElement(i,i);

		// Matrice per Gauss Jordan
		// La matrice Ld e' identica alla matrice Ls, con la differenza che non ha pero'
		// elementi sulla diagonale principale
		ReplaceBzzMatrixSparseWithMinusBzzMatrixSparseLocked(&Ms,&Ld);
	}

	// Memory allocation
	ChangeDimensions(numComponents, &cReactor);
	ChangeDimensions(numComponents, &omegaReactor);
	ChangeDimensions(numComponents, &xReactor);
	ChangeDimensions(numComponents, &RReactor);

	// Inizializzazione della matrice delle frazioni massive
	if(iaia == 0)
	{
		cout << " - 05 - Reading data from FirstGuess.bzz file..." << endl;

		int numC;

		BzzLoad load(firstGuessFile);
		load >> numC >> iSpec;
		load.End();

		if(numC != numComponents)
			BzzError("numC != numComponents");
		ChangeDimensions(numCSTRReactors,numComponents,&massFractionsInReactors);

		for(i = 1;i <= numCSTRReactors;i++)
			for(k = 1;k <= iSpec.Size();k++)
				massFractionsInReactors[i][iSpec[k]] = omegaReduced_Cluster[i][k];

		ChangeDimensions(numCSTRReactors,&cstrSequence);
		for(i = 1;i <= numCSTRReactors;i++)
			cstrSequence[i] = i;
	}

	else
	{
		cout << " - 05 - Recovering data from mass.tmp file..." << endl;

		int numC,jC;
		BzzVectorInt aui;

		BzzLoad load('*', "Temp/mass.tmp");
		BzzMatrix aux;
		load >> numC >> aui >> aux;
		load.End();

		if(aui.Size() != originalNumCSTRReactors)
			BzzError("Data in mass.tmp file are not congruent!");
		ChangeDimensions(numCSTRReactors,numComponents,&massFractionsInReactors);

		for(i = 1;i <= originalNumCSTRReactors;i++)
		{
			j = giveClusterIndexFromCellIndex[i];
			jC = aui[i];

			for(k = 1;k <= numComponents;k++)
				massFractionsInReactors[j][k] = aux[jC][k];
		}

		ChangeDimensions(numCSTRReactors,&cstrSequence);
		for(i = 1;i <= numCSTRReactors;i++)
			cstrSequence[i] = i;
	}


	// Write cluster network file
	{
		cout << " - 06 - Writing cluster network file..." << endl;
		ofstream fClusterNetwork;
		openOutputFileAndControl(fClusterNetwork, "Temp/ClusterNetwork.tmp");
		fClusterNetwork.setf(ios::scientific);
		for (i=1; i<=originalNumCSTRReactors; i++)
			fClusterNetwork << giveClusterIndexFromCellIndex[i] << endl;
		fClusterNetwork.close();
	}

	// Write temperature network file
	{
		cout << " - 07 - Writing temperature network file..." << endl;
		ofstream fTemperature;
		openOutputFileAndControl(fTemperature, "Output/maps/Temperature.map");
		fTemperature.setf(ios::scientific);
		for (i=1; i<=originalNumCSTRReactors; i++)
			fTemperature << temperature[giveClusterIndexFromCellIndex[i]] << endl;
		fTemperature.close();
	}

	// Write sum cells network file
	{
		cout << " - 08 - Writing sum cells network file..." << endl;
		ofstream fSumCells;
		openOutputFileAndControl(fSumCells, "Output/maps/SumCells.map");
		fSumCells.setf(ios::scientific);
		BzzVectorInt sum(numCSTRReactors);
		for (i=1; i<=originalNumCSTRReactors; i++)
			sum[giveClusterIndexFromCellIndex[i]]++;
		for (i=1; i<=originalNumCSTRReactors; i++)
			fSumCells << double(sum[giveClusterIndexFromCellIndex[i]]) << endl;
		fSumCells.close();
	}

	// Write random network file
	{
		cout << " - 09 - Writing random network file..." << endl;
		ofstream fRandom;
		openOutputFileAndControl(fRandom, "Output/maps/Random.map");
		fRandom.setf(ios::scientific);
		BzzVector vRandom(giveClusterIndexFromCellIndex.Size());
		srand( (unsigned)time(NULL) );
		for(i=1;i<=giveClusterIndexFromCellIndex.Size();i++)
			vRandom[i] = rand() % 100 + 1;
		for (i=1; i<=originalNumCSTRReactors; i++)
			fRandom << vRandom[giveClusterIndexFromCellIndex[i]] << endl;
		fRandom.close();
	}

	// Write file clustering
	if (OnlyClustering == true)
	{
		{
			cout << "I am writing Clustering.FirstGuess.bzz" << endl;

			ofstream fClustering;
			fClustering.open("Clustering.FirstGuess.bzz", ios::out);
			fClustering.setf(ios::scientific);
		
			fClustering << numComponents << endl;
			fClustering << iSpec.Size() << endl;
			for(k = 1;k <= iSpec.Size();k++)
				fClustering << iSpec[k] << endl;
			fClustering << numCSTRReactors << " " << iSpec.Size() << endl;
			for(int i=1;i<=numCSTRReactors;i++)
			{
				for(k = 1;k <= iSpec.Size();k++)
					fClustering << omegaReduced_Cluster[i][k] << endl;	
			}
			fClustering.close();
		}

		{
			cout << "I am writing Clustering.CFDNetwork.bzz" << endl;
			ofstream fClustering;
			fClustering.open("Clustering.CFDNetwork.bzz", ios::out);
			fClustering.setf(ios::scientific);
			fClustering << numCSTRReactors << endl;
			for(int i=1;i<=numCSTRReactors;i++)
			{
				fClustering << i << "\t" << temperature[i] << "\t" << pressure[i] << "\t" << volume[i] << "\t" << qanm[i] << endl;

				// External inputs
				double sum = 0.;
				for(int j=1;j<=iSpec.Size();j++)
					sum+=feed[i][iSpec[j]];
				
				if (massInput[i] > 0.)
				{	
					if (fabs(massInput[i]-sum)/sum > 1.e-4)
					{
						cout << massInput[i] << " " << sum << endl;
						ErrorMessage("Clsuter  Error in sum of feeds!");
					}

					fClustering << 1 << endl;
					fClustering << sum << endl;
					fClustering << iSpec.Size() << endl;
					for(int j=1;j<=iSpec.Size();j++)
						fClustering << iSpec[j] << "\t" << feed[i][iSpec[j]]/sum << endl;
					
				}
				else
					fClustering << 0 << endl;

			
				
				// Communication
				
				if (FromClusterToCluster_Index[i].Size() > 0)
				{
					BzzVectorInt indices = FromClusterToCluster_Index[i];
					BzzVector masses     = FromClusterToCluster_MassFlowRate[i];
					BzzVector diffusions = FromClusterToCluster_DiffusionFlowRate[i];
					BzzVectorInt sort(indices.Size());
					Sort(&indices, &sort);
					Reorder(&masses, sort);
					Reorder(&diffusions, sort);
					
					BzzVectorInt count(1,1);
					for(int j=1;j<=sort.Size()-1;j++)
						if (indices[j+1] != indices[j])	count.Append(j+1);
					
					
					fClustering << count.Size() << endl;
					
					count.Append(sort.Size()+1);
					for(int j=1;j<=count.Size()-1;j++)
					{	
						double sum_masses = 0.;
						double sum_diffusions = 0.;
						for (int kk=count[j];kk<count[j+1];kk++)
						{
							sum_masses += masses[kk];
							sum_diffusions += diffusions[kk];
						}
						fClustering << indices[count[j]] << "\t" << sum_masses << "\t" << sum_diffusions << endl;
					}

				//	fClustering << FromClusterToCluster_Index[i].Size() << endl;
				//	for(int j=1;j<=FromClusterToCluster_Index[i].Size();j++)
				//		fClustering << FromClusterToCluster_Index[i][j] << "\t" << FromClusterToCluster_MassFlowRate[i][j] << "\t" << FromClusterToCluster_DiffusionFlowRate[i][j] << endl;
				}
			}
			
			fClustering.close();	
		}

		cout << "Clustering successfull!" << endl;
		exit(0);
	}

	// Open fResidualFile
	openOutputFileAndControl(fHistory,		"History.out");
	fHistory.setf(ios::scientific);
	openOutputFileAndControl(fResiduals_1,	"Output/residuals/Residuals_1.out");
	fResiduals_1.setf(ios::scientific);
	fResiduals_1 << "#Iter.(total)[1]  " << "\t";
	fResiduals_1 << "#Iter.(partial)[2]" << "\t";
	fResiduals_1 << "Res.(mean)[3]     " << "\t";
	fResiduals_1 << "Res.(max)[4]      " << "\t";
	fResiduals_1 << "CSTRNoOK[5]       " << "\t";
	fResiduals_1 << "EqsNoOK[6]        " << "\t";
	fResiduals_1 << "Time(global)[7]   " << "\t";
	fResiduals_1 << "Time(local)[8]    " << "\t";
	fResiduals_1 << endl << endl;
	
	openOutputFileAndControl(fResiduals_2,	"Output/residuals/Residuals_2.out");
	fResiduals_2.setf(ios::scientific);
	fResiduals_2 << "#Iter.(total)[1]  " << "\t";
	fResiduals_2 << "F0(Before)[2]     " << "\t";
	fResiduals_2 << "F0Max(Before)[3]  " << "\t";
	fResiduals_2 << "F1(Before)[4]     " << "\t";
	fResiduals_2 << "F1Max(Before)[5]  " << "\t";
	fResiduals_2 << endl << endl;

	openOutputFileAndControl(fResiduals_3,	"Output/residuals/Residuals_3.out");
	fResiduals_3.setf(ios::scientific);
	fResiduals_3 << "#Iter.(total)[1]  " << "\t";
	fResiduals_3 << "#Iter.(partial)[2]" << "\t";
	fResiduals_3 << "F0(Before)[3]     " << "\t";
	fResiduals_3 << "F0Max(Before)[4]  " << "\t";
	fResiduals_3 << "F1(Before)[5]     " << "\t";
	fResiduals_3 << "F1Max(Before)[6]  " << "\t";
	fResiduals_3 << endl << endl;

	if (OnlyODE == true)
	{
	openOutputFileAndControl(fResiduals_4,	"Output/residuals/Residuals_4.out");
	fResiduals_4.setf(ios::scientific);
	fResiduals_4 << "Iteration   " << "\t";
	fResiduals_4 << "time        " << "\t";
	fResiduals_4 << "Fmean[-]    " << "\t";
	fResiduals_4 << "Fmax[5]     " << "\t";
	fResiduals_4 << endl << endl;
	}

	countThirdTotal  = 0;
	countSecondTotal = 0;
	countFourthTotal = 0;
	countTuttoTotal	 = 0;
}

void OpenSMOKE_CSTRNetwork::operator()(const string fileNetwork, const string first, const string fileTolerances,
								  const double cicloCluster, const int cicloDiffusion, const int iaia, 
								  const int relaxation, const int _iAnalyticalJacobian, const SymbolicKinetics _analyticalJacobian, const int _iKindOfCorrection)
{
	int countInput;
	int numCluster;
	BzzVectorInt cstrClusterSize;
	BzzMatrix omegaReduced_Cluster;
	BzzMatrix omegaReduced_CFD;

	firstGuessFile		= first;
	networkFile			= fileNetwork;
	iAnalyticalJacobian = _iAnalyticalJacobian;
	analyticalJacobian	= _analyticalJacobian;
	iKindOfCorrection	= _iKindOfCorrection;

	read_first_guess(first);
	reading_the_network_file(fileNetwork, countInput);
	read_tolerances(fileTolerances, numCSTRReactors, cicloCluster);
	build_cluster(cicloCluster, numCSTRReactors, countInput, numCluster, cstrClusterSize);
	calculate_initial_massfractions_in_clusters(numCluster, omegaReduced_Cluster, omegaReduced_CFD);
	calculate_clustering(fileNetwork, cicloDiffusion, numCluster, cstrClusterSize, omegaReduced_Cluster, omegaReduced_CFD);
	complete_clustering(countInput, cicloDiffusion, relaxation, iaia, omegaReduced_Cluster);
	start_solving_network(iaia);
}

void OpenSMOKE_CSTRNetwork::start_solving_network(const int from_cfd_results)
{
	cout << endl;
	cout << " - 10 - Temperature dependences creation and memorization..."		<< endl;

	// 0=None 1=Sin expansion 2= Delta Dirac  (3=Beta PDF 4=ClippedGaussian)
	MemoTemperatureFunctions(iKindOfCorrection);

	cout << "-----------------------------------------------------------------" << endl;
	cout << endl;

	cout << "-----------------------------------------------------------------" << endl;
	cout << " START COMPUTATIONS" << endl;
	cout << "-----------------------------------------------------------------" << endl;
	cout << endl;

	// Initial Solution
	massFractionsInReactorsSolution = massFractionsInReactors;

	// Print Temperature File
	SaveTemperature('*', "Temp/temperature.tmp");

	// Connection Matrix
	//ConnectionMatrix();

	// Enthalpy analysis
	if (iEnergyAnalysis == true)
	{
	//	if (originalNumCSTRReactors == numCSTRReactors && from_cfd_results == 0)
	//	{
	//		cout << "Energy Analysis on the original network ... " << endl;
	//		string fileNameEnthalpy				= folderCase + "Temp/Enthalpy/OriginalEnthalpy.map";
	//		string fileNameDeltaEnthalpy		= folderCase + "Temp/Enthalpy/OriginalDeltaEnthalpy.map";
	//		string fileNameDifferenceEnthalpy	= folderCase + "Temp/Enthalpy/OriginalDifferenceEnthalpy.map";
	//		string fileNameInletEnthalpy		= folderCase + "Temp/Enthalpy/OriginalInletEnthalpy.map";
	//		string fileNameErrorEnthalpy		= folderCase + "Temp/Enthalpy/OriginalErrorEnthalpy.map";
	//		EnthalpyAnalysis(	fileNameEnthalpy, fileNameDeltaEnthalpy, fileNameDifferenceEnthalpy, 
	//							fileNameInletEnthalpy, fileNameErrorEnthalpy);
	//		cout << "DONE!" << endl;
	//	}
	//	else
		{
			cout << "Energy Analysis on the reactor network ... " << endl;
			string fileNameEnthalpy				= "Output/enthalpy/Enthalpy.map";
			string fileNameInletEnthalpy		= "Output/enthalpy/InletEnthalpy.map";
			string fileNameDeltaEnthalpy		= "Output/enthalpy/DeltaEnthalpy.map";
			string fileNameDifferenceEnthalpy	= "Output/enthalpy/DifferenceEnthalpy.map";
			string fileNameErrorEnthalpy		= "Output/enthalpy/ErrorEnthalpy.map";
			EnthalpyAnalysis(	fileNameEnthalpy, fileNameDeltaEnthalpy, fileNameDifferenceEnthalpy, 
								fileNameInletEnthalpy, fileNameErrorEnthalpy);
			cout << "DONE!" << endl;
		}
	}

	// Global Algorithm
	if (OnlyODE == false)	globalAlgorithm();
	if (OnlyODE == true)	odeOnlyAlgorithm();
}

void OpenSMOKE_CSTRNetwork::globalAlgorithm()
{
	double start;
	double end;

	const int MaxNumberOfIterations = 4;

	// Convergence Criteria
	maxNumberOfIterations_Global = Max(numCSTRReactors / 3,30);

	int okok;
	int countCicle = 1;
	for(int iter = 1; iter<=MaxNumberOfIterations;iter++)
	{
		fHistory << "Cycle:  " << iter << "." << countCicle       << endl;
		fHistory << "TolRel: " << tolRel << "\t" << "TolAbs: " << tolAbs << endl;
		fHistory << "--------------------------------------------------" << endl;

		// 1. Vengono risolti singolarmente tutti i reattori CSTR
		//    L'ordine di risoluzione non e' importante e viene praticamente mantenuto
		//    quello della simulazione originale CFD
		//    0 = se la convergenza non e' stata raggiunta in maniera appropriata
		//	  1 = se la convergenza e' stata ottenuta in maniera corretta
		fHistory << "First: " << endl;
		start = BzzGetCpuTime();
		
			okok = GetFirst();
		
		end = BzzGetCpuTime();
		fHistory << " * cpu time (min): " << (end-start)/60. << endl;
		fHistory << " * exit status:    " << okok << endl;

		// 2. Stampa le informazioni a video
		OutputPrint(countCicle++);

		// Se le iterazioni sui singoli reattori sono andate a buon fine e si sono fatte
		// almeno tre iterazioni globali del metodo, allora la rete di reattori puo'
		// considerarsi risolta ed e' possibile uscire dal programma
		if(okok == 1 && iter == 3)
			return;

		// 3. Viene risolta l'intera rete di reattori
		fHistory << "Third: " << endl;
		start = BzzGetCpuTime();
			
			okok = GetThird();
			
		end = BzzGetCpuTime();
		fHistory << " * cpu time (min): " << (end-start)/60. << endl;
		fHistory << " * exit status:    " << okok << endl;

		// 4. Stampa le informazioni a video
		OutputPrint(countCicle++);

		// Se il sistema globale ha raggiunto la convergenza e' possibile uscire
		if(okok == 1)
			return;

		// 5. Viene chiamata solo a partire dalla seconda iterazione in poi
		if(iter >= 2)
		{
			fHistory << "Second: " << endl;
			start = BzzGetCpuTime();

				GetSecond();

			end = BzzGetCpuTime();
			fHistory << " * cpu time (min): " << (end-start)/60. << endl;

			fHistory << "Third: " << endl;
			start = BzzGetCpuTime();

				okok = GetThird();

			end = BzzGetCpuTime();
			fHistory << " * cpu time (min): " << (end-start)/60. << endl;
			fHistory << " * exit status:    " << okok << endl;
			OutputPrint(countCicle++);
			if(okok == 1)
				return;
		}

		fHistory << endl;

		// New strict tolerances
		tolRel *= .1;
		tolAbs *= .1;
		if(tolRel < MaxTolRel)	tolRel = MaxTolRel;
		if(tolAbs < MaxTolAbs)	tolAbs = MaxTolAbs;
	}
}

void OpenSMOKE_CSTRNetwork::odeOnlyAlgorithm()
{
	double start;
	double end;

	const int MaxNumberOfIterations = 6;

	// Convergence Criteria
	maxNumberOfIterations_Global = Max(numCSTRReactors / 3,30);

	int okok;
	int countCicle = 1;
	double tEnd0 = 10.;
	for(int iter = 1; iter<=MaxNumberOfIterations;iter++)
	{
		fHistory << "Cycle:  " << iter << "." << countCicle       << endl;
		fHistory << "TolRel: " << tolRel << "\t" << "TolAbs: " << tolAbs << endl;
		fHistory << "--------------------------------------------------" << endl;

			{
				fHistory << "Fourth: " << endl;
			
				start = BzzGetCpuTime();

				if (iter<=3)	GetFourth(tEnd0);
				else			GetFourth(tEnd0*pow(3., double(iter-3.)));

				end = BzzGetCpuTime();
			
				fHistory << " * cpu time (min): " << (end-start)/60. << endl;

				OutputPrint(countCicle++);
			}

			{
				fHistory << "First: " << endl;
			
				start = BzzGetCpuTime();
			
					okok = GetFirst();
			
				end = BzzGetCpuTime();
				
				fHistory << " * cpu time (min): " << (end-start)/60. << endl;
				fHistory << " * exit status:    " << okok << endl;

				OutputPrint(countCicle++);
			}

			if (iter >= 3)
			{
				fHistory << "Third: " << endl;
				
				start = BzzGetCpuTime();

					okok = GetThird();

				end = BzzGetCpuTime();
				fHistory << " * cpu time (min): " << (end-start)/60. << endl;

				OutputPrint(countCicle++);
				
				if(okok==1 && iter>=4)	return;
			}


		fHistory << endl;

		// New strict tolerances
	//	tolRel *= .1;
	//	tolAbs *= .1;
	//	if(tolRel < MaxTolRel)	tolRel = MaxTolRel;
	//	if(tolAbs < MaxTolAbs)	tolAbs = MaxTolAbs;
	}
}

void OpenSMOKE_CSTRNetwork::GetResiduals(BzzMatrix &massFractionsInReactors, BzzMatrix &residuals)
{
	Product(Ls,massFractionsInReactors,&residuals);
	GetReactionsRateInAllReactorsFromMassFractions(massFractionsInReactors, dummyMatrix);
	ElementByElementProduct(volumeMolecularWeightT,dummyMatrix,&dummyMatrix);
	Difference(dummyMatrix,&residuals);
	Sum(feed,&residuals);
}

void OpenSMOKE_CSTRNetwork::GetAllResiduals(BzzVector &omega,BzzVector &res)
{
	Product(Ls,omega,&res);
	GetReactionsRateInAllReactorsFromMassFractions(omega,dummyVector);
	ElementByElementProduct(volumeMolecularWeightT_Vector,dummyVector,&dummyVector);
	Difference(dummyVector,&res);
	Sum(feed_Vector,&res);
}

void OpenSMOKE_CSTRNetwork::GetReactionRatesInSingleReactor(int iReactor,
							BzzVector &omega, BzzVector &molefractions,
							BzzVector &concentrations, BzzVector &ReactionRates)
{
	wM = Reactions->GetMWFromMassFractions(omega);
	Reactions->GetMoleFractionsFromMassFractionsAndMW(molefractions, omega, wM);
	concentrations = cTot*molefractions;

	// Reaction Rates
	#if SYMBOLIC_KINETICS==1
		if (iAnalyticalJacobian == 0)
			Reactions->ComputeFromConcentrations_map(iReactor, concentrations, &ReactionRates);
		else
			reactor[iReactor]->giveReactionRates(cTot, concentrations, ReactionRates);
	#else
			Reactions->ComputeFromConcentrations_map(iReactor, concentrations, &ReactionRates);
	#endif
}

void OpenSMOKE_CSTRNetwork::GetResiduals(BzzVector &m,BzzVector &residuals)
{
	GetReactionRatesInSingleReactor(kReactor, m, xReactor, cReactor, RReactor);

	volumeMolecularWeightT.GetRow(kReactor,&cReactor);
	ElementByElementProduct(RReactor,cReactor,&xReactor);
	Product(sM[kReactor],m,&cReactor);
	Difference(&cReactor,feedInkRreactor);
	Difference(cReactor,xReactor,&residuals);
}

void OpenSMOKE_CSTRNetwork::GetJacobian(BzzVector &omega,BzzMatrix &JJ)
{
	GetReactionRatesInSingleReactor(kReactor, omega, xReactor, cReactor, RReactor);
	GetJacobianForSingleReactor(kReactor, cReactor, RReactor, JJ);
}

void OpenSMOKE_CSTRNetwork::GetJacobianForSingleReactor(int iReactor, BzzVector &cReactor, BzzVector &R, BzzMatrix &dRC)
{
	int i, j;
	double wc;

	// Calculation of the Jacobian
	if (iAnalyticalJacobian == 0)
	{
		Reactions->kinetics.mc		= cReactor;
		Reactions->kinetics.mR		= R;
		Reactions->kinetics.mr		= Reactions->r;
		Reactions->kinetics.mrDirT	= Reactions->coeffFallOff;
		Reactions->kinetics.mrDirC	= Reactions->rDirC;
		Reactions->kinetics.mrInvC	= Reactions->rInvC;

		Reactions->kinetics.GetDerivativesC(T, cTot, &dRC, cReactor, R);
	}

	else
	{
		#if SYMBOLIC_KINETICS==1
			ChangeDimensions(numComponents, numComponents, &dRC);
			reactor[iReactor]->giveJacobian(cReactor, dRC);
		#endif
	}

	// 1. Contributi dai termini di reazione
	wc = wM * cTot;
	for(i = 1;i <= numComponents;i++)
		xReactor[i] = wc / Reactions->M[i];
	dRC.ColumnsElementsProduct(xReactor);
	for(i = 1;i <= numComponents;i++)
		for(j = 1;j <= numComponents;j++)
			dRC[i][j] *= volumeMolecularWeightT[iReactor][i];

	// 2. Contributo derivante dal termine convettivo
	for(i = 1;i <= numComponents;i++)
		dRC[i][i] -= sM[iReactor];
}


int OpenSMOKE_CSTRNetwork::GetFirst(void)
{
	// Indici
	bool iPrintReactor;
	int i,j;
	BzzLoad load;

	// Contatori per il numero di reattori/specie non soddisfatti
	int minNoOK;					// numero minimo di reattori non OK
	int minNoOkTot;					// numero minimo delle specie (totale) non OK

	// Residui per tutte le variabili
	double resMean;					// Residuo medio
	double resMax;					// Residuo massimo
	double resMeanOpt;				// Residuo medio ottimale
	double maxResOpt;				// Residuo massimo ottimale

	// Residui corrispondenti al metodo di Newton vero e proprio
	double F0;						// Residuo medio corrispondente al punto di first guess
	double F1;						// Residuo medio corrispondente al punto aggiornato
	double F1Max;					// Residuo massimo
	double FIntegral;				// Residuo medio relativo al termine di accumulo

	// Informazioni sul peggior reattore integrato temporalmente
	int kMaxIntegration;			// Indice peggior reattore in seguito integrazione temporale
	double F1MaxIntegration;		// Massimo residuo di accumulo
	double TMaxIntegration;			// Temperatura peggior reattore

	// Fattori di normalizzazione
	double uvariables;
	double ucomponents;

	// Contatori
	int countTutto;					// Contatore delle iterazioni globali
	int countReactorsOK_FirstGuess;	// Contatore di reattori OK al first guess
	int countIntegrations;			// Contatore locale di ODE systems risolti
	int countNewtonIterations;		// Contatore del numero di iterazioni Newton sul singolo CSTR
	int countNoMaxRes;				// Contatore delle iterazioni senza aggiornamento o miglioramento della soluzione
	int countNoIntegration;			// Contatore iterazioni globali senza nessun reattore risolto con il transitorio

	// Varie
	int memoCount;					// Costante
	double maxOdeSum;				// Criterio di stop per il sistema ODE sul singolo CSTR
	int numF16;						// Numero di reattori che non soddifano
									// un valore sufficientemente basso del residuo in seguito a Newton
	double h;
	double hMin;


	// Inizializzazione del sistema ODE per il singolo CSTR
	BzzOdeStiffObject ode;
	OpenSMOKE_CSTRNetwork_MyOdeSystemObjectOneCSTR cstrMono(this);

//	cstrMono(this);

	// Allocazione di memoria
	BzzVector x0,x1,d1;
	BzzVector Ras;
	BzzFactorizedGauss G;
	BzzMatrix dRC(numComponents,numComponents);
	BzzMatrix Res(numCSTRReactors,numComponents);
	BzzVector f0(numComponents);
	BzzVector	f1(numComponents);
	BzzVector R(numComponents);
	BzzVector Rvw(numComponents);
	BzzVectorInt vectorReactorsNotOK(numCSTRReactors);	// flag per individuare reattori non OK

	// Inizializzazione variabili
	memoCount			= 30;
	minNoOK				= numCSTRReactors + 1;
	minNoOkTot			= (numCSTRReactors*numComponents) + 1;
	maxOdeSum			= .1 * F1Stop * double(numComponents);
	countNoIntegration	= 0;
	countNoMaxRes		= 0;
	countTutto			= 0;
	uvariables			= 1. / double(numComponents * numCSTRReactors);
	ucomponents			= 1. / double(numComponents);
	
	const double	odeTolAbs		= 1.e-15;	// Default 1e-15
	const int		odeMaxJacobian	= 3;		// Default 3
	const double	MaxTau			= 100.;		// Default 100.

	// Start times
	double start = BzzGetCpuTime();
	double memoStart = start;

	// START COMPUTATIONS
	GetResiduals(massFractionsInReactorsSolution,Res);
	resMeanOpt = uvariables * Res.GetSumAbsElements();
	maxResOpt = Res.MaxAbs();

	cout << "STARTING GET FIRST APPROACH!" << endl;
	cout << "  ** Optimal Mean Residual: " << resMeanOpt << endl;
	cout << "  ** Optimal Max  Residual: " << maxResOpt << endl;
	

	// ---------------------------------------------------------------------------------------//
	//							GLOBAL ITERATION
	// ---------------------------------------------------------------------------------------//
	int iSequence=1;

	ripetiTutto:

	// Reset Counts
	countTutto++;
	countTuttoTotal++;
	countIntegrations			= 0;
	countReactorsOK_FirstGuess	= 0;
	F1Max						= 0.;
	numF16						= 0;
	F1MaxIntegration			= 0.;
	kMaxIntegration				= 0;
	maxOdeSum = .1 * F1Stop * double(numComponents);

	// Loading MemoTemperatureFunctions.tmp
	if(memoTemperature == 1)
		load('*', "Temp/MemoTemperatureFunctions.tmp");

	BzzVectorInt sequence(numCSTRReactors);
	if (iSequence==1)
		for(kReactor = 1;kReactor <= numCSTRReactors;kReactor++)
			sequence[kReactor] = kReactor;
	else
		for(kReactor = 1;kReactor <= numCSTRReactors;kReactor++)
			sequence[kReactor] = numCSTRReactors+1-kReactor;

	// Cycle on every reactor
	//for(kReactor = 1;kReactor <= numCSTRReactors;kReactor++)
	//for(kReactor = numCSTRReactors;kReactor >= 1;kReactor--)
	for(int index=1;index<=numCSTRReactors;index++)
	{
		iSequence *= -1;
		kReactor = sequence[index];

		if (kReactor%(int(0.05*numCSTRReactors)+1) == 0)	iPrintReactor = true;
		else											iPrintReactor = false;
		if(iPrintReactor==true)
		{
			cout << " -- R#: "		<< kReactor	<< "/"		<< numCSTRReactors << "\t";
			cout << " -- Tour# "	<< countTutto	<< "-"	<< countTuttoTotal << "\t";
			cout << endl;
		}

		countNewtonIterations = 0;

		// Recovering Physical Information
		T = temperature[kReactor];
		P = pressure[kReactor];
		logT = logTm[kReactor];
		loguRT = loguRTm[kReactor];
		cTot = cTotm[kReactor];

		// Recovering thermodynamic properties
		if(memoTemperature == 1)
		{
			load >> Reactions->uKeq;
			load >> Reactions->k1 >> Reactions->k2;
			load >> Reactions->logFcent;
		}

		// Recovering Actual solution in reactor
		massFractionsInReactorsSolution.GetRow(kReactor,&x0);

		// Checking on Sum mass fractions
		double su = x0.GetSumAbsElements();
		if(su > 1.0001 || su < .9999)
		{
			cout << endl << "WARNING: Reactor " << kReactor << " - Sum: " << su << endl;
			if(su == 0.)
				cout << "WARNING: Wrong composition in reactor " << kReactor;

			// Correction
			su = 1. / su;
			for(i = 1;i <= numComponents;i++)
				x0[i] *= su;
		}

		// Residual calculations
		omegaReactor = x0;
		feed.GetRow(kReactor,&feedInkRreactor);
		ProductForSelectedRow(Ls,massFractionsInReactorsSolution,kReactor,&cReactor);

		// Questa riga viene aggiunta perche' serve per i successivi usi di feedInkReactor
		for(i = 1;i <= numComponents;i++)
			cReactor[i] -= massFractionsInReactorsSolution[kReactor][i] * sM[kReactor];
		Difference(&feedInkRreactor,cReactor);

		// Reaction Rates Calculation
		GetReactionRatesInSingleReactor(kReactor, omegaReactor, xReactor, cReactor, R);
		for(i = 1;i <= numComponents;i++)
			Rvw[i] = R[i] * volumeMolecularWeightT[kReactor][i];

		// Residual calculations
		// Questa riga viene inserita per compensare quella scritta sopra, utilizzata
		// per la definizione di feedInkReactor
		Product(sM[kReactor],x0,&Ras);
		Difference(&Ras,feedInkRreactor);
		Difference(Ras,Rvw,&f0);

		// Mean Residual for Single Reactor
		F0 = f0.GetSumAbsElements() * ucomponents;

		// If the mean Residual is low enough the cycle goes to the next reactor
		if(F0 < 1.e-6 && F0 < resMeanOpt && vectorReactorsNotOK[kReactor] == 0)
		{
			countReactorsOK_FirstGuess++;
			continue;
		}

		// Calcolo dello jacobiano
		GetJacobianForSingleReactor(kReactor, cReactor, R, dRC);
		G = dRC;

		// Newton-Raphson method
		// The following equation is solved to get the new solution
		// J * (xNew-xOld) = -F
		ripetiNewton:

			countNewtonIterations++;

			// In the expression above J is the jacobian and F is the vector of residuals
			Solve(&G,f0,&d1);

			// Viene ottenuto il nuovo punto predetto dal metodo di Newton; tuttavia
			// successivamente questo nuovo punto potrebbe essere addolcito nel caso in
			// vengano incontrate delle opportune condizioni
			Sum(x0,d1,&x1);


			hMin = 1.;
			for(i = 1;i <= numComponents;i++)
			{
				// In questo caso le frazioni massive potrebbero diventare negative
				if(d1[i] < 0.)  //-1.e-5 && x0[i] > 1.e-12)
				{
					// A. Se il punto di partenza era gia' zero ci sono due possibilita'
					if(x0[i] == 0.)
					{
						// A1. Se la variazione non e' piccolissima allora e' meglio abbandonare
						// da subito il metodo di Newton e quindi settare hMin=0, il che
						// obbliga successivamente a eseguire l'integrazione temporale
						if(d1[i] < -1.e-15)
						{
							hMin = 0.;
							break;
						}
						// A2. Se invece la variazione e' piccolissima e' sufficiente porre tale
						// variazione pari a zero econtinuare
						else
							d1[i] = 0.;
					}
					else
					{
						// E' una misura relativa di quanto ho modificato il punto ed e' sempre
						// positiva; tanto maggiore e', tanto piu' grande e' lo spostamento del punto
						h = -x0[i] / d1[i];
						if(h < hMin)
							hMin = h;
					}
				}
			}

			// Se vengono incontrate le condizioni opportune viene risolto il transitorio sul CSTR
			// in modo da avere una maggiore robustezza per l'ottenimento della soluzione
			// Il transitorio viene risolto o quando e' stato incontrato il punto A1, oppure
			// se la soluzione non e' stata spostata di molto, ma solo la prima volta
			double tau = .01;
			if(hMin <= tau && countNewtonIterations == 1)
			{
				ode.SetInitialConditions(x0, 0., &cstrMono, &dRC);
				ode.SetAnalyticalJacobian();									// TODO
				ode.SetTollAbs(odeTolAbs);
				ode.StopIntegrationBeforeRecalcuatingJacobian(odeMaxJacobian);	// TODO
				ode.StopIntegrationWhenSumAbsY1IsLessThan(maxOdeSum);
				BzzVector yMin(numComponents);
				ode.SetMinimumConstraints(&yMin);
				omegaReactor = ode(MaxTau);

				if(ode.GetOdeCalculationState() < 0)
				{
					omegaReactor.BzzPrint("omegaReactor");
					BzzWarning("WARNING: The solution of the ODE system on reactor %d met some problems: %d",
								kReactor, ode.GetOdeCalculationState());
				}
				else
					massFractionsInReactorsSolution.SetRow(kReactor,omegaReactor);

				// Residui (termini di accumulo)
				f1 = ode.GetY1InMeshPoint();

				// Viene individuato il reattore che presenta i maggiori residui temporali
				FIntegral = f1.GetSumAbsElements() * ucomponents;
				if(FIntegral > F1MaxIntegration)
				{
					TMaxIntegration = T;
					kMaxIntegration = kReactor;
					F1MaxIntegration = FIntegral;
				}

				// Nel caso in cui stato condotto il transitorio si passa comunque al
				// successivo reattore e si salta completamente la parte seguente
				countIntegrations++;
				continue;
			}

			// Questa parte del codice viene presa in considerazione solo se non e' stato
			// effettuato il transitorio sul CSTR, segno che quindi il metodo di Newton
			// aveva dato buoni risultati.
			// Non viene preso in considerazione l'intero spostamento predetto dal metodo di Newton
			// ma soltanto una frazione che tiene conto di 0.90*hMin e quindi evita che possano
			// esserci valori negtivi delle frazioni massive
			// Viene quindi ricalcolato il nuovo punto x1
			double sumOmega;
			if(hMin != 1.)
			{
				hMin *= .9;
				Product(hMin,&d1);
				Sum(x0,d1,&x1);
			}

			// Normalizzazione delle frazioni massive del nuovo punto predetto dal metodo di Newton
			sumOmega = x1.GetSumElements();
			if(sumOmega < .99999999 || sumOmega > 1.000000001)
			{
				if(sumOmega == 0.)	ErrorMessage("The sum of mass fractions is ouside the ranges!");
				sumOmega = 1. / sumOmega;
				for(i = 1;i <= numComponents;i++)
					x1[i] *= sumOmega;
			}

			// Calcolo dei residui relativi alla nuova soluzione
			GetResiduals(x1,f1);

			// Calcolo del nuovo residuo medio corrispondente al nuovo punto
			F1 = f1.GetSumAbsElements() * ucomponents;
			double OldF0_Corrected = F0 * .9;

			// Se il nuovo residuo medio e' migliore del precedente, la nuova soluzione viene
			// considerata tale e pertanto il nuovo punto x0 diventa quello appena calcolato;
			// di conseguenza vengono modificate anche le altre informazioni
			if(F1 < F0)
			{
				// Aggiornamento del massimo residuo incontrato nel corso delle iterazioni
				if(F1 > F1Max)
					F1Max = F1;

				omegaReactor = x1;
				massFractionsInReactorsSolution.SetRow(kReactor,omegaReactor);
				F0 = F1;
				x0 = x1;
				f0 = f1;

				// Se il nuovo residuo non e' comunque piu' basso del valore ottimale richiesto,
				// viene aggiornato opportunamente il corrispondente contatore
				if(F1 > resMeanOpt)
					numF16++;

				// Se la soluzione non e' ancora soddisfacente il metodo di Newton viene
				// ripetuto; lo jacobiano non viene pero' aggiornato
				if(countNewtonIterations < MaxCountNewtonIterations && F1<OldF0_Corrected && F1>MAX_SUM_F)
					goto ripetiNewton;
			}

			// Se al contrario il nuovo residuo medio non e' migliore del precedente, non e'
			// detto necessariamente che il nuovo punto non possa essere accettato; bisogna
			// condurre degli ulteriori controlli

			// Se si tratta della prima iterazione, o comunque il vecchio residuo e' alto, viene
			// risolto il transitorio sul CSTR, anche se precedentemente era stato scartato
			else if(countNewtonIterations == 1 || F0 > resMeanOpt)
			{
				ode.SetInitialConditions(x0, 0., &cstrMono, &dRC);
				ode.SetAnalyticalJacobian();						// TODO
				ode.SetTollAbs(odeTolAbs); 
				ode.StopIntegrationBeforeRecalcuatingJacobian(odeMaxJacobian);	// TODO
				ode.StopIntegrationWhenSumAbsY1IsLessThan(maxOdeSum);
				BzzVector yMin(numComponents);
				ode.SetMinimumConstraints(&yMin);
				omegaReactor = ode(MaxTau);
				massFractionsInReactorsSolution.SetRow(kReactor,omegaReactor);

				if(ode.GetOdeCalculationState() < 0)
				{
					omegaReactor.BzzPrint("omegaReactor");
					BzzWarning("WARNING: The solution of the ODE system on cluster %d met some problems: %d",
								kReactor, ode.GetOdeCalculationState());
				}

				// Eventuale aggiornamento del peggior reattore
				f1 = ode.GetY1InMeshPoint();
				FIntegral = f1.GetSumAbsElements() * ucomponents;
				if(FIntegral > F1MaxIntegration)
				{
					TMaxIntegration = T;
					kMaxIntegration = kReactor;
					F1MaxIntegration = FIntegral;
				}

				// Una volta fatto il transitorio, si passa comunque in ogni
				// caso al reattore successivo
				countIntegrations++;
				continue;
			}
	}
	if(memoTemperature == 1)	load.End();
	
	// Time counters
	double global_time = BzzGetCpuTime() - start;
	double local_time  = BzzGetCpuTime() - memoStart;

	// Tutti i reattori sono stati risolti; bisogna adesso individuare
	// degli opportuni criteri di convergenza

	// Informazioni a video
	cout	<< endl
			<< "Total number of global calls: "			<< countTutto					<< endl
			<< "Number of CSTR already ok: "			<< countReactorsOK_FirstGuess	<< endl
			<< "Total number of ODE system solved: "	<< countIntegrations			<< endl
			<< "FStop Value (Assigned): "				<< F1Stop						<< endl
			<< "Maximum residual (Newton): "			<< F1Max						<< endl
			<< "Maximum Residual (Accumulation): "		<< F1MaxIntegration				<< endl
			<< "CPU Time (Global): "					<< BzzGetCpuTime() - start		<< endl
			<< "CPU Time (Local): "						<< BzzGetCpuTime() - memoStart	<< endl
			<< endl;

	// Se sono state fatte piu' di 30 chiamate e il residuo (di accumulo) e' ancora
	// alto cio' vuol dire che potrebbero esserci problemi: il tutto viene segnalato a video
	if(countTutto > 30 && F1MaxIntegration > 1.e-5)
	{
		BzzWarning("WARNING: The number of calls is higher than 30 and the integral residual is still high!");
		cout << "WARNING: The number of calls is higher than 30 and the integral residual is still high!" << endl;
		cout << "Maximum residual (Accumulation): " << F1MaxIntegration << endl;
	}

	// Stampa a video le informazioni relative all'eventuale reattore che ha dato i peggiori
	// risultati dal punto di vista dell'integrazione temporale
	if(kMaxIntegration != 0)
		cout	<< "  ** Reactor number: "			<< kMaxIntegration					<< endl
				<< "  ** Index in the Sequence: "	<< cstrSequence[kMaxIntegration]	<< endl
				<< "  ** Temperature: "				<< TMaxIntegration					<< endl
				<< endl;

	// Stampa a video dei reattori per i quali il metodo di Newton non ha dato un residuo
	// piu' basso di quello ottimale
	if (numF16 != 0)
		cout << "  ** Number of reactors whose max residual (Newton) is higher than optimal value: "
			 << numF16 << endl;

	// Viene aggiornato il residuo globale con le soluzioni appena ottenute
	// su tutti i reattori; vengono quindi aggiornati il massimo e il valore medio
	GetResiduals(massFractionsInReactorsSolution,Res);
	resMean = uvariables * Res.GetSumAbsElements();
	resMax  = Res.MaxAbs();
	cout << "  ** Mean Residual: "		<< resMean << endl;
	cout << "  ** Maximum Residual: "	<< resMax << endl;

	// Viene eventualmente aggiornato il criterio di convergenza
	F1Stop = Min(F1Stop,resMean);

	// Viene aggiornata la variabile che tiene conto del numero di iterazioni globali
	// che non hanno portato a nessun aggiornamento apprezzabile della soluzione
	countNoMaxRes++;
	if (countTutto < 100) countNoMaxRes = 0;

	// Vengono esaminati i residui sulle singole specie per ciascun reattore
	int nReactorsNotOK	= 0;	// e' il numero totale di reattori che non sono ok
	int nVariablesNotOK = 0;	// e' il numero totale di specie che non sono ok
	vectorReactorsNotOK = 0;
	for(i = 1;i <= numCSTRReactors;i++)
	{
		nReactorsNotOK++;
		for(j = 1;j <= numComponents;j++)
		{
			if(fabs(Res[i][j]) >  tolRel * massFractionsInReactorsSolution[i][j] + tolAbs)
			{
				nVariablesNotOK++;
				vectorReactorsNotOK[nReactorsNotOK] = 1;
			}
		}
	}
	nReactorsNotOK = vectorReactorsNotOK.GetSumElements();
	cout << "Total number of reactors which are not OK:          " << nReactorsNotOK << endl;
	cout << "Total number of balance equations which are not OK: " << nVariablesNotOK << endl;

	memoStart = BzzGetCpuTime();

	// Write on file
	fResiduals_1	<< countTuttoTotal													<< "\t" 
					<< countTutto														<< "\t"		
					<< resMean															<< "\t"		
					<< resMax															<< "\t" 
					<< double(nReactorsNotOK)/double(numCSTRReactors)					<< "\t" 
					<< double(nVariablesNotOK)/double(numCSTRReactors*numComponents)	<< "\t" 
					<< global_time														<< "\t" 
					<< local_time														<< "\t" 
					<< endl;

	// Convergence (1)
	// Non si tratta di una convergenza veramente piena
	if(	nReactorsNotOK    < ::Max(10, numCSTRReactors / 20)		&&
		countIntegrations < ::Max(5, numCSTRReactors / 1000)	&&
		countTutto > memoCount + 300							&&
		resMean <= F1Stop	)
	{
		cout	<< "Get First Completed!"								<< endl
				<< "  ** Convergence crterium: 1"						<< endl
				<< "  ** Mean residual: "				<< resMeanOpt	<< endl
				<< "  ** Mean residual: "				<< maxResOpt	<< endl
				<< "  ** FStop Value: "					<< F1Stop		<< endl
				<< "  ** Number of Calls: "				<< countTutto	<< endl;

		return 0;
	}

	// Convergence (2): convergenza piena
	// La convergenza e' stata raggiunta con successo dal momento che i residui per
	// tutte le variabili e quindi per tutti i reattori sono sufficientemente bassi
	if(nVariablesNotOK == 0)
	{
		cout	<< "Get First Succesfully Completed!"					<< endl
				<< "  ** Convergence criterium: 2"						<< endl
				<< "  ** Mean residual: "				<< resMeanOpt	<< endl
				<< "  ** Mean residual: "				<< maxResOpt	<< endl
				<< "  ** FStop Value: "					<< F1Stop		<< endl
				<< "  ** Number of Calls: "				<< countTutto	<< endl;

		return 1;
	}

	// Viene aggiornato il contatore che tiene conto di quante volte il metodo globale
	// e' stato applicato senza che sia stato necessario ricorrere all'integrazione
	// temporale su un reattore
	if(countTutto >= memoCount + 10 && countIntegrations == 0)
		countNoIntegration++;
	else
		countNoIntegration = 0;

	// Convergence (3)
	// Questo criterio viene applicato se da almeno 4 iterazioni non viene chiamato il
	// sistema ODE per nessun reattore e comunque il residuo medio e' basso
	if(countNoIntegration > 4 && resMean <= F1Stop)
	{
		cout	<< "Get First Completed!"								<< endl
				<< "  ** Convergence criterium: 3"						<< endl
				<< "  ** Mean residual: "				<< resMeanOpt	<< endl
				<< "  ** Mean residual: "				<< maxResOpt	<< endl
				<< "  ** FStop Value: "					<< F1Stop		<< endl
				<< "  ** Number of Calls: "				<< countTutto	<< endl;

		return 0;
	}

	// Updating parameters
	if(countTutto > memoCount)
	{
		// Se le specie o i reattori non OK diminuiti rispetto ai valori precedenti
		if(nVariablesNotOK < minNoOkTot || nReactorsNotOK < minNoOK)
			massFractionsInReactors = massFractionsInReactorsSolution;

		// Updating the minimum number of not OK species
		if(nVariablesNotOK < minNoOkTot)
		{
			countNoMaxRes = 0;
			minNoOkTot = nVariablesNotOK;
			cout	<< "  ** Updating the Number of minimum number of species not OK: " << endl
					<< "    ** New Value: " << minNoOkTot << endl;
		}

		// Updating the minimum number of not OK reactors
		if(nReactorsNotOK < minNoOK)
		{
			countNoMaxRes = 0;
			minNoOK = nReactorsNotOK;
			cout	<< "  ** Updating the Number of minimum number of reactors not OK: " << endl
					<< "    **  New Value: " << minNoOK << endl;
		}

		// Updating the Optimal mean residual
		if(resMean < resMeanOpt)
		{
			countNoMaxRes = 0;
			if(fabs(resMeanOpt - resMean) > .01 * resMean)
				countNoIntegration = 0;

			cout	<< "  ** Updating the Optimal Mean Residual: " << endl
					<< "    ** Difference: " << fabs(resMeanOpt - resMean)/resMean << endl
					<< "    ** Number of iterations without ODE solutions: " << countNoIntegration << endl;

			resMeanOpt = resMean;

			cout	<< "    ** New value: " << resMeanOpt << endl;
		}

		// Updating the Optimal max residual
		if(resMax < maxResOpt)
		{
			countNoMaxRes = 0;
			if(fabs(maxResOpt - resMax) > .01 * resMax)
				countNoIntegration = 0;
			cout	<< "  ** Updating the Optimal Max Residual: " << endl
					<< "    ** Difference: " << fabs(maxResOpt - resMax)/resMax << endl
					<< "    ** Number of iterations without ODE solutions: " << countNoIntegration << endl;

			maxResOpt = resMax;

			cout	<< "    ** New value: " << maxResOpt << endl;
		}
	}

	// Convergence (4)
	// Se il numero di iterazioni globali ha superato il massimo valore consentito e
	// comunque il residuo medio e' basso
	if(countTutto >= maxNumberOfIterations_Global && resMean <= F1Stop)
	{
		cout	<< "Get First Completed!"								<< endl
				<< "  ** Convergence criterium: 4"						<< endl
				<< "  ** Mean residual: "				<< resMeanOpt	<< endl
				<< "  ** Mean residual: "				<< maxResOpt	<< endl
				<< "  ** FStop Value: "					<< F1Stop		<< endl
				<< "  ** Number of Calls: "				<< countTutto	<< endl;

		return 0;
	}

	// Convergence (5)
	// Se non viene aggiornato nessun valore nei criteri di convergenza da piu'
	// di 50 iterazioni globali, esci
	if(countNoMaxRes > 50)
	{
		cout	<< "Get First Completed!"								<< endl
				<< "  ** Convergence criterium: 5"						<< endl
				<< "  ** Mean residual: "				<< resMeanOpt	<< endl
				<< "  ** Mean residual: "				<< maxResOpt	<< endl
				<< "  ** FStop Value: "					<< F1Stop		<< endl
				<< "  ** Number of Calls: "				<< countTutto	<< endl;

		massFractionsInReactorsSolution = massFractionsInReactors;
		return 0;
	}

	if (OnlyODE==true && countTutto>30)
		return 0;

	// If none of these convergence criteria is satisfied, repeat the global iteration
	goto ripetiTutto;
}

// La rete di reattori viene risolta GLOBALMENTE attraverso un sistema ODE
void OpenSMOKE_CSTRNetwork::GetSecond(void)
{
	double F0,F1,maxF0,maxF1;
	double uvariables;
	BzzVectorInt blockDimensions;
	BzzVector initialValues;
	BzzVector residual_Vector;

	countSecondTotal++;
	fResiduals_2 << countSecondTotal << "\t";


	// Memory Allocation
	ChangeDimensions(numCSTRReactors, &blockDimensions);
	ChangeDimensions(numCSTRReactors * numComponents, &residual_Vector);
	ChangeDimensions(numCSTRReactors * numComponents, &dummyVector);

	// Initializzazioni
	blockDimensions = numComponents;
	uvariables = 1. / (numComponents * numCSTRReactors);

	// Initializzazioni
	initialValues = massFractionsInReactorsSolution;
	volumeMolecularWeightT_Vector = volumeMolecularWeightT;
	feed_Vector = feed;

	// Starting Computations
	double start = BzzGetCpuTime();

	// Restituisce i residui per tutte le equazioni in tutti i reattori
	GetAllResiduals(initialValues, residual_Vector);

	// Calcola il residuo medio e massimo
	F0 = residual_Vector.GetSumAbsElements() * uvariables;
	maxF0 = residual_Vector.MaxAbs();
	cout << "Solving Linear System in Second Method (ODE)" << endl;
	cout << "  ** Mean Residual Before Second Method" << F0 << endl;
	cout << "  ** Max  Residual Before Second Method" << maxF0 << endl;
	fResiduals_2 << F0 << "\t" << maxF0 << "\t";

	// Solving the ODE System
	{
		double tStart,tEnd;
		double stopODE;
		BzzVector yMin(initialValues.Size());
		BzzVector yMax(initialValues.Size());
		OpenSMOKE_CSTRNetwork_MyOdeSystemObjectAllCSTR obj(this);

		tStart = 0.;
		tEnd = 10.;
		stopODE = tolAbs * uvariables;
		BzzMatrixSparseLockedByRows Q = Ld;
	//	obj(this);

        string fileSparseMatrix = "Temp/A.tmp";
        char fileSM[200];
        strcpy(fileSM, fileSparseMatrix.c_str());
		BzzOdeSparseStiffObject o(initialValues, tStart, &obj, blockDimensions, &Q, fileSM);
		yMin = 0.;
		yMax = 1.;
		o.SetMinimumConstraints(yMin);
		o.SetMaximumConstraints(yMax);
		o.SetTollAbs(1.e-15);
		o.StopIntegrationBeforeRecalcuatingJacobian(30);
		o.StopIntegrationWhenSumAbsY1IsLessThan(stopODE);
		massFractionsInReactorsSolution_Vector = o(tEnd);

		cout	<< "CPU Time for Solving Linear System in Second Method (ODE): "
				<< BzzGetCpuTime() - start << endl;
	}

	// Calcolo dei nuovi residui dopo l'integrazione del sistema ODE
	GetAllResiduals(massFractionsInReactorsSolution_Vector, residual_Vector);
	F1 = residual_Vector.GetSumAbsElements() * uvariables;
	maxF1 = residual_Vector.MaxAbs();
	cout << "  ** Mean Residual After Second Method (ODE): " << F1 << endl;
	cout << "  ** Max  Residual After Second Method (ODE): " << maxF1 << endl;
	fResiduals_2 << F1 << "\t" << maxF1 << "\t";

	// Accetta la nuova soluzione che sicuramente e' buona in quanto ottenuta
	// attraverso la soluzione di un transitorio reale
	CopyDataFromVector(&massFractionsInReactorsSolution,massFractionsInReactorsSolution_Vector);
	fResiduals_2 << endl;
}

// Viene risolto il sistema lineare GLOBALE associato a tutta la rete di reattori.
// Per fare questa operazione si sfrutta la classe BzzFactoredDiagonalBlockGaussAndMatrixLocked
// che prevede che la matrice del sistema abbia sulla diagonale blocchi anche di diverse
// dimensioni e i termini extradiagonali piccoli se confrontati con i corrispondenti
// sulla diagonale principale. L'algoritmo e' un po' complicato perche' e' necessario
// salvare la fattorizzazione su file onde evitare problemi di allocazioni di memoria
int OpenSMOKE_CSTRNetwork::GetThird(void)
{
	// Indici
	int i,kr;
	BzzLoad load;

	// Numero di iterazioni
	int countThird;

	// Residui
	double F0;		// Residuo medio relativo al vecchio punto
	double F1;		// Residuo medio relativo al nuovo punto
	double maxF0;	// Residuo massimo relativo al vecchio punto
	double maxF1;	// Residuo massimo relativo al nuovo punto

	// Norme delle correzioni di Newton
	double F1dx;	// norma della nuova correzione
	double F0dx;	// norma della vecchia correzione

	// Fattore di normalizzazione dei residui
	double uvariables;

	// Variabili ausiliarie
	double start;
	double hMin,h,sum;

	// Allocazione di memoria
	BzzFactorizedDiagonalBlocksGaussAndMatrixLocked	GlobalMatrix;

	BzzMatrixSparseLockedByRows	ExtraDiagonalTermsMatrix;
	BzzVector ma;
	BzzMatrix reactionsRateInReactors(numCSTRReactors,numComponents);
	BzzMatrix residuals;
	BzzMatrix dRC(numComponents,numComponents);
	BzzVector R(numComponents);
	BzzFactorizedGauss G;
	BzzVector d1;
	BzzVector x0;
	BzzVector x1;

	// Inizializzazione variabili
	countThird = 0;
	uvariables = 1. / (numComponents * numCSTRReactors);
	start = BzzGetCpuTime();

	ripeti_third:

	// Aggiornamento numero di iterazioni
	countThird++;
	countThirdTotal++;
	cout << "Solving Linear System in Third Method - Iteration: " << countThird << endl;
	fResiduals_3 << countThirdTotal << "\t" << countThird << "\t";
	if(memoTemperature == 1)
		load('*', "Temp/MemoTemperatureFunctions.tmp");

	// Viene aperto il file su cui verra' scritta la fattorizzaizone delle matrici
	// che costituiscono i blocchi della matrice del sistema lineare
	// In particolare viene scritto prima di tutto il numero di blocchi, ovvero il
	// numero di reattori della rete
	BzzSave save('*', "Temp/MemoJacobian.tmp");
	save << numCSTRReactors;

	// Ciclo sui singoli reattori
	for(kr = 1;kr <= numCSTRReactors;kr++)
	{
		// Recupero informazioni per il singolo reattore CSTR
		T = temperature[kr];
		P = pressure[kr];
		logT = logTm[kr];
		loguRT = loguRTm[kr];
		cTot = cTotm[kr];

		if(memoTemperature == 1)
		{
			load >> Reactions->uKeq;
			load >> Reactions->k1 >> Reactions->k2;
			load >> Reactions->logFcent;
		}

		// Calcolo delle velocita' di reazione
		massFractionsInReactorsSolution.UseMatrixRowAsVector(kr,&ma);
		GetReactionRatesInSingleReactor(kr, ma, xReactor, cReactor, R);
		reactionsRateInReactors.SetRow(kr,R);

		// Calcolo dello Jacobiano del sistema di equazioni
		// In questo caso lo jacobiano del sistema rappresenta ovviamente la il blocco
		// sulla diagonale della matrice del sistema lineare dal momento che stiamo applicando
		// il metodo di Newton
		GetJacobianForSingleReactor(kr, cReactor, R, dRC);

		// Fattorizzazione dello Jacobiano (Secondo Gauss)
		ReplaceBzzMatrixWithBzzFactorized(&dRC,&G);
		Factorize(&G);

		// Scrittura dello Jacobiano su file
		save << G;
	}

	save.End();
	if(memoTemperature == 1) load.End();


	// La matrice Ld e' una BzzMatrixSparseLockedByRows che di fatto corrisponde
	// esattamente alla matrice Ls, ma senza i termini sulla diagonale principale; questo perche'
	// la definizione della matrice Pcstr prevede che vengano fornite separamente le matrici
	// dei blocchi della diagonale (che in questo caso e' presa dal file) e gli elementi che
	// si trovano fuori
	ExtraDiagonalTermsMatrix = Ld;
	string fileJacobian = "Temp/MemoJacobian.tmp";
	char fileJ[200];
	strcpy(fileJ, fileJacobian.c_str());
	GlobalMatrix(fileJ, &ExtraDiagonalTermsMatrix);

	// Scrittura delle equazioni per tutta la rete di reattori: calcolo dei residui
	Product(Ls,massFractionsInReactorsSolution,&residuals);
	ElementByElementProduct(volumeMolecularWeightT,reactionsRateInReactors,&dummyMatrix);
	Difference(&residuals,dummyMatrix);
	Difference(&residuals,feed);

	// Calcolo del residuo medio e del residuo massimo
	F0 = uvariables * residuals.GetSumAbsElements();
	maxF0 = residuals.MaxAbs();
	cout << " ** Mean Residuals Before Global Newton "<< F0 << endl;
	cout << " ** Max  Residuals Before Global Newton "<< maxF0 << endl;
	fResiduals_3 << F0 << "\t" << maxF0 << "\t";

	double startLocal = BzzGetCpuTime();

	// Viene ricavato il vettore dei termini noti del sistema lineare per poter applicare
	// il metodo di Newton
	residuals.UseMatrixAsVector(&d1);

	// Viene applicato il metodo di Newton per ottenere la correzione
	int fail = GlobalMatrix.GaussIterativeSolve(&d1);
	cout	<< "Flag for Solving Linear System in Third Method: " << fail << endl;
	cout	<< "CPU Time for Solving Linear System in Third Method: "
			<< BzzGetCpuTime() - startLocal << endl;

	if(fail > 0) // TODO
		return 0;

	// Calcolo della norma della correzione
	if(countThird == 1)
		F1dx = F0dx = d1.GetSumAbsElements();
	else
		F1dx = d1.GetSumAbsElements();
	cout << "Newton Correction: F0 Norm " << F0dx << " F1 Norm " << F1dx << endl;

	if(memoTemperature == 1)
		load('*', "Temp/MemoTemperatureFunctions.tmp");
	save('*', "Temp/MemoJacobian.tmp");

	// Ciclo su tutti i singoli reattori
	for(kReactor = 1;kReactor <= numCSTRReactors;kReactor++)
	{
		residuals.UseMatrixRowAsVector(kReactor,&d1);
		massFractionsInReactorsSolution.UseMatrixRowAsVector(kReactor,&x0);

		// La correzione di Newton non puo' essere applicata direttamente ma bisogna
		// controllare che questa non porti ad avere delle frazioni massive negative
		hMin = 1.;
		for(i = 1;i <= numComponents;i++)
		{
			// Se la correzione e' negativa
			if(d1[i] < -1.e-5)
			{
				// Se la specie aveva una frazione massiva negativa
				if(x0[i] == 0.)
				{
					// Se la correzione non e' piccolissima e' meglio lasciare stare
					if(fabs(d1[i]) > 1.e-5)
					{
						hMin = 0.;
						break;
					}
					// Altrimenti togliamo la correzione per questo punto
					else
						d1[i] = 0.;
				}

				// Se invece la frazione massiva non era negativa possiamo sempre stimare
				// un valore percentiuale della correzione, che mi dia una misura di quanto
				// sto spostando il nuovo punto
				// In particolare in questo modo vado a individuare il minimo spostamento
				// pecentuale tra tutti i punti
				else
				{
					h = -x0[i] / d1[i];
					if(h < hMin)
						hMin = h;
				}
			}
		}

		// In questo caso e' meglio abbandonare il metodo
		if(hMin == 0.)
		{
			cout << "TODO: Solving Linear System in Third Method - hMin = 0" << endl;
			BzzWarning("TODO: Solving Linear System in Third Method - hMin = 0");
			return 0;
		}

		else if(hMin <= .999 || hMin >= 1.001)
		{
			BzzWarning("kReactor %d Third hMin %e", kReactor, hMin);
			hMin *= .99;
			Product(hMin,&d1);
		}
		else
			Product(hMin,&d1);

		// Viene calcolato il nuovo punto servendosi della correzione di Newton opportunamente
		// scalata
		massFractionsInReactors.UseMatrixRowAsVector(kReactor,&x1);
		Sum(x0,d1,&x1);

		for(i = 1;i <= numComponents;i++)
			if(x1[i] < 0.)
				x1[i] = 0.;

		// Normalizzaizone delle frazioni massive
		sum = x1.GetSumElements();
		if(sum < .99999999 || sum > 1.000000001)
		{
			if(sum == 0.)
				BzzError("TODO: Solving Linear System in Third Method - Sum = 0.");
			sum = 1. / sum;
			for(i = 1;i <= numComponents;i++)
				x1[i] *= sum;
		}

		// Recupero informazioni per il singolo reattore
		T = temperature[kReactor];
		P = pressure[kReactor];
		logT = logTm[kReactor];
		loguRT = loguRTm[kReactor];
		cTot = cTotm[kReactor];
		if(memoTemperature == 1)
		{
			load >> Reactions->uKeq;
			load >> Reactions->k1 >> Reactions->k2;
			load >> Reactions->logFcent;
		}

		// Calcolo delle velocita' di reazione con il nuovo punto
		massFractionsInReactors.UseMatrixRowAsVector(kReactor,&ma);
		GetReactionRatesInSingleReactor(kReactor, ma, xReactor, cReactor, R);
		reactionsRateInReactors.SetRow(kReactor,R);
	}

	if(memoTemperature == 1) load.End();
	save.End();

	// Calcola i residui delle equazioni per tutti i reattori
	Product(Ls,massFractionsInReactors,&residuals);
	ElementByElementProduct(volumeMolecularWeightT,reactionsRateInReactors,&dummyMatrix);
	Difference(&residuals,dummyMatrix);
	Difference(&residuals,feed);

	// Calcolo del residuo medio e del residuo massimo per il nuovo punto
	F1 = uvariables * residuals.GetSumAbsElements();
	maxF1 = residuals.MaxAbs();
	cout	<< " ** Mean Residuals After Global Newton "<< F1 << endl;
	cout	<< " ** Max  Residuals After Global Newton "<< maxF1 << endl;
	cout	<< " ** CPU Time for Solving Linear System in Third Method: "
			<< BzzGetCpuTime() - start << endl;
	fResiduals_3 << F1 << "\t" << maxF1 << "\t";
	fResiduals_3 << endl;

	// Se la norma della correzione e' piu' piccola di quella
	// del punto precedente e non ho fatto molte iterazioni conviene accettare il
	// nuovo punto e ripetere l'iterazione
	if(F1dx < F0dx && countThird < MaxCountNewtonIterations)
	{
		massFractionsInReactorsSolution = massFractionsInReactors;
		F0dx = F1dx;
		goto ripeti_third;
	}

	// Se il residuo medio oppure quello massimo sono migliorati conviene accettare la
	// soluzione; bisogna fare dei controlli
	// per vedere se questo si e' verificato perche' abbiamo comunque raggiunto la
	// soluzione o se conviene uscire
	if(F1 < F0 || maxF1 < maxF0)
	{
		massFractionsInReactorsSolution = massFractionsInReactors;

		// Se il miglioramento non e' cosi' grande bisogna eseguire degli ulteriori controlli
		if(F1 > .9 * F0 && maxF1 > .9 * maxF0)
		{
			// Se pero' i valori sono sufficientemente bassi vuol dire che la convergenza
			// e' stata raggiunta  che quindi possiamo uscire con successo
			if(F1 < MAX_SUM_F || maxF1 < tolAbs)
				return 1;

			// Se non si tratta della prima iterazioni ci potrebbero essere delle condizioni
			// tali da consentire l'uscita con successo o meno
			else if(countThird > 1)
			{
				cout << "Mean residual F1 "<< F1 << " Mean Residual F1Stop " << F1Stop << endl;
				BzzWarning("TODO: Solving Linear System in Third Method - F1 > MAX_SUM_F");

				if(F1 < F1Stop)
					F1Stop = F1;

				// Se il numero di iterazioni e' eccessivo e il residuo medio non e'
				// ancora basso conviene uscire senza successo; se invece il numero
				// di iterazioni non e' ancora alto conviene tentare ancora una volta
				if(F1 > 1.e-8 && countThird > 5)
					return 0;
				else if(F1 > 1.e-8 && countThird <= 5)
					goto ripeti_third;

				// Se invece il residuo medio e' sufficientemente basso possiamo ritenere
				// comunque di aver raggiunto la soluzione
				return 1;
			}
		}

		goto ripeti_third;
	}

	// Se i residui del nuovo punto non sono milgiori del precedente
	else
	{
		// Se si tratta della prima iterazione e' possibile fare una nuova iterazione
		if(countThird == 1)
		{
			massFractionsInReactorsSolution = massFractionsInReactors;
			goto ripeti_third;
		}
		else
		{
			// Se il residuo e' davvero cattivo conviene uscire, altrimenti se il residuo
			// e' molto buono vuol dire che la soluzione e' accettabile e quindi
			// possiamo uscire con successo
			if(F1 > MAX_SUM_F)
				return 0;
			else
				return 1;
		}
	}
}


void OpenSMOKE_CSTRNetwork::ObjectBzzPrint(void)
{
	int i,j,k;

	BzzVector ma(numComponents);
	BzzMatrix res(numCSTRReactors,numComponents);
	BzzVector sum(numComponents);

	GetResiduals(massFractionsInReactorsSolution,res);
	massFractionsInReactorsSolution.GetRowsSum(&sum);

	::BzzPrint("\nBzzReaction No.%d",whoAmI);
	::BzzPrint("\n\nComponents Number: %d\n",numComponents);
	::BzzPrint("\n\n\nReactions Number: %d\n",numReactions);

	for(j = 1;j <= numReactions;j++)
		::BzzPrint("\n%5d   %s",j,Reactions->strReaction[j]);

	for(i = 1;i <= numCSTRReactors;i++)
	{
		k = cstrSequence[i];

		cTot = pressure[i] / (temperature[i] * R_CTOT);
		massFractionsInReactorsSolution.GetRow(i,&ma);
		wM = Reactions->GetMWFromMassFractions(ma);
		Reactions->GetMoleFractionsFromMassFractionsAndMW(xReactor,ma,wM);
		cReactor = cTot*xReactor;

		cout << endl;
		cout << "========================================" << endl;
		cout << "          Reactor number " << endl;
		cout << "========================================" << endl;

		cout << "T[K]\tP[atm]\tV[m3]\tRate[kg/s]" << endl;
		cout << endl;
		cout << temperature[i] << "\t"
			 << pressure[i] << "\t"
			 << volume[i] << "\t"
			 << massRate[i] << "\t"
			 << endl
			 << endl;

		cout << "Species\tConc\tx\tomega\tResidual" << endl;
		cout << endl;
		for(j = 1;j <= numComponents;j++)
		{
			cout << j << " " << Reactions->names[j]	<< "\t"
				 << cReactor[j]						<< "\t"
				 << xReactor[j]						<< "\t"
				 << ma[j]							<< "\t"
				 << res[i][j]						<< "\t"
				 << endl;
		}
	}

	OutputPrint(0);
}

void OpenSMOKE_CSTRNetwork::OutputPrint(int cicle)
{
	int i,k;
	double flowrate;
	double totalOutFlowRate = 0.;
	BzzVector outFlowRate(numComponents);

	cout.setf(ios::scientific);
	Save('*', "Temp/mass.tmp");
	for(k = 1;k <= numCSTRReactors;k++)
	{
		for(i = 1;i <= numComponents;i++)
		{
			flowrate = massFractionsInReactorsSolution[k][i] * massOutput[k];
			outFlowRate[i] += flowrate;
			totalOutFlowRate += flowrate;
		}
	}
	flowrate = massOutput.GetSumElements();

	if (iEnergyAnalysis == true)
	{
		string fileNameEnthalpy				= "Output/enthalpy/Enthalpy.map";
		string fileNameInletEnthalpy		= "Output/enthalpy/InletEnthalpy.map";
		string fileNameDeltaEnthalpy		= "Output/enthalpy/DeltaEnthalpy.map";
		string fileNameDifferenceEnthalpy	= "Output/enthalpy/DifferenceEnthalpy.map";
		string fileNameErrorEnthalpy		= "Output/enthalpy/ErrorEnthalpy.map";

		EnthalpyAnalysis(	fileNameEnthalpy, fileNameDeltaEnthalpy, fileNameDifferenceEnthalpy, 
							fileNameInletEnthalpy, fileNameErrorEnthalpy);
	}

	cout << endl << endl;
	cout << "===================================================="	<< endl;
	if(cicle == 0)
		cout << "                Final Cicle                     "	<< endl;
	else
		cout << "                Cicle Number " << cicle			<< endl;
		cout << "                CSTR Number "  << numCSTRReactors	<< endl;
	cout << "===================================================="	<< endl;

	cout << "\t\tRate [kg/s]\tMass Fract." << endl << endl;


	for(i = 1;i <= numComponents;i++)
		cout	<< i							<< "\t"
				<< Reactions->names[i]			<< "\t"
				<< outFlowRate[i]				<< "\t"
				<< outFlowRate[i] / flowrate	<< endl;

	cout	<< endl							<< "\t"
			<< "TOTAL"						<< "\t"
			<< totalOutFlowRate				<< "\t"
			<< totalOutFlowRate / flowrate	<< endl
			<< endl;

	cout << "=====================================================" << endl;

	OutputFinalFile("FinalSummary.out");
}

void OpenSMOKE_CSTRNetwork::EnthalpyAnalysis(const string fileNameEnthalpy, const string fileNameDeltaEnthalpy,
										const string fileNameDifferenceEnthalpy, const string fileNameInletEnthalpy,
										const string fileNameErrorEnthalpy
										)
{
	int i, k;

	// Entalpie di miscela in ciascun reattore
	BzzVector Hmix(numCSTRReactors);
	for(k = 1; k<=numCSTRReactors; k++)
	{
		double temperature_cluster;
		BzzVector omega_cluster(numComponents);
		BzzVector h_cluster(numComponents);

		// Temperature and composition of each cluster
		temperature_cluster = temperature[k];							// [K]
		omega_cluster = massFractionsInReactorsSolution.GetRow(k);		// [-]

		// Enthalpy
		Reactions->SpeciesEnthalpy(temperature_cluster);									// [-]
		Product((Constants::R_J_kmol*temperature_cluster), Reactions->h,  &h_cluster);	// [J/kmol]
		ElementByElementProduct(h_cluster,  Reactions->uM, &h_cluster);					// [J/kg]

		// Mixture enthalpy and enthalpy flux
		Hmix[k] = Dot(h_cluster, omega_cluster);	// [J/kg]
	}


	// Calcolo della potenza uscente dal singolo reattore (contributo convettivo)
	Hout_tot = 0.;
	BzzVector Hout(numCSTRReactors);
	for(i = 1;i <= numCSTRReactors;i++)
	{
		double outputMassFlowRate = 0.;
		for(int j=1;j<=FromClusterToCluster_Index[i].Size();j++)
		{
			outputMassFlowRate += FromClusterToCluster_MassFlowRate[i][j];
			if (FromClusterToCluster_Index[i][j] == 0)	
				Hout_tot += FromClusterToCluster_MassFlowRate[i][j]*Hmix[i];
		}
		Hout[i] = outputMassFlowRate*Hmix[i];
	}

	// Calcolo della potenza entrante nel singolo reattore (contributo convettivo)
	BzzVector Hin(numCSTRReactors);
	for(i = 1;i <= numCSTRReactors;i++)
	{
		for(int k = 1; k<=numCSTRReactors; k++)
			for(int j=1;j<=FromClusterToCluster_Index[k].Size();j++)
				if (FromClusterToCluster_Index[k][j] == i)
					Hin[i] += Hmix[k]*FromClusterToCluster_MassFlowRate[k][j];
		Hin[i] += Hinput[i];
	}

	// Calcolo delle differenze di potenza (contributo convettivo)
	BzzVector DH(numCSTRReactors);
	for(k = 1; k<=numCSTRReactors; k++)
		DH[k] = Hout[k] - Hin[k];

	
	// Potenza totale entrante nel sistema
	Hin_tot = Hinput.GetSumElements();	// [J/s]


	cout << "Energy analysis" << endl;
	cout << "------------------------------------------------------"		<< endl;
	cout << "Htot Input:  " << Hin_tot*1.e-3				<< " [kW]"		<< endl;
	cout << "Htot Output: " << Hout_tot*1.e-3				<< " [kW]"		<< endl;
	cout << "Hloss:       " << (Hin_tot-Hout_tot)*1.e-3		<< " [kW]"		<< endl;

	// Enthalpy of each reactor
	{
		ofstream fEnthalpy;
		openOutputFileAndControl(fEnthalpy, fileNameEnthalpy);
		fEnthalpy.setf(ios::scientific);
		for (i=1; i<=originalNumCSTRReactors; i++)
			fEnthalpy << Hmix[giveClusterIndexFromCellIndex[i]] << endl;					// [J/kg]
		fEnthalpy.close();
	}
}

void OpenSMOKE_CSTRNetwork::Maps()
{
    int i,j,k, kSpecies;
    double MWmix;

	cout << "Writing maps on file... " << endl;

    vector<string> listOfSpecies;
	for(k=1;k<=Reactions->NumberOfSpecies();k++)
		listOfSpecies.push_back(Reactions->names[k]);

    BzzVector mass_fractions(massFractionsInReactors.Columns());
    BzzVector mole_fractions(massFractionsInReactors.Columns());
    BzzMatrix moleFractionsInReactors(massFractionsInReactors.Rows(),massFractionsInReactors.Columns());

    for(k=1;k<=massFractionsInReactors.Rows();k++)
    {
        mass_fractions = massFractionsInReactors.GetRow(k);
        Reactions->GetMWAndMoleFractionsFromMassFractions(MWmix, mole_fractions, mass_fractions);
        moleFractionsInReactors.SetRow(k, mole_fractions);
    }

    for(k=0;k<int(listOfSpecies.size());k++)
    {
        ofstream fOutputMass;
        ofstream fOutputMole;

        openOutputFileAndControl(fOutputMass, "Output/maps/MassFraction_" + listOfSpecies[k] + ".map");
        openOutputFileAndControl(fOutputMole, "Output/maps/MoleFraction_" + listOfSpecies[k] + ".map");
        fOutputMass.setf(ios::scientific);
        fOutputMole.setf(ios::scientific);

        for(i = 1;i <= originalNumCSTRReactors;i++)
        {
            j = giveClusterIndexFromCellIndex[i];
            kSpecies = Reactions->recognize_species_without_exit(listOfSpecies[k]);
			if (kSpecies>0)
			{
				fOutputMass << massFractionsInReactors[j][kSpecies] << endl;
				fOutputMole << moleFractionsInReactors[j][kSpecies] << endl;
			}
			else
			{
				fOutputMass << 0. << endl;
				fOutputMole << 0. << endl;
			}
        }

        fOutputMass.close();
        fOutputMole.close();
    }

	// Mass Fractions
	/*
	{
		ofstream fOutputRoryMassFractions;
		openOutputFileAndControl(fOutputRoryMassFractions, "Output/maps/MassFractions.rory");
        fOutputRoryMassFractions.setf(ios::scientific);

		fOutputRoryMassFractions << "#\t";
		fOutputRoryMassFractions << "H2\t";
		fOutputRoryMassFractions << "O2\t";
		fOutputRoryMassFractions << "OH\t";
		fOutputRoryMassFractions << "H2O\t";
		fOutputRoryMassFractions << "CH4\t";
		fOutputRoryMassFractions << "CO\t";
		fOutputRoryMassFractions << "CO2\t";
		fOutputRoryMassFractions << "N2\t";
		fOutputRoryMassFractions << endl;

		for(i = 1;i <= originalNumCSTRReactors;i++)
		{
				int j = giveClusterIndexFromCellIndex[i];

				fOutputRoryMassFractions << i << "\t"; 
				fOutputRoryMassFractions << massFractionsInReactors[j][Reactions->recognize_species("H2")]  << "\t";
				fOutputRoryMassFractions << massFractionsInReactors[j][Reactions->recognize_species("O2")]  << "\t";
				fOutputRoryMassFractions << massFractionsInReactors[j][Reactions->recognize_species("OH")]  << "\t";
				fOutputRoryMassFractions << massFractionsInReactors[j][Reactions->recognize_species("H2O")] << "\t";
				fOutputRoryMassFractions << massFractionsInReactors[j][Reactions->recognize_species("CH4")] << "\t";
				fOutputRoryMassFractions << massFractionsInReactors[j][Reactions->recognize_species("CO")]  << "\t";
				fOutputRoryMassFractions << massFractionsInReactors[j][Reactions->recognize_species("CO2")] << "\t";
				fOutputRoryMassFractions << massFractionsInReactors[j][Reactions->recognize_species("N2")]  << "\t";
				fOutputRoryMassFractions << endl;
		}

		fOutputRoryMassFractions.close();
	}

	// Mole Fractions
	{
		ofstream fOutputRoryMoleFractions;
		openOutputFileAndControl(fOutputRoryMoleFractions, "Output/maps/MoleFractions.rory");
        fOutputRoryMoleFractions.setf(ios::scientific);

		fOutputRoryMoleFractions << "#\t";
		fOutputRoryMoleFractions << "H2\t";
		fOutputRoryMoleFractions << "O2\t";
		fOutputRoryMoleFractions << "OH\t";
		fOutputRoryMoleFractions << "H2O\t";
		fOutputRoryMoleFractions << "CH4\t";
		fOutputRoryMoleFractions << "CO\t";
		fOutputRoryMoleFractions << "CO2\t";
		fOutputRoryMoleFractions << "N2\t";
		fOutputRoryMoleFractions << endl;

		for(i = 1;i <= originalNumCSTRReactors;i++)
		{
				int j = giveClusterIndexFromCellIndex[i];

				fOutputRoryMoleFractions << i << "\t"; 
				fOutputRoryMoleFractions << moleFractionsInReactors[j][Reactions->recognize_species("H2")]  << "\t";
				fOutputRoryMoleFractions << moleFractionsInReactors[j][Reactions->recognize_species("O2")]  << "\t";
				fOutputRoryMoleFractions << moleFractionsInReactors[j][Reactions->recognize_species("OH")]  << "\t";
				fOutputRoryMoleFractions << moleFractionsInReactors[j][Reactions->recognize_species("H2O")] << "\t";
				fOutputRoryMoleFractions << moleFractionsInReactors[j][Reactions->recognize_species("CH4")] << "\t";
				fOutputRoryMoleFractions << moleFractionsInReactors[j][Reactions->recognize_species("CO")]  << "\t";
				fOutputRoryMoleFractions << moleFractionsInReactors[j][Reactions->recognize_species("CO2")] << "\t";
				fOutputRoryMoleFractions << moleFractionsInReactors[j][Reactions->recognize_species("N2")]  << "\t";
				fOutputRoryMoleFractions << endl;
		}

		fOutputRoryMoleFractions.close();
	}
	*/
}

void OpenSMOKE_CSTRNetwork::OutputFinalFile(string fileName)
{
	int i,k;
	double flowrate;
	ofstream fFinal;
	openOutputFileAndControl(fFinal, fileName );
	fFinal.setf(ios::scientific);

	double massFlowRate = 0.;
	double moleFlowRate = 0.;
	BzzVector outMassFlowRate(numComponents);
	BzzVector outMoleFlowRate(numComponents);
	BzzVector outW(numComponents);
	BzzVector outX(numComponents);

	for(k = 1;k <= numCSTRReactors;k++)
	{
		for(i = 1;i <= numComponents;i++)
		{
			flowrate = massFractionsInReactorsSolution[k][i] * massOutput[k];
			outMassFlowRate[i]	+= flowrate;
			massFlowRate		+= flowrate;

			flowrate			/= Reactions->M[i];
			outMoleFlowRate[i]	+= flowrate;
			moleFlowRate		+= flowrate;
		}
	}

	outW = outMassFlowRate / massFlowRate;
	outX = outMoleFlowRate / moleFlowRate;

	cout << "\t\tRate [kg/s]\tMass Fract." << endl << endl;

	fFinal << "Number of CFD Cells:     "	<< originalNumCSTRReactors	<<endl;
	fFinal << "Number of Reactor Cells: "	<< numCSTRReactors			<<endl;
	fFinal << "Clustering ratio:        "	<< double(numCSTRReactors)/double(originalNumCSTRReactors)*100.	<< " %" <<endl;
	fFinal << endl << endl;

	fFinal << "#     "			        << "\t";
	fFinal << "Name         "	        << "\t\t";
	fFinal << "MassFraction "	        << "\t\t";
	fFinal << "MoleFraction "	        << "\t";
	fFinal << "MassFlowRate[kg/s] "	    << "\t";
	fFinal << "MoleFlowRate[kmol/s] "	<< "\t";
	fFinal << endl;
	fFinal << "------------------------------------------------------------------------------" << endl;
	for(i = 1;i <= numComponents;i++)
	{
		fFinal	<< i							<< "\t\t";
		fFinal	<< Reactions->names[i]			<< "\t\t";
		fFinal	<< outW[i]						<< "\t\t";
		fFinal	<< outX[i]						<< "\t\t";
		fFinal	<< outMassFlowRate[i]			<< "\t\t";
		fFinal	<< outMoleFlowRate[i]			<< "\t\t";
		fFinal  << endl;
	}
	fFinal << endl;
	fFinal << endl;

	int iO2  = Reactions->recognize_species("O2");
	int iCO2 = Reactions->recognize_species("CO2");
	int iH2O = Reactions->recognize_species("H2O");
	int iCO  = Reactions->recognize_species("CO");


		int nMain = 4;
		BzzVectorInt indexMain(nMain);

		k=1;
		indexMain(k++) = iO2;
		indexMain(k++) = iCO2;
		indexMain(k++) = iH2O;
		indexMain(k++) = iCO;

		BzzVector Main_W(nMain);
		BzzVector Main_X(nMain);

		for(i=1;i<=nMain;i++)
		{
			Main_W[i]	= outW[indexMain[i]];
			Main_X[i]	= outX[indexMain[i]];
		}

		fFinal << "#        O2             CO2           H2O           CO" << endl;

		fFinal << "Mass:   ";
		for (i=1;i<=nMain;i++)	fFinal << Main_W[i] << "\t";
		fFinal << endl;

		fFinal << "Mole:   ";
		for (i=1;i<=nMain;i++)	fFinal << Main_X[i] << "\t";
		fFinal << endl;

		fFinal << endl << endl;


		int iNO  = Reactions->recognize_species_without_exit("NO");
		int iNO2 = Reactions->recognize_species_without_exit("NO2");
		int iN2O = Reactions->recognize_species_without_exit("N2O");
		int iHCN = Reactions->recognize_species_without_exit("HCN");
		int iNCO = Reactions->recognize_species_without_exit("NCO");

		int nNOX = 5;
		BzzVectorInt indexNOX(nNOX);

		k=1;
		indexNOX(k++) = iNO;
		indexNOX(k++) = iNO2;
		indexNOX(k++) = iN2O;
		indexNOX(k++) = iHCN;
		indexNOX(k++) = iNCO;

		BzzVector NOX_ppm(nNOX);
		BzzVector NOX_dry(nNOX);
		BzzVector NOX_3O2(nNOX);
		BzzVector NOX_6O2(nNOX);

		for(i=1;i<=nNOX;i++)
		{
			if (indexNOX[i] > 0)
			{
				NOX_ppm[i] = 1.e6*outX[indexNOX[i]];
				NOX_dry[i] = NOX_ppm[i] / (1.-outX[iH2O]);
				NOX_3O2[i] = NOX_ppm[i] * (0.21-outX[iO2]) / 0.18;
				NOX_6O2[i] = NOX_ppm[i] * (0.21-outX[iO2]) / 0.15;
			}
			else
			{
				NOX_ppm[i] = 0.;
				NOX_dry[i] = 0.;
				NOX_3O2[i] = 0.;
				NOX_6O2[i] = 0.;
			}
		}

		fFinal << "#[mole]       NO             NO2           N2O           HCN           NCO" << endl;

        fFinal << "Out[ppm]:  ";
		for (i=1;i<=nNOX;i++)	fFinal << NOX_ppm[i] << "\t";
		fFinal<< endl;

		fFinal << "Dry[ppm]:  ";
		for (i=1;i<=nNOX;i++)	fFinal << NOX_dry[i] << "\t";
		fFinal<< endl;

		fFinal << "3%O2[ppm]: ";
		for (i=1;i<=nNOX;i++)	fFinal << NOX_3O2[i] << "\t";
		fFinal << endl;

		fFinal << "6%O2[ppm]: ";
		for (i=1;i<=nNOX;i++)	fFinal << NOX_6O2[i] << "\t";
		fFinal << endl;

		fFinal << endl << endl;

		fFinal << "Energy analysis" << endl;
		fFinal << "------------------------------------------------------"		<< endl;
		fFinal << "htot Input:  " << Hin_tot/massFlowRate*1.e-3		<< " [kJ/kg]"	<< endl;
		fFinal << "htot Output: " << Hout_tot/massFlowRate*1.e-3	<< " [kJ/kg]"	<< endl;
		fFinal << "Htot Input:  " << Hin_tot*1.e-3					<< " [kW]"		<< endl;
		fFinal << "Htot Output: " << Hout_tot*1.e-3					<< " [kW]"		<< endl;
		fFinal << "Hloss:       " << (Hin_tot-Hout_tot)*1.e-3		<< " [kW]"		<< endl;

	fFinal.close();
}

void OpenSMOKE_CSTRNetwork::GiveMeOutputLabel(ofstream &fOutput)
{
	int i;
	fOutput.setf(ios::scientific);

	fOutput << "Index(1)     " << "\t";
	fOutput << "Reactors(2)  " << "\t";
	fOutput << "TotReact(3)  " << "\t";
	fOutput << "Ratio(4)     " << "\t";
	fOutput << "NO_ppm(5)    " << "\t"
			<< "NO2_ppm(6)   " << "\t"
			<< "N2O_ppm(7)   " << "\t"
			<< "HCN_ppm(8)   " << "\t"
			<< "NCO_ppm(9)   " << "\t"
			<< "NO_dry(10)    " << "\t"
			<< "NO2_dry(11)   " << "\t"
			<< "N2O_dry(12)  " << "\t"
			<< "HCN_dry(13)  " << "\t"
			<< "NCO_dry(14)  " << "\t"
			<< "NO_3O2(15)   " << "\t"
			<< "NO2_3O2(16)  " << "\t"
			<< "N2O_3O2(17)  " << "\t"
			<< "HCN_3O2(18)  " << "\t"
			<< "NCO_3O2(19)  " << "\t"
			<< "NO_6O2(20)   " << "\t"
			<< "NO2_6O2(21)  " << "\t"
			<< "N2O_6O2(22)  " << "\t"
			<< "HCN_6O2(23)  " << "\t"
			<< "NCO_6O2(24)  " << "\t";
	
	int count = 25;
	for(i = 1;i <= numComponents;i++)	fOutput << Reactions->names[i] << "(x" << count++ << ")\t";
	for(i = 1;i <= numComponents;i++)	fOutput << Reactions->names[i] << "(w" << count++ << ")\t";
	fOutput << endl << endl;
}

void OpenSMOKE_CSTRNetwork::GiveMeOutput(ofstream &fOutput, const int index)
{
	int i,k;
	double flowrate;
	double massFlowRate = 0.;
	double moleFlowRate = 0.;
	BzzVector outMassFlowRate(numComponents);
	BzzVector outMoleFlowRate(numComponents);
	BzzVector outW(numComponents);
	BzzVector outX(numComponents);

	for(k = 1;k <= numCSTRReactors;k++)
	{
		for(i = 1;i <= numComponents;i++)
		{
			flowrate = massFractionsInReactorsSolution[k][i] * massOutput[k];
			outMassFlowRate[i]	+= flowrate;
			massFlowRate		+= flowrate;

			flowrate			/= Reactions->M[i];
			outMoleFlowRate[i]	+= flowrate;
			moleFlowRate		+= flowrate;
		}
	}

	outW = outMassFlowRate / massFlowRate;
	outX = outMoleFlowRate / moleFlowRate;

	int iO2  = Reactions->recognize_species("O2");
	int iCO2 = Reactions->recognize_species("CO2");
	int iH2O = Reactions->recognize_species("H2O");
	int iCO  = Reactions->recognize_species("CO");

	int nMain = 4;
	BzzVectorInt indexMain(nMain);
	BzzVector Main_W(nMain);
	BzzVector Main_X(nMain);

	k=1;
	indexMain(k++) = iO2;
	indexMain(k++) = iCO2;
	indexMain(k++) = iH2O;
	indexMain(k++) = iCO;

	for(i=1;i<=nMain;i++)
	{
		Main_W[i]	= outW[indexMain[i]];
		Main_X[i]	= outX[indexMain[i]];
	}

	int nNOX = 5;
	int iNO  = Reactions->recognize_species_without_exit("NO");
	int iNO2 = Reactions->recognize_species_without_exit("NO2");
	int iN2O = Reactions->recognize_species_without_exit("N2O");
	int iHCN = Reactions->recognize_species_without_exit("HCN");
	int iNCO = Reactions->recognize_species_without_exit("NCO");

	BzzVectorInt	indexNOX(nNOX);
	BzzVector		NOX_ppm(nNOX);
	BzzVector		NOX_dry(nNOX);
	BzzVector		NOX_3O2(nNOX);
	BzzVector		NOX_6O2(nNOX);

	k=1;
	indexNOX(k++) = iNO;
	indexNOX(k++) = iNO2;
	indexNOX(k++) = iN2O;
	indexNOX(k++) = iHCN;
	indexNOX(k++) = iNCO;

	for(i=1;i<=nNOX;i++)
	{
		if(indexNOX[i]>0)
		{
			NOX_ppm[i] = 1.e6*outX[indexNOX[i]];
			NOX_dry[i] = NOX_ppm[i] / (1.-outX[iH2O]);
			NOX_3O2[i] = NOX_ppm[i] * (0.21-outX[iO2]) / 0.18;
			NOX_6O2[i] = NOX_ppm[i] * (0.21-outX[iO2]) / 0.15;
		}
		else
		{
			NOX_ppm[i] = 0.;
			NOX_dry[i] = 0.;
			NOX_3O2[i] = 0.;
			NOX_6O2[i] = 0.;
		}
	}

	// ------------------------------------------------------------------------------ //
	//								Write on file                                     //
	// ------------------------------------------------------------------------------ //

	fOutput << index					<< "\t";
	fOutput << numCSTRReactors			<< "\t";
	fOutput << originalNumCSTRReactors	<< "\t";
	fOutput << double(numCSTRReactors)/double(originalNumCSTRReactors)	<< "\t";
	
	for (i=1;i<=nNOX;i++)	fOutput << NOX_ppm[i]		<< "\t";
	for (i=1;i<=nNOX;i++)	fOutput << NOX_dry[i]		<< "\t";
	for (i=1;i<=nNOX;i++)	fOutput << NOX_3O2[i]		<< "\t";
	for (i=1;i<=nNOX;i++)	fOutput << NOX_6O2[i]		<< "\t";
	
	for(i = 1;i <= numComponents;i++)	fOutput << outX[i]		<< "\t";
	for(i = 1;i <= numComponents;i++)	fOutput << outW[i]		<< "\t";
	
	fOutput << endl;
}

void OpenSMOKE_CSTRNetwork::CleanTemperatureVariance(const int kind)
{
	if (kind == 0)			// None
	{	
	}
	else if (kind == 1)		// Sin
	{
		double min_csi =  1.e16;
		double max_csi = -1.e16;
		for(int k=1;k<=numCSTRReactors;k++)
		{
			double T = temperature[k];
			double csi = (T-Tmin)/DeltaT[k];
			double qanm_max = 0.99 * csi*csi*DeltaT[k]*DeltaT[k]/T/T;
			if (qanm[k]>qanm_max)	qanm[k] = qanm_max;

			if (csi < min_csi)	min_csi = csi;
			if (csi > max_csi)	max_csi = csi;
		}

		cout << "Minimum csi: " << min_csi << endl;
		cout << "Maximum csi: " << max_csi << endl;
	}
	else if (kind >= 2)		// Beta || Gaussian
	{
		int iCountMin = 0;
		int iCountMax = 0;
		double min_csi =  1.e16;
		double max_csi = -1.e16;
		double min_variance_ratio =  1.e16;
		double max_variance_ratio = -1.e16;

		for(int k=1;k<=numCSTRReactors;k++)
		{
			double T = temperature[k];
			double csi = (T-Tmin)/DeltaT[k];
			double max_g = csi*(1.-csi);
			double max_allowed_g = Constants::Kg_max*max_g;
			double min_allowed_g = Constants::Kg_min*max_g;

			double qanm_max = max_allowed_g * 2.*BzzPow2(DeltaT[k]/T);
			double qanm_min = min_allowed_g * 2.*BzzPow2(DeltaT[k]/T);

			if (qanm[k]>qanm_max)	{	qanm[k] = qanm_max;	iCountMax++; }
			if (qanm[k]<qanm_min)	{	qanm[k] = qanm_min;	iCountMin++; }
		
			if (csi < min_csi)	min_csi = csi;
			if (csi > max_csi)	max_csi = csi;
		}

		cout << "Minimum csi: " << min_csi << endl;
		cout << "Maximum csi: " << max_csi << endl;
		cout << "Number of corrections of minimum variance: " << iCountMin << endl;
		cout << "Number of corrections of maximum variance: " << iCountMax << endl;

		if (min_csi < Constants::csi_min)	ErrorMessage("Csi minimum is too small...");
		if (max_csi > Constants::csi_max)	ErrorMessage("Csi maximum is too large...");
	}
}


void OpenSMOKE_CSTRNetwork::MemoTemperatureFunctions(const int kind)
{
	int kCSTR;
	double uRT;

	cout << "        Memory allocation... " << endl;
	BzzVector u_temperature(numCSTRReactors);
	BzzVector log_temperature(numCSTRReactors);

	BzzMatrix matrix_correction_k1(numCSTRReactors,numReactions);
	BzzMatrix matrix_correction_uKeq(numCSTRReactors,numReactions);
	BzzMatrix matrix_correction_k2(numCSTRReactors,numReactions);

	cout << "        Preparation... " << endl;
	for(kCSTR = 1;kCSTR <= numCSTRReactors;kCSTR++)
	{
		u_temperature[kCSTR] = 1. / temperature[kCSTR];
		log_temperature[kCSTR] = log(temperature[kCSTR]);

		uRT = u_temperature[kCSTR] * UR_CTOT;
		loguRT = log(uRT);

		logTm[kCSTR] = log_temperature[kCSTR];
		loguRTm[kCSTR] = loguRT;
		cTotm[kCSTR] = pressure[kCSTR] * uRT;
	}

	cout << "        Step A: Kinetic constants... "		<< endl;
	cout << "                Map setup... "				<< endl;
	Reactions->InitializeMap(numCSTRReactors);
	cout << "                Kinetic parameters... "	<< endl;
	
	if (pressure.Max()*101325. > 101325.*100.)
		ErrorMessage("Pressure is over 100 atm");
	BzzVector pressure_Pa = pressure;
	pressure_Pa *= 101325.;

	Reactions->ComputeKineticParameters_map(temperature, cTotm, pressure_Pa);

	ChangeDimensions(numCSTRReactors, &Tmax);
	ChangeDimensions(numCSTRReactors, &DeltaT);

	Tmin   = TminGlobal-20.;
	Tmax   = TmaxGlobal+Fluctuations_DeltaMax;

	// User defined maximum temperature
	if (Fluctuations_TMax != 0.)
		Tmax = Fluctuations_TMax;
	
	// User defined Correction coefficient (local)
	if (Fluctuations_CcMax != 0.)
		for(kCSTR = 1;kCSTR <= numCSTRReactors;kCSTR++)
			Tmax[kCSTR] = temperature[kCSTR]*Fluctuations_CcMax;

	// Max-Min temperature difference
	for(kCSTR = 1;kCSTR <= numCSTRReactors;kCSTR++)
		DeltaT[kCSTR] = (Tmax[kCSTR]-Tmin);
	
	cout << "                Clean Temperature Variance... " << endl;
	CleanTemperatureVariance(kind);

	cout << " * Minimum Global Temperature:  " << TminGlobal << " K" << endl;
	cout << " * Maximum Global Temperature:  " << TmaxGlobal << " K" << endl;
	cout << " * Minimum Allowed Temperature: " << Tmin		 << " K" << endl;
	cout << " * Maximum Global Temperature:  " << Tmax.Max() << " K" << endl;

	// Fluctuating reactions
	{
		BzzVectorInt fluctuating_reactions;
		ChangeDimensions(numReactions, &switchReactions);
		if (list_of_fluctuating_species.size() > 0)
		{
			fluctuating_reactions = Reactions->kinetics.ReactionIndices(list_of_fluctuating_species);
			
			for(int j=1;j<=numReactions;j++)
				for(int k=1;k<=fluctuating_reactions.Size();k++)
					if (fluctuating_reactions[k] == j)
					{
						switchReactions[j] = 1;
						break;
					}
		}
		else
		{
			switchReactions = 1;
		}	
	}

	openOutputFileAndControl(fWarning, "WarningCorrections.log");
	fWarning.setf(ios::scientific);

	if (kind == 0)
	{
		cout << "        Step B: Correction coefficients [none]... " << endl;
		CorrectionCoefficient_None(	u_temperature, log_temperature,
									matrix_correction_k1, matrix_correction_uKeq, matrix_correction_k2);
	}
	else if (kind == 1)
	{
		cout << "        Step B: Correction coefficients [sin expansion]... " << endl;
		CorrectionCoefficient_SinExpansion(	u_temperature, log_temperature,
										matrix_correction_k1, matrix_correction_uKeq, matrix_correction_k2);
	}
	else if (kind == 2)
	{
		cout << "        Step B: Correction coefficients [double delta dirac]... " << endl;
		CorrectionCoefficient_DeltaDirac(	u_temperature, log_temperature,
										matrix_correction_k1, matrix_correction_uKeq, matrix_correction_k2);
	}
	else if (kind == 3)
	{
		cout << "        Step B: Correction coefficients [beta PDF]... " << endl;
		CorrectionCoefficient_BetaPDF(u_temperature, log_temperature,
										matrix_correction_k1, matrix_correction_uKeq, matrix_correction_k2);
	}
	else if (kind == 4)
	{
		cout << "        Step B: Correction coefficients [clipped gaussian PDF]... " << endl;
		CorrectionCoefficient_ClippedGaussianPDF(	u_temperature, log_temperature,
													matrix_correction_k1, matrix_correction_uKeq, matrix_correction_k2);
	}

	fWarning.close();

	cout << "        Step C: Correction mapping... " << endl;
	Reactions->CorrectKineticParameters_map(matrix_correction_k1, matrix_correction_uKeq, matrix_correction_k2);
	
	cout << "        Step D: Kinetic temperature maps... " << endl;
	{
		{
			ofstream fOutput;
			openOutputFileAndControl(fOutput, "Output/maps/Tk_20000cal_mol.map");
			fOutput.setf(ios::scientific);
			for (int i=1; i<=originalNumCSTRReactors; i++)
				fOutput << Tk_20000[giveClusterIndexFromCellIndex[i]] << endl;
			fOutput.close();
		}
		{
			ofstream fOutput;
			openOutputFileAndControl(fOutput, "Output/maps/Tk_40000cal_mol.map");
			fOutput.setf(ios::scientific);
			for (int i=1; i<=originalNumCSTRReactors; i++)
				fOutput << Tk_40000[giveClusterIndexFromCellIndex[i]] << endl;
			fOutput.close();
		}
		{
			ofstream fOutput;
			openOutputFileAndControl(fOutput, "Output/maps/Tk_60000cal_mol.map");
			fOutput.setf(ios::scientific);
			for (int i=1; i<=originalNumCSTRReactors; i++)
				fOutput << Tk_60000[giveClusterIndexFromCellIndex[i]] << endl;
			fOutput.close();
		}
	}

	// Writing on File
	if(memoTemperature == 1)
	{
		BzzSave save('*', "Temp/MemoTemperatureFunctions.tmp");
		for(kCSTR = 1;kCSTR <= numCSTRReactors;kCSTR++)
		{
			BzzVector auxiliar;

			auxiliar = Reactions->uKeq_map.GetRow(kCSTR);			save << auxiliar;
			auxiliar = Reactions->k1_map.GetRow(kCSTR);				save << auxiliar;
			auxiliar = Reactions->k2_map.GetRow(kCSTR);				save << auxiliar;
			auxiliar = Reactions->logFcent_map.GetRow(kCSTR);	save << auxiliar;
		}
		// Closing File
		save.End();
	}

	#if SYMBOLIC_KINETICS==1

	if (iAnalyticalJacobian == 1)
	{
		cout << "        Step E: Initialize Analytical Jacobian... " << endl; 
		
		if (analyticalJacobian == GRI12)
			for (int kCSTR=1;kCSTR<=numCSTRReactors;kCSTR++)
				reactor[kCSTR] = new OpenSMOKE_SymbolicKinetics_GRI12();
		else
			ErrorMessage("No Symbolic Jacobian available!");
		/*
		if (analyticalJacobian == POLIMI_C1C3HTNOX_0810)
			for (int kCSTR=1;kCSTR<=numCSTRReactors;kCSTR++)
				reactor[kCSTR] = new OpenSMOKE_SymbolicKinetics_Polimi_C1C3HTNOX_0810();
		else if (analyticalJacobian == POLIMI_C1C3HTNOX_AVIO)
			for (int kCSTR=1;kCSTR<=numCSTRReactors;kCSTR++)
				reactor[kCSTR] = new OpenSMOKE_SymbolicKinetics_Polimi_C1C3HTNOX_AVIO();
		else if (analyticalJacobian == FLUENT_GLARBORG152)
			for (int kCSTR=1;kCSTR<=numCSTRReactors;kCSTR++)
				reactor[kCSTR] = new OpenSMOKE_SymbolicKinetics_Fluent_Glarborg152();

		else if (analyticalJacobian == FLUENT_DRM22_POLIMI)
			for (int kCSTR=1;kCSTR<=numCSTRReactors;kCSTR++)
				reactor[kCSTR] = new OpenSMOKE_SymbolicKinetics_Fluent_DRM22_Polimi();
		else if (analyticalJacobian == FLUENT_DRM22_POLIMI_NOX)
			for (int kCSTR=1;kCSTR<=numCSTRReactors;kCSTR++)
				reactor[kCSTR] = new OpenSMOKE_SymbolicKinetics_Fluent_DRM22_Polimi_NOX();
		else if (analyticalJacobian == FLUENT_DRM22_POLIMI_THERMALNOX)
			for (int kCSTR=1;kCSTR<=numCSTRReactors;kCSTR++)
				reactor[kCSTR] = new OpenSMOKE_SymbolicKinetics_Fluent_DRM22_Polimi_ThermalNOX();

		else if (analyticalJacobian == GRI30)
			for (int kCSTR=1;kCSTR<=numCSTRReactors;kCSTR++)
				reactor[kCSTR] = new OpenSMOKE_SymbolicKinetics_GRI30();
		else if (analyticalJacobian == SANDIEGO_AVIO)
			for (int kCSTR=1;kCSTR<=numCSTRReactors;kCSTR++)
				reactor[kCSTR] = new OpenSMOKE_SymbolicKinetics_SanDiego_AVIO();
		else
			ErrorMessage("No Symbolic Jacobian available!");
			*/
		if (reactor[1]->NC != Reactions->NumberOfSpecies())		ErrorMessage("The number of species in the Symbolic Kinetics is wrong");
		if (reactor[1]->NR != Reactions->NumberOfReactions())	ErrorMessage("The number of reactions in the Symbolic Kinetics is wrong");
		
		for(kCSTR = 1;kCSTR <= numCSTRReactors;kCSTR++)
		{
			BzzVector k1_map			= Reactions->k1_map.GetRow(kCSTR);
			BzzVector uKeq_map			= Reactions->uKeq_map.GetRow(kCSTR);
			BzzVector logFcent_map		= Reactions->logFcent_map.GetRow(kCSTR);
			BzzVector k2_map			= Reactions->k2_map.GetRow(kCSTR);
			reactor[kCSTR]->assignKineticConstants(k1_map, uKeq_map, logFcent_map, k2_map);
		}
	}
	#endif
}

void OpenSMOKE_CSTRNetwork::WriteCorrectionMap(const int kReaction, BzzMatrix &matrix_correction_coefficient, const string fileName)
{
	ofstream fOutputCorrection;
	string name = "Temp/" + fileName;
	openOutputFileAndControl(fOutputCorrection, name);
	fOutputCorrection.setf(ios::scientific);

	for(int i = 1;i <= originalNumCSTRReactors;i++)
	{
		int j = giveClusterIndexFromCellIndex[i];
        fOutputCorrection << matrix_correction_coefficient[j][kReaction] << endl;
    }
	fOutputCorrection.close();
}

void OpenSMOKE_CSTRNetwork::GetReactionsRateInAllReactorsFromMassFractions
		(
			BzzMatrix &massFractionsInReactors, // reactors as rows
			BzzMatrix &reactionsRateInReactors  // reactors as rows
		)
{
	int kr;
	BzzLoad load;

	BzzVector R;
	BzzVector ma;

	ChangeDimensions(numComponents, &R);
	ChangeDimensions(numCSTRReactors,numComponents, &reactionsRateInReactors);

	if(memoTemperature == 1)
		load('*', "Temp/MemoTemperatureFunctions.tmp");

	for(kr = 1;kr <= numCSTRReactors;kr++)
	{
		T = temperature[kr];
		P = pressure[kr];

		logT = logTm[kr];
		loguRT = loguRTm[kr];
		cTot = cTotm[kr];

		if(memoTemperature == 1)
		{
			load >> Reactions->uKeq;
			load >> Reactions->k1 >> Reactions->k2;
			load >> Reactions->logFcent;
		}

		massFractionsInReactors.UseMatrixRowAsVector(kr,&ma);

		// Concentrations in the reactor
		wM = Reactions->GetMWFromMassFractions(ma);
		Reactions->GetMoleFractionsFromMassFractionsAndMW(xReactor,ma,wM);
		cReactor = cTot*xReactor;
		reactionsRateInReactors.UseMatrixRowAsVector(kr,&R);

		// Reaction Rates
		#if SYMBOLIC_KINETICS==1
			if (iAnalyticalJacobian == 0)
				Reactions->ComputeFromConcentrations_map(kr, cReactor, &R);
			else
				reactor[kr]->giveReactionRates(cTot, cReactor, R);
		#else
			Reactions->ComputeFromConcentrations_map(kr, cReactor, &R);
		#endif

	}
	if(memoTemperature == 1)
		load.End();
}

void OpenSMOKE_CSTRNetwork::GetReactionsRateInAllReactorsFromMassFractions
		(
			BzzVector &massFractionsInReactorsV, // reactors as rows
			BzzVector &reactionsRateInReactorsV  // reactors as rows
		)
{
	int kr;
	BzzLoad load;
	BzzVector ma;
	BzzVector R(numComponents);

	ChangeDimensions(numCSTRReactors*numComponents,&reactionsRateInReactorsV);

	if(memoTemperature == 1)
		load('*', "Temp/MemoTemperatureFunctions.tmp");

	for(kr = 1;kr <= numCSTRReactors;kr++)
	{
		T = temperature[kr];
		P = pressure[kr];
		logT = logTm[kr];
		loguRT = loguRTm[kr];
		cTot = cTotm[kr];

		if(memoTemperature == 1)
		{
			load >> Reactions->uKeq;
			load >> Reactions->k1 >> Reactions->k2;
			load >> Reactions->logFcent;
		}

		ma.UseSubVectorAsVector(numComponents, (kr - 1)*numComponents + 1, massFractionsInReactorsV);
		wM = Reactions->GetMWFromMassFractions(ma);
		Reactions->GetMoleFractionsFromMassFractionsAndMW(xReactor,ma,wM);
		cReactor = cTot*xReactor;
		R.UseSubVectorAsVector(numComponents,(kr - 1)*numComponents + 1,reactionsRateInReactorsV);
		 
		// Reaction Rates
		#if SYMBOLIC_KINETICS==1
			if (iAnalyticalJacobian == 0)
				Reactions->ComputeFromConcentrations_map(kr, cReactor, &R);
			else
				reactor[kr]->giveReactionRates(cTot, cReactor, R);
		#else
			Reactions->ComputeFromConcentrations_map(kr, cReactor, &R);
		#endif
		
	}

	if(memoTemperature == 1)	load.End();
}


void  OpenSMOKE_CSTRNetwork::GetDiagonalFactoredMatricesForLinearizedSistem
		(BzzMatrix &massFractionsInReactors)
{
	int kr;
	BzzLoad load;
	BzzMatrix dRC(numComponents,numComponents);
	BzzVector R(numComponents);
	BzzFactorizedGauss G;

	if(memoTemperature == 1)
		load('*', "Temp/MemoTemperatureFunctions.tmp");
	BzzSave save('*', "Temp/MemoJacobian.tmp");

	BzzVector ma;

	double start = BzzGetCpuTime();
	save << numCSTRReactors;

	for(kr = 1;kr <= numCSTRReactors;kr++)
	{
		T = temperature[kr];
		P = pressure[kr];
		logT = logTm[kr];
		loguRT = loguRTm[kr];
		cTot = cTotm[kr];

		if(memoTemperature == 1)
		{
			load >> Reactions->uKeq;
			load >> Reactions->k1 >> Reactions->k2;
			load >> Reactions->logFcent;
		}

		// Reaction Rates
		massFractionsInReactors.UseMatrixRowAsVector(kr,&ma);
		GetReactionRatesInSingleReactor(kr, ma, xReactor, cReactor, R);

		// Jacobian for the single reactor
		GetJacobianForSingleReactor(kr, cReactor, R, dRC);
		ReplaceBzzMatrixWithBzzFactorized(&dRC,&G);
		Factorize(&G);
		save << G;
	}

	if(memoTemperature == 1)	load.End();

	::BzzPrint("\nSeconds for CSTR: %e",BzzGetCpuTime() - start);
}

void  OpenSMOKE_CSTRNetwork::GetDiagonalMatricesForLinearizedSistem(string file,BzzVector &mfV)
{
	int kr;
	BzzLoad load;

	BzzMatrix dRC(numComponents,numComponents);
	BzzVector R(numComponents);
	BzzFactorizedGauss G;
	BzzVector ma;

	if(memoTemperature == 1)
		load('*', "Temp/MemoTemperatureFunctions.tmp");

	// The file called fileDiagonal in the BzzODE object is opened and the number of blocks
	// is written on this file
	BzzSave save('*',file);
	save << numCSTRReactors;

	double start = BzzGetCpuTime();

	for(kr = 1;kr <= numCSTRReactors;kr++)
	{
		T = temperature[kr];
		P = pressure[kr];

		logT = logTm[kr];
		loguRT = loguRTm[kr];
		cTot = cTotm[kr];
		if(memoTemperature == 1)
		{
			load >> Reactions->uKeq;
			load >> Reactions->k1 >> Reactions->k2;
			load >> Reactions->logFcent;
		}

		// Reaction Rates
		ma.UseSubVectorAsVector(numComponents, (kr - 1)*numComponents + 1, mfV);
		GetReactionRatesInSingleReactor(kr, ma, xReactor, cReactor, R);

		// Jacobian for the single reactor
		GetJacobianForSingleReactor(kr, cReactor, R, dRC);

		// Viene salvato su file il blocco diagonale della matrice del sistema, senza pero'
		// che ne sia stata fatta la fattorizzazione
		save << dRC;
	}
	if(memoTemperature == 1)	load.End();

	::BzzPrint("\nSeconds for CSTR: %e",BzzGetCpuTime() - start);
}


void OpenSMOKE_CSTRNetwork::SwapMassFractionsInReactors(BzzMatrix *massFractionsInReactors)
{
	ReorderByRows(massFractionsInReactors,cstrSequence);
}

void OpenSMOKE_CSTRNetwork::Save(string file) // formatted
{
	BzzSave save(file);
	save << numCSTRReactors;
	save << giveClusterIndexFromCellIndex;
	save << massFractionsInReactorsSolution;
	save.End();
}

void OpenSMOKE_CSTRNetwork::Save(char, string file)// binary
{
	BzzSave save('*',file);
	save << numCSTRReactors;
	save << giveClusterIndexFromCellIndex;
	save << massFractionsInReactorsSolution;
	save.End();
}

void OpenSMOKE_CSTRNetwork::SaveTemperature(string file) // formatted
{
	BzzSave save(file);
	save << temperature;
	save << qanm;
	save.End();
}

void OpenSMOKE_CSTRNetwork::SaveTemperature(char, string file)// binary
{
	BzzSave save('*',file);
	save << temperature;
	save << qanm;
	save.End();
}

void OpenSMOKE_CSTRNetwork::CorrectionCoefficient_None(	BzzVector &u_temperature, BzzVector &log_temperature,
													BzzMatrix &matrix_correction_k1,
													BzzMatrix &matrix_correction_uKeq,
													BzzMatrix &matrix_correction_k2)
{
	int		j, k, kCSTR;
	double	uT;
	double	uRT;

	for(kCSTR = 1;kCSTR <= numCSTRReactors; kCSTR++)
	{
		T = temperature[kCSTR];
		P = pressure[kCSTR];
		cTot = cTotm[kCSTR];
		uT = u_temperature[kCSTR];
		logT = log_temperature[kCSTR];
		uRT = uT * UR_CTOT;
		loguRT = log(uRT);

		for(k = 1;k <= Reactions->kinetics.numEquilibrium;k++)
		{
			j = Reactions->kinetics.reactionWithEquil[k];
			Reactions->uKeq[j] = exp(-Reactions->reactionDS_map[kCSTR][j] + Reactions->reactionDH_map[kCSTR][j] - loguRT * Reactions->kinetics.sumNuij[j]);

			matrix_correction_uKeq[kCSTR][j] = 1.00;
		}

		for(j = 1;j <= numReactions;j++)
		{
			Reactions->k1[j] = exp(Reactions->kinetics.k01[j] + Reactions->kinetics.beta1[j] * logT + Reactions->kinetics.E1[j] * uT);

			matrix_correction_k1[kCSTR][j] = 1.;

			if(Reactions->kinetics.jThirdBody[j] >= 2)
			{
				Reactions->k2[j] = exp(Reactions->kinetics.k02[j] + Reactions->kinetics.beta2[j] * logT + Reactions->kinetics.E2[j] * uT);
				
				matrix_correction_k2[kCSTR][j] = 1.00;
			}
		}
	}
}

void OpenSMOKE_CSTRNetwork::CorrectionCoefficient_SinExpansion(BzzVector &u_temperature, BzzVector &log_temperature,
														  BzzMatrix &matrix_correction_k1,
														  BzzMatrix &matrix_correction_uKeq,
														  BzzMatrix &matrix_correction_k2)
{
	int		j, k, kCSTR;
	double	uT;
	double	uRT;

	OpenSMOKE_SinIntegralDistribution SinIntegralDistribution;

	for(kCSTR = 1;kCSTR <= numCSTRReactors; kCSTR++)
	{
		double CoeffCorr;

		T = temperature[kCSTR];
		P = pressure[kCSTR];
		cTot = cTotm[kCSTR];
		uT = u_temperature[kCSTR];
		logT = log_temperature[kCSTR];
		uRT = uT * UR_CTOT;
		loguRT = log(uRT);

		double csiMean = (T-Tmin)/(Tmax[kCSTR]-Tmin);
		double g = 0.50*qanm[kCSTR]*T*T/DeltaT[kCSTR]/DeltaT[kCSTR];
		SinIntegralDistribution.Set(csiMean, g);

		for(k = 1;k <= Reactions->kinetics.numEquilibrium;k++)
		{
			j = Reactions->kinetics.reactionWithEquil[k];

			CoeffCorr = 1.;
			if (switchReactions[j] == 1)
			{
				Reactions->uKeq[j] = exp(-Reactions->reactionDS_map[kCSTR][j] + Reactions->reactionDH_map[kCSTR][j] - loguRT * Reactions->kinetics.sumNuij[j]);
				{
					double	EsuR = -(Reactions->kinetics.E1[j] + T*Reactions->reactionDH_map[kCSTR][j]);
					double	n = Reactions->kinetics.beta1[j] + Reactions->kinetics.sumNuij[j];

					CoeffCorr = SinIntegralDistribution.CorrectionCoefficient(n, EsuR, T);
					if (CoeffCorr > MaxCoeffCorr)	WarningLargeCorrectionCoefficient(kCSTR, T, g, j, Reactions->strReaction[j], EsuR, n, CoeffCorr);
					if (CoeffCorr < 0.)				WarningSmallCorrectionCoefficient(kCSTR, T, g, j, Reactions->strReaction[j], EsuR, n, CoeffCorr);
				}
			}

			matrix_correction_uKeq[kCSTR][j] = CoeffCorr;
		}

		for(j = 1;j <= numReactions;j++)
		{
			CoeffCorr = 1.;
			if (switchReactions[j] == 1)
			{
				Reactions->k1[j] = exp(Reactions->kinetics.k01[j] + Reactions->kinetics.beta1[j] * logT + Reactions->kinetics.E1[j] * uT);
				{
					double	EsuR = -Reactions->kinetics.E1[j];
					double	n = Reactions->kinetics.beta1[j];

					CoeffCorr = SinIntegralDistribution.CorrectionCoefficient(n, EsuR, T);
					if (CoeffCorr > MaxCoeffCorr)	WarningLargeCorrectionCoefficient(kCSTR, T, g, j, Reactions->strReaction[j], EsuR, n, CoeffCorr);
					if (CoeffCorr < 0.)				WarningSmallCorrectionCoefficient(kCSTR, T, g, j, Reactions->strReaction[j], EsuR, n, CoeffCorr);
				}
			}
			matrix_correction_k1[kCSTR][j] = CoeffCorr;

			if(Reactions->kinetics.jThirdBody[j] >= 2)
			{
				CoeffCorr = 1.;
				if (switchReactions[j] == 1)
				{
					Reactions->k2[j] = exp(Reactions->kinetics.k02[j] + Reactions->kinetics.beta2[j] * logT + Reactions->kinetics.E2[j] * uT);
					{
						double	EsuR = -Reactions->kinetics.E2[j];
						double	n = Reactions->kinetics.beta2[j];

						CoeffCorr = SinIntegralDistribution.CorrectionCoefficient(n, EsuR, T);
						
						if (CoeffCorr > MaxCoeffCorr)	WarningLargeCorrectionCoefficient(kCSTR, T, g, j, Reactions->strReaction[j], EsuR, n, CoeffCorr);
						if (CoeffCorr < 0.)				WarningSmallCorrectionCoefficient(kCSTR, T, g, j, Reactions->strReaction[j], EsuR, n, CoeffCorr);
					}
				}
				matrix_correction_k2[kCSTR][j] = CoeffCorr;
			}
		}

		{
			Tk_20000[kCSTR] = KineticTemperature(SinIntegralDistribution.CorrectionCoefficient(0., 20000./1.987, T), T, 20000./1.987, 0.);
			Tk_40000[kCSTR] = KineticTemperature(SinIntegralDistribution.CorrectionCoefficient(0., 40000./1.987, T), T, 40000./1.987, 0.);
			Tk_60000[kCSTR] = KineticTemperature(SinIntegralDistribution.CorrectionCoefficient(0., 60000./1.987, T), T, 60000./1.987, 0.);
		}

	}	// End Cycle On each Reactor
}

void OpenSMOKE_CSTRNetwork::CorrectionCoefficient_DeltaDirac(	BzzVector &u_temperature, BzzVector &log_temperature,
																BzzMatrix &matrix_correction_k1,
																BzzMatrix &matrix_correction_uKeq,
																BzzMatrix &matrix_correction_k2)
{
	int		j, k, kCSTR;
	double	uT;
	double	uRT;

	OpenSMOKE_DoubleDeltaDiracDistribution DoubleDeltaDirac;

	for(kCSTR = 1;kCSTR <= numCSTRReactors;kCSTR++)
	{
		double CoeffCorr;

		T = temperature[kCSTR];
		P = pressure[kCSTR];
		cTot = cTotm[kCSTR];
		uT = u_temperature[kCSTR];
		logT = log_temperature[kCSTR];
		uRT = uT * UR_CTOT;
		loguRT = log(uRT);


		double csiMean = (T-Tmin)/(Tmax[kCSTR]-Tmin);
		double g = 0.50*qanm[kCSTR]*T*T/DeltaT[kCSTR]/DeltaT[kCSTR];
		DoubleDeltaDirac.Set(csiMean, g);

		for(k = 1;k <= Reactions->kinetics.numEquilibrium;k++)
		{
			j = Reactions->kinetics.reactionWithEquil[k];

			CoeffCorr = 1.;
			if (switchReactions[j] == 1)
			{		
				Reactions->uKeq[j] = exp(-Reactions->reactionDS_map[kCSTR][j] + Reactions->reactionDH_map[kCSTR][j] - loguRT * Reactions->kinetics.sumNuij[j]);
				{
					double	EsuR = -(Reactions->kinetics.E1[j] + T*Reactions->reactionDH_map[kCSTR][j]);
					double	n = Reactions->kinetics.beta1[j] + Reactions->kinetics.sumNuij[j];

					CoeffCorr = DoubleDeltaDirac.ReactionCorrectionCoefficient(n, EsuR, T, Tmin, Tmax[kCSTR]);
					if (CoeffCorr > MaxCoeffCorr)	WarningLargeCorrectionCoefficient(kCSTR, T, g, j, Reactions->strReaction[j], EsuR, n, CoeffCorr);
					if (CoeffCorr < 0.)				WarningSmallCorrectionCoefficient(kCSTR, T, g, j, Reactions->strReaction[j], EsuR, n, CoeffCorr);
				}
			}
			matrix_correction_uKeq[kCSTR][j] = CoeffCorr;
		}

		for(j = 1;j <= numReactions;j++)
		{
			CoeffCorr = 1.;
			if (switchReactions[j] == 1)
			{
				Reactions->k1[j] = exp(Reactions->kinetics.k01[j] + Reactions->kinetics.beta1[j] * logT + Reactions->kinetics.E1[j] * uT);
				{
					double	EsuR = -Reactions->kinetics.E1[j];
					double	n = Reactions->kinetics.beta1[j];

					CoeffCorr = DoubleDeltaDirac.ReactionCorrectionCoefficient(n, EsuR, T, Tmin, Tmax[kCSTR]);
					if (CoeffCorr > MaxCoeffCorr)	WarningLargeCorrectionCoefficient(kCSTR, T, g, j, Reactions->strReaction[j], EsuR, n, CoeffCorr);
					if (CoeffCorr < 0.)				WarningSmallCorrectionCoefficient(kCSTR, T, g, j, Reactions->strReaction[j], EsuR, n, CoeffCorr);
				}
			}
			matrix_correction_k1[kCSTR][j] = CoeffCorr;

			if(Reactions->kinetics.jThirdBody[j] >= 2)
			{
				CoeffCorr = 1.;
				if (switchReactions[j] == 1)
				{
					Reactions->k2[j] = exp(Reactions->kinetics.k02[j] + Reactions->kinetics.beta2[j] * logT + Reactions->kinetics.E2[j] * uT);
					{
						double	EsuR = -Reactions->kinetics.E2[j];
						double	n = Reactions->kinetics.beta2[j];

						CoeffCorr = DoubleDeltaDirac.ReactionCorrectionCoefficient(n, EsuR, T, Tmin, Tmax[kCSTR]);
						if (CoeffCorr > MaxCoeffCorr)	WarningLargeCorrectionCoefficient(kCSTR, T, g, j, Reactions->strReaction[j], EsuR, n, CoeffCorr);
						if (CoeffCorr < 0.)				WarningSmallCorrectionCoefficient(kCSTR, T, g, j, Reactions->strReaction[j], EsuR, n, CoeffCorr);
					}
				}
				matrix_correction_k2[kCSTR][j] = CoeffCorr;
			}
		}

		{
			Tk_20000[kCSTR] = KineticTemperature(DoubleDeltaDirac.ReactionCorrectionCoefficient(0., 20000./1.987, T, Tmin, Tmax[kCSTR]), T, 20000./1.987, 0.);
			Tk_40000[kCSTR] = KineticTemperature(DoubleDeltaDirac.ReactionCorrectionCoefficient(0., 40000./1.987, T, Tmin, Tmax[kCSTR]), T, 40000./1.987, 0.);
			Tk_60000[kCSTR] = KineticTemperature(DoubleDeltaDirac.ReactionCorrectionCoefficient(0., 60000./1.987, T, Tmin, Tmax[kCSTR]), T, 60000./1.987, 0.);
		}
		
		/*BzzVector veq = matrix_correction_uKeq.GetRow(kCSTR);
		for(int j=1;j<=veq.Size();j++)
			cout << T << " " << j << " " << veq[j] << " " << Tmax[kCSTR] << " " << qanm[kCSTR] << " " << DeltaT[kCSTR] << " " << csiMean << " " << g << " " << Tmin << endl;
		getchar();*/
	}	// End Cycle On each Reactor
}


void OpenSMOKE_CSTRNetwork::CorrectionCoefficient_BetaPDF(	BzzVector &u_temperature, BzzVector &log_temperature,
															BzzMatrix &matrix_correction_k1, BzzMatrix &matrix_correction_uKeq, BzzMatrix &matrix_correction_k2)
{
	int		j, k, kCSTR;
	double	uT;
	double	uRT;

	OpenSMOKE_BetaDistribution BetaDistribution;

	for(kCSTR = 1;kCSTR <= numCSTRReactors;kCSTR++)
	{
		double CoeffCorr;

		T = temperature[kCSTR];
		P = pressure[kCSTR];
		cTot = cTotm[kCSTR];
		uT = u_temperature[kCSTR];
		logT = log_temperature[kCSTR];
		uRT = uT * UR_CTOT;
		loguRT = log(uRT);

		// Beta PDF parameters
		double csiMean = (T-Tmin)/(Tmax[kCSTR]-Tmin);
		double g = 0.50*qanm[kCSTR]*T*T/DeltaT[kCSTR]/DeltaT[kCSTR];
		double alfa = csiMean*(csiMean*(1.-csiMean)/g-1.);
		double beta = (1.-csiMean)*(csiMean*(1.-csiMean)/g-1.);

		BetaDistribution.Set(alfa, beta);
		for(k = 1;k <= Reactions->kinetics.numEquilibrium;k++)
		{
			j = Reactions->kinetics.reactionWithEquil[k];
			CoeffCorr = 1.;
			if (switchReactions[j] == 1)
			{
				Reactions->uKeq[j] = exp(-Reactions->reactionDS_map[kCSTR][j] + Reactions->reactionDH_map[kCSTR][j] - loguRT * Reactions->kinetics.sumNuij[j]);
				{
					double	EsuR = -(Reactions->kinetics.E1[j] + T*Reactions->reactionDH_map[kCSTR][j]);
					double	n = Reactions->kinetics.beta1[j] + Reactions->kinetics.sumNuij[j];

					CoeffCorr = BetaDistribution.ReactionCorrectionCoefficient(n, EsuR, T, Tmin, Tmax[kCSTR]);
					if (CoeffCorr > MaxCoeffCorr)	WarningLargeCorrectionCoefficient(kCSTR, T, g, j, Reactions->strReaction[j], EsuR, n, CoeffCorr);
					if (CoeffCorr < 0.)				WarningSmallCorrectionCoefficient(kCSTR, T, g, j, Reactions->strReaction[j], EsuR, n, CoeffCorr);
				}
			}
			matrix_correction_uKeq[kCSTR][j] = CoeffCorr;
		}

		for(j = 1;j <= numReactions;j++)
		{
			CoeffCorr = 1.;
			if (switchReactions[j] == 1)
			{
				Reactions->k1[j] = exp(Reactions->kinetics.k01[j] + Reactions->kinetics.beta1[j] * logT + Reactions->kinetics.E1[j] * uT);
				{
					double	EsuR = -Reactions->kinetics.E1[j];
					double	n = Reactions->kinetics.beta1[j];

					CoeffCorr = BetaDistribution.ReactionCorrectionCoefficient(n, EsuR, T, Tmin, Tmax[kCSTR]);
					if (CoeffCorr > MaxCoeffCorr)	WarningLargeCorrectionCoefficient(kCSTR, T, g, j, Reactions->strReaction[j], EsuR, n, CoeffCorr);
					if (CoeffCorr < 0.)				WarningSmallCorrectionCoefficient(kCSTR, T, g, j, Reactions->strReaction[j], EsuR, n, CoeffCorr);
				}
			}	
			matrix_correction_k1[kCSTR][j] = CoeffCorr;

			if(Reactions->kinetics.jThirdBody[j] >= 2)
			{
				CoeffCorr = 1.;
				if (switchReactions[j] == 1)
				{
					Reactions->k2[j] = exp(Reactions->kinetics.k02[j] + Reactions->kinetics.beta2[j] * logT + Reactions->kinetics.E2[j] * uT);
					{
						double	EsuR = -Reactions->kinetics.E2[j];
						double	n = Reactions->kinetics.beta2[j];

						CoeffCorr = BetaDistribution.ReactionCorrectionCoefficient(n, EsuR, T, Tmin, Tmax[kCSTR]);
						if (CoeffCorr > MaxCoeffCorr)	WarningLargeCorrectionCoefficient(kCSTR, T, g, j, Reactions->strReaction[j], EsuR, n, CoeffCorr);
						if (CoeffCorr < 0.)				WarningSmallCorrectionCoefficient(kCSTR, T, g, j, Reactions->strReaction[j], EsuR, n, CoeffCorr);
					}				
				}
				matrix_correction_k2[kCSTR][j] = CoeffCorr;
			}
		}
		
		{
			Tk_20000[kCSTR] = KineticTemperature(BetaDistribution.ReactionCorrectionCoefficient(0., 20000./1.987, T, Tmin, Tmax[kCSTR]), T, 20000./1.987, 0.);
			Tk_40000[kCSTR] = KineticTemperature(BetaDistribution.ReactionCorrectionCoefficient(0., 40000./1.987, T, Tmin, Tmax[kCSTR]), T, 40000./1.987, 0.);
			Tk_60000[kCSTR] = KineticTemperature(BetaDistribution.ReactionCorrectionCoefficient(0., 60000./1.987, T, Tmin, Tmax[kCSTR]), T, 60000./1.987, 0.);
		}
	}	// End Cycle On each Reactor
}

void OpenSMOKE_CSTRNetwork::CorrectionCoefficient_ClippedGaussianPDF(	BzzVector &u_temperature, 
																		BzzVector &log_temperature, BzzMatrix &matrix_correction_k1,
																		BzzMatrix &matrix_correction_uKeq, BzzMatrix &matrix_correction_k2)
{
	int		j, k, kCSTR;
	double	uT;
	double	uRT;

	OpenSMOKE_ClippedGaussianAccurateDistribution ClippedGaussianDistribution;

	for(kCSTR = 1;kCSTR <= numCSTRReactors;kCSTR++)
	{
		double CoeffCorr;
		
		T = temperature[kCSTR];
		P = pressure[kCSTR];
		cTot = cTotm[kCSTR];
		uT = u_temperature[kCSTR];
		logT = log_temperature[kCSTR];
		uRT = uT * UR_CTOT;
		loguRT = log(uRT);

		double csiMean = (T-Tmin)/(Tmax[kCSTR]-Tmin);
		double g = 0.50*qanm[kCSTR]*T*T/DeltaT[kCSTR]/DeltaT[kCSTR];

		ClippedGaussianDistribution.Set(csiMean, g);
		
		for(k = 1;k <= Reactions->kinetics.numEquilibrium;k++)
		{
			j = Reactions->kinetics.reactionWithEquil[k];
			CoeffCorr = 1.;
			if (switchReactions[j] == 1)
			{
				Reactions->uKeq[j] = exp(-Reactions->reactionDS_map[kCSTR][j] + Reactions->reactionDH_map[kCSTR][j] - loguRT * Reactions->kinetics.sumNuij[j]);
				{
					double	EsuR = -(Reactions->kinetics.E1[j] + T*Reactions->reactionDH_map[kCSTR][j]);
					double	n = Reactions->kinetics.beta1[j] + Reactions->kinetics.sumNuij[j];

					CoeffCorr = ClippedGaussianDistribution.ReactionCorrectionCoefficient(n, EsuR, T, Tmin, Tmax[kCSTR]);
					if (CoeffCorr > MaxCoeffCorr)	WarningLargeCorrectionCoefficient(kCSTR, T, g, j, Reactions->strReaction[j], EsuR, n, CoeffCorr);
					if (CoeffCorr < 0.)				WarningSmallCorrectionCoefficient(kCSTR, T, g, j, Reactions->strReaction[j], EsuR, n, CoeffCorr);
				}	
			}
			matrix_correction_uKeq[kCSTR][j] = CoeffCorr;
		}

		for(j = 1;j <= numReactions;j++)
		{
			CoeffCorr = 1.;
			if (switchReactions[j] == 1)
			{
				Reactions->k1[j] = exp(Reactions->kinetics.k01[j] + Reactions->kinetics.beta1[j] * logT + Reactions->kinetics.E1[j] * uT);
				{
					double	EsuR = -Reactions->kinetics.E1[j];
					double	n = Reactions->kinetics.beta1[j];

					CoeffCorr = ClippedGaussianDistribution.ReactionCorrectionCoefficient(n, EsuR, T, Tmin, Tmax[kCSTR]);
					if (CoeffCorr > MaxCoeffCorr)	WarningLargeCorrectionCoefficient(kCSTR, T, g, j, Reactions->strReaction[j], EsuR, n, CoeffCorr);
					if (CoeffCorr < 0.)				WarningSmallCorrectionCoefficient(kCSTR, T, g, j, Reactions->strReaction[j], EsuR, n, CoeffCorr);
				}
			}
			matrix_correction_k1[kCSTR][j] = CoeffCorr;
						
			if(Reactions->kinetics.jThirdBody[j] >= 2)
			{
				CoeffCorr = 1.;
				if (switchReactions[j] == 1)
				{
					Reactions->k2[j] = exp(Reactions->kinetics.k02[j] + Reactions->kinetics.beta2[j] * logT + Reactions->kinetics.E2[j] * uT);
					{
						double	EsuR = -Reactions->kinetics.E2[j];
						double	n = Reactions->kinetics.beta2[j];

						CoeffCorr = ClippedGaussianDistribution.ReactionCorrectionCoefficient(n, EsuR, T, Tmin, Tmax[kCSTR]);
						if (CoeffCorr > MaxCoeffCorr)	WarningLargeCorrectionCoefficient(kCSTR, T, g, j, Reactions->strReaction[j], EsuR, n, CoeffCorr);
						if (CoeffCorr < 0.)				WarningSmallCorrectionCoefficient(kCSTR, T, g, j, Reactions->strReaction[j], EsuR, n, CoeffCorr);
					}	
				}
				matrix_correction_k2[kCSTR][j] = CoeffCorr;
			}
		}
	
		{
		
			bool iPrintReactor = false;
			if (kCSTR%(int(0.25*numCSTRReactors)+1) == 0)	iPrintReactor = true;
			if(iPrintReactor==true)
			{
				cout << " -- R#:     " << kCSTR	<< "/"		<< numCSTRReactors << endl;
				cout << "    csi:    " << csiMean << endl;
				cout << "    g/gmax: " << g/(csiMean*(1.-csiMean)) << endl;
			}

			Tk_20000[kCSTR] = KineticTemperature(ClippedGaussianDistribution.ReactionCorrectionCoefficient(0., 20000./1.987, T, Tmin, Tmax[kCSTR]), T, 20000./1.987, 0.);
			Tk_40000[kCSTR] = KineticTemperature(ClippedGaussianDistribution.ReactionCorrectionCoefficient(0., 40000./1.987, T, Tmin, Tmax[kCSTR]), T, 40000./1.987, 0.);
			Tk_60000[kCSTR] = KineticTemperature(ClippedGaussianDistribution.ReactionCorrectionCoefficient(0., 60000./1.987, T, Tmin, Tmax[kCSTR]), T, 60000./1.987, 0.);
		}

	}	// End Cycle On each Reactor
}

void OpenSMOKE_CSTRNetwork::ConnectionMatrix()
{
	int i;

	BzzVectorInt *connection;
	connection = new BzzVectorInt[numCSTRReactors+1];

	// Reattori in comunicazione con il reattore i
	for(i = 1;i <= numCSTRReactors;i++)
	{
		connection[i].Append(i);
		for(int j=1;j<=FromClusterToCluster_Index[i].Size();j++)
		{
			if (FromClusterToCluster_Index[i][j] != 0)
			{
				connection[i].Append(FromClusterToCluster_Index[i][j]);
				connection[FromClusterToCluster_Index[i][j]].Append(i);
			}
		}
	}

	ofstream fOutput;
	openOutputFileAndControl(fOutput, "Temp/ConnectionMatrix.out");

	for(i = 1;i <= numCSTRReactors;i++)
	{
		for(int j=1;j<=connection[i].Size();j++)
			fOutput << connection[i][j] << "\t";
		fOutput << endl;
	}
	fOutput.close();

	openOutputFileAndControl(fOutput, "Temp/ConnectionMatrix_MATLAB.out");
	for(i = 1;i <= numCSTRReactors;i++)
		for(int j=1;j<=connection[i].Size();j++)
			fOutput << i << "\t" << connection[i][j] << "\t" << "1" << endl;
	fOutput.close();
}

void OpenSMOKE_CSTRNetwork::WarningLargeCorrectionCoefficient(const int k, const double T, const double g, const int iReaction, 
							const string stringReaction, const double EsuR, const double n, double &CoeffCorr)
{
	double csi = (T-Tmin)/DeltaT[k];
	double Tk = KineticTemperature(CoeffCorr, T, EsuR, n);

	if (CoeffCorr >= 1.e48)
	{
		fWarning << "Correction coefficient too large..." << endl;
		fWarning << "  Cluster:      " << k << endl;
		fWarning << "  Temperature:  " << T << " K" << endl;
		fWarning << "  g/gmax:       " << g/(1.-csi)/csi << endl;
		fWarning << "  Reaction:     " << iReaction	<<  "\t" << stringReaction << endl;
		fWarning << "  E:            " << EsuR*1.987 << " cal/mol"	<< endl;
		fWarning << "  n:            " << n			<< endl;
		fWarning << "  Correction:   " << CoeffCorr	<< endl;
		fWarning << "  Tk:           " << Tk << " K"<< endl;
		fWarning << "  Tk/T:         " << Tk/T << endl;
	}
	CoeffCorr = MaxCoeffCorr;

	if (Tk > Tmax[k])
	{
		fWarning << "Kinetic equivalent temperature too large..." << endl;
		fWarning << "  Cluster:      " << k << endl;
		fWarning << "  Temperature:  " << T << " K" << endl;
		fWarning << "  g/gmax:       " << g/(1.-csi)/csi << endl;
		fWarning << "  Reaction:     " << iReaction	<<  "\t" << stringReaction << endl;
		fWarning << "  E:            " << EsuR*1.987 << " cal/mol"	<< endl;
		fWarning << "  n:            " << n			<< endl;
		fWarning << "  Correction:   " << CoeffCorr	<< endl;
		fWarning << "  Tk:           " << Tk << " K"<< endl;
		fWarning << "  Tmax:         " << Tmax[k] << " K"<< endl;
		fWarning << "  Tk/T:         " << Tk/T << endl;
	}
}

void OpenSMOKE_CSTRNetwork::WarningSmallCorrectionCoefficient(const int k, const double T, const double g, const int iReaction, 
							const string stringReaction, const double EsuR, const double n, double &CoeffCorr)
{
	double csi = (T-Tmin)/DeltaT[k];
	fWarning << "Negative correction coefficient..." << endl;
	fWarning << "  Cluster:      " << k << endl;
	fWarning << "  Temperature:  " << T << " K" << endl;
	fWarning << "  g/gmax:       " << g/(1.-csi)/csi << endl;
	fWarning << "  Reaction:     " << iReaction	<<  "\t" << stringReaction << endl;
	fWarning << "  E:            " << EsuR*1.987 << " cal/mol"	<< endl;
	fWarning << "  n:            " << n			<< endl;
	fWarning << "  Correction:   " << CoeffCorr	<< endl;

	CoeffCorr = 1.;

	ErrorMessage("Negative correction coefficient...\nSee Warning.log file...");
}

// La rete di reattori viene risolta GLOBALMENTE attraverso un sistema ODE
void OpenSMOKE_CSTRNetwork::GetFourth(const double tEnd)
{
	double F0,F1,maxF0,maxF1;
	double uvariables;
	BzzVectorInt blockDimensions;
	BzzVector initialValues;
	BzzVector residual_Vector;

	// Memory Allocation
	ChangeDimensions(numCSTRReactors, &blockDimensions);
	ChangeDimensions(numCSTRReactors * numComponents, &residual_Vector);
	ChangeDimensions(numCSTRReactors * numComponents, &dummyVector);

	// Initializzazioni
	blockDimensions = numComponents;
	uvariables = 1. / (numComponents * numCSTRReactors);

	// Initializzazioni
	initialValues = massFractionsInReactorsSolution;
	volumeMolecularWeightT_Vector = volumeMolecularWeightT;
	feed_Vector = feed;

	// Starting Computations
	double start = BzzGetCpuTime();

	// Restituisce i residui per tutte le equazioni in tutti i reattori
	GetAllResiduals(initialValues, residual_Vector);

	// Calcola il residuo medio e massimo
	F0 = residual_Vector.GetSumAbsElements() * uvariables;
	maxF0 = residual_Vector.MaxAbs();
	cout << "Solving Linear System in Fourth Method (ODE)" << endl;
	cout << "  ** Mean Residual Before Fourth Method" << F0 << endl;
	cout << "  ** Max  Residual Before Fourth Method" << maxF0 << endl;

	// Solving the ODE System
	{
		double tStart;
		double stopODE;
		BzzVector yMin(initialValues.Size());
		BzzVector yMax(initialValues.Size());
		OpenSMOKE_CSTRNetwork_MyOdeSystemObjectAllCSTR obj(this);

		tStart = 0.;
		stopODE = tolAbs * uvariables;
		BzzMatrixSparseLockedByRows Q = Ld;
		//obj(this);

        string fileSparseMatrix = "Temp/A.tmp";
        char fileSM[200];
        strcpy(fileSM, fileSparseMatrix.c_str());
		BzzOdeSparseStiffObject o(initialValues, tStart, &obj, blockDimensions, &Q, fileSM);
		yMin = 0.;
		yMax = 1.;
		o.SetMinimumConstraints(yMin);
		o.SetMaximumConstraints(yMax);
		
		// Default values: (A) 1e-15      (R) 100*MachEps()
		//o.SetTollRel(1.e2*MachEpsFloat());
		o.SetTollAbs(1.e-15);

		o.StopIntegrationBeforeRecalcuatingJacobian(30);
		o.StopIntegrationWhenSumAbsY1IsLessThan(stopODE);

		if		(F0 > 1.e-6)	massFractionsInReactorsSolution_Vector = o(tEnd);
		else					massFractionsInReactorsSolution_Vector = o(1.);

		cout	<< "CPU Time for Solving Linear System in Fourth Method (ODE): "
				<< BzzGetCpuTime() - start << endl;
	}

	// Calcolo dei nuovi residui dopo l'integrazione del sistema ODE
	GetAllResiduals(massFractionsInReactorsSolution_Vector, residual_Vector);
	F1 = residual_Vector.GetSumAbsElements() * uvariables;
	maxF1 = residual_Vector.MaxAbs();
	cout << "  ** Mean Residual After Fourth Method (ODE): " << F1 << endl;
	cout << "  ** Max  Residual After Fourth Method (ODE): " << maxF1 << endl;

	// Accetta la nuova soluzione che sicuramente e' buona in quanto ottenuta
	// attraverso la soluzione di un transitorio reale
	CopyDataFromVector(&massFractionsInReactorsSolution,massFractionsInReactorsSolution_Vector);
}
