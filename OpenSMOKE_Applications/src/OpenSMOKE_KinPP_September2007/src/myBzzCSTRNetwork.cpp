#include "BzzMath.hpp"
#include "myBzzCSTRNetwork.hpp"
#include "MyOdeSystemObjectOneCSTR.h"
#include "MyOdeSystemObjectAllCSTR.h"


const double MAX_SUM_F = 1.e-10;

const double myBzzCSTRNetwork::R_CTOT = 0.0820578337034;
const double myBzzCSTRNetwork::UR_CTOT = 1./R_CTOT;
const char *const myBzzCSTRNetwork::BZZ_ERROR = "myBzzCSTRNetwork Class\n";


int myBzzCSTRNetwork::count = 0;
int myBzzCSTRNetwork::countInScope = 0;


void myBzzCSTRNetwork::SetSolution(BzzMatrix &omega)
{
	massFractionsInReactorsSolution = omega;
}

int myBzzCSTRNetwork::WhoAmI(void) const 
{
	return whoAmI;
}

void myBzzCSTRNetwork::SetTolRel(double tolr)
{
	tolRel = tolr;
}
	
void myBzzCSTRNetwork::SetTolAbs(double tola)
{
	tolAbs = tola;
}
	
void myBzzCSTRNetwork::SetMemoTemperature(void)
{
	memoTemperature = 0;	// 1 salva su file 0 memorizza
} 
	
int myBzzCSTRNetwork::GetNumComponents(void)
{
	return Reactions.NumberOfSpecies();
}
	
int myBzzCSTRNetwork::GetNumReactions(void)
{
	return Reactions.NumberOfReactions();
}
	
int myBzzCSTRNetwork::GetNumCSTRReactors(void)
{
	return numCSTRReactors;
}

void myBzzCSTRNetwork::SetTasksPrint(void)
{
	printTasks = 1;
}

void myBzzCSTRNetwork::SetSubTasksPrint(void)
{
	printSubTasks = 1;
}

myBzzCSTRNetwork::~myBzzCSTRNetwork(void)
{
	if(Reactions.NumberOfReactions() == 0 || Reactions.NumberOfSpecies() == 0) 
		return;
}

myBzzCSTRNetwork::myBzzCSTRNetwork(void)
{
	tolRel = 1.e-2;
	tolAbs = 1.e-8;
	count++;
	countInScope++;
	whoAmI = count;
	numCSTRReactors = 1;
}

void myBzzCSTRNetwork::readKineticSheme()
{
//	Reactions.setOptions(1, "ONLY_THERMODYNAMICS");
//	Reactions.setupGlobal("Input/KineticScheme", 10);

	Reactions.SetupBinary("Input/KineticScheme");

	numReactions = Reactions.NumberOfReactions();
	numComponents = Reactions.NumberOfSpecies();
}

void openInputFile(ifstream &file, char* fileName)
{
	file.open(fileName, ios::in);
	if(!file.good())
	{
		cerr << "FATAL ERROR: The file " << fileName << " cannot be opened!" << endl;
		cerr << endl;
		exit(-1);
	}
}

void openOutputFile(ofstream &file, char* fileName)
{
	file.open(fileName, ios::out);
	if(!file.good())
	{
		cerr << "FATAL ERROR: The file " << fileName << " cannot be opened!" << endl;
		cerr << endl;
		exit(-1);
	}
}

void checkInputfile(ifstream &file, char* fileName)
{
	if (file.fail())
	{
		cerr << "FATAL ERROR: It was impossible to correctly read some values in the file " << fileName << endl;
		cerr << endl;
		exit(-1);
	}
}

void myBzzCSTRNetwork::reading_the_network_file(char *fileNetwork, int &countInput)
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
	openInputFile(fInput, fileNetwork);

	// Reading the number of original CFD cells
	fInput >> numCSTRReactors;
	checkInputfile(fInput, fileNetwork);
	
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
		fInput >> indexCell;
		fInput >> temperatureCell;
		fInput >> pressureCell;
		fInput >> volumeCell;
		fInput >> TvarianceCell;
		checkInputfile(fInput, fileNetwork);		
		
		// The cells must be ordered
		if(k != indexCell)	BzzError("Wrong CSTR sequence");
			
		// Checking if the cell numbers are correct
		if(indexCell < 1 || indexCell > numCSTRReactors)
			BzzError("%s%sReactor out of range", myBzzCSTRNetwork::BZZ_ERROR,fileNetwork);

		temperature[indexCell] = temperatureCell;
		pressure[indexCell] = pressureCell;
		volume[indexCell] = volumeCell;
		qanm[indexCell] = TvarianceCell;
			
		// Number of external input in the CFD cell
		// ------------------------------------------------------------------------------
		fInput >> nInput;
		checkInputfile(fInput, fileNetwork);	
		if(nInput != 0)
		{
			countInput++;

			for(j=1;j<=nInput;j++)
			{
				fInput >> inputMassFlowRate;
				checkInputfile(fInput, fileNetwork);	
					
				fInput >> numberOfSpeciesCFD;
				checkInputfile(fInput, fileNetwork);	
					
				for(i = 1;i <= numberOfSpeciesCFD;i++)
				{				
					fInput >> indexSpecies;
					fInput >> omegaSpecies;
					checkInputfile(fInput, fileNetwork);	
				}
			}
		}
	
		// Number of output (internal) in the CFD cell
		// ------------------------------------------------------------------------------
		fInput >> nOutput;
		checkInputfile(fInput, fileNetwork);	
		for(i = 1;i <= nOutput;i++)
		{
			fInput >> destination_cell;
			fInput >> outputMassFlowRate;
			fInput >> outputDiffusionRate;
			checkInputfile(fInput, fileNetwork);

			// Connection number of each CFD cell
			if(destination_cell > 0)
			{
				countCSTRConnect[indexCell]++;
				countCSTRConnect[destination_cell]++;
			}
		}
	}
	fInput.close();

	// Memory allocation for the sparse matrix of connections
	cstrConnect(countCSTRConnect);

	// Saving information on CFD Network
	BzzSave saveNetworkInput("Temp/CFDNetworkInput.bzz");
	saveNetworkInput << countInput;

	// Opening the network file
	openInputFile(fInput, fileNetwork);

	// Reading the number of original CFD cells
	fInput >> numCSTRReactors;
	checkInputfile(fInput, fileNetwork);
		
	
	// Cycle on every original CFD cell
	for(k = 1;k <= numCSTRReactors;k++)
	{		
		fInput >> indexCell;
		fInput >> temperatureCell;
		fInput >> pressureCell;
		fInput >> volumeCell;
		fInput >> TvarianceCell;
		checkInputfile(fInput, fileNetwork);	
		
		// Number of external input in the CFD cell
		// ------------------------------------------------------------------------------
		fInput >> nInput;
		checkInputfile(fInput, fileNetwork);	
		if(nInput != 0)
		{
			double totalInputMassFlowRate = 0.;
			BzzVector speciesInletFeed(numComponents);
			
			saveNetworkInput << indexCell;
			
			for(j=1; j<=nInput; j++)
			{
				fInput >> inputMassFlowRate;
				checkInputfile(fInput, fileNetwork);		
			
				fInput >> numberOfSpeciesCFD;
				checkInputfile(fInput, fileNetwork);	
			
				totalInputMassFlowRate += inputMassFlowRate;
				if(j == nInput)
					saveNetworkInput << totalInputMassFlowRate << numberOfSpeciesCFD; 
					
				double omegaSum = 0.;
				for(i = 1;i <= numberOfSpeciesCFD;i++)
				{
					fInput >> indexSpecies;
					fInput >> omegaSpecies;
					checkInputfile(fInput, fileNetwork);

					omegaSum += omegaSpecies;
					speciesInletFeed[indexSpecies] += inputMassFlowRate * omegaSpecies;
					
					if(j == nInput)
						saveNetworkInput << indexSpecies << speciesInletFeed[indexSpecies];
				}
					
				if(omegaSum < .9999999 || omegaSum > 1.0000001)
					BzzError("Wrong fraction in %d rector: %22.14e", indexCell, omegaSum);
			}

			massExternalInputCFDCells[indexCell] = totalInputMassFlowRate;
		}

		// Number of output (internal) in the CFD cell
		// ------------------------------------------------------------------------------
		fInput >> nOutput;
		checkInputfile(fInput, fileNetwork);
		for(i = 1;i <= nOutput;i++)
		{
			fInput >> destination_cell;
			fInput >> outputMassFlowRate;
			fInput >> outputDiffusionRate;
			checkInputfile(fInput, fileNetwork);

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
	openOutputFile(fOutput, "Temp/MassFluxesCFDNetwork.out");
	for(k = 1;k <= numCSTRReactors;k++)
	{
		double mIN  = massExternalInputCFDCells[k] + massInputCFDCells[k];
		double mOUT = massExternalOutputCFDCells[k] + massOutputCFDCells[k];
		double relative_error = fabs((mOUT-mIN)/(mIN+mOUT));
		meanError+=relative_error;
		if (relative_error>maxError)
			maxError = relative_error;

		if (relative_error > 5.e-3)
		{
			cerr << "FATAL ERROR: Wrong mass balance on CFD cell " << k << endl;
			cerr << endl;
			exit(-1);
		}

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
	cout << " ** Mean relative error in mass balances on CFD cells: " << meanError << endl;
	cout << " ** Max  relative error in mass balances on CFD cells: " << maxError << endl;

}

void myBzzCSTRNetwork::read_tolerances(int numCSTRReactors, int cicloCluster)
{
	int i;

	BzzLoad loadTol("Input/Tolerances.bzz");
	loadTol >> tolRelCBase >> tolAbsCBase >> numTolT;

	ChangeDimensions(numTolT, &tInf);
	ChangeDimensions(numTolT, &tSup);
	ChangeDimensions(numTolT, &dt);
	ChangeDimensions(numTolT, &dtBase);
	
	for(i = 1;i <= numTolT;i++)
		loadTol >> tInf(i) >> tSup(i) >> dtBase(i);
	
	loadTol.End();

	tolRelC = double(cicloCluster) * tolRelCBase;
	tolAbsC = double(cicloCluster) * tolAbsCBase;
	for(i = 1;i <= numTolT;i++)
		dt[i] = double(cicloCluster) * dtBase(i);
}

void myBzzCSTRNetwork::read_first_guess(char *first)
{
	BzzLoad load(first);
	load >> numComponents >> iSpec >> massFractionsInReactors;
	load.End();
	numComponentsReduced = iSpec.Size();
}

void myBzzCSTRNetwork::build_cluster(int cicloCluster, int numCSTRReactors, int countInput, int &numCluster, BzzVectorInt &cstrClusterSize)
{
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
	BzzSave saveCluster("Temp/NetworkCluster.bzz");
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
}

void myBzzCSTRNetwork::calculate_initial_massfractions_in_clusters(int numCluster, BzzMatrix &omegaReduced_Cluster, BzzMatrix &omegaReduced_CFD)
{
	int i, k, j;
	int ispecSize;
	int numC;

	BzzVector clusterVolume(numCluster);

	
	BzzLoad load("Input/FirstGuess.bzz");	
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
}

void myBzzCSTRNetwork::calculate_clustering(char *fileNetwork, int cicloDiffusion, int numCluster, 
											BzzVectorInt &cstrClusterSize, 
											BzzMatrix &omegaReduced_Cluster, 
											BzzMatrix &omegaReduced_CFD)
{
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

	openInputFile(fInput, fileNetwork);
	fInput >> numCSTRReactors;
	checkInputfile(fInput, fileNetwork);
	
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

	Mg.SetDiagonal(0,1.);
	

	// Ciclo su tutte le celle della simulazione CFD originale
	for(k = 1;k <= originalNumCSTRReactors;k++)
	{		
		fInput >> indexCell;
		fInput >> temperatureCell;
		fInput >> pressureCell;
		fInput >> volumeCell;
		fInput >> TvarianceCell;
		checkInputfile(fInput, fileNetwork);

		// Recupero indice cluster a cui appartiene la cella
		int indexCluster = giveClusterIndexFromCellIndex[indexCell];

		// Calcolo delle proprieta' medie del cluster
		temperature[indexCluster] += temperatureCell;
		pressure[indexCluster] += pressureCell;
		volume[indexCluster] += volumeCell;
		qanm[indexCluster] += TvarianceCell;
		
		fInput >> nInput;
		checkInputfile(fInput, fileNetwork);
		
		// Nel caso in cui ci siano degli input esterni alla cella di calcolo
		if(nInput != 0)
		{
			for(j = 1;j <= nInput;j++)
			{				
				fInput >> inputMassFlowRate;
				checkInputfile(fInput, fileNetwork);
		
				fInput >> numberOfSpeciesCFD;
				checkInputfile(fInput, fileNetwork);
		
				massInput[indexCluster] += inputMassFlowRate;
								
				double sumOmega = 0.;
				for(i = 1;i <= numberOfSpeciesCFD;i++)
				{
					fInput >> indexSpecies;
					fInput >> omegaSpecies;
					checkInputfile(fInput, fileNetwork);

					sumOmega += omegaSpecies;
				}
				
				if(sumOmega < .9999999 || sumOmega > 1.0000001)
					BzzError("Wrong fraction in %d rector: %22.14e", indexCluster, sumOmega);
			}
		}

		// Numero di celle collegate alla cella in esame
		fInput >> nOutput;
		checkInputfile(fInput, fileNetwork);
		
		// Ciclo su tutte le celle vicine
		for(i = 1;i <= nOutput;i++)
		{
			fInput >> destinationCell;
			fInput >> outputMassFlowRate;
			fInput >> outputDiffusionRate;
			checkInputfile(fInput, fileNetwork);
			
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
		}
	}
	fInput.close();

	// Mean Values for each cluster
	double uNumberOfCFDCells;
	for(i = 1;i <= numCSTRReactors;i++)
	{
		uNumberOfCFDCells	 = 1. / double(cstrClusterSize[i]);
		temperature[i]		*= uNumberOfCFDCells;
		qanm[i]				*= uNumberOfCFDCells;
		pressure[i]			*= uNumberOfCFDCells;
	}
	
	ProductT(volume, Reactions.M, &volumeMolecularWeightT);
}

void myBzzCSTRNetwork::calculate_massflowrate_in_reactors(BzzVector &inputMassFlowRate,
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
		double relative_error = fabs((mOUT-mIN)/(mIN+mOUT));
		meanError+=relative_error;
		if (relative_error>maxError)
			maxError = relative_error;
	}

	meanError /= numCSTRReactors;
	cout << " ** Mean relative error in mass balances on Cluster: " << meanError << endl;
	cout << " ** Max  relative error in mass balances on Cluster: " << maxError << endl;
}

void myBzzCSTRNetwork::complete_clustering(int countInput, int cicloDiffusion, int relaxation, int iaia, BzzMatrix &omegaReduced_Cluster)
{
	int i, j, k;
	double *ptrVal;
	double val;
	ElementBzzMatrixSparse *elem;
	BzzVector provisionalInputMassRate;
	BzzVector provisionalOutputMassRate;

	// Calculate provisional mass flow rate
	cout << "Relative Errors before normalization: " << endl;
	calculate_massflowrate_in_reactors(provisionalInputMassRate, provisionalOutputMassRate);

	// Vengono normalizzate le portate massive della matrice Mg e ne viene cambiato il segno e viene rifatta la trasposizione;
	// i termini sulla diagonale principale vengono settati pari al valore unitario
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
	for(i = 1;i <= numCSTRReactors;i++)
	{
		if(provisionalOutputMassRate[i] != 0.)
			massOutput[i] /=  provisionalOutputMassRate[i];
	}
	
	// Vengono calcolati i nuovi flussi convettivi che interessano ciscun cluster
	// vengono inseriti nel vettore massRate
	if(relaxation == 0 || relaxation == 2 || (relaxation == 1 && cicloDiffusion == 4))
	{
		BzzSave massRateFile('*', "Temp/MassRate.tmp");
		Fg = Mg;
		massRate = massInput;
		Solve(&Fg,&massRate);
		massRateFile << massRate;
		massRateFile.End();
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
	cout << "Relative Errors after normalization: " << endl;
	Transpose(&Mg);
	calculate_massflowrate_in_reactors(provisionalInputMassRate, provisionalOutputMassRate);
	Delete(&Mg);


	// Updating External Feeds
	{
		int indexCell;
		double massIn;
		int numOriginalSpecies, iCluster;
	
		BzzVector massFeed(numCSTRReactors);
		ChangeDimensions(numCSTRReactors,numComponents,&feed);
		BzzLoad loadClusterInput("Temp/CFDNetworkInput.bzz");
	
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

	// Cleaning Ms (sparse)
	// The elements very small are erased to increase the sparsity pattern of this
	// matrix without affecting the accuracy
	Ms.CleanMatrix(1.e-100);

	// Costruzione della matrice Ls contenente tutte le informazioni per la risoluzione
	// del sistema non lineare; questa matrice contiene anche gli elementi sulla diagonale
	// principale
	Ls = Ms;

	//Ms.BzzPrint("Ms");
	//Ls.BzzPrint("Ls");
	//sM.BzzPrint("sM");
	//exit(1);


	// Eliminazione della diagonale principale dalla matrice Ms, in modo da aumentarne
	// ulteriormente la struttura di sparsita'
	for(i = 1;i <= numCSTRReactors;i++)
		Ms.RemoveElement(i,i);

	// Matrice per Gauss Jordan
	// La matrice Ld e' identica alla matrice Ls, con la differenza che non ha pero' 
	// elementi sulla diagonale principale
	ReplaceBzzMatrixSparseWithMinusBzzMatrixSparseLocked(&Ms,&Ld);

	// Memory allocation
	ChangeDimensions(numComponents, &cReactor);
	ChangeDimensions(numComponents, &omegaReactor);
	ChangeDimensions(numComponents, &xReactor);
	ChangeDimensions(numComponents, &RReactor);

	// Inizializzazione della matrice delle frazioni massive
	if(iaia == 0)
	{
		int numC;

		BzzLoad load("Input/FirstGuess.bzz");
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
		int numC,jC;
		BzzVectorInt aui;

		BzzLoad load('*',"Temp/mass.tmp");
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
}

void myBzzCSTRNetwork::operator()(char *fileNetwork,char *first, int cicloCluster,int cicloDiffusion,int iaia,int relaxation)
{
	int countInput;
	int numCluster;
	BzzVectorInt cstrClusterSize;
	BzzMatrix omegaReduced_Cluster;
	BzzMatrix omegaReduced_CFD;

	readKineticSheme();
	read_first_guess(first);
	reading_the_network_file(fileNetwork, countInput);
	read_tolerances(numCSTRReactors, cicloCluster);
	build_cluster(cicloCluster, numCSTRReactors, countInput, numCluster, cstrClusterSize);
	calculate_initial_massfractions_in_clusters(numCluster, omegaReduced_Cluster, omegaReduced_CFD);
	calculate_clustering(fileNetwork, cicloDiffusion, numCluster, cstrClusterSize, omegaReduced_Cluster, omegaReduced_CFD);
	complete_clustering(countInput, cicloDiffusion, relaxation, iaia, omegaReduced_Cluster);

	cout << "Total number of CSTRs: " << numCSTRReactors << endl;

	start_solving_network();
}

void myBzzCSTRNetwork::start_solving_network()
{
	cout << "-----------------------------------------------------------------" << endl;
	cout << " Temperature dependences creation and memorization"				<< endl;
	
			MemoTemperatureFunctions();
	
	cout << "-----------------------------------------------------------------" << endl;
	cout << endl;

	cout << "-----------------------------------------------------------------" << endl;
	cout << " START COMPUTATIONS" << endl;
	cout << "-----------------------------------------------------------------" << endl;
	cout << endl;
	
	// Initial Solution
	massFractionsInReactorsSolution = massFractionsInReactors;

	// Global Algorithm
	globalAlgorithm();
}

void myBzzCSTRNetwork::globalAlgorithm()
{
	// Convergence Criteria
	F1Stop = 1.e-6;
//	F1Stop = 1.e-7;
//	F1Stop = 1.e-15;
	maxNumberOfIterations_Global = Max(numCSTRReactors / 3,30);
	
	int okok;
	int countCicle = 1;
	for(int iter = 1;iter <= 4;iter++)
	{
		// 1. Vengono risolti singolarmente tutti i reattori CSTR
		//    L'ordine di risoluzione non e' importante e viene praticamente mantenuto 
		//    quello della simulazione originale CFD
		//    0 = se la convergenza non e' stata raggiunta in maniera appropriata
		//	  1 = se la convergenza e' stata ottenuta in maniera corretta
		okok = GetFirst();

		// 2. Stampa le informazioni a video
		OutputPrint(countCicle++);
		
		// Se le iterazioni sui singoli reattori sono andate a buon fine e si sono fatte 
		// almeno tre iterazioni globali del metodo, allora la rete di reattori puo'
		// considerarsi risolta ed e' possibile uscire dal programma
		if(okok == 1 && iter == 3)
			return;
		
		// 3. Viene risolta l'intera rete di reattori
		cout << " ITERATION GET THIRD " << endl;
		okok = GetThird();

		// 4. Stampa le informazioni a video
		OutputPrint(countCicle++);
		
		// Se il sistema globale ha raggiunto la convergenza e' possibile uscire
		if(okok == 1)
			return;
		
		// 5. Viene chiamata solo a partire dalla seconda iterazione in poi
		if(iter >= 2)
		{
			GetSecond();
			okok = GetThird();
			OutputPrint(countCicle++);
			if(okok == 1)
				return;
		}
		
		// New strict tolerances
		tolRel *= .1;
		tolAbs *= .1;
		if(tolRel < 1.e-5)
			tolRel = 1.e-5;
		if(tolAbs < 1.e-9)
			tolAbs = 1.e-9;
	}
}

void myBzzCSTRNetwork::GetResiduals(BzzMatrix &massFractionsInReactors, BzzMatrix &residuals)
{	
	Product(Ls,massFractionsInReactors,&residuals);
	GetReactionsRateInAllReactorsFromMassFractions(massFractionsInReactors, dummyMatrix);	
	ElementByElementProduct(volumeMolecularWeightT,dummyMatrix,&dummyMatrix);
	Difference(dummyMatrix,&residuals);
	Sum(feed,&residuals);
}

void myBzzCSTRNetwork::GetAllResiduals(BzzVector &omega,BzzVector &res)
{
	Product(Ls,omega,&res);
	GetReactionsRateInAllReactorsFromMassFractions(omega,dummyVector);
	ElementByElementProduct(volumeMolecularWeightT_Vector,dummyVector,&dummyVector);
	Difference(dummyVector,&res);
	Sum(feed_Vector,&res);
}

void myBzzCSTRNetwork::GetReactionRatesInSingleReactor(int iReactor,
							BzzVector &omega, BzzVector &molefractions,
							BzzVector &concentrations, BzzVector &ReactionRates)
{
	wM = Reactions.GetMWFromMassFractions(omega);
	Reactions.GetMoleFractionsFromMassFractionsAndMW(molefractions,omega,wM);
	concentrations = cTot*molefractions;
	Reactions.ComputeFromConcentrations_map(iReactor, concentrations, &ReactionRates);
}

void myBzzCSTRNetwork::GetResiduals(BzzVector &m,BzzVector &residuals)
{
	GetReactionRatesInSingleReactor(kReactor, m, xReactor, cReactor, RReactor);
	
	volumeMolecularWeightT.GetRow(kReactor,&cReactor);
	ElementByElementProduct(RReactor,cReactor,&xReactor);
	Product(sM[kReactor],m,&cReactor);
	Difference(&cReactor,feedInkRreactor);
	Difference(cReactor,xReactor,&residuals);
}

void myBzzCSTRNetwork::GetJacobian(BzzVector &omega,BzzMatrix &JJ)
{
	GetReactionRatesInSingleReactor(kReactor, omega, xReactor, cReactor, RReactor);
	GetJacobianForSingleReactor(kReactor, cReactor, RReactor, JJ);
}

void myBzzCSTRNetwork::GetJacobianForSingleReactor(int iReactor, BzzVector &cReactor, BzzVector &R, BzzMatrix &dRC)
{
	int i, j;
	double wc;

	// Calculation of the Jacobian
	Reactions.kinetics.mc		= cReactor;
	Reactions.kinetics.mR		= R;
	Reactions.kinetics.mr		= Reactions.r;
	Reactions.kinetics.mrDirT	= Reactions.coeffFallOff;
	Reactions.kinetics.mrDirC	= Reactions.rDirC;
	Reactions.kinetics.mrInvC	= Reactions.rInvC;

	// 1. Contributi dai termini di reazione
	Reactions.kinetics.GetDerivativesC(T, cTot, &dRC, cReactor, R);
	wc = wM * cTot;
	for(i = 1;i <= numComponents;i++)
		xReactor[i] = wc / Reactions.M[i];
	dRC.ColumnsElementsProduct(xReactor);
	for(i = 1;i <= numComponents;i++)
		for(j = 1;j <= numComponents;j++)
			dRC[i][j] *= volumeMolecularWeightT[iReactor][i];
		
	// 2. Contributo derivante dal termine convettivo
	for(i = 1;i <= numComponents;i++)
		dRC[i][i] -= sM[iReactor];
}


int myBzzCSTRNetwork::GetFirst(void)
{
	// Indici
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
	MyOdeSystemObjectOneCSTR cstrMono;
	cstrMono(this);
	
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
	ripetiTutto:

	// Reset Counts
	countTutto++;
	countIntegrations			= 0;
	countReactorsOK_FirstGuess	= 0;
	F1Max						= 0.;
	numF16						= 0;
	F1MaxIntegration			= 0.;
	kMaxIntegration				= 0;
	maxOdeSum = .1 * F1Stop * double(numComponents);

	// Loading MemoTemperatureFunctions.tmp
	if(memoTemperature == 1)
		load('*',"Temp/MemoTemperatureFunctions.tmp");
	
	// Cycle on every reactor
	for(kReactor = 1;kReactor <= numCSTRReactors;kReactor++)
	{
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
			load >> Reactions.uKeq;
			load >> Reactions.k1 >> Reactions.k2;
			load >> Reactions.logFcent;
		}

		// Recovering Actual solution in reactor
		massFractionsInReactorsSolution.GetRow(kReactor,&x0);

		// Checking on Sum mass fractions
		double su = x0.GetSumAbsElements();
		if(su > 1.0001 || su < .9999)
		{
			cout << endl << "WARNING: Reactor " << kReactor << " - Sum: " << su << endl;
			if(su == 0.)
				BzzError("FATAL ERROR: Wrong Composition in Reactor %d", kReactor);
			
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
					// A1. Se l variazione non e' piccolissima allora e' meglio abbandonare
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
			ode.SetAnalyticalJacobian();						// TODO
			ode.SetTollAbs(1.e-15);
			ode.StopIntegrationBeforeRecalcuatingJacobian(3);	// TODO 
			ode.StopIntegrationWhenSumAbsY1IsLessThan(maxOdeSum);
			BzzVector yMin(numComponents);
			ode.SetMinimumConstraints(&yMin);
			omegaReactor = ode(10000.*tau);
			
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
			if(sumOmega == 0.)	BzzError("FATAL ERROR: CODE 0001");
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
			if(countNewtonIterations < 10 && F1<OldF0_Corrected && F1>MAX_SUM_F)
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
			ode.SetTollAbs(1.e-15);
			ode.StopIntegrationBeforeRecalcuatingJacobian(3);	// TODO
			ode.StopIntegrationWhenSumAbsY1IsLessThan(maxOdeSum);
			BzzVector yMin(numComponents);
			ode.SetMinimumConstraints(&yMin);
			omegaReactor = ode(10000.*tau);
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
	cout << "Total number of reactors which are not OK: " << nReactorsNotOK << endl;
	cout << "Total number of species  which are not OK: " << nVariablesNotOK << endl;

	memoStart = BzzGetCpuTime();
	
	// Convergence (1)

	if(	nReactorsNotOK < ::Max(10, numCSTRReactors / 20)		&& 
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

	// Convergence (2)
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
	// e' stato applicato senza che sia stao necessario ricorrere all'integrazione
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

	// If none of these convergence criteria is satisfied, repeat the global iteration
	goto ripetiTutto;
}

// La rete di reattori viene risolta GLOBALMENTE attraverso un sistema ODE
void myBzzCSTRNetwork::GetSecond(void)
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
	cout << "Solving Linear System in Second Method" << endl;
	cout << "  ** Mean Residual Before Second Method" << F0 << endl;
	cout << "  ** Max  Residual Before Second Method" << maxF0 << endl;

	
	// Solving the ODE System
	{
		double tStart,tEnd;
		double stopODE;
		BzzVector yMin(initialValues.Size());
		BzzVector yMax(initialValues.Size());
		MyOdeSystemObjectAllCSTR obj;

		tStart = 0.;
		tEnd = 10.;
		stopODE = tolAbs * uvariables;
	
		BzzMatrixSparseLockedByRows Q = Ld;
		obj(this);

		BzzOdeSparseStiffObject o(initialValues, tStart, &obj, blockDimensions, &Q, "A.tmp");
		yMin = 0.;
		yMax = 1.;
		o.SetMinimumConstraints(yMin);
		o.SetMaximumConstraints(yMax);
		o.SetTollAbs(1.e-15);
		o.StopIntegrationBeforeRecalcuatingJacobian(30);
		o.StopIntegrationWhenSumAbsY1IsLessThan(stopODE);
		
		massFractionsInReactorsSolution_Vector = o(tEnd);
	
		cout	<< "CPU Time for Solving Linear System in Second Method: " 
				<< BzzGetCpuTime() - start << endl;
	}

	// Calcolo dei nuovi residui dopo l'integrazione del sistema ODE
	GetAllResiduals(massFractionsInReactorsSolution_Vector, residual_Vector);
	F1 = residual_Vector.GetSumAbsElements() * uvariables;
	maxF1 = residual_Vector.MaxAbs();
	cout << "  ** Mean Residual After Second Method" << F1 << endl;
	cout << "  ** Max  Residual After Second Method" << maxF1 << endl;

	// Accetta la nuova soluzione che sicuramente e' buona in quanto ottenuta
	// attraverso la soluzione di un transitorio reale
	CopyDataFromVector(&massFractionsInReactorsSolution,massFractionsInReactorsSolution_Vector);
}

// Viene risolto il sistema lineare GLOBALE associato a tutta la rete di reattori.
// Per fare questa operazione si sfrutta la classe BzzFactoredDoubleDiagonalBlockGaussAndMatrixLocked
// che prevede che la matrice del sistema abbia sulla diagonale blocchi anche di diverse
// dimensioni e i termini extradiagonali piccoli se confrontati con i corrispondenti 
// sulla diagonale principale. L'algoritmo e' un po' complicato perche' e' necessario
// salvare la fattorizzazione su file onde evitare problemi di allocazioni di memoria
int myBzzCSTRNetwork::GetThird(void)
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
	BzzMatrixSparseLockedByRows					ExtraDiagonalTermsMatrix;
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
	cout << "Solving Linear System in Third Method - Iteration: " << countThird << endl;
	
	if(memoTemperature == 1)
		load('*',"Temp/MemoTemperatureFunctions.tmp");
	
	// Viene aperto il file su cui verra' scritta la fattorizzaizone delle matrici
	// che costituiscono i blocchi della matrice del sistema lineare
	// In particolare viene scritto prima di tutto il numero di blocchi, ovvero il
	// numero di reattori della rete
	BzzSave save('*',"Temp/MemoJacobian.tmp");
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
			load >> Reactions.uKeq;
			load >> Reactions.k1 >> Reactions.k2;
			load >> Reactions.logFcent;
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


	// La matrice Ld e' una BzzMatrixDoubleSparseLockedByRows che di fatto corrisponde
	// esattamente alla matrice Ls, ma senza i termini sulla diagonale principale; questo perche'
	// la definizione della matrice Pcstr prevede che vengano fornite separamente le matrici
	// dei blocchi della diagonale (che in questo caso e' presa dal file) e gli elementi che
	// si trovano fuori
	ExtraDiagonalTermsMatrix = Ld;
	GlobalMatrix("Temp/MemoJacobian.tmp", &ExtraDiagonalTermsMatrix);
	
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


	double startLocal = BzzGetCpuTime();

	// Viene ricavato il vettore dei termini noti del sistema lineare per poter applicare
	// il metodo di Newton
	residuals.UseMatrixAsVector(&d1);

	// Viene applicato il metodo di Newton per ottenere la correzione 
//	int fail = GlobalMatrix.GaussJordanSolve(&d1);	
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
		load('*',"Temp/MemoTemperatureFunctions.tmp");
	save('*',"Temp/MemoJacobian.tmp");
	
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
			load >> Reactions.uKeq;
			load >> Reactions.k1 >> Reactions.k2;
			load >> Reactions.logFcent;
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

	// Se la norma della correzione e' piu' piccola di quella 
	// del punto precedente e non ho fatto molte iterazioni conviene accettare il
	// nuovo punto e ripetere l'iterazione
	if(F1dx < F0dx && countThird < 10)
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
			BzzWarning("TODO: Solving Linear System in Third Method - CODE 0001");
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


void myBzzCSTRNetwork::ObjectBzzPrint(void)
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
		::BzzPrint("\n%5d   %s",j,Reactions.strReaction[j]);
	
	for(i = 1;i <= numCSTRReactors;i++)
	{
		k = cstrSequence[i];

		cTot = pressure[i] / (temperature[i] * R_CTOT);
		massFractionsInReactorsSolution.GetRow(i,&ma);
		wM = Reactions.GetMWFromMassFractions(ma);
		Reactions.GetMoleFractionsFromMassFractionsAndMW(xReactor,ma,wM);
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
			cout << j << " " << Reactions.names[j]	<< "\t"
				 << cReactor[j]							<< "\t"
				 << xReactor[j]						<< "\t"
				 << ma[j]							<< "\t"
				 << res[i][j]						<< "\t"
				 << endl;
		}
	}

	OutputPrint(0);
}

void myBzzCSTRNetwork::OutputPrint(int cicle)
{
	int i,k;
	double flowrate;
	double totalOutFlowRate = 0.;
	BzzVector outFlowRate(numComponents);

	cout.setf(ios::scientific);
	Save('*',"Temp/mass.tmp");
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
				<< Reactions.names[i]			<< "\t" 
				<< outFlowRate[i]				<< "\t" 
				<< outFlowRate[i] / flowrate	<< endl;

	cout	<< endl							<< "\t"
			<< "TOTAL"						<< "\t" 
			<< totalOutFlowRate				<< "\t" 
			<< totalOutFlowRate / flowrate	<< endl
			<< endl;
	
	cout << "=====================================================" << endl;
}

void myBzzCSTRNetwork::MemoTemperatureFunctions(void)
{
	int j,k,kCSTR;
	double uT;
	double uRT;

	BzzVector u_temperature(numCSTRReactors);
	BzzVector log_temperature(numCSTRReactors);
	
	BzzMatrix matrix_correction_k1(numCSTRReactors,numReactions);
	BzzMatrix matrix_correction_uKeq(numCSTRReactors,numReactions);
	BzzMatrix matrix_correction_k2(numCSTRReactors,numReactions);
	

	// INIZIALIZZAZIONE DELLE VARIABILI DIPENDENTI SOLO DALLA TEMPERATURA
	cout << "Inizializzazione delle costanti cinetiche... ";

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

	Reactions.InitializeMap(numCSTRReactors);
	Reactions.ComputeKineticParameters_map(temperature, cTotm, pressure);	

	for(kCSTR = 1;kCSTR <= numCSTRReactors;kCSTR++)
	{
		T = temperature[kCSTR];
		P = pressure[kCSTR];
		cTot = cTotm[kCSTR];
		uT = u_temperature[kCSTR];
		logT = log_temperature[kCSTR];
		uRT = uT * UR_CTOT;
		loguRT = log(uRT);

		double CoeffCorr;
		for(k = 1;k <= Reactions.kinetics.numEquilibrium;k++)
		{
			j = Reactions.kinetics.reactionWithEquil[k];
			Reactions.uKeq[j] = exp(-Reactions.reactionDS_map[kCSTR][j] + Reactions.reactionDH_map[kCSTR][j] - loguRT * Reactions.kinetics.sumNuij[j]);
			{
				double n,n2,n3,n4,n5,n6;
				double EsuRT,EsuRT2,EsuRT3,EsuRT4,EsuRT5,EsuRT6;
				double sommaTer2,sommaTer4,sommaTer6;
				double qan, qan2, qan3;
				EsuRT = -(Reactions.kinetics.E1[j] + T*Reactions.reactionDH_map[kCSTR][j]) * uT;
				n = Reactions.kinetics.beta1[j] + Reactions.kinetics.sumNuij[j];
				qan = qanm[kCSTR];
				n2=n*n;
				n3=n2*n;
				n4=n3*n;
				n5=n4*n;
				n6=n5*n;

				EsuRT2=EsuRT*EsuRT;
				EsuRT3=EsuRT2*EsuRT;
				EsuRT4=EsuRT3*EsuRT;
				EsuRT5=EsuRT4*EsuRT;
				EsuRT6=EsuRT5*EsuRT;
				qan2=qan*qan;
				qan3=qan2*qan;
				sommaTer2=-n+n2+EsuRT*(-2+2*n)+EsuRT2;
				sommaTer4=-6*n+11*n2-6*n3+n4+EsuRT*(-24+44*n-24*n2+4*n3)+EsuRT2*(36-30*n+6*n2)+EsuRT3*(-12+4*n)+EsuRT4;
				sommaTer6=-120*n+274*n2-225*n3+85*n4-15*n5+n6+EsuRT*(-720+1644*n-1350*n2+510*n3-90*n4+6*n5);
				sommaTer6+=EsuRT2*(1800-2310*n+1065*n2-210*n3+15*n4)+EsuRT3*(-1200+940*n-240*n2+20*n3);
				sommaTer6+=EsuRT4*(300-135*n+15*n2)+EsuRT5*(-30+6*n)+EsuRT6;
				CoeffCorr=1.+0.25*sommaTer2*qan+0.015625*sommaTer4*qan2+4.340278E-4*sommaTer6*qan3;
			}

			matrix_correction_uKeq[kCSTR][j] = CoeffCorr;
		}

		for(j = 1;j <= numReactions;j++)
		{
			Reactions.k1[j] = exp(Reactions.kinetics.k01[j] + Reactions.kinetics.beta1[j] * logT + Reactions.kinetics.E1[j] * uT);
			{
				double n,n2,n3,n4,n5,n6;
				double EsuRT,EsuRT2,EsuRT3,EsuRT4,EsuRT5,EsuRT6;
				double sommaTer2,sommaTer4,sommaTer6;
				double qan, qan2, qan3;
				EsuRT = -Reactions.kinetics.E1[j] * uT;
				n = Reactions.kinetics.beta1[j];
				qan = qanm[kCSTR];
				n2=n*n;
				n3=n2*n;
				n4=n3*n;
				n5=n4*n;
				n6=n5*n;

				EsuRT2=EsuRT*EsuRT;
				EsuRT3=EsuRT2*EsuRT;
				EsuRT4=EsuRT3*EsuRT;
				EsuRT5=EsuRT4*EsuRT;
				EsuRT6=EsuRT5*EsuRT;
				qan2=qan*qan;
				qan3=qan2*qan;
				sommaTer2=-n+n2+EsuRT*(-2+2*n)+EsuRT2;
				sommaTer4=-6*n+11*n2-6*n3+n4+EsuRT*(-24+44*n-24*n2+4*n3)+EsuRT2*(36-30*n+6*n2)+EsuRT3*(-12+4*n)+EsuRT4;
				sommaTer6=-120*n+274*n2-225*n3+85*n4-15*n5+n6+EsuRT*(-720+1644*n-1350*n2+510*n3-90*n4+6*n5);
				sommaTer6+=EsuRT2*(1800-2310*n+1065*n2-210*n3+15*n4)+EsuRT3*(-1200+940*n-240*n2+20*n3);
				sommaTer6+=EsuRT4*(300-135*n+15*n2)+EsuRT5*(-30+6*n)+EsuRT6;
				CoeffCorr=1+0.25*sommaTer2*qan+0.015625*sommaTer4*qan2+4.340278E-4*sommaTer6*qan3;
			}
			
			matrix_correction_k1[kCSTR][j] = CoeffCorr;
			
			if(Reactions.kinetics.jThirdBody[j] >= 2)
			{
				Reactions.k2[j] = exp(Reactions.kinetics.k02[j] + Reactions.kinetics.beta2[j] * logT + Reactions.kinetics.E2[j] * uT);
				{
					double n,n2,n3,n4,n5,n6;
					double EsuRT,EsuRT2,EsuRT3,EsuRT4,EsuRT5,EsuRT6;
					double sommaTer2,sommaTer4,sommaTer6;
					double qan, qan2, qan3;
					EsuRT = -Reactions.kinetics.E2[j] * uT;
					n = Reactions.kinetics.beta2[j];
					qan = qanm[kCSTR];
					n2=n*n;
					n3=n2*n;
					n4=n3*n;
					n5=n4*n;
					n6=n5*n;

					EsuRT2=EsuRT*EsuRT;
					EsuRT3=EsuRT2*EsuRT;
					EsuRT4=EsuRT3*EsuRT;
					EsuRT5=EsuRT4*EsuRT;
					EsuRT6=EsuRT5*EsuRT;
					qan2=qan*qan;
					qan3=qan2*qan;
					sommaTer2=-n+n2+EsuRT*(-2+2*n)+EsuRT2;
					sommaTer4=-6*n+11*n2-6*n3+n4+EsuRT*(-24+44*n-24*n2+4*n3)+EsuRT2*(36-30*n+6*n2)+EsuRT3*(-12+4*n)+EsuRT4;
					sommaTer6=-120*n+274*n2-225*n3+85*n4-15*n5+n6+EsuRT*(-720+1644*n-1350*n2+510*n3-90*n4+6*n5);
					sommaTer6+=EsuRT2*(1800-2310*n+1065*n2-210*n3+15*n4)+EsuRT3*(-1200+940*n-240*n2+20*n3);
					sommaTer6+=EsuRT4*(300-135*n+15*n2)+EsuRT5*(-30+6*n)+EsuRT6;
					CoeffCorr=1+0.25*sommaTer2*qan+0.015625*sommaTer4*qan2+4.340278E-4*sommaTer6*qan3;
				}
				
				matrix_correction_k2[kCSTR][j] = CoeffCorr;
			}
		}

	}	// End Cycle On each Reactor


	Reactions.CorrectKineticParameters_map(matrix_correction_k1, matrix_correction_uKeq, matrix_correction_k2);
	cout << "DONE" << endl;

	// Writing on File
	if(memoTemperature == 1)
	{
		#if WINDOWS
		
		BzzSave save('*',"Temp/MemoTemperatureFunctions.tmp");
		for(kCSTR = 1;kCSTR <= numCSTRReactors;kCSTR++)
		{
			save << Reactions.uKeq_map.GetRow(kCSTR);
			save << Reactions.k1_map.GetRow(kCSTR) << Reactions.k2_map.GetRow(kCSTR);
			save << Reactions.logFcent_map.GetRow(kCSTR);
		}
		// Closing File
		save.End();

		#else

		BzzSave save('*',"Temp/MemoTemperatureFunctions.tmp");
		for(kCSTR = 1;kCSTR <= numCSTRReactors;kCSTR++)
		{
			BzzVector auxiliar;

			auxiliar = Reactions.uKeq_map.GetRow(kCSTR);
			save << auxiliar;
		
			auxiliar = Reactions.k1_map.GetRow(kCSTR);	
			save << auxiliar;

			auxiliar = Reactions.k2_map.GetRow(kCSTR);
			save << auxiliar;

			auxiliar = Reactions.logFcent_map.GetRow(kCSTR);
			save << auxiliar;
		}
		// Closing File
		save.End();
		#endif
	}
}

void myBzzCSTRNetwork::GetReactionsRateInAllReactorsFromMassFractions
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
		load('*',"Temp/MemoTemperatureFunctions.tmp");

	for(kr = 1;kr <= numCSTRReactors;kr++)
	{
		T = temperature[kr];
		P = pressure[kr];

		logT = logTm[kr];
		loguRT = loguRTm[kr];
		cTot = cTotm[kr];
		
		if(memoTemperature == 1)
		{
			load >> Reactions.uKeq;
			load >> Reactions.k1 >> Reactions.k2;
			load >> Reactions.logFcent;
		}

		massFractionsInReactors.UseMatrixRowAsVector(kr,&ma);

		wM = Reactions.GetMWFromMassFractions(ma);
		Reactions.GetMoleFractionsFromMassFractionsAndMW(xReactor,ma,wM);
		cReactor = cTot*xReactor;
		reactionsRateInReactors.UseMatrixRowAsVector(kr,&R);
		Reactions.ComputeFromConcentrations_map(kr, cReactor, &R);
	}
	if(memoTemperature == 1)
		load.End();
}	

void myBzzCSTRNetwork::GetReactionsRateInAllReactorsFromMassFractions
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
		load('*',"Temp/MemoTemperatureFunctions.tmp");

	for(kr = 1;kr <= numCSTRReactors;kr++)
	{
		T = temperature[kr];
		P = pressure[kr];
		logT = logTm[kr];
		loguRT = loguRTm[kr];
		cTot = cTotm[kr];

		if(memoTemperature == 1)
		{
			load >> Reactions.uKeq;
			load >> Reactions.k1 >> Reactions.k2;
			load >> Reactions.logFcent;
		}

		ma.UseSubVectorAsVector(numComponents, (kr - 1)*numComponents + 1, massFractionsInReactorsV);
		wM = Reactions.GetMWFromMassFractions(ma);
		Reactions.GetMoleFractionsFromMassFractionsAndMW(xReactor,ma,wM);	
		cReactor = cTot*xReactor;
		R.UseSubVectorAsVector(numComponents,(kr - 1)*numComponents + 1,reactionsRateInReactorsV);
		Reactions.ComputeFromConcentrations_map(kr, cReactor, &R);
	}

	if(memoTemperature == 1)	load.End();
}	


void  myBzzCSTRNetwork::GetDiagonalFactoredMatricesForLinearizedSistem
		(BzzMatrix &massFractionsInReactors)
{
	int kr;
	BzzLoad load;
	BzzMatrix dRC(numComponents,numComponents);
	BzzVector R(numComponents);
	BzzFactorizedGauss G;

	if(memoTemperature == 1)
		load('*',"Temp/MemoTemperatureFunctions.tmp");
	BzzSave save('*',"Temp/MemoJacobian.tmp");

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
			load >> Reactions.uKeq;
			load >> Reactions.k1 >> Reactions.k2;
			load >> Reactions.logFcent;
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

void  myBzzCSTRNetwork::GetDiagonalMatricesForLinearizedSistem(char *file,BzzVector &mfV)
{
	int kr;
	BzzLoad load;

	BzzMatrix dRC(numComponents,numComponents);
	BzzVector R(numComponents);
	BzzFactorizedGauss G;
	BzzVector ma;

	if(memoTemperature == 1)
		load('*',"Temp/MemoTemperatureFunctions.tmp");
	
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
			load >> Reactions.uKeq;
			load >> Reactions.k1 >> Reactions.k2;
			load >> Reactions.logFcent;
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


void myBzzCSTRNetwork::SwapMassFractionsInReactors(BzzMatrix *massFractionsInReactors)
{
	ReorderByRows(massFractionsInReactors,cstrSequence);
}

void myBzzCSTRNetwork::Save(char *file) // formatted
{
	BzzSave save(file);
	save << numCSTRReactors;
	save << giveClusterIndexFromCellIndex;
	save << massFractionsInReactorsSolution;
	save.End();
}

void myBzzCSTRNetwork::Save(char,char *file)// binary
{
	BzzSave save('*',file);
	save << numCSTRReactors;
	save << giveClusterIndexFromCellIndex;
	save << massFractionsInReactorsSolution;
	save.End();
}
