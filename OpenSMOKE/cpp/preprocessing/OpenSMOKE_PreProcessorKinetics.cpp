/***************************************************************************
 *   Copyright (C) 2006 by Alberto Cuoci   	   *
 *   alberto.cuoci@polimi.it   						   *
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

#include "preprocessing/OpenSMOKE_PreProcessorReactingGas.h"
#include "preprocessing/OpenSMOKE_PreProcessorKinetics.h"
#include "basic/OpenSMOKE_Constants.h"
#include <sstream>

const int		OpenSMOKE_PreProcessorKinetics::NUM_MAX	= 40000;

void OpenSMOKE_PreProcessorKinetics::ErrorMessage(const std::string message)
{
    cout << endl;
    cout << "Class:  OpenSMOKE_PreProcessorKinetics"	<< endl;
    cout << "Object: " << name_object		<< endl;
    cout << "Error:  " << message           << endl;
    cout << "Press enter to continue... "   << endl;
    getchar();
    exit(-1);
}

void OpenSMOKE_PreProcessorKinetics::WarningMessage(const std::string message)
{
    cout << endl;
    cout << "Class:	  OpenSMOKE_PreProcessorKinetics"	<< endl;
    cout << "Object:  " << name_object		<< endl;
    cout << "Warning: " << message          << endl;
	cout << endl;
}

void OpenSMOKE_PreProcessorKinetics::readFromFile(const std::string fileSt, const std::string fileKin)
{
	int i,j,k;
	
	// Recupero del numero di reazioni e del numero di specie
	// -------------------------------------------------------------------------------
	NR = reactionRates->NumberOfReactions();
	NC = reactionRates->NumberOfSpecies();
	
	// Inizializzazione dei contatori per l'analisi dei tipi di reazioni che
	// compaiono nello schema cinetico
	// -------------------------------------------------------------------------------
	numMortoCheParla		= 0;
	numThirdBodyOnly		= 0;
	numFallOff				= 0;
	numThirdBody			= 0;
	numConventional			= 0;
	numIrreversible			= 0;
	numEquilibrium			= 0;
	numCABR					= 0;
	numLandauTeller			= 0;
	numPowerSeries			= 0;
	numJanevLanger			= 0;
	numChebishev			= 0;
	numLogarithmicPressure	= 0;
	numCollisionEfficiency	= 0;
	numTAR					= 0;

	// Dimensionamento dei reattori
	// -------------------------------------------------------------------------------
	iEfficiency	 =	 new BzzVectorInt [NR + 1];
	efficiency	 =	 new BzzVector [NR + 1];

	reactionRates->strReaction  =   new char* [NR + 1];
	for(j = 1;j <= NR;j++)
		reactionRates->strReaction[j] = new char[Constants::REACTION_NAME_SIZE];


	ChangeDimensions(NR,&sumNuij);
	ChangeDimensions(NR,&forwardOrders);
	ChangeDimensions(NR,&backwardOrders);
	ChangeDimensions(NR,&jEquil);
	ChangeDimensions(NR,&jThirdBody);
	ChangeDimensions(NR,&jNumEfficiency);
	ChangeDimensions(NR,&beta1);
	ChangeDimensions(NR,&beta2);
	ChangeDimensions(NR,&E1);
	ChangeDimensions(NR,&E2);
	ChangeDimensions(NR,&k01);
	ChangeDimensions(NR,&exp_k01);
	ChangeDimensions(NR,&k02);
	ChangeDimensions(NR,&aPressure);
	ChangeDimensions(NR,&bPressure);
	ChangeDimensions(NR,&cPressure);
	ChangeDimensions(NR,&dPressure);
	ChangeDimensions(NR,&ePressure);
	ChangeDimensions(NC,&numDir1);
	ChangeDimensions(NC,&numDir2);
	ChangeDimensions(NC,&numDir3);
	ChangeDimensions(NC,&numDir4);
	ChangeDimensions(NC,&numDir5);
	ChangeDimensions(NC,&numInvTot1);
	ChangeDimensions(NC,&numInvTot2);
	ChangeDimensions(NC,&numInvTot3);
	ChangeDimensions(NC,&numInvTot4);
	ChangeDimensions(NC,&numInvTot5);
	ChangeDimensions(NC,&numInvEq1);
	ChangeDimensions(NC,&numInvEq2);
	ChangeDimensions(NC,&numInvEq3);
	ChangeDimensions(NC,&numInvEq4);
	ChangeDimensions(NC,&numInvEq5);

	ChangeDimensions(NUM_MAX,&iFallOff);
	ChangeDimensions(NUM_MAX,&iCABR);
	ChangeDimensions(NUM_MAX,&iThirdBody);
	ChangeDimensions(NUM_MAX,&iThirdBodyOnly);
	ChangeDimensions(NUM_MAX,&iLandauTeller);
	ChangeDimensions(NUM_MAX,&iPowerSeries);
	ChangeDimensions(NUM_MAX,&iJanevLanger);
	ChangeDimensions(NUM_MAX,&iChebishev);
	ChangeDimensions(NUM_MAX,&iLogarithmicPressure);
	ChangeDimensions(NUM_MAX,&iCollisionEfficiency);
	ChangeDimensions(NUM_MAX,&iTAR);


	// -------------------------------------------------------------------------------
	// 3. Lettura del file REAZ.BZZ
	// -------------------------------------------------------------------------------
	ifstream inputFile;
	openInputFileAndControl(inputFile, fileKin);
	
	BzzSave outputFile;
	BzzSave asciiFile;
	outputFile('*', "reactions.bin");
	asciiFile("reactions.ascii");
	outputFile << NC;
	outputFile << NR;
	asciiFile << NC;
	asciiFile << NR;

	
	for(j = 1;j <= NR;j++)
	{
		inputFile >> reactionRates->strReaction[j];
		if (inputFile.eof()!=0) ErrorMessage("End of file");

		inputFile  >> jEquil[j] >> jThirdBody[j] >> jNumEfficiency[j];
		outputFile << jEquil[j] << jThirdBody[j] << jNumEfficiency[j];
		asciiFile << jEquil[j] << jThirdBody[j] << jNumEfficiency[j];
		if (inputFile.eof()!=0) ErrorMessage("End of file");

		if(jNumEfficiency[j] == -1)
		{
			BzzWarning("Morto che parla");
			ChangeDimensions(1,&iEfficiency[j]);
			ChangeDimensions(1,&efficiency[j]);
			inputFile >> iEfficiency[j][1] >> efficiency[j][1];
			outputFile << iEfficiency[j][1] << efficiency[j][1];
			asciiFile << iEfficiency[j][1] << efficiency[j][1];
			
			numMortoCheParla++;
			if (inputFile.eof()!=0) ErrorMessage("End of file");
		}
		
		// Reazioni con il terzo corpo: lettura delle efficienze
		else if(jNumEfficiency[j]  > 0)
		{
			ChangeDimensions(jNumEfficiency[j],&iEfficiency[j]);
			ChangeDimensions(jNumEfficiency[j],&efficiency[j]);
			for(int ke = 1;ke <= jNumEfficiency[j];ke++)
			{
				inputFile  >> iEfficiency[j][ke] >> efficiency[j][ke];
				outputFile << iEfficiency[j][ke] << efficiency[j][ke];
				asciiFile << iEfficiency[j][ke] << efficiency[j][ke];
				efficiency[j][ke]-=1.;
				if (inputFile.eof()!=0) ErrorMessage("End of file");
			}
		}

		// Aggiornamento delle variabili per il conteggio e l'analisi dei tipi di
		// reazione coinvolte nello schema cinetico
		if(jThirdBody[j]==0 && jNumEfficiency[j]==0) numConventional++;
		else
		{
			// - tutte le reazioni di FallOff sono automaticamente anche di terzo corpo
			//   e dunque si rende comunque necessario il calcolo del coefficiente
			//   correttivo coeffM

			// la reazione solo di terzo corpo
			if(jThirdBody[j]==1  || jThirdBody[j]==10 || 
			   jThirdBody[j]==11 || jThirdBody[j]==12 ||
			   jThirdBody[j]==13 || jThirdBody[j]==14  )
			{
				numThirdBodyOnly++;
				iThirdBodyOnly[numThirdBodyOnly]=j;
			}

			// reazioni con terzo corpo (sia quelle pure che di falloff che cabr)
			if(jThirdBody[j]<100) 
			{  
				numThirdBody++;
				iThirdBody[numThirdBody]=j;
			}

			// in questo caso la reazione e' sicuramente di FallOff e quindi in realta
			// anche di terzo corpo
			if(jThirdBody[j]>=2 && jThirdBody[j]<=4) 
			{
				numFallOff++;
				iFallOff[numFallOff]=j;
			}

			if(jThirdBody[j]>=5 && jThirdBody[j]<=7) 
			{
				numCABR++;
				iCABR[numCABR]=j;
			}

			if(jThirdBody[j]==10 || jThirdBody[j]==100) 
			{
				numLandauTeller++;
				iLandauTeller[numLandauTeller]=j;
			}

			if(jThirdBody[j]==11 || jThirdBody[j]==110) 
			{
				numJanevLanger++;
				iJanevLanger[numJanevLanger]=j;
			}

			if(jThirdBody[j]==12 || jThirdBody[j]==120) 
			{
				numPowerSeries++;
				iPowerSeries[numPowerSeries]=j;
			}

			if(jThirdBody[j]==13 || jThirdBody[j]==130) 
			{
				numChebishev++;
				iChebishev[numChebishev]=j;
			}

			if(jThirdBody[j]==14 || jThirdBody[j]==140) 
			{
				numLogarithmicPressure++;
				iLogarithmicPressure[numLogarithmicPressure]=j;
			}

			if(jThirdBody[j]==15 || jThirdBody[j]==150) 
			{
				numCollisionEfficiency++;
				iCollisionEfficiency[numCollisionEfficiency]=j;
			}

			if(jThirdBody[j]==160) 
			{
				numTAR++;
				iTAR[numTAR]=j;
			}
		}

		inputFile  >> k01[j] >> beta1[j] >> E1[j];
		outputFile << k01[j] << beta1[j] << E1[j];
		asciiFile << k01[j] << beta1[j] << E1[j];
		
		if (inputFile.eof()!=0) ErrorMessage("End of file");

		exp_k01[j] = k01[j];
		
		if(k01[j] <= 0.)
		{
			stringstream number;
			number << j;
			std::string message = "Reaction #" + number.str() +"\n";
			cout << k01[j] << "\t" << beta1[j] << "\t" << E1[j];
			message += reactionRates->strReaction[j];
			message += "\nReaction with negative pre-exponential factor!";
			WarningMessage(message);
			k01[j] *= -1.;
			negativeSigns.Append(j);
		}

		k01[j] = log(k01[j]);
		E1[j] = -E1[j] / Constants::R_cal_mol;

		if(jThirdBody[j] >= 2 && jThirdBody[j] <= 7)
		{
			inputFile  >> k02[j] >> beta2[j] >> E2[j];
			outputFile << k02[j] << beta2[j] << E2[j];
			asciiFile << k02[j] << beta2[j] << E2[j];
			if (inputFile.eof()!=0) ErrorMessage("End of file");
			if(k02[j] <= 0.)
				ErrorMessage("Reaction with negative pre-exponential (third body) factor!");
			k02[j] = log(k02[j]);
			E2[j] = -E2[j] / Constants::R_cal_mol;
			
			if (inputFile.eof()!=0) ErrorMessage("End of file");

		
			inputFile  >> aPressure[j] >> bPressure[j] >> cPressure[j] >> dPressure[j] >> ePressure[j];
			outputFile << aPressure[j] << bPressure[j] << cPressure[j] << dPressure[j] << ePressure[j];
			asciiFile << aPressure[j] << bPressure[j] << cPressure[j] << dPressure[j] << ePressure[j];
			if (inputFile.eof()!=0) ErrorMessage("End of file");
		}

		if(jThirdBody[j] == 10 || jThirdBody[j] == 100)
		{
			double landau_teller;
			inputFile   >> landau_teller;	outputFile  << landau_teller; asciiFile  << landau_teller;
			inputFile   >> landau_teller;	outputFile  << landau_teller; asciiFile  << landau_teller;
				if (inputFile.eof()!=0) ErrorMessage("End of file");
		}
		
		if(jThirdBody[j] == 11 || jThirdBody[j] == 110)
		{
			double janev_langer;
			for(int k=1;k<=9;k++)
			{
				inputFile   >> janev_langer;
				outputFile  << janev_langer;
				asciiFile  << janev_langer;
				if (inputFile.eof()!=0) ErrorMessage("End of file");
			}
		}

		if(jThirdBody[j] == 12 || jThirdBody[j] == 120)
		{
			double power_series;
			for(int k=1;k<=4;k++)
			{
				inputFile   >> power_series;
				outputFile  << power_series;
				asciiFile  << power_series;
				if (inputFile.eof()!=0) ErrorMessage("End of file");
			}
		}

		if(jThirdBody[j] == 160)
		{
			double TAR_series;
			for(int k=1;k<=6;k++)
			{
				inputFile   >> TAR_series;
				outputFile  << TAR_series;
				asciiFile  << TAR_series;
				if (inputFile.eof()!=0) ErrorMessage("End of file");
			}
		}

		if(jThirdBody[j] == 13 || jThirdBody[j] == 130)
		{
			double chebishev;
			int N, M;
			inputFile   >> N >> M;
			outputFile  << N << M;
			asciiFile   << N << M;
			for(int k=1;k<=N*M;k++)
			{
				inputFile   >> chebishev;
				outputFile  << chebishev;
				asciiFile  << chebishev;
				if (inputFile.eof()!=0) ErrorMessage("End of file");
			}
			inputFile   >> chebishev;	outputFile << chebishev;	asciiFile << chebishev;// conversion
			inputFile   >> chebishev;	outputFile << chebishev;	asciiFile << chebishev;// Tmin
			inputFile   >> chebishev;	outputFile << chebishev;	asciiFile << chebishev;// Tmax
			inputFile   >> chebishev;	outputFile << chebishev;	asciiFile << chebishev;// Pmin
			inputFile   >> chebishev;	outputFile << chebishev;	asciiFile << chebishev;// Pmax
		}

		if(jThirdBody[j] == 14 || jThirdBody[j] == 140)
		{
			double log_pressure;
			int N;
			inputFile  >> N;
			outputFile << N;
			asciiFile << N;
			for(int k=1;k<=N*4;k++)
			{
				inputFile   >> log_pressure;
				outputFile  << log_pressure;
				asciiFile  << log_pressure;
				if (inputFile.eof()!=0) ErrorMessage("End of file");
			}
		}

		if(jThirdBody[j] == 15 || jThirdBody[j] == 150)
		{
			double collision_efficiency;
			inputFile   >> collision_efficiency;
			outputFile  << collision_efficiency;
			asciiFile  << collision_efficiency;
			if (inputFile.eof()!=0) ErrorMessage("End of file");
		}
	}

	inputFile.close();
	
	iFallOff.DeleteLastNElements(NUM_MAX-numFallOff);
	iThirdBody.DeleteLastNElements(NUM_MAX-numThirdBody);
	iThirdBodyOnly.DeleteLastNElements(NUM_MAX-numThirdBodyOnly);
	iLandauTeller.DeleteLastNElements(NUM_MAX-numLandauTeller);
	iCollisionEfficiency.DeleteLastNElements(NUM_MAX-numCollisionEfficiency);
	iJanevLanger.DeleteLastNElements(NUM_MAX-numJanevLanger);
	iPowerSeries.DeleteLastNElements(NUM_MAX-numPowerSeries);
	iTAR.DeleteLastNElements(NUM_MAX-numTAR);
	iChebishev.DeleteLastNElements(NUM_MAX-numChebishev);
	iLogarithmicPressure.DeleteLastNElements(NUM_MAX-numLogarithmicPressure);


	numEquilibrium = 0;
	for(j = 1;j <= NR;j++)
	{
		if(jEquil[j] == 1)	numEquilibrium++;
		if(jEquil[j] == 0)	numIrreversible++;
	}
	ChangeDimensions(numEquilibrium,&reactionWithEquil);	
	int neq = 0;
	for(j = 1;j <= NR;j++)
	{
		if(jEquil[j] == 1)
			reactionWithEquil[++neq] = j;
	}

	// -------------------------------------------------------------------------------
	// 4. Lettura del file COEFSTEC.BZZ
	// -------------------------------------------------------------------------------
	readStoichiometricFile(fileSt);

 
	ChangeDimensions(numDir1.Sum(),&jDir1);
	ChangeDimensions(numDir2.Sum(),&jDir2);
	ChangeDimensions(numDir3.Sum(),&jDir3);
	ChangeDimensions(numDir4.Sum(),&jDir4);
	ChangeDimensions(numDir5.Sum(),&jDir5);
	ChangeDimensions(numDir5.Sum(),&valDir5);

	ChangeDimensions(numInvTot1.Sum(),&jInvTot1);
	ChangeDimensions(numInvTot2.Sum(),&jInvTot2);
	ChangeDimensions(numInvTot3.Sum(),&jInvTot3);
	ChangeDimensions(numInvTot4.Sum(),&jInvTot4);
	ChangeDimensions(numInvTot5.Sum(),&jInvTot5);
	ChangeDimensions(numInvTot5.Sum(),&valInvTot5);

	ChangeDimensions(numInvEq1.Sum(),&jInvEq1);
	ChangeDimensions(numInvEq2.Sum(),&jInvEq2);
	ChangeDimensions(numInvEq3.Sum(),&jInvEq3);
	ChangeDimensions(numInvEq4.Sum(),&jInvEq4);
	ChangeDimensions(numInvEq5.Sum(),&jInvEq5);
	ChangeDimensions(numInvEq5.Sum(),&valInvEq5);


	readStoichiometricFile_2(fileSt);

	// 5. Costruzione del vettore somma algebrica dei coefficienti stechiometrici
	//    di reazione
	// -------------------------------------------------------------------------------
	{
		int *jD1 = jDir1.GetHandle();
		int *jD2 = jDir2.GetHandle();
		int *jD3 = jDir3.GetHandle();
		int *jD4 = jDir4.GetHandle();
		int *jD5 = jDir5.GetHandle();
		double *vD5 = valDir5.GetHandle();

		int *jIT1 = jInvTot1.GetHandle();
		int *jIT2 = jInvTot2.GetHandle();
		int *jIT3 = jInvTot3.GetHandle();
		int *jIT4 = jInvTot4.GetHandle();
		int *jIT5 = jInvTot5.GetHandle();
		double *vIT5 = valInvTot5.GetHandle();

		for(i = 1;i <= NC;i++)
		{
			for(k = 1;k <= numDir1[i];k++)
				sumNuij[*jD1++] -= 1.;
			for(k = 1;k <= numDir2[i];k++)
				sumNuij[*jD2++] -= 2.;
			for(k = 1;k <= numDir3[i];k++)
				sumNuij[*jD3++] -= 3.;
			for(k = 1;k <= numDir4[i];k++)
				sumNuij[*jD4++] -= 0.5;
			for(k = 1;k <= numDir5[i];k++)
				sumNuij[*jD5++] -= *vD5++;

			for(k = 1;k <= numInvTot1[i];k++)
				sumNuij[*jIT1++] += 1.;
			for(k = 1;k <= numInvTot2[i];k++)
				sumNuij[*jIT2++] += 2.;
			for(k = 1;k <= numInvTot3[i];k++)
				sumNuij[*jIT3++] += 3.;
			for(k = 1;k <= numInvTot4[i];k++)
				sumNuij[*jIT4++] += 0.5;
			for(k = 1;k <= numInvTot5[i];k++)
				sumNuij[*jIT5++] += *vIT5++;
		}
	}

	// 5. Reaction orders
	// -------------------------------------------------------------------------------
	{
		int *jD1 = jDir1.GetHandle();
		int *jD2 = jDir2.GetHandle();
		int *jD3 = jDir3.GetHandle();
		int *jD4 = jDir4.GetHandle();
		int *jD5 = jDir5.GetHandle();
		double *vD5 = valDir5.GetHandle();

		int *jIT1 = jInvTot1.GetHandle();
		int *jIT2 = jInvTot2.GetHandle();
		int *jIT3 = jInvTot3.GetHandle();
		int *jIT4 = jInvTot4.GetHandle();
		int *jIT5 = jInvTot5.GetHandle();
		double *vIT5 = valInvTot5.GetHandle();

		for(i = 1;i <= NC;i++)
		{
			for(k = 1;k <= numDir1[i];k++)
				forwardOrders[*jD1++] += 1.;
			for(k = 1;k <= numDir2[i];k++)
				forwardOrders[*jD2++] += 2.;
			for(k = 1;k <= numDir3[i];k++)
				forwardOrders[*jD3++] += 3.;
			for(k = 1;k <= numDir4[i];k++)
				forwardOrders[*jD4++] += 0.5;
			for(k = 1;k <= numDir5[i];k++)
				forwardOrders[*jD5++] += *vD5++;

			for(k = 1;k <= numInvTot1[i];k++)
				backwardOrders[*jIT1++] += 1.;
			for(k = 1;k <= numInvTot2[i];k++)
				backwardOrders[*jIT2++] += 2.;
			for(k = 1;k <= numInvTot3[i];k++)
				backwardOrders[*jIT3++] += 3.;
			for(k = 1;k <= numInvTot4[i];k++)
				backwardOrders[*jIT4++] += 0.5;
			for(k = 1;k <= numInvTot5[i];k++)
				backwardOrders[*jIT5++] += *vIT5++;
		}
	}
	
	// 5. Riepilogo delle informazioni a video
	// -------------------------------------------------------------------------------
	cout << "--------------------------------------------------------------------------" << endl;
	cout << "                       KINETIC SCHEME SUMMARY                             " << endl;
	cout << "--------------------------------------------------------------------------" << endl;
	cout << "  Total number of species                   = " << NC			 << endl;
	cout << "  Total number of reactions                 = " << NR			 << endl;
	cout << "    Number of irreversible reactions        = " << numIrreversible			 << endl;
	cout << "    Number of equilibrium reactions         = " << numEquilibrium			 << endl;
	cout << "    Number of conventional reactions        = " << numConventional			 << endl;
	cout << "    Number of third body reactions (total)  = " << numThirdBody		 	 << endl;
	cout << "    Number of third body reactions          = " << numThirdBodyOnly		 << endl;
	cout << "    Number of Fall-Off reactions            = " << numFallOff				 << endl;
	cout << "    Number of C.A.B. reactions              = " << numCABR					 << endl;
	cout << "    Number of Landau-Teller reactions       = " << numLandauTeller			 << endl;
	cout << "    Number of Janev-Langer reactions        = " << numJanevLanger			 << endl;
	cout << "    Number of Power-Series reactions        = " << numPowerSeries			 << endl;
	cout << "    Number of Chebishev-Pol. reactions      = " << numChebishev			 << endl;
	cout << "    Number of Log-Pressure reactions        = " << numLogarithmicPressure	 << endl;
	cout << "    Number of Bimolecular Coll. reactions   = " << numCollisionEfficiency	 << endl;
	cout << "    Number of TAR reactions                 = " << numTAR					 << endl;
	cout << "--------------------------------------------------------------------------" << endl;
	cout << endl;

	// 9. Write Binary File
	// -------------------------------------------------------------------------------
	{
		outputFile << numDir1;
		outputFile << numDir2;
		outputFile << numDir3;
		outputFile << numDir4;
		outputFile << numDir5;

		outputFile << numInvTot1;
		outputFile << numInvTot2;
		outputFile << numInvTot3;
		outputFile << numInvTot4;
		outputFile << numInvTot5;

		outputFile << numInvEq1;
		outputFile << numInvEq2;
		outputFile << numInvEq3;
		outputFile << numInvEq4;
		outputFile << numInvEq5;

		outputFile << jDir1;
		outputFile << jDir2;
		outputFile << jDir3;
		outputFile << jDir4;
		outputFile << jDir5;
		outputFile << valDir5;

		outputFile << jInvTot1;
		outputFile << jInvTot2;
		outputFile << jInvTot3;
		outputFile << jInvTot4;
		outputFile << jInvTot5;
		outputFile << valInvTot5;

		outputFile << jInvEq1;
		outputFile << jInvEq2;
		outputFile << jInvEq3;
		outputFile << jInvEq4;
		outputFile << jInvEq5;
		outputFile << valInvEq5;

		outputFile << sumNuij;
	}

	{
		asciiFile << numDir1;
		asciiFile << numDir2;
		asciiFile << numDir3;
		asciiFile << numDir4;
		asciiFile << numDir5;

		asciiFile << numInvTot1;
		asciiFile << numInvTot2;
		asciiFile << numInvTot3;
		asciiFile << numInvTot4;
		asciiFile << numInvTot5;

		asciiFile << numInvEq1;
		asciiFile << numInvEq2;
		asciiFile << numInvEq3;
		asciiFile << numInvEq4;
		asciiFile << numInvEq5;

		asciiFile << jDir1;
		asciiFile << jDir2;
		asciiFile << jDir3;
		asciiFile << jDir4;
		asciiFile << jDir5;
		asciiFile << valDir5;

		asciiFile << jInvTot1;
		asciiFile << jInvTot2;
		asciiFile << jInvTot3;
		asciiFile << jInvTot4;
		asciiFile << jInvTot5;
		asciiFile << valInvTot5;

		asciiFile << jInvEq1;
		asciiFile << jInvEq2;
		asciiFile << jInvEq3;
		asciiFile << jInvEq4;
		asciiFile << jInvEq5;
		asciiFile << valInvEq5;

		asciiFile << sumNuij;
	}

	// Writing reaction names on file
	for(j = 1;j <= NR;j++)
	{
		char name[Constants::REACTION_NAME_SIZE];
		strcpy(name, reactionRates->strReaction[j]);
		outputFile.fileSave.write((char*) name, sizeof(name));
		asciiFile << name;
	}

	// Writing additional info on stoichiometry (number of global reactions)
	outputFile << forwardOrders;		// forward orders
	outputFile << backwardOrders;		// backward orders
	outputFile << 0;					// no global reactions
	outputFile << 0;					// no soot mode

	asciiFile << forwardOrders;		// forward orders
	asciiFile << backwardOrders;		// backward orders
	asciiFile << 0;					// no global reactions
	asciiFile << 0;					// no soot mode

	outputFile.End();
	asciiFile.End();
};

void OpenSMOKE_PreProcessorKinetics::CheckingStoichiometry()
{
	// Checking stoichiometry
	cout << "  ** Checking Stoichiometry..." << endl;
	{
		int i, k;

		int *jD1 = jDir1.GetHandle();
		int *jD2 = jDir2.GetHandle();
		int *jD3 = jDir3.GetHandle();
		int *jD4 = jDir4.GetHandle();
		int *jD5 = jDir5.GetHandle();
		double *vD5 = valDir5.GetHandle();

		int *jIT1 = jInvTot1.GetHandle();
		int *jIT2 = jInvTot2.GetHandle();
		int *jIT3 = jInvTot3.GetHandle();
		int *jIT4 = jInvTot4.GetHandle();
		int *jIT5 = jInvTot5.GetHandle();
		double *vIT5 = valInvTot5.GetHandle();

		BzzVector sum_stoichiometry_r(NR);
		BzzVector sum_stoichiometry_p(NR);
		for(i = 1;i <= NC;i++)
		{
			for(k = 1;k <= numDir1[i];k++)
				sum_stoichiometry_r[*jD1++] += 1.*reactionRates->M[i];
			for(k = 1;k <= numDir2[i];k++)
				sum_stoichiometry_r[*jD2++] += 2.*reactionRates->M[i];
			for(k = 1;k <= numDir3[i];k++)
				sum_stoichiometry_r[*jD3++] += 3.*reactionRates->M[i];
			for(k = 1;k <= numDir4[i];k++)
				sum_stoichiometry_r[*jD4++] += 0.5*reactionRates->M[i];
			for(k = 1;k <= numDir5[i];k++)
				sum_stoichiometry_r[*jD5++] += (*vD5++)*reactionRates->M[i];;

			for(k = 1;k <= numInvTot1[i];k++)
				sum_stoichiometry_p[*jIT1++] += 1.*reactionRates->M[i];
			for(k = 1;k <= numInvTot2[i];k++)
				sum_stoichiometry_p[*jIT2++] += 2.*reactionRates->M[i];
			for(k = 1;k <= numInvTot3[i];k++)
				sum_stoichiometry_p[*jIT3++] += 3.*reactionRates->M[i];
			for(k = 1;k <= numInvTot4[i];k++)
				sum_stoichiometry_p[*jIT4++] += 0.5*reactionRates->M[i];
			for(k = 1;k <= numInvTot5[i];k++)
				sum_stoichiometry_p[*jIT5++] += (*vIT5++)*reactionRates->M[i];;
		}

		BzzVectorInt indices_wrong;
		BzzVector	 reactants_wrong;
		BzzVector    products_wrong;
		vector<string> names_wrong; names_wrong.push_back("names");

		for(i = 1;i <= NR;i++)
			if ( (sum_stoichiometry_p[i]-sum_stoichiometry_r[i]) >  0.0001 || 
				 (sum_stoichiometry_p[i]-sum_stoichiometry_r[i]) < -0.0001  )
			{
				cout << "Reaction: " << i << endl;
				cout << reactionRates->strReaction[i] << endl;
				cout << " React.:  " << sum_stoichiometry_r[i] << endl;
				cout << " Prod.:   " << sum_stoichiometry_p[i] << endl;
				indices_wrong.Append(i);
				reactants_wrong.Append(sum_stoichiometry_r[i]);
				products_wrong.Append(sum_stoichiometry_p[i]);
				names_wrong.push_back(reactionRates->strReaction[i]);
			}
		
		if (indices_wrong.Size()>0)
		{
			ofstream fWrong;
			openOutputFileAndControl(fWrong, "wrong.out");
			fWrong.setf(ios::scientific);
			for(int j=1;j<=indices_wrong.Size();j++)
				fWrong	<< j << "\t"
						<< indices_wrong[j] << "\t"
						<< reactants_wrong[j] << "\t"
						<< products_wrong[j] << "\t"
						<< names_wrong[j] << "\t"
						<< endl;
		}
	}
}

void OpenSMOKE_PreProcessorKinetics::readStoichiometricFile(const std::string fileSt)
{
	int i,j;
	double val;
	int NC, NR;

	ifstream inputFile;
	openInputFileAndControl(inputFile, fileSt);
	inputFile  >> NC >> NR;
	
	while(1)
	{
		inputFile >> i >> j >> val;
		if (inputFile.eof()!=0) break;

		if(val == -1.)
			numDir1[i]++;
		else if(val == -2.)
			numDir2[i]++;
		else if(val == -3.)
			numDir3[i]++;
		else if(val == -0.50)
			numDir4[i]++;
		else if(val < 0.)
			numDir5[i]++;
		else if(val == 1.)
			{
			numInvTot1[i]++;
			if(jEquil[j] != 0)
				numInvEq1[i]++;
			}
		else if(val == 2.)
			{
			numInvTot2[i]++;
			if(jEquil[j] != 0)
				numInvEq2[i]++;
			}
		else if(val == 3.)
			{
			numInvTot3[i]++;
			if(jEquil[j] != 0)
				numInvEq3[i]++;
			}
		else if(val == 0.50)
			{
			numInvTot4[i]++;
			if(jEquil[j] != 0)
				numInvEq4[i]++;
			}
		else if(val > 0.)
			{
			numInvTot5[i]++;
			if(jEquil[j] != 0)
				numInvEq5[i]++;
			}
	}

	inputFile.close();
}

void OpenSMOKE_PreProcessorKinetics::SparsityStructures(const std::string fileSt)
{
	int i,j;
	double val;
	int numComponents, NR;

	ifstream inputFile;
	openInputFileAndControl(inputFile, fileSt);
	inputFile  >> numComponents >> NR;

	
	BzzMatrix dRoverdC, dRStaroverdC, dRStarStaroverdC, nu;
	ChangeDimensions(numComponents, NR,  &dRoverdC);
	ChangeDimensions(numComponents, NR,  &dRStaroverdC);
	ChangeDimensions(numComponents, NR,  &dRStarStaroverdC);
	ChangeDimensions(NR,  numComponents, &nu);

	while(1)
	{
		inputFile >> i >> j >> val;
		if (inputFile.eof()!=0) break;

		// Matrice della struttura di sparsita
		if (val>0.) dRoverdC[i][j] = 1;			
		else if (val < 0. && jEquil[j]==1) dRoverdC[i][j] = 1;
		
		int ke; 
		if (jThirdBody[j]==1)
		{
			for( ke=1;ke<=numComponents;ke++)
				dRStaroverdC[ke][j] = 1;
			for( ke = 1;ke <= jNumEfficiency[j];ke++)
				if (efficiency[j][ke] == -1.) dRStaroverdC[iEfficiency[j][ke]][j] = 0;
		}

		if (jThirdBody[j]>=2 && jThirdBody[j]<=7)
		{
			for( ke=1;ke<=numComponents;ke++)
				dRStarStaroverdC[ke][j] = 1;
			for( ke = 1;ke <= jNumEfficiency[j];ke++)
				if (efficiency[j][ke] == -1.) dRStarStaroverdC[iEfficiency[j][ke]][j] = 0;
		}
		
		nu[j][i] = 1;

		// ------------------------------------
	}
 
	inputFile.close();

	// Strutture di sparsita
	// --------------------------------------------------------------------
	ofstream outputNu("Sparsity/nuStructure.dat", ios::out);
	for(i=1;i<=numComponents;i++)
	{
		for(j=1;j<=NR;j++)
			outputNu << nu[j][i] << "\t";
		outputNu << endl;
	}
	outputNu.close();

	ofstream outputR("Sparsity/RStructure.dat", ios::out);
	for(i=1;i<=numComponents;i++)
	{
		for(j=1;j<=NR;j++)
			outputR << dRoverdC[i][j] << "\t";
		outputR << endl;
	}
	outputR.close();


	ofstream outputRStar("Sparsity/RStarStructure.dat", ios::out);
	for(i=1;i<=numComponents;i++)
	{
		for(j=1;j<=NR;j++)
			outputRStar << dRStaroverdC[i][j] << "\t";
		outputRStar << endl;
	}
	outputRStar.close();

	ofstream outputRStarStar("Sparsity/RStarStarStructure.dat", ios::out);
	for(i=1;i<=numComponents;i++)
	{
		for(j=1;j<=NR;j++)
			outputRStarStar << dRStarStaroverdC[i][j] << "\t";
		outputRStarStar << endl;
	}
	outputRStarStar.close();

	cout << "Sparsity structures correctly written!" << endl;
}

void OpenSMOKE_PreProcessorKinetics::readStoichiometricFile_2(const std::string fileSt)
{
	int i,j;
	double val;

	ifstream inputFile;
	openInputFileAndControl(inputFile, fileSt);
	
	int ro, co;
	inputFile >> ro >> co ;

	int *jD1 = jDir1.GetHandle();
	int *jD2 = jDir2.GetHandle();
	int *jD3 = jDir3.GetHandle();
	int *jD4 = jDir4.GetHandle();
	int *jD5 = jDir5.GetHandle();
	double *vD5 = valDir5.GetHandle();

	int *jIT1 = jInvTot1.GetHandle();
	int *jIT2 = jInvTot2.GetHandle();
	int *jIT3 = jInvTot3.GetHandle();
	int *jIT4 = jInvTot4.GetHandle();
	int *jIT5 = jInvTot5.GetHandle();
	double *vIT5 = valInvTot5.GetHandle();
	
	int *jIE1 = jInvEq1.GetHandle();
	int *jIE2 = jInvEq2.GetHandle();
	int *jIE3 = jInvEq3.GetHandle();
	int *jIE4 = jInvEq4.GetHandle();
	int *jIE5 = jInvEq5.GetHandle();
	double *vIE5 = valInvEq5.GetHandle();

	while(1)
	{
		inputFile >> i >> j >> val;
		if (inputFile.eof()!=0) break;
		
		if(val == -1.)
			*jD1++ = j;
		else if(val == -2.)
			*jD2++ = j;
		else if(val == -3.)
			*jD3++ = j;
		else if(val == -0.50)
			*jD4++ = j;
		else if(val < 0.)
			{
			*jD5++ = j;
			*vD5++ = fabs(val);
			}
		
		else if(val == 1.)
			{
			*jIT1++ = j;
			if(jEquil[j] != 0)
				*jIE1++ = j;
			}
		else if(val == 2.)
			{
			*jIT2++ = j;
			if(jEquil[j] != 0)
				*jIE2++ = j;
			}
		else if(val == 3.)
			{
			*jIT3++ = j;
			if(jEquil[j] != 0)
				*jIE3++ = j;
			}
		else if(val == 0.50)
			{*jIT4++ = j;
			if(jEquil[j] != 0)
				*jIE4++ = j;
			}
		else if(val > 0.)
			{
			*jIT5++ = j;
			*vIT5++ = val;
			if(jEquil[j] != 0)
				{
				*jIE5++ = j;
				*vIE5++ = val;
				}
			}	
	}

	inputFile.close();
}
