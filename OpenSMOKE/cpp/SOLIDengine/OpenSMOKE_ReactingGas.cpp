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

#include "engine/OpenSMOKE_ReactingGas.h"
#include "basic/OpenSMOKE_Constants.h"

void OpenSMOKE_ReactingGas::ErrorMessage(const string message)
{
    cout << endl;
    cout << "Class:  OpenSMOKE_ReactingGas"		<< endl;
    cout << "Object: " << name_object			<< endl;
    cout << "Error:  " << message				<< endl;
    cout << "Press enter to continue... "		<< endl;
    getchar();
    exit(-1);
}

void OpenSMOKE_ReactingGas::WarningMessage(const string message)
{
    cout << endl;
    cout << "Class:	  OpenSMOKE_ReactingGas"	<< endl;
    cout << "Object:  " << name_object			<< endl;
    cout << "Warning: " << message				<< endl;
	cout << endl;
}

OpenSMOKE_ReactingGas::OpenSMOKE_ReactingGas()
{
	name_object	= "[Name not assigned]";
}

void OpenSMOKE_ReactingGas::SetName(const string name)
{
	name_object = name;
}

int	OpenSMOKE_ReactingGas::NumberOfReactions()
{
	return NR;
}

void OpenSMOKE_ReactingGas::SetupBinary(const string pathName)
{
	string fileNameReactions		= pathName + "/reactions.bin";
	string fileNameThermodynamics	= pathName + "/thermodynamics.bin";
	string fileNameFitting			= pathName + "/fitting.bin";

	kinetics.reactionRates = this;

	// 2. Lettura del numero di componenti e del numero di reazioni	
	BzzLoad inputFile;
	inputFile('*', fileNameReactions);
	inputFile >> NC >> NR;
	inputFile.End();

	// 3. Lettura dello schema cinetico e dei parametri cinetici
	kinetics.readFromFileBinary(fileNameReactions);

	// 4. Lettura delle informazioni termodinamiche	
	BzzLoad binaryFile;
	binaryFile('*', pathName + "/idealgas.bin");
	Setup(binaryFile);
	Fitting(fileNameFitting, binaryFile);
	binaryFile.End();

	// Kinetic Scheme name
	size_t found;
	found = pathName.find_last_of("/\\");
	folder_path			= pathName.substr(0,found);
	name_kinetic_scheme	= pathName.substr(found+1);

	// 5. Dimensionamento dei vettori per il calcolo delle Reactions Rates
	ChangeDimensions(NR,&reactionDH);
	ChangeDimensions(NR,&reactionDS);
	ChangeDimensions(NR,&r);
	ChangeDimensions(NR,&rDirC);
	ChangeDimensions(NR,&rInvC);
	ChangeDimensions(NR,&coeffM);
	ChangeDimensions(NR,&coeffFallOff);
	ChangeDimensions(NR,&uKeq);
	ChangeDimensions(NR,&k1);
	ChangeDimensions(NR,&k2);
	ChangeDimensions(NR,&logFcent);
}

// --------------------------------------------------------------------------
//								MAIN FUNCTIONS
// --------------------------------------------------------------------------

void OpenSMOKE_ReactingGas::ComputeKineticParameters(double T, double logT, double uT)
{
	double loguRT = logCATM - logT;

	SpeciesEnthalpyAndEntropy(T);

	kinetics.ComputeDHandDSreaction(h, s, reactionDH, reactionDS);

	kinetics.ComputeKineticParameters(reactionDH, reactionDS, T, logT, uT, loguRT, uKeq, k1, k2, logFcent);
}

void OpenSMOKE_ReactingGas::ComputeFromConcentrations(double T, BzzVector &c, double cTot, BzzVector *R)
{
	// Parte 1: Costruzione dei contributi alla velocita di reazione		
	// -----------------------------------------------------------------------------
		// Calcolo delle velocita di reazione dirette e inverse (in realta si 
		// tratta soltanto della produttoria delle concentrazioni, non ancora 
		// moltiplicata per alcuna costante cinetica)
		kinetics.ComputeDirectAndInverse(c,rDirC,rInvC);

		// Calcolo dei coefficienti correttivi per TUTTE le reazioni con terzo corpo
		// (anche quelle di FallOff)
		kinetics.thirdBody(c, cTot, coeffM);

		// Corregge la Costante di Arrhenius per le reazioni di FallOff
		kinetics.fallOff(T, coeffM, k1, k2, logFcent, coeffFallOff);
		
	// Parte 2: Assemblaggio dei contributi per la costruzione delle velocita di reazione
	// ----------------------------------------------------------------------------------
		kinetics.reactionsWithEquilibrium(rDirC, rInvC, uKeq, r);
		kinetics.reactionsWithThirdBody(coeffM,r);
		//reactionsWithFallOff(coeffFallOff, r);

		for(int j = 1;j <= NR;j++)
			r[j] *= coeffFallOff[j];

	// Parte 3: Costruzione delle velocita di reazione per ciascun componente (1...NC)
	// -------------------------------------------------------------------------------
		kinetics.compositionReactionRates(r,R);	
}

double OpenSMOKE_ReactingGas::ComputeQReaction(double T)
{
	double QReaction;
	QReaction = Dot(r, reactionDH);				// [kmol/m3s] [-]
	QReaction *= (Constants::R_J_kmol * T);		// [J/m3.s]
	
	return QReaction;
}

// --------------------------------------------------------------------------
//								MAPS
// --------------------------------------------------------------------------
void OpenSMOKE_ReactingGas::InitializeMap(int numPoints)
{
	ChangeDimensions(numPoints, NR, &k1_map);	
	ChangeDimensions(numPoints, NR, &k2_map);	
	ChangeDimensions(numPoints, NR, &uKeq_map);	
	ChangeDimensions(numPoints, NR, &logFcent_map);	
	ChangeDimensions(numPoints, NR, &reactionDH_map);	
	ChangeDimensions(numPoints, NR, &reactionDS_map);	

	ChangeDimensions(numPoints,  &T_map);	
	ChangeDimensions(numPoints,  &cTot_map);	
	ChangeDimensions(NC,  &RVector);	
}

void OpenSMOKE_ReactingGas::ComputeKineticParameters_map(BzzVector &T, BzzVector &cTot)
{
	int i;

	T_map = T;
	cTot_map = cTot;

	for(i=1;i<=T.Size();i++)
	{
		ComputeKineticParameters( T[i], log(T[i]), 1./T[i]);
		k1_map.SetRow(i,k1);
		k2_map.SetRow(i,k2);
		uKeq_map.SetRow(i,uKeq);
		logFcent_map.SetRow(i,logFcent);
		reactionDH_map.SetRow(i,reactionDH);
		reactionDS_map.SetRow(i,reactionDS);
	}
}

void OpenSMOKE_ReactingGas::ComputeFromConcentrations_map(int i, BzzVector &c, BzzVector *RVector)
{
	// Recupero Dati dalle tabelle di lookup
	k1_map.GetRow(i,&k1);
	k2_map.GetRow(i,&k2);
	uKeq_map.GetRow(i,&uKeq);
	logFcent_map.GetRow(i,&logFcent);
	reactionDH_map.GetRow(i,&reactionDH);
	reactionDS_map.GetRow(i,&reactionDS);
	ComputeFromConcentrations( T_map[i], c, cTot_map[i], RVector);// [kmol/m3/s]
//	ElementByElementProduct(RVector, mix.M, &RVector);
//	R.SetRow(i, RVector);			// [kg/m3/s]
}

// Creazione della tabella per le velocita di reazione
void OpenSMOKE_ReactingGas::ComputeFromConcentrations_map(BzzVector &c, BzzMatrix &R, BzzVector &QReaction)
{
	int i;
	for(i=1;i<=T_map.Size();i++)
	{
		ComputeFromConcentrations_map(i, c, &RVector);
		ElementByElementProduct(RVector, M, &RVector);
		R.SetRow(i, RVector);							// [kg/m3/s]
		QReaction[i] = - ComputeQReaction(T_map[i]);	// [J/m3.s]	
	}
}

void OpenSMOKE_ReactingGas::CorrectKineticParameters_map(BzzMatrix &correction_k1, BzzMatrix &correction_uKeq, BzzMatrix &correction_k2)
{
	int j,k,s;
	for (k=1;k<=T_map.Size();k++)
	{
		for(j = 1;j <= NR;j++)
			k1_map[k][j] *= correction_k1[k][j];

		for(s = 1; s <= kinetics.numEquilibrium; s++)
		{
			j = kinetics.reactionWithEquil[s];
			uKeq_map[k][j] *= correction_uKeq[k][j];
		}
		
		for(s = 1;s <= kinetics.numFallOff;s++)
		{
			j=kinetics.iFallOff[s];
			k2_map[k][j] *= correction_k2[k][j];
		}	
	}	
}

// --------------------------------------------------------------------------
//							SENSITIVITY FUNCTIONS
// --------------------------------------------------------------------------

void OpenSMOKE_ReactingGas::JalfaFromConcentrations(char kind, int index, double T, BzzVector &c, double cTot, BzzVector *Jalfa)
{
	int j, k;
	// kind = 'A': only pre-exponential factor
	// kind = 'k': the whole kinetic constant

	// Parte 1: Costruzione dei contributi alla velocita di reazione		
	// -----------------------------------------------------------------------------
		// Calcolo delle velocita di reazione dirette e inverse (in realta si 
		// tratta soltanto della produttoria delle concentrazioni, non ancora 
		// moltiplicata per alcuna costante cinetica)
		kinetics.ComputeDirectAndInverse(c,rDirC,rInvC);

		// Calcolo dei coefficienti correttivi per TUTTE le reazioni con terzo corpo
		// (anche quelle di FallOff)
		kinetics.thirdBody(c, cTot, coeffM);

		// Corregge la Costante di Arrhenius per le reazioni di FallOff
		kinetics.fallOff(T, coeffM, k1, k2, logFcent, coeffFallOff);
		

	// Parte 2: Assemblaggio dei contributi per la costruzione delle velocita di reazione
	// ----------------------------------------------------------------------------------
		kinetics.reactionsWithEquilibrium(rDirC, rInvC, uKeq, r);
		kinetics.reactionsWithThirdBody(coeffM,r);
		ElementByElementProduct(r, coeffFallOff, &r);


	// Parte 2 bis: Costruzione delle derivate rispetto alla costante cinetica
	// ----------------------------------------------------------------------------------
		BzzVector denominator(NR);
		BzzVector alfaParameter(NR);
		if (kind == 'A')	alfaParameter = kinetics.exp_k01;
		if (kind == 'k')	alfaParameter = k1;
		denominator = alfaParameter;

		for(k=1; k<=kinetics.numFallOff; k++)
		{
			j=kinetics.iFallOff[k];
			denominator[j] = alfaParameter[j]/k2[j]*(k2[j]+coeffM[j]*k1[j]);
		}

		for(j=1;j<=NR;j++)
			if (j==index)	r[j] = r[j] / denominator[j];
			else			r[j] = 0.;


	// Parte 3: Costruzione delle velocita di reazione per ciascun componente (1...NC)
	// -------------------------------------------------------------------------------
		kinetics.compositionReactionRates(r, Jalfa);	
}

void OpenSMOKE_ReactingGas::GiveMe_Jalfa_Species(BzzMatrix &JAlfa, OpenSMOKE_NuManager *nu, const int indexSpecies)
{
	BzzVector denominator;

	// Scaling factor for numerical computation
	denominator = kinetics.exp_k01;
	for(int k=1; k<=kinetics.numFallOff; k++)
	{
		int j=kinetics.iFallOff[k];
//		denominator[j] = kinetics.exp_k01[j]/k2[j]*(k2[j]+coeffM[j]*k1[j]);
	}

	// Jalfa for species
	for(int i=indexSpecies;i<=indexSpecies+NC-1;i++)
		for (int j=1;j<=nu[i].nReactions;j++)
		{
			int jReaction = nu[i].iReactions[j];
			JAlfa[i][jReaction] = nu[i].nuReactions[j]*r[jReaction] * M[i] / denominator[jReaction];
		}
}

void OpenSMOKE_ReactingGas::GiveMe_Jalfa_Temperature(BzzMatrix &JAlfa, const double T, const int indexTemperature)
{
	double coefficient = -Constants::R_J_kmol*T;
	for(int i=1;i<=NR;i++)
		JAlfa[indexTemperature][i] = r[i]*reactionDH[i] / kinetics.exp_k01[i] * coefficient;
}

void OpenSMOKE_ReactingGas::GiveMeIndexOfSpeciesInReactions(const string fileName, BzzVectorIntArray &indices)
{
	char comment[Constants::COMMENT_SIZE];
	string pathName;

	ifstream fInput;
	openInputFileAndControl(fInput, fileName);
	fInput >> pathName;				fInput.getline(comment, Constants::COMMENT_SIZE);
	fInput.close();

	string fileNameStoichiometry = pathName + "/stoichiometry.bzz";

	indices(NR);
	kinetics.GiveMeIndexOfSpeciesInEachReaction(fileNameStoichiometry, indices);
}

// --------------------------------------------------------------------------
//								GLOBAL REACTIONS
// --------------------------------------------------------------------------

void OpenSMOKE_ReactingGas::ComputeKEq(const double T, BzzMatrix &nu, BzzVector &sumNuij,BzzVector &uKeq)
{
	int i;
	double logT   = log(T);
	double loguRT = logCATM - logT;
	BzzVector DHReaction(nu.Rows());
	BzzVector DSReaction(nu.Rows());
 
	SpeciesEnthalpyAndEntropy(T);

	DHReaction = 0.;
	DSReaction  = 0.;
	for(i=1;i<=nu.Rows();i++)
		for(int j=1;j<=nu.Columns();j++)
		{
			DHReaction[i] += nu[i][j]*h[j];
			DSReaction[i] += nu[i][j]*s[j];
		}

	for(i=1;i<=nu.Rows();i++)
		uKeq[i] = exp(-DSReaction[i]+DHReaction[i]-loguRT*sumNuij[i]);
}

