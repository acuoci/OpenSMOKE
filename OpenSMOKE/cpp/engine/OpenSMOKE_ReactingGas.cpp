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
#include <iomanip>

void inverse_kinetics_non_linear_regression_complete(int model, int ex, BzzVector &b, BzzVector &x, BzzVector &y);
void inverse_kinetics_non_linear_regression_partial(int model, int ex, BzzVector &b, BzzVector &x, BzzVector &y);

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
	iSootMode = false;
}



void OpenSMOKE_ReactingGas::SetName(const string name)
{
	name_object = name;
}

int	OpenSMOKE_ReactingGas::NumberOfReactions()
{
	return NR;
}

bool OpenSMOKE_ReactingGas::IsTransportModeAvailable()
{
	return iTransportMode;
}

void OpenSMOKE_ReactingGas::SetupBinary(const string pathName)
{
	string fileNameReactions;
	string fileNameIdealGas;
	if (binary_version_ == true)
	{
		cout << "BINARY " << pathName << endl;
		fileNameReactions		= pathName + "/reactions.bin";
		fileNameIdealGas		= pathName + "/idealgas.bin";
		cout << "BINARY " << fileNameReactions << endl;
		cout << "BINARY " << fileNameIdealGas << endl;
	}
	else
	{
		cout << "ASCII " << pathName << endl;
		fileNameReactions		= pathName + "/reactions.oldstyle.ascii";
		fileNameIdealGas		= pathName + "/idealgas.oldstyle.ascii";
		cout << "ASCII " << fileNameReactions << endl;
		cout << "ASCII " << fileNameIdealGas << endl;
	}

	kinetics.reactionRates = this;

	// 1. Kinetics version
	BzzLoad binaryFile;
	if (binary_version_ == true)	binaryFile('*', fileNameIdealGas);
	else							binaryFile(fileNameIdealGas);

	cout << fileNameIdealGas << endl;
	ReadVersion(binaryFile);

	// 2. Lettura del numero di componenti e del numero di reazioni	
	BzzLoad inputFile;
	if (binary_version_ == true)	inputFile('*', fileNameReactions);
	else							inputFile(fileNameReactions);
	inputFile >> NC >> NR;
	inputFile.End();

	// 3. Lettura dello schema cinetico e dei parametri cinetici
	kinetics.readFromFileBinary(fileNameReactions);

	// 4. Lettura delle informazioni termodinamiche	
	if (binary_version_ == true)	binaryFile('*', fileNameIdealGas);
	else							binaryFile(fileNameIdealGas);
	Setup(binaryFile);
	Fitting(binaryFile);
	binaryFile.End();

	// Kinetic Scheme name
	size_t found;
	found = pathName.find_last_of("/\\");
	path_complete		= pathName;
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

void OpenSMOKE_ReactingGas::ComputeKineticParameters(double T, double logT, double uT, double P_Pa)
{
	double loguRT = logCATM - logT;

	SpeciesEnthalpyAndEntropy(T);

	kinetics.ComputeDHandDSreaction(h, s, reactionDH, reactionDS);

	kinetics.ComputeKineticParameters(reactionDH, reactionDS, T, P_Pa, logT, uT, loguRT, uKeq, k1, k2, logFcent);
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
		kinetics.ThirdBody(c, cTot, coeffM);

		// Corregge la Costante di Arrhenius per le reazioni di FallOff
		kinetics.FallOff(T, coeffM, k1, k2, logFcent, coeffFallOff);

		// Corregge la Costante di Arrhenius per le reazioni CABR
		kinetics.ChemicallyActivatedBimolecularReactions(T, coeffM, k1, k2, logFcent, coeffFallOff);
	
	// Parte 2: Assemblaggio dei contributi per la costruzione delle velocita di reazione
	// ----------------------------------------------------------------------------------
		kinetics.reactionsWithEquilibrium(rDirC, rInvC, uKeq, r);
		kinetics.reactionsWithThirdBody(coeffM,r);

		for(int j = 1;j <= NR;j++)
			r[j] *= coeffFallOff[j];

	// Parte 3: Costruzione delle velocita di reazione per ciascun componente (1...NC)
	// -------------------------------------------------------------------------------
		kinetics.compositionReactionRates(r,R);	
}

void OpenSMOKE_ReactingGas::ComputeFromConcentrations(double T, BzzVector &c, double cTot, BzzVector &rForward, BzzVector &rBackward)
{
	// Parte 1: Costruzione dei contributi alla velocita di reazione		
	// -----------------------------------------------------------------------------
		// Calcolo delle velocita di reazione dirette e inverse (in realta si 
		// tratta soltanto della produttoria delle concentrazioni, non ancora 
		// moltiplicata per alcuna costante cinetica)
		kinetics.ComputeDirectAndInverse(c,rForward,rBackward);

		// Calcolo dei coefficienti correttivi per TUTTE le reazioni con terzo corpo
		// (anche quelle di FallOff)
		kinetics.ThirdBody(c, cTot, coeffM);

		// Corregge la Costante di Arrhenius per le reazioni di FallOff
		kinetics.FallOff(T, coeffM, k1, k2, logFcent, coeffFallOff);

		// Corregge la Costante di Arrhenius per le reazioni CABR
		kinetics.ChemicallyActivatedBimolecularReactions(T, coeffM, k1, k2, logFcent, coeffFallOff);
	
	// Parte 2: Assemblaggio dei contributi per la costruzione delle velocita di reazione
	// ----------------------------------------------------------------------------------
		kinetics.reactionsWithEquilibriumForwardAndBackward(rBackward, uKeq);
		kinetics.reactionsWithThirdBodyForwardAndBackward(coeffM,rForward,rBackward);

		for(int j=1;j<=NR;j++)
		{
			rForward[j]  *= coeffFallOff[j];
			rBackward[j] *= coeffFallOff[j];
		}
		
		for(int j=1;j<=NR;j++)
			rBackward[j] *= -1.;
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
	ChangeDimensions(numPoints,  &P_Pa_map);	
	ChangeDimensions(NC,  &RVector);	
}

void OpenSMOKE_ReactingGas::ComputeKineticParameters_map(BzzVector &T, BzzVector &cTot, BzzVector &P_Pa)
{
	int i;

	T_map = T;
	cTot_map = cTot;
	P_Pa_map = P_Pa;
	for(i=1;i<=T.Size();i++)
	{
		ComputeKineticParameters( T[i], log(T[i]), 1./T[i], P_Pa[i]);
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
		
		for(s = 1;s <= kinetics.numCABR;s++)
		{
			j=kinetics.iCABR[s];
			k2_map[k][j] *= correction_k2[k][j];
		}
	}	
}

void OpenSMOKE_ReactingGas::CorrectKineticParameters(BzzVector &correction_k1,BzzVector &correction_uKeq, BzzVector &correction_k2)
{
	for(int j=1;j<=NR;j++)
		k1[j] *= correction_k1[j];

	for(int s=1;s<=kinetics.numEquilibrium;s++)
	{
		int j = kinetics.reactionWithEquil[s];
		uKeq[j] *= correction_uKeq[j];
	}

	for(int s = 1;s <= kinetics.numFallOff;s++)
	{
		int j=kinetics.iFallOff[s];
		k2[j] *= correction_k2[j];
	}

	for(int s=1;s<=kinetics.numCABR;s++)
	{
		int j=kinetics.iCABR[s];
		k2[j] *= correction_k2[j];
	}	
}

// --------------------------------------------------------------------------
//							SENSITIVITY FUNCTIONS
// --------------------------------------------------------------------------
void OpenSMOKE_ReactingGas::GiveMe_Jalfa_A(	BzzMatrix &JAlfa, OpenSMOKE_NuManager *nu, const double T, 
											const int indexSpecies, const int indexTemperature, BzzVector &parameters)
{
	int j,k;		
	// TODO: CABR	reactions
	// TODO: LT		reactions
	
	BzzVector denominator(NR+kinetics.numFallOff);

	// Scaling factor for numerical computation
	for(j=1;j<=NR;j++)
		denominator[j] = kinetics.exp_k01[j];
	for(k=1; k<=kinetics.numFallOff; k++)
		denominator[NR+k] = 0.;

	for(k=1; k<=kinetics.numFallOff; k++)
	{
		int j=kinetics.iFallOff[k];
		double Pr = coeffM[j]*k1[j]/k2[j];

		denominator[j]		= kinetics.exp_k01[j]*(1.+Pr);
		denominator[NR+k]	= kinetics.exp_k02[j]*(1.+Pr)/Pr;

		double Cc1, Cc2;
		if (kinetics.jThirdBody[j] == 2)		// LINDEMANN
		{
			Cc1 = 1.;
			Cc2 = 1.;
		}
		else if (kinetics.jThirdBody[j] == 3)	// TROE
		{
			double logPr = log10(Pr);
			double nTroe = 0.75 - 1.27 * logFcent[j];		// m
			double cTroe = -0.4 - 0.67 * logFcent[j];		// c
			double sTroe = logPr + cTroe;					// csi + c

			double BetaTroe = sTroe/(nTroe-0.14*sTroe);
			double AlfaTroe = logFcent[j]/(1.+BetaTroe*BetaTroe);

			double coeff =  2.*AlfaTroe*BetaTroe/(1.+BetaTroe*BetaTroe) * 
						    nTroe/BzzPow2(0.14*sTroe-nTroe);
			Cc1 = 1. - (1.+Pr)*coeff;
			Cc2 = 1. + (1.+Pr)/Pr*coeff;
		}
		else if (kinetics.jThirdBody[j] == 4)	// SRI
		{
			double logPr = log10(Pr);
			double xSRI = 1. / (1. + BzzPow2(logPr));

			double coeff = log(logFcent[j])/log(10.)*2.*logPr*xSRI/(1.+logPr*logPr);
			Cc1 = 1. - (1+Pr)*coeff;
			Cc2 = 1. + (1+Pr)/Pr*coeff;
		}
		else	
		{
			Cc1 = 1.;
			Cc2 = 1.;
		}

		if (Cc1 == 0.)	ErrorMessage("Correction coefficient Cc1 equal to zero");
		if (Cc2 == 0.)	ErrorMessage("Correction coefficient Cc2 equal to zero");
		denominator[j]		/= Cc1;
		denominator[NR+k]	/= Cc2;
	}

	// Jalfa for species
	for(int i=1;i<=NC;i++)
		for (int j=1;j<=nu[i].nReactions;j++)
		{
			int jReaction = nu[i].iReactions[j];
			JAlfa[indexSpecies+i-1][jReaction] = nu[i].nuReactions[j]*r[jReaction] * M[i] / denominator[jReaction];

			int kFallOff = kinetics.IsAFallOffReaction(jReaction);
			if (kFallOff != 0)
				JAlfa[indexSpecies+i-1][NR+kFallOff] = nu[i].nuReactions[j]*r[jReaction] * M[i] / denominator[NR+kFallOff];
		}

	// Jalfa for temperature
	if (indexTemperature > 0)
	{
		double coefficient = -Constants::R_J_kmol*T;
		for(j=1;j<=NR;j++)
			JAlfa[indexTemperature][j] = r[j]*reactionDH[j] / denominator[j] * coefficient;
		for(k=1; k<=kinetics.numFallOff; k++)
		{
			int j=kinetics.iFallOff[k];
			JAlfa[indexTemperature][NR+k] = r[j]*reactionDH[j] / denominator[NR+k] * coefficient;
		}
	}

	// Parameters (Pre-exponential factors)
	for(j=1;j<=NR;j++)
		parameters[j] = kinetics.exp_k01[j];
	for(k=1; k<=kinetics.numFallOff; k++)
	{
		int j=kinetics.iFallOff[k];
		parameters[NR+k] = kinetics.exp_k02[j];
	}
}

void OpenSMOKE_ReactingGas::GiveMe_Jalfa_Beta(	BzzMatrix &JAlfa, OpenSMOKE_NuManager *nu, const double T, 
												const int indexSpecies, const int indexTemperature, BzzVector &parameters)
{
}

void OpenSMOKE_ReactingGas::GiveMe_Jalfa_Eatt(	BzzMatrix &JAlfa, OpenSMOKE_NuManager *nu, const double T, 
												const int indexSpecies, const int indexTemperature, BzzVector &parameters)
{
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

string MirrorString(const string s)
{
	string reactants;
	string products;
	int count = 0;
	for (int i=0;i<s.size();i++)
	{
		count++;
		if (s.at(i)!='=')
		{
			reactants += s.at(i);
		}
		else
			break;
	}
	for (int i=count;i<s.size();i++)
		products += s.at(i);

	string reaction = products + "=>" + reactants;

	return reaction;
}

void OpenSMOKE_ReactingGas::VerboseDataKinetics(ofstream &fOutput, ofstream &fOutputInverseKinetics)
{
	int nIntervals = int((TMAX-TMIN)/100.);

	BzzVector T_Vector(nIntervals);
	BzzMatrix kappa_Matrix(nIntervals, NR);
	BzzMatrix uKappaEquilibrium_Matrix(nIntervals, NR);
	BzzMatrix DH_Matrix(nIntervals,NR);
	BzzMatrix DS_Matrix(nIntervals,NR);

	BzzVector DH_Tref(NR);
	BzzVector DS_Tref(NR);

	ComputeKineticParameters(Constants::T_Reference, log(Constants::T_Reference), 1./Constants::T_Reference, 101325.);

	DH_Tref = reactionDH;
	DS_Tref = reactionDS;

	for(int i=1;i<=nIntervals;i++)
	{
		T_Vector[i] = TMIN+(i-1)*100.;

		ComputeKineticParameters(T_Vector[i], log(T_Vector[i]), 1./T_Vector[i], 101325.);

		kappa_Matrix.SetRow(i, k1);
		uKappaEquilibrium_Matrix.SetRow(i, uKeq);

		DH_Matrix.SetRow(i, reactionDH);
		DS_Matrix.SetRow(i, reactionDS);
	}

	VerboseDataSummary(fOutput, DH_Tref, DS_Tref);

	for(int k=1;k<=NR;k++)
		VerboseDataKinetics(fOutput, k, T_Vector, kappa_Matrix, uKappaEquilibrium_Matrix, DH_Matrix, DS_Matrix);

	// Inverse kinetics
	{
		int N = T_Vector.Size();
		BzzMatrix F(N,3);
		BzzVector Y(N);
		for(int i=1;i<=N;i++)
		{
			F[i][1] = 1.;
			F[i][2] = log(T_Vector[i]);
			F[i][3] = -1./T_Vector[i];
		}

		fOutput << endl;
		fOutput << " --------------------------------------------------------------------------------------------------------------------------------" << endl;
		fOutput << "    Reaction        A           Beta           E         Reaction                                                                        " << endl; 
        fOutput << "     index     [kmol,m3,s]                 [cal/mol]                                                                                     " << endl;
		fOutput << " --------------------------------------------------------------------------------------------------------------------------------" << endl;

		for(int k=1;k<=NR;k++)
		{
			if (kinetics.IsAReversibleReaction(k) != 0)
			{
				for(int i=1;i<=N;i++)
					Y[i] = log(kappa_Matrix[i][k]*uKappaEquilibrium_Matrix[i][k]);

				// Arrhenius modified (A T^Beta exp(-Tatt/T))
				{
					BzzLinearRegression linReg(F,Y);
					BzzVector b = linReg.GetParameters();

					fOutput << setw(8)  << right << k;
					fOutput << setw(18) << right << scientific << setprecision(4) << exp(b[1]);
					fOutput << setw(10) << right << fixed << setprecision(3) << b[2];
					fOutput << setw(16) << right << fixed << setprecision(2) << b[3]*Constants::R_cal_mol;
					fOutput << setw(5)  << "";
					fOutput << setw(60) << left << strReaction[k];
					fOutput << endl;
				}
			}
		}

		// TODO: Veronica
		for(int k=1;k<=NR;k++)
		{
			if (kinetics.IsAReversibleReaction(k) != 0)
			{
				for(int i=1;i<=N;i++)
					Y[i] = log(kappa_Matrix[i][k]*uKappaEquilibrium_Matrix[i][k]);

				// Arrhenius modified (A T^Beta exp(-Tatt/T))
				{
					BzzLinearRegression linReg(F,Y);
					BzzVector b = linReg.GetParameters();

					fOutputInverseKinetics << setw(60) << left << MirrorString(strReaction[k]);
					fOutputInverseKinetics << setw(18) << right << scientific << setprecision(4) << exp(b[1]);
					fOutputInverseKinetics << setw(10) << right << fixed << setprecision(3) << b[2];
					fOutputInverseKinetics << setw(16) << right << fixed << setprecision(2) << b[3]*Constants::R_cal_mol;
					fOutputInverseKinetics << endl;
				}
			}
		}
	}

	// Inverse kinetics
	{
		int N = T_Vector.Size();
		BzzMatrix F(N,2);
		BzzVector Y(N);
		for(int i=1;i<=N;i++)
		{
			F[i][1] = 1.;
			F[i][2] = -1./T_Vector[i];
		}

		fOutput << endl;
		fOutput << " --------------------------------------------------------------------------------------------------------------------------------" << endl;
		fOutput << "    Reaction        A                E         Reaction                                                                        " << endl; 
        fOutput << "     index     [kmol,m3,s]       [cal/mol]                                                                                     " << endl;
		fOutput << " --------------------------------------------------------------------------------------------------------------------------------" << endl;

		for(int k=1;k<=NR;k++)
		{
			if (kinetics.IsAReversibleReaction(k) != 0)
			{
				for(int i=1;i<=N;i++)
					Y[i] = log(kappa_Matrix[i][k]*uKappaEquilibrium_Matrix[i][k]);

				// Arrhenius (A exp(-Tatt/T))
				{
					BzzLinearRegression linReg(F,Y);
					BzzVector b = linReg.GetParameters();

					fOutput << setw(8)  << right << k;
					fOutput << setw(18) << right << scientific << setprecision(4) << exp(b[1]);
					fOutput << setw(16) << right << fixed << setprecision(2) << b[2]*Constants::R_cal_mol;
					fOutput << setw(5)  << "";
					fOutput << setw(60) << left << strReaction[k];
					fOutput << endl;
				}
			}
		}
	}
}

void OpenSMOKE_ReactingGas::VerboseDataSummary(ofstream &fOutput, BzzVector &DH_Tref, BzzVector &DS_Tref)
{
	double conversion  = Constants::T_Reference*Constants::R_cal_mol/1000.;

	fOutput << endl;
	fOutput << " ================================================================================================================================" << endl;
	fOutput << "   KINETIC DATA SUMMARY  " << endl;
	fOutput << " ================================================================================================================================" << endl;
	fOutput << "   T = " << Constants::T_Reference << " K" << endl;
	fOutput << " ================================================================================================================================" << endl;
	fOutput << "        Reaction			                                      DG            DH           DS            " << endl;          
    fOutput << "					                                              [kcal/mol]    [kcal/mol]   [cal/mol/K]   " << endl;        
	fOutput << " ================================================================================================================================" << endl;

	for(int k=1;k<=NR;k++)
	{
		fOutput << " " << setw(6) << left << k << setw(50) << strReaction[k];
		fOutput << setw(14) << left << setprecision(4) << fixed << -(DS_Tref[k]-DH_Tref[k])*conversion;	// [kcal/mol]
		fOutput << setw(14) << left << setprecision(4) << fixed << DH_Tref[k]*conversion;						// [kcal/mol]
		fOutput << setw(14) << left << setprecision(4) << fixed << DS_Tref[k]*Constants::R_cal_mol;			// [cal/mol/K]									// [cal/mol/K]
		fOutput << endl;
	}

	fOutput << " ================================================================================================================================" << endl;
	fOutput << endl;
	fOutput << endl;
}

void OpenSMOKE_ReactingGas::VerboseDataKinetics(ofstream &fOutput, const int k, BzzVector &T_Vector, BzzMatrix &kappa_Matrix, BzzMatrix &uKappaEquilibrium_Matrix, BzzMatrix &DH_Matrix, BzzMatrix &DS_Matrix)
{
	double conversion_forward  = pow(1.e3, kinetics.forwardOrders[k]-1.);
	double conversion_backward = pow(1.e3, kinetics.backwardOrders[k]-1.);

	fOutput << " ================================================================================================================================" << endl;
	fOutput << "   KINETIC DATA - REACTION  " << k << endl;
	fOutput << "     " << strReaction[k]	<< endl;
	fOutput << " ================================================================================================================================" << endl;
	kinetics.KineticExpressionString(fOutput, k);

	fOutput << " --------------------------------------------------------------------------------------------------------------------------------" << endl;
	fOutput << "    Temperature         kF           Keq           kR            DG           DH           DS             kF            kR       " << endl;          
    fOutput << "       [K]          [kmol,m3,s]      [-]       [kmol,m3,s]    [kcal/mol]  [kcal/mol]   [cal/mol/K]    [mol,cm3,s]   [mol,cm3,s]  " << endl;        
	fOutput << " --------------------------------------------------------------------------------------------------------------------------------" << endl;

	for(int i=1;i<=T_Vector.Size();i++)
	{
		fOutput << "     " << setw(14) << left << setprecision(4) << fixed << T_Vector[i];
		fOutput            << setw(14) << left << setprecision(4) << scientific << kappa_Matrix[i][k];
		if (kinetics.IsAReversibleReaction(k) != 0)
		{
			fOutput << setw(14) << left << setprecision(4) << scientific << 1./uKappaEquilibrium_Matrix[i][k];
			fOutput << setw(14) << left << setprecision(4) << scientific << kappa_Matrix[i][k]*uKappaEquilibrium_Matrix[i][k];
		}
		else
		{
			fOutput << setw(14) << left << setprecision(4) << scientific << 0.;
			fOutput << setw(14) << left << setprecision(4) << scientific << 0.;
		}
		fOutput            << setw(14) << left << setprecision(4) << fixed << -(DS_Matrix[i][k]-DH_Matrix[i][k])*T_Vector[i]*Constants::R_cal_mol/1000.;	// [kcal/mol]
		fOutput            << setw(14) << left << setprecision(4) << fixed << DH_Matrix[i][k]*T_Vector[i]*Constants::R_cal_mol/1000.;					// [kcal/mol]
		fOutput            << setw(14) << left << setprecision(4) << fixed << DS_Matrix[i][k]*Constants::R_cal_mol;										// [cal/mol/K]

		fOutput            << setw(14) << left << setprecision(4) << scientific << kappa_Matrix[i][k]*conversion_forward;
		fOutput			   << setw(14) << left << setprecision(4) << scientific << kappa_Matrix[i][k]*uKappaEquilibrium_Matrix[i][k]*conversion_backward;	
		fOutput			   << endl;
	}
	fOutput << " --------------------------------------------------------------------------------------------------------------------------------" << endl;
	fOutput << endl;
}

void OpenSMOKE_ReactingGas::SaveOnBinaryFile(BzzSave &fSave)
{
	int j;
	string dummy;
	char name[Constants::NAME_SIZE];
	char name_reaction[Constants::REACTION_NAME_SIZE];

	dummy = name_kinetic_scheme;
	strcpy(name, dummy.c_str());
	fSave.fileSave.write((char*) name, sizeof(name));

	dummy = "NC";
	strcpy(name, dummy.c_str());
	fSave.fileSave.write((char*) name, sizeof(name));
	fSave << NumberOfSpecies();

	dummy = "NR";
	strcpy(name, dummy.c_str());
	fSave.fileSave.write((char*) name, sizeof(name));
	fSave << NumberOfReactions();

	dummy = "M";
	strcpy(name, dummy.c_str());
	fSave.fileSave.write((char*) name, sizeof(name));
	fSave << M;

	dummy = "SPECIES";
	strcpy(name, dummy.c_str());
	fSave.fileSave.write((char*) name, sizeof(name));
	for(j=1;j<=NC;j++)
	{
		char name[Constants::NAME_SIZE];
		strcpy(name, names[j].c_str());
		fSave.fileSave.write((char*) name, sizeof(name));
	}

	dummy = "REACTIONS";
	strcpy(name, dummy.c_str());
	fSave.fileSave.write((char*) name, sizeof(name));
	for(j=1;j<=NR;j++)
	{
		strcpy(name_reaction, kinetics.reactionRates->strReaction[j]);
		fSave.fileSave.write((char*) name_reaction, sizeof(name_reaction));
	}
}


void OpenSMOKE_ReactingGas::ComputeFromConcentrationsAndCorrect(double T, BzzVector &c, double cTot, BzzVector *R, BzzVector &correction)
{
	// Parte 1: Costruzione dei contributi alla velocita di reazione		
	// -----------------------------------------------------------------------------
		// Calcolo delle velocita di reazione dirette e inverse (in realta si 
		// tratta soltanto della produttoria delle concentrazioni, non ancora 
		// moltiplicata per alcuna costante cinetica)
		kinetics.ComputeDirectAndInverse(c,rDirC,rInvC);

		// Calcolo dei coefficienti correttivi per TUTTE le reazioni con terzo corpo
		// (anche quelle di FallOff)
		kinetics.ThirdBody(c, cTot, coeffM);

		// Corregge la Costante di Arrhenius per le reazioni di FallOff
		kinetics.FallOff(T, coeffM, k1, k2, logFcent, coeffFallOff);

		// Corregge la Costante di Arrhenius per le reazioni CABR
		kinetics.ChemicallyActivatedBimolecularReactions(T, coeffM, k1, k2, logFcent, coeffFallOff);
	
	// Parte 2: Assemblaggio dei contributi per la costruzione delle velocita di reazione
	// ----------------------------------------------------------------------------------
		kinetics.reactionsWithEquilibrium(rDirC, rInvC, uKeq, r);
		kinetics.reactionsWithThirdBody(coeffM,r);

		for(int j = 1;j <= NR;j++)
			r[j] *= coeffFallOff[j];

	
		for(int j = 1;j <= NR;j++)
			r[j] *= correction[j];

	// Parte 3: Costruzione delle velocita di reazione per ciascun componente (1...NC)
	// -------------------------------------------------------------------------------
		kinetics.compositionReactionRates(r,R);	
}

void OpenSMOKE_ReactingGas::WriteElementsFile(const string fileName)
{
	ofstream fOutput;
	openOutputFileAndControl(fOutput, fileName);

	fOutput << "ELEMENTS" << endl;
	for(int i=1;i<=NumberOfElements();i++)
	{
		fOutput << setw(12) << left << StringToLower(list_of_elements[i-1]);
		fOutput << setprecision(16) << m_elements[i] << endl;
	}
	fOutput << "END" << endl;


	for(int j=1;j<=NC;j++)
	{
		fOutput << setw(20) << left << names[j];
		for(int i=1;i<=NumberOfElements();i++)
			if (elements[i][j] != 0)
				fOutput << StringToLower(list_of_elements[i-1]) << " " << elements[i][j] << " ";
		fOutput << "/" << endl;
	}
	fOutput << "END" << endl;

	fOutput.close();
}

void inverse_kinetics_non_linear_regression_complete(int model, int ex, BzzVector &b, BzzVector &x, BzzVector &y)
{
	y[1] = b[1] +b[2]*log(x[1]) -b[3]/x[1];
}

void inverse_kinetics_non_linear_regression_partial(int model, int ex, BzzVector &b, BzzVector &x, BzzVector &y)
{
	y[1] = b[1]-b[2]/x[1];
}
