/***************************************************************************
 *   Copyright (C) 2007 by Alberto Cuoci   	                               *
 *   alberto.cuoci@polimi.it   						                       *
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

#include "basic/OpenSMOKE_Utilities.h"
#include "engine/OpenSMOKE_ReactingGas.h"
#include "basic/OpenSMOKE_Constants.h"
#include "engine/OpenSMOKE_GlobalKinetics.h"

const int OpenSMOKE_GlobalKinetics::MAXREACTIONS = 30;

void OpenSMOKE_GlobalKinetics::ErrorMessage(string message)
{
    cout << "FATAL ERROR: OpenSMOKE_GlobalKinetics" << endl;
    cout << message << endl;
    cout << "Press a key to continue... " << endl;
    getchar();
    exit(-1);
}

void OpenSMOKE_GlobalKinetics::assign_mix(OpenSMOKE_ReactingGas *mixture)
{
	mix = mixture;

	NC = mix->NumberOfSpecies();
	
	ChangeDimensions(MAXREACTIONS, NC, &lambda);
	ChangeDimensions(MAXREACTIONS, NC, &nu);

	// Default
	default_set();
}

void OpenSMOKE_GlobalKinetics::default_set()
{
	rotation_index = 0;
	Cstar = 1.e-8;
	ALFA  = 1.e-5;
	delta = 1.e9;
}

void OpenSMOKE_GlobalKinetics::change_set()
{
	rotation_index++;
	if (rotation_index == 1)
	{
		Cstar = 1.e-8;
		ALFA  = 1.e-4;
		delta = 1.e9;
	}
	else if (rotation_index == 2)
	{
		Cstar = 1.e-8;
		ALFA  = 1.e-3;
		delta = 1.e10;
	}
	else if (rotation_index == 3)
	{
		Cstar = 1.e-8;
		ALFA  = 1.e-2;
		delta = 1.e10;
	}
	else if (rotation_index == 4)
	{
		Cstar = 5.e-8;
		ALFA  = 1.e-2;
		delta = 1.e9;
	}
	else if (rotation_index == 5)
	{
		Cstar = 5.e-9;
		ALFA  = 1.e-3;
		delta = 1.e10;
	}
	else if (rotation_index == 6)
	{
		Cstar = 5.e-9;
		ALFA  = 1.e-1;
		delta = 1.e10;
	}
	else
		ErrorMessage("Rotation tolerances not available...");
}

void OpenSMOKE_GlobalKinetics::read_from_file(string fileName)
{
	ifstream fInput;
	int iReaction;
	int int_number;
	double double_number;
	char label[Constants::NAME_SIZE];
	char name_kinetic_scheme[Constants::NAME_SIZE];
	
	openInputFileAndControl(fInput, fileName.c_str());

	fInput >> label >> name_kinetic_scheme;
	checkInputFile(label, "NAME");

	nReactions = 1;
	fInput >> label;
	for(;;)
	{
		fInput >> iReaction;
		checkInputFile(label, "Reaction");
		checkInputFile(iReaction, nReactions);

		fInput >> label >> double_number;
		checkInputFile(label, "A");
		A.Append(double_number);

		fInput >> label >> double_number;
		checkInputFile(label, "Beta");
		Beta.Append(double_number);

		fInput >> label >> double_number;
		checkInputFile(label, "Tatt");
		Tatt.Append(double_number);

		string dummy;
		fInput >> label >> dummy;
		checkInputFile(label, "Equilibrium");
		if		(dummy == "Y")	iEquilibrium.Append(1);
		else if (dummy == "N")	iEquilibrium.Append(0);
		else	ErrorMessage("Equilibrium must be Y || N");

		fInput >> label;
		checkInputFile(label, "Stoichiometry");
		
		for(;;)
		{
			fInput >> label;
			if (!strcmp(label, "Kinetics"))
				break;
			
			fInput >> double_number;
			int_number = mix->recognize_species(label);
			nu[iReaction][int_number] = double_number;

		}

		for(;;)
		{
			fInput >> label;
			if (!strcmp(label, "Reaction") || !strcmp(label, "END"))
				break;
			
			fInput >> double_number;
			int_number = mix->recognize_species(label);
			lambda[iReaction][int_number] = double_number;
		}


		if (!strcmp(label, "END"))
			break;

		nReactions++;
	}
	
	int i;
	for(i=MAXREACTIONS;i>=nReactions+1;i--)
	{
		nu.DeleteRow(i);
		lambda.DeleteRow(i);
	}

	ChangeDimensions(nReactions, &kappa);
	ChangeDimensions(nReactions, &r);
	ChangeDimensions(nReactions, &rInverse);
	ChangeDimensions(nReactions, &DHReaction);
	ChangeDimensions(nReactions, &DSReaction);
	ChangeDimensions(nReactions, &sumNuij);
	ChangeDimensions(nReactions, &uKeq);
	ChangeDimensions(NC, &h);
	ChangeDimensions(nReactions, NC, &lambdaInverse);

	for(i=1;i<=nReactions;i++)
	{
		BzzVector aux = nu.GetRow(i);
		sumNuij[i] = aux.GetSumElements();
	}

	// Reaction orders: inverse reactions
	Sum(lambda, nu, &lambdaInverse);

	// Check Stoichiometry
	for(i=1;i<=nReactions;i++)
	{
		double sum = 0.;
		for(int j=1;j<=NC;j++)
			sum += nu[i][j]*mix->M[j];
		if (fabs(sum) >= 1.e-4)
		{	
			cout << "Reaction: " << i << endl;
			cout << "Sum:      " << sum << endl;
 			ErrorMessage("Wrong stoichiometry...");
		}
	}
}

void OpenSMOKE_GlobalKinetics::checkInputFile(const char *found, const char *expected)
{
	if (strcmp(found, expected))
	{
		cout << "FATAL ERROR in GlobalKinetics.inp File!" << endl;
		cout << "  Expected: " << expected << " - " << "Found: " << found << endl;
		exit(1);
	}
}

void OpenSMOKE_GlobalKinetics::checkInputFile(const int found, const int expected)
{
	if (found != expected)
	{
		cout << "FATAL ERROR in GlobalKinetics.inp File!" << endl;
		cout << "  Expected: " << expected << " - " << "Found: " << found << endl;
		exit(1);
	}
}

void OpenSMOKE_GlobalKinetics::KineticConstants(const double T)
{
	for(int i=1;i<=nReactions;i++)
		kappa[i] = A[i] * pow(T, Beta[i]) * exp(-Tatt[i]/T);
}

void OpenSMOKE_GlobalKinetics::KineticConstants(const int index, const double T)
{
	kappa = 0.;
	kappa[index] = A[index] * pow(T, Beta[index]) * exp(-Tatt[index]/T);
}

void OpenSMOKE_GlobalKinetics::ReactionRates(const double T, BzzVector &c)
{
	int i;

	// Default
	// Cstar = 1e-8
	// ALFA  = 1e-5
	// delta = 1e+9

	r			= kappa;
	for(i=1;i<=nReactions;i++)
		rInverse[i]	= kappa[i]*uKeq[i];

	double H = 1.50*log(ALFA/(1.-ALFA));
	double K = 2.00*log((1.-ALFA)/ALFA) / Cstar;

	// Direct reactions
	for(i=1;i<=nReactions;i++)
		for(int j=1;j<=NC;j++)
			if (lambda[i][j] != 0.)
			{
				if (lambda[i][j]>=1.)
				{
					r[i] *= pow(c[j], lambda[i][j]);
				}
				else
				{
					double m = (tanh(K*c[j] + H)+1.)/2.;	//con questi coefficienti inizia a decrescere per c<1.5e-5 e finisce per 9.e-6
					double gamma =	m*pow(c[j]+m/delta,lambda[i][j]) + 
									(1-m)*pow(Cstar,lambda[i][j]-1.)*c[j];
					r[i] *= gamma;
				}
			}

	// Inverse reactions
	for(i=1;i<=nReactions;i++)
	{
		if (iEquilibrium[i] == 1)
		{
			for(int j=1;j<=NC;j++)
				if (lambdaInverse[i][j] != 0.)
				{
					if (lambdaInverse[i][j]>=1.)
					{
						rInverse[i] *= pow(c[j], lambdaInverse[i][j]);
					}
					else
					{
						double delta = 1.e9;
						double m = (tanh(K*c[j] + H)+1.)/2.;	//con questi coefficienti inizia a decrescere per c<1.5e-5 e finisce per 9.e-6
						double gamma =	m*pow(c[j]+m/delta,lambdaInverse[i][j]) + 
										(1-m)*pow(Cstar,lambdaInverse[i][j]-1.)*c[j];
						rInverse[i] *= gamma;
					}
				}

			r[i] += -rInverse[i];
		}
	}
}

void OpenSMOKE_GlobalKinetics::FormationRates(BzzVector &R)
{
	R=0.;
	for(int j=1;j<=NC;j++)
		for(int i=1;i<=nReactions;i++)
				R[j] += nu[i][j]*r[i];
	
	ElementByElementProduct(R, mix->M, &R);
}

void OpenSMOKE_GlobalKinetics::GiveMeReactionHeat(const double T, BzzVector &R, double &QReaction)
{
	mix->GetMixAveragedEnthalpy_Mass(h, T);		// [J/kg]
	QReaction = -Dot(h,R);						// [J/kg] * [kg/m3.s] = [J/m3.s] 
}

void OpenSMOKE_GlobalKinetics::GiveMeSpecificEnthalpy(const double T, BzzVector &omega, double &specificEnthalpy)
{
	mix->GetMixAveragedEnthalpy_Mass(h, T);		// [J/kg]
	specificEnthalpy = Dot(omega,h);			// [J/kg]
}

void OpenSMOKE_GlobalKinetics::GiveMeKEquilibrium(const double T)
{
	mix->ComputeKEq(T, nu, sumNuij, uKeq);
}

void OpenSMOKE_GlobalKinetics::GiveMeFormationRates(const double T, BzzVector &c, BzzVector &R)
{
	KineticConstants(T);
	GiveMeKEquilibrium(T);
	ReactionRates(T, c);
	FormationRates(R);
}

void OpenSMOKE_GlobalKinetics::GiveMedRdk(char kind, const int index, const double T, BzzVector &c, BzzVector &dRdk)
{
	KineticConstants(index, T);
	ReactionRates(T, c);
	FormationRates(dRdk);

	if (kind == 'A')	dRdk *= (1./A[index]);
	if (kind == 'k')	dRdk *= (1./kappa[index]);
}


void OpenSMOKE_GlobalKinetics::SetupOptimization(BzzVectorInt &_iKind, BzzVectorInt &_iReaction, BzzVectorInt &_iSpecies)
{
	iKind		= _iKind;
	iReaction	= _iReaction;
	iSpecies	= _iSpecies;
}

void OpenSMOKE_GlobalKinetics::ChangeKineticParameters(BzzVector &newParameters)
{
	for(int i=1;i<=iKind.Size();i++)
	{	
		if (iKind[i] == 1)	// A
			A[iReaction[i]]				= newParameters[i];

		else if (iKind[i] == 2)	// Tatt
			Tatt[iReaction[i]]			= newParameters[i];

		else if (iKind[i] == 3)	// Beta
			Beta[iReaction[i]]			= newParameters[i];

		else if (iKind[i] == 4)	// lambda
			lambda[iReaction[i]][iSpecies[i]]	= newParameters[i];
	}

	// Reaction orders: inverse reactions
	Sum(lambda, nu, &lambdaInverse);
}

int OpenSMOKE_GlobalKinetics::IsTattAKineticParameter(const int index_reaction)
{
	for(int i=1;i<=iKind.Size();i++)
		if (iKind[i] == 2)
			if (iReaction[i] == index_reaction)
				return i;
	return 0;
}

void OpenSMOKE_GlobalKinetics::IsLambdaAKineticParameter(const int index_reaction, BzzVectorInt &index_parameters, BzzVectorInt &index_species)
{
	for(int i=1;i<=iKind.Size();i++)
	{	
		if (iKind[i] == 4)	// lambda
			if (iReaction[i] == index_reaction)
			{
				index_parameters.Append(i);
				index_species.Append(iSpecies[i]);
			}
	}
}

void OpenSMOKE_GlobalKinetics::SensitivityCoefficients(const string fileName, BzzVector &grid, BzzVector &T, BzzMatrix &C)
{
	int i, j, k;

	BzzVector *S_A;
	BzzVector *S_Beta;
	BzzVector *S_Tatt;
	BzzMatrix *S_lambda;
	BzzVector R(NC);

	S_A			= new BzzVector[NC+1];
	S_Beta		= new BzzVector[NC+1];
	S_Tatt		= new BzzVector[NC+1];
	S_lambda	= new BzzMatrix[NC+1];

	for(i=1;i<=NC;i++)	ChangeDimensions(nReactions, &S_A[i]);
	for(i=1;i<=NC;i++)	ChangeDimensions(nReactions, &S_Beta[i]);
	for(i=1;i<=NC;i++)	ChangeDimensions(nReactions, &S_Tatt[i]);
	for(i=1;i<=NC;i++)	ChangeDimensions(nReactions, NC, &S_lambda[i]);

	ofstream fOutput;
	openOutputFileAndControl(fOutput, fileName);
	fOutput.setf(ios::scientific);

	fOutput << "x(1)\t\t";
	int count = 2;
	for(j=1;j<=nReactions;j++)	fOutput << "r_"	<< j				<< "[kmol/m3/s](" << count++ << ")\t";
	for(i=1;i<=NC;i++)			fOutput << "R_" << mix->names[i]	<< "[kmol/m3/s](" << count++ << ")\t";

	for(i=1;i<=NC;i++)
	{
		for(j=1;j<=nReactions;j++)	fOutput << mix->names[i] << "_A"	<< j << "(" << count++ << ")\t";
		for(j=1;j<=nReactions;j++)	fOutput << mix->names[i] << "_Beta" << j << "(" << count++ << ")\t";
		for(j=1;j<=nReactions;j++)	fOutput << mix->names[i] << "_Tatt" << j << "(" << count++ << ")\t";
		for(j=1;j<=nReactions;j++)
			for(k=1;k<=NC;k++)
				if (lambda[j][k] != 0)	fOutput << mix->names[i] << "_lambda_" << mix->names[k] << "_r_" << j << "(" << count++ << ")\t";
	}
	fOutput << endl << endl;

	for(k=1;k<=grid.Size();k++)
	{
		for(i=1;i<=NC;i++)
		{
			S_A[i]		= 0.;
			S_Beta[i]	= 0.;
			S_Tatt[i]	= 0.;
			S_lambda[i] = 0.;
		}
		
		BzzVector aux = C.GetRow(k);
		GiveMeFormationRates(T[k], aux, R);
		ElementByElementProduct(R, mix->uM, &R);

		for(i=1;i<=NC;i++)
		{
			for(j=1;j<=nReactions;j++)
				if (R[i] != 0.)	
				{
				//	S_A[i][j]		=  nu[j][i]*r[j]/R[i]; 
					S_A[i][j]		=  nu[j][i]*r[j]; 
					S_Beta[i][j]	=  S_A[i][j]*Beta[j]*log(T[k]); 
					S_Tatt[i][j]	= -S_A[i][j]*Tatt[j]/T[k]; 
					for(int w=1;w<=NC;w++)
						S_lambda[i][j][w]	= S_A[i][j]*lambda[j][w]*log(C[k][w]+1.e-12); 
				}
		}

		fOutput	<< grid[k] << "\t";
		
		for(j=1;j<=nReactions;j++)		fOutput << r[j]	<< "\t";
		for(i=1;i<=NC;i++)				fOutput << R[i]	<< "\t";

		for(i=1;i<=NC;i++)
		{
			for(j=1;j<=nReactions;j++)	fOutput << S_A[i][j]	<< "\t";
			for(j=1;j<=nReactions;j++)	fOutput << S_Beta[i][j] << "\t";
			for(j=1;j<=nReactions;j++)	fOutput << S_Tatt[i][j] << "\t";
			for(j=1;j<=nReactions;j++)
				for(int w=1;w<=NC;w++)		
					if (lambda[j][w] != 0)	fOutput << S_lambda[i][j][w] << "\t";
		}
		fOutput << endl;
	}

	fOutput.close();
}
