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

	r			= kappa;
	for(i=1;i<=nReactions;i++)
		rInverse[i]	= kappa[i]*uKeq[i];
	
	double Cstar = 1.e-8;
	double ALFA  = 1.e-5;
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
					double delta = 1.e9;
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






/*
void OpenSMOKE_GlobalKinetics::SetupKineticParametersOptimization(string fileName)
{	
	int    idummy;	
	string dummy;
	double minValue;
	double maxValue;

	ifstream fInput;
	openInputFileAndControl(fInput, fileName.c_str());

	nParameters = 0;
	for(;;)
	{
		fInput >> dummy;
		if(dummy == "END")
			break;

		if(dummy != "REACTION")
			ErrorMessage("Expected REACTION key word!");

		fInput >> idummy;
		iReaction.Append(idummy);

		fInput >> dummy;
		if(dummy == "A")
		{
			iKind.Append(1);
			iSpecies.Append(0);
		
		}
		else if(dummy == "Tatt")
		{
			iKind.Append(2);
			iSpecies.Append(0);
		}
		else if(dummy == "Beta")
		{
			iKind.Append(3);
			iSpecies.Append(0);
		}
		else if(dummy == "lambda")
		{
			iKind.Append(4);
			
			fInput >> dummy;
			iSpecies.Append(mix->recognize_species(dummy));
		}
		else
			ErrorMessage("Error code 0001");

		fInput >> minValue;
		fInput >> maxValue;
		if(minValue >= maxValue) ErrorMessage("Error code 0000");
		
		minValues.Append(minValue);
		maxValues.Append(maxValue);

		nParameters++;
	}	

	fInput.close();	

	CheckKineticParametersOptimization();

	InitializeStartingPoint();
}

void OpenSMOKE_GlobalKinetics::CheckKineticParametersOptimization()
{
	int i;
	
	for(i=1;i<=nParameters;i++)
	{
		stringstream i_string;
		stringstream iReaction_string;

		i_string 		<< i;
		iReaction_string 	<< iReaction[i];


		if (iReaction[i] > nReactions)
		{
			string message 	 = "The Reaction " + iReaction_string.str() + " is not included in the Global Kinetics!\n";
			message 	+= "Check the Optimization Table file!";
			ErrorMessage(message);
		}
		
		if (iKind[i] == 1)	// A
		{
			string message = "Please check A min and max values in Reaction " 	+ iReaction_string.str() + 
				         " in Optimization Table file at request #" 		+ i_string.str();

			if (minValues[i] < 0.)		ErrorMessage(message);
			if (maxValues[i] > 1.e32)	ErrorMessage(message);
		}
		
		if (iKind[i] == 2)	// Tatt
		{
			string message = "Please check Tatt min and max values in Reaction " 	+ iReaction_string.str() + 
				         " in Optimization Table file at request #" 		+ i_string.str();

			if (minValues[i] < -1e8)	ErrorMessage(message);
			if (maxValues[i] >  1e8)	ErrorMessage(message);
		}

		if (iKind[i] == 3)	// Beta
		{
			string message = "Please check Beta min and max values in Reaction " 	+ iReaction_string.str() + 
				         " in Optimization Table file at request #" 		+ i_string.str();

			if (minValues[i] < -1e1)	ErrorMessage(message);
			if (maxValues[i] >  1e1)	ErrorMessage(message);
		}
		
		if (iKind[i] == 4)	// lambda
		{
			string message = "Please check lambda min and max values in Reaction " 	+ iReaction_string.str() + 
				         " in Optimization Table file at request #" 		+ i_string.str();

			if (minValues[i] < -1e1)	ErrorMessage(message);
			if (maxValues[i] >  1e1)	ErrorMessage(message);
		}		
	}
	
	// Additional checks 
	for(i=1;i<=nParameters;i++)
	{
		
	}
}

void OpenSMOKE_GlobalKinetics::ChangeKineticParameters(BzzVector &newParameters)
{
	for(int i=1;i<=nParameters;i++)
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
}

void OpenSMOKE_GlobalKinetics::InitializeStartingPoint()
{
	ChangeDimensions(nParameters, &startingPoint);
	for(int i=1;i<=nParameters;i++)
	{	
		if (iKind[i] == 1)	// A
			startingPoint[i] = A[iReaction[i]];

		else if (iKind[i] == 2)	// Tatt
			startingPoint[i] = Tatt[iReaction[i]];

		else if (iKind[i] == 3)	// Beta
			startingPoint[i] = Beta[iReaction[i]];

		else if (iKind[i] == 4)	// lambda
			startingPoint[i] = lambda[iReaction[i]][iSpecies[i]];
	}
}

*/

void OpenSMOKE_GlobalKinetics::SensitivityAnalysis(const string fileName, int indexT, int dimBlock,
										 BzzFactorizedTridiagonalBlocksGauss &J,
										 BzzVector &grid, 
										 BzzVector &T, BzzMatrix &C,
										 BzzVector &rho, BzzVector &Cp
										 )
{
/*	int i, j, k;

	const double	Cmin		= 1.e-12;
	const int		Npoints		= grid.Size();
	const int		Nequations	= Npoints*dimBlock;
	const int		Nparameters	= (NC+3)*nReactions;

	BzzVector *S_A;
	BzzVector *S_Beta;
	BzzVector *S_Tatt;
	BzzMatrix *S_lambda;
	BzzVector R(NC);

	S_A			= new BzzVector[NC+2];
	S_Beta		= new BzzVector[NC+2];
	S_Tatt		= new BzzVector[NC+2];
	S_lambda	= new BzzMatrix[NC+2];

	for(i=1;i<=NC+1;i++)	ChangeDimensions(nReactions,		&S_A[i]);
	for(i=1;i<=NC+1;i++)	ChangeDimensions(nReactions,		&S_Beta[i]);
	for(i=1;i<=NC+1;i++)	ChangeDimensions(nReactions,		&S_Tatt[i]);
	for(i=1;i<=NC+1;i++)	ChangeDimensions(nReactions, NC,	&S_lambda[i]);

	BzzMatrix Jalfa(Nequations, Nparameters);

	for(k=1;k<=Npoints;k++)
	{
		for(i=1;i<=NC+1;i++)
		{
			S_A[i]		= 0.;
			S_Beta[i]	= 0.;
			S_Tatt[i]	= 0.;
			S_lambda[i] = 0.;
		}
		
		BzzVector aux = C.GetRow(k);
		GiveMeFormationRates(T[k], aux, R);			// !WARNING: r [kmol/m3/s] - R [kg/m3/s]	

		// Sensitivity coefficients for formation rates
		for(i=1;i<=NC;i++)
		{
			for(j=1;j<=nReactions;j++)
			{
				S_A[i][j]		=  (nu[j][i]*r[j]*mix->M[i])/A[j]/rho[k]; 
				S_Beta[i][j]	=  S_A[i][j]*A[j]*log(T[k]); 
				S_Tatt[i][j]	= -S_A[i][j]*A[j]/T[k]; 
				for(int w=1; w<=NC; w++)
					S_lambda[i][j][w]	= S_A[i][j]*A[j]*log(C[k][w]+Cmin); 
			}		
		}

		mix->SpeciesEnthalpy(T[k]);						// [-]
		Product((Constants::R_J_kmol*T[k]), mix->h,  &h);	// [J/kmol]
		ElementByElementProduct(h,  mix->uM, &h);		// [J/kg]

		// Sensitivity coefficients for reaction heat
		for(j=1;j<=nReactions;j++)
			for(i=1;i<=NC;i++)
			{
				S_A[NC+1][j]	-= h[i]*S_A[i][j]/Cp[k];
				S_Beta[NC+1][j] -= h[i]*S_Beta[i][j]/Cp[k];
				S_Tatt[NC+1][j] -= h[i]*S_Tatt[i][j]/Cp[k];
				for(int w=1; w<=NC; w++)
					S_lambda[NC+1][j][w]	-= h[i]*S_lambda[i][j][w]/Cp[k];
			}


		int indexRow = (k-1)*dimBlock;
		for(j=1;j<=nReactions;j++)
		{
			int indexColumn = (j-1)*(NC+3);
			Jalfa[indexRow + indexT][indexColumn+1]		= S_A[NC+1][j];
			Jalfa[indexRow + indexT][indexColumn+2]		= S_Beta[NC+1][j];
			Jalfa[indexRow + indexT][indexColumn+3]		= S_Tatt[NC+1][j];
			for(int w=1; w<=NC; w++)
				Jalfa[indexRow + indexT][indexColumn+3+w] = S_lambda[NC+1][j][w];
		}

		for(i=1;i<=NC;i++)
		{
			for(j=1;j<=nReactions;j++)
			{
				int indexColumn = (j-1)*(NC+3);
				Jalfa[indexRow + indexT + i][indexColumn+1]		= S_A[i][j];
				Jalfa[indexRow + indexT + i][indexColumn+2]		= S_Beta[i][j];
				Jalfa[indexRow + indexT + i][indexColumn+3]		= S_Tatt[i][j];
				for(int w=1; w<=NC; w++)
					Jalfa[indexRow + indexT + i][indexColumn+3+w]	= S_lambda[i][j][w];
			}
		}
	}


	cout << " - Solving linear system for sensitivity coefficients..." << endl;

	double startTime = BzzGetCpuTime();
	BzzMatrix E(Nequations, Nparameters);
	Jalfa *= -1.;
	Solve(&J, Jalfa, &E);

	E.BzzPrint();
	
	cout << "   CPU time for solving: " << BzzGetCpuTime() - startTime << " s" << endl;

	ofstream fOutput;
	openOutputFileAndControl(fOutput, fileName);
	fOutput.setf(ios::scientific);
	for(k=1;k<=Npoints;k++)
	{
		fOutput << grid[k] << "\t";
		for(j=1;j<=nReactions;j++)	// Sens di T rispetto a Tatt delle 4 reazioni
			fOutput << E[(k-1)*dimBlock+indexT][(NC+3)*(j-1)+3] << "\t";
		fOutput << endl;
	}
	fOutput.close();
*/}
