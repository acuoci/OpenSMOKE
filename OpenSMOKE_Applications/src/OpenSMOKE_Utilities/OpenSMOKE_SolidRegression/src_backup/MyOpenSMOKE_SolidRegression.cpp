/***************************************************************************
 *   Copyright (C) 2009 by Tiziano Maffei and Alberto Cuoci  	           *
 *   tiziano.maffei@chem.polimi.it   				                       *
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

#include <sstream>
#include <iomanip>
#include "basic/OpenSMOKE_Utilities.h"
#include "basic/OpenSMOKE_Constants.h"
#include "MyOpenSMOKE_SolidRegression.h"
#include "OpenSMOKE_SolidExperiment.h"

MyOpenSMOKE_SolidRegression* ptRegression;
void MyBzzModelOdeRegression(int model, int ex, BzzVector &b, BzzVector &x, BzzVector &y);

const int iEquationModel = 1;
const double MyOpenSMOKE_SolidRegression::MWchar	= 12.0;	// [kg/kmol]
const double MyOpenSMOKE_SolidRegression::psi		= 3.0;	// [-]			// TODO

string GetLabelIndex(const int count);
string GetNumber(const int count);
void ODE_Print(BzzVector &x, double t)
{
	ptRegression->MyODE_Print(x, t);
}

void MyOpenSMOKE_SolidRegression::ErrorMessage(const string message)
{
    cout << endl;
    cout << "Class:  OpenSMOKE_SolidRegression"		<< endl;
    cout << "Object: " << name_object				<< endl;
    cout << "Error:  " << message					<< endl;
    cout << "Press enter to continue... "			<< endl;
    getchar();
    exit(-1);
}

void MyOpenSMOKE_SolidRegression::WarningMessage(const string message)
{
    cout << endl;
    cout << "Class:	  OpenSMOKE_SolidRegression"	<< endl;
    cout << "Object:  " << name_object				<< endl;
    cout << "Warning: " << message					<< endl;
	cout << endl;
}

MyOpenSMOKE_SolidRegression::MyOpenSMOKE_SolidRegression()
{
	name_object = "[Name not assigned]";
	count_global = 1;
	iWriteOnFile = false;
}

void MyOpenSMOKE_SolidRegression::SetName(const string name)
{
	name_object = name;
}

void MyOpenSMOKE_SolidRegression::SetKinetics(OpenSMOKE_CharKineticScheme *_kinetics)
{
	kinetics = _kinetics;
}

void MyOpenSMOKE_SolidRegression::Setup(const string filename)
{
	PrepareKineticScheme();

	ptRegression = this;

	int i;
	string dummy;
	ifstream fInput;
	openInputFileAndControl(fInput, filename);

	while (fInput.eof()==0)
	{
		fInput >> dummy;
		list_of_names_of_files.push_back(dummy);
	}

	fInput.close();

	nCases = list_of_names_of_files.size();
	ChangeDimensions(nCases,		&indices);
	ChangeDimensions(nCases,		&initialconditions_x);
	ChangeDimensions(nCases, 2+Nadsorbed,	&initialconditions_y);
	ChangeDimensions(2+Nadsorbed,			&yMin);
	ChangeDimensions(2+Nadsorbed,			&yMax);
	

	// Set each experiment
	experiments = new OpenSMOKE_SolidExperiment[nCases+1];
	for(i=1;i<=nCases;i++)
		experiments[i].ReadFromFile(list_of_names_of_files[i-1]);

	// Assembling experiments
	numTotal = 0;
	for(i=1;i<=nCases;i++)
	{
		indices[i]				= numTotal+1;

		initialconditions_x[i]	= experiments[i].x[1]; 
		
		initialconditions_y[i][1]	= 0.; 
		initialconditions_y[i][2]	= 1.0; 
		for(int j=1;j<=Nadsorbed;j++)
			initialconditions_y[i][2+j]	= 0.0; 
		
		numTotal  += experiments[i].nPoints;
	}

	yMax[1] = 1.e0;
	yMax[2] = 1.e16;
	for(int j=1;j<=Nadsorbed;j++)
		yMax[2+j] = 1.e0;

	yMin[1] = 0.0;
	yMin[2] = 0.0;
	for(int j=1;j<=Nadsorbed;j++)
		yMin[2+j] = 0.0;
	
	ChangeDimensions(numTotal, 1, &YY);

	for(i=1;i<=nCases;i++)
	{
		cout << "-----------------------------------------"		<< endl;
		cout << "Experiment #" << i << endl;
		cout << "-----------------------------------------"		<< endl;
		cout << "  * Index:  " << indices[i] << endl;
		cout << "  * x0:     " << initialconditions_x[i] << endl;

		cout << "  * y1:     " << initialconditions_y[i][1] << endl;
		cout << "  * y2:     " << initialconditions_y[i][2] << endl;
		for(int j=1;j<=Nadsorbed;j++)
			cout << "  * y" << 2+j << " :     " << initialconditions_y[i][2+j] << endl;
		cout << endl;
	}

	Prepare();

	vector<string> list_of_gas_names;
	list_of_gas_names.push_back("list");
	list_of_gas_names.push_back("CO");
	list_of_gas_names.push_back("CO2");
	list_of_gas_names.push_back("O2");
	for(i=1;i<=nCases;i++)
		experiments[i].PrepareConcentrations(list_of_gas_names);
}

void MyOpenSMOKE_SolidRegression::Prepare()
{
	ChangeDimensions(numTotal, &temperatures_global);
	ChangeDimensions(numTotal, &utemperatures_global);
	ChangeDimensions(numTotal, &pressures_Gas_global);

	int k=1;
	for(int i=1;i<=nCases;i++)
	{
		for(int j=1;j<=experiments[i].nPoints;j++)
		{
			temperatures_global[k]   = experiments[i].temperature;
			utemperatures_global[k]  = 1./temperatures_global[k];
			pressures_Gas_global[k]  = experiments[i].pressure;
			k++;
		}
	}
	
	TMeanGlobal		= Mean(temperatures_global);
	PGasMeanGlobal   = Mean(pressures_Gas_global);

	C2 = -1./TMeanGlobal;
	C3 = pressures_Gas_global.Norm2();

	BzzVector auxiliar_vector = utemperatures_global;
	for(k=1;k<=numTotal;k++) auxiliar_vector[k] += C2;

	C1 = 1./auxiliar_vector.Norm2();

	cout << "Tmin:    " << temperatures_global.Min()	<< " K"		<< endl;
	cout << "Tmax:    " << temperatures_global.Max()	<< " K"		<< endl;
	cout << "Tmean:   " << TMeanGlobal					<< " K"		<< endl;
	cout << "PGasmin:  " << pressures_Gas_global.Min()	<< " atm"		<< endl;
	cout << "PGasmax:  " << pressures_Gas_global.Max()	<< " atm"		<< endl;
	cout << "PGasmean: " << PGasMeanGlobal				<< " atm"		<< endl;

	cout << "C1:    " << C1							<< " K"		<< endl;
	cout << "C2:    " << C2							<< " 1/K"	<< endl;
	cout << "C3:    " << C3							<< " atm"	<< endl;
}

void MyOpenSMOKE_SolidRegression::PrepareKineticScheme()
{
	// TODO
	NR			= 4;		// Number of reactions
	Nadsorbed	= 1;		// Number of surface species
	Nbulk		= 1;		// Number of bulk species
	Ngas		= 3;		// Number of bulk species
	
	ChangeDimensions(Nadsorbed, &teta);
	ChangeDimensions(Nadsorbed, &Cadsorbed);
	ChangeDimensions(Nbulk,		&Cbulk);
	ChangeDimensions(Nadsorbed, &eta);

	ChangeDimensions(NR, &kappa0);
	ChangeDimensions(NR, &Ea);
	
	ChangeDimensions(NR, &kappa);
	ChangeDimensions(NR, &rr);

	ChangeDimensions(Nadsorbed, NR, &nuAdsorbed);
	ChangeDimensions(Nadsorbed, NR, &lambdaAdsorbed);
	ChangeDimensions(Nbulk, NR, &nuBulk);
	ChangeDimensions(Nbulk, NR, &lambdaBulk);

	ChangeDimensions(NR, &nuChar);
	ChangeDimensions(Ngas, NR, &nuGas);
	ChangeDimensions(Ngas, NR, &lambdaGas);
	ChangeDimensions(NR, &lambdaCarbon);

	ChangeDimensions(Nadsorbed, &R);
	ChangeDimensions(Nadsorbed, &dteta_over_dt);

	// Number of char surface sites occupied by an adsorbed molecule
	eta[1] = 1.;

	
	// Stoichiometric coefficients: Surface
	nuAdsorbed[1][1] = 2.;
	nuAdsorbed[1][2] = 0.;
	nuAdsorbed[1][3] = 1.;
	nuAdsorbed[1][4] = -1.;
	
	// Reaction orders: Surface
	lambdaAdsorbed[1][1] = 0.;
	lambdaAdsorbed[1][2] = 1.;
	lambdaAdsorbed[1][3] = 1.;
	lambdaAdsorbed[1][4] = 1.;

	// Stoichiometric coefficients: Bulk
	nuBulk[1][1] = 0.;
	nuBulk[1][2] = -1.;
	nuBulk[1][3] = -1.;
	nuBulk[1][4] = -1.;

	// Reaction orders: Bulk
	lambdaBulk[1][1] = 0.;
	lambdaBulk[1][2] = 1.;
	lambdaBulk[1][3] = 1.;
	lambdaBulk[1][4] = 1.;


	// Stoichiometric coefficients: Char
	nuChar[1] = 0.;
	nuChar[2] = 1.;
	nuChar[3] = 1.;
	nuChar[4] = 1.;

	// Stoichiometric coefficients: Gas Species
	nuGas[1][1] = -1.;
	nuGas[1][2] = -1.;
	nuGas[1][3] = -1.;
	nuGas[2][4] =  0.;
	nuGas[2][1] = -1.;
	nuGas[2][2] = -1.;
	nuGas[2][3] = -1.;
	nuGas[3][4] =  0.;
	nuGas[3][1] = -1.;
	nuGas[3][2] = -1.;
	nuGas[3][3] = -1.;
	nuGas[3][4] =  0.;

	// Reaction orders: Gas Species
	lambdaGas[1][1] = 0.;
	lambdaGas[1][2] = 0.;
	lambdaGas[1][3] = 0.;
	lambdaGas[1][4] = 0.;
	lambdaGas[2][1] = 0.;
	lambdaGas[2][2] = 0.;
	lambdaGas[2][3] = 0.;
	lambdaGas[2][4] = 0.;
	lambdaGas[3][1] = 1.;
	lambdaGas[3][2] = 1.;
	lambdaGas[3][3] = 1.;
	lambdaGas[3][4] = 0.;

	// Reaction orders: Adsorbed Species
	lambdaCarbon[1] = 2.;
	lambdaCarbon[2] = 1.;
	lambdaCarbon[3] = 1.;
	lambdaCarbon[4] = 0.;
}


void MyOpenSMOKE_SolidRegression::Run(const int iModel, BzzVector &bFirstGuess)
{
	int numModels	= 1;
	int numX		= 1;
	int numY		= 1;

	BzzMatrix X(numTotal, numX);
	BzzMatrix Y(numTotal, numY);

	int i,j;
	int k=1;
	for(i=1;i<=nCases;i++)
	{
		for(int j=1;j<=experiments[i].nPoints;j++)
		{
			X[k][1]		= experiments[i].x[j]; 
			Y[k++][1]	= experiments[i].y[j]; 
		}
	}

	if (iEquationModel == 0)
	{

	}
	else if (iEquationModel == 1)
	{
		for(j=1;j<=NR;j++)	bFirstGuess[j]		= log(bFirstGuess[j]);
		for(j=1;j<=NR;j++)	bFirstGuess[NR+j]  /= Constants::R_J_kmol;
	}
	else if (iEquationModel == 2)
	{
		for(j=1;j<=NR;j++)	bFirstGuess[j]		= log(bFirstGuess[j]) + C2*bFirstGuess[NR+j]/Constants::R_J_kmol;
		for(j=1;j<=NR;j++)	bFirstGuess[NR+j]   = -bFirstGuess[NR+j]/Constants::R_J_kmol/C1;
	}


	BzzVector b = bFirstGuess;
	BzzVector bMin = bFirstGuess;	bMin *= 0.1;
	BzzVector bMax = bFirstGuess;	bMax *= 10.;


	openOutputFileAndControl(fLog, "Log.out");
	fLog.setf(ios::scientific);

	// Non linear regression
	BzzNonLinearRegression nonLinReg(numModels, X, Y, MyBzzModelOdeRegression);
	ode.assign(this);
	nonLinReg.InitializeModel(1,b, bMin, bMax);
//	nonLinReg.InitializeModel(1,b);
//	BzzVector v(1, 0.1);
//	nonLinReg.SetVariance(1, v);
//	nonLinReg.RobustAnalysis();

	double start_time = BzzGetCpuTime();
	nonLinReg.LeastSquaresAnalysis();
	double end_time = BzzGetCpuTime();
	fLog.close();

	nonLinReg.GetSolution(b);

	cout << "--------------------------" << endl;
	cout << "Solution" << endl;
	cout << "--------------------------" << endl;

	for(j=1;j<=NR;j++)
		cout << "k0[" << j << "] = " << b[j] << endl;

	for(j=1;j<=NR;j++)
		cout << "Ea[" << j << "] = " << b[NR+j] << endl;

	// --------------------------------------------------------------------------------------------
	// Post-Processing
	// --------------------------------------------------------------------------------------------
	for(i=1;i<=nCases;i++)
	{
		// Fitted data
		iWriteOnFile = true;
		openOutputFileAndControl(fOutput, experiments[i].name_of_file + ".fitted");
		fOutput.setf(ios::scientific);
		LabelODE_File();

		// Buzzi model 00
		if (iEquationModel == 0)
		{
			for(j=1;j<=NR;j++)	kappa0[j]	= b[j];
			for(j=1;j<=NR;j++)	Ea[j]		= b[NR+j];
		}
		else if (iEquationModel == 1)
		{
			for(j=1;j<=NR;j++)	kappa0[j]	= exp(b[j]);
			for(j=1;j<=NR;j++)	Ea[j]		= b[NR+j]*Constants::R_J_kmol;
		}
		else if (iEquationModel == 2)
		{
			for(j=1;j<=NR;j++)	kappa0[j]	= exp(b[j]+C1*C2*b[j+NR]);
			for(j=1;j<=NR;j++)	Ea[j]		= -C1*Constants::R_J_kmol*b[NR+j];
		}


		ptT			= experiments[i].temperature;		// [K]
		ptSg0		= experiments[i].Sg0;				// []
		ptSigma		= experiments[i].Sigma;				// []
		ptGas		= experiments[i].gas_c;				// [kmol/m3]		


		// Kinetic constants
		for(j=1;j<=NR;j++)
			kappa[j] = kappa0[j]*exp(-Ea[j]/Constants::R_J_kmol/ptT);		// [kmol, m, s]
//		Update(ptT);

		o.Deinitialize();
		o.SetInitialConditions(initialconditions_y.GetRow(i), initialconditions_x[i], &ode);
		o.SetMinimumConstraints(yMin);
		o.SetMaximumConstraints(yMax);
		o.StepPrint(ODE_Print);
		o(experiments[i].x[experiments[i].nPoints], experiments[i].x[experiments[i].nPoints]);
		fOutput.close();

		// Experimental data
		{
			ofstream fExp;
			openOutputFileAndControl(fExp, experiments[i].name_of_file + ".exp");
			fExp.setf(ios::scientific);
			for(int j=1;j<=experiments[i].nPoints;j++)
				fExp << experiments[i].x[j] << "\t" << experiments[i].y[j] << endl;
			fExp.close();
		}
	}	

	cout << "CPU time: " << end_time - start_time << " s" << endl;
}

void MyOpenSMOKE_SolidRegression::MyODE_Print(BzzVector &x, double t)
{
	if (iWriteOnFile == true)
	{
		int i;

		fOutput << setw(18) << left << t;
		
		fOutput << setw(18) << left << xchar;
		fOutput << setw(18) << left << dxchar_over_dt;
		fOutput << setw(18) << left << Sg;
		fOutput << setw(18) << left << dSg_over_dt;
		
		for(i=1;i<=Nadsorbed;i++)
			fOutput << setw(18) << left << teta[i];

		for(i=1;i<=Nadsorbed;i++)
			fOutput << setw(18) << left << dteta_over_dt[i];

		for(i=1;i<=NR;i++)
			fOutput << setw(18) << left << rr[i];

		fOutput << setw(18) << left << Rchar;
		for(i=1;i<=Nadsorbed;i++)
			fOutput << setw(18) << left << R[i];

		fOutput << endl;
	}
}

void MyOpenSMOKE_SolidRegression::LabelODE_File()
{
	if (iWriteOnFile == true)
	{
		int i;
		int count = 1;

		// Time
		fOutput << setw(18) << left << "t[s]" + GetLabelIndex(count++);
		
		// Char
		fOutput << setw(18) << left << "xc[-]" + GetLabelIndex(count++);
		fOutput << setw(18) << left << "dxc_dt[1/s]" + GetLabelIndex(count++);
		
		// Surface
		fOutput << setw(18) << left << "Sg[-]" + GetLabelIndex(count++);
		fOutput << setw(18) << left << "dSg_dt[1/s]" + GetLabelIndex(count++);
		
		// Adsorbed species
		for(i=1;i<=Nadsorbed;i++)
			fOutput << setw(18) << left << "Teta" + GetNumber(i) + "[-]" + GetLabelIndex(count++);
		for(i=1;i<=Nadsorbed;i++)
			fOutput << setw(18) << left << "dTetadt" + GetNumber(i) + "[-]" + GetLabelIndex(count++);

		// Reaction rates
		for(i=1;i<=NR;i++)
			fOutput << setw(18) << left << "RR" + GetNumber(i) + "[kmol/m2/s]" + GetLabelIndex(count++);

		// Formation rates
		fOutput << setw(16) << left << "FRChar[kmol/m2/s]" + GetLabelIndex(count++);
		for(i=1;i<=Nadsorbed;i++)
			fOutput << setw(18) << left << "FR" + GetNumber(i) + "[kmol/m2/s]" + GetLabelIndex(count++);

		fOutput << endl;
	}
}

void MyBzzModelOdeRegression(int model, int ex, BzzVector &b, BzzVector &x, BzzVector &y)
{
	ptRegression->ModelOdeRegression(model, ex, b, x, y);
}

void MyOpenSMOKE_SolidRegression::ModelOdeRegression(int model, int ex, BzzVector &b, BzzVector &x, BzzVector &y)
{
	int i;

	if(ex == 1)
	{
		// From Regression parameters to Kinetic parameters
		int j;
		if (iEquationModel == 0)
		{
			for(j=1;j<=NR;j++)	kappa0[j]	= b[j];
			for(j=1;j<=NR;j++)	Ea[j]		= b[NR+j];
		}
		else if (iEquationModel == 1)
		{
			for(j=1;j<=NR;j++)	kappa0[j]	= exp(b[j]);
			for(j=1;j<=NR;j++)	Ea[j]		= b[NR+j]*Constants::R_J_kmol;
		}
		else if (iEquationModel == 2)
		{
			for(j=1;j<=NR;j++)	kappa0[j]	= exp(b[j]+C1*C2*b[j+NR]);
			for(j=1;j<=NR;j++)	Ea[j]		= -C1*Constants::R_J_kmol*b[NR+j];
		}

		// Write on file
		{
			fLog << setw(8) << left << count_global	<< "\t";
			for(j=1;j<=NR;j++)	fLog <<kappa0[j]	<< "\t";
			for(j=1;j<=NR;j++)	fLog <<Ea[j]		<< "\t";
			fLog << endl;
		}

		// Write on video
		if ((count_global-1)%100 == 0)	cout << count_global << endl;

		// Update global counter
		count_global++;

		for(i=1;i<=nCases;i++)
		{
			int j;

			ptT			= experiments[i].temperature;
			ptGas		= experiments[i].gas_c;				// [kmol/m3]
			ptSg0		= experiments[i].Sg0;
			ptSigma		= experiments[i].Sigma;


			// Kinetic constants
			for(j=1;j<=NR;j++)
				kappa[j] = kappa0[j]*exp(-Ea[j]/Constants::R_J_kmol/ptT);		// [kmol, m, s]
//			Update(ptT);

			// Solve ODE System
			o.Deinitialize();
			o.SetInitialConditions(initialconditions_y.GetRow(i), initialconditions_x[i], &ode);
			o.SetMinimumConstraints(yMin);
			o.SetMaximumConstraints(yMax);

			for(j=1;j<=experiments[i].nPoints;j++)
			{
				yy = o(experiments[i].x[j]);
				YY.SetRow(indices[i]+j-1, yy[1]); 
			}
		}
	}

	// Recover regression data
	YY.GetRow(ex, &y);
}

void MyOpenSMOKE_SolidRegression::GetModel_01(BzzVector &y, double t, BzzVector &f)
{
	int i, j;

	// Recover unknowns
	xchar = y[1];									// [-]
	Sg = y[2];										// [m2/kg]
	for(i=1;i<=Nadsorbed;i++)
		teta[i] = y[2+i];							// [-]
	
	double etaf  = 1.;
	double tetaf = 1. - teta.GetSumElements();		// [-]

	// Surface concentrations
	for(i=1;i<=Nadsorbed;i++)
		Cadsorbed[i] = teta[i]*ptSigma/Constants::Nav_kmol/eta[i];				// [kmol/m2]
	double Cf = tetaf*ptSigma/Constants::Nav_kmol/etaf;					// [kmol/m2]


	// Bulk concentrations
	for(i=1;i<=Nbulk;i++)
		Cbulk[i] = 1.;														// [kmol/m2]

	// Reaction rates
	rr = kappa;							

	for(j=1;j<=NR;j++)	// S
		for(i=1;i<=Nadsorbed;i++)
			rr[j] *= pow(Cadsorbed[i], lambdaAdsorbed[i][j]);			// [kmol/m2/s]
	
	for(j=1;j<=NR;j++)	// B
		for(i=1;i<=Nbulk;i++)
			rr[j] *= pow(Cbulk[i], lambdaBulk[i][j]);			// [kmol/m2/s]

	for(j=1;j<=NR;j++)	// Gas
		for(i=1;i<=Ngas;i++)
			rr[j] *= pow(ptGas[i], lambdaGas[i][j]);			// [kmol/m2/s]

	for(j=1;j<=NR;j++)	// Cf
		rr[j] *= pow(Cf, lambdaCarbon[j]);				// [kmol/m2/s]

	// Char reaction rates
	Rchar = Dot(rr, nuChar);							// [kmol/m2/s]
	
	// Formation rates
	for(i=1;i<=Nadsorbed;i++)
		R[i] = Dot(rr, nuAdsorbed.GetRow(i));

	// Equation 1: Char conversion
	dxchar_over_dt = ptSg0*Sg * MWchar * Rchar;				// [1/s]

	// Equation 2: Surface
	dSg_over_dt = ( psi*(1.-xchar)/2./Sg - Sg/(1.-xchar+1.e-12) ) * dxchar_over_dt;		// [m2/kg/s]

	// Equation 3...Nadsorbed+2: 
	for(i=1;i<=Nadsorbed;i++)
		dteta_over_dt[i] = - teta[i]/Sg*dSg_over_dt + Constants::Nav_kmol*eta[i]/ptSigma*R[i] ;	

	// Recover residuals
	f[1] = dxchar_over_dt;
	f[2] = dSg_over_dt;

	for(i=1;i<=Nadsorbed;i++)
		f[2+i] = dteta_over_dt[i];

	if (xchar > 0.999)	f = 0.;

//	cout << t << " " << xchar << endl;
}


void MyOdeSystem::GetSystemFunctions(BzzVector &x,double t,BzzVector &f)
{
	ptMyRegression->GetModel_01(x, t, f);
}

void MyOdeSystem::assign(MyOpenSMOKE_SolidRegression *myregression)
{
	ptMyRegression = myregression;
}

string GetLabelIndex(const int count)
{
	stringstream number;
	number << count;
	string label = "(" + number.str() + ")";
	return label;
}

string GetNumber(const int count)
{
	stringstream number;
	number << count;
	return number.str();
}
