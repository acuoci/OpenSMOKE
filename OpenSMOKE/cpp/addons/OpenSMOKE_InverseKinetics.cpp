/***************************************************************************
 *   Copyright (C) 2006-2008 by Alberto Cuoci   	                       *
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
#include "engine/OpenSMOKE_ReactingGas.h"
#include "basic/OpenSMOKE_Utilities.h"
#include "basic/OpenSMOKE_Constants.h"
#include "addons/OpenSMOKE_InverseKinetics.h"
#include "math.h"

using namespace std;

void inverse_kinetics_non_linear_regression(int model, int ex, BzzVector &b, BzzVector &x, BzzVector &y);

void OpenSMOKE_InverseKinetics::Setup(const int iReaction, double Patm, double Tmin, double Tmax, double deltaT)
{
	int i;
	double P = Patm*101325.;
	double sumLambdaInverse = 0.;

	cout << "Reaction #" << iReaction << endl;
	cout << "SumNu: "    << global->sumNuij[iReaction] << endl;


	Npoints = int((Tmax-Tmin)/deltaT) + 1;
	ChangeDimensions(Npoints, &T);
	ChangeDimensions(Npoints, &kappa);
	ChangeDimensions(Npoints, &Keq);
	ChangeDimensions(Npoints, &kappainv);
	ChangeDimensions(Npoints, &logkappainv);
	for(i=1;i<=Npoints;i++)
		T[i]= Tmin+(i-1)*deltaT;

	for(i=1;i<=Npoints;i++)
	{
		global->KineticConstants(T[i]);
		global->GiveMeKEquilibrium(T[i]);

		kappa[i]	= global->kappa[iReaction];
		Keq[i]		= 1./global->uKeq[iReaction];

		kappainv[i]	= kappa[i]/Keq[i];

		logkappainv[i] = log(kappainv[i]);
	}

	for(i=1;i<=global->NC;i++)
		if (global->lambdaInverse[iReaction][i] != 0.)
				sumLambdaInverse += global->lambdaInverse[iReaction][i];
	
	BzzMatrix X(Npoints,1);
	BzzMatrix Y(Npoints,1);
	X.SetColumn(1,T);
	Y.SetColumn(1,logkappainv);

	BzzVector b(3, log(1.e10), 0., 30000.);
	BzzNonLinearRegression o(1, X, Y, inverse_kinetics_non_linear_regression);
	o.InitializeModel(1,b);
	o.LeastSquaresAnalysis();

	BzzVector xSolution(3);
	o.GetSolution(xSolution);

	cout << endl;
	cout << "----------------------------------------------"    << endl;
	cout << "SOLUTION"											<< endl;
	cout << "----------------------------------------------"    << endl;
	cout << "Tmin   =  " << Tmin					<< " [K]"				<< endl;
	cout << "Tmax   =  " << Tmax					<< " [K]"				<< endl;
	cout << "deltaT =  " << deltaT					<< " [K]"				<< endl;
	cout << endl;
	cout << "A    =  " << exp(xSolution[1])		<< " 1/s.(m3/kmol)^"		<< sumLambdaInverse-1.	<< endl; 
	cout << "Beta =  " << xSolution[2]			<< " [-]"					<< endl; 
	cout << "Tatt =  " << xSolution[3]			<< " [K]"					<< endl; 
	cout << endl;
	cout << "A    =  " << exp(xSolution[1])*pow(1000.,sumLambdaInverse-1.)	<< " 1/s.(cm3/mol)^" << sumLambdaInverse-1.	<< endl; 
	cout << "Beta =  " << xSolution[2]			<< " [-]"					<< endl; 
	cout << "Tatt =  " << xSolution[3]			<< " [K]"					<< endl; 
	cout << endl;
	for(i=1;i<=global->NC;i++)
		if (global->lambdaInverse[iReaction][i] != 0.)
			cout << mix->names[i].c_str() << "\t" << global->lambdaInverse[iReaction][i] << endl;

	ofstream fOutput;
	openOutputFileAndControl(fOutput, "InverseKinetics.out");
	fOutput.setf(ios::scientific);
	 
	fOutput << "T[K]    "		<< "\t";
	fOutput << "kForward"		<< "\t";
	fOutput << "KEq     "		<< "\t";
	fOutput << "kInverse"		<< "\t";
	fOutput << "kInverseFitted" << "\t";
	fOutput << endl;

	for(i=1;i<=Npoints;i++)
	{
		fOutput << T[i]			<< "\t";
		fOutput << kappa[i]		<< "\t";
		fOutput << Keq[i]																<< "\t";
		fOutput << kappainv[i]															<< "\t";
		fOutput << exp(xSolution[1])*pow(T[i], xSolution[2])*exp(-xSolution[3]/T[i])	<< "\t";
		fOutput << endl;
	}

	fOutput.close();
}

void OpenSMOKE_InverseKinetics::Setup(const int iReaction, const double Patm, const double Tmin, const double Tmax, const double deltaT,
									  double &A, double &Beta, double &Tatt, BzzVector &lambdaInverse)
{
	int i;
	double P = Patm*101325.;
	double sumLambdaInverse = 0.;

	cout << "Reaction #" << iReaction << endl;
	cout << "SumNu: "    << global->sumNuij[iReaction] << endl;


	Npoints = int((Tmax-Tmin)/deltaT) + 1;
	ChangeDimensions(Npoints, &T);
	ChangeDimensions(Npoints, &kappa);
	ChangeDimensions(Npoints, &Keq);
	ChangeDimensions(Npoints, &kappainv);
	ChangeDimensions(Npoints, &logkappainv);
	for(i=1;i<=Npoints;i++)
		T[i]= Tmin+(i-1)*deltaT;

	for(i=1;i<=Npoints;i++)
	{
		global->KineticConstants(T[i]);
		global->GiveMeKEquilibrium(T[i]);

		kappa[i]	= global->kappa[iReaction];
		Keq[i]		= 1./global->uKeq[iReaction];

		kappainv[i]	= kappa[i]/Keq[i];

		logkappainv[i] = log(kappainv[i]);
	}

	for(i=1;i<=global->NC;i++)
		if (global->lambdaInverse[iReaction][i] != 0.)
				sumLambdaInverse += global->lambdaInverse[iReaction][i];
	
	BzzMatrix X(Npoints,1);
	BzzMatrix Y(Npoints,1);
	X.SetColumn(1,T);
	Y.SetColumn(1,logkappainv);

	BzzVector b(3, log(1.e10), 0., 30000.);
	BzzNonLinearRegression o(1, X, Y, inverse_kinetics_non_linear_regression);
	o.InitializeModel(1,b);
	o.LeastSquaresAnalysis();

	BzzVector xSolution(3);
	o.GetSolution(xSolution);

	A    = exp(xSolution[1]);
	Beta = xSolution[2];
	Tatt = xSolution[3];
	lambdaInverse = global->lambdaInverse.GetRow(iReaction);
}

void OpenSMOKE_InverseKinetics::AssignKineticScheme(OpenSMOKE_ReactingGas &_mix)
{
    mix     = &_mix;
}

void OpenSMOKE_InverseKinetics::AssignGlobalKineticScheme(OpenSMOKE_GlobalKinetics &_global)
{
    global = &_global;
}

void inverse_kinetics_non_linear_regression(int model, int ex, BzzVector &b, BzzVector &x, BzzVector &y)
{
	y[1] = b[1] +b[2]*log(x[1]) -b[3]/x[1];
}


