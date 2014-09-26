/***************************************************************************
 *   Copyright (C) 2009 by Alberto Cuoci								   *
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

#include "basic/OpenSMOKE_Constants.h"
#include "basic/OpenSMOKE_Utilities.h"
#include "distributions/OpenSMOKE_SinIntegralDistribution.h"

void OpenSMOKE_SinIntegralDistribution::ErrorMessage(const std::string message)
{
    cout << endl;
    cout << "Class:  OpenSMOKE_SinIntegralDistribution"			<< endl;
    cout << "Error:  " << message           << endl;
    cout << "Press a key to continue... "   << endl;
    getchar();
    exit(-1);
}

void OpenSMOKE_SinIntegralDistribution::WarningMessage(const std::string message)
{
    cout << endl;
    cout << "Class:    OpenSMOKE_SinIntegralDistribution"		<< endl;
    cout << "Warning:  "	<< message      << endl;
    cout << "Press a key to continue... "   << endl;
    getchar();
}

OpenSMOKE_SinIntegralDistribution::OpenSMOKE_SinIntegralDistribution()
{
	int j;

	Npoints = 10000;
	ChangeDimensions(Npoints, &x);
	ChangeDimensions(Npoints, &sinx);
	ChangeDimensions(Npoints, &f);
	
	dx = 2.*Constants::pi/double(Npoints-1.);

	x[1] = 0.;
	for(j=2;j<=Npoints;j++)
		x[j] = x[j-1]+dx;

	for(j=1;j<=Npoints;j++)
		sinx[j] = sin(x[j]);
}

OpenSMOKE_SinIntegralDistribution::~OpenSMOKE_SinIntegralDistribution()
{
}

void OpenSMOKE_SinIntegralDistribution::SetBoundaries(const double _Tmin, const double _Tmax)
{
	Tmin = _Tmin;
	Tmax = _Tmax;
	dT   = Tmax-Tmin;
}

void OpenSMOKE_SinIntegralDistribution::Set(const double csi, const double g)
{
	double Tmean  = Tmin + dT*csi;
	double qanm   = 2.*dT*dT/Tmean/Tmean * g;

	dTsin = sqrt(qanm)*Tmean;

	sigmaT2 = g*dT*dT;
	Tsigned = Tmean;

	/*
	BzzVector fResiduals(2);
	BzzVector xMin(2, Tmin, 0.);
	BzzVector xMax(2, Tmax, 1e16);
	
	BzzVector xFirstGuess(2, Tmean, sqrt(qanm)*Tmean);

	OpenSMOKE_SinIntegralDistributionMyNonLinearSystem nls;
	nls.AssignSinIntegralDistribution(this);

	BzzNonLinearSystemObject o(xFirstGuess, &nls);
	o.SetMinimumConstraints(xMin);
	o.SetMaximumConstraints(xMax);

	cout << "dTsin(O):     " << sqrt(qanm)*Tmean << endl;
	cout << "Tmin(O):      " << Tmean-sqrt(qanm)*Tmean << endl;
	cout << "Tmax(O):      " << Tmean+sqrt(qanm)*Tmean << endl;
	cout << "g/gmax:       " << g/(csi*(1.-csi)) << endl;
	cout << "Tmean:        " << Tmean << endl;
	
	char control = o();

	if (control<1 || control>7)
		ErrorMessage("It was impossible to obtain Tsigned and dTsin...");

	o.GetSolution(&xFirstGuess, &fResiduals);

	Tsigned = xFirstGuess[1];
	dTsin   = xFirstGuess[2];
	

	cout << "Tsigned:      " << Tsigned << endl;
	cout << "alfa:         " << alfa/2./Constants::pi*360. << " degrees" << endl;
	cout << "beta:         " << beta/2./Constants::pi*360. << " degrees" << endl;
	cout << "dTsin(N):     " << dTsin << endl;
	cout << "Tmin(N):      " << Tmean-dTsin << endl;
	cout << "Tmax(N):      " << Tmean+dTsin << endl;
	cout << "dT(N)/dT(O):  " << dTsin/(sqrt(qanm)*Tmean) << endl;
	getchar();
	*/
}

double OpenSMOKE_SinIntegralDistribution::CorrectionCoefficient(const double n, const double Tatt, const double Tmean)
{
	int j;

	if (n!=0.)
	{
		for(j=1;j<=Npoints;j++)
		{
			double T = Tsigned + dTsin*sinx[j];
			if (T<=Tmin) T=Tmin;
			
			f[j] = exp(-Tatt*(1./T-1./Tmean)) * pow(T/Tmean, n);
		}
	}
	else
	{
		for(j=1;j<=Npoints;j++)
		{
			double T = Tsigned + dTsin*sinx[j];
			if (T<=Tmin) T=Tmin;

			f[j] = exp(-Tatt*(1./T-1./Tmean));
		}
	}

	double sum = 0.;
	for(j=1;j<=Npoints-1;j++)
		sum += f[j]+f[j+1];
	sum *= (0.50*dx);

	return sum/(2.*Constants::pi);
}

void OpenSMOKE_SinIntegralDistributionMyNonLinearSystem::AssignSinIntegralDistribution(OpenSMOKE_SinIntegralDistribution *sin)
{
	ptSin	=  sin;
}

void OpenSMOKE_SinIntegralDistributionMyNonLinearSystem::GetResiduals(BzzVector &x,BzzVector &f)
{
/*	double pi = Constants::pi;
	double iLow  = false;
	double iHigh = false;

	double a = Constants::pi_over_2;
	double b = Constants::pi_over_2;
	
	double T  = x[1];
	double dT = x[2];
	double Tm = ptSin->Tmean;

	if ((ptSin->Tmax-T)/dT < 1.0)	
	{
		iHigh = true;
		a = asin((ptSin->Tmax-T)/dT);
	}
	
	if ((T-ptSin->Tmin)/dT < 1.0)	
	{
		iLow = true;
		b = asin((T-ptSin->Tmin)/dT);
	}

	double dTmax = ptSin->Tmax-Tm;
	double dTmin = Tm-ptSin->Tmin;
	double gamma = 2.*( T*T+Tm*Tm ) - 4.*T*Tm + dT*dT;
	
	ptSin->alfa = a;
	ptSin->beta = b;
	
	double f1 = 2.*pi*T;
	if (iHigh == true)
	{
		f1 += (pi-2.*a)*ptSin->Tmax;
		f1 -= 2.*cos(a)*dT-(2.*a-pi)*T;
	}
	if (iLow == true)
	{
		f1 += (pi-2.*b)*ptSin->Tmin;
		f1 -= -2.*cos(b)*dT+(pi-2.*b)*T;
	}

	double f2 = pi*gamma;
	if (iHigh == true)
	{
		f2 += (pi-2.*a)*dTmax*dTmax;
		f2 -= cos(a)*dT*(sin(a)*dT+4.*(T-Tm)) - 0.50*(2.*a-pi)*gamma;
	}
	if (iLow == true)
	{
		f2 += (pi-2.*b)*dTmin*dTmin;
		f2 -= cos(b)*dT*(sin(b)*dT+4.*(Tm-T)) - 0.50*(2.*b-pi)*gamma;
	}

	f[1] = f1 - 2.*pi*Tm;
	f[2] = f2 - 2.*pi*ptSin->sigmaT2;

	cout << "T: " << T << " " << Tm << " " << dT << endl;
	cout << "ang: " << a/2./Constants::pi*360. << " " << b/2./Constants::pi*360. << endl;
	cout << Tm << " " << f1/2./pi << endl;
*/
}