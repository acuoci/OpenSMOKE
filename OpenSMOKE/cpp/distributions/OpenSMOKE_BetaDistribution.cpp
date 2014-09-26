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
#include "distributions/OpenSMOKE_BetaDistribution.h"

void OpenSMOKE_BetaFunction::ErrorMessage(const std::string message)
{
    cout << endl;
    cout << "Class:  OpenSMOKE_BetaFunction"			<< endl;
    cout << "Error:  " << message           << endl;
    cout << "Press a key to continue... "   << endl;
    getchar();
    exit(-1);
}

void OpenSMOKE_BetaFunction::WarningMessage(const std::string message)
{
    cout << endl;
    cout << "Class:    OpenSMOKE_BetaFunction"		<< endl;
    cout << "Warning:  "	<< message      << endl;
    cout << "Press a key to continue... "   << endl;
    getchar();
}

OpenSMOKE_BetaFunction::OpenSMOKE_BetaFunction()
{
	iIntegralFormulation = false;
	PrepareBrizuelaIntegral(1000);
}

OpenSMOKE_BetaFunction::~OpenSMOKE_BetaFunction()
{
}

void OpenSMOKE_BetaFunction::SetIntegralFormulation(const double N)
{
	iIntegralFormulation = true;
	PrepareBrizuelaIntegral(N);
}

void OpenSMOKE_BetaFunction::PrepareBrizuelaIntegral(const double N)
{
	int j;
	
	Nb = N;
	ChangeDimensions(Nb+1, &xb);
	ChangeDimensions(Nb+1, &uxb);
	ChangeDimensions(Nb+1, &fb);
	dh = 1./double(Nb);
	
	xb[1] = 0.;
	for(j=2;j<=Nb+1;j++)
		xb[j] = xb[j-1] + dh;

	for(j=1;j<=Nb+1;j++)
		uxb[j] = 1.-xb[j];

	fb[1]  = 0.;
	fb[Nb+1] = 0.;
}

double OpenSMOKE_BetaFunction::BrizuelaIntegral(const double a, const double b)
{
	int j;

	double a_plus_one = a+1.;
	double b_plus_one = b+1.;
	
	for(j=2;j<=Nb;j++)
		fb[j] = pow(xb[j],a_plus_one)*pow(uxb[j], b_plus_one);

	double sum = 0.;
	for(j=1;j<=Nb;j++)
		sum += 0.50*dh*(fb[j]+fb[j+1]);

	return sum;
}

double OpenSMOKE_BetaFunction::at(const double a, const double b)
{
	if (a<=0.)		ErrorMessage("The a coefficient must be positive...");
	if (b<=0.)		ErrorMessage("The b coefficient must be positive...");

	if (iIntegralFormulation == true || (a+b>170.))
	{
		double coeff = (a+b)*(a+b+1.)*(a+b+2.)*(a+b+3.)/(a*b)/(a+1.)/(b+1.);
		return coeff*BrizuelaIntegral(a, b);
	}
	else
	{
		return gamma.at(a)*gamma.at(b)/gamma.at(a+b);
	}
}

void OpenSMOKE_BetaDistribution::ErrorMessage(const std::string message)
{
    cout << endl;
    cout << "Class:  OpenSMOKE_BetaDistribution"			<< endl;
    cout << "Error:  " << message           << endl;
    cout << "Press a key to continue... "   << endl;
    getchar();
    exit(-1);
}

void OpenSMOKE_BetaDistribution::WarningMessage(const std::string message)
{
    cout << endl;
    cout << "Class:    OpenSMOKE_BetaDistribution"		<< endl;
    cout << "Warning:  "	<< message      << endl;
    cout << "Press a key to continue... "   << endl;
    getchar();
}

OpenSMOKE_BetaDistribution::OpenSMOKE_BetaDistribution()
{
	int j, k;

	epsilon = 1.e-6;
	Nsub = log10(0.1/epsilon);
	N	 = 20;
	M	 = 50;
	
	ChangeDimensions(Nsub, &dhb);
	ChangeDimensions(Nsub, &dhf);

	ChangeDimensions(Nsub*N, &xb);
	ChangeDimensions(Nsub*N, &yb);
	ChangeDimensions(Nsub*N, &xf);
	ChangeDimensions(Nsub*N, &yf);

	ChangeDimensions(M, &xcenter);
	ChangeDimensions(M, &ycenter);

	ChangeDimensions(Nsub*N, &ff);
	ChangeDimensions(Nsub*N, &fb);
	ChangeDimensions(M, &fcenter);

	// 1. Center interval
	dhcenter = (0.9-0.1)/double(M);
	xcenter[1] = 0.1 + dhcenter*0.50;
	for(j=2;j<=M;j++)
		xcenter[j] = xcenter[j-1] + dhcenter;

	// 2. Backward Interval
	for(k=1;k<=Nsub;k++)
		dhb[k] = epsilon*(pow(10.,k)-pow(10.,k-1))/double(N);

	for(k=1;k<=Nsub;k++)
		xb[(k-1)*N+1]	= epsilon*pow(10.,k-1.)	+ dhb[k]*0.50;
	
	for(k=1;k<=Nsub;k++)
		for(j=(k-1)*N+2;j<=k*N;j++)	xb[j] = xb[j-1] + dhb[k];

	// 3. Forward Interval 
	dhf=dhb;
	for(j=1;j<=Nsub*N;j++)
		xf[j] = 1.-xb[j];
}

OpenSMOKE_BetaDistribution::~OpenSMOKE_BetaDistribution()
{
}

void OpenSMOKE_BetaDistribution::Set(const double _a, const double _b)
{
	a = _a;
	b = _b;

	if (a<=0.)		ErrorMessage("The a coefficient must be positive...");
	if (b<=0.)		ErrorMessage("The b coefficient must be positive...");

	if (a>1. && b>1.)
	{
		double fmax = 1./(1.+(b-1.)/(a-1.));
		
		if (a>500. && b>500.)	ErrorMessage("The a and b coefficient are too large...");

		if (a>500.)
		{
			a = 500.;
			b = (a-1.-fmax*(a-2.))/fmax;
		}
		else if(b>500.)
		{
			b = 500.;
			a = (1.+fmax*(b-2.))/(1.-fmax);
		}
	}

	double a_minus_one = a-1.;
	double b_minus_one = b-1.;
	
	int j;
	for(j=1;j<=Nsub*N;j++)	yb[j]		= pow(xb[j],a_minus_one)*pow(xf[j],b_minus_one);
	for(j=1;j<=Nsub*N;j++)	yf[j]		= pow(xf[j],a_minus_one)*pow(xb[j],b_minus_one);

	for(j=1;j<=M;j++)	ycenter[j]	= pow(xcenter[j],a_minus_one)*pow(1.-xcenter[j],b_minus_one);

	extreme_a = pow(epsilon,a)/a;
	extreme_b = pow(epsilon,b)/b;

	BetaInf = FlatIntegral();
}

double OpenSMOKE_BetaDistribution::FlatIntegral()
{
	int j;
	double sum_partial;
	double sum = 0.;

	int k;
	for(k=1;k<=Nsub;k++)
	{
		sum_partial = 0.;	
		for(j=(k-1)*N+1;j<=k*N;j++)	
			sum_partial += yb[j];	
		sum_partial*=dhb[k];	
		sum+=sum_partial;
	}

	for(k=1;k<=Nsub;k++)
	{
		sum_partial = 0.;	
		for(j=(k-1)*N+1;j<=k*N;j++)	
			sum_partial += yf[j];	
		sum_partial*=dhf[k];	
		sum+=sum_partial;
	}

	sum_partial = 0.;	
	for(j=1;j<=M;j++)		
		sum_partial += ycenter[j];	
	sum_partial*=dhcenter;	
	sum+=sum_partial;

	sum += extreme_a + extreme_b;

	return sum;
}

double OpenSMOKE_BetaDistribution::IntegralNormalized()
{
	int j;
	double sum_partial;
	double sum = 0.;

	int k;
	for(k=1;k<=Nsub;k++)
	{
		sum_partial = 0.;	
		for(j=(k-1)*N+1;j<=k*N;j++)	
			sum_partial += yb[j]*fb[j];	
		sum_partial*=dhb[k];	
		sum+=sum_partial;
	}

	for(k=1;k<=Nsub;k++)
	{
		sum_partial = 0.;	
		for(j=(k-1)*N+1;j<=k*N;j++)	
			sum_partial += yf[j]*ff[j];	
		sum_partial*=dhf[k];	
		sum+=sum_partial;
	}

	sum_partial = 0.;	
	for(j=1;j<=M;j++)		
		sum_partial += ycenter[j]*fcenter[j];	
	sum_partial*=dhcenter;	
	sum+=sum_partial;

	sum += extreme_a*f_a + extreme_b*f_b;

	return sum/BetaInf;
}

double OpenSMOKE_BetaDistribution::ReactionCorrectionCoefficient(const double Tatt, const double Tmean, const double Tmin, const double Tmax)
{
	int j;

	double T;
	double DeltaT = Tmax-Tmin;
	
	for(j=1;j<=Nsub*N;j++)	
	{
		T = Tmin + xb[j]*DeltaT;
		fb[j] = exp(-Tatt*(1./T-1./Tmean));
	}

	for(j=1;j<=Nsub*N;j++)	
	{
		T = Tmin + xf[j]*DeltaT;
		ff[j] = exp(-Tatt*(1./T-1./Tmean));
	}

	for(j=1;j<=M;j++)	
	{
		T = Tmin + xcenter[j]*DeltaT;
		fcenter[j]	= exp(-Tatt*(1./T-1./Tmean));
	}

	f_a = exp(-Tatt*(1./Tmin-1./Tmean));
	f_b = exp(-Tatt*(1./Tmax-1./Tmean));

	return IntegralNormalized();
}

double OpenSMOKE_BetaDistribution::ReactionCorrectionCoefficient(const double n, const double Tatt, const double Tmean, const double Tmin, const double Tmax)
{
	int j;

	double T;
	double DeltaT = Tmax-Tmin;
	
	for(j=1;j<=Nsub*N;j++)	
	{
		T = Tmin + xb[j]*DeltaT;
		fb[j] = exp(-Tatt*(1./T-1./Tmean)) * pow(T/Tmean, n);
	}

	for(j=1;j<=Nsub*N;j++)	
	{
		T = Tmin + xf[j]*DeltaT;
		ff[j] = exp(-Tatt*(1./T-1./Tmean)) * pow(T/Tmean, n);
	}

	for(j=1;j<=M;j++)	
	{
		T = Tmin + xcenter[j]*DeltaT;
		fcenter[j]	= exp(-Tatt*(1./T-1./Tmean)) * pow(T/Tmean, n);
	}

	f_a = exp(-Tatt*(1./Tmin-1./Tmean)) * pow(Tmin/Tmean, n);
	f_b = exp(-Tatt*(1./Tmax-1./Tmean)) * pow(Tmax/Tmean, n);

	return IntegralNormalized();
}

double OpenSMOKE_BetaDistribution::GetValue(LinearInterpolation &interpolation)
{
	int j;
	double csi;

	for(j=1;j<=Nsub*N;j++)	
	{
		csi = xb[j];
		fb[j] = interpolation(csi);
	}

	for(j=1;j<=Nsub*N;j++)	
	{
		csi = xf[j];
		ff[j] = interpolation(csi);
	}

	for(j=1;j<=M;j++)	
	{
		csi = xcenter[j];
		fcenter[j]	= interpolation(csi);
	}

	f_a = interpolation(0.);
	f_b = interpolation(1.);

	return IntegralNormalized();
}

