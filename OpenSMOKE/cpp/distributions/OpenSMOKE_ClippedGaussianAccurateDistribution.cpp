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
#include "distributions/OpenSMOKE_ClippedGaussianAccurateDistribution.h"

const double OpenSMOKE_ClippedGaussianAccurateDistribution::sqrt_2pi = sqrt(2.*acos(-1.));
const double OpenSMOKE_ClippedGaussianAccurateDistribution::sqrt_pi  = sqrt(acos(-1.));
const double OpenSMOKE_ClippedGaussianAccurateDistribution::sqrt_2   = sqrt(2.);

void OpenSMOKE_ClippedGaussianAccurateDistribution::ErrorMessage(const string message)
{
    cout << endl;
    cout << "Class:  OpenSMOKE_ClippedGaussianAccurateDistribution"			<< endl;
    cout << "Error:  " << message           << endl;
    cout << "Press a key to continue... "   << endl;
    getchar();
    exit(-1);
}

void OpenSMOKE_ClippedGaussianAccurateDistribution::WarningMessage(const string message)
{
    cout << endl;
    cout << "Class:    OpenSMOKE_ClippedGaussianAccurateDistribution"		<< endl;
    cout << "Warning:  "	<< message      << endl;
    cout << "Press a key to continue... "   << endl;
    getchar();
}

OpenSMOKE_ClippedGaussianAccurateDistribution::OpenSMOKE_ClippedGaussianAccurateDistribution()
{
}

OpenSMOKE_ClippedGaussianAccurateDistribution::~OpenSMOKE_ClippedGaussianAccurateDistribution()
{
}

void OpenSMOKE_ClippedGaussianAccurateDistribution::Set(const double _csi, const double _g)
{
	int j;
	iDeltaDirac = false;

	if (_csi<=0.)	ErrorMessage("The mean value must be stricly positive...");
	if (_g<=0.)		ErrorMessage("The variance value must be stricly positive...");

	if (_g >= 0.99*_csi*(1.-_csi))
		ErrorMessage("The variance is too large to be accepted...");

	csi = _csi;
	g	= _g;

	if (g <= csi*(1.-csi)/1000.)
		iDeltaDirac = true;
	else
	{
	N = 200;

	ChangeDimensions(N+1, &x);
	ChangeDimensions(N+1, &y);
	ChangeDimensions(N+1, &f1);
	ChangeDimensions(N+1, &f2);
	ChangeDimensions(N+1, &f);

	dx = 1./double(N);
	x[1] = 0.;
	for(j=2;j<=N+1;j++)		x[j]  = x[j-1]+dx;
	for(j=1;j<=N+1;j++)		f1[j] = x[j];
	for(j=1;j<=N+1;j++)		f2[j] = (x[j]-csi)*(x[j]-csi);

	BzzVector fResiduals(2);
	BzzVector xMin(2, -1e16, 1.e-16);
	BzzVector xMax(2,  1e16, 1e16);
		
	BzzVector xFirstGuess(2, csi, sqrt(g));

	OpenSMOKE_ClippedGaussianAccurateDistributionMyNonLinearSystem nls;
	nls.AssignGaussian(this);

	BzzNonLinearSystemObject o(xFirstGuess, &nls);
	o.SetMinimumConstraints(xMin);
	o.SetMaximumConstraints(xMax);

	char control = o();

	if (control<=1 || control>=7)
		ErrorMessage("It was impossible to obtain mu and sigma...");

	o.GetSolution(&xFirstGuess, &fResiduals);
	
	mu		= xFirstGuess[1];
	sigma	= xFirstGuess[2];
	Prepare_yFunction(mu,sigma);

	double constant = sqrt(2.)*sigma;
	double C0	= BzzErf((0.-mu)/constant);
	double C1	= BzzErf((1.-mu)/constant);
	double Icentral	= fabs( 0.50 * (C1-C0) );
	alfa1		= fabs( 0.50 * (C0-BzzErf((-1.e9-mu)/constant)));
	alfa2		= fabs( 0.50 * (BzzErf((1.e9-mu)/constant)-C1) );
	}
}

void OpenSMOKE_ClippedGaussianAccurateDistribution::Prepare_yFunction(const double mu, const double sigma)
{
	double coeff = 1./(2.*sigma*sigma);

	for(int j=1;j<=N+1;j++)
		y[j] = exp(-coeff*(x[j]-mu)*(x[j]-mu));
}

double OpenSMOKE_ClippedGaussianAccurateDistribution::Integral1(const double mu, const double sigma)
{
	const double C1 = sqrt_2*sigma*exp(-0.50*mu*mu/sigma/sigma);
	const double C2 = mu*BzzErf(0.50*sqrt_2*mu/sigma)*sqrt_pi;
	const double C3 = sqrt_2*sigma*exp(-0.50*(mu-1.)*(mu-1.)/sigma/sigma);
	const double C4 = mu*BzzErf(0.50*sqrt_2*(mu-1.)/sigma)*sqrt_pi;

	return 1./2./sqrt_pi*(C1+C2-C3-C4);
}

double OpenSMOKE_ClippedGaussianAccurateDistribution::Integral2(const double mu, const double sigma)
{
	const double A1 = BzzErf(0.50*sqrt_2*(mu-1.)/sigma);
	const double A2 = BzzErf(0.50*sqrt_2*mu/sigma);

	const double B1 = A1*sqrt_pi;
	const double B2 = A2*sqrt_pi;
	
	const double E1 = exp(-0.50*mu*mu/sigma/sigma);
	const double E2 = exp(-0.50*(mu-1.)*(mu-1.)/sigma/sigma);

	const double C1  = sqrt_2*sigma*mu*E1;
	const double C2  = mu*mu*B2;
	const double C3  = sigma*sigma*B2;
	const double C4  = 2*sqrt_2*csi*sigma*E1;
	const double C5  = 2.*csi*mu*B2;
	const double C6  = csi*csi*B2;
	const double C7  = csi*csi*B1;
	const double C8  = sqrt_2*sigma*E2;
	const double C9  = sqrt_2*mu*sigma*E2;
	const double C10 = mu*mu*B1;
	const double C11 = sigma*sigma*B1;
	const double C12 = 2.*sqrt_2*csi*sigma*E2;
	const double C13 = 2.*csi*mu*B1;

	return 1./2./sqrt_pi*(C1+C2+C3-C4-C5+C6-C7-C8-C9-C10-C11+C12+C13);
}

double OpenSMOKE_ClippedGaussianAccurateDistribution::CentralIntegral(const double sigma)
{
	double sum = 0.;
	for(int j=1;j<=N;j++)
		sum += (y[j]+y[j+1]);
	sum *= (0.50*dx/sigma/sqrt_2pi);

	return sum;
}

double OpenSMOKE_ClippedGaussianAccurateDistribution::IntegralNormalized()
{
	double central = 0.;
	for(int j=1;j<=N;j++)	
		central += (y[j]*f[j]+y[j+1]*f[j+1]);
	central *= (0.50*dx/sigma/sqrt_2pi);

	return alfa1*f_a + central +alfa2*f_b;
}

double OpenSMOKE_ClippedGaussianAccurateDistribution::ReactionCorrectionCoefficient(const double Tatt, const double Tmean, const double Tmin, const double Tmax)
{
	if (iDeltaDirac == true)
		return 1.;
	else
	{
		double T;
		double DeltaT = Tmax-Tmin;
		
		for(int j=1;j<=N+1;j++)	
		{
			T = Tmin + x[j]*DeltaT;
			f[j] = exp(-Tatt*(1./T-1./Tmean));
		}

		f_a = exp(-Tatt*(1./Tmin-1./Tmean));
		f_b = exp(-Tatt*(1./Tmax-1./Tmean));

		return IntegralNormalized();
	}
}

double OpenSMOKE_ClippedGaussianAccurateDistribution::ReactionCorrectionCoefficient(const double n, const double Tatt, const double Tmean, const double Tmin, const double Tmax)
{
	if (iDeltaDirac == true)
		return 1.;
	else
	{
		double T;
		double DeltaT = Tmax-Tmin;
		
		for(int j=1;j<=N+1;j++)	
		{
			T = Tmin + x[j]*DeltaT;
			f[j] = exp(-Tatt*(1./T-1./Tmean))*pow(T/Tmean, n);
		}

		f_a = exp(-Tatt*(1./Tmin-1./Tmean))*pow(Tmin/Tmean, n);
		f_b = exp(-Tatt*(1./Tmax-1./Tmean))*pow(Tmax/Tmean, n);

		return IntegralNormalized();
	}
}

void OpenSMOKE_ClippedGaussianAccurateDistributionMyNonLinearSystem::AssignGaussian(OpenSMOKE_ClippedGaussianAccurateDistribution *gaussian)
{
	ptGaussian	=  gaussian;
}

void OpenSMOKE_ClippedGaussianAccurateDistributionMyNonLinearSystem::GetResiduals(BzzVector &x,BzzVector &f)
{
	double constant = sqrt(2.)*x[2];
	double C0	= BzzErf((0.-x[1])/constant);
	double C1	= BzzErf((1.-x[1])/constant);
	double Icentral	= fabs( 0.50 * (C1-C0) );
	double alfa1	= fabs( 0.50 * (C0-BzzErf((-1000000.-x[1])/constant)));
	double alfa2	= fabs( 0.50 * (BzzErf((1000000.-x[1])/constant)-C1) );

	f[1] = ptGaussian->Integral1(x[1], x[2])+alfa1*0+alfa2*1. - ptGaussian->csi;
	f[2] = ptGaussian->Integral2(x[1],x[2])+alfa1*ptGaussian->csi*ptGaussian->csi+alfa2*BzzPow2(1.-ptGaussian->csi) - ptGaussian->g;
}

double OpenSMOKE_ClippedGaussianAccurateDistribution::GetValue(LinearInterpolation &interpolation)
{
	if (iDeltaDirac == true)
	{
		ErrorMessage("iDeltaDirac = 1");
		return interpolation(csi);
	}
	else
	{
		double csi;
		
		for(int j=1;j<=N+1;j++)	
		{
			csi  = x[j];
			f[j] = interpolation(csi);
		}

		f_a = interpolation(csi);
		f_b = interpolation(csi);

		return IntegralNormalized();
	}
}