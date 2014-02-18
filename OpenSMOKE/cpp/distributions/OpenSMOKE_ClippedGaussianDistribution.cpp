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
#include "distributions/OpenSMOKE_ClippedGaussianDistribution.h"


void OpenSMOKE_ClippedGaussianDistribution::ErrorMessage(const string message)
{
    cout << endl;
    cout << "Class:  OpenSMOKE_ClippedGaussianDistribution"			<< endl;
    cout << "Error:  " << message           << endl;
    cout << "Press a key to continue... "   << endl;
    getchar();
    exit(-1);
}

void OpenSMOKE_ClippedGaussianDistribution::WarningMessage(const string message)
{
    cout << endl;
    cout << "Class:    OpenSMOKE_ClippedGaussianDistribution"		<< endl;
    cout << "Warning:  "	<< message      << endl;
    cout << "Press a key to continue... "   << endl;
    getchar();
}

OpenSMOKE_ClippedGaussianDistribution::OpenSMOKE_ClippedGaussianDistribution()
{
	N = 200;
}

OpenSMOKE_ClippedGaussianDistribution::~OpenSMOKE_ClippedGaussianDistribution()
{
}

void OpenSMOKE_ClippedGaussianDistribution::Set(const double _csi, const double _csiV2)
{
	int j;

	if (_csi<=0.)	ErrorMessage("The mean value must be stricly positive...");
	if (_csiV2<=0.)	ErrorMessage("The variance value must be stricly positive...");

	if (_csiV2 > 0.99*_csi*(1.-_csi))
		ErrorMessage("The variance is too large to be accepted...");

	csi		= _csi;
	csiV2	= _csiV2;

	g = csiV2;
	if (g>BzzPow2(0.50*csi))		g = BzzPow2(0.50*csi);
	if (g>BzzPow2(0.50*(1.-csi)))	g = BzzPow2(0.50*(1.-csi));

	double gamma = 1.-(csi*(1.-csi)-csiV2)/(csi*(1.-csi)-g);
	
	alfa1 = gamma*(1.-csi);
	alfa2 = gamma*csi;
	x1 = csi - 2.*sqrt(g);
	x2 = csi + 2.*sqrt(g);

	ChangeDimensions(N+1, &x);
	ChangeDimensions(N+1, &y);
	ChangeDimensions(N+1, &f);

	dx = (x2-x1)/double(N);
	x[1] = x1;
	for(j=2;j<=N+1;j++)	x[j] = x[j-1]+dx;
	for(j=2;j<=N+1;j++)	y[j] = exp(-BzzPow2(x[j]-csi)/(2.*g));
}

double OpenSMOKE_ClippedGaussianDistribution::ReactionCorrectionCoefficient(const double Tatt, const double Tmean, const double Tmin, const double Tmax)
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

double OpenSMOKE_ClippedGaussianDistribution::ReactionCorrectionCoefficient(const double n, const double Tatt, const double Tmean, const double Tmin, const double Tmax)
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

double OpenSMOKE_ClippedGaussianDistribution::IntegralNormalized()
{
	double central = 0.;

	for(int j=1;j<=N;j++)	
		central += (y[j]*f[j]+y[j+1]*f[j+1]);
	central *= dx*0.50;
	central /= sqrt(2.*Constants::pi*g);

	double constant = sqrt(2.*g);
	double C0	= BzzErf((0.-csi)/constant);
	double C1	= BzzErf((1.-csi)/constant);
	double Icentral	= fabs( 0.50 * (C1-C0) );

	return alfa1*f_a + (1.-alfa1-alfa2)*central/Icentral +alfa2*f_b;
}

double OpenSMOKE_ClippedGaussianDistribution::FlatCentralIntegral()
{
	double central = 0.;
	for(int j=1;j<=N;j++)	
		central += (y[j]+y[j+1]);
	central *= dx*0.50;
	central /= sqrt(2.*Constants::pi*g);

	return central;
}

double OpenSMOKE_ClippedGaussianDistribution::GetValue(LinearInterpolation &interpolation)
{
	double csi;
	
	for(int j=1;j<=N+1;j++)	
	{
		csi = x[j];
		f[j] = interpolation(csi);
	}

	f_a = interpolation(0.);
	f_b = interpolation(1.);

	return IntegralNormalized();
}
