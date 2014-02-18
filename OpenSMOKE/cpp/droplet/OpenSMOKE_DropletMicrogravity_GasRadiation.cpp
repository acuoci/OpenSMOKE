/***************************************************************************
 *   Copyright (C) 2010 by Alberto Cuoci								   *
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

#include <iomanip>
#include "droplet/OpenSMOKE_DropletMicrogravity_GasRadiation.h"


void OpenSMOKE_DropletMicrogravity_GasRadiation::ErrorMessage(const string message)
{
    cout << endl;
    cout << "Class:  OpenSMOKE_Droplet"	<< endl;
    cout << "Object: " << name_object			<< endl;
    cout << "Error:  " << message				<< endl;
    cout << "Press a key to continue... "		<< endl;
    getchar();
    exit(-1);
}

void OpenSMOKE_DropletMicrogravity_GasRadiation::WarningMessage(const string message)
{
    cout << endl;
    cout << "Class:  OpenSMOKE_DropletMicrogravity"	<< endl;
    cout << "Object: "		<< name_object		<< endl;
    cout << "Warning:  "	<< message			<< endl;
    cout << "Press a key to continue... "		<< endl;
    getchar();
}

OpenSMOKE_DropletMicrogravity_GasRadiation::OpenSMOKE_DropletMicrogravity_GasRadiation()
{
	mode		= VISKANTA_MODE_ANALYTICAL;

	N_K			= 20;
	N_subK		= 20;
	N_Psi		= N_K*3;
	N_subPsi	= N_subK*3;

	G3 = 0.;
	qSurface = 0.;
}

void OpenSMOKE_DropletMicrogravity_GasRadiation::Initialize(const int N_)
{
	N = N_;
	ChangeDimensions(N, &divq);
	ChangeDimensions(N, &G0);
	ChangeDimensions(N, &G1);
	ChangeDimensions(N, &G2);
	ChangeDimensions(N, &rTilde);
	ChangeDimensions(N, &TTilde);
	ChangeDimensions(N, &KpTilde);
	ChangeDimensions(N, &y_aux);
	ChangeDimensions(N, N, &Kernel);
}

void OpenSMOKE_DropletMicrogravity_GasRadiation::SetOptions(const vector<string> options)
{
	if (options[0] == "Analytical")			mode = VISKANTA_MODE_ANALYTICAL;
	else if (options[0] == "Numerical")		mode = VISKANTA_MODE_NUMERICAL;
	else ErrorMessage("Wrong options!");
	 
	if ( atoi(options[1].c_str()) >= 3 && atoi(options[1].c_str()) <= 60)
	{
		N_K			= atoi(options[1].c_str());
		N_Psi		= N_K*3;
	}
	else ErrorMessage("Wrong options!");

	if (atoi(options[2].c_str()) >= 3 && atoi(options[2].c_str()) <= 60)
	{
		N_subK			= atoi(options[2].c_str());
		N_subPsi		= N_subK*3;
	}
}

void OpenSMOKE_DropletMicrogravity_GasRadiation::Calculate(const BzzVector &x_, const BzzVector &T_, const BzzVector &Kp_)
{
	if (N != x_.Size())
		ErrorMessage("Vector sizes do not match during the calculation!");

	T  = T_;
	Kp = Kp_;

	for (int j=1;j<=N;j++)
	{
		rTilde[j]  = x_[j]/x_[1];
		TTilde[j]  = T[j]/T[1];
		KpTilde[j] = Kp[j]*x_[1];
	}

	KpTilde_interpolation(rTilde, KpTilde);
	rhoTilde = rTilde;

	// Fill kernel
	FillKernel();
	
	// Calculate divq	[W/m3]
	divq[1] = 0.;
	for (int j=2;j<=N;j++)
		divq[j] = div_q(j);

	// Calculate surface radiation	[W/m2]
	G3 = g3();
	qSurface = Constants::sigma*BzzPow4(T[1])*(1.-G3);
}

double OpenSMOKE_DropletMicrogravity_GasRadiation::Phi(const double r, const double rho, const double z)
{
	return (r+rho-z)*(r+rho+z)*(r-rho+z)*(r-rho-z)/(4.*z*z);
}

void OpenSMOKE_DropletMicrogravity_GasRadiation::FillKernel()
{
	for(int k=1;k<=N;k++)
		for(int j=k;j<=N;j++)
			Kernel[k][j] = K(rTilde[k], rhoTilde[j]);
}

double OpenSMOKE_DropletMicrogravity_GasRadiation::div_q( const int j )
{
	G0[j] = 2.*BzzPow4(TTilde[j]);
	G1[j] = 0.; 
	G2[j] = g2(j);
	return 2.*Constants::sigma*BzzPow4(T[1]) * Kp[j] * (G0[j]-G1[j]-G2[j]);
}

double OpenSMOKE_DropletMicrogravity_GasRadiation::g2( const int k )
{
	// Symmetric lower part (indices are switched)
	for(int j=1;j<k;j++)
		y_aux[j] = K(rTilde[j], rhoTilde[k]) * rhoTilde[j] * BzzPow4(TTilde[j]) * KpTilde[j] ;
	
	// Symmetric upper part (indecies are taken in the original order)
	for(int j=k;j<=N;j++)
		y_aux[j] = K(rTilde[k], rhoTilde[j]) * rhoTilde[j] * BzzPow4(TTilde[j]) * KpTilde[j] ;

	double sum = 0.;
	for(int j=1;j<=N-1;j++)
		sum += 0.50*(y_aux[j]+y_aux[j+1]) * (rhoTilde[j+1]-rhoTilde[j]);

	return sum/rTilde[k];
}

double OpenSMOKE_DropletMicrogravity_GasRadiation::g3()
{
	for(int j=1;j<=N;j++)
		y_aux[j] = Psi(rhoTilde[j]) * rhoTilde[j] * BzzPow4(TTilde[j]) * KpTilde[j] ;

	double sum = 0.;
	for(int j=1;j<=N-1;j++)
		sum += 0.50*(y_aux[j]+y_aux[j+1]) * (rhoTilde[j+1]-rhoTilde[j]);

	return 2.*sum;
}

double OpenSMOKE_DropletMicrogravity_GasRadiation::K(const double r, const double rho)
{
	const double z1 = fabs(r-rho);
	const double z2 = sqrt(r*r-1.) + sqrt(rho*rho-1.);

	const double deltaz = (z2-z1)/double(N_K);

	double sum = 0.;
	
	if (mode == VISKANTA_MODE_NUMERICAL)
	{
		for(int j=1;j<=N_K;j++)
		{
			double z = (z1+deltaz/2.) + deltaz*(j-1);
			sum += exp ( -SubIntegral_K(r, rho, z) ) / z  * deltaz;
		}
	}
	else if (mode == VISKANTA_MODE_ANALYTICAL)
	{
		for(int j=1;j<=N_K;j++)
		{
			double z = (z1+deltaz/2.) + deltaz*(j-1);
			sum += exp ( -SubIntegral_K(r, rho, z, KpTilde_interpolation(r)) ) / z  * deltaz;
		}
	}

	return sum;
}

double OpenSMOKE_DropletMicrogravity_GasRadiation::SubIntegral_K(const double r, const double rho, const double z)
{
	double phi = Phi(r, rho, z);

	double a = fabs(r*r-rho*rho)/(2.*z);
	double b = z/2.;
	double u1 = a-b;
	double u2 = a+b;
	
	const double deltau = (u2-u1)/double(N_subK);

	double sum = 0.;
	for(int j=1;j<=N_subK;j++)
	{
		double u = (u1+deltau/2.) + deltau*(j-1);
		sum += KpTilde_interpolation( sqrt(u*u-phi) ) * deltau;
	}

	return sum;
}

double OpenSMOKE_DropletMicrogravity_GasRadiation::SubIntegral_K(const double r, const double rho, const double z, const double kptilde)
{
	double phi = Phi(r, rho, z);

	double a = fabs(r*r-rho*rho)/(2.*z);
	double b = z/2.;
	double c = a-b;
	double d = a+b;
	double e = sqrt(c*c-phi); 
	double f = sqrt(d*d-phi); 

	return kptilde/2. * (log(e+c)*phi-c*e - log(f+d)*phi+d*f);
}

double OpenSMOKE_DropletMicrogravity_GasRadiation::Psi(const double rho)
{
	const double z1 = sqrt(rho*rho-1.);
	const double z2 = rho;

	const double deltaz = (z2-z1)/double(N_Psi);

	double sum = 0.;

	if (mode == VISKANTA_MODE_NUMERICAL)
	{
		for(int j=1;j<=N_K;j++)
		{
			double z = (z1+deltaz/2.) + deltaz*(j-1);
			sum += exp ( -SubIntegral_Psi(rho, z) ) * deltaz;
		}
	}
	else if (mode == VISKANTA_MODE_ANALYTICAL)
	{
		for(int j=1;j<=N_K;j++)
		{
			double z = (z1+deltaz/2.) + deltaz*(j-1);
			sum += exp ( -SubIntegral_Psi(rho, z, KpTilde_interpolation(rho)) ) * deltaz;
		}
	}

	return sum;
}

double OpenSMOKE_DropletMicrogravity_GasRadiation::SubIntegral_Psi(const double rho, const double z)
{
	double u1 = sqrt(z*z-rho*rho+1.);
	double u2 = z;
	
	const double deltau = (u2-u1)/double(N_subPsi);

	double sum = 0.;
	for(int j=1;j<=N_subPsi;j++)
	{
		double u = (u1+deltau/2.) + deltau*(j-1);
		sum += KpTilde_interpolation( sqrt(u*u+rho*rho-z*z) ) * deltau;
	}

	return sum;
}

double OpenSMOKE_DropletMicrogravity_GasRadiation::SubIntegral_Psi(const double rho, const double z, const double kptilde)
{
	double a = z*z-rho*rho;
	double b = sqrt(z*z-a);
	double c = sqrt(a+1.);
	
	return kptilde/2. * ( a*log(b+z)-z*b -a*log(c+1)+c);
}