/***************************************************************************
 *   Copyright (C) 2003-2008 by                                            *
 *   Guido Buzzi-Ferraris, Alessio Frassoldati and Alberto Cuoci		   *						   *
 *                                                                         *
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
#include "kinpp/OpenSMOKE_CSTRNetwork_CorrectionCoefficients.h"

const double	OpenSMOKE_CSTRNetwork_CorrectionCoefficients::epsilon		= 1.e-6;
const int		OpenSMOKE_CSTRNetwork_CorrectionCoefficients::NSTEPSPDF	= 101;

void OpenSMOKE_CSTRNetwork_CorrectionCoefficients::Setup(const double _Tmin, const double _Tmax)
{
	Tmin = _Tmin;
	Tmax = _Tmax;
	deltaT = Tmax-Tmin;

	ChangeDimensions(NSTEPSPDF, &x);
	ChangeDimensions(NSTEPSPDF, &functionBeta);
	ChangeDimensions(NSTEPSPDF, &functionCc);
	ChangeDimensions(NSTEPSPDF, &functionI);
}


void OpenSMOKE_CSTRNetwork_CorrectionCoefficients::GiveMeNormalizedVariables(const double T, const double TVN, double &csiMean, double &csiV)
{
	csiMean	= (T-Tmin)/(Tmax-Tmin);
	csiV	= TVN*T*T/(Tmax-Tmin)/(Tmax-Tmin);

	if (csiMean >= 1. || csiMean <= 0.)
	{
		cout << "The boundary values are not appropriate for this simulation" << endl;
		cout << "Tmax:  " << Tmax	<< endl;
		cout << "Tmin:  " << Tmin	<< endl;
		cout << "Tmean: " << T		<< endl;
		cout << "Press enter to continue... " << endl;
		getchar();
		exit(-1);
	}
}

double OpenSMOKE_CSTRNetwork_CorrectionCoefficients::GiveMeCorrectionCoefficient(const double csi, const double csiMean, const double Tatt, const double n)
{
	double T		= deltaT*csi		+	Tmin;
	double TMean	= deltaT*csiMean	+	Tmin;

	if (n != 0.)	return pow(T/TMean, n)*exp(-Tatt*(1./T-1./TMean));
	else			return exp(-Tatt*(1./T-1./TMean));
}

double OpenSMOKE_CSTRNetwork_CorrectionCoefficients::GiveMeCorrectionDoubleDirac(const double TMean, const double sigmaFluent, const double Tatt, const double n)
{
	int iRegion;
	double deltaT	= sqrt(sigmaFluent/2.)*TMean;
	double fTilde	= (TMean-Tmin)/(Tmax-Tmin);
	double g		= BzzPow2((deltaT)/(Tmax-Tmin));

	if (g > (1.-fTilde)*fTilde)	g = (1.-fTilde)*fTilde;

	double sqrt_g = sqrt(g);

	if		(fTilde <=0.50 && sqrt_g<=fTilde)		iRegion = 1;
	else if (fTilde >=0.50 && sqrt_g<=1.-fTilde)	iRegion = 1;
	else if (fTilde <=0.50 && sqrt_g>fTilde)		iRegion = 2;
	else if (fTilde >=0.50 && sqrt_g>1.-fTilde)		iRegion = 3;
	else											iRegion = 4;

	double fMinus;
	double fPlus;
	double wMinus;
	double wPlus;

	if (iRegion == 1)
	{
		fMinus = fTilde - sqrt_g;
		fPlus  = fTilde + sqrt_g;
		wMinus = 0.50;
		wPlus  = 0.50;
	}
	else if (iRegion == 2)
	{
		fMinus = 0.;
		fPlus  = fTilde + sqrt_g;
		wMinus = g/(fTilde*fTilde+g);
		wPlus  = fTilde/(fTilde+g/fTilde); 
	}
	else if (iRegion == 3)
	{
		fMinus = fTilde - sqrt_g;
		fPlus  = 1.;
		wMinus = (1.-fTilde)/(1.-fTilde+g/(1.-fTilde));
		wPlus  = g/((1.-fTilde)*(1.-fTilde)+g); 
	}
	else
	{
		fMinus = 0.;
		fPlus  = 1.;
		wMinus = 1.-fTilde;
		wPlus  = fTilde;
	}

	double TMinus = Tmin + fMinus*(Tmax-Tmin);
	double TPlus  = Tmin + fPlus*(Tmax-Tmin);

				
	double CoeffCorr  = (	wMinus*pow(TMinus,n) * exp(-Tatt/TMinus) + wPlus*pow(TPlus,n) * exp(-Tatt/TPlus) )	/
										(	pow(TMean,n) * exp(-Tatt/TMean) );

	return CoeffCorr;
}

double OpenSMOKE_CSTRNetwork_CorrectionCoefficients::Beta(const double a, const double b)
{
	int j;

	// Building function
	double step = (1.-2.*epsilon)/double(NSTEPSPDF-1.);
	for(j=1;j<=NSTEPSPDF;j++)
	{
		x[j]			= epsilon + step*(j-1);
		functionBeta[j]	= pow(x[j], a-1.) * pow((1.-x[j]), b-1.);
	}

	// Building Beta Function
	double beta = 0.;
	for(j=1;j<=NSTEPSPDF-1;j++)
	{
		double f = 0.50*(functionBeta[j]+functionBeta[j+1]);
		beta += f*(x[j+1]-x[j]);
	}

	beta += pow(epsilon, a)/a + pow(epsilon, b)/b;

	return beta;
}

double OpenSMOKE_CSTRNetwork_CorrectionCoefficients::GiveMeCorrectionBetaPDF( const double csiMean, 
													    const double a, const double b, const double beta,
														const double Tatt, const double n)
{
	int j;

	// Building PDF
	double step = (1.-2.*epsilon)/double(NSTEPSPDF-1.);
	for(j=1;j<=NSTEPSPDF;j++)
	{
		x[j]			= epsilon + step*(j-1);
		functionCc[j]	= GiveMeCorrectionCoefficient(x[j], csiMean, Tatt, n) * pow(x[j], a-1.) * pow((1.-x[j]), b-1.);
	}

	// Building Correction Coefficient
	double CoeffCorr = 0.;
	for(j=1;j<=NSTEPSPDF-1;j++)
	{
		double f = 0.50*(functionCc[j]+functionCc[j+1]);
		CoeffCorr += f*(x[j+1]-x[j]);
	}

	// Adding two Dirac delta at the boundaries
	double d0 = GiveMeCorrectionCoefficient(0. , csiMean, Tatt, n);
	double d1 = GiveMeCorrectionCoefficient(1. , csiMean, Tatt, n);
	CoeffCorr += pow(epsilon, a)/a * d0 + pow(epsilon, b)/b * d1;

	// Renormalizing
	CoeffCorr /= beta;

	// Return Correction Coefficient
	return CoeffCorr;
}

double OpenSMOKE_CSTRNetwork_CorrectionCoefficients::GiveMeCorrectionClippedGaussianPDF(	
								const double csiMean, const double x1,  const double x2, 
								const double alfa1, const double alfa2, const double g, const double Tatt, const double n)
{
	int j;

	// Building PDF
	double step = (x2-x1)/double(NSTEPSPDF-1);
	for(j=1;j<=NSTEPSPDF;j++)
		functionCc[j]	= GiveMeCorrectionCoefficient(x[j], csiMean, Tatt, n) * functionI[j];

	// Building Correction Coefficient
	double CoeffCorr = 0.;
	for(j=1;j<=NSTEPSPDF-1;j++)
	{
		double f = 0.50*(functionCc[j]+functionCc[j+1]);
		CoeffCorr += f*step;
	}
	CoeffCorr /= sqrt(2.*Constants::pi*g);


	// Adding two Dirac delta at the boundaries
	double d0 = GiveMeCorrectionCoefficient(0. , csiMean, Tatt, n);
	double d1 = GiveMeCorrectionCoefficient(1. , csiMean, Tatt, n);
	CoeffCorr = alfa1*d0 + alfa2*d1 + (1.-alfa1-alfa2)*CoeffCorr/I;

	// Return Correction Coefficient
	return CoeffCorr;
}

void OpenSMOKE_CSTRNetwork_CorrectionCoefficients::SetClippedGaussianPDF(const double csiMean, const double x1,  
												   const double x2, const double g)
{
	int j;

	// Building PDF
	double step = (x2-x1)/double(NSTEPSPDF-1);
	for(j=1;j<=NSTEPSPDF;j++)
	{
		x[j]			= x1 + step*(j-1);
		functionI[j]	= exp( - BzzPow2(x[j]-csiMean)/(2.*g) );
	}

	// Building Correction Coefficient (I)
//	I = 0.;
//	for(j=1;j<=NSTEPSPDF-1;j++)
//	{
//		double f = 0.50*(functionI[j]+functionI[j+1]);
//		I += f*step;
//	}
//	I /= sqrt(2.*Constants::pi*g);

	double constant = sqrt(2.*g);
	I  = fabs( 0.50 * (BzzErf((1.0-csiMean)/constant)-BzzErf((0.-csiMean)/constant)) );
}