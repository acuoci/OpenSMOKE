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

#include "qmom/OpenSMOKE_PhysicalModels.h"

void OpenSMOKE_PhysicalModels::setup(int N)
{
	_N = N;
	_2N = 2*_N;
	_2N_minus_1 = _2N - 1;

	ChangeDimensions(_N, &w);
	ChangeDimensions(_N, &csi);

	ChangeDimensions(_2N, &localSource);
	ChangeDimensions(_2N, &Source);
	ChangeDimensions(_2N, &coeffA1);
	ChangeDimensions(_2N, &coeffA2);
	ChangeDimensions(_2N, &coeffA3);

	ChangeDimensions(_2N, _N,  &A3);
	ChangeDimensions(_2N, _N, &powerOfCsi);

	ChangeDimensions(_2N, N, &auxiliary_2NxN);

	ChangeDimensions(_N, _N,	&betaKernel);
	ChangeDimensions(_N,		&aKernel);
	ChangeDimensions(_2N, _N,	&bKernel);

	ChangeDimensions(_N, &Correction);

	prepareCoefficients();
	prepareA3();
}

void OpenSMOKE_PhysicalModels::update(BzzVector &_w, BzzVector &_csi)
{
	w	= _w;
	csi = _csi;

	computePowerOfCsi();
	prepareA3();
	assemblingA3();

	Source = 0.;
}

void OpenSMOKE_PhysicalModels::prepareCoefficients()
{
	int i;

	// -----------------------------------------------------
	// Coefficients
	// -----------------------------------------------------
	for(i=2;i<=_2N_minus_1;i++)
	{
		coeffA1[i+1] = (1-i);
		coeffA2[i+1] = i;
	}
	for(i=3;i<=_2N_minus_1;i++)
		coeffA3[i+1] = i*(i-1);
}

void OpenSMOKE_PhysicalModels::prepareA3()
{
	int i;
	// -----------------------------------------------------
	// First three rows of matrix A3
	// -----------------------------------------------------
	for(i=1;i<=_N;i++)
	{
		A3[1][i] = 0.;
		A3[2][i] = 0.;
		A3[3][i] = 2.;
	}
}

void OpenSMOKE_PhysicalModels::computePowerOfCsi()
{
	int i, j;

	// -----------------------------------------------------
	// Matrix powerOfCsi
	// -----------------------------------------------------
	powerOfCsi.SetRow(1, 1.);
	for(i=2;i<=_2N;i++)
		for(j=1;j<=_N;j++)
			powerOfCsi[i][j] = powerOfCsi[i-1][j]*csi[j];

}

void OpenSMOKE_PhysicalModels::assemblingA3()
{
	int i, j;
	int row;

	// The first three rows of A3 are constant and therefore it is not
	// necessary to calculate them
	for(i=4;i<=_2N;i++)
	{
		row = i-3;		
		for(j=1;j<=_N;j++)
			A3[i][j] = coeffA3[i]*powerOfCsi[row+1][j];
	}
}

void OpenSMOKE_PhysicalModels::powerLawGrowth(double alfa, double rGrowth)
{
	int k, i;

	localSource = 0.;
	for(k=1;k<=_2N_minus_1;k++)
	{
		for(i=1;i<=_N;i++)
			localSource[k+1]+=w[i]*pow(csi[i], k+rGrowth-1);
		localSource[k+1] *= k;
	}
	
	localSource *= alfa;

	Source += localSource;
}

void OpenSMOKE_PhysicalModels::homogeneousDispersion(BzzVector &D)
{
	// Diagonal Matrix
	BzzMatrixDiagonal P;
	P = w;

	// First Step: a = A3 x P   :: we obtain a 2N x N matrix
	Product(A3, P, &auxiliary_2NxN);

	// Second Step: localSource = A x D  
	Product(auxiliary_2NxN, D, &localSource);

	// Third Step: Assembling source  
	Source += localSource;
}

void OpenSMOKE_PhysicalModels::nucleation(double epsilon, double J, BzzVector &L)
{
	int j,k;
	double sum;

//  For uniform distributions
//	localSource = J;
//	for(k=1;k<=_2N_minus_1;k++)
//		localSource[k+1] *= pow(epsilon, k) / double(k+1.);

//	General Case
	localSource = J;
	for(k=1;k<=_2N_minus_1;k++)
	{
		sum = 0.;
		for(j=1;j<=_N;j++)
			sum += pow(L[j], k);
		localSource[k+1] *= pow(epsilon, k)*sum / double(_N);
	}

	Source += localSource;	
}

void OpenSMOKE_PhysicalModels::homogeneousAggregation()
{
	int i, j, k;
	double sum, sumB, sumD;

	for(k=0;k<=_2N_minus_1;k++)		// k=0..2N-1
	{
		// Birth rate for aggregation
		sumB = 0.;
		for(i=1;i<=_N;i++)
		{
			sum = 0.;
			for(j=1;j<=_N;j++)
				sum += ( pow((powerOfCsi[4][i]+powerOfCsi[4][j]),k/3.) *
				         betaKernel[i][j] * w[j] );
			sumB += (sum*w[i]);
		}
		
		// Death rate for aggregation
		sumD = 0.;
		for(i=1;i<=_N;i++)
		{
			sum = 0.;
			for(j=1;j<=_N;j++)
				sum += ( betaKernel[i][j] * w[j] );
			sumD += (sum*w[i]*powerOfCsi[k+1][i]);
		}

		localSource[k+1] = (0.50*sumB-sumD);
	}

	Source += localSource;
}

void OpenSMOKE_PhysicalModels::homogeneousBreakage()
{
	int i, k;

	double sumB, sumD;

	for(k=0;k<=_2N_minus_1;k++)		// k=0..2N-1
	{
		// Birth rate for aggregation
		sumB = 0.;
		for(i=1;i<=_N;i++)
			sumB += (bKernel[k+1][i]*w[i]*aKernel[i]);

		// Death rate for aggregation
		sumD = 0.;
		for(i=1;i<=_N;i++)
			sumD += (powerOfCsi[k+1][i]*w[i]*aKernel[i]);

		localSource[k+1] = sumB-sumD;
	}

	Source += localSource;
}

void OpenSMOKE_PhysicalModels::diffusiveTermCorrectionForDQMOM(BzzVector &dcsi_over_dx, 
													 double Dmix)
{	
	int i;

	// First Step: C calculation
	for(i=1;i<=_N;i++)
		Correction[i] = w[i]*Dmix*dcsi_over_dx[i]*dcsi_over_dx[i];

	// Second Step: localSource = A x D  
	Product(A3, Correction, &localSource);

	// Third Step: Assembling source  
	Source += localSource;

//	localSource.BzzPrint("localSource");
//	dcsi_over_dx.BzzPrint("dcsi");	
}

void OpenSMOKE_PhysicalModels::diffusiveTermCorrectionForDQMOM(BzzVector &Cvector)
{
	// First Step: localSource = A x C  
	Product(A3, Cvector, &localSource);

	// Second Step: Assembling source  
	Source += localSource;
}

void OpenSMOKE_PhysicalModels::aggregationKernel(int kind, double parameter)
{
	int i, j;

	switch (kind)
	{
		case 0:
			// Case 0: Constant (Marchisio 2005, p. 56)
			betaKernel = parameter;
		break;

		case 1:
			// Case 1: Hydrodynamic aggregation (Marchisio 2005, p. 56)
			for(i=1;i<=_N;i++)
				for(j=1;j<=_N;j++)
					betaKernel[i][j] = BzzPow3(csi[i]+csi[j]);
		break;

		case 2:
			// Case 2: Sum aggregation (Marchisio 2003, p. 325)
			for(i=1;i<=_N;i++)
				for(j=1;j<=_N;j++)
					betaKernel[i][j] = powerOfCsi[4][i]+powerOfCsi[4][j];
		break;

		case 3:
			// Case 3: Brownian aggregation (Marchisio 2003, p. 325)
			for(i=1;i<=_N;i++)
				for(j=1;j<=_N;j++)
					betaKernel[i][j] = BzzPow2(csi[i]+csi[j])/(csi[i]*csi[j]);
		break;

		case 4:
			// Case 4: Differential force aggregation (Marchisio 2003, p. 325)
			for(i=1;i<=_N;i++)
				for(j=1;j<=_N;j++)
					betaKernel[i][j] = BzzPow2(csi[i]+csi[j]) * 
					             fabs(powerOfCsi[3][i]-powerOfCsi[3][j]);
		break;
		
		default:
			cout << "ERROR: Aggregation kernel undefined!!";
			exit(1);
	}
}

void OpenSMOKE_PhysicalModels::breakageKernel(int kind, double parameterA, double parameterB, double L0, double coefficient)
{	
	int i;
	double dimension = coefficient*L0;

	switch (kind)
	{
		case 0:
			// Case 0: Constant (Marchisio 2003, p. 325)
			aKernel = parameterA;
		break;

		case 1:
			// Case 1: Power Law (Marchisio 2003, p. 325)
			for(i=1;i<=_N;i++)
				aKernel[i] = 
				     (csi[i]<=dimension) ? 0. : parameterA * pow(csi[i], parameterB);

		break;
		case 2:
			// Case 2: Exponential Law (Marchisio 2003, p. 325)
			for(i=1;i<=_N;i++)
				aKernel[i] = 
				   (csi[i]<=dimension) ? 0 : parameterA * exp(parameterB*powerOfCsi[4][i]);
		break;
		
		default:
			cout << "ERROR: Breakage Kernel undefined!!";
			exit(1);
	}
}

void OpenSMOKE_PhysicalModels::fragmentationKernel(int kind, double parameter)
{	
	int i, k;
	double coefficientFragment;

	switch (kind)
	{
		case 1:
			// Case 1: Symmetric Fragmentation (Marchisio 2003, p. 325)
			for(k=0;k<=_2N_minus_1;k++)
			{
				coefficientFragment = pow(2., (3.-k)/3.);
				for(i=1;i<=_N;i++)
					bKernel[k+1][i] = coefficientFragment * powerOfCsi[k+1][i];
			}
		break;

		case 2:
			// Case 2: Mass Ratio 1:4 (Marchisio 2003, p. 325)
			for(k=0;k<=_2N_minus_1;k++)
			{
				coefficientFragment = ( pow(4.,k/3.) + 1.) / pow(5.,k/3.);
				for(i=1;i<=_N;i++)
					bKernel[k+1][i] = coefficientFragment * powerOfCsi[k+1][i];
			}
		break;

		case 3:
			// Case 3: Uniform (Marchisio 2003, p. 325)
			for(k=0;k<=_2N_minus_1;k++)
			{
				coefficientFragment = 6./(k+3.);
				for(i=1;i<=_N;i++)
					bKernel[k+1][i] = coefficientFragment * powerOfCsi[k+1][i];
			}
		break;
		
		case 4:
			// Case 4: Parabolic (Marchisio 2003, p. 325)
			for(k=0;k<=_2N_minus_1;k++)
			{
				coefficientFragment = 3.*parameter/(3+k) + (1.-0.50*parameter) * 
					                  (72./(9.+k) - 72./(6.+k) + 18./(3.+k));
					
				for(i=1;i<=_N;i++)
					bKernel[k+1][i] = coefficientFragment * powerOfCsi[k+1][i];
			}
		break;

		case 5:
			// Case 5: Erosion (Marchisio 2003, p. 325)
			for(k=0;k<=_2N_minus_1;k++)
				for(i=1;i<=_N;i++)
				{
					double a = powerOfCsi[4][i]-1.;
					if (a>=0.)
						bKernel[k+1][i] = 1. + pow(a, k/3.);
					else
						bKernel[k+1][i] = 1. - pow(-a, k/3.);
				}
		break;

		default:
			cout << "ERROR: Fragment Function undefined!!";
			exit(1);
	}	
}

void OpenSMOKE_PhysicalModels::sources()
{
	// ---------------------------------------------------------------
	// Parameters
	// ---------------------------------------------------------------	
	double L0 = 1.;
	BzzVector L(_N);
	L[1] = 0.2113;
	L[2] = 0.7887;

	// ---------------------------------------------------------------
	// Sources
	// ---------------------------------------------------------------
	int iCase	= 0;
	int iFigure = 0;
	int iMcGraw = 0;
	int iCaseFLUENT = 0;
	int iAcetylene = 0;
	int iMoss = 0;
	
	if (iMcGraw==1)						// Figure 2
	{
		double kappa = 0.78;
		double rGrowth = -1.;

		powerLawGrowth(kappa, rGrowth);
	}
	
	if (iFigure==1)						// Figure 1
	{
		double kappa = 1.;
		double rGrowth = 1.;

		powerLawGrowth(kappa, rGrowth);
	}
	
	if (iFigure==2)						// Figure 2
	{
		BzzVector D(_N); D = 1.;
		homogeneousDispersion(D);
	}

	if (iFigure==3)						// Figure 3-5
	{
		double breakageConstant = 0.02;

		aggregationKernel(0, 1.);								// Constant
		breakageKernel(0, breakageConstant, 0., L0, 1.);		// Constant
		fragmentationKernel(1, 0.);								// Symmetric
		homogeneousAggregation();
		homogeneousBreakage();
	}

	if (iFigure==4)						// Figure 4-6
	{
		double breakageConstant = 0.02;

		aggregationKernel(2, 0.);							// Hydrodynamic
		breakageKernel(1, breakageConstant, 3., L0, 1);		// Power-Law
		fragmentationKernel(1, 0.);							// Symmetric
		homogeneousAggregation();
		homogeneousBreakage();
	}
	
	if (iCase==1)						// Case 1 - Vanni: ?-OK
	{
		double breakageConstant = 0.02;

		aggregationKernel(0, 1.);								// Constant
		breakageKernel(0, breakageConstant, 0., L0, 1.);		// Constant
		fragmentationKernel(1, 0.);								// Symmetric	
		homogeneousAggregation();
		homogeneousBreakage();
	}

	if (iCase==2)						// Case 2 - Vanni: OK-OK
	{
		double breakageConstant = 0.02;

		aggregationKernel(2, 0.);								// Sum
		breakageKernel(1, breakageConstant, 3., L0, 1.);		// Power-Law
		fragmentationKernel(1, 0.);								// Symmetric
		homogeneousAggregation();
		homogeneousBreakage();
	}

	if (iCase==3)						// Case 3 - Vanni: OK-OK
	{
		double breakageConstant = 0.1;

		aggregationKernel(1, 0.);							// Hydrodynamic
		breakageKernel(2, breakageConstant, 0.01, L0, 1.);	// Exponential
		fragmentationKernel(1, 0.);							// Symmetric
		homogeneousAggregation();
		homogeneousBreakage();
	}

	if (iCase==4)						// Case 4 - Vanni: OK-OK
	{
		double breakageConstant = 0.1;

		aggregationKernel(1, 0.);										// Hydrodynamic
		breakageKernel(2, breakageConstant, 0.01, L0, pow(5., .333));	// Exponential
		fragmentationKernel(2, 0.);										// Mass ratio 1:4
		homogeneousAggregation();
		homogeneousBreakage();
	}

	if (iCase==5)						// Case 5 - Vanni: OK-OK
	{
		double breakageConstant = 1.0;

		aggregationKernel(1, 0.);											// Hydrodynamic
		breakageKernel(2, breakageConstant, 0.01, L0, pow(3., .333));		// Exponential
		fragmentationKernel(5, 0.);											// Erosion
		homogeneousAggregation();
		homogeneousBreakage();
	}

	if (iCase==6)						// Case 6 - Vanni: OK-OK
	{
		double breakageConstant = 0.1;

		aggregationKernel(1, 0.);								// Hydrodynamic
		breakageKernel(2, breakageConstant, 0.01, L0, 1.);		// Exponential
		fragmentationKernel(3, 0.);								// Uniform
		homogeneousAggregation();
		homogeneousBreakage();
	}

	if (iCase==8)						// Case 8 - Vanni: OK-OK
	{
		double breakageConstant = 0.01;

		aggregationKernel(1, 0.);										// Hydrodynamic
		breakageKernel(1, breakageConstant, 6., L0, pow(3., .333));		// Power-Law
		fragmentationKernel(5, 0.);										// Erosion
		homogeneousAggregation();
		homogeneousBreakage();
	}

	if (iCase==9)						// Case 9 - Vanni: OK-OK
	{
		double breakageConstant = 2.0;

		aggregationKernel(1, 0.);								// Hydrodynamic
		breakageKernel(1, breakageConstant, 1.5, L0, 1.0);		// Pawer-Law
		fragmentationKernel(4, 0.5);							// Parabolic C=0.50
		homogeneousAggregation();
		homogeneousBreakage();
	}

	if (iCase==10)						// Case 10 - Vanni: OK-OK
	{
		double breakageConstant = 0.01;

		aggregationKernel(4, 0.);								// Differential force
		breakageKernel(1, breakageConstant, 6., L0, 1.);		// Power-Law
		fragmentationKernel(3, 0.);								// Uniform
		homogeneousAggregation();
		homogeneousBreakage();
	}

	if (iCaseFLUENT==11)						// Case 1.1
	{
		double kappa = 4.00;
		double rGrowth = -1.;

		powerLawGrowth(kappa, rGrowth);
	}

	if (iCaseFLUENT==21)						// Case 2.1
	{
		BzzVector D(_N); D = 1.;
		homogeneousDispersion(D);
	}

	if (iCaseFLUENT==31)						// Case 3.1 (Vanni case 1)
	{
		double breakageConstant = 0.02;
		double L0 = 1.;
		aggregationKernel(0, 1.);								// Constant
		breakageKernel(0, breakageConstant, 0., L0, 1.);		// Constant
		fragmentationKernel(1, 0.);								// Symmetric	
		homogeneousAggregation();
		homogeneousBreakage();
	}

	if (iCaseFLUENT==32)						// Case 3.1 (Vanni case 2)
	{
		double breakageConstant = 0.02;
		double L0 =1.;

		aggregationKernel(2, 0.);								// Sum
		breakageKernel(1, breakageConstant, 3., L0, 1.);		// Power-Law
		fragmentationKernel(1, 0.);								// Symmetric
		homogeneousAggregation();
		homogeneousBreakage();
	}

	if (iCaseFLUENT==33)						// Case 3.3
	{
		aggregationKernel(0, 1.);								// Constant
		homogeneousAggregation();
	}

	if (iCaseFLUENT==34)						// Case 3.4
	{
		double breakageConstant = 2.0;
		double L0 = 1.0;

		breakageKernel(0, breakageConstant, 0., L0, 1.);		// Constant
		fragmentationKernel(1, 0.);								// Symmetric	
		homogeneousBreakage();
	}

	if (iCaseFLUENT==41)						// Only nucleation
	{
		if (_N!=2)
		{
			cout << "Error: Nucleation just for N=2" << endl;
			exit(1);
		}

		double epsilon = 1.;
		double J = 1.;
		nucleation(epsilon, J, L);
	}

	if (iAcetylene==1)						// Only nucleation
	{
		if (_N!=2)
		{
			cout << "Error: Nucleation just for N=2" << endl;
			exit(1);
		}

		double epsilon = 1.;
		double J = 1.;
		nucleation(epsilon, J, L);
	}

	if (iMoss==1)						// Figure 2
	{
		double beta = 1.;
		double kappa = 2*beta/pow(3.14159,1./3.)/pow(6.,2./3.)/pow(1800.,1./3.);
		double rGrowth = 0.;

		powerLawGrowth(kappa, rGrowth);
	}
}

void OpenSMOKE_PhysicalModels::aggregationKernel_Soot(BzzMatrix &Kernel)
{
	betaKernel = Kernel;
}
