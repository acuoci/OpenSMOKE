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

#include "qmom/OpenSMOKE_Distributions.h"
using namespace std;

void OpenSMOKE_Distributions::printOnVideo()
{
	cout << endl;
	cout << "Order\t" << "moment" << endl;
	cout << "--------------------------------" << endl;
	for(int i=1;i<=moments.Size();i++)
		cout << i-1 << "\t" << moments[i] << endl;
	cout << endl;
}

double OpenSMOKE_Distributions::variance()
{
	return moments[3] - moments[2]*moments[2];
}

// ------------------------------------------------------------------------//
//						   1. DIRAC DISTRIBUTION						   //
// ------------------------------------------------------------------------//

void OpenSMOKE_DiracDistribution::setup()
{
}
void OpenSMOKE_DiracDistribution::update(double _csiN)
{
	csiN = _csiN;
}

void OpenSMOKE_DiracDistribution::control()
{
}

void OpenSMOKE_DiracDistribution::GetMoment(int orderOfMoment, double &m)
{
	m = pow(csiN, orderOfMoment);
}

void OpenSMOKE_DiracDistribution::GetMomentOrderOne(double &m)
{
	m = csiN;
}

void OpenSMOKE_DiracDistribution::GetMomentOrderTwo(double &m)
{
	m = csiN*csiN;
}

void OpenSMOKE_DiracDistribution::GetMomentOrderThree(double &m)
{
	m = csiN*csiN*csiN;
}

void OpenSMOKE_DiracDistribution::GetMomentOrderFour(double &m)
{
	m = csiN*csiN*csiN*csiN;
}

void OpenSMOKE_DiracDistribution::GetMomentRatioThreeOverTwo(double &d32)
{
	d32 = csiN;
}

void OpenSMOKE_DiracDistribution::GetMomentRatioFourOverThree(double &d43)
{
	d43 = csiN;
}

void OpenSMOKE_DiracDistribution::GetMoments(int maxOrderOfMoments)
{
	int i;

	ChangeDimensions(maxOrderOfMoments+1, &moments);

	moments[1] = 1.;							// order 0
	for(i=1;i<=maxOrderOfMoments;i++)			// order 1... N
		moments[i+1] = moments[i]*csiN;
}


// ------------------------------------------------------------------------//
//					 2. MULTI DIRAC DISTRIBUTION						   //
// ------------------------------------------------------------------------//

void OpenSMOKE_MultiDiracDistribution::setup(int _N)
{
	N = _N;

	epsilon = 1.e-4;
	_1_plus_epsilon = 1.+epsilon;
	_1_minus_epsilon = 1.-epsilon;
}

void OpenSMOKE_MultiDiracDistribution::update(BzzVector &_csiN, BzzVector &_pN)
{
	csiN = _csiN;
	pN = _pN;
//	control();
}

void OpenSMOKE_MultiDiracDistribution::updateWithoutControl(BzzVector &_csiN, BzzVector &_pN)
{
	csiN = _csiN;
	pN = _pN;
}

void OpenSMOKE_MultiDiracDistribution::control()
{
	// Error messages
	// ---------------------------------------------------------------------
	double sum = pN.GetSumElements(); 

	if (sum<=_1_minus_epsilon || sum>=_1_plus_epsilon)
	{
		cout << "ERROR! MultiDiracDistribution!" << endl;
		cout << "The sum of probabilities (" << sum 
			 << ") is different than 1." << endl;
		exit(1);
	}

	if (csiN.Size()!=pN.Size())
	{
		cout << "ERROR! MultiDiracDistribution!" << endl;
		cout << "The dimensions of abscissas and probabilities must agree!" << endl;
		exit(1);
	}
}

void OpenSMOKE_MultiDiracDistribution::GetMoment(int numberOfMoment, double &m)
{
	m=0.;
	for(int j=1;j<=N;j++)
		m += pow(csiN[j],numberOfMoment)*pN[j];
}

void OpenSMOKE_MultiDiracDistribution::GetMomentOrderOne(double &m)
{
	m=0.;
	for(int j=1;j<=N;j++)
		m += csiN[j]*pN[j];
}

void OpenSMOKE_MultiDiracDistribution::GetMomentOrderTwo(double &m)
{
	m=0.;
	for(int j=1;j<=N;j++)
		m += csiN[j]*csiN[j]*pN[j];
}

void OpenSMOKE_MultiDiracDistribution::GetMomentOrderThree(double &m)
{
	m=0.;
	for(int j=1;j<=N;j++)
		m += csiN[j]*csiN[j]*csiN[j]*pN[j];
}

void OpenSMOKE_MultiDiracDistribution::GetMomentOrderFour(double &m)
{
	m=0.;
	for(int j=1;j<=N;j++)
		m += csiN[j]*csiN[j]*csiN[j]*csiN[j]*pN[j];
}

void OpenSMOKE_MultiDiracDistribution::GetMomentRatioThreeOverTwo(double &d32)
{
	double a, m;

	m = 0.;
	d32 = 0.;
	for(int j=1;j<=N;j++)
	{
		a = csiN[j]*csiN[j]*pN[j];
		m += a;
		d32 += a*csiN[j];
	}
	d32 /= m; 
}

void OpenSMOKE_MultiDiracDistribution::GetMomentRatioFourOverThree(double &d43)
{
	double a, m;

	m = 0.;
	d43 = 0.;
	for(int j=1;j<=N;j++)
	{
		a = csiN[j]*csiN[j]*csiN[j]*pN[j];
		m += a;
		d43 += m*csiN[j];
	}
	d43 /= m; 
}

void OpenSMOKE_MultiDiracDistribution::GetMoments(int maxOrderOfMoments)
{
	int i, j;
	BzzVector csiStar(N);
	csiStar = 1.;

	ChangeDimensions(maxOrderOfMoments+1, &moments);

	for(i=0;i<=maxOrderOfMoments;i++)			// order 0... N
	{
		for(j=1;j<=N;j++)
			moments[i+1] += csiStar[j]*pN[j];
		ElementByElementProduct(csiStar, csiN, &csiStar);
	}
}


// ------------------------------------------------------------------------//
//					 3. CONSTANT DISTRIBUTION						   //
// ------------------------------------------------------------------------//

void OpenSMOKE_ConstantDistribution::setup()
{
}

void OpenSMOKE_ConstantDistribution::update(double _csiA, double _csiB, double _H)
{
	csiA = _csiA;
	csiB = _csiB;
	H = _H;

	control();
}

void OpenSMOKE_ConstantDistribution::control()
{
	// Error messages
	// ---------------------------------------------------------------------
	if (csiB<=csiA)
	{
		cout << "ERROR! ConstantDistribution!" << endl;
		cout << "The extremes of interval of definition are wrong!" << endl;
		exit(1);
	}

	if (H*(csiB-csiA)!=1.)
	{
		cout << "ERROR! ConstantDistribution!" << endl;
		cout << "The distribution is not normal!" << endl;
		exit(1);
	}
}

void OpenSMOKE_ConstantDistribution::GetMoment(int numberOfMoment, double &m)
{
	m = H * (pow(csiB, numberOfMoment) - pow(csiA, numberOfMoment)) 
		                   / (numberOfMoment+1);
}

void OpenSMOKE_ConstantDistribution::GetMomentOrderOne(double &m)
{
	m = H * (csiB - csiA) / 2.;
}

void OpenSMOKE_ConstantDistribution::GetMomentOrderTwo(double &m)
{
	m = H * (csiB*csiB - csiA*csiA) / 3.; 
}

void OpenSMOKE_ConstantDistribution::GetMomentOrderThree(double &m)
{
	m = H * (csiB*csiB*csiB - csiA*csiA*csiA) / 4.; 
}

void OpenSMOKE_ConstantDistribution::GetMomentOrderFour(double &m)
{
	m = H * (csiB*csiB*csiB*csiB - csiA*csiA*csiA*csiA) / 5.; 
}

void OpenSMOKE_ConstantDistribution::GetMomentRatioThreeOverTwo(double &d32)
{
	double a = csiB*csiB;
	double b = csiA*csiA;
	d32 = 0.75 * (b*csiB - a*csiA) / (b - a); 
}

void OpenSMOKE_ConstantDistribution::GetMomentRatioFourOverThree(double &d43)
{
	double a = csiB*csiB*csiB;
	double b = csiA*csiA*csiA;
	d43 = 0.80 * (b*csiB - a*csiA) / (b - a);
}

void OpenSMOKE_ConstantDistribution::GetMoments(int maxOrderOfMoments)
{
	int i;
	double csiAStar = csiA;
	double csiBStar = csiB;

	ChangeDimensions(maxOrderOfMoments+1, &moments);

	for(i=0;i<=maxOrderOfMoments;i++)			// order 0... N
	{
		moments[i+1] = H * (csiBStar - csiAStar) 
		                   / (i+1);
		csiAStar*=csiA;
		csiBStar*=csiB;
	}
}

// ------------------------------------------------------------------------//
//					 3. KM CLOUD DISTRIBUTION							   //
// ------------------------------------------------------------------------//
int factorial(int n)
{
	int i;
	int sum;

	sum = 1;
	for(i=2;i<=n;i++) sum *= i;
	return sum;
}

void OpenSMOKE_KM_Cloud_Distribution::setup()
{
}

void OpenSMOKE_KM_Cloud_Distribution::update(double _a, double _b)
{
	a = _a;
	b = _b;

	control();
}

void OpenSMOKE_KM_Cloud_Distribution::control()
{
	// Error messages
	// ---------------------------------------------------------------------
	if (a<=0.)
	{
		cout << "ERROR! KM_Cloud_Distribution!" << endl;
		cout << "The parameters a and b must be positive!" << endl;
		exit(1);
	}
}

void OpenSMOKE_KM_Cloud_Distribution::GetMoment(int numberOfMoment, double &m)
{
	m = a * factorial(numberOfMoment+2) * pow(b, -3-numberOfMoment);
}

void OpenSMOKE_KM_Cloud_Distribution::GetMomentOrderOne(double &m)
{
	m = 6.*a/BzzPow4(b);
}

void OpenSMOKE_KM_Cloud_Distribution::GetMomentOrderTwo(double &m)
{
	m = 24.*a/BzzPow5(b);
}

void OpenSMOKE_KM_Cloud_Distribution::GetMomentOrderThree(double &m)
{
	m = 120.*a/BzzPow6(b);
}

void OpenSMOKE_KM_Cloud_Distribution::GetMomentOrderFour(double &m)
{
	m = 720.*a/BzzPow7(b);
}

void OpenSMOKE_KM_Cloud_Distribution::GetMomentRatioThreeOverTwo(double &d32)
{
	d32 = 5./b;
}

void OpenSMOKE_KM_Cloud_Distribution::GetMomentRatioFourOverThree(double &d43)
{
	d43 = 6/b;
}

void OpenSMOKE_KM_Cloud_Distribution::GetMoments(int maxOrderOfMoments)
{
	int i;

	ChangeDimensions(maxOrderOfMoments+1, &moments);

	for(i=0;i<=maxOrderOfMoments;i++)			// order 0... N
		moments[i+1] = a * factorial(i+2) * pow(b, -3-i);
}

