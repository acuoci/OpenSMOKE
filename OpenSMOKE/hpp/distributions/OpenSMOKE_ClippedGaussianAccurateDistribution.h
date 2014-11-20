/***************************************************************************
 *   Copyright (C) 2007 by Alberto Cuoci								   *
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

#if !defined(OPENSMOKE_CLIPPEDGAUSSIANACCURATEDISTRIBUTION_H)
#define OPENSMOKE_CLIPPEDGAUSSIANACCURATEDISTRIBUTION_H

#include "BzzMath.hpp"
#include "basic/OpenSMOKE_Utilities.h"

class OpenSMOKE_ClippeGaussianAccurateDistributionMyNonLinearSystem;

class OpenSMOKE_ClippedGaussianAccurateDistribution 
{

friend class OpenSMOKE_ClippedGaussianAccurateDistributionMyNonLinearSystem;

public:
	OpenSMOKE_ClippedGaussianAccurateDistribution();
	virtual ~OpenSMOKE_ClippedGaussianAccurateDistribution();

	void Set(const double _csi, const double _csiV2);
	double IntegralNormalized();

	double ReactionCorrectionCoefficient(const double Tatt, const double Tmean, const double Tmin, const double Tmax);
	double ReactionCorrectionCoefficient(const double n, const double Tatt, const double Tmean, const double Tmin, const double Tmax);
	double GetValue(LinearInterpolation &interpolation);

private:

	void ErrorMessage(const std::string message);
	void WarningMessage(const std::string message);

private:

	double mu;
	double sigma;
	double csi;
	double g;

	void Prepare_yFunction(const double mu, const double sigma);
	double Integral1(const double mu, const double sigma);
	double Integral2(const double mu, const double sigma);
	double CentralIntegral(const double sigma);

	bool iDeltaDirac;

	int N;
	double dx;
	BzzVector x;
	BzzVector y;
	BzzVector f1;
	BzzVector f2;
	BzzVector f;
	double alfa1;
	double alfa2;
	double f_a;
	double f_b;

	static const double sqrt_2pi;
	static const double sqrt_2;
	static const double sqrt_pi;
};

class OpenSMOKE_ClippedGaussianAccurateDistributionMyNonLinearSystem: public BzzMyNonLinearSystemObject
{
public:
	void AssignGaussian(OpenSMOKE_ClippedGaussianAccurateDistribution *gaussian);
	OpenSMOKE_ClippedGaussianAccurateDistribution *ptGaussian;

public:
	virtual void GetResiduals(BzzVector &x, BzzVector &f);
	virtual void ObjectBzzPrint(void) {};
};

#endif // !defined(OPENSMOKE_CLIPPEDGAUSSIANACCURATEDISTRIBUTION_H)
