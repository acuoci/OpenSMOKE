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

#if !defined(OPENSMOKE_BETADISTRIBUTION_H)
#define OPENSMOKE_BETADISTRIBUTION_H

#include "BzzMath.hpp"
#include "distributions/OpenSMOKE_GammaFunction.h"
#include "basic/OpenSMOKE_Utilities.h"

class OpenSMOKE_BetaFunction  
{
public:

	OpenSMOKE_BetaFunction();
	virtual ~OpenSMOKE_BetaFunction();

	double at(const double a, const double b);
	
	void SetIntegralFormulation(const double N);

private:

	bool iIntegralFormulation;

	void ErrorMessage(const string message);
	void WarningMessage(const string message);

private:
	
	OpenSMOKE_GammaFunction gamma;

private:

	int Nb;
	double dh;
	BzzVector xb;
	BzzVector uxb;
	BzzVector fb;

	double BrizuelaIntegral(const double a, const double b);
	void PrepareBrizuelaIntegral(const double N);
};

class OpenSMOKE_BetaDistribution 
{
public:
	OpenSMOKE_BetaDistribution();
	virtual ~OpenSMOKE_BetaDistribution();

	void Set(const double _a, const double _b);

	double FlatIntegral();
	double IntegralNormalized();

	double ReactionCorrectionCoefficient(const double Tatt, const double Tmean, const double Tmin, const double Tmax);
	double ReactionCorrectionCoefficient(const double n, const double Tatt, const double Tmean, const double Tmin, const double Tmax);

	double GetValue(LinearInterpolation &interpolation);

private:

	void ErrorMessage(const string message);
	void WarningMessage(const string message);

private:
	
	double a;
	double b;
	double BetaInf;
	double extreme_a;
	double extreme_b;
	OpenSMOKE_BetaFunction betaFunction;

	int Nsub;
	int N;
	int M;
	double epsilon;

	double dhcenter;
	BzzVector dhb;
	BzzVector dhf;

	BzzVector xb;
	BzzVector xf;
	BzzVector yb;
	BzzVector yf;

	BzzVector xcenter;
	BzzVector ycenter;

	double f_a;
	double f_b;
	BzzVector fb;
	BzzVector ff;
	BzzVector fcenter;
};

#endif // !defined(OPENSMOKE_BETADISTRIBUTION_H)
