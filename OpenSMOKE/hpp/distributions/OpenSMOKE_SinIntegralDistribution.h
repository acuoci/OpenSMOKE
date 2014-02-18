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

#if !defined(OPENSMOKE_SININTEGRALDISTRIBUTION_H)
#define OPENSMOKE_SININTEGRALDISTRIBUTION_H

#include "BzzMath.hpp"

class OpenSMOKE_SinIntegralDistributionMyNonLinearSystem;

class OpenSMOKE_SinIntegralDistribution 
{
friend class OpenSMOKE_SinIntegralDistributionMyNonLinearSystem;

public:

	OpenSMOKE_SinIntegralDistribution();
	virtual ~OpenSMOKE_SinIntegralDistribution();

	void SetBoundaries(const double _Tmin, const double _Tmax);
	void Set(const double csi, const double g);
	double CorrectionCoefficient(const double n, const double Tatt, const double T);

private:

	void ErrorMessage(const string message);
	void WarningMessage(const string message);

private:

	int Npoints;
	double dx;
	BzzVector x;
	BzzVector sinx;
	BzzVector f;

	double Tmin;
	double Tmax;
	double dT;

	double alfa;
	double beta;
	double dTsin;

	double Tsigned;
	double sigmaT2;
};

class OpenSMOKE_SinIntegralDistributionMyNonLinearSystem: public BzzMyNonLinearSystemObject
{
public:
	void AssignSinIntegralDistribution(OpenSMOKE_SinIntegralDistribution *sin);
	OpenSMOKE_SinIntegralDistribution *ptSin;

public:
	virtual void GetResiduals(BzzVector &x, BzzVector &f);
	virtual void ObjectBzzPrint(void) {};
};

#endif // !defined(OPENSMOKE_SININTEGRALDISTRIBUTION_H)