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

#if !defined(OPENSMOKE_CLIPPEDGAUSSIANDISTRIBUTION_H)
#define OPENSMOKE_CLIPPEDGAUSSIANDISTRIBUTION_H

#include "BzzMath.hpp"
#include "basic/OpenSMOKE_Utilities.h"

class OpenSMOKE_ClippedGaussianDistribution 
{
public:
	OpenSMOKE_ClippedGaussianDistribution();
	virtual ~OpenSMOKE_ClippedGaussianDistribution();

	void Set(const double _csi, const double _csiV2);

	double FlatIntegral();
	double IntegralNormalized();
	double FlatCentralIntegral();

	double ReactionCorrectionCoefficient(const double Tatt, const double Tmean, const double Tmin, const double Tmax);
	double ReactionCorrectionCoefficient(const double n, const double Tatt, const double Tmean, const double Tmin, const double Tmax);
	double GetValue(LinearInterpolation &interpolation);

private:

	void ErrorMessage(const string message);
	void WarningMessage(const string message);

private:

	double csi;
	double csiV2;

	double alfa1;
	double alfa2;
	double x1;
	double x2;
	double g;

	int N;
	double dx;
	BzzVector x;
	BzzVector y;
	BzzVector f;
	double f_a;
	double f_b;
};

#endif // !defined(OPENSMOKE_CLIPPEDGAUSSIANDISTRIBUTION)
