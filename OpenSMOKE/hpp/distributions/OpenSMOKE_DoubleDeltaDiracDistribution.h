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

#if !defined(OPENSMOKE_DOUBLEDELTADIRACDISTRIBUTION_H)
#define OPENSMOKE_DOUBLEDELTADIRACDISTRIBUTION_H

#include "BzzMath.hpp"

class OpenSMOKE_DoubleDeltaDiracDistribution 
{
public:

	OpenSMOKE_DoubleDeltaDiracDistribution();
	virtual ~OpenSMOKE_DoubleDeltaDiracDistribution();

	void Set(const double _csi, const double _csiV2);

	double ReactionCorrectionCoefficient(const double Tatt, const double Tmean, const double Tmin, const double Tmax);
	double ReactionCorrectionCoefficient(const double n, const double Tatt, const double Tmean, const double Tmin, const double Tmax);

private:

	void ErrorMessage(const std::string message);
	void WarningMessage(const std::string message);

private:

	double csi;
	double g;
	double fMinus;
	double fPlus;
	double wMinus;
	double wPlus;
};

#endif // !defined(OPENSMOKE_DOUBLEDELTADIRACDISTRIBUTION_H)
