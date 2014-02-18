/***************************************************************************
 *   Copyright (C) 2006-2008 by Alberto Cuoci   	                       *
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

#ifndef OPENSMOKE_INVERSEKINETICS
#define OPENSMOKE_INVERSEKINETICS

#include "BzzMath.hpp"

class OpenSMOKE_ReactingGas;
class OpenSMOKE_GlobalKinetics;


class OpenSMOKE_InverseKinetics
{
private:
	OpenSMOKE_ReactingGas		*mix;
	OpenSMOKE_GlobalKinetics	*global;

public:

	void AssignKineticScheme(OpenSMOKE_ReactingGas &_mix);
	void AssignGlobalKineticScheme(OpenSMOKE_GlobalKinetics &_global);
	void Setup(const int iReaction, double Patm, double Tmin, double Tmax, double deltaT);
	void Setup(const int iReaction, const double Patm, const double Tmin, const double Tmax, const double deltaT,
    		   double &A, double &Beta, double &Tatt, BzzVector &lambdaInverse);

	BzzVector T;
	BzzVector kappa;
	BzzVector Keq;
	BzzVector kappainv;
	BzzVector logkappainv;
	int		Npoints;
	double	P;
};

#endif // OPENSMOKE_INVERSEKINETICS