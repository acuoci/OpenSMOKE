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

#if !defined(OPENSMOKE_SINDISTRIBUTION_H)
#define OPENSMOKE_SINDISTRIBUTION_H

#include "BzzMath.hpp"

class OpenSMOKE_SinDistribution 
{
public:

	OpenSMOKE_SinDistribution();
	virtual ~OpenSMOKE_SinDistribution();

	double CorrectionCoefficient(const double qT2, const double EsuR, const double n, const double T);

private:

	void ErrorMessage(const string message);
	void WarningMessage(const string message);

private:

	int iSinExpansion;

	static const double c2;
	static const double c4;
	static const double c6;
	static const double c8;
	static const double c10;
	static const double c12;
	static const double c14;
	static const double c16;
};

#endif // !defined(OPENSMOKE_SINDISTRIBUTION_H)
