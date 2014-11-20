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

#if !defined(OPENSMOKE_GAMMAFUNCTION_H)
#define OPENSMOKE_GAMMAFUNCTION_H

#include "BzzMath.hpp"

class OpenSMOKE_GammaFunction  
{
public:
	OpenSMOKE_GammaFunction();
	virtual ~OpenSMOKE_GammaFunction();

	double at(const double y);

private:

	void ErrorMessage(const std::string message);
	void WarningMessage(const std::string message);

	double PolynomialExpansion(const double y);
	double RecursionFormula(const double y);

	static const double c1;
	static const double c2;
	static const double c3;
	static const double c4;
	static const double c5;
	static const double c6;
	static const double c7;
	static const double c8;
};

#endif // !defined(OPENSMOKE_GAMMA_FUNCTION)
