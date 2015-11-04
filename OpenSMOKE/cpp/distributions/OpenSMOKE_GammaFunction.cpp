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

#include <string>
#include "basic/OpenSMOKE_Constants.h"
#include "distributions/OpenSMOKE_GammaFunction.h"
using namespace std;

const double OpenSMOKE_GammaFunction::c1 = -0.577191652;
const double OpenSMOKE_GammaFunction::c2 =  0.988205891;
const double OpenSMOKE_GammaFunction::c3 = -0.897056937;
const double OpenSMOKE_GammaFunction::c4 =  0.918206857;
const double OpenSMOKE_GammaFunction::c5 = -0.756704078;
const double OpenSMOKE_GammaFunction::c6 =  0.482199394;
const double OpenSMOKE_GammaFunction::c7 = -0.193527818;
const double OpenSMOKE_GammaFunction::c8 =  0.035868343;

void OpenSMOKE_GammaFunction::ErrorMessage(const std::string message)
{
    cout << endl;
    cout << "Class:  OpenSMOKE_GammaFunction"			<< endl;
    cout << "Error:  " << message           << endl;
    cout << "Press a key to continue... "   << endl;
    getchar();
    exit(-1);
}

void OpenSMOKE_GammaFunction::WarningMessage(const std::string message)
{
    cout << endl;
    cout << "Class:    OpenSMOKE_GammaFunction"		<< endl;
    cout << "Warning:  "	<< message      << endl;
    cout << "Press a key to continue... "   << endl;
    getchar();
}

OpenSMOKE_GammaFunction::OpenSMOKE_GammaFunction()
{
}

OpenSMOKE_GammaFunction::~OpenSMOKE_GammaFunction()
{
}

double OpenSMOKE_GammaFunction::at(const double y)
{
	if (y<0.)
		ErrorMessage("The argument must be positive...");

	if (y>=1. && y<=2.)	
		return PolynomialExpansion(y);
	else if (y<1.)
		return PolynomialExpansion(y+1.)/y;
	else
		return RecursionFormula(y);
}

double OpenSMOKE_GammaFunction::PolynomialExpansion(const double y)
{
	double x = y-1.;
	return 1.+x*(c1+x*(c2+x*(c3+x*(c4+x*(c5+x*(c6+x*(c7+x*c8)))))));
}

double OpenSMOKE_GammaFunction::RecursionFormula(const double y)
{
	int Y = y/1;
	double z = y-Y;
	
	double product=1.;
	for(int i=1;i<=Y-1;i++)
		product*=(z+i);

	return PolynomialExpansion(z+1.)*product;
}


