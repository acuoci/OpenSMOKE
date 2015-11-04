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
#include "distributions/OpenSMOKE_DoubleDeltaDiracDistribution.h"
using namespace std;

void OpenSMOKE_DoubleDeltaDiracDistribution::ErrorMessage(const std::string message)
{
    cout << endl;
    cout << "Class:  OpenSMOKE_DoubleDeltaDiracDistribution"			<< endl;
    cout << "Error:  " << message           << endl;
    cout << "Press a key to continue... "   << endl;
    getchar();
    exit(-1);
}

void OpenSMOKE_DoubleDeltaDiracDistribution::WarningMessage(const std::string message)
{
    cout << endl;
    cout << "Class:    OpenSMOKE_DoubleDeltaDiracDistribution"		<< endl;
    cout << "Warning:  "	<< message      << endl;
    cout << "Press a key to continue... "   << endl;
    getchar();
}

OpenSMOKE_DoubleDeltaDiracDistribution::OpenSMOKE_DoubleDeltaDiracDistribution()
{
}

OpenSMOKE_DoubleDeltaDiracDistribution::~OpenSMOKE_DoubleDeltaDiracDistribution()
{
}

void OpenSMOKE_DoubleDeltaDiracDistribution::Set(const double _csi, const double _g)
{
	int iRegion;

	if (_csi<=0.)	ErrorMessage("The mean value must be stricly positive...");
	if (_g<=0.)		ErrorMessage("The variance value must be stricly positive...");

	double gmax = 0.99*_csi*(1.-_csi);
	if (_g > gmax)
	{
		if (_g > 1.25*gmax)
		{
			cout << "csi:    " << _csi << endl;
			cout << "g:      " << _g << endl;
			cout << "gmax:   " << 0.99*_csi*(1.-_csi) << endl;
			ErrorMessage("The variance is too large to be accepted...");
		}
	}

	csi = _csi;
	g	= _g;
	if (_g > gmax)
		g = gmax;

	double sqrt_g = sqrt(g);

	if		(csi <=0.50 && sqrt_g<=csi)		iRegion = 1;
	else if (csi >=0.50 && sqrt_g<=1.-csi)	iRegion = 1;
	else if (csi <=0.50 && sqrt_g>csi)		iRegion = 2;
	else if (csi >=0.50 && sqrt_g>1.-csi)	iRegion = 3;
	else									iRegion = 4;

	if (iRegion == 1)
	{
		fMinus = csi - sqrt_g;
		fPlus  = csi + sqrt_g;
		wMinus = 0.50;
		wPlus  = 0.50;
	}
	else if (iRegion == 2)
	{
		fMinus = 0.;
		fPlus  = csi + sqrt_g;
		wMinus = g/(csi*csi+g);
		wPlus  = csi/(csi+g/csi); 
	}
	else if (iRegion == 3)
	{
		fMinus = csi - sqrt_g;
		fPlus  = 1.;
		wMinus = (1.-csi)/(1.-csi+g/(1.-csi));
		wPlus  = g/((1.-csi)*(1.-csi)+g); 
	}
	else
	{
		fMinus = 0.;
		fPlus  = 1.;
		wMinus = 1.-csi;
		wPlus  = csi;
	}
}

double OpenSMOKE_DoubleDeltaDiracDistribution::ReactionCorrectionCoefficient(const double Tatt, const double Tmean, const double Tmin, const double Tmax)
{
	double TMinus = Tmin + fMinus*(Tmax-Tmin);
	double TPlus  = Tmin + fPlus*(Tmax-Tmin);
	
	return ( wMinus*exp(-Tatt/TMinus) + wPlus*exp(-Tatt/TPlus) )	/
							( exp(-Tatt/Tmean) );
}

double OpenSMOKE_DoubleDeltaDiracDistribution::ReactionCorrectionCoefficient(const double n, const double Tatt, const double Tmean, const double Tmin, const double Tmax)
{
	double TMinus = Tmin + fMinus*(Tmax-Tmin);
	double TPlus  = Tmin + fPlus*(Tmax-Tmin);
	
	return ( wMinus*pow(TMinus,n)*exp(-Tatt/TMinus) + wPlus*pow(TPlus,n)*exp(-Tatt/TPlus) )	/
							    ( pow(Tmean,n)*exp(-Tatt/Tmean) );
}


