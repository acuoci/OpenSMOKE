/***************************************************************************
 *   Copyright (C) 2006-2008 by Alberto Cuoci							   *
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


#ifndef OPENSMOKE_CSTR_GEOMETRY_H
#define OPENSMOKE_CSTR_GEOMETRY_H

#include "BzzMath.hpp"
#include "basic/OpenSMOKE_Utilities.h"
#include "idealreactors/OpenSMOKE_UD_Profile.h"

class OpenSMOKE_CSTR_Geometry
{
public:

    OpenSMOKE_CSTR_Geometry();

    void Setup(const string fileName);
    void SetName(const string name);
	void Update(const double x, double &D, double &Volume, double &Area);
	void Update(const double x, double &Tau);

private:

    static const double  pi_over_6;
    static const double _6_over_pi;

	string  name_object;
    int     iKind;
    
	double  D0;
    double  Tau0;
    double  Volume0;

    double slope_diameter;
    double slope_volume;
    double slope_tau;

	OpenSMOKE_UD_Profile profile;

    LinearInterpolation interpolation_diameter;
    LinearInterpolation interpolation_volume;
    LinearInterpolation interpolation_tau;

    LinearInterpolation interpolation_dVolume;
    LinearInterpolation interpolation_dDiameter;
    LinearInterpolation interpolation_dTau;

    void ErrorMessage(const string message);
};

#endif // OPENSMOKE_CSTR_GEOMETRY_H
