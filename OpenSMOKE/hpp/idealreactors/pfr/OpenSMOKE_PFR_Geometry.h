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


#ifndef OPENSMOKE_PFR_GEOMETRY_H
#define OPENSMOKE_PFR_GEOMETRY_H

#include "BzzMath.hpp"
#include "basic/OpenSMOKE_Utilities.h"
#include "idealreactors/OpenSMOKE_UD_Profile.h"

class OpenSMOKE_PFR_Geometry
{
public:

    OpenSMOKE_PFR_Geometry();

    void Setup(const std::string fileName);
    void SetName(const std::string name);
	void Update(const double x, double &D, double &Area);

private:

    static const double  pi_over_4;
    static const double _4_over_pi;

	std::string  name_object;
    int     iKind;
    double  D0;
    double  Area0;

    double slope_diameter;
    double slope_area;

	OpenSMOKE_UD_Profile profile;

    LinearInterpolation interpolation_diameter;
    LinearInterpolation interpolation_area;
    LinearInterpolation interpolation_dArea;
    LinearInterpolation interpolation_dDiameter;

    void ErrorMessage(const std::string message);
};

#endif // OPENSMOKE_PFR_GEOMETRY_H
