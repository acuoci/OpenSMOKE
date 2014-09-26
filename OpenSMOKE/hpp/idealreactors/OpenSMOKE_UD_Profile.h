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

#ifndef OPENSMOKE_UD_PROFILE_H
#define OPENSMOKE_UD_PROFILE_H

#include "BzzMath.hpp"
#include "basic/OpenSMOKE_Utilities.h"

class OpenSMOKE_UD_Profile
{
friend class OpenSMOKE_PFR_Geometry;
friend class OpenSMOKE_CSTR_Geometry;
friend class OpenSMOKE_ShockTube_Geometry;

public:

    OpenSMOKE_UD_Profile();

    void    AssignFromFile(const std::string fileName, const std::string name);
    void    SetName(const std::string name);
    void	Check(const double abscissa, const double expected_value);
    double  GiveMeValue(const double t, const double x);
	void Setup(const std::string units_x, const std::string units_y, const double x0, const double x1, const double y0, const double y1);


public:

	std::string				name_object;
    bool				timeSupport;
	bool				iLinear;
    LinearInterpolation interpolation;

	double				yZero;
	double				xMax;
	double				slope;

	BzzVector		x;
    BzzVector		y;

	bool CheckForComment(const std::string read_word);
	bool CheckForKeyWord(const std::string read_word);

    void ErrorMessage(const std::string message);
};

#endif // OPENSMOKE_UD_PROFILE_H
