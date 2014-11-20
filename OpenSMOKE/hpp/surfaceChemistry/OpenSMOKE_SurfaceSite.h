/***************************************************************************
 *   Copyright (C) 2011 by Alberto Cuoci   	   *
 *   alberto.cuoci@polimi.it   						   *
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

#ifndef OPENSMOKE_SURFACESITE_H
#define OPENSMOKE_SURFACESITE_H

#include "BzzMath.hpp"

class OpenSMOKE_SurfaceSite
{
public:

	OpenSMOKE_SurfaceSite();
	void Setup(OpenSMOKE_SurfaceMaterial* material, const std::string name_object, const double density);

	inline std::string name()	{ return name_object; }
	inline double density() { return density_; }		// [kmol/m2]

private:
	
	OpenSMOKE_SurfaceMaterial* ptMaterial;
	double density_;									// [kmol/m2]

private:

	std::string name_object;
	void ErrorMessage(const std::string message);
	void WarningMessage(const std::string message);
};

#endif	// OPENSMOKE_SURFACESITE_H