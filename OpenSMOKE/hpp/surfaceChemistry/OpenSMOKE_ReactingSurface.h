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

#ifndef OPENSMOKE_REACTINGSURFACE_H
#define OPENSMOKE_REACTINGSURFACE_H

#include "BzzMath.hpp"
#include <vector>

class OpenSMOKE_SurfaceKinetics;
class OpenSMOKE_SurfaceMaterial;
class OpenSMOKE_SurfaceThermodynamics;

enum version_surface_file { Surface110514 };

class OpenSMOKE_ReactingSurface
{
public:

	void ReadFromBinaryFile(const string fileName);

	inline int NumberMaterials()   { return numberMaterials_; }
	inline int NumberSiteSpecies() { return numberSiteSpecies_; }
	inline int NumberBulkSpecies() { return numberBulkSpecies_; }

	inline string NameSiteSpecies(const int i)	{ return namesSiteSpecies_[i]; }
	inline string NameBulkSpecies(const int i)	{ return namesBulkSpecies_[i]; }

	inline OpenSMOKE_SurfaceMaterial*			material()		 { return material_; };
	inline OpenSMOKE_SurfaceThermodynamics*		thermodynamics() { return thermodynamics_; };

private:

	OpenSMOKE_SurfaceMaterial		*material_;
	OpenSMOKE_SurfaceThermodynamics	*thermodynamics_;

	int numberMaterials_;
	int numberSiteSpecies_;
	int numberBulkSpecies_;

	vector<string> namesSiteSpecies_;
	vector<string> namesBulkSpecies_;

	// Elements
	BzzMatrix site_elements_;
	BzzVector site_m_elements_;
	BzzMatrix bulk_elements_;
	BzzVector bulk_m_elements_;

	void SummaryOnVideo();

private:

	version_surface_file iVersion;

	string name_object;
	void ErrorMessage(const string message);
	void WarningMessage(const string message);
};

#endif	// OPENSMOKE_REACTINGSURFACE_H


