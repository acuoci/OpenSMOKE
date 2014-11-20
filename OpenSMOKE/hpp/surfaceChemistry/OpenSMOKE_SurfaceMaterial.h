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

#ifndef OPENSMOKE_SURFACEMATERIAL_H
#define OPENSMOKE_SURFACEMATERIAL_H

#include "BzzMath.hpp"

class OpenSMOKE_SurfaceSite;
class OpenSMOKE_SurfaceBulk;
class OpenSMOKE_SurfaceKinetics;
class OpenSMOKE_ReactingSurface;

class OpenSMOKE_SurfaceMaterial
{
public:

	OpenSMOKE_SurfaceMaterial();

	void Setup(OpenSMOKE_ReactingSurface *surface, const std::string name_);
	void ReadFromBinaryFile(BzzLoad &binaryFile);
	void ReadKineticsFromBinaryFile(BzzLoad &binaryFile);
	void ReactionRates( const double T, const BzzVector &cGas, BzzVector &cSurface, BzzVector &aBulk,
						BzzVector &RGas, BzzVector &RSurface, BzzVector &RBulk, double &QReaction);

	inline std::string name()	   { return name_object;  }
	inline int NumberOfSites() { return numberSites_; }
	inline int NumberOfBulks() { return numberBulks_; }

	inline OpenSMOKE_SurfaceSite*		sites()		{return sites_; };
	inline OpenSMOKE_SurfaceBulk*		bulks()		{return bulks_; };
	inline OpenSMOKE_SurfaceKinetics*	kinetics()	{return kinetics_; };


private:
	
	OpenSMOKE_ReactingSurface *ptSurface;
	OpenSMOKE_SurfaceSite     *sites_;
	OpenSMOKE_SurfaceBulk	  *bulks_;
	OpenSMOKE_SurfaceKinetics *kinetics_;

	int numberSites_;
	int numberBulks_;

	// Matrices
	BzzMatrix site_occupancy_matrix_;
	BzzMatrix bulk_density_matrix_;

private:

	std::string name_object;
	void ErrorMessage(const std::string message);
	void WarningMessage(const std::string message);
};

#endif	// OPENSMOKE_SURFACEMATERIAL_H


