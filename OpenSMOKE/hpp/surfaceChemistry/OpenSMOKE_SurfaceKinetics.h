/***************************************************************************
 *   Copyright (C) 2011 by Alberto Cuoci   	                               *
 *   alberto.cuoci@polimi.it   						                       *
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

#ifndef OPENSMOKE_SURFACEKINETICS_H
#define OPENSMOKE_SURFACEKINETICS_H

#include "BzzMath.hpp"
#include "surfaceChemistry/OpenSMOKE_SurfaceReaction.h"
#include "basic/OpenSMOKE_Constants.h"

class OpenSMOKE_ReactingSurface;


class OpenSMOKE_SurfaceKinetics
{
public:

	OpenSMOKE_SurfaceKinetics();
	void Setup(OpenSMOKE_ReactingSurface *surface, const string name_);
	void ReadFromBinaryFile(BzzLoad &binaryFile);
	void ReactionRates(	const double T, const double lnT, const BzzVector &cGas, BzzVector &cSurface, BzzVector &zSurface, BzzVector &aBulk,
				BzzVector &RGas, BzzVector &RSurface, BzzVector &RBulk, double &QReaction);

	inline int NumberOfReactions() { return numberReactions_; }
	inline OpenSMOKE_SurfaceReaction *reaction() { return reaction_; }
	inline double ActivationEnergy(const int i) const { return reaction_[i].E_over_R*Constants::R_J_kmol; }		// [J/kmol]
	inline double ActivationTemperature(const int i) const { return reaction_[i].E_over_R; }					// [K]

private:

	void Allocate();

private:

	OpenSMOKE_ReactingSurface *ptSurface;
	OpenSMOKE_SurfaceReaction *reaction_;

	int numberReactions_;
	BzzVector dhrs;

private:

	string name_object;
    void coupleDirectAndInverseReactions();
	void ErrorMessage(const string message);
	void WarningMessage(const string message);
};	

#endif	// OPENSMOKE_SURFACEKINETICS_H


