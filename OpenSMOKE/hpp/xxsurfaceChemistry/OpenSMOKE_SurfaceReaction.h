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

#ifndef OPENSMOKE_SURFACEREACTION_H
#define OPENSMOKE_SURFACEREACTION_H

#include "BzzMath.hpp"

class OpenSMOKE_SurfaceReaction
{
public:

	void ReadFromBinaryFile(BzzLoad &binaryFile);
	void ReactionRate(	const double T, const double lnT, const BzzVector &cGas, BzzVector &cSurface, BzzVector &aBulk,
						BzzVector &RGas, BzzVector &RSurface, BzzVector &RBulk);

public:
	
	int nReactants;
	int nProducts;
	int nGasReactants;
	int nSiteReactants;
	int nBulkReactants;
	int nGasProducts;
	int nSiteProducts;
	int nBulkProducts;
	double sumNuGas;
	double sumNuSite;
	double sumNuBulk;

	double A;
	double lnA;
	double Beta;
	double E;
	double E_over_R;

	BzzVectorInt indexDirect;
	BzzVector    nuDirect;
	BzzVector    lambdaDirect;

	BzzVectorInt indexInverse;
	BzzVector    nuInverse;
	BzzVector    lambdaInverse;

	int iReversible;
	int tag;

	string reaction_string;

private:

	string name_object;
	void ErrorMessage(const string message);
	void WarningMessage(const string message);
};

#endif	// OPENSMOKE_SURFACEREACTION_H


