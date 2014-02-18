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

#ifndef OPENSMOKE_SURFACETHERMODYNAMICS_H
#define OPENSMOKE_SURFACETHERMODYNAMICS_H

#include "BzzMath.hpp"

class OpenSMOKE_ReactingSurface;

class OpenSMOKE_SurfaceThermodynamics
{
public:

	void ReadFromBinaryFile(BzzLoad &binaryFile, OpenSMOKE_ReactingSurface *surface);

private:

	void Allocate();

private:

	OpenSMOKE_ReactingSurface *ptSurface;

	BzzVector *site_aDH; 
	BzzVector *site_bDH; 
	BzzVector *site_aDS; 
	BzzVector *site_bDS; 
	BzzVector site_T1;
	BzzVector site_T2;
	BzzVector site_T3;
	BzzVector site_M;
	BzzVector *site_CpHT;
	BzzVector *site_CpLT;


	BzzVector *bulk_aDH; 
	BzzVector *bulk_bDH; 
	BzzVector *bulk_aDS; 
	BzzVector *bulk_bDS; 
	BzzVector bulk_T1;
	BzzVector bulk_T2;
	BzzVector bulk_T3;
	BzzVector bulk_M;
	BzzVector *bulk_CpHT;
	BzzVector *bulk_CpLT;

private:

	string name_object;
	void ErrorMessage(const string message);
	void WarningMessage(const string message);
};	

#endif	// OPENSMOKE_SURFACETHERMODYNAMICS_H


