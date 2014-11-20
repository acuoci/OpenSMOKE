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

#include <string>
#include "surfaceChemistry/OpenSMOKE_ReactingSurface.h"
#include "surfaceChemistry/OpenSMOKE_SurfaceThermodynamics.h"

void OpenSMOKE_SurfaceThermodynamics::ErrorMessage(const std::string message)
{
    cout << endl;
    cout << "Class:  OpenSMOKE_SurfaceThermodynamics"	<< endl;
    cout << "Object: " << name_object		<< endl;
    cout << "Error:  " << message           << endl;
    cout << "Press enter to continue... "   << endl;
    getchar();
    exit(-1);
}

void OpenSMOKE_SurfaceThermodynamics::WarningMessage(const std::string message)
{
    cout << endl;
    cout << "Class:	  OpenSMOKE_SurfaceThermodynamics"	<< endl;
    cout << "Object:  " << name_object		<< endl;
    cout << "Warning: " << message          << endl;
	cout << endl;
}

void OpenSMOKE_SurfaceThermodynamics::Allocate()
{
	// Site species
	{
		int n = ptSurface->NumberSiteSpecies();

		site_aDH	=	new BzzVector [n+1];
		site_bDH	=   new BzzVector [n+1];
		site_aDS	=   new BzzVector [n+1];
		site_bDS	=   new BzzVector [n+1];

		site_CpHT	=	new BzzVector [n+1];
		site_CpLT	=	new BzzVector [n+1];

		for(int k=1;k<=n;k++)
		{
			ChangeDimensions(6,&site_aDH[k]);
			ChangeDimensions(6,&site_bDH[k]);
			ChangeDimensions(6,&site_aDS[k]);
			ChangeDimensions(6,&site_bDS[k]);

			ChangeDimensions(5,&site_CpHT[k]);
			ChangeDimensions(5,&site_CpLT[k]);
		}
	}

	// Bulk species
	{
		int n = ptSurface->NumberBulkSpecies();
		if (n>0)
		{
			bulk_aDH	=	new BzzVector [n+1];
			bulk_bDH	=   new BzzVector [n+1];
			bulk_aDS	=   new BzzVector [n+1];
			bulk_bDS	=   new BzzVector [n+1];

			bulk_CpHT	=	new BzzVector [n+1];
			bulk_CpLT	=	new BzzVector [n+1];

			for(int k=1;k<=n;k++)
			{
				ChangeDimensions(6,&bulk_aDH[k]);
				ChangeDimensions(6,&bulk_bDH[k]);
				ChangeDimensions(6,&bulk_aDS[k]);
				ChangeDimensions(6,&bulk_bDS[k]);

				ChangeDimensions(5,&bulk_CpHT[k]);
				ChangeDimensions(5,&bulk_CpLT[k]);
			}
		}
	}
}


void OpenSMOKE_SurfaceThermodynamics::ReadFromBinaryFile(BzzLoad &binaryFile, OpenSMOKE_ReactingSurface *surface)
{
	ptSurface = surface;

	// Memory allocation
	Allocate();

	// Site thermodynamics
	{
		for(int k=1;k<=ptSurface->NumberSiteSpecies();k++)
		{
			binaryFile >> site_CpHT[k];
			binaryFile >> site_CpLT[k];
			binaryFile >> site_aDH[k];
			binaryFile >> site_bDH[k];
			binaryFile >> site_aDS[k];
			binaryFile >> site_bDS[k];
		}
		binaryFile >> site_M;
		binaryFile >> site_T1;
		binaryFile >> site_T2;
		binaryFile >> site_T3;
	}

	// Bulk thermodynamics
	if (ptSurface->NumberBulkSpecies()>0)
	{
		for(int k=1;k<=ptSurface->NumberBulkSpecies();k++)
		{
			binaryFile >> bulk_CpHT[k];
			binaryFile >> bulk_CpLT[k];
			binaryFile >> bulk_aDH[k];
			binaryFile >> bulk_bDH[k];
			binaryFile >> bulk_aDS[k];
			binaryFile >> bulk_bDS[k];
		}
		binaryFile >> bulk_M;
		binaryFile >> bulk_T1;
		binaryFile >> bulk_T2;
		binaryFile >> bulk_T3;
	}
}
