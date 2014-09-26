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
#include "basic/OpenSMOKE_Constants.h"
#include "surfaceChemistry/OpenSMOKE_ReactingSurface.h"
#include "surfaceChemistry/OpenSMOKE_SurfaceSite.h"
#include "surfaceChemistry/OpenSMOKE_SurfaceBulk.h"
#include "surfaceChemistry/OpenSMOKE_SurfaceKinetics.h"
#include "surfaceChemistry/OpenSMOKE_SurfaceMaterial.h"



void OpenSMOKE_SurfaceMaterial::ErrorMessage(const std::string message)
{
    cout << endl;
    cout << "Class:  OpenSMOKE_SurfaceSiteMaterial"	<< endl;
    cout << "Object: " << name_object		<< endl;
    cout << "Error:  " << message           << endl;
    cout << "Press enter to continue... "   << endl;
    getchar();
    exit(-1);
}

void OpenSMOKE_SurfaceMaterial::WarningMessage(const std::string message)
{
    cout << endl;
    cout << "Class:	  OpenSMOKE_SurfaceMaterial"	<< endl;
    cout << "Object:  " << name_object		<< endl;
    cout << "Warning: " << message          << endl;
	cout << endl;
}

OpenSMOKE_SurfaceMaterial::OpenSMOKE_SurfaceMaterial()
{
	name_object = "[not assigned]";
}

void OpenSMOKE_SurfaceMaterial::Setup(OpenSMOKE_ReactingSurface *surface, const std::string name_)
{
	ptSurface  = surface;
	name_object = name_;
}

void OpenSMOKE_SurfaceMaterial::ReadFromBinaryFile(BzzLoad &binaryFile)
{
	char dummy[Constants::NAME_SIZE];

	{
		binaryFile >> numberSites_;
		sites_ = new OpenSMOKE_SurfaceSite[numberSites_+1];

		for(int j=1;j<=numberSites_;j++)
		{
			double density;
			binaryFile.fileLoad.read((char*) dummy, sizeof(dummy));
			binaryFile >> density;												// [mol/cm2]
			sites_[j].Setup(this, dummy, density*10.);	// [kmol/m2]
		}
	}

	binaryFile >> site_occupancy_matrix_;

	// Bulks
	{
		binaryFile >> numberBulks_;
		if (numberBulks_ > 0)
		{
			bulks_ = new OpenSMOKE_SurfaceBulk[numberBulks_+1];

			for(int j=1;j<=numberBulks_;j++)
			{
				binaryFile.fileLoad.read((char*) dummy, sizeof(dummy));
				bulks_[j].Setup(this, dummy);
			}

			binaryFile >> bulk_density_matrix_;
		}
	}
}

void OpenSMOKE_SurfaceMaterial::ReadKineticsFromBinaryFile(BzzLoad &binaryFile)
{
	kinetics_ = new OpenSMOKE_SurfaceKinetics;

	kinetics_->Setup(ptSurface, name());
	kinetics_->ReadFromBinaryFile(binaryFile);
}

void OpenSMOKE_SurfaceMaterial::ReactionRates(const double T, const BzzVector &cGas, BzzVector &ZSurface, BzzVector &aBulk,
											   BzzVector &RGas, BzzVector &RSurface, BzzVector &RBulk, double &QReaction)
{
	const double lnT = log(T);
	
	// TODO // Only One site
	int k=1;

	BzzVector cSurface(ptSurface->NumberSiteSpecies());
	for(int j=1;j<=ptSurface->NumberSiteSpecies();j++)
	{
		cSurface[j] = ZSurface[j]*sites_[k].density()/site_occupancy_matrix_[j][k];
	
		// TODO
	//	cout << "Zsurf " << endl;
	//	 cout << ptSurface->NameSiteSpecies(j) << " " << ZSurface[j] << " " << sites_[k].density() << " " << site_occupancy_matrix_[j][k] << endl;
	}

	RGas	 = 0.;
	RSurface = 0.;
	RBulk	 = 0.;

	kinetics_->ReactionRates(T, lnT, cGas, cSurface, ZSurface, aBulk, RGas, RSurface, RBulk, QReaction);

	// TODO
	cout << "React rate" << endl;
	for(int j=1;j<=ptSurface->NumberSiteSpecies();j++)
		cout << ptSurface->NameSiteSpecies(j) << " " << RSurface[j]*47.5/1000. << endl;
	for(int j=1;j<=RGas.Size();j++)
		cout << j << " " << RGas[j]*47.5/1000. << endl;
}

