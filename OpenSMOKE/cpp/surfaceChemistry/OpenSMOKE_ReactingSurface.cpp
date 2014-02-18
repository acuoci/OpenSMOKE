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

#include "basic/OpenSMOKE_Constants.h"
#include "surfaceChemistry/OpenSMOKE_SurfaceMaterial.h"
#include "surfaceChemistry/OpenSMOKE_SurfaceSite.h"
#include "surfaceChemistry/OpenSMOKE_SurfaceBulk.h"
#include "surfaceChemistry/OpenSMOKE_ReactingSurface.h"
#include "surfaceChemistry/OpenSMOKE_SurfaceKinetics.h"
#include "surfaceChemistry/OpenSMOKE_SurfaceThermodynamics.h"


void OpenSMOKE_ReactingSurface::ErrorMessage(const string message)
{
    cout << endl;
    cout << "Class:  OpenSMOKE_ReactingSurface"	<< endl;
    cout << "Object: " << name_object		<< endl;
    cout << "Error:  " << message           << endl;
    cout << "Press enter to continue... "   << endl;
    getchar();
    exit(-1);
}

void OpenSMOKE_ReactingSurface::WarningMessage(const string message)
{
    cout << endl;
    cout << "Class:	  OpenSMOKE_ReactingSurface"	<< endl;
    cout << "Object:  " << name_object		<< endl;
    cout << "Warning: " << message          << endl;
	cout << endl;
}

void OpenSMOKE_ReactingSurface::ReadFromBinaryFile(const string fileName)
{
	char dummy[Constants::NAME_SIZE];

	BzzLoad binaryFile;
	binaryFile('*', fileName);

	cout << " * Reading surface kinetic scheme..." << endl;

	binaryFile.fileLoad.read((char*) dummy, sizeof(dummy));
	string version = dummy;
	if		(version == "Surface110514")	iVersion = Surface110514;
	else ErrorMessage("This version is not supported: " + version);

	binaryFile >> numberSiteSpecies_;
	binaryFile >> numberBulkSpecies_;

	cout << " * Reading site species..." << endl;
	namesSiteSpecies_.push_back("site_species");
	for(int k=1;k<=numberSiteSpecies_;k++)
	{
		binaryFile.fileLoad.read((char*) dummy, sizeof(dummy));
		namesSiteSpecies_.push_back(dummy);
	}

	cout << " * Reading bulk species..." << endl;
	namesBulkSpecies_.push_back("bulk_species");
	for(int k=1;k<=numberBulkSpecies_;k++)
	{
		binaryFile.fileLoad.read((char*) dummy, sizeof(dummy));
		namesBulkSpecies_.push_back(dummy);
	}

	cout << " * Reading materials..." << endl;
	binaryFile >> numberMaterials_;
	material_ = new OpenSMOKE_SurfaceMaterial[numberMaterials_+1];

	for(int k=1;k<=numberMaterials_;k++)
	{
		binaryFile.fileLoad.read((char*) dummy, sizeof(dummy));
		material_[k].Setup(this, dummy);
		material_[k].ReadFromBinaryFile(binaryFile);
	}

	// Thermodynamics
	cout << " * Reading thermodynamics..." << endl;
	thermodynamics_ = new OpenSMOKE_SurfaceThermodynamics;
	thermodynamics_->ReadFromBinaryFile(binaryFile, this);

	// Read elements
	cout << " * Reading atomic elements..." << endl;
	binaryFile >> site_elements_;
	binaryFile >> site_m_elements_;
	if (numberBulkSpecies_ > 0)
	{
		binaryFile >> bulk_elements_;
		binaryFile >> bulk_m_elements_;
	}

	// Kinetics
	cout << " * Reading reactions..." << endl;
	for(int k=1;k<=numberMaterials_;k++)
		material_[k].ReadKineticsFromBinaryFile(binaryFile);
	
	binaryFile.End();

	SummaryOnVideo();
}

void OpenSMOKE_ReactingSurface::SummaryOnVideo()
{
	if (numberMaterials_ > 1)
		ErrorMessage("Only 1 Material is currently available");

	for (int k=1;k<=numberMaterials_;k++)
	{
		if (material_[k].NumberOfBulks() > 0)
			ErrorMessage("Only 0 Bulks are currently available");

		if (material_[k].NumberOfSites() > 1)
			ErrorMessage("Only 1 Site is currently available");
	}
}
