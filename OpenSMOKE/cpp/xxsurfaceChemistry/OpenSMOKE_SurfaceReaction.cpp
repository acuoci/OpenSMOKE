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
#include "surfaceChemistry/OpenSMOKE_SurfaceReaction.h"

void OpenSMOKE_SurfaceReaction::ErrorMessage(const string message)
{
    cout << endl;
    cout << "Class:  OpenSMOKE_SurfaceReaction"	<< endl;
    cout << "Object: " << name_object		<< endl;
    cout << "Error:  " << message           << endl;
    cout << "Press enter to continue... "   << endl;
    getchar();
    exit(-1);
}

void OpenSMOKE_SurfaceReaction::WarningMessage(const string message)
{
    cout << endl;
    cout << "Class:	  OpenSMOKE_SurfaceReaction"	<< endl;
    cout << "Object:  " << name_object		<< endl;
    cout << "Warning: " << message          << endl;
	cout << endl;
}

void OpenSMOKE_SurfaceReaction::ReadFromBinaryFile(BzzLoad &binaryFile)
{
	binaryFile >> iReversible;
	binaryFile >> tag;

	binaryFile >> nReactants;
	binaryFile >> nProducts;
	binaryFile >> nGasReactants;
	binaryFile >> nSiteReactants;
	binaryFile >> nBulkReactants;
	binaryFile >> nGasProducts;
	binaryFile >> nSiteProducts;
	binaryFile >> nBulkProducts;
	binaryFile >> sumNuGas;
	binaryFile >> sumNuSite;
	binaryFile >> sumNuBulk;

	cout << nGasReactants <<endl;
	cout << nSiteReactants <<endl;
	cout << nBulkReactants <<endl;

	ChangeDimensions(nReactants, &indexDirect);
	ChangeDimensions(nReactants, &nuDirect);
	ChangeDimensions(nReactants, &lambdaDirect);
	ChangeDimensions(nProducts,  &indexInverse);
	ChangeDimensions(nProducts,  &nuInverse);
	ChangeDimensions(nProducts,  &lambdaInverse);

	for(int j=1;j<=nReactants;j++)
		binaryFile >> indexDirect[j] >> nuDirect[j] >> lambdaDirect[j];

	for(int j=1;j<=nProducts;j++)
		binaryFile >> indexInverse[j] >> nuInverse[j] >> lambdaInverse[j];

	// Name
	char dummy[Constants::REACTION_NAME_SIZE];
	binaryFile.fileLoad.read((char*) dummy, sizeof(dummy));
	reaction_string = dummy;


	// Kinetic parameters
	binaryFile >> A;
	binaryFile >> Beta;
	binaryFile >> E;

	// Conversions
	lnA = log(A);
	E_over_R = E/Constants::R_cal_mol;
}

void OpenSMOKE_SurfaceReaction::ReactionRate(const double T, const double lnT, const BzzVector &cGas, BzzVector &cSurface, BzzVector &aBulk,
											 BzzVector &RGas, BzzVector &RSurface, BzzVector &RBulk)
{
	double r;		// [kmol/m2/s]

//  TODO
//	cout << A << " " << Beta << " " << E << " " << E_over_R << endl;

	r = lnA + Beta*lnT - E_over_R/T;
	
	r = exp(r);
	
//  TODO
//	cout << " " << r << endl;

	for (int j=1;j<=nGasReactants;j++)
	{
		if (lambdaDirect[j] == 1.)	r *= cGas[indexDirect[j]];
		else						r *= pow(cGas[indexDirect[j]], lambdaDirect[j]);
		//  TODO
		// cout << j << " " << " G " << indexDirect[j] << " " <<  lambdaDirect[j] << " " << cGas[indexDirect[j]] << " " << endl;
	}

	for (int j=nGasReactants+1;j<=nGasReactants+nSiteReactants;j++)
	{	
		if (lambdaDirect[j] == 1.)	r *= cSurface[indexDirect[j]];
		else						r *= pow(cSurface[indexDirect[j]], lambdaDirect[j]);
		//  TODO
		// cout << j << " " << " S " << indexDirect[j] << " " <<  lambdaDirect[j] << " " << cSurface[indexDirect[j]] << endl;
	}
	for (int j=nGasReactants+nSiteReactants+1;j<=nGasReactants+nSiteReactants+nBulkReactants;j++)
	{	if (lambdaDirect[j] == 1.)	r *= aBulk[indexDirect[j]];
		else						r *= pow(aBulk[indexDirect[j]], lambdaDirect[j]);
		//  TODO
		//ErrorMessage("No bulk species are available (yet)");
	}

	// Consumption rates
	for (int j=1;j<=nGasReactants;j++)
	{	
		// TODO
		//cout << j << " " << indexDirect[j] << " " << nuDirect[j] << " " << RGas.Size() << endl;
		RGas[indexDirect[j]] -= nuDirect[j]*r;
	}

	for (int j=nGasReactants+1;j<=nGasReactants+nSiteReactants;j++)
		RSurface[indexDirect[j]] -= nuDirect[j]*r;

	for (int j=nGasReactants+nSiteReactants+1;j<=nGasReactants+nSiteReactants+nBulkReactants;j++)
		RBulk[indexDirect[j]] -= nuDirect[j]*r;

	// Formation rates
	for (int j=1;j<=nGasProducts;j++)
		RGas[indexInverse[j]] += nuInverse[j]*r;

	for (int j=nGasProducts+1;j<=nGasProducts+nSiteProducts;j++)
		RSurface[indexInverse[j]] += nuInverse[j]*r;

	for (int j=nGasProducts+nSiteProducts+1;j<=nGasProducts+nSiteProducts+nBulkProducts;j++)
		RBulk[indexInverse[j]] += nuInverse[j]*r;

// TODO
//	cout << r << endl;
//	cout << r*47.5 << endl;
//	cout << r*47.5/1000. << endl;
}