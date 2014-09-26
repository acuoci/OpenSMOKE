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
#include <algorithm>
#include "basic/OpenSMOKE_Constants.h"
#include "surfaceChemistry/OpenSMOKE_SurfaceReaction.h"

void OpenSMOKE_SurfaceReaction::ErrorMessage(const std::string message)
{
    cout << endl;
    cout << "Class:  OpenSMOKE_SurfaceReaction"	<< endl;
    cout << "Object: " << name_object		<< endl;
    cout << "Error:  " << message           << endl;
    cout << "Press enter to continue... "   << endl;
    getchar();
    exit(-1);
}

void OpenSMOKE_SurfaceReaction::WarningMessage(const std::string message)
{
    cout << endl;
    cout << "Class:	  OpenSMOKE_SurfaceReaction"	<< endl;
    cout << "Object:  " << name_object		<< endl;
    cout << "Warning: " << message          << endl;
	cout << endl;
}

void OpenSMOKE_SurfaceReaction::ReadFromBinaryFile(BzzLoad &binaryFile)
{
	int tag;

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

//	cout << nGasReactants <<endl;
//	cout << nSiteReactants <<endl;
//	cout << nBulkReactants <<endl;

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

	if (iReversible == 1)
	{
		binaryFile >> powerSiteDensity >> powerOccupancySites;
	}

	// Tag
	iStickyReaction = false;
	iConvenctionalReaction = false;
	iCoverageDependentReaction=false;

	if (tag == 0)
	{
		iConvenctionalReaction = true;
	}
	
	if (tag == 1)
	{
		iStickyReaction = true;
		binaryFile >> stickyGasSpeciesMW >> stickyPowerOccupancy;
	}

	if (tag == 2)
	{
		cout << "Read 1LH" << endl;
		iLangmuirHinshelwoodReaction = true;
		binaryFile >> langmuirHinshelwoodReactionSpecies;
		binaryFile >> langmuirHinshelwoodReactionParameters;
		cout << "Read 2LH" << endl;
		binaryFile >> langmuirHinshelwoodReactantSpecies;
		binaryFile >> langmuirHinshelwoodReactantExponent;
		cout << "Read 3LH" << endl;
		binaryFile >> langmuirHinshelwoodDenominatorExponentParameter;
		binaryFile >> langmuirHinshelwoodEquilibriumPressure ;
		cout << "Read LH" << endl;
	}

	if (tag == 3)
	{
		iCoverageDependentReaction = true;
		binaryFile >> coverageDependentReactionSpecies >> coverageDependentReactionParameters;
	}
}

double OpenSMOKE_SurfaceReaction::ReactionRate(const double T, const double lnT, const BzzVector &cGas, BzzVector &cSurface, BzzVector &ZSurface, BzzVector &aBulk,
											 BzzVector &RGas, BzzVector &RSurface, BzzVector &RBulk)
{
	double r = 0.0;
	double Kforward = 0.;
	double Kbackward = 0.;

	if (iConvenctionalReaction == true)
	{
		Kforward = lnA + Beta*lnT - E_over_R/T;
		Kforward = exp(Kforward);
	}

	else if (iStickyReaction == true)
	{

		double gamma = min(1., exp(lnA + Beta*lnT - E_over_R/T));
		double coeff = sqrt(Constants::R_J_kmol*T/2./Constants::pi/stickyGasSpeciesMW);
		Kforward = gamma * stickyPowerOccupancy * coeff;
	}

	else if (iCoverageDependentReaction == true)
	{
		cout << "Coverage " << endl;
		Kforward = lnA + Beta*lnT - E_over_R/T;
		Kforward = exp(Kforward);

		for(int k=1;k<=coverageDependentReactionSpecies.Size();k++)
		{
			int j=(k-1)*3;
			Kforward *=
				pow(10., coverageDependentReactionParameters[j+1]*ZSurface[coverageDependentReactionSpecies[k]]) *
				pow(ZSurface[coverageDependentReactionSpecies[k]], coverageDependentReactionParameters[j+2])*
				exp(-coverageDependentReactionParameters[j+3]*ZSurface[coverageDependentReactionSpecies[k]]/(Constants::R_cal_mol*T));
		}
	}

	else if (iLangmuirHinshelwoodReaction == true)
	{
		//BzzVector K(langmuirHinshelwoodReactionSpecies.Size());
		cout << "Langmuir Hinshelwood Reaction";
		cout << "Size " << langmuirHinshelwoodReactionSpecies.Size() << endl;
		Kforward = lnA + Beta*lnT - E_over_R/T;
		Kforward = exp(Kforward);
		double denominator = 1.;

		
		for(int k=1;k<=langmuirHinshelwoodReactionSpecies.Size();k++)
		{
			int j=(k-1)*4;
			cout << "1 " << langmuirHinshelwoodReactionParameters[j+1] << endl; 
			cout << "2 " << langmuirHinshelwoodReactionParameters[j+2] << endl; 
			cout << "3 " << langmuirHinshelwoodReactionParameters[j+3] << endl; 
			cout << "4 " << langmuirHinshelwoodReactionParameters[j+4] << endl; 

			double K = langmuirHinshelwoodReactionParameters[j+1] * 
				       pow(T, langmuirHinshelwoodReactionParameters[j+2]) * 
				       exp(-langmuirHinshelwoodReactionParameters[j+3]/T);///(Constants::R_cal_mol*T));
			cout << k << " " << langmuirHinshelwoodReactionSpecies[k] << " " << K << endl;
			cout << langmuirHinshelwoodReactionSpecies[k] << " " << cGas[langmuirHinshelwoodReactionSpecies[k]] << endl;
			denominator += K * pow(cGas[langmuirHinshelwoodReactionSpecies[k]], langmuirHinshelwoodReactionParameters[j+4]);
			cout << denominator << endl;
		}

		Kforward /= pow(denominator, langmuirHinshelwoodDenominatorExponentParameter);
	}

	
	// Reversible reactions
	if (iReversible == 1)
	{
		double Kequilibrium =	pow(101325./Constants::R_J_kmol/T, sumNuGas) * 
								powerSiteDensity * powerOccupancySites;
		Kbackward = Kequilibrium * Kforward;

//		cout << "Reversible " << powerSiteDensity << " " << powerOccupancySites << " " << sumNuGas << " " << Kequilibrium << " " << Kbackward << endl;
		ErrorMessage("No equilibrium reactions");
//		getchar();
	}

	double productoryReactants = 1.;
	double productoryProducts = 1.;

	// Reactants
	{
		for (int j=1;j<=nGasReactants;j++)
		{
			if (lambdaDirect[j] == 1.)	productoryReactants *= cGas[indexDirect[j]];
			else						productoryReactants *= pow(cGas[indexDirect[j]], lambdaDirect[j]);
		}
		for (int j=nGasReactants+1;j<=nGasReactants+nSiteReactants;j++)
		{	
			if (lambdaDirect[j] == 1.)	productoryReactants *= cSurface[indexDirect[j]];
			else						productoryReactants *= pow(cSurface[indexDirect[j]], lambdaDirect[j]);
		}
		for (int j=nGasReactants+nSiteReactants+1;j<=nGasReactants+nSiteReactants+nBulkReactants;j++)
		{	if (lambdaDirect[j] == 1.)	productoryReactants *= aBulk[indexDirect[j]];
			else						productoryReactants *= pow(aBulk[indexDirect[j]], lambdaDirect[j]);
		}
	}

	// Products
	if (iReversible == 1)
	{
		for (int j=1;j<=nGasProducts;j++)
		{
			if (lambdaInverse[j] == 1.)	productoryProducts *= cGas[indexInverse[j]];
			else						productoryProducts *= pow(cGas[indexInverse[j]], lambdaInverse[j]);
		}
		for (int j=nGasProducts+1;j<=nGasProducts+nSiteProducts;j++)
		{	
			if (lambdaInverse[j] == 1.)	productoryProducts *= cSurface[indexInverse[j]];
			else						productoryProducts *= pow(cSurface[indexInverse[j]], lambdaInverse[j]);
		}
		for (int j=nGasProducts+nSiteProducts+1;j<=nGasProducts+nSiteProducts+nBulkProducts;j++)
		{	if (lambdaInverse[j] == 1.)	productoryProducts *= aBulk[indexInverse[j]];
			else						productoryProducts *= pow(aBulk[indexInverse[j]], lambdaInverse[j]);
		}

//		cout << productoryProducts << endl;
//		getchar();
	}

	r = Kforward*productoryReactants;			// reaction rate!
	if (iReversible == 1)
		r -= Kbackward*productoryProducts;


	// Consumption rates
	for (int j=1;j<=nGasReactants;j++)
		RGas[indexDirect[j]] -= nuDirect[j]*r;

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

	return r;
}
