/***************************************************************************
 *   Copyright (C) 2011 by Alberto Cuoci								   *
 *   alberto.cuoci@polimi.it                                               *
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

#ifndef OpenSMOKE_KPP_SingleReactor_KineticsManager_H
#define OpenSMOKE_KPP_SingleReactor_KineticsManager_H

#include "BzzMath.hpp"
#include "engine/OpenSMOKE_ReactingGas.h"
#include "OpenSMOKE_KPP_Definitions.h"

class OpenSMOKE_ReactingGas;
class OpenSMOKE_SymbolicKinetics;

class OpenSMOKE_KPP_SingleReactor_KineticsManager
{
public:

	OpenSMOKE_KPP_SingleReactor_KineticsManager() {} ;
	void SetMinMax(const double Tmin, const double Tmax) { Tmin_ = Tmin; Tmax_ = Tmax;}
	void Setup(OpenSMOKE_ReactingGas* mixture);
	void Setup(const KPP_Correction correction, OpenSMOKE_ReactingGas* mixture, const double temperature, const double pressure, const double variance);

	void UpdateProperties(BzzVector& omega, const double temperature, const double pressure, BzzVector& R);
	void GetFormationRatesDerivatives(BzzVector& omega, const double temperature, const double pressure, BzzMatrix &JJ);
	double GetMWFromMassFractions(BzzVector &y);

private:

	void RecoverKineticParameters(const double temperature, const double pressure);
	void GetMoleFractionsFromMassFractions(BzzVector &x, BzzVector &y);
	void GetMWAndMoleFractionsFromMassFractions(double &MWmix, BzzVector &x, BzzVector &y);
	
private:

	// Symbolic Kinetics
   	OpenSMOKE_SymbolicKinetics* symbolicKinetics;
	OpenSMOKE_ReactingGas* mix_;
	
	// Auxiliary variables
	bool iSavedKineticConstants_;
	int numberOfSpecies;
	BzzVector x_;
	BzzVector c_;
	BzzVector R_;
	double cTot_;
	BzzVector M;
	BzzVector uM;
	
	double Tmin_;
	double Tmax_;

private:

	void ErrorMessage(const string message_);
	void WarningMessage(const string message_);
};


#endif	// OpenSMOKE_KPP_SingleReactor_KineticsManager_H
