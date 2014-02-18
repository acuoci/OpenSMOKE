/***************************************************************************
 *   Copyright (C) 2010 by Alberto Cuoci								   *
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

#if !defined(OPENSMOKE_LIQUIDSPECIES_H)
#define OPENSMOKE_LIQUIDSPECIES_H

#include "OpenSMOKE.hpp"

#include "liquid/OpenSMOKE_LiquidDensity_Dictionary.h"
#include "liquid/OpenSMOKE_LiquidVaporizationHeat_Dictionary.h"
#include "liquid/OpenSMOKE_LiquidVaporPressure_Dictionary.h"
#include "liquid/OpenSMOKE_LiquidSpecificHeat_Dictionary.h"
#include "liquid/OpenSMOKE_LiquidCriticalConstants_Dictionary.h"
#include "liquid/OpenSMOKE_LiquidThermalConductivity_Dictionary.h"

class OpenSMOKE_LiquidProperties_Database;

class OpenSMOKE_LiquidSpecies
{

public:

	OpenSMOKE_LiquidSpecies();
	void SetName(const string name);
	void SetProperties(OpenSMOKE_LiquidProperties_Database &database);
	void Summary();

	double rho(const double T);				// [kg/m3]
	double Pv(const double T);				// [Pa]
	double Hv(const double T);				// [J/kg]
	double Cp(const double T);				// [J/kg/K]
	double lambda(const double T);			// [W/m/K]
	double TNormalBoiling(const double P_Pa);	// [K]
	double MW;

private:

	int CAS;

	string name_extended;
	string name;

	
	double Tc;
	double Pc;
	double Vc;
	double Zc;
	double omega;

	liquid_vaporpressure_equation::equation Pv_equation;
	double Pv_C1;
	double Pv_C2;
	double Pv_C3;
	double Pv_C4;
	double Pv_C5;
	double Pv_Tmin;
	double Pv_Tmax;

	liquid_specificheat_equation::equation Cp_equation;
	double Cp_C1;
	double Cp_C2;
	double Cp_C3;
	double Cp_C4;
	double Cp_C5;
	double Cp_Tmin;
	double Cp_Tmax;

	liquid_vaporizationheat_equation::equation Hv_equation;
	double Hv_C1;
	double Hv_C2;
	double Hv_C3;
	double Hv_C4;
	double Hv_Tmin;
	double Hv_Tmax;

	liquid_density_equation::equation Rho_equation;
	double Rho_C1;
	double Rho_C2;
	double Rho_C3;
	double Rho_C4;
	double Rho_Tmin;
	double Rho_Tmax;

	liquid_thermalconductivity_equation::equation Lambda_equation;
	double Lambda_A;
	double Lambda_C1;
	double Lambda_C2;
	double Lambda_C3;
	double Lambda_C4;
	double Lambda_C5;
	double Lambda_Tmin;
	double Lambda_Tmax;

private:

	string name_object;
	void ErrorMessage(const string message);
	void WarningMessage(const string message);
};

#endif // OPENSMOKE_LIQUIDSPECIES_H
