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

#if !defined(OPENSMOKE_LIQUIDPROPERTIES_DATABASE_H)
#define OPENSMOKE_LIQUIDPROPERTIES_DATABASE_H

#include "OpenSMOKE.hpp"

#include "liquid/OpenSMOKE_LiquidDensity_Dictionary.h"
#include "liquid/OpenSMOKE_LiquidVaporizationHeat_Dictionary.h"
#include "liquid/OpenSMOKE_LiquidVaporPressure_Dictionary.h"
#include "liquid/OpenSMOKE_LiquidSpecificHeat_Dictionary.h"
#include "liquid/OpenSMOKE_LiquidCriticalConstants_Dictionary.h"
#include "liquid/OpenSMOKE_LiquidThermalConductivity_Dictionary.h"

class OpenSMOKE_LiquidProperties_Database
{
friend class OpenSMOKE_LiquidSpecies;

public:

	OpenSMOKE_LiquidProperties_Database();
	void ReadFromFolder(const string folder_name);
	void SaveOnFile(const string file_name);
	void LoadFromFile(const string file_name);
	int  RecognizeSpecies(const string name);

private:

	OpenSMOKE_LiquidDensity_Dictionary				*dictionary_density;
	OpenSMOKE_LiquidVaporizationHeat_Dictionary		*dictionary_vaporization_heat;
	OpenSMOKE_LiquidVaporPressure_Dictionary		*dictionary_vapor_pressure;
	OpenSMOKE_LiquidSpecificHeat_Dictionary			*dictionary_specific_heat;
	OpenSMOKE_LiquidCriticalConstants_Dictionary	*dictionary_critical_constants;
	OpenSMOKE_LiquidThermalConductivity_Dictionary	*dictionary_thermal_conductivity;

	void CheckConsistency();

private:

	int N;
	vector<string> name_extended;
	vector<string> name_first;
	vector<string> name_second;

	BzzVectorInt CAS;
	BzzVector MW;
	BzzVector Tc;
	BzzVector Pc;
	BzzVector Vc;
	BzzVector Zc;
	BzzVector omega;

	vector<liquid_vaporpressure_equation::equation> Pv_equation;
	BzzVector Pv_C1;
	BzzVector Pv_C2;
	BzzVector Pv_C3;
	BzzVector Pv_C4;
	BzzVector Pv_C5;
	BzzVector Pv_Tmin;
	BzzVector Pv_Tmax;

	vector<liquid_specificheat_equation::equation> Cp_equation;
	BzzVector Cp_C1;
	BzzVector Cp_C2;
	BzzVector Cp_C3;
	BzzVector Cp_C4;
	BzzVector Cp_C5;
	BzzVector Cp_Tmin;
	BzzVector Cp_Tmax;

	vector<liquid_vaporizationheat_equation::equation> Hv_equation;
	BzzVector Hv_C1;
	BzzVector Hv_C2;
	BzzVector Hv_C3;
	BzzVector Hv_C4;
	BzzVector Hv_Tmin;
	BzzVector Hv_Tmax;

	vector<liquid_density_equation::equation> Rho_equation;
	BzzVector Rho_C1;
	BzzVector Rho_C2;
	BzzVector Rho_C3;
	BzzVector Rho_C4;
	BzzVector Rho_Tmin;
	BzzVector Rho_Tmax;

	vector<liquid_thermalconductivity_equation::equation> Lambda_equation;
	BzzVector Lambda_C1;
	BzzVector Lambda_C2;
	BzzVector Lambda_C3;
	BzzVector Lambda_C4;
	BzzVector Lambda_C5;
	BzzVector Lambda_Tmin;
	BzzVector Lambda_Tmax;

private:

	string name_object;
	void ErrorMessage(const string message);
	void WarningMessage(const string message);
};

#endif // OPENSMOKE_LIQUIDPROPERTIESDATABASE_H
