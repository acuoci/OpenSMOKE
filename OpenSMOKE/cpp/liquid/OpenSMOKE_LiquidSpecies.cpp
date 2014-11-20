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

#include <sstream>
#include "liquid/OpenSMOKE_LiquidSpecies.h"
#include "liquid/OpenSMOKE_LiquidProperties_Database.h"
#include "liquid/OpenSMOKE_LiquidVaporizationHeat_Dictionary.h"
#include "liquid/OpenSMOKE_LiquidVaporPressure_Dictionary.h"
#include "liquid/OpenSMOKE_LiquidSpecificHeat_Dictionary.h"
#include "liquid/OpenSMOKE_LiquidCriticalConstants_Dictionary.h"
#include "liquid/OpenSMOKE_LiquidThermalConductivity_Dictionary.h"

void OpenSMOKE_LiquidSpecies::ErrorMessage(const std::string message)
{
    cout << endl;
    cout << "Class:  OpenSMOKE_LiquidSpecies"	<< endl;
    cout << "Object: " << name_object			<< endl;
    cout << "Error:  " << message				<< endl;
    cout << "Press a key to continue... "		<< endl;
    getchar();
    exit(-1);
}

void OpenSMOKE_LiquidSpecies::WarningMessage(const std::string message)
{
    cout << endl;
    cout << "Class:  OpenSMOKE_LiquidSpecies"	<< endl;
    cout << "Object: "		<< name_object		<< endl;
    cout << "Warning:  "	<< message			<< endl;
    cout << "Press a key to continue... "		<< endl;
    getchar();
}

OpenSMOKE_LiquidSpecies::OpenSMOKE_LiquidSpecies()
{
	name_object = "[not assigned]";
}

void OpenSMOKE_LiquidSpecies::SetName(const std::string _name)
{
	name = _name;
}

void OpenSMOKE_LiquidSpecies::SetProperties(OpenSMOKE_LiquidProperties_Database &database)
{
	int index = database.RecognizeSpecies(name);

	// Extended name
	name_extended = database.name_extended[index];

	// Constants
	CAS = database.CAS[index];
	MW = database.MW[index];
	Tc = database.Tc[index];
	Pc = database.Pc[index];
	Vc = database.Vc[index];
	Zc = database.Zc[index];
	omega = database.omega[index];

	// Density
	Rho_equation = database.Rho_equation[index];
	Rho_C1 = database.Rho_C1[index];
	Rho_C2 = database.Rho_C2[index];
	Rho_C3 = database.Rho_C3[index];
	Rho_C4 = database.Rho_C4[index];
	Rho_Tmin = database.Rho_Tmin[index];
	Rho_Tmax = database.Rho_Tmax[index];

	// Vapor pressure
	Pv_equation = database.Pv_equation[index];
	Pv_C1 = database.Pv_C1[index];
	Pv_C2 = database.Pv_C2[index];
	Pv_C3 = database.Pv_C3[index];
	Pv_C4 = database.Pv_C4[index];
	Pv_C5 = database.Pv_C5[index];
	Pv_Tmin = database.Pv_Tmin[index];
	Pv_Tmax = database.Pv_Tmax[index];

	// Vaporization heat
	Hv_equation = database.Hv_equation[index];
	Hv_C1 = database.Hv_C1[index];
	Hv_C2 = database.Hv_C2[index];
	Hv_C3 = database.Hv_C3[index];
	Hv_C4 = database.Hv_C4[index];
	Hv_Tmin = database.Hv_Tmin[index];
	Hv_Tmax = database.Hv_Tmax[index];

	// Specific heat
	Cp_equation = database.Cp_equation[index];
	Cp_C1 = database.Cp_C1[index];
	Cp_C2 = database.Cp_C2[index];
	Cp_C3 = database.Cp_C3[index];
	Cp_C4 = database.Cp_C4[index];
	Cp_C5 = database.Cp_C5[index];
	Cp_Tmin = database.Cp_Tmin[index];
	Cp_Tmax = database.Cp_Tmax[index];

	// Thermal conductivity
	Lambda_equation = database.Lambda_equation[index];
	Lambda_C1 = database.Lambda_C1[index];
	Lambda_C2 = database.Lambda_C2[index];
	Lambda_C3 = database.Lambda_C3[index];
	Lambda_C4 = database.Lambda_C4[index];
	Lambda_C5 = database.Lambda_C5[index];
	Lambda_Tmin = database.Lambda_Tmin[index];
	Lambda_Tmax = database.Lambda_Tmax[index];
	if (Lambda_equation == liquid_thermalconductivity_equation::EQ1)
		Lambda_A = Lambda_C1*pow(Lambda_C5, Lambda_C2)/pow(MW, Lambda_C3)/pow(Tc, Lambda_C4);
	if (Lambda_equation == liquid_thermalconductivity_equation::EQ2)
		Lambda_A = Lambda_C3*Lambda_C1*pow(MW, Lambda_C2);
}

void OpenSMOKE_LiquidSpecies::Summary()
{
	double TebN = TNormalBoiling(101325.);

	cout << "-------------------------------------------------" << endl;
	cout << "Liquid properties" << endl;
	cout << "-------------------------------------------------" << endl;
	cout << "Name:        " << name_extended	<< endl;
	cout << "Formula:     " << name				<< endl;
	cout << "CAS:         " << CAS				<< endl;
	cout << "MW[kg/kmol]: " << MW				<< endl;
	cout << "Tc[K]:       " << Tc				<< endl;
	cout << "Pc[bar]:     " << Pc/1.e5			<< endl;
	cout << "Vc[m3/kmol]: " << Vc				<< endl;
	cout << "Zc[-]:       " << Zc				<< endl;
	cout << "omega[-]:    " << omega			<< endl;
	cout << endl;
	cout << "rho@298[kg/m3]:  " << rho(298)			<< endl;
	cout << "Pv@298[Pa]:      " << Pv(298)			<< endl;
	cout << "Hv@298[J/kg]:    " << Hv(298)			<< endl;
	cout << "Cp@298[J/kg/K]:  " << Cp(298)			<< endl;
	cout << "k@298[W/m/K]:    " << lambda(298)		<< endl;
	cout << "TN@1atm[K]:      " << TebN				<< endl;
	cout << endl;
	cout << "rho@TNeb[kg/m3]:  " << rho(TebN)		<< endl;
	cout << "Pv@TNeb[Pa]:      " << Pv(TebN)		<< endl;
	cout << "Hv@TNeb[J/kg]:    " << Hv(TebN)		<< endl;
	cout << "Cp@TNeb[J/kg/K]:  " << Cp(TebN)		<< endl;
	cout << "k@TNeb[W/m/K]:    " << lambda(TebN)	<< endl;
	cout << endl;
}

double OpenSMOKE_LiquidSpecies::rho(const double T)		// [kg/m3]
{
	if (Rho_equation == liquid_density_equation::EQ0)
		return Rho_C1/pow(Rho_C2, (1.+pow(1.-T/Rho_C3,Rho_C4)))*MW;

	else if (Rho_equation == liquid_density_equation::EQ1)
		return 1./(Constants::R_J_kmol*Tc/Pc * pow(Zc, 1.+pow(1.-T/Tc, 2./7.))/MW);

	else if (Rho_equation == liquid_density_equation::EQ7)
	{
		if (T <= 333.15)					return Rho_C1/pow(Rho_C2, (1.+pow(1.-T/Rho_C3,Rho_C4)))*MW;
		else if (T > 333.15 && T<=403.15)	return 4.96838/pow(2.7788e-1, (1.+pow(1.-T/647.13,1.874e-1)))*MW;
		else if (T>403.15)					return 4.358/pow(2.487e-1, (1.+pow(1.-T/647.13,2.534e-1)))*MW;
	}

	return 0;
}

double OpenSMOKE_LiquidSpecies::Pv(const double T)		// [Pa]
{
	if (Pv_equation == liquid_vaporpressure_equation::EQ1)
		return exp(Pv_C1+Pv_C2/T+Pv_C3*log(T)+Pv_C4*pow(T, Pv_C5));

	return 0;
}

double OpenSMOKE_LiquidSpecies::Hv(const double T)		// [J/kg]
{
	if (Hv_equation == liquid_vaporizationheat_equation::EQ1)
	{
		double Tr = T/Tc;
		return Hv_C1*pow(1.-Tr, Hv_C2+Hv_C3*Tr+Hv_C4*Tr*Tr)/MW*1.e7;	
	}

	return 0;
}

double OpenSMOKE_LiquidSpecies::Cp(const double T)		// [J/kg/K]
{
	if (Cp_equation == liquid_specificheat_equation::EQ1)
		return (Cp_C1+T*(Cp_C2+T*(Cp_C3+T*(Cp_C4+Cp_C5*T))))/MW;
	else if (Cp_equation == liquid_specificheat_equation::EQ2)
		ErrorMessage("Specific heat equation 2 not yet implemented...");

	return 0;
}

double OpenSMOKE_LiquidSpecies::lambda(const double T)	// [W/m/K]
{
	if (Lambda_equation == liquid_thermalconductivity_equation::EQ1)
		return Lambda_A*pow(1.-T/Tc, 0.38)/pow(T/Tc,1./6.);
	else if (Lambda_equation == liquid_thermalconductivity_equation::EQ2)
		return Lambda_A*((3.+20.*pow(1.-T/Tc, 0.6666667))/(3.+20.*pow(1-293.15/Tc, 0.6666667)));
	else 
		ErrorMessage("Thermal conductivity equation 0 not yet implemented...");

	return 0;
}

double OpenSMOKE_LiquidSpecies::TNormalBoiling(const double P_Pa)		// [K]
{
	double eps = 1e-2;	// [K]
	int nMax = 18;

	double TA = Pv_Tmin+0.01; 
	double TB = Pv_Tmax-0.01;
	double fA = Pv(TA)-P_Pa;
	double fB = Pv(TB)-P_Pa;

	if (fA*fB>0.)
		ErrorMessage("Impossible to find boiling temperature...");

	for(int i=1;i<=nMax;i++)
	{
		double TC = (TA+TB)*0.50;
		double fC = Pv(TC)-P_Pa;
	
		if (fA*fC<=0.)
		{
			TB = TC;
			fB = fC;
		}
		else
		{
			TA = TC;
			fA = fC;
		}

		if ( (TB-TA)<=eps)
			return TC;
	}

	ErrorMessage("Boiling temperature: Maximum number of iteration...");

	return 0;
}