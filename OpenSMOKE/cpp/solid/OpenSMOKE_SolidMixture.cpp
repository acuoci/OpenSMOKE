/***************************************************************************
 *   Copyright (C) 2008 by Alberto Cuoci								   *
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

#include "engine/OpenSMOKE_ReactingGas.h"
#include "solid/OpenSMOKE_SolidDatabase.h"
#include "solid/OpenSMOKE_SolidMixture.h"

#include <iostream>
using namespace std;

const double OpenSMOKE_SolidMixture::Cstar	= 1.e-8;
const double OpenSMOKE_SolidMixture::ALFA	= 1.e-5;
const double OpenSMOKE_SolidMixture::H		= 1.50*log(ALFA/(1.-ALFA));
const double OpenSMOKE_SolidMixture::K		= 2.00*log((1.-ALFA)/ALFA) / Cstar;
const double OpenSMOKE_SolidMixture::delta	= 1.e9;

void OpenSMOKE_SolidMixture::ErrorMessage(const string message)
{
    cout << endl;
    cout << "FATAL ERROR" << endl;
    cout << "Class:  OpenSMOKE_SolidMixture" << endl;
    cout << "Object: " << nameObject << endl;
    cout << "Error:  " << message << endl;
    cout << endl;
    cout << "Press a key to continue... " << endl;
    getchar();
    exit(-1);
}

void OpenSMOKE_SolidMixture::WarningMessage(const string message)
{
    cout << endl;
    cout << "WARNING ERROR" << endl;
    cout << "Class:  OpenSMOKE_SolidMixture" << endl;
    cout << "Object:  " << nameObject << endl;
    cout << "Warning: " << message << endl;
    cout << endl;
}

OpenSMOKE_SolidMixture::OpenSMOKE_SolidMixture()
{
    nameObject = "[name not assigned]";
}

double OpenSMOKE_SolidMixture::GetDensityApparent(BzzVectorDouble &omega)
{
	double rho = 0.;
	for(int i=1; i<=NC_solid; i++)
		rho +=  omega[i] / (rho_solid[i]*(1. - epsilon_solid[i]));
	return 1. / rho;
}

double OpenSMOKE_SolidMixture::GetEpsilon(BzzVectorDouble &omega, const double rho_apparent)
{
	double epsilon_total = 0.;
	for(int i=1; i<=NC_solid; i++)
		epsilon_total +=  omega[i] * epsilon_solid[i] / (rho_solid[i]*(1. - epsilon_solid[i]));
	return epsilon_total*rho_apparent;
}

double OpenSMOKE_SolidMixture::GetSpecificHeat(BzzVectorDouble &omega, const double T)
{
	GiveMeSpecificHeat(T);

	return Dot(omega, Cp_solid);
}

double OpenSMOKE_SolidMixture::GetSpecificEnthalpy(BzzVectorDouble &omega, const double T)
{
	GiveMeSpecificEnthalpy(T);

	return Dot(omega, h_solid);
}

double OpenSMOKE_SolidMixture::GetThermalConductivity(BzzVectorDouble &omega, const double rho_apparent)
{
    double lambda = 0.;
	for(int i=1; i<=NC_solid; i++)
		lambda += omega[i] * lambda_solid[i] / rho_solid[i];
		lambda *= rho_apparent;
	return lambda;
}

void OpenSMOKE_SolidMixture::GiveMeSpecificHeat(const double T)
{
    // Cp_solid = 112.  + 4.85*(T-273);		// Wood (Koufopanos, J. Chem. Eng. 1989, 67-75)
    // Cp_solid = 1003. + 2.09*(T-273);		// Char (Koufopanos, J. Chem. Eng. 1989, 67-75)
    // Cp_solid = (-200.+ 80. * sqrt(T));	// J/kg/K  PE
	// Cp_solid = 1255.8;					// J/kg/K	approssimato biomasse
}

void OpenSMOKE_SolidMixture::GiveMeSpecificEnthalpy(const double T)
{
	int i;

	for(i=1; i<=NC_solid; i++)
		//h_solid[i] = h_solid_298[i] - 200.*(T-298.) + 160./3.*(T*sqrt(T) - 5144.277597);  // J/kg  PE
		//h_solid[i] = h_solid_298[i] + 1255.8 * (T - 298.);					// J/kg   approssimato biomasse
		h_solid[i] = h_solid_298[i] + Cp_solid[i] * (T - 298.);					// J/kg
		//h_solid[i] = h_solid_298[i] + 1390. * T + 0.18 * T^2 - 430204;		// J/kg  con Cp Di Blasi
}

void OpenSMOKE_SolidMixture::GiveMekR(double T)
{
    // A    = [kmol, m3, s]
    // Eatt = [J/mol]

    for(int j=1; j<=NR; j++)
    {
        kR[j] = A[j]*exp(-Eatt[j]/(Constants::R_J_mol*T));  // [kmol, m3, s]
        if(Beta[j]!=0)  kR[j] *= pow(T, Beta[j]);           // [kmol, m3, s]
        if(Beta[j]!=0)  kR[j] *= pow(T, Beta[j]);           // [kmol, m3, s]
    }
}

void OpenSMOKE_SolidMixture::GetReactionRates(BzzVectorDouble &Csolid, BzzVectorDouble &xgas, const double T)
{
    int i, j;
    int index;

	GiveMekR(T);

    r = kR;     // reaction rates [kmol/m3/s]

    index = 0;
    for(j=1; j<=NR; j++)
    {
        for(i=1; i<=nEtaSolid[j]; i++)
        {
            index++;
			r[j] *= GiveMeGamma(compact_solid_eta[index], Csolid[compact_solid_eta_index[index]]);
            //r[j] *= pow(Csolid[compact_solid_eta_index[index]], compact_solid_eta[index]);
        }
    }

    index = 0;
    for(j=1; j<=NR; j++)
    {
        for(i=1; i<=nEtaGas[j]; i++)
        {
            index++;
            r[j] *= GiveMeGamma(compact_gas_eta[index], xgas[compact_gas_eta_index[index]]);
            //r[j] *= pow(xgas[compact_gas_eta_index[index]], compact_gas_eta[index]);
        }
    }
}

void OpenSMOKE_SolidMixture::GetFormationRates(BzzVectorDouble &Rsolid, BzzVectorDouble &Rgas)
{
    int i, j;
    int index;

    Rsolid = 0.;    // formation rates [kmol/m3/s]
    Rgas   = 0.;    // formation rates [kmol/m3/s]

    index = 0;
    for(j=1; j<=NR; j++)
    {
        for(i=1; i<=nNuSolid[j]; i++)
        {
            index++;
            Rsolid[compact_solid_nu_index[index]] +=  r[j]*compact_solid_nu[index];
        }
    }

    index = 0;
    for(j=1; j<=NR; j++)
    {
        for(i=1; i<=nNuGas[j]; i++)
        {
            index++;
            Rgas[compact_gas_nu_index[index]] +=  r[j]*compact_gas_nu[index];
        }
    }

    ElementByElementProduct(Rsolid, MW_solid, &Rsolid);     // [kg/m3/s]
    ElementByElementProduct(Rgas, MW_gas, &Rgas);           // [kg/m3/s]
}

void OpenSMOKE_SolidMixture::GetReactionEnthalpies(const double T, BzzVectorDouble &h_gas)
{
    int i, j;
    int index;

	GiveMeSpecificEnthalpy(T);

	DH = 0.;	// [j/kmol]

	// Reaction enthalpies: solid contribution
	index = 0.;
	for(j=1; j<=NR; j++)
    {
        for(i=1; i<=nNuSolid[j]; i++)
        {
            index++;
            DH[j] +=  compact_solid_nu[index]*h_solid[compact_solid_nu_index[index]]*MW_solid[compact_solid_nu_index[index]];
        }
    }

	// Reaction enthalpies: gas contribution
	index = 0.;
	for(j=1; j<=NR; j++)
    {
        for(i=1; i<=nNuGas[j]; i++)
        {
            index++;
		//	cout << "  "  << i << " " << MW_gas[compact_gas_nu_index[index]] << endl;
            DH[j] +=  compact_gas_nu[index] * (h_gas[compact_gas_nu_index[index]]*MW_gas[compact_gas_nu_index[index]]);		// [J/kmol]
        }
    }
	//getchar();
}

void OpenSMOKE_SolidMixture::setup(	string path,
									const nModules,			const string *names_modules,
									OpenSMOKE_ReactingGas	*gas)
{
    int i;

    // 0. Variables declaration
    BzzMatrixDouble nu_solid;   // stoichiometric coefficients
    BzzMatrixDouble nu_gas;     // stoichiometric coefficients
    BzzMatrixDouble eta_solid;  // reaction order
    BzzMatrixDouble eta_gas;    // reaction order

    // 1. Gas Setup
    NC_gas    = gas->NumberOfSpecies();
    gas_names = new string[NC_gas + 1];
    for(i=1;i<=NC_gas;i++)
    {
        gas_names[i] = gas->names[i];
        MW_gas.Append(gas->M[i]);
        uMW_gas.Append(1./gas->M[i]);
    }

    // 2. Biomass database
	{
		OpenSMOKE_SolidDatabase database;
		database.setup(path, NC_gas, gas_names, MW_gas, nModules, names_modules);

		// 3. Transfer info from database to solid mixture
		database.transfer_info(		NC_biomass, biomass_names, A, Beta, Eatt,
									nu_solid, nu_gas, eta_solid, eta_gas);
		database.extract_properties(rho_solid, Cp_solid, h_solid_298, epsilon_solid, lambda_solid, MW_solid, uMW_solid);

	
		// 4. Dimensions and memory allocation
		NR           = A.Size();
		NC_solid     = database.NC_solid;
		solid_names  = new string[NC_solid + 1];
		for(i=1;i<=NC_solid;i++)
			solid_names[i] = database.names_solid_database[i];
	}

    ChangeDimensions(NC_solid, &h_solid); // TODO
    ChangeDimensions(NR, &kR);
    ChangeDimensions(NR, &DH);
    ChangeDimensions(NR, &r);
    ChangeDimensions(NR, &nEtaGas);
    ChangeDimensions(NR, &nEtaSolid);
    ChangeDimensions(NR, &nNuGas);
    ChangeDimensions(NR, &nNuSolid);

    // 5. Compact kinetics
    compact(nu_solid, nu_gas, eta_solid, eta_gas);

}

void OpenSMOKE_SolidMixture::compact( BzzMatrixDouble &nu_solid, BzzMatrixDouble &nu_gas,
                            BzzMatrixDouble &eta_solid, BzzMatrixDouble &eta_gas)
{
    int i,j;

    nEtaGas     = 0;
    nEtaSolid   = 0;
    nNuGas      = 0;
    nNuSolid    = 0;

    for (j=1;j<=NR;j++)
    {
        for(i=1;i<=NC_gas;i++)
            if(eta_gas[j][i] != 0.)
            {
                nEtaGas[j]++;
                compact_gas_eta.Append(eta_gas[j][i]);
                compact_gas_eta_index.Append(i);
            }

        for(i=1;i<=NC_solid;i++)
            if(eta_solid[j][i] != 0.)
            {
                nEtaSolid[j]++;
                compact_solid_eta.Append(eta_solid[j][i]);
                compact_solid_eta_index.Append(i);
            }

        for(i=1;i<=NC_gas;i++)
            if(nu_gas[j][i] != 0.)
            {
                nNuGas[j]++;
                compact_gas_nu.Append(nu_gas[j][i]);
                compact_gas_nu_index.Append(i);
            }

        for(i=1;i<=NC_solid;i++)
            if(nu_solid[j][i] != 0.)
            {
                nNuSolid[j]++;
                compact_solid_nu.Append(nu_solid[j][i]);
                compact_solid_nu_index.Append(i);
            }
    }
}

int OpenSMOKE_SolidMixture::recognize_species(const string name)
{
    for (int i=1;i<=NC_solid;i++)
        if (name == solid_names[i])
            return i;

	ErrorMessage("The following solid species is not included in any solid database: " + name);
    return 0;
}

int OpenSMOKE_SolidMixture::recognize_species_from_biomass(const string name)
{
    int index = 0;

    index = recognize_species(name);

    if (index > 0)
        return index;
    else
	{
		ErrorMessage("The following solid species is not included in any solid database: " + name);
		return 0;
	}
}

double OpenSMOKE_SolidMixture::GiveMeGamma(const double lambda, const double C)
{
	if (lambda>=1.)
		return pow(C, lambda);
	else
	{
		double m = (tanh(K*C + H)+1.)/2.;
		double gamma =	m*pow(C+m/delta,lambda) + 
						(1-m)*pow(Cstar,lambda-1.)*C;
		return gamma;
	}
}