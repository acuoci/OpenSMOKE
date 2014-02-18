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

#ifndef OPENSMOKE_SOLIDMIXTURE_H
#define OPENSMOKE_SOLIDMIXTURE_H

class OpenSMOKE_ReactingGas;

class OpenSMOKE_SolidMixture
{
public:

	OpenSMOKE_SolidMixture();				// Default constructor

    int NC_solid;

    string *gas_names;          // gas names
    string *biomass_names;      // biomass names
    string *solid_names;        // solid names


private:

	string nameObject;

    void compact( BzzMatrixDouble &nu_solid, BzzMatrixDouble &nu_gas,
                  BzzMatrixDouble &eta_solid, BzzMatrixDouble &eta_gas);

    
    int NC_biomass;
    int NC_gas;


    BzzVectorDouble h_solid_298;    // solid specific enthalpy (298K)

    BzzVectorDouble Cp_solid;       // solid specific heat
    BzzVectorDouble lambda_solid;   // solid thermal conductivity
    BzzVectorDouble rho_solid;      // solid densities
    BzzVectorDouble epsilon_solid;  // solid void coefficient

    BzzVectorDouble MW_gas;     // gas molecular weights
    BzzVectorDouble uMW_gas;    // gas molecular weights

    BzzVectorDouble A;          // pre-exponential factor
    BzzVectorDouble Eatt;       // activation energy
    BzzVectorDouble Beta;       // reaction heat

    BzzVectorInt    nEtaGas;
    BzzVectorInt    nEtaSolid;
    BzzVectorInt    compact_gas_eta_index;
    BzzVectorInt    compact_solid_eta_index;
    BzzVectorDouble compact_solid_eta;
    BzzVectorDouble compact_gas_eta;

    BzzVectorInt    nNuGas;
    BzzVectorInt    nNuSolid;
    BzzVectorInt    compact_solid_nu_index;
    BzzVectorInt    compact_gas_nu_index;
    BzzVectorDouble compact_solid_nu;
    BzzVectorDouble compact_gas_nu;
	

	void	GiveMeSpecificHeat(const double T);
    void	GiveMeSpecificEnthalpy(const double T);
	void	GiveMekR(const double T);

	void    ErrorMessage(string message);
    void    WarningMessage(string message);

public:

    BzzVectorDouble uMW_solid;      // solid molecular weights
	   BzzVectorDouble MW_solid;       // solid molecular weights

    void    setup(  string path,
					const nModules,			const string *names_modules,
                    OpenSMOKE_ReactingGas	*gas);

    double  GetDensityApparent(BzzVectorDouble &omega);
	double	GetEpsilon(BzzVectorDouble &omega, const double rho_apparent);
    double  GetSpecificHeat(BzzVectorDouble &omega, const double T);
    double  GetSpecificEnthalpy(BzzVectorDouble &omega, const double T);
    double  GetThermalConductivity(BzzVectorDouble &omega, const double rho_apparent);

    void	GetReactionRates(BzzVectorDouble &Csolid, BzzVectorDouble &xgas, const double T);
    void	GetFormationRates(BzzVectorDouble &Rsolid, BzzVectorDouble &Rgas);
	void	GetReactionEnthalpies(const double T, BzzVectorDouble &h_gas);

    int     recognize_species(const string name);
    int     recognize_species_from_biomass(const string name);

	int NR;
    BzzVectorDouble r;          // kinetic constant
	BzzVectorDouble kR;			// kinetic constant
	BzzVectorDouble DH;			// kinetic constant
	BzzVectorDouble h_solid;        // solid specific enthalpy

private:
	static const double Cstar;
	static const double ALFA;
	static const double H;
	static const double K;
	static const double delta;
	double GiveMeGamma(const double lambda, const double C);

};

#endif // OPENSMOKE_SOLIDMIXTURE_H
