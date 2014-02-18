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

#ifndef OPENSMOKE_SOLIDDATABASE_H
#define OPENSMOKE_SOLIDBIOMASSDATABASE_H

#include "basic/OpenSMOKE_Utilities.h"

class OpenSMOKE_SolidDatabaseModule;

class aliasClass
{
friend class OpenSMOKE_SolidDatabaseModule;
friend class OpenSMOKE_SolidDatabase;

public:

	aliasClass();

protected:

    int NC;
    string *names;
    BzzVectorDouble x;
    BzzVectorInt    index;
    void setup(const int _NC);
    void check(const string _name);

private:

	string nameObject;
    void ErrorMessage(const string message);
    void WarningMessage(const string message);
};

class OpenSMOKE_SolidDatabase
{
public:

        OpenSMOKE_SolidDatabase();

        void setup(	const string path, const int _NC_gas, string *_names_gas, BzzVectorDouble &_MW_gas, 
					const nModules, const string *names_modules);

		void transfer_info( int &_NC_biomass, string *_biomass_names,
							BzzVectorDouble &_A, BzzVectorDouble &_Beta, BzzVectorDouble &_Eatt,
							BzzMatrixDouble &_nu_solid, BzzMatrixDouble &_nu_gas,
							BzzMatrixDouble &_eta_solid, BzzMatrixDouble &_eta_gas);

        void extract_properties( BzzVectorDouble &_rho_solid, BzzVectorDouble &_Cp_solid,
                                 BzzVectorDouble &_h_solid_298,
                                 BzzVectorDouble &_epsilon_solid, BzzVectorDouble &_lambda_solid,
                                 BzzVectorDouble &_MW_solid, BzzVectorDouble &_uMW_solid);

        int NC_gas;
        int NC_solid;
        string *names_solid_database;


 private:

    int NC_alias;
	int NC_biomass;
	string *names_biomass_database;

    BzzVectorDouble *A;                 // pre-exponential factor
    BzzVectorDouble *Beta;              // temperature exponent
    BzzVectorDouble *Eatt;              // activation energy
    BzzVectorInt     nReactionsBiomass; //

    void    ErrorMessage(string message);
    void    WarningMessage(string message);

    void print_summary();

public:
	int     recognize_gas_from_database(string name);
	int     recognize_alias_from_database(string name);
    int     recognize_solid_from_database(string name);

    aliasClass *alias;
    string *names_gas_database;
    string *names_alias_database;

	
    BzzVectorDouble MW_solid;
    BzzVectorDouble MW_gas;

    void extract_info(  const int _NC_biomass, string *_biomass_names,
						BzzVectorDouble &_A, BzzVectorDouble &_Beta, BzzVectorDouble &_Eatt,
						BzzMatrixDouble &_nu_solid, BzzMatrixDouble &_nu_gas,
						BzzMatrixDouble &_eta_solid, BzzMatrixDouble &_eta_gas);


private:
	int     recognize_biomass_from_database(const string name);
    
	void    load_solid_database(const string path);
    void    load_alias_database(const string path);
    
	void	Merge(const int nModules, OpenSMOKE_SolidDatabaseModule *modules);
	void	MergeNames(const int nModules, OpenSMOKE_SolidDatabaseModule *modules);

    BzzMatrixDouble *nu_solid_forward;
    BzzMatrixDouble *nu_solid_backward;
    BzzMatrixDouble *nu_gas_forward;
    BzzMatrixDouble *nu_gas_backward;
    BzzMatrixDouble *eta_solid_forward;
    BzzMatrixDouble *eta_gas_forward;

    BzzVectorDouble rho_solid;
    BzzVectorDouble epsilon_solid;

    BzzVectorDouble lambda_solid;
    BzzVectorDouble h_solid_298;
    BzzVectorDouble Cp_solid;

    string nameObject;
};

#endif // OPENSMOKE_SOLIDDATABASE_H
