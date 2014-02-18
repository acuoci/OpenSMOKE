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

#ifndef OPENSMOKE_SOLIDDATABASEMODULE_H
#define OPENSMOKE_SOLIDDATABASEMODULE_H

class OpenSMOKE_SolidDatabase;

class OpenSMOKE_SolidDatabaseModule
{
public:

        OpenSMOKE_SolidDatabaseModule();

        void setup(const string path, const string moduleName, OpenSMOKE_SolidDatabase &database);
		void load_biomass_module_only_names(const string path, const string moduleName);
		int		NC_biomass;
		string *names_biomass_module;

	BzzVectorDouble *A;                 // pre-exponential factor
    BzzVectorDouble *Beta;              // temperature exponent
    BzzVectorDouble *Eatt;              // activation energy
	BzzMatrixDouble *nu_solid_forward;
    BzzMatrixDouble *nu_solid_backward;
    BzzMatrixDouble *nu_gas_forward;
    BzzMatrixDouble *nu_gas_backward;
    BzzMatrixDouble *eta_solid_forward;
    BzzMatrixDouble *eta_gas_forward;

	BzzVectorInt     nReactionsBiomass; // number of reactions for each component


protected:

private:

	OpenSMOKE_SolidDatabase *database;



    void    ErrorMessage(const string message);
    void    WarningMessage(const string message);
    void    read_reaction(ifstream &iFile, const int indexBiomass, const int indexReaction);
    double  read_number(ifstream &iFile, const int Biomass, const int indexReaction, const string symbol);
    void    check_reaction(const int indexBiomass, const int indexReaction);

    void print_summary();

    void    load_biomass_module(string moduleName);
	int     recognize_biomass_from_module(string name);



    string nameObject;
};

#endif // OPENSMOKE_SOLIDDATABASEMODULE_H