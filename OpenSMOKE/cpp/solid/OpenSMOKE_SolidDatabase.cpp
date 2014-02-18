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

#include <sstream>
#include "solid/OpenSMOKE_SolidDatabase.h"
#include "solid/OpenSMOKE_SolidDatabaseModule.h"

char comments[Constants::COMMENT_SIZE];

void OpenSMOKE_SolidDatabase::ErrorMessage(const string message)
{
    cout << endl;
    cout << "FATAL ERROR" << endl;
    cout << "Class:  OpenSMOKE_SolidDatabase" << endl;
    cout << "Object: " << nameObject << endl;
    cout << "Error:  " << message << endl;
    cout << endl;
    cout << "Press a key to continue... " << endl;
    getchar();
    exit(-1);
}

void OpenSMOKE_SolidDatabase::WarningMessage(const string message)
{
    cout << endl;
    cout << "WARNING ERROR" << endl;
    cout << "Class:  OpenSMOKE_SolidDatabase" << endl;
    cout << "Object:  " << nameObject << endl;
    cout << "Warning: " << message << endl;
    cout << endl;
}

OpenSMOKE_SolidDatabase::OpenSMOKE_SolidDatabase()
{
    nameObject = "[name not assigned]";
}


int OpenSMOKE_SolidDatabase::recognize_solid_from_database(const string name)
{
    for(int i=1;i<=NC_solid;i++)
        if (name == names_solid_database[i])
            return i;

    return 0;
}

int OpenSMOKE_SolidDatabase::recognize_gas_from_database(const string name)
{
    for(int i=1;i<=NC_gas;i++)
        if (name == names_gas_database[i])
            return i;

    return 0;
}

int OpenSMOKE_SolidDatabase::recognize_alias_from_database(const string name)
{
    for(int i=1;i<=NC_alias;i++)
        if (name == names_alias_database[i])
            return i;

    return 0;
}

int OpenSMOKE_SolidDatabase::recognize_biomass_from_database(const string name)
{
    for(int i=1;i<=NC_biomass;i++)
        if (name == names_biomass_database[i])
            return i;

    string  message = "The following biomass was not found in the database: ";
    message += name;
    ErrorMessage(message);

    return 0;
}

void OpenSMOKE_SolidDatabase::setup(const string path, const int _NC_gas, string *_names_gas, BzzVectorDouble &_MW_gas,
							const nModules, const string *names_modules)
{
    int i;
    string name;
    ifstream iFile;


    NC_gas = _NC_gas;
    MW_gas = _MW_gas;

    names_gas_database = new string[NC_gas+1];
    for(i=1;i<=NC_gas;i++)
        names_gas_database[i] = _names_gas[i];

    load_alias_database(path);
    load_solid_database(path);

	// Merge Kinetic modules
	{
		OpenSMOKE_SolidDatabaseModule *modules;
		modules = new OpenSMOKE_SolidDatabaseModule[1+nModules];

		for(i=1;i<=nModules;i++)
			modules[i].load_biomass_module_only_names(path, "KineticModules/" + names_modules[i]);
		MergeNames(nModules, modules);

		for(i=1;i<=nModules;i++)
			modules[i].setup(path, "KineticModules/" + names_modules[i], *this);
		Merge(nModules, modules);
	}

    print_summary();
}

void OpenSMOKE_SolidDatabase::load_solid_database(const string path)
{
    ifstream iFile;
    int iEnd;
    string name;

    string fileName = path + "/solid_database.inp";
    openInputFileAndControl(iFile, fileName.c_str());

    NC_solid = 0;

    iEnd = 0;
    while (iEnd == 0)
    {
        iFile >> name;  iFile.getline(comments, Constants::COMMENT_SIZE);
        if (name == "END")
            iEnd = 1;
        else
        {
            NC_solid++;
            iFile.getline(comments, Constants::COMMENT_SIZE);
            iFile.getline(comments, Constants::COMMENT_SIZE);
            iFile.getline(comments, Constants::COMMENT_SIZE);
            iFile.getline(comments, Constants::COMMENT_SIZE);
            iFile.getline(comments, Constants::COMMENT_SIZE);
			iFile.getline(comments, Constants::COMMENT_SIZE);
        }
    }

    iFile.seekg (0, ios::beg);

    names_solid_database = new string[NC_solid+1];
    ChangeDimensions(NC_solid, &rho_solid);
    ChangeDimensions(NC_solid, &epsilon_solid);
    ChangeDimensions(NC_solid, &MW_solid);
    ChangeDimensions(NC_solid, &lambda_solid);
    ChangeDimensions(NC_solid, &h_solid_298);
    ChangeDimensions(NC_solid, &Cp_solid);

    for(int i=1;i<=NC_solid;i++)
    {
		double rho_apparent;
        iFile >> names_solid_database[i];
        iFile >> rho_apparent;	            iFile.getline(comments, Constants::COMMENT_SIZE);
        iFile >> epsilon_solid[i];          iFile.getline(comments, Constants::COMMENT_SIZE);
        iFile >> MW_solid[i];               iFile.getline(comments, Constants::COMMENT_SIZE);
        iFile >> lambda_solid[i];           iFile.getline(comments, Constants::COMMENT_SIZE);
        iFile >> h_solid_298[i];            iFile.getline(comments, Constants::COMMENT_SIZE);
		iFile >> Cp_solid[i];	            iFile.getline(comments, Constants::COMMENT_SIZE);

		rho_solid[i] = rho_apparent/(1.-epsilon_solid[i]);
    }

    iFile.close();
}

void OpenSMOKE_SolidDatabase::load_alias_database(const string path)
{
    int i, j;
    ifstream iFile;
    int iEnd;
    string name;
    int number;

    string fileName = path + "/alias.inp";
    openInputFileAndControl(iFile, fileName.c_str());

    NC_alias = 0;
    iEnd = 0;

    while (iEnd == 0)
    {
        iFile >> name;
        if (name == "END")
            iEnd = 1;
        else
        {
            iFile >> number;    iFile.getline(comments, Constants::COMMENT_SIZE);
            for(i=1;i<=number;i++)
                iFile.getline(comments, Constants::COMMENT_SIZE);
            NC_alias++;
        }
    }
    iFile.seekg(0, ios::beg);

    names_alias_database = new string[NC_alias + 1];
    alias                = new aliasClass[NC_alias + 1];

    for(i=1;i<=NC_alias;i++)
    {
        iFile >> names_alias_database[i];
        iFile >> number;
        alias[i].setup(number);
        for(j=1;j<=number;j++)
        {
            iFile >> alias[i].names[j];
            iFile >> alias[i].x[j];
            alias[i].index[j] = recognize_gas_from_database(alias[i].names[j]);
        }

        alias[i].check(names_alias_database[i]);
    }

    iFile.close();
}

void OpenSMOKE_SolidDatabase::print_summary()
{
    ofstream oFile;
    openOutputFileAndControl(oFile, "ReactionSummary.out");
    oFile.setf(ios::scientific);

    int NR_total = 0;
    for(int indexBiomass=1;indexBiomass<=NC_biomass;indexBiomass++)
    {
        oFile << "-----------------------------------------------------" << endl;
        oFile << names_biomass_database[indexBiomass]                    << endl;
        oFile << "-----------------------------------------------------" << endl;
        for(int indexReaction=1;indexReaction<=A[indexBiomass].Size();indexReaction++)
        {
            NR_total++;
            oFile << "Reaction: "   << indexReaction << " (" << NR_total << ")" << endl;
            oFile << "A    = "      << A[indexBiomass][indexReaction]    << endl;
            oFile << "Beta = "      << Beta[indexBiomass][indexReaction] << endl;
            oFile << "Eatt = "      << Eatt[indexBiomass][indexReaction] << endl;
            oFile << endl;
        }
        oFile << endl;
    }
    oFile.close();
}

void OpenSMOKE_SolidDatabase::transfer_info( int &_NC_biomass, string *_biomass_names,
                                     BzzVectorDouble &_A, BzzVectorDouble &_Beta, BzzVectorDouble &_Eatt,
                                     BzzMatrixDouble &_nu_solid, BzzMatrixDouble &_nu_gas,
                                     BzzMatrixDouble &_eta_solid, BzzMatrixDouble &_eta_gas)
{
    int j;

	_NC_biomass = NC_biomass;
	_biomass_names = new string[_NC_biomass+1];
	for(j=1;j<=NC_biomass;j++)
		_biomass_names[j] = names_biomass_database[j];

	extract_info(_NC_biomass, _biomass_names, _A, _Beta, _Eatt, _nu_solid, _nu_gas, _eta_solid, _eta_gas);
}

void OpenSMOKE_SolidDatabase::extract_info( const int _NC_biomass, string *_biomass_names,
                                    BzzVectorDouble &_A, BzzVectorDouble &_Beta, BzzVectorDouble &_Eatt,
                                    BzzMatrixDouble &_nu_solid, BzzMatrixDouble &_nu_gas,
                                    BzzMatrixDouble &_eta_solid, BzzMatrixDouble &_eta_gas)
{
    int i, j, k;

	cout << "-----------------------------------"	<< endl;
	cout << "Index" << "\t" << "Name" << "\t\t"		<< endl;
	cout << "-----------------------------------"	<< endl;

    for(i=1;i<=_NC_biomass;i++)
    {
        int iFound = 0;
        for(j=1;j<=NC_biomass;j++)
            if(_biomass_names[i] == names_biomass_database[j])
            {
				cout << j << "\t" << names_biomass_database[j] << endl;
                
				for(k=1;k<=A[j].Size();k++)
                {
                    _A.Append(A[j][k]);
                    _Beta.Append(Beta[j][k]);
                    _Eatt.Append(Eatt[j][k]);

                    BzzVectorDouble auxA = nu_solid_backward[j].GetRow(k)   - nu_solid_forward[j].GetRow(k);
                    _nu_solid.AppendRow(auxA);

                    BzzVectorDouble auxB = nu_gas_backward[j].GetRow(k)     - nu_gas_forward[j].GetRow(k);
                    _nu_gas.AppendRow(auxB);

                    BzzVectorDouble auxC = eta_solid_forward[j].GetRow(k);
                    _eta_solid.AppendRow(auxC);

                    BzzVectorDouble auxD = eta_gas_forward[j].GetRow(k);
                    _eta_gas.AppendRow(auxD);

				}

                iFound = 1;
                break;
            }

        if (iFound == 0)
            ErrorMessage("This species cannot be extracted from database: " + _biomass_names[i]);
    }

	cout << "-----------------------------------"	<< endl;
	cout << endl;
}

void OpenSMOKE_SolidDatabase::extract_properties(   BzzVectorDouble &_rho_solid, BzzVectorDouble &_Cp_solid,
                                            BzzVectorDouble &_h_solid_298,
                                            BzzVectorDouble &_epsilon_solid, BzzVectorDouble &_lambda_solid,
                                            BzzVectorDouble &_MW_solid, BzzVectorDouble &_uMW_solid)
{
    _rho_solid      = rho_solid;
    _Cp_solid       = Cp_solid;
    _epsilon_solid  = epsilon_solid;
    _lambda_solid   = lambda_solid;
    _h_solid_298    = h_solid_298;
    _MW_solid       = MW_solid;

    ChangeDimensions(NC_solid, &_uMW_solid);
    for(int i=1;i<=NC_solid;i++)
        _uMW_solid[i] = 1./MW_solid[i];
}

void OpenSMOKE_SolidDatabase::Merge(const int nModules, OpenSMOKE_SolidDatabaseModule *modules)
{
	int i,k,j;

	for(k=1;k<=nModules;k++)
	{
		for(i=1;i<=modules[k].NC_biomass;i++)
		{	
			int indexBiomassGlobal = recognize_biomass_from_database(modules[k].names_biomass_module[i]);
			for(j=1;j<=modules[k].nReactionsBiomass[i];j++)
			{
				A[indexBiomassGlobal].Append(modules[k].A[i][j]);
				Beta[indexBiomassGlobal].Append(modules[k].Beta[i][j]);
				Eatt[indexBiomassGlobal].Append(modules[k].Eatt[i][j]);

				BzzVectorDouble aux;

				aux = modules[k].nu_solid_forward[i].GetRow(j);
				nu_solid_forward[indexBiomassGlobal].AppendRow(aux);
				aux = modules[k].nu_solid_backward[i].GetRow(j);
				nu_solid_backward[indexBiomassGlobal].AppendRow(aux);
				
				aux = modules[k].nu_gas_forward[i].GetRow(j);
				nu_gas_forward[indexBiomassGlobal].AppendRow(aux);
				aux = modules[k].nu_gas_backward[i].GetRow(j);
				nu_gas_backward[indexBiomassGlobal].AppendRow(aux);
				
				aux = modules[k].eta_solid_forward[i].GetRow(j);
				eta_solid_forward[indexBiomassGlobal].AppendRow(aux);
				aux = modules[k].eta_gas_forward[i].GetRow(j);
				eta_gas_forward[indexBiomassGlobal].AppendRow(aux);
			}
		}
	}

	for(i=1;i<=NC_biomass;i++)
		nReactionsBiomass[i] = A[i].Size();

	// Check on double reactions // TODO
}

void OpenSMOKE_SolidDatabase::MergeNames(const int nModules, OpenSMOKE_SolidDatabaseModule *modules)
{
	int i,j,k;

	string *names_biomass_database_provisional;
	names_biomass_database_provisional = new string[1000];

	NC_biomass = 0;
	for(k=1;k<=nModules;k++)
	{
		for(i=1;i<=modules[k].NC_biomass;i++)
		{
			bool iAdd = true;
			for(j=1;j<=NC_biomass;j++)
				if (modules[k].names_biomass_module[i] == names_biomass_database_provisional[j]) 
				{
					iAdd = false;
					break;
				}
			if (iAdd == true)
			{
				NC_biomass++;
				names_biomass_database_provisional[NC_biomass] = modules[k].names_biomass_module[i];
			}
		}
	}

    names_biomass_database = new string[NC_biomass+1];

	for(k=1;k<=NC_biomass;k++)
		names_biomass_database[k] = names_biomass_database_provisional[k];

    A					= new BzzVectorDouble[NC_biomass+1];
    Beta				= new BzzVectorDouble[NC_biomass+1];
    Eatt				= new BzzVectorDouble[NC_biomass+1];
    nu_solid_forward	= new BzzMatrixDouble[NC_biomass+1];
    nu_solid_backward	= new BzzMatrixDouble[NC_biomass+1];
    nu_gas_forward		= new BzzMatrixDouble[NC_biomass+1];
    nu_gas_backward		= new BzzMatrixDouble[NC_biomass+1];
    eta_solid_forward	= new BzzMatrixDouble[NC_biomass+1];
    eta_gas_forward		= new BzzMatrixDouble[NC_biomass+1];
    ChangeDimensions(NC_biomass, &nReactionsBiomass);
}


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

void aliasClass::setup(const int _NC)
{
    NC = _NC;
    names = new string[NC+1];
    ChangeDimensions(NC, &x);
    ChangeDimensions(NC, &index);
}

void aliasClass::check(const string _name)
{
    int j;

    for(j=1;j<=NC;j++)
        if (index[j] <= 0)
        {
            string message = "The " + _name + " alias contains the species " + names[j] + " which was not found in the Gas Kinetic Mechanism";
            ErrorMessage(message);
        }
}

void aliasClass::ErrorMessage(const string message)
{
    cout << endl;
    cout << "FATAL ERROR" << endl;
    cout << "Class:  Alias" << endl;
    cout << "Object: " << nameObject << endl;
    cout << "Error:  " << message << endl;
    cout << endl;
    cout << "Press a key to continue... " << endl;
    getchar();
    exit(-1);
}

void aliasClass::WarningMessage(const string message)
{
    cout << endl;
    cout << "WARNING ERROR" << endl;
    cout << "Class:  Alias" << endl;
    cout << "Object:  " << nameObject << endl;
    cout << "Warning: " << message << endl;
    cout << endl;
}

aliasClass::aliasClass()
{
    nameObject = "[name not assigned]";
}

