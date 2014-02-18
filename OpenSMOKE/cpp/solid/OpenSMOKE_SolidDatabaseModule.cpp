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

#include <iostream>
using namespace std;

void OpenSMOKE_SolidDatabaseModule::ErrorMessage(const string message)
{
    cout << endl;
    cout << "FATAL ERROR" << endl;
    cout << "Class:  OpenSMOKE_SolidDatabaseModule" << endl;
    cout << "Object: " << nameObject << endl;
    cout << "Error:  " << message << endl;
    cout << endl;
    cout << "Press a key to continue... " << endl;
    getchar();
    exit(-1);
}

void OpenSMOKE_SolidDatabaseModule::WarningMessage(const string message)
{
    cout << endl;
    cout << "WARNING ERROR" << endl;
    cout << "Class:  OpenSMOKE_SolidDatabaseModule" << endl;
    cout << "Object:  " << nameObject << endl;
    cout << "Warning: " << message << endl;
    cout << endl;
}

OpenSMOKE_SolidDatabaseModule::OpenSMOKE_SolidDatabaseModule()
{
    nameObject = "[name not assigned]";
}

int OpenSMOKE_SolidDatabaseModule::recognize_biomass_from_module(string name)
{
    for(int i=1;i<=NC_biomass;i++)
        if (name == names_biomass_module[i])
            return i;

    string  message = "The following biomass was not found in the module: ";
    message += name;
    ErrorMessage(message);

    return 0;
}

void OpenSMOKE_SolidDatabaseModule::read_reaction(ifstream &iFile, const int indexBiomass, const int indexReaction)
{
    string dummy;
    double number;
    int indexSpecies;
    string nameReactants[20];

    int iEnd       = 0;
    int iProducts  = 0;
    int nReactants = 0;

    while (iEnd == 0)
    {
        iFile >> dummy;
		

        if (dummy == "/")
        {
            number = read_number(iFile, indexBiomass, indexReaction, "/");
            A[indexBiomass][indexReaction] = number;

            number = read_number(iFile, indexBiomass, indexReaction, "/");
            Beta[indexBiomass][indexReaction] = number;

            number = read_number(iFile, indexBiomass, indexReaction, "//");
            Eatt[indexBiomass][indexReaction] = number;

            for(int j=1;j<=nReactants;j++)
            {
                number = read_number(iFile, indexBiomass, indexReaction, "/");
                indexSpecies = database->recognize_solid_from_database(nameReactants[j]);
                if (indexSpecies > 0)
                    eta_solid_forward[indexBiomass][indexReaction][indexSpecies] = number;
                else
                {
                    indexSpecies = database->recognize_gas_from_database(nameReactants[j]);
                    if (indexSpecies > 0)
                        eta_gas_forward[indexBiomass][indexReaction][indexSpecies] = number;
                    else
                    {
                        stringstream s;
                        s << indexReaction;
                        string message = "Error in reading reaction order: ";
                        message += "Biomass  " + names_biomass_module[indexBiomass] + " - " +
                                   "Reaction " + s.str() + " - " +
                                   "Species  " + nameReactants[j];
                        ErrorMessage(message);
                    }
                }
            }

            iEnd = 1;
        }
        else if (dummy == "=>")
        {
            iProducts = 1;
        }
        else if (dummy == "+")
        {
        }
        else
        {
			
            number = atof(dummy.c_str());
			
            if (number <=0)
            {
                stringstream sA, sB;
                sA << indexReaction;
                sB << number;
                string message = "Error in reading stoichiometric coefficient - Biomass "
                                    + names_biomass_module[indexBiomass]
                                            + "  -  Reaction: "    + sA.str()
                                            + "  -  Read number: " + sB.str();
                ErrorMessage(message);
            }
            iFile >> dummy;

            indexSpecies = database->recognize_solid_from_database(dummy);
            if (indexSpecies > 0)
            {
                if (iProducts == 0)
                {
                    nReactants++;
                    nameReactants[nReactants] = dummy;
                    nu_solid_forward[indexBiomass][indexReaction][indexSpecies] = number;
                }
                else
                    nu_solid_backward[indexBiomass][indexReaction][indexSpecies] = number;
            }
            else
            {
                indexSpecies = database->recognize_gas_from_database(dummy);
                if (indexSpecies > 0)
                {
                    if (iProducts == 0)
                    {
                        nReactants++;
                        nameReactants[nReactants] = dummy;
                        nu_gas_forward[indexBiomass][indexReaction][indexSpecies] = number;
                    }
                    else
                        nu_gas_backward[indexBiomass][indexReaction][indexSpecies] = number;
                }
                else
                {
                    int indexAlias = database->recognize_alias_from_database(dummy);
                    if (indexAlias > 0)
                    {
                        for(int j=1;j<=database->alias[indexAlias].NC;j++)
                        {
                            if (iProducts == 0)
                                nu_gas_forward[indexBiomass][indexReaction][database->alias[indexAlias].index[j]]
                                                    = number * database->alias[indexAlias].x[j];
                            else
                                nu_gas_backward[indexBiomass][indexReaction][database->alias[indexAlias].index[j]]
                                                    = number * database->alias[indexAlias].x[j];
                        }
                    }
                    else
                    {
                        string message = "The following species was not found: ";
                        message += dummy;
                        ErrorMessage(message);
                    }
                }
            }
        }
    }

    check_reaction(indexBiomass, indexReaction);
}

void OpenSMOKE_SolidDatabaseModule::check_reaction(const int indexBiomass, const int indexReaction)
{
    int i;
    double sum_Forward;
    double sum_Backward;

    sum_Forward =0.;
    for(i=1;i<=database->NC_solid;i++)
        sum_Forward += nu_solid_forward[indexBiomass][indexReaction][i] * database->MW_solid[i];
    for(i=1;i<=database->NC_gas;i++)
        sum_Forward += nu_gas_forward[indexBiomass][indexReaction][i] * database->MW_gas[i];

    sum_Backward =0.;
    for(i=1;i<=database->NC_solid;i++)
        sum_Backward += nu_solid_backward[indexBiomass][indexReaction][i] * database->MW_solid[i];
    for(i=1;i<=database->NC_gas;i++)
        sum_Backward += nu_gas_backward[indexBiomass][indexReaction][i] * database->MW_gas[i];

    double error_relative = fabs(sum_Forward-sum_Backward)/max(sum_Forward, sum_Backward);
    if (error_relative > 1.e-4)
    {
        stringstream sA;
        stringstream sB;
        stringstream sC;

        sA << indexReaction;
        string message = "The following reaction is not balanced: " + sA.str();

        message += "  -  Biomass: " + names_biomass_module[indexBiomass] + "\n";

        sB << sum_Forward;
        message += "  MW_forward  = " + sB.str() + " kg/kmol  -  ";

        sC << sum_Backward;
        message += "  MW_backward = " + sC.str() + " kg/kmol";

        ErrorMessage(message);
    }
}

double OpenSMOKE_SolidDatabaseModule::read_number(ifstream &iFile, const int indexBiomass, const int indexReaction, const string symbol)
{
    string dummy;
    double number;
    stringstream s;

    string message = "Error in reaction number ";
    s << indexReaction;
    message += s.str();
    message += "  -  Biomass: " + names_biomass_module[indexBiomass];


    iFile >> dummy;
    number = atof(dummy.c_str());
    if (number <=-1e16 || number >=1e32)
    {
        message += " - Found " + dummy;
        ErrorMessage(message);
    }

    iFile >> dummy;
    if (dummy != symbol)
    {
        message += " - Expected " + symbol + " Found " + dummy;
        ErrorMessage(message);
    }
    return number;
}


void OpenSMOKE_SolidDatabaseModule::setup(string path, string moduleName, OpenSMOKE_SolidDatabase &_database)
{
	database = &_database;

    int i;
    string name;
    int indexReaction;
    int indexBiomass;
    ifstream iFile;
	char comments[Constants::COMMENT_SIZE];

	// Load biomass module
    string fileName = path + "/" + moduleName;
	load_biomass_module(fileName);

	// Read reactions
	openInputFileAndControl(iFile, fileName.c_str());
    for(i=1;i<=NC_biomass+2;i++)
        iFile.getline(comments, Constants::COMMENT_SIZE);

    for(i=1;i<=NC_biomass;i++)
    {
        iFile >> name;
        indexBiomass = recognize_biomass_from_module(name);
        iFile >> nReactionsBiomass[indexBiomass];

        ChangeDimensions(nReactionsBiomass[indexBiomass], &A[indexBiomass]);
        ChangeDimensions(nReactionsBiomass[indexBiomass], &Beta[indexBiomass]);
        ChangeDimensions(nReactionsBiomass[indexBiomass], &Eatt[indexBiomass]);

        ChangeDimensions(nReactionsBiomass[indexBiomass], database->NC_solid,	&nu_solid_forward[indexBiomass]);
        ChangeDimensions(nReactionsBiomass[indexBiomass], database->NC_solid,	&nu_solid_backward[indexBiomass]);
        ChangeDimensions(nReactionsBiomass[indexBiomass], database->NC_gas,		&nu_gas_forward[indexBiomass]);
        ChangeDimensions(nReactionsBiomass[indexBiomass], database->NC_gas,		&nu_gas_backward[indexBiomass]);

        ChangeDimensions(nReactionsBiomass[indexBiomass], database->NC_solid,	&eta_solid_forward[indexBiomass]);
        ChangeDimensions(nReactionsBiomass[indexBiomass], database->NC_gas,		&eta_gas_forward[indexBiomass]);


        for(indexReaction=1;indexReaction<=nReactionsBiomass[indexBiomass];indexReaction++)
            read_reaction(iFile, indexBiomass, indexReaction);
    }

    print_summary();
}

void OpenSMOKE_SolidDatabaseModule::load_biomass_module(string moduleName)
{
    ifstream iFile;
    int iEnd;
    string name;

	char comments[Constants::COMMENT_SIZE];

    openInputFileAndControl(iFile, moduleName.c_str());

    iFile >> name;
    if (name != "SPECIES")
        ErrorMessage("The SPECIES keyword was not found in the solid kinetic database!");

    NC_biomass = 0;
    iEnd = 0;
    while(iEnd == 0)
    {
        iFile >> name;
        if (name == "END") iEnd = 1;
        else               NC_biomass++;
    }
    iFile.seekg(0, ios::beg);

    names_biomass_module = new string[NC_biomass+1];

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

    iFile.getline(comments, Constants::COMMENT_SIZE);
    for(int i=1;i<=NC_biomass;i++)
        iFile >> names_biomass_module[i];

    iFile.close();
}


void OpenSMOKE_SolidDatabaseModule::print_summary()
{
    ofstream oFile;
    openOutputFileAndControl(oFile, "ReactionSummary.out");
    oFile.setf(ios::scientific);

    int NR_total = 0;
    for(int indexBiomass=1;indexBiomass<=NC_biomass;indexBiomass++)
    {
        oFile << "-----------------------------------------------------" << endl;
        oFile << names_biomass_module[indexBiomass]                    << endl;
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

void OpenSMOKE_SolidDatabaseModule::load_biomass_module_only_names(const string path, const string moduleName)
{
    ifstream iFile;
    int iEnd;
    string name;

	char comments[Constants::COMMENT_SIZE];

	string nameFile = path+ "/" + moduleName;
    openInputFileAndControl(iFile, nameFile.c_str());

    iFile >> name;
    if (name != "SPECIES")
        ErrorMessage("The SPECIES keyword was not found in the solid kinetic database!");

    NC_biomass = 0;
    iEnd = 0;
    while(iEnd == 0)
    {
        iFile >> name;
        if (name == "END") iEnd = 1;
        else               NC_biomass++;
    }
    iFile.seekg(0, ios::beg);

    names_biomass_module = new string[NC_biomass+1];

    iFile.getline(comments, Constants::COMMENT_SIZE);
    for(int i=1;i<=NC_biomass;i++)
        iFile >> names_biomass_module[i];

    iFile.close();
}