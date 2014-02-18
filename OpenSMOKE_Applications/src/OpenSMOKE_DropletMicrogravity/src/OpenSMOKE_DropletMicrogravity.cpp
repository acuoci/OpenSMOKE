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
#include "OpenSMOKE.hpp"
#include "droplet/OpenSMOKE_DropletMicrogravity.h"
#include "droplet/OpenSMOKE_DropletMicrogravity_DataManager.h"

#include "liquid/OpenSMOKE_LiquidProperties_Database.h"
#include "liquid/OpenSMOKE_LiquidSpecies.h"

using namespace std;

void ErrorMessage(const string message)
{
    cout << endl;
    cout << "Executable:  OpenSMOKE_DropletMicrogravity"	<< endl;
    cout << "Error:       "  << message			<< endl;
    cout << "Press a key to continue... "		<< endl;
    getchar();
    exit(-1);
}

// define the ID values to indentify the option
enum { OPT_HELP, OPT_FLAG, OPT_ARG };

CSimpleOpt::SOption parser_options[] =
{
    { OPT_FLAG, _T("-noverbose"),		 SO_NONE    }, // "-noverbose"
	{ OPT_FLAG, _T("--2e"),				 SO_REQ_SEP }, // "--2e"
    { OPT_ARG,  _T("--backup"),          SO_REQ_SEP }, // "--input    ARG"
    { OPT_ARG,  _T("--liquid"),          SO_REQ_SEP }, // "--input    ARG"
    { OPT_ARG,  _T("--input"),           SO_REQ_SEP }, // "--input    ARG"
    { OPT_ARG,  _T("--output"),          SO_REQ_SEP }, // "--output   ARG"
    { OPT_ARG,  _T("--kinetics"),        SO_REQ_SEP }, // "--kinetic  ARG"
    { OPT_ARG,  _T("--global"),          SO_REQ_SEP }, // "--global   ARG"
    { OPT_ARG,  _T("--schedule"),        SO_REQ_SEP }, // "--schedule ARG"

    { OPT_HELP, _T("-?"),     SO_NONE    }, // "-?"
    { OPT_HELP, _T("-help"),  SO_NONE    }, // "-help"
    SO_END_OF_OPTIONS                       // END
};


// ---------------------------------------------------------------------------------
// MAIN
// ---------------------------------------------------------------------------------
int main(int argc, char* argv[])
{
	string								argument;
    ParserClass							parser;

	OpenSMOKE_DropletMicrogravity				droplet;
	OpenSMOKE_DropletMicrogravity_DataManager	data;
	
	OpenSMOKE_GlobalKinetics			global;
	OpenSMOKE_ReactingGas				mix;

	// Prepare
	droplet.Assign(&mix);
	droplet.Assign(&global);
	droplet.Assign(&data);
	data.Assign(&mix);
	data.Assign(&droplet);

	// 0. Parser setup
    parser.setup(argc, argv, parser_options);

	OpenSMOKE_LiquidProperties_Database liquid_database;
    if (parser.parse("--liquid", argument))
	{
		liquid_database.ReadFromFolder(argument);
		data.Assign(&liquid_database);
	}

//	liquid_database.LoadFromFile("LiquidProperties/Database.liq");
	
    // 1. Detailed kinetic scheme setup
    if (parser.parse("--kinetics", argument))
		mix.SetupBinary(argument);
    else ErrorMessage("The kinetics must be specified: --kinetics FOLDERNAME");

	// 2. Global kinetic scheme setup
    global.assign_mix(&mix);
    if (parser.parse("--global", argument))
	{
		global.read_from_file(argument);
		data.iGlobalKinetics = true;
	}

	if (parser.parse("--input", argument))
		data.ReadFromFile(argument);
	else ErrorMessage("The input file must be specified: --input FILENAME");

	if (parser.parse("--output", argument))
		data.SetOutputFolderName(argument);

/*
	// 4. Input data
    if (parser.parse("--backup", argument))
	{
		data.iBackUp = true;
		flamelet.FoldersAndFilesManager(argument);
		flamelet.RecoverFromBackUp(flamelet.nameFileBackupInputData);
	}

	if (parser.parse("--2e", argument))
	{
		data.twoEquation_file_name = argument;
		data.i2E = true;
	}
	else data.i2E = false;
*/

	// Solution

	double tStart = BzzGetCpuTime();

	droplet.Setup();
	droplet.Run();

	cout << "Total time: " << BzzGetCpuTime() - tStart << " s" << endl;

	OpenSMOKE_logo("OpenSMOKE_DropletMicrogravity", "0.12", "March 2011");

	return(0);
}

void ParserClass::ShowUsage()
{
    cout << "Usage: OpenSMOKE_DropletMicrogravity [-option] [--option FILE] [-?] [-help]" << endl;

    cout << "Options: "																					<< endl;
    cout << "   -?                    this help"														<< endl;
    cout << "   -help				  this help"														<< endl;
    cout << "   -noverbose            no file nor video output"											<< endl;
	cout << "   --2e FILENAME		  2e simulation"													<< endl;
	cout << "   --liquid FILENAME     liquid database file name"        										<< endl;
	cout << "   --backup FILENAME     start from backup"        										<< endl;
	cout << "   --global FILENAME     global kinetics file name"										<< endl;
    cout << "   --input FILENAME      input file for setting up"										<< endl;
    cout << "   --kinetics FOLDER     detailed kinetic scheme folder"									<< endl;
    cout << "   --output FOLDER	      output folder for setting up (default 'Output')"					<< endl;
    cout << "   --schedule FILENAME   input file for operations (default 'Input/Operations.inp')"		<< endl;
	cout << endl;
}
