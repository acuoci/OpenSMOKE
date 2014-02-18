/***************************************************************************
 *   Copyright (C) 2009 by Alberto Cuoci							   *
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
#include "interfaces/SimpleOpt.h"
#include "idealreactors/pfr/OpenSMOKE_PFR.h"
#include "preprocessing_chemkin/OpenSMOKE_CHEMKINInterpreter.h"

void ErrorMessage(const string message)
{ 
    cout << endl;
    cout << "Executable:  OpenSMOKE_ChemkinInterpreter"		<< endl;
    cout << "Error:       "  << message			<< endl;
    cout << "Press a key to continue... "		<< endl;
    getchar();
    exit(-1);
}

// define the ID values to indentify the option
enum { OPT_HELP, OPT_FLAG, OPT_ARG };

CSimpleOpt::SOption parser_options[] =
{
    { OPT_FLAG, _T("-no-verbose"),		 SO_NONE    },	// "-no-verbose"
    { OPT_FLAG, _T("-bzz-output"),		 SO_NONE    },	// "-bzz-output"
    { OPT_FLAG, _T("-soot-mode"),					SO_NONE    },	// "-soot-mode"
    { OPT_FLAG, _T("-lennard-jones-sensitivity"),	SO_NONE    },	// "-lennard-jones-sensitivity"
    
	{ OPT_ARG,  _T("--output"),          SO_REQ_SEP },	// "-output			ARG"
    { OPT_ARG,  _T("--kinetics"),        SO_REQ_SEP },	// "-kinetic		ARG"
    { OPT_ARG,  _T("--transport"),       SO_REQ_SEP },	// "-transport		ARG"
    { OPT_ARG,  _T("--thermo"),          SO_REQ_SEP },	// "-thermo			ARG"
    { OPT_ARG,  _T("--surface"),         SO_REQ_SEP },	// "-surface			ARG"
    { OPT_ARG,  _T("--reduced-list"),    SO_REQ_SEP },	// "-surface			ARG"

    { OPT_FLAG, _T("--name"),			 SO_REQ_SEP },	// "--input				ARG"
    { OPT_FLAG, _T("--author"),			 SO_REQ_SEP },	// "--author			ARG"
    { OPT_ARG,  _T("--tmin"),            SO_REQ_SEP },	// "--tmin				ARG"
    { OPT_ARG,  _T("--tmax"),            SO_REQ_SEP },	// "--tmin				ARG"
    { OPT_ARG,  _T("--fitting-points"),  SO_REQ_SEP },	// "--fitting-points	ARG"

    { OPT_HELP, _T("-?"),				 SO_NONE    },	// "-?"
    { OPT_HELP, _T("-help"),			 SO_NONE    },	// "--help"
    SO_END_OF_OPTIONS									// END
};

int main(int argc, char* argv[])
{
	OpenSMOKE_logo("OpenSMOKE_ChemkinInterpreter", "0.2", "November 2010");

    string              argument;
    ParserClass         parser;

	OpenSMOKE_CHEMKINInterpreter interpreter;

    // 0. Parser setup
    parser.setup(argc, argv, parser_options);


    // 1. Kinetic scheme setup
    if (parser.parse("--kinetics", argument))
            interpreter.SetKineticFileName(argument);
    else    ErrorMessage("The --kinetics option is compulsory");

    // 2. Thermodynamic file setup
    if (parser.parse("--thermo", argument))
            interpreter.SetThermodynamicFileName(argument);
    else    ErrorMessage("The --thermo option is compulsory");

    // 3. Transport properties file setup
    if (parser.parse("--transport", argument))
            interpreter.SetTransportFileName(argument);

    // 2. Thermodynamic file setup
    if (parser.parse("--surface", argument))
		interpreter.SetSurfaceKineticFileName(argument);

    // 4. Output folder name
    if (parser.parse("--output", argument))
            interpreter.SetOutputFolderName(argument);
    else    ErrorMessage("The --output option is compulsory");

    // 4. Output folder name
    if (parser.parse("--reduced-list", argument))
            interpreter.SetReducedListFileName(argument);

	// 5. Minimum fitting temperature
	if (parser.parse("--tmin", argument))		
	{
		double			double_number;
		stringstream	number;
		number << argument;
		number >> double_number;
		interpreter.SetMinimumTemperature(double_number);
	}

	// 6. Maximum fitting temperature
	if (parser.parse("--tmax", argument))		
	{
		double			double_number;
		stringstream	number;
		number << argument;
		number >> double_number;
		interpreter.SetMaximumTemperature(double_number);
	}

	// 6. Number of fitting points
	if (parser.parse("--fitting-points", argument))		
	{
		int			int_number;
		stringstream	number;
		number << argument;
		number >> int_number;
		interpreter.SetFittingPoints(int_number);
	}

	// 7. Kinetic Scheme name
	if (parser.parse("--name", argument))
		interpreter.SetName(argument);

	// 8. Author name
	if (parser.parse("--author", argument))
		interpreter.SetAuthorName(argument);
	
	// 10. Bzz output
	if (parser.parse("-bzz-output"))
		interpreter.SetBuzziMode();

	// 11. Verbose mode
	if (parser.parse("-no-verbose"))
		interpreter.SetNoVerboseMode();
	else
		interpreter.SetVerboseMode();

	// 12. Verbose mode
	if (parser.parse("-soot-mode"))
		interpreter.SetSootMode();

	// 13. Lennard-Jones parameters
	if (parser.parse("-lennard-jones-sensitivity"))
		interpreter.SetBuzziMode();

	interpreter.Run();

    OpenSMOKE_logo("OpenSMOKE_CHEMKINInterpreter", "0.3", "January 2014");

	// Option: sensitivity to lennard-jones parameters
	if (parser.parse("-lennard-jones-sensitivity"))
	{
		OpenSMOKE_PreProcessorReactingGas preprocessor;

		double epsilon = 0.02;
		if (parser.parse("--epsilon", argument))	epsilon = atof(argument.c_str());

		preprocessor.PreProcessingTransportSensitivity(interpreter.nameOutputFolder(), interpreter.nameTransportFile(), interpreter.nameOutputFolder() + "/Elements.out", epsilon);

		OpenSMOKE_logo("OpenSMOKE_PreProcessorTransportSensitivity", "0.13", "July 2010");
	}

	cout << "Press enter to continue..." << endl;
	getchar();

    return 0;
}

void ParserClass::ShowUsage()
{
    cout << "Usage: OpenSMOKE_ChemkinInterpreter.exe [-option] [--option STRING] [-?] [-help]" << endl;

    cout << "Options: "																						<< endl;
    cout << "    -?                           this help"													<< endl;
    cout << "    -help                        this help"													<< endl;
    cout << "    -bzz-output                  bzz format"													<< endl;
    cout << "    -no-verbose                  no additional files"											<< endl;
    cout << "    -soot-mode                   soot mode"													<< endl;
    cout << "    -lennard-jones-sensitivity   additional data for sensitivity to lennard jones parameters "	<< endl;
 
	cout << "   --kinetics FILENAME      detailed kinetic scheme file          "	<< endl;
    cout << "   --thermo FILENAME        thermodynamic data file               "	<< endl;
    cout << "   --transport FILENAME     transport properties file             "	<< endl;
    cout << "   --surface FILENAME       surface kinetic scheme file             "	<< endl;
    cout << "   --output FOLDERNAME      output folder                         "	<< endl;
    cout << "   --reduced-list FILENAME  output folder                         "	<< endl;
	cout << "   --author STRING	         author name                           "	<< endl;
	cout << "   --name STRING	         kinetic scheme name                   "	<< endl;
    cout << "   --tmax DOUBLE            maximum temperature (default 3600 K)  "	<< endl;
    cout << "   --tmin DOUBLE            minimum temperature (default 300 K)   "	<< endl;
    cout << "   --fitting-points INT     number of fitting points (default 50) "	<< endl;

    cout << endl;
}
