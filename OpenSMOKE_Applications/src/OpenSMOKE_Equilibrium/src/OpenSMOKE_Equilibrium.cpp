/***************************************************************************
 *   Copyright (C) 2009 by Alberto Cuoci								   *
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

#include "OpenSMOKE.hpp"
#include "OpenSMOKE_EquilibriumModule.h"

void ErrorMessage(const string message)
{
    cout << endl;
    cout << "Executable:  OpenSMOKE_Equilibrium"	<< endl;
    cout << "Error:       "  << message				<< endl;
    cout << "Press a key to continue... "			<< endl;
    getchar();
    exit(-1);
}

// define the ID values to indentify the option
enum { OPT_HELP, OPT_FLAG, OPT_ARG };

CSimpleOpt::SOption parser_options[] =
{
    { OPT_FLAG, _T("-noverbose"),		 SO_NONE    }, // "-noverbose"
    { OPT_ARG,  _T("--input"),           SO_REQ_SEP }, // "--input		ARG"
    { OPT_ARG,  _T("--kinetics"),        SO_REQ_SEP }, // "--kinetics	ARG"
    { OPT_ARG,  _T("--output"),		     SO_REQ_SEP }, // "--output		ARG"

    { OPT_HELP, _T("-?"),     SO_NONE    }, // "-?"
    { OPT_HELP, _T("-help"),  SO_NONE    }, // "-help"
    SO_END_OF_OPTIONS                       // END
};


// ---------------------------------------------------------------------------------
// MAIN
// ---------------------------------------------------------------------------------
int main(int argc, char* argv[])
{
	RemoveWarningWindow();

	string								argument;
    ParserClass							parser;

	OpenSMOKE_EquilibriumModule			equilibrium;
	OpenSMOKE_ReactingGas				mix;


	// Prepare
	equilibrium.Assign(&mix);

	// 0. Parser setup
    parser.setup(argc, argv, parser_options);

    // 1. Detailed kinetic scheme setup
    if (parser.parse("--kinetics", argument))
		mix.SetupBinary(argument);
    else ErrorMessage("The kinetic mechanism must be specified: --kinetics NAMEFOLDER");

	// 2. Data Manager Setup
	if (parser.parse("--input", argument))
		equilibrium.SetInputFileName(argument);
	else ErrorMessage("The input file must be specified: --input NAMEFILE");


	double tStart = BzzGetCpuTime();

		equilibrium.Run();
	
	double tEnd = BzzGetCpuTime();

	// Finalize
	cout << "Total time for solution: " << tEnd - tStart << " s" << endl;
	OpenSMOKE_logo("OpenSMOKE_Equilibrium", "0.2", "November 2010");

	return 1;
}

void ParserClass::ShowUsage()
{
    cout << "Usage: OpenSMOKE_Equilibrium [-option] [--option FILE] [-?] [-help]" << endl;

    cout << "Options: "																					<< endl;
    cout << "    -?									this help"													<< endl;
    cout << "    -help								this help"													<< endl;
    cout << "    -noverbose							no file nor video output"									<< endl;
    cout << "   --input FILENAME					input file for setting up (default 'Input/Flamelet.inp')"	<< endl;
    cout << "   --kinetics FOLDER					detailed kinetic scheme folder (default 'CRECK')"			<< endl;
    cout << "   --output FOLDER						output folder for setting up (default 'Output')"				<< endl;
	cout << endl;
}