/***************************************************************************
 *   Copyright (C) 2006-2008 by Alberto Cuoci							   *
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
#include "interfaces/SimpleOpt.h"
#include "addons/OpenSMOKE_InverseKinetics.h"
 
using namespace std;

void ErrorMessage(const string message)
{
    cout << endl;
    cout << "Executable:  OpenSMOKE_InverseKinetics"		<< endl;
    cout << "Error:       "  << message			<< endl;
    cout << "Press a key to continue... "		<< endl;
    getchar();
    exit(-1);
}

// define the ID values to indentify the option
enum { OPT_HELP, OPT_FLAG, OPT_ARG };

CSimpleOpt::SOption parser_options[] =
{
    { OPT_ARG,  _T("--reaction"),        SO_REQ_SEP }, // "-reaction		ARG"
    { OPT_ARG,  _T("--pressure"),        SO_REQ_SEP }, // "-reaction		ARG"
    { OPT_ARG,  _T("--kinetics"),        SO_REQ_SEP }, // "-kinetic			ARG"
    { OPT_ARG,  _T("--global"),          SO_REQ_SEP }, // "-global			ARG"
    { OPT_ARG,  _T("--tmin"),            SO_REQ_SEP }, // "-reaction		ARG"
    { OPT_ARG,  _T("--tmax"),            SO_REQ_SEP }, // "-reaction		ARG"
    { OPT_ARG,  _T("--deltat"),          SO_REQ_SEP }, // "-reaction		ARG"
 
    { OPT_HELP, _T("-?"),     SO_NONE    }, // "-?"
    { OPT_HELP, _T("-help"),  SO_NONE    }, // "-help"
    SO_END_OF_OPTIONS                       // END
};

int main(int argc, char* argv[])
{
	int					iReaction;
	double				pressure;
    string              argument;
    ParserClass         parser;
    
	OpenSMOKE_ReactingGas		mix;
    OpenSMOKE_GlobalKinetics    global;
    OpenSMOKE_InverseKinetics	inversekinetics;


    // 0. Parser setup
    parser.setup(argc, argv, parser_options);


    // 1. Detailed kinetic scheme setup
    if (parser.parse("--kinetics", argument))
            mix.SetupBinary(argument);
    else    ErrorMessage("The --kinetics option is compulsory");


    // 2a. Global kinetic scheme setup
    global.assign_mix(&mix);
    if (parser.parse("--global", argument))
		global.read_from_file(argument);
	else    ErrorMessage("The --global option is compulsory");

	if (parser.parse("--reaction", argument))
		iReaction = atoi(argument.c_str());
	else    ErrorMessage("The --reaction option is compulsory");

	if (parser.parse("--pressure", argument))
		pressure = atof(argument.c_str());
	else    ErrorMessage("The --pressure option is compulsory");

	double tmin   = 300.;
	double tmax   = 3000.;
	double deltat = 100.;

	if (parser.parse("--tmin",   argument))	tmin   = atof(argument.c_str());
	if (parser.parse("--tmax",   argument))	tmax   = atof(argument.c_str());
	if (parser.parse("--deltat", argument))	deltat = atof(argument.c_str());

    // 4. Setup
    inversekinetics.AssignKineticScheme(mix);
    inversekinetics.AssignGlobalKineticScheme(global);

	// 5. Reaction
	inversekinetics.Setup(iReaction, pressure, tmin, tmax, deltat);

    OpenSMOKE_logo("OpenSMOKE_InverseKinetics", "0.1", "May 2009");
	getchar();

    return 0;
}

void ParserClass::ShowUsage()
{
    cout << "Usage: OpenSMOKE_InverseKinetics [-option] [--option STRING] [-?] [--help]" << endl;

    cout << "Options: "																				<< endl;
    cout << "   -?                      this help"													<< endl;
    cout << "   --pressure DOUBLE       pressure in atm"											<< endl;
    cout << "   --tmin DOUBLE           minimum temperature in K"									<< endl;
    cout << "   --tmax DOUBLE           maximum temperature in K"									<< endl;
    cout << "   --deltat DOUBLE         delta temperature in K"										<< endl;
    cout << "   --reaction INTEGER      reaction index"												<< endl;
    cout << "   --global FILENAME       global kinetics file name"									<< endl;
    cout << "   --kinetics FOLDERNAME   detailed kinetic scheme folder (default 'CRECK')"			<< endl;
    cout << "   -help                   this help"													<< endl;
    cout << endl;
}
