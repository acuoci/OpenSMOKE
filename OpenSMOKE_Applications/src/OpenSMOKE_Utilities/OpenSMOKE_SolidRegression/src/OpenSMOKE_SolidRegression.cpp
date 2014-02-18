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
#include "SimpleOpt.h"
#include "MyOpenSMOKE_SolidRegression.h"
#include "MyOpenSMOKE_CharKineticScheme.h"

using namespace std;

void ErrorMessage(const string message)
{
    cout << endl;
    cout << "Executable:  OpenSMOKE_SolidRegression"		<< endl;
    cout << "Error:       "  << message			<< endl;
    cout << "Press a key to continue... "		<< endl;
    getchar();
    exit(-1);
}

// define the ID values to indentify the option
enum { OPT_HELP, OPT_FLAG, OPT_ARG };

CSimpleOpt::SOption parser_options[] =
{
    { OPT_ARG,  _T("--input"),        SO_REQ_SEP }, // "--input		ARG"
    { OPT_ARG,  _T("--kinetics"),     SO_REQ_SEP }, // "--kinetics	ARG"
    { OPT_ARG,  _T("--database"),     SO_REQ_SEP }, // "--database	ARG"

	{ OPT_ARG,  _T("-char"),		  SO_NONE }, // "--database	ARG"
	{ OPT_ARG,  _T("-bio"),			  SO_NONE }, // "--database	ARG"
 
    { OPT_HELP, _T("-?"),     SO_NONE    }, // "-?"
    { OPT_HELP, _T("--help"), SO_NONE    }, // "--help"
    
	SO_END_OF_OPTIONS                       // END
};

int main(int argc, char* argv[])
{
	string input_name;
	string kinetics_name;
	string database_name;

	string argument;
	int iModel;

    MyOpenSMOKE_SolidRegression	regression;
	OpenSMOKE_CharKineticScheme char_kinetics;
	OpenSMOKE_ReactingGas		bio_kinetics;
    ParserClass					parser;

    // 0. Parser setup
	parser.setup(argc, argv, parser_options);

	if (parser.parse("--kinetics", argument))
		kinetics_name = argument;
	else    ErrorMessage("The --kinetics option is compulsory");

	if (parser.parse("--input", argument))
		input_name = argument;
	else    ErrorMessage("The --input option is compulsory");

	iModel = 0;
	if (parser.parse("-char"))		iModel = 1;
	else if (parser.parse("-bio"))	iModel = 2;
	else ErrorMessage("The -char || -bio options are compulsory...");
	
	
	// Kinetic scheme setup
	if (iModel == 1)
	{
		if (parser.parse("--database", argument))
			database_name = argument;
		else ErrorMessage("The --database option is compulsory...");

		regression.AssignKinetics(&char_kinetics);
		char_kinetics.ReadDatabase(database_name);
		char_kinetics.ReadKineticScheme(kinetics_name);
	}
	
	if (iModel == 2)
	{
		bio_kinetics.SetupBinary(kinetics_name);
		regression.AssignKinetics(&bio_kinetics);

		if (parser.parse("--database", argument))
			ErrorMessage("The --database option is not required...");
	}

	// Regression setup
	regression.Setup(input_name);

	// Regression running
	double start_time = BzzGetCpuTime();
		if (iModel == 1)	regression.Run_Char();
		if (iModel == 2)	regression.Run_DSmoke_Bio();
	double end_time = BzzGetCpuTime();
	cout << "Total cpu time: " << end_time - start_time << " s" << endl;

	// Logo
    OpenSMOKE_logo("OpenSMOKE_SolidRegression", "0.4", "December 2009");
	getchar();

    return 0;
}

void ParserClass::ShowUsage()
{
    cout << "Usage: OpenSMOKE_SolidRegression [-option] [--option STRING] [-?] [--help]" << endl;

    cout << "Options: "																				<< endl;
    cout << "   -?                      this help"													<< endl;

    cout << "   --input    FILENAME        file name of input options"						<< endl;
    cout << "   --kinetics FILENAME        file name containing the kinetic scheme"					<< endl;
    cout << "   --database FILENAME        file name containing the species database"					<< endl;
 
	cout << "   -char      NONE            char regression"					<< endl;
	cout << "   -bio       NONE            bio regression"					<< endl;

	cout << "   -help                   this help"													<< endl;
    cout << endl;
}
