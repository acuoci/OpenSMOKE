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
    { OPT_ARG,  _T("--model"),        SO_REQ_SEP }, // "--model		ARG"
    { OPT_ARG,  _T("--list"),         SO_REQ_SEP }, // "--list		ARG"
    { OPT_ARG,  _T("--kinetics"),     SO_REQ_SEP }, // "--kinetics	ARG"
    { OPT_ARG,  _T("--database"),     SO_REQ_SEP }, // "--database	ARG"
 
    { OPT_HELP, _T("-?"),     SO_NONE    }, // "-?"
    { OPT_HELP, _T("--help"), SO_NONE    }, // "--help"
    
	SO_END_OF_OPTIONS                       // END
};

int main(int argc, char* argv[])
{
	int    iModel;
	string list_name;
	string kinetics_name;
	string database_name;

	string argument;

    MyOpenSMOKE_SolidRegression	regression;
	OpenSMOKE_CharKineticScheme char_kinetics;
    ParserClass					parser;

    // 0. Parser setup
	parser.setup(argc, argv, parser_options);

	if (parser.parse("--model", argument))
		iModel = atoi(argument.c_str());
	else    ErrorMessage("The --model option is compulsory");

	if (parser.parse("--list", argument))
		list_name = argument;
	else    ErrorMessage("The --list option is compulsory");

	if (parser.parse("--kinetics", argument))
		kinetics_name = argument;
	else    ErrorMessage("The --kinetics option is compulsory");

	if (parser.parse("--database", argument))
		database_name = argument;
	else    ErrorMessage("The --database option is compulsory");


	// Kinetic scheme setup
	char_kinetics.ReadDatabase(database_name);
	char_kinetics.ReadKineticScheme(kinetics_name);
	

	// Regression setup
	regression.SetKinetics(&char_kinetics);
	regression.Setup(list_name);


	BzzVector kappa0(4);
	BzzVector Ea(4);

	kappa0[1] = 2.52e3;		// [m5/mol2/s]
	kappa0[2] = 3.64e10;	// [m5/mol2/s]
	kappa0[3] = 1.84e15;	// [m5/mol2/s]
	kappa0[4] = 2.53e11;	// [1/s]

	Ea[1] = 42.e6;			// [J/kmol]
	Ea[2] = 134.e6;			// [J/kmol]
	Ea[3] = 226.e6;			// [J/kmol]
	Ea[4] = 251.e6;			// [J/kmol]

	kappa0[1] *= 1e6;		// [m5/kmol2/s]
	kappa0[2] *= 1e6;		// [m5/kmol2/s]
	kappa0[3] *= 1e6;		// [m5/kmol2/s]
	kappa0[4] *= 1.0;		// [1/s]

/*
kappa0[1] = 3.376306e+009;
kappa0[2] = 3.198205e+016;
kappa0[3] = 1.592472e+021;
kappa0[4] = 5.869724e+011;	
Ea[1] = 5.875736e+007;
Ea[2] = 1.622735e+008;
Ea[3] = 2.154851e+008;
Ea[4] = 2.162609e+008;
*/
	int i;
	BzzVector bFirstGuess(8);
	for(i=1;i<=4;i++)	bFirstGuess[i] = kappa0[i];
	for(i=1;i<=4;i++)	bFirstGuess[4+i] = Ea[i];
	

	for(i=1;i<=8;i++)
		cout << bFirstGuess[i] << endl;

	getchar();

	regression.Run(iModel, bFirstGuess);

	// Logo
    OpenSMOKE_logo("OpenSMOKE_SolidRegression", "0.1", "May 2009");
	getchar();

    return 0;
}

void ParserClass::ShowUsage()
{
    cout << "Usage: OpenSMOKE_SolidRegression [-option] [--option STRING] [-?] [--help]" << endl;

    cout << "Options: "																				<< endl;
    cout << "   -?                      this help"													<< endl;

	cout << "   --model    INTEGER         model index"     											<< endl;
    cout << "   --list     FILENAME        file name containing the list of files"						<< endl;
    cout << "   --kinetics FILENAME        file name containing the kinetic scheme"					<< endl;
    cout << "   --database FILENAME        file name containing the species database"					<< endl;

	cout << "   -help                   this help"													<< endl;
    cout << endl;
}
