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
#include "SimpleOpt.h"
#include "preprocessing/OpenSMOKE_PreProcessorReactingGas.h"

void ErrorMessage(const string message)
{ 
    cout << endl;
    cout << "Executable:  OpenSMOKE_PreProcessorTransportSensitivity"		<< endl;
    cout << "Error:       "  << message			<< endl;
    cout << "Press a key to continue... "		<< endl;
    getchar();
    exit(-1);
}

// define the ID values to indentify the option
enum { OPT_HELP, OPT_FLAG, OPT_ARG };

CSimpleOpt::SOption parser_options[] =
{
    { OPT_ARG,  _T("--kinetics"),        SO_REQ_SEP }, // "-kinetic			ARG"
    { OPT_ARG,  _T("--transport"),       SO_REQ_SEP }, // "-transport		ARG"
    { OPT_ARG,  _T("--elements"),        SO_REQ_SEP }, // "-thermo			ARG"
    { OPT_ARG,  _T("--epsilon"),         SO_REQ_SEP }, // "-thermo			ARG"

    { OPT_HELP, _T("-?"),				 SO_NONE    }, // "-?"
    { OPT_HELP, _T("-help"),			 SO_NONE    }, // "--help"
    SO_END_OF_OPTIONS                       // END
};

int main(int argc, char* argv[])
{
	OpenSMOKE_logo("OpenSMOKE_PreProcessorTransportSensitivity", "0.13", "July 2010");

    string              argument;
    ParserClass         parser;

	string pathName;
	string fileNameElements;
	string fileNameTransport;
	string epsilon;

	OpenSMOKE_PreProcessorReactingGas preprocessor;

    // 0. Parser setup
    parser.setup(argc, argv, parser_options);


    // 1. Kinetic scheme setup
    if (parser.parse("--kinetics", argument))
            pathName = argument;
    else    ErrorMessage("The --kinetics option is compulsory");

    // 2. Thermodynamic file setup
    if (parser.parse("--elements", argument))
            fileNameElements = argument;
    else    ErrorMessage("The --elements option is compulsory");

    // 3. Transport properties file setup
    if (parser.parse("--transport", argument))
            fileNameTransport = argument;
    else    ErrorMessage("The --transport option is compulsory");

    // 3. Transport properties file setup
    if (parser.parse("--epsilon", argument))
		epsilon = argument;
    else    ErrorMessage("The --epsilon option is compulsory");

	preprocessor.PreProcessingTransportSensitivity(pathName, fileNameTransport, fileNameElements, atof(epsilon.c_str()));

    OpenSMOKE_logo("OpenSMOKE_PreProcessorTransportSensitivity", "0.13", "July 2010");

	OpenSMOKE_LennardJonesSensitivityCoefficientsManager os;
	os.Setup("TransportSensitivity.out");

	cout << "Press enter to continue..." << endl;
	getchar();

    return 0;
}

void ParserClass::ShowUsage()
{
    cout << "Usage: OpenSMOKE_PreProcessorTransportSensitivity [-option] [--option STRING] [-?] [--help]" << endl;

    cout << "Options: "																				<< endl;
    cout << "    -?                     this help"													<< endl;
    cout << "    -help                  this help"													<< endl;
 
	cout << "   --kinetics FILENAME     detailed kinetic scheme file        "	<< endl;
    cout << "   --elements FILENAME     thermodynamic data file             "	<< endl;
    cout << "   --transport FILENAME    transport properties file           "	<< endl;

    cout << endl;
}
