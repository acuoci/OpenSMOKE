/***************************************************************************
 *   Copyright (C) 2006-2008 by Alberto Cuoci						       *
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

#include "BzzMath.hpp"
#include "preprocessing/OpenSMOKE_PreProcessorReactingGas.h"
#include <sstream>


// define the ID values to indentify the option
enum { OPT_HELP, OPT_FLAG, OPT_ARG };

CSimpleOpt::SOption parser_options[] =
{
    { OPT_FLAG, _T("--name"),			 SO_REQ_SEP }, // "--input		ARG"
    { OPT_FLAG, _T("--author"),			 SO_REQ_SEP }, // "--author		ARG"
    { OPT_FLAG, _T("--place"),			 SO_REQ_SEP }, // "--place		ARG"
    { OPT_FLAG, _T("--input"),			 SO_REQ_SEP }, // "--input		ARG"
    { OPT_ARG,  _T("--transport"),       SO_REQ_SEP }, // "--transport  ARG"
    { OPT_ARG,  _T("--elements"),        SO_REQ_SEP }, // "--elements   ARG"
    { OPT_ARG,  _T("--Tmin"),            SO_REQ_SEP }, // "--Tmin		ARG"
    { OPT_ARG,  _T("--Tmax"),            SO_REQ_SEP }, // "--Tmin		ARG"

    { OPT_HELP, _T("-?"),    SO_NONE    }, // "-?"
    { OPT_HELP, _T("-help"), SO_NONE    }, // "-help"
    SO_END_OF_OPTIONS                       // END
};

int main(int argc, char* argv[])
{
	ParserClass	parser;
	string      argument;
	string		input_folder_name;
	string		transport_file_name;
	string		elements_file_name;
	OpenSMOKE_PreProcessorReactingGas mix;

    // 0. Parser setup
    parser.setup(argc, argv, parser_options);

    if (parser.parse("--transport", argument))	transport_file_name = argument;
    else										transport_file_name = "tran.bzz";
    if (parser.parse("--elements", argument))	elements_file_name	= argument;
    else										elements_file_name	= "elements.bzz";
    if (parser.parse("--input", argument))		input_folder_name	= argument;
    else										input_folder_name	= "BzzKineticScheme";
    
	if (parser.parse("--Tmin", argument))		
	{
		double			double_number;
		stringstream	number;
		number << argument;
		number >> double_number;
		mix.SetMinimumTemperature(double_number);
	}

	if (parser.parse("--Tmax", argument))		
	{
		double			double_number;
		stringstream	number;
		number << argument;
		number >> double_number;
		mix.SetMaximumTemperature(double_number);
	}

	if (parser.parse("--name", argument))
		mix.SetName(argument);

	if (parser.parse("--author", argument))
		mix.SetAuthorName(argument);

	cout << "" << endl;
	cout << "---------------------------------------------------------------" << endl;
	cout << "                      OpenSMOKE_Converter                      " << endl;
	cout << "                  Version 0.14 - August 2009                   " << endl;
	cout << "             Alberto Cuoci - Politecnico di Milano             " << endl;
	cout << "                    alberto.cuoci@polimi.it                    " << endl;
	cout << "---------------------------------------------------------------" << endl;
	cout << "" << endl;

	mix.PreProcessing("BzzKineticScheme", transport_file_name, elements_file_name, "DEFAULT");

	cout << "Press enter to continue... " << endl;
	getchar();
	return 0;
}

void ParserClass::ShowUsage()
{
    cout << "Usage: OpenSMOKE_Converter [-option] [--option VALUE] [-?] [-help]" << endl;

    cout << "Options: "																						<< endl;
    cout << "   -?                    this help"															<< endl;
    cout << "   -help                 this help"															<< endl;
	cout << "   --author STRING	      author name (default '[not assigned]')"								<< endl;
    cout << "   --elements FILENAME   elements file name (default 'elements.bzz')"							<< endl;
	cout << "   --input FOLDERNAME    folder containing the kinetic scheme (default 'BzzKineticScheme')"	<< endl;
	cout << "   --name STRING	      kinetic scheme name (default '[not assigned]')"						<< endl;
	cout << "   --place STRING	      place name (default '[not assigned]')"								<< endl;
    cout << "   --Tmax DOUBLE         maximum temperature (default 3700 K)"									<< endl;
    cout << "   --Tmin DOUBLE         minimum temperature (default 240 K)"									<< endl;
    cout << "   --transport FILENAME  transport properties file name (default 'tran.bzz')"					<< endl;
    cout << endl;
}
