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

#include "BzzMathAdvanced.hpp"
#include "OpenSMOKE.hpp"
#include "solid/OpenSMOKE_SolidMixture.h"
//#include "idealreactors/batch_solid/OpenSMOKE_SolidBatch.h"
#include "XOpenSMOKE_SolidBatch.h"

using namespace std;


// define the ID values to indentify the option
enum { OPT_HELP, OPT_FLAG, OPT_ARG };

CSimpleOpt::SOption parser_options[] =
{
    { OPT_FLAG, _T("--noverbose"),       SO_NONE    }, // "-verbose"
    { OPT_ARG,  _T("--input"),           SO_REQ_SEP }, // "-input   ARG"
    { OPT_ARG,  _T("--output"),          SO_REQ_SEP }, // "-output  ARG"
    { OPT_ARG,  _T("--kinetics"),        SO_REQ_SEP }, // "-kinetic ARG"
    { OPT_ARG,  _T("--global"),          SO_REQ_SEP }, // "-global  ARG"
    { OPT_ARG,  _T("--inlet"),           SO_REQ_SEP }, // "-inlet   ARG"
	{ OPT_ARG,  _T("--soot2EModel"),     SO_REQ_SEP }, // "-soot2EModel		ARG"

    { OPT_HELP, _T("-?"),     SO_NONE    }, // "-?"
    { OPT_HELP, _T("--help"), SO_NONE    }, // "--help"
    SO_END_OF_OPTIONS                       // END
};

int main(int argc, char* argv[])
{
	string              argument;
    ParserClass         parser;
    
	OpenSMOKE_ReactingGas		mix;
    OpenSMOKE_GlobalKinetics    global;
	OpenSMOKE_SolidBatch		solid_batch;
	OpenSMOKE_GasStream			inlet;
	OpenSMOKE_GasStream			outlet;
	OpenSMOKE_2EModel			soot2EModel;
	OpenSMOKE_SolidMixture		solid;
	

    // 0. Parser setup
    parser.setup(argc, argv, parser_options);


    // 1a. Detailed kinetic scheme setup
	mix.SetLookUpTables();
	mix.SetMixViscositySimplified();
    if (parser.parse("--kinetics", argument))
            mix.Setup(argument);
    else    mix.Setup("CRECK");

    // 1b. Solid kinetic scheme setup
	int nModules = 2;
	string *names_modules = new string[nModules+1];
	names_modules[1] = "Biomass_0810.inp";
	names_modules[2] = "Char_0810.inp";
	solid.setup("SolidDatabase", nModules, names_modules, &mix);
	

    // 2a. Global kinetic scheme setup
    global.assign_mix(&mix);
    if (parser.parse("--global", argument))
		global.read_from_file(argument);

    // 2b. Soot 2E Model
//	soot2EModel.assign_mixture(mix);
//	if (parser.parse("--soot2EModel", argument))
//		soot2EModel.setupFromFile(argument);
	

    // 3. Inlet stream setup
    inlet.AssignKineticScheme(mix);
    if (parser.parse("--inlet", argument))
            inlet.DefineFromFile(argument);
    else    inlet.DefineFromFile("InletStream.inp");
	inlet.VideoSummary();


    // 4. Batch Setup
    solid_batch.AssignKineticScheme(mix);
    solid_batch.AssignSolidKineticScheme(solid);
    solid_batch.AssignGlobalKineticScheme(global);
	solid_batch.AssignSoot2EModel(soot2EModel);
    solid_batch.AssignInletFlows(inlet);
    if (parser.parse("--input", argument))
            solid_batch.DefineFromFile(argument);
    else    solid_batch.DefineFromFile("SolidBatch.inp");


    // 5. Batch Video Summary
    solid_batch.VideoSummary();

    // 6. Batch Solution
    solid_batch.Solve();

    // 7. Video Solution
    solid_batch.VideoFinalResult();

	// 8. Outlet Stream
	solid_batch.OutletStream(outlet);
	outlet.VideoSummary();

	// 9. Mass and Energy balances
	solid_batch.SummaryOnFile("Output/Summary.out");
	solid_batch.MassAnalysis(outlet);
	solid_batch.EnergyAnalysis(outlet);

    OpenSMOKE_logo("OpenSMOKE_SolidBatch", "0.0.12", "December 2008");

	cout << "Last update: 20081226 6pm" << endl;

    return 0;
}

void ParserClass::ShowUsage()
{
    cout << "Usage: OpenSMOKE_SolidBatch [-option] [--option STRING] [-?] [--help]"					<< endl;

    cout << "Options: "																				<< endl;
    cout << "   -?                      this help"													<< endl;
    cout << "   --global FILENAME       global kinetics file name"									<< endl;
    cout << "   --inlet FILENAME        stream definition file name (default 'InletStream.inp')"	<< endl;
    cout << "   --input FILENAME        input file for reactor setup (default 'SolidBatch.inp')"	<< endl;
    cout << "   --kinetics FOLDERNAME   detailed kinetic scheme folder (default 'CRECK')"			<< endl;
    cout << "   -help                   this help"													<< endl;
    cout << "   --output FOLDERNAME     output folder for setting up (default 'Output')"			<< endl;
    cout << "   -noverbose              no file nor video output"									<< endl;
    cout << "   --soot2EModel FILENAME  global kinetics file name"									<< endl;
    cout << endl;
}
