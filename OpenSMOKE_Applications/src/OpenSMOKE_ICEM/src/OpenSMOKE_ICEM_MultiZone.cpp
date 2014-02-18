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

#include "BzzMath.hpp"
#include "OpenSMOKE.hpp"
#include "idealreactors/icem/OpenSMOKE_ICEM_MultiZone.h"

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
    { OPT_HELP, _T("-help"),  SO_NONE    }, // "-help"
    SO_END_OF_OPTIONS                       // END
};

int main(int argc, char* argv[])
{
//	bzzOpenMP = 0;

	string              argument;
    ParserClass         parser;
    
	OpenSMOKE_ReactingGas		mix;
    OpenSMOKE_GlobalKinetics    global;
	OpenSMOKE_ICEM_MultiZone	icem;
	OpenSMOKE_GasStream			inlet;
	OpenSMOKE_GasStream			outlet;
	OpenSMOKE_2EModel			soot2EModel;
	

    // 0. Parser setup
    parser.setup(argc, argv, parser_options);


    // 1. Detailed kinetic scheme setup
    if (parser.parse("--kinetics", argument))
            mix.SetupBinary(argument);
//    else    ErrorMessage("--kinetics option is compulsory");


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

	// Output folder
    if (parser.parse("--output", argument))
		icem.SetOutputFolder(argument);

    // 4. Batch Setup
    icem.AssignKineticScheme(mix);
    icem.AssignGlobalKineticScheme(global);
	icem.AssignSoot2EModel(soot2EModel);
    icem.AssignInletFlows(inlet);
    if (parser.parse("--input", argument))
            icem.DefineFromFile(argument);
    else    icem.DefineFromFile("Batch.inp");


    // 5. Batch Video Summary
    icem.VideoSummary();

    // 6. Batch Solution
	double timeStart = BzzGetCpuTime();
    icem.Solve();
	double timeEnd = BzzGetCpuTime();

    // 7. Video Solution
    icem.VideoFinalResult();

	// 8. Outlet Stream
//	icem.OutletStream(outlet);
//	outlet.VideoSummary();

	// 9. Mass and Energy balances
//	icem.SummaryOnFile();
//	icem.MassAnalysis(outlet);
//	icem.EnergyAnalysis(outlet);

	cout << "Total time for engine simulation: " << timeEnd-timeStart << " " << endl;

    OpenSMOKE_logo("OpenSMOKE_ICEM", "0.2", "November 2010");

	cout << "Last update: 20101117 6pm" << endl;
	cout << endl;

	cout << "Press enter to continue..." << endl;
	getchar();

    return 0;
}

void ParserClass::ShowUsage()
{
    cout << "Usage: OpenSMOKE_ICEM [-option] [--option STRING] [-?] [--help]"						<< endl;

    cout << "Options: "																				<< endl;
    cout << "   -?                      this help"													<< endl;
    cout << "   --global FILENAME       global kinetics file name"									<< endl;
    cout << "   --inlet FILENAME        stream definition file name (default 'InletStream.inp')"	<< endl;
    cout << "   --input FILENAME        input file for reactor setup (default 'PFR.inp')"			<< endl;
    cout << "   --kinetics FOLDERNAME   detailed kinetic scheme folder (default 'CRECK')"			<< endl;
    cout << "   -help                   this help"													<< endl;
    cout << "   --output FOLDERNAME     output folder for setting up (default 'Output')"			<< endl;
    cout << "   -noverbose              no file nor video output"									<< endl;
    cout << "   --soot2EModel FILENAME  global kinetics file name"									<< endl;
    cout << endl;
}
