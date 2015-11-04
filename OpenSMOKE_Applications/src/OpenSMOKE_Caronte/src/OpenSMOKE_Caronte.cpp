/***************************************************************************
 *   Copyright (C) 2006 by Alberto Cuoci								   *
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
#include "idealreactors/flame1d/OpenSMOKE_Flame1D.h"

void ErrorMessage(const string message)
{
    cout << endl;
    cout << "Executable:  OpenSMOKE_Flame1D"	<< endl;
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
	{ OPT_FLAG, _T("-premix"),			 SO_NONE    }, // "-premix"
	{ OPT_FLAG, _T("-opposed"),			 SO_NONE    }, // "-opposed"
	{ OPT_FLAG, _T("-twin"),			 SO_NONE    }, // "-opposed"
	
	{ OPT_FLAG, _T("--qmom"),			 SO_REQ_SEP    }, // "-premix"
	{ OPT_FLAG, _T("--2e"),				 SO_REQ_SEP    }, // "-premix"
    { OPT_ARG,  _T("--input"),           SO_REQ_SEP }, // "-input    ARG"
    { OPT_ARG,  _T("--output"),          SO_REQ_SEP }, // "-output   ARG"
    { OPT_ARG,  _T("--flame"),           SO_REQ_SEP }, // "-flame    ARG"
    { OPT_ARG,  _T("--kinetics"),        SO_REQ_SEP }, // "-kinetics ARG"
	{ OPT_ARG,  _T("--kinetics-ascii"),  SO_REQ_SEP }, // "-kinetics-ascii ARG"
    { OPT_ARG,  _T("--global"),          SO_REQ_SEP }, // "-global   ARG"
    { OPT_ARG,  _T("--schedule"),        SO_REQ_SEP }, // "-schedule ARG"
    { OPT_ARG,  _T("--backup"),          SO_REQ_SEP }, // "-backup   ARG"
    { OPT_ARG,  _T("--unsteady"),        SO_REQ_SEP }, // "-unsteady ARG"
    { OPT_ARG,  _T("--flame-speed-analysis"),   SO_REQ_SEP }, // "-unsteady ARG"
	{ OPT_ARG, _T("--opposed-flame-analysis-regular"), SO_REQ_SEP }, // "-unsteady ARG"
	{ OPT_ARG, _T("--opposed-flame-analysis-extinction"), SO_REQ_SEP }, // "-unsteady ARG"
	{ OPT_ARG, _T("--opposed-flame-analysis-ignition"), SO_REQ_SEP }, // "-unsteady ARG"

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

	OpenSMOKE_Flame1D									flame;
	OpenSMOKE_Flame1D_DataManager		data;
	OpenSMOKE_Flame1D_ScheduleClass		operations;
	
	OpenSMOKE_GlobalKinetics			global;
	OpenSMOKE_ReactingGas				mix;


	// Prepare
	flame.Assign(&mix);
	flame.Assign(&global);
	flame.Assign(&data);
	flame.Assign(&operations);
	data.Assign(&mix);
	data.Assign(&flame);

	// 0. Parser setup
    parser.setup(argc, argv, parser_options);

    // 1. Detailed kinetic scheme setup
    if (parser.parse("--kinetics", argument))
	{
		mix.SetupBinary(argument);
	}
    else if (parser.parse("--kinetics-ascii", argument))
	{
		mix.SetASCII();
		mix.SetupBinary(argument);
	}
    else ErrorMessage("The kinetic mechanism must be specified: --kinetics NAMEFOLDER || --kinetics-ascii NAMEFOLDER");

	// 2a. Global kinetic scheme setup
    global.assign_mix(&mix);
    if (parser.parse("--global", argument))
	{
		global.read_from_file(argument);
		flame.data->iGlobalKinetics = true;
	}

	// 2a. Kind of flame
	if (parser.parse("--qmom", argument))
	{
		data.qmom_file_name = argument;
		data.iQMOM = true;
	}
	else data.iQMOM = false;
	if (parser.parse("--2e", argument))
	{
		data.twoEquation_file_name = argument;
		data.i2E = true;
	}
	else data.i2E = false;

		 if (parser.parse("-premix")  && data.iQMOM == false && data.i2E == false) data.Setup("PREMIXED");
	else if (parser.parse("-opposed") && data.iQMOM == false && data.i2E == false) data.Setup("OPPOSED");
	else if (parser.parse("-twin")	  && data.iQMOM == false && data.i2E == false) data.Setup("TWIN");
	else if (parser.parse("-premix")  && data.iQMOM == true && data.i2E == false) data.Setup("PREMIXED_QMOM");
	else if (parser.parse("-opposed") && data.iQMOM == true && data.i2E == false) data.Setup("OPPOSED_QMOM");
	else if (parser.parse("-twin")    && data.iQMOM == true && data.i2E == false) data.Setup("TWIN_QMOM");
	else if (parser.parse("-premix")  && data.iQMOM == false && data.i2E == true) data.Setup("PREMIXED_SOOT");
	else if (parser.parse("-opposed") && data.iQMOM == false && data.i2E == true) data.Setup("OPPOSED_SOOT");
	else if (parser.parse("-twin")    && data.iQMOM == false && data.i2E == true) data.Setup("TWIN_SOOT");
	else    ErrorMessage("The kind of 1D Flame must be specified: -premix || -opposed || -twin || -premix-2e || -opposed-2e || -twin-2e || -premix-qmom || -opposed-qmom || -twin-qmom");

	// Output folder
    if (parser.parse("--output", argument))
	{
		data.SetOutputFolder(argument);
		if (parser.parse("--backup", argument))
			ErrorMessage("--output and --backup options cannot be used together...");
	}			

	// 4. Schedule Class Setup
	if (parser.parse("--schedule", argument))
            operations.readOperations(argument);
    else    operations.readOperations("Input/Operations.inp");

	// 5. Unsteady calculations
    if (parser.parse("--unsteady", argument))
	{
		data.iUnsteady = true;						// TODO
		data.unsteady_flame_file_name = argument;	// TODO
		flame.iUnsteadyFromBackUp = 0;				// TODO
	}
	else data.iUnsteady = false;

	// Backup
    if (parser.parse("--backup", argument))
	{
		data.iBackUp = true;
		flame.FoldersAndFilesManager(argument);
		flame.recoverFromBackUp(flame.nameFileBackupInputData);
	}
	else
	{	
		// 3a. Data Manager Setup
		if (parser.parse("--input", argument))
		{
			if (data.kind_of_flame == FLAME1D_PHYSICS_OPPOSED)	data.readFromFileForOpposed(argument);
			if (data.kind_of_flame == FLAME1D_PHYSICS_TWIN)		data.readFromFileForOpposed(argument);
			if (data.kind_of_flame == FLAME1D_PHYSICS_PREMIXED)	data.readFromFileForPremixed(argument);
		}
		else ErrorMessage("The input file must be specified: --input NAMEFILE");

		flame.setup();
		flame.FoldersAndFilesManager();
	}

	// 6. Flame speed analysis
    if (parser.parse("--flame-speed-analysis", argument))
	{
		data.iFlameSpeedAnalysis = true;			
		data.flameSpeedAnalysisFileName = argument;	
	}
	else data.iFlameSpeedAnalysis = false;	


	// 6. Flame speed analysis
    if (parser.parse("--opposed-flame-analysis-regular", argument))
	{
		data.iOpposedFlameAnalysis = true;	
		data.OpposedAnalysisType = "Regular";
		data.opposedFlameAnalysisFileName = argument;	
	}
	else if (parser.parse("--opposed-flame-analysis-ignition", argument))
	{
		data.iOpposedFlameAnalysis = true;
		data.OpposedAnalysisType = "Ignition";
		data.opposedFlameAnalysisFileName = argument;
	}
	else if (parser.parse("--opposed-flame-analysis-extinction", argument))
	{
		data.iOpposedFlameAnalysis = true;
		data.OpposedAnalysisType = "Extinction";
		data.opposedFlameAnalysisFileName = argument;
	}
	else data.iOpposedFlameAnalysis = false;	

	double tStart = BzzGetCpuTime();

		flame.Run();
	
	double tEnd = BzzGetCpuTime();

	// Finalize
	cout << "Total time for solution: " << tEnd - tStart << " s" << endl;
	OpenSMOKE_logo("OpenSMOKE_Flame1D", "0.5", "November 2015");

	return 1;
}

void ParserClass::ShowUsage()
{
    cout << "Usage: OpenSMOKE_Flame1D [-option] [--option FILE] [-?] [-help]" << endl;

    cout << "Options: "																					<< endl;
    cout << "    -?									this help"													<< endl;
    cout << "    -help								this help"													<< endl;
    cout << "    -noverbose							no file nor video output"									<< endl;
	cout << "    -premix							premixed flame"												<< endl;
	cout << "    -opposed							opposed flame"												<< endl;
	cout << "    -twin								twin flame"												<< endl;
	cout << "   --qmom FILENAME						qmom simulation"									<< endl;
	cout << "   --2e FILENAME						2e simulation"									<< endl;
	cout << "   --global FILENAME					global kinetics file name"									<< endl;
    cout << "   --input FILENAME					input file"	<< endl;
    cout << "   --kinetics FOLDER					detailed kinetic scheme folder (Binary version)"			<< endl;
	cout << "   --kinetics-ascii FOLDER				detailed kinetic scheme folder (ASCII version)"			<< endl;
	cout << "   --output FOLDER						output folder"				<< endl;
    cout << "   --schedule FILENAME					input file for operations"	<< endl;
    cout << "   --flame-speed-analysis FILENAME		input file for flame speed analysis"					<< endl;
	cout << "   --opposed-flame-analysis-regular FILENAME	 input file for opposed flame analysis" << endl;
	cout << "   --opposed-flame-analysis-extinction FILENAME input file for opposed flame analysis" << endl;
	cout << "   --opposed-flame-analysis-ignition FILENAME	 input file for opposed flame analysis" << endl;
	cout << endl;
}