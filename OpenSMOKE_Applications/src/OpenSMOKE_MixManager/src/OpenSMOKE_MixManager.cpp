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
#include "SimpleOpt.h"
#include "mixmanager/OpenSMOKE_MixManager.h"
#include "engine/OpenSMOKE_EquilibriumStanjan.h"
#include "basic/OpenSMOKE_Utilities.h"

void ErrorMessage(const string message)
{
    cout << endl;
    cout << "Executable:  OpenSMOKE_MixManager"		<< endl;
    cout << "Error:       "  << message				<< endl;
    cout << "Press a key to continue... "			<< endl;
    getchar();
    exit(-1);
}

// define the ID values to indentify the option
enum { OPT_HELP, OPT_FLAG, OPT_ARG };

CSimpleOpt::SOption parser_options[] =
{
    { OPT_FLAG, _T("--noverbose"),       SO_NONE    }, // "-verbose"
    { OPT_ARG,  _T("--input"),           SO_REQ_SEP }, // "-input			ARG"
    { OPT_ARG,  _T("--output"),          SO_REQ_SEP }, // "-output			ARG"
    { OPT_ARG,  _T("--kinetics"),        SO_REQ_SEP }, // "-kinetic			ARG"
    { OPT_ARG,  _T("--global"),          SO_REQ_SEP }, // "-global			ARG"
    { OPT_ARG,  _T("--inlet"),           SO_REQ_SEP }, // "-inlet			ARG"
	{ OPT_ARG,  _T("--soot2EModel"),     SO_REQ_SEP }, // "-soot2EModel		ARG"

    { OPT_HELP, _T("-?"),     SO_NONE    }, // "-?"
    { OPT_HELP, _T("--help"), SO_NONE    }, // "--help"
    SO_END_OF_OPTIONS                       // END
};

int main(int argc, char* argv[])
{
	cout.setf(ios::scientific);

    string						argument;
    ParserClass					parser;

	OpenSMOKE_ReactingGas		mix;
    OpenSMOKE_GlobalKinetics    global;
    OpenSMOKE_MixManager		mixture;
	OpenSMOKE_GasStream			inlet;
	OpenSMOKE_GasStream			outlet;
	OpenSMOKE_2EModel			soot2EModel;

    // 0. Parser setup
    parser.setup(argc, argv, parser_options);


    // 1. Detailed kinetic scheme setup
    if (parser.parse("--kinetics", argument))
            mix.SetupBinary(argument);
    else    ErrorMessage("The --kinetics option is compulsory");


    // 2a. Global kinetic scheme setup
    global.assign_mix(&mix);
    if (parser.parse("--global", argument))
    {
		global.read_from_file(argument);
		mixture.SetGlobalKinetics();
	}

	// 2b. Soot 2E Model
	if (parser.parse("--soot2EModel", argument))
	{
		soot2EModel.assign_mixture(mix);
		soot2EModel.setupFromFile(argument);
	}

    // 3. Inlet stream setup
    inlet.AssignKineticScheme(mix);
    if (parser.parse("--inlet", argument))
		inlet.DefineFromFile(argument);
    else    ErrorMessage("The --inlet option is compulsory");


    // 4. PFR Setup
    mixture.AssignKineticScheme(mix);
    mixture.AssignGlobalKineticScheme(global);
	mixture.AssignSoot2EModel(soot2EModel);
    mixture.AssignInletFlows(inlet);
//    if (parser.parse("--input", argument))
//		mixture.DefineFromFile(argument);
//    else    ErrorMessage("The --input option is compulsory");


	mixture.Initialize();
	mixture.Properties();

    // 5. PFR Video Summary
    mixture.VideoSummary();

	// 8. Outlet Stream
//	mixture.OutletStream(outlet);
//	outlet.VideoSummary();

/*PROVV */
/*
	vector<string> fuel_names;
	vector<string> oxidizer_names;
	fuel_names.push_back("fuels");
	fuel_names.push_back("CH4");
	fuel_names.push_back("H2");
	fuel_names.push_back("CO");
	fuel_names.push_back("N2");
	fuel_names.push_back("O2");
	oxidizer_names.push_back("oxidizer");
	oxidizer_names.push_back("O2");
	oxidizer_names.push_back("N2");
	oxidizer_names.push_back("AR");
	BzzVector moles_fuel(5, 1., 1., 0.5, 0, 0.25);
	BzzVector moles_oxidizer(3, 1., 2., 3.);
	inlet.AssignEquivalenceRatio(1.0,fuel_names, moles_fuel,	oxidizer_names, moles_oxidizer);
	getchar();

*/
	BzzVector xInitial = inlet.x;
	BzzVector xFinal(mix.NumberOfSpecies());


/*	ofstream fOut;
	openOutputFileAndControl(fOut, "Cp.txt");
	BzzVector wVector(mix.NumberOfSpecies());
	wVector = 1./double(mix.NumberOfSpecies());
	for(int i=1;i<=1400;i++)
	{
		double Tlocal = 373+double(i);
		mix.ComputeKineticParameters( Tlocal, log(Tlocal), 1./Tlocal);

		mix.SpeciesCp(373+double(i));
		mix.MixCp_FromMassFractions(wVector);
		fOut << 373+double(i) << "\t";
		for(int j=1;j<=mix.NumberOfSpecies();j++)
			fOut << mix.Cp[j]/1000.*mix.M[j] << "\t";

		for(j=1;j<=mix.NumberOfReactions();j++)
			fOut << mix.reactionDH[j]*Constants::R_J_kmol*Tlocal/1000. << "\t";

		fOut << endl;

	}
	fOut.close();
*/
/*	{
		ifstream fIn;
		openInputFileAndControl(fIn, "ListOfFiles.txt");
		ofstream fOut;
		openOutputFileAndControl(fOut, "ListOfFiles.out");

		int i;
		string dummy;
		vector<string> list;
		for(i=1;i<=342;i++)
		{
			fIn >> dummy;
			list.push_back(dummy);
		}
		for(i=1;i<=342;i++)
		{
			fOut << "\t\t$(PATH_OBJECTS)/" << list[i-1] << ".o \\" << endl;
		}
		fOut << endl << endl;
		for(i=1;i<=342;i++)
		{
			fOut << "$(PATH_OBJECTS)/" << list[i-1] << ".o : $(PATH_CPP)/" << list[i-1]<< ".cpp" << endl;
			fOut << "\t\t$(CCP) $(PATH_CPP)/" << list[i-1] << ".cpp -o $(PATH_OBJECTS)/" << list[i-1] << ".o" << endl;
			fOut << endl;
		}

		fIn.close();
		fOut.close();
	}
	*/
/*	
	double tStart = BzzGetCpuTime();
	for(int k=1;k<=1;k++)
	{	
		OpenSMOKE_EquilibriumStanjan equilibrium;
		equilibrium.Setup(&mix);
		equilibrium.SetElementalCompositionFromSpeciesMoleFractions(xInitial);
		equilibrium.Equilibrate(inlet.T, inlet.P, xFinal);
	}
	double tEnd = BzzGetCpuTime();
	cout << endl << "Equilibrium time: " << tEnd - tStart << endl << endl;	
*/

//	mix.FlameTemperature(inlet.T, inlet.P, xInitial);

/*	OpenSMOKE_EquilibriumStanjan equilibrium;
		equilibrium.Setup(&mix);
		equilibrium.SetElementalCompositionFromSpeciesMoleFractions(xInitial);
		equilibrium.Equilibrate(inlet.T, inlet.P, xFinal);
*/
    OpenSMOKE_logo("OpenSMOKE_MixManager", "0.12", "October 2009");
	getchar();
    return 0;
}

void ParserClass::ShowUsage()
{
    cout << "Usage: OpenSMOKE_MixManager [-option] [--option STRING] [-?] [--help]" << endl;

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
