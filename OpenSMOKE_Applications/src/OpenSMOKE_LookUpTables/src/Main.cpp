/***************************************************************************
 *   Copyright (C) 2006-2010 by Alberto Cuoci							   *
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
#include "flamelet_group.h"
#include "basic/OpenSMOKE_Utilities.h"
#include "interfaces/SimpleOpt.h"
#include "engine/OpenSMOKE_ReactingGas.h"
#include "OpenSMOKE_LookUp_Table_Manager.h"
#include "OpenSMOKE_LookUp_Table_Flame.h"
#include "OpenSMOKE_LookUp_Table_Executables.h"


// define the ID values to indentify the option
enum { OPT_HELP, OPT_FLAG, OPT_ARG };

void ErrorMessage(const string message)
{
    cout << endl;
    cout << "Executable:  OpenSMOKE_LookUpTable"	<< endl;
    cout << "Error:       "  << message			<< endl;
    cout << "Press a key to continue... "		<< endl;
    getchar();
    exit(-1);
}

CSimpleOpt::SOption parser_options[] =
{
    { OPT_HELP, _T("-?"),                SO_NONE    }, // "-?"
    { OPT_HELP, _T("-help"),             SO_NONE    }, // "-help"

	{ OPT_ARG,  _T("-build-library"),								SO_NONE    }, // "-build-library"
	
	{ OPT_ARG,  _T("-fluent-lut-exe"),								SO_NONE    }, // "-fluent-lut-exe"
	{ OPT_ARG,  _T("-fluent-lut-uncorrelated-adiabatic-exe"),		SO_NONE    }, // "-fluent-lut-uncorrelated-adiabatic-exe"
	{ OPT_ARG,  _T("-fluent-lut-uncorrelated-non-adiabatic-exe"),	SO_NONE    }, // "-fluent-lut-uncorrelated-non-adiabatic-exe"
	{ OPT_ARG,  _T("-fluent-lut-correlated-adiabatic-exe"),			SO_NONE    }, // "-fluent-lut-correlated-adiabatic-exe"
	{ OPT_ARG,  _T("-fluent-lut-correlated-non-adiabatic-exe"),		SO_NONE    }, // "-fluent-lut-correlated-non-adiabatic-exe"

	{ OPT_ARG,  _T("--kinetics"),									SO_REQ_SEP },	// "--kinetics   ARG"
	{ OPT_ARG,  _T("--input"),										SO_REQ_SEP },	// "--input		ARG"
	{ OPT_ARG,  _T("-verbose-formation-rates"),					SO_NONE },	// "--input		ARG"
    SO_END_OF_OPTIONS																// END
};

// Executables
void FLUENT_LookUpTable();


int main(int argc, char* argv[])
{
	string          argument;

	OpenSMOKE_logo("OpenSMOKE_LookUpTable", "0.1", "February 2010");

	//	0. Parser setup
	ParserClass     parser;
	parser.setup(argc, argv, parser_options);

	// Build library
	if (parser.parse("-build-library"))
	{
		string kineticsFolder;
		string inputFile;

		if (parser.parse("--kinetics", argument))	kineticsFolder = argument.c_str();
		else										ErrorMessage("The --kinetics option is compulsory");

		if (parser.parse("--input", argument))		inputFile = argument.c_str();
		else										ErrorMessage("The --input option is compulsory");

		OpenSMOKE_ReactingGas mix;
		mix.SetupBinary(kineticsFolder);
		
		OpenSMOKE_LookUp_Table_Manager manager;

		if (parser.parse("-verbose-formation-rates"))
				manager.SetVerboseFormationRates(true);

		manager.SetKineticScheme(mix);
		manager.DefineFromFile(inputFile);
		

		manager.Run();

		exit(-1);
	}

	else if (parser.parse("-fluent-lut-exe"))
	{
		FLUENT_LookUpTable();
		exit(1);
	}

	else if (parser.parse("-fluent-lut-uncorrelated-adiabatic-exe"))
	{
		OpenSMOKE_LookUp_Table_Executables exe;
		exe.SetSootClosureModel("UNCORRELATED");
		exe.SetSootSourcesModel("ADIABATIC");
		exe.Run();
	}

	else if (parser.parse("-fluent-lut-uncorrelated-non-adiabatic-exe"))
	{
		OpenSMOKE_LookUp_Table_Executables exe;
		exe.SetSootClosureModel("UNCORRELATED");
		exe.SetSootSourcesModel("NON_ADIABATIC");
		exe.Run();
	}

	else if (parser.parse("-fluent-lut-correlated-adiabatic-exe"))
	{
		OpenSMOKE_LookUp_Table_Executables exe;
		exe.SetSootClosureModel("CORRELATED");
		exe.SetSootSourcesModel("ADIABATIC");
		exe.Run();
	}

	else if (parser.parse("-fluent-lut-correlated-non-adiabatic-exe"))
	{
		OpenSMOKE_LookUp_Table_Executables exe;
		exe.SetSootClosureModel("CORRELATED");
		exe.SetSootSourcesModel("NON_ADIABATIC");
		exe.Run();
	}

	return 0;
}

void ParserClass::ShowUsage()
{
    cout << "Usage: OpenSMOKE_LookUpTable [-option] [--option FILE] [-?] [--help]" << endl;

    cout << "Options: "                                                                          << endl;
    cout << "   -?                              this help"                                                << endl;
    cout << "   --kinetics FOLDER               folder (kinetic scheme)"                                  << endl;
    cout << "   --input FILENAME                input file name"                                          << endl;
	cout << "   -help                           this help"                                                << endl;
	cout << "   -build-library                  build pdf library"                                        << endl;
	cout << "   -fluent-lut-exe                 build fluent look up table executable"                    << endl;
	cout << "   -fluent-lut-uncorrelated-exe    build fluent look up table executable: soot sources"      << endl;
	cout << "   -fluent-lut-correlated-exe      build fluent look up table executable: soot sources"      << endl;
    cout << endl;
}


void choose_flames(BzzVector &list_of_chi, const double chi, int &aChi, int &bChi)	
{
	if (chi <= list_of_chi[1])
	{
		aChi = 1;					bChi = aChi;
	}

	else if (chi >= list_of_chi[list_of_chi.Size()])
	{
		aChi = list_of_chi.Size();	bChi = aChi;
	}

	else
		for(int j=2;j<=list_of_chi.Size();j++)
			if (chi <= list_of_chi[j])
			{
				aChi = j-1;		bChi = j;
				break;
			}
}


void extract_square_coordinates(const double xcsi, const double ycsiV, BzzVector &csi, BzzVector &csiV, int &aCsi, int &bCsi, int &aCsiV, int &bCsiV)
{
	int j;

	for(j=2;j<=csi.Size();j++)
		if (xcsi <= csi[j])
		{
			aCsi = j-1;
			bCsi = j;
			break;
		}

	for(j=2;j<=csiV.Size();j++)
		if (ycsiV <= csiV[j])
		{
			aCsiV = j-1;
			bCsiV = j;
			break;
		}
}

double extract_square_value(const double xcsi, const double ycsiV, BzzMatrix &y, BzzVector &csi, BzzVector &csiV, const int aCsi, const int bCsi, const int aCsiV, const int bCsiV)
{
	double a = y[aCsi][aCsiV] + (y[bCsi][aCsiV]-y[aCsi][aCsiV])/(csi[bCsi]-csi[aCsi])*(xcsi-csi[aCsi]);
	double b = y[aCsi][bCsiV] + (y[bCsi][bCsiV]-y[aCsi][bCsiV])/(csi[bCsi]-csi[aCsi])*(xcsi-csi[aCsi]);

	return  a + (b-a)/(csiV[bCsiV]-csiV[aCsiV])*(ycsiV-csiV[aCsiV]);
}


void extract_values(const double csi, const double csiV, OpenSMOKE_LookUp_Table_Flame &flame, 
					double &temperature, double &density, double &cp, double &as)
{
	int aCsi,	bCsi;
	int aCsiV,	bCsiV;

	extract_square_coordinates(csi, csiV, flame.csi, flame.csiV, aCsi, bCsi, aCsiV, bCsiV);
	
	temperature = extract_square_value(csi, csiV, flame.temperature,	flame.csi, flame.csiV, aCsi, bCsi, aCsiV, bCsiV);
	density		= extract_square_value(csi, csiV, flame.density,		flame.csi, flame.csiV, aCsi, bCsi, aCsiV, bCsiV);
	cp			= extract_square_value(csi, csiV, flame.cp,				flame.csi, flame.csiV, aCsi, bCsi, aCsiV, bCsiV);
	as			= extract_square_value(csi, csiV, flame.as,				flame.csi, flame.csiV, aCsi, bCsi, aCsiV, bCsiV);
}

void FLUENT_LookUpTable()
{
	int nc;
	int nChi;
	int NCells;
	BzzVector list_of_chi;
	BzzVector csi;
	BzzVector csiV;
	BzzVector chi;
	BzzVector temperature;
	BzzVector density;
	BzzVector cp;
	BzzVector as;

	// Read head file
	{
		ifstream fLUT;
		openInputFileAndControl(fLUT, "LUT/LookUpTable.out");

		string title;
		string dummy;

		fLUT >> title;
		fLUT >> dummy;
		fLUT >> nc;
		fLUT >> dummy;
		fLUT >> nChi;

		ChangeDimensions(nChi, &list_of_chi);
		for(int k=1;k<=nChi;k++)
			fLUT >> list_of_chi[k];
		fLUT.close();
	}

	// Read look up tables
	OpenSMOKE_LookUp_Table_Flame *flames;
	flames = new OpenSMOKE_LookUp_Table_Flame[nChi+1];
	for(int k=1;k<=nChi;k++)
	{
		stringstream kstring;
		kstring << k;
		flames[k].read_from_file("LUT/SR_" + kstring.str() + ".out");
	}

	// Read file from fluent
	{
		ifstream fFLUENT;
		openInputFileAndControl(fFLUENT, "FromFLUENT.out");

		fFLUENT >> NCells;
		ChangeDimensions(NCells, &csi);
		ChangeDimensions(NCells, &csiV);
		ChangeDimensions(NCells, &chi);

		for(int k=1;k<=NCells;k++)
		{
			fFLUENT >> csi[k];
			fFLUENT >> csiV[k];
			fFLUENT >> chi[k];
		}
		fFLUENT.close();
	}

	// Extract data
	{
		int aChi, bChi;
		double temperature_a,	temperature_b;
		double density_a,		density_b;
		double cp_a,			cp_b;
		double as_a,			as_b;

		ChangeDimensions(NCells, &temperature);
		ChangeDimensions(NCells, &density);
		ChangeDimensions(NCells, &cp);
		ChangeDimensions(NCells, &as);
		for(int k=1;k<=NCells;k++)
		{		
			choose_flames(list_of_chi, chi[k], aChi, bChi);
			if (aChi != bChi)
			{
				extract_values(csi[k], csiV[k], flames[aChi], temperature_a, density_a, cp_a, as_a);
				extract_values(csi[k], csiV[k], flames[bChi], temperature_b, density_b, cp_b, as_b);
		
				double coefficient = (chi[k]-list_of_chi[aChi])/(list_of_chi[bChi]-list_of_chi[aChi]);
				
				temperature[k]	= temperature_a + (temperature_b-temperature_a)*coefficient;
				density[k]		= density_a + (density_b-density_a)*coefficient;
				cp[k]			= cp_a + (cp_b-cp_a)*coefficient;
				as[k]			= as_a + (as_b-as_a)*coefficient;
			}
			else
				extract_values(csi[k], csiV[k], flames[aChi], temperature[k], density[k], cp[k], as[k]);
		}
	}

	// Extract files to FLUENT
	{
		ofstream fFLUENT;
		openOutputFileAndControl(fFLUENT, "ToFLUENT.out");
		fFLUENT.setf(ios::scientific);

		fFLUENT << NCells << endl;

		for(int k=1;k<=NCells;k++)
		{
			fFLUENT << temperature[k]	<< "\t";
			fFLUENT << density[k]		<< "\t";
			fFLUENT << cp[k]			<< "\t";
			fFLUENT << as[k]			<< "\t";
			fFLUENT << endl;
		}
		fFLUENT.close();
	}
}

