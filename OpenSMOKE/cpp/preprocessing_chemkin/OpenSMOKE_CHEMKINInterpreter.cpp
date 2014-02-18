/***************************************************************************
 *   Copyright (C) 2009 by Alberto Cuoci								   *
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

#include "basic/OpenSMOKE_Utilities.h"
#include "preprocessing_chemkin/OpenSMOKE_CHEMKINInterpreter.h"
#include "preprocessing/OpenSMOKE_PreProcessorIdealGas.h"
#include "engine/OpenSMOKE_ReactingGas.h"
#include <sstream>

void OpenSMOKE_CHEMKINInterpreter::ErrorMessage(const int iLine, const string message)
{
    cout << endl;
    cout << "Class:  OpenSMOKE_CHEMKINInterpreter "	<< endl;
    cout << "Object: " << name_object		<< endl;
	cout << "Line:   " << iLine				<< endl;
    cout << "Error:  " << message           << endl;
    cout << "Press a key to continue... "   << endl;
    getchar();
    exit(-1);
}

void OpenSMOKE_CHEMKINInterpreter::WarningMessage(const int iLine, const string message)
{
    cout << endl;
    cout << "Class:    OpenSMOKE_CHEMKINInterpreter"	<< endl;
    cout << "Object:   " << name_object	<< endl;
	cout << "Line:     " << iLine				<< endl;
    cout << "Warning:  " << message      << endl;
    cout << "Press a key to continue... "   << endl;
    getchar();
}

OpenSMOKE_CHEMKINInterpreter::OpenSMOKE_CHEMKINInterpreter()
{
	TMIN = 300.;
	TMAX = 3600.;
	FITTINGPOINTS = 50;

	name_kinetics_file	       = "[null]";
	name_surface_kinetics_file = "[null]";
	name_thermo_file	       = "[null]";
	name_transport_file        = "[null]";
	name_reduced_list_file     = "[null]";

	name_object   = "[not assigned]";
	name_author   = "[not assigned]";
	building_date = GiveMeTimeAndDate();

	iSurfaceKinetics = false;
	iBuzziMode		 = false;
	iVerboseMode	 = false;
	iSootMode		 = false;

	startElements	= 0;
	startSpecies	= 0;
	startReactions	= 0;
}

string OpenSMOKE_CHEMKINInterpreter::nameKineticsFile()
{
	return name_kinetics_file;
}

string OpenSMOKE_CHEMKINInterpreter::nameSurfaceKineticsFile()
{
	return name_surface_kinetics_file;
}

string OpenSMOKE_CHEMKINInterpreter::nameThermoFile()
{
	return name_thermo_file;
}

string OpenSMOKE_CHEMKINInterpreter::nameTransportFile()
{
	return name_transport_file;
}

string OpenSMOKE_CHEMKINInterpreter::nameOutputFolder()
{
	return name_output_folder;
}



void OpenSMOKE_CHEMKINInterpreter::SetMinimumTemperature(const double tmin)
{
	TMIN		= tmin;
	if (TMIN >= TMAX || TMIN <= 0.)
		ErrorMessage(0, "Please check your minimum temperature!");
}

void OpenSMOKE_CHEMKINInterpreter::SetMaximumTemperature(const double tmax)
{
	TMAX		= tmax;
	if (TMIN >= TMAX || TMAX >= 6000.)
		ErrorMessage(0, "Please check your maximum temperature!");
}

void OpenSMOKE_CHEMKINInterpreter::SetFittingPoints(const int fittingpoints)
{
	FITTINGPOINTS	= fittingpoints;
	if (FITTINGPOINTS > 100 || FITTINGPOINTS <= 4)
		ErrorMessage(0, "Please check your number of fitting points!");
}

void OpenSMOKE_CHEMKINInterpreter::SetName(const string name)
{
	name_object = name;
}

void OpenSMOKE_CHEMKINInterpreter::SetAuthorName(const string name)
{
	name_author = name;
}

void OpenSMOKE_CHEMKINInterpreter::SetKineticFileName(const string name)
{
	name_kinetics_file = name;
}

void OpenSMOKE_CHEMKINInterpreter::SetSurfaceKineticFileName(const string name)
{
	iSurfaceKinetics = true;
	name_surface_kinetics_file = name;
}

void OpenSMOKE_CHEMKINInterpreter::SetThermodynamicFileName(const string name)
{
	name_thermo_file = name;
}

void OpenSMOKE_CHEMKINInterpreter::SetTransportFileName(const string name)
{
	name_transport_file = name;
}

void OpenSMOKE_CHEMKINInterpreter::SetReducedListFileName(const string name)
{
	name_reduced_list_file = name;
}

void OpenSMOKE_CHEMKINInterpreter::SetOutputFolderName(const string name)
{
	name_output_folder = name;

	string MSDOScommand = "mkdir " + name_output_folder;
	system(MSDOScommand.c_str());

	name_ascii_folder = name_output_folder;
	#if LINUX_SO==1
		name_ascii_folder     += "/ascii";
	#else
		name_ascii_folder     += "\\ascii";
	#endif
}

void OpenSMOKE_CHEMKINInterpreter::SetBuzziMode()
{
	iBuzziMode = true;
}

void OpenSMOKE_CHEMKINInterpreter::SetVerboseMode()
{
	iVerboseMode = true;
}

void OpenSMOKE_CHEMKINInterpreter::SetNoVerboseMode()
{
	iVerboseMode = false;
}

void OpenSMOKE_CHEMKINInterpreter::SetSootMode()
{
	iSootMode = true;
}

void OpenSMOKE_CHEMKINInterpreter::Run()
{
	double timeStart = BzzGetCpuTime();

	fWarning.Setup("Warning.log");
	openOutputFileAndControl(fLog, "Log.log");

	kinetics.Setup(&fLog, &fWarning);
	species.Setup(&fLog, &fWarning);
	units.Setup(&fLog, &fWarning);
	elements.Setup(&fLog, &fWarning);
	
	cout << " ** Reading thermodynamic database..." << endl;
	thermo.ReadThermoData(name_thermo_file, &fLog);

	if (name_transport_file != "[null]")
	{
		cout << " ** Reading transport properties database..." << endl;
		transport.ReadTransportData(name_transport_file, &fLog);
	}

	cout << " ** Reading kinetic scheme..." << endl;
	ReadKineticScheme();

	// ---------------------------------------------------------------
	// 1. Parsing sections
	// ---------------------------------------------------------------
	cout << "    a. Parsing sections..." << endl;
	ParsingSections();

	// ---------------------------------------------------------------
	// 2. Syntax corrections
	// ---------------------------------------------------------------
	cout << "    b. Syntax corrections..." << endl;	
	SyntaxCorrections();

	// ---------------------------------------------------------------
	// 3. Parsing units
	// ---------------------------------------------------------------
	cout << "    c. Parsing units..." << endl;	
	ParsingUnits();

	// ---------------------------------------------------------------
	// 4. Cleaning lines
	// ---------------------------------------------------------------
	cout << "    d. Cleaning lines..." << endl;	
	CleaningLines();

	// ---------------------------------------------------------------
	// 5. Parsing elements
	// ---------------------------------------------------------------
	cout << "    e. Parsing elements..." << endl;	
	ParsingElements();

	// ---------------------------------------------------------------
	// 6. Parsing species
	// ---------------------------------------------------------------
	cout << "    f. Parsing species..." << endl;	
	ParsingSpecies();

	// ---------------------------------------------------------------
	// 7. Checking thermodynamic data
	// ---------------------------------------------------------------
	cout << "    g. Checking thermodynamic data..." << endl;	
	indices_species_thermo_database = thermo.GiveMeSpeciesIndices(kinetics.species_list);
	thermo.GiveMeMolecularWeights(indices_species_thermo_database, elements);

	// ---------------------------------------------------------------
	// 8. Checking transport properties
	// ---------------------------------------------------------------
	if (transport.IsActivated() == true)
	{
		cout << "    h. Checking transport properties..." << endl;	
		indices_species_transport_database = transport.GiveMeSpeciesIndices(kinetics.species_list);
	}

	// ---------------------------------------------------------------
	// 9. Parsing reactions
	// ---------------------------------------------------------------
	cout << "    i. Parsing reactions..." << endl;	
	ParsingReactions();

	// ---------------------------------------------------------------
	// 10. Parsing reactions (additional data)
	// ---------------------------------------------------------------
	cout << "    j. Parsing reactions (additional data)..." << endl;	
	ParsingReactionsAdditionalData();

	if (name_reduced_list_file != "[null]")
	{
		cout << "    k  Reducing the scheme..." << endl;
		WriteReducedScheme(name_reduced_list_file);
		cout << "Press enter to exit..." << endl;
		getchar();
		exit(-1);
	}

	// ---------------------------------------------------------------
	// 11. Checking thermodinamyc data
	// ---------------------------------------------------------------
	cout << "    k. Final Checks..." << endl << endl;	
	kinetics.FinalChecks();

	// ---------------------------------------------------------------
	// 12. Reading Transport Properties
	// ---------------------------------------------------------------
	if (transport.IsActivated() == true)
	{
		cout << " ** Preparing transport properties..." << endl;
		transport.ProcessTransportData(indices_species_transport_database);
	}

	// ---------------------------------------------------------------
	// 13. Preparing Thermodynamics
	// ---------------------------------------------------------------
	cout << " ** Preparing elements..." << endl;
	thermo.ProcessElementData(indices_species_thermo_database, elements);

	// ---------------------------------------------------------------
	// 14. Preparing Thermodynamics
	// ---------------------------------------------------------------
	cout << " ** Preparing thermodynamics..." << endl;
	thermo.ProcessThermoData(indices_species_thermo_database);

	// ---------------------------------------------------------------
	// 15. Preparing Stoichiometry
	// ---------------------------------------------------------------
	cout << " ** Preparing stoichiometry..." << endl;
	kinetics.PrepareStoichiometry();

	// ---------------------------------------------------------------
	// 16. Statistics
	// ---------------------------------------------------------------
	cout << " ** Statistics..." << endl;
	kinetics.Statistics();

	// ---------------------------------------------------------------
	// 19. Writing Binary Files
	// ---------------------------------------------------------------
	cout << " ** Writing Ideal Gas Binary Files..." << endl;
	PrintIdealGasBinaryFile(name_output_folder + "/idealgas.bin");

	// ---------------------------------------------------------------
	// 20. Writing Binary Files
	// ---------------------------------------------------------------
	cout << " ** Writing Kinetics Binary Files (gas)..." << endl;
	kinetics.PrintBinaryFile(name_output_folder + "/reactions.bin", thermo, transport, *this);

	// ---------------------------------------------------------------
	// 17. Writing Binary Files
	// ---------------------------------------------------------------
	if (iBuzziMode == true)
	{
		cout << " ** Writing ASCII Files (Buzzi Mode)..." << endl;
		PrintBuzziFiles();
	}
	double timeEnd = BzzGetCpuTime();
		
	cout << " ** Kinetic Scheme succesfully interpreted!" << endl;
	cout << " ** Total time: " << timeEnd - timeStart << "s";
	cout << endl;

	// Surface
	if (iSurfaceKinetics == true)
	{
		cout << "    Surface kinetics..." << endl;	
		ReadSurfaceKineticScheme();

		cout << "    b. Syntax corrections..." << endl;	
		SurfaceSyntaxCorrections();

		cout << "    a. Parsing sections..." << endl;	
		SurfaceParsingSections();
	
		//cout << "    c. Parsing units..." << endl;	
		//SurfaceParsingUnits();

		cout << "    d. Cleaning lines..." << endl;	
		SurfaceCleaningLines();

		cout << "    e. Parsing materials..." << endl;	
		SurfaceParsingMaterials();

		cout << "    e. Parsing sites..." << endl;	
		SurfaceParsingSites();

		cout << "    f. Parsing bulks..." << endl;	
		SurfaceParsingBulk();

		cout << "    g. Checking thermodynamic data..." << endl;	
		indices_site_thermo_database = thermo.GiveMeSpeciesIndices(site_species_global_list);
		indices_bulk_thermo_database = thermo.GiveMeSpeciesIndices(bulk_species_global_list);
		thermo.GiveMeMolecularWeights(indices_site_thermo_database, elements);
		thermo.GiveMeMolecularWeights(indices_bulk_thermo_database, elements);

		cout << " ** Preparing thermodynamics..." << endl;
		thermo.ProcessThermoData(indices_site_thermo_database, indices_bulk_thermo_database);

		cout << "    g. Surface Parsing reactions..." << endl;	
		SurfaceParsingReactions();

		cout << "    h. Surface Parsing reactions (additional data)..." << endl;	
		SurfaceParsingReactionsAdditionalData();

		cout << "    k. Surface Final Checks..." << endl << endl;	
		for (int k=1;k<=nSurfaceMaterials;k++)
			material[k].surfaceKinetics.FinalChecks();

		cout << " ** Preparing elements..." << endl;
		thermo.ProcessElementData(indices_site_thermo_database, indices_bulk_thermo_database, elements);

		cout << " ** Statistics..." << endl;
		for (int k=1;k<=nSurfaceMaterials;k++)
			material[k].surfaceKinetics.Statistics();

		cout << " ** Writing Kinetics Binary Files (surface +++)..." << endl;
		{
			PrintSurfaceBinaryFile(name_output_folder + "/surface.bin");
		}

		cout << " ** Writing Summary on Files..." << endl;
		SurfaceSummaryOnFile();
	}

	// Print Additional data
	if (iVerboseMode == true)
	{
		cout << " ** Writing GasPhaseSummary.out file..." << endl;
		PrintSummaryFile(name_output_folder + "/GasPhaseSummary.out");
	}

	if (iVerboseMode == true)
	{
		WriteAdditionalFiles();

		OpenSMOKE_ReactingGas mix;
		mix.SetupBinary(name_output_folder);

		cout << " ** Writing Analysis.out file..." << endl;
		ofstream fAnalyzer;
		openOutputFileAndControl(fAnalyzer, name_output_folder + "/Analysis.out");
		ofstream fOutputInverseKinetics;
		openOutputFileAndControl(fOutputInverseKinetics, name_output_folder + "/InverseKinetics.out");
		mix.VerboseDataSpecies(fAnalyzer);
		mix.VerboseDataKinetics(fAnalyzer, fOutputInverseKinetics);
		mix.VerboseDataSpeciesCheckConsistency(fAnalyzer);
		fAnalyzer.close();

		cout << " ** Writing Elements.out file..." << endl;
		mix.WriteElementsFile(name_output_folder + "/Elements.out");		
	}

}

void OpenSMOKE_CHEMKINInterpreter::ParsingElements(const int iLine, vector<string> instructions)
{
	int number_of_instructions = instructions.size()-1;
	BzzVectorInt indexSlash;
	BzzVectorInt isotope_position;
	BzzVectorInt element_position;
	BzzVector	 isotope_mw;

	int i;
	for(i=1;i<=number_of_instructions;i++)
	{
		if ( caseInsCompare(instructions[i],"END") && i!= number_of_instructions)
			if (caseInsCompare(instructions[i+1],"ELEMENTS")==false && caseInsCompare(instructions[i+1],"ELEM")==false)
				ErrorMessage(iLine, "Found END");
		
		if (instructions[i] == "/")	indexSlash.Append(i);
	}

	if (indexSlash.Size()%2 != 0)
		ErrorMessage(iLine, "Wrong number of slashes");

	for(i=1;i<=indexSlash.Size()/2;i++)
	{
		string parsed_string;
		int iBegin	= indexSlash[(i-1)*2+1];
		int iEnd	= indexSlash[(i-1)*2+2];
		if (iEnd != iBegin+2)
			ErrorMessage(iLine, "Error in slashes");

		double mw = atof(instructions[iBegin+1].c_str());
		isotope_position.Append(iBegin-1);
		isotope_mw.Append(mw);
	}

	i=0;
	for(;;)
	{	
		i++;
		if (caseInsCompare(instructions[i],"END") || caseInsCompare(instructions[i],"ELEM") || caseInsCompare(instructions[i],"ELEMENTS")) 
		{
			if (i>=number_of_instructions) break;
			continue;
		}

		bool iIsotope = false;
		for(int j=1;j<=isotope_position.Size();j++)
			if (isotope_position[j] == i)
			{
				i += 3;
				iIsotope = true;
				break;
			}

		if (iIsotope == false)
		{
			element_position.Append(i);
		}
		
		if (i>=number_of_instructions) break;
	}


	for(i=1;i<=isotope_position.Size();i++)
		elements.Parse_Isotope_Name(instructions[isotope_position[i]], isotope_mw[i]);

	for(i=1;i<=element_position.Size();i++)
	{
		bool success = elements.Parse_Element_Name(instructions[element_position[i]]);
		if (success == false)
			ErrorMessage(iLine, "The " + instructions[element_position[i]] + " element is not included in the database");
	}
}

void OpenSMOKE_CHEMKINInterpreter::ParsingSpecies(const int iLine, vector<string> instructions)
{
	int number_of_instructions = instructions.size()-1;

	int i;
	for(i=1;i<=number_of_instructions;i++)
	{
		if ( caseInsCompare(instructions[i],"END") && i!= number_of_instructions)
			if (caseInsCompare(instructions[i+1],"SPECIES")==false && caseInsCompare(instructions[i+1],"SPEC")==false)
				ErrorMessage(iLine, "Found END");
		
		if (instructions[i] == "/")	ErrorMessage(iLine, "Wrong number of slashes");
	}

	i=0;
	for(;;)
	{	
		i++;
		if (caseInsCompare(instructions[i],"END") || caseInsCompare(instructions[i],"SPEC") || caseInsCompare(instructions[i],"SPECIES")) 
		{
			if (i>=number_of_instructions) break;
			continue;
		}

		bool success;
		success = species.Parse_Species_Name(instructions[i]);

		if (i>=number_of_instructions) break;
	}

	if (iSootMode == true)	species.Parse_Species_Name("CSOOT");
}


void OpenSMOKE_CHEMKINInterpreter::ParsingReactions(const int iReaction, const int iLine, vector<string> instructions)
{
	int number_of_instructions = instructions.size()-1;
	{
		kinetics.reaction[iReaction].A		= atof(instructions[number_of_instructions-2].c_str());
		kinetics.reaction[iReaction].Beta	= atof(instructions[number_of_instructions-1].c_str());
		kinetics.reaction[iReaction].E		= atof(instructions[number_of_instructions].c_str());
	}

	int i;
	string whole_reaction;
	for(i=1;i<=number_of_instructions-3;i++)
		whole_reaction+= (instructions[i] + " ");

	AddSpaces(whole_reaction, '+');

	vector<string> stoichiometry;
	stoichiometry.push_back("stoichiometry");
		
	string dummy;
	string dummy_complete;
	stringstream parsed_string(whole_reaction);

	int index;
	kinetics.reaction[iReaction].thirdbodysingle = "null";
	for(;;)
	{
		parsed_string >> dummy;
		if (parsed_string.fail())
			break;
		stoichiometry.push_back(dummy);
		
		if (dummy == "#REVERSIBLE#")
		{
			kinetics.reaction[iReaction].iReversible = true;
			index = stoichiometry.size()-1;
		}	
		if (dummy == "#DIRECT#")
		{
			kinetics.reaction[iReaction].iReversible = false;
			index = stoichiometry.size()-1;
		}
		if (dummy == "#THIRDBODYSINGLE#")
		{
			stoichiometry.pop_back();
			parsed_string >> dummy;
			kinetics.reaction[iReaction].thirdbodysingle = dummy;
		}	
	}

	// Reactants
	{
		bool iWarning = false;
		for(i=1;i<=index-1;i++)
		{
			string name;
			double number;

			SeparateNumberFromString(stoichiometry[i], name, number);
			if (name != ""  && name != "+" && name != "M")			kinetics.reaction[iReaction].nameDirect.push_back(name);
			if (name != "+" && iWarning == false && name != "M")	kinetics.reaction[iReaction].nuDirect.Append(number);

			if (name == "")		iWarning = true;
			else				iWarning = false;

			if (name == "M")	kinetics.reaction[iReaction].iThirdBody = true;
		}

		if (kinetics.reaction[iReaction].nameDirect.size()-1 != kinetics.reaction[iReaction].nuDirect.Size())
			ErrorMessage(iLine, "Error in reaction stoichiometry (reactant side)");

		kinetics.reaction[iReaction].nReactants = kinetics.reaction[iReaction].nameDirect.size()-1;

		//for(i=1;i<=kinetics.reaction[iReaction].nReactants;i++)
		//	cout << "R " << kinetics.reaction[iReaction].nameDirect[i] << " " << kinetics.reaction[iReaction].nuDirect[i] << endl;
	
		for(i=1;i<=kinetics.reaction[iReaction].nReactants;i++)
			kinetics.reaction[iReaction].indexDirect.Append(kinetics.RecognizeSpecies(kinetics.reaction[iReaction].nameDirect[i]));

	}

	// Products
	{
		bool iWarning = false;
		for(i=index+1;i<=int(stoichiometry.size())-1;i++)
		{
			string name;
			double number;

			SeparateNumberFromString(stoichiometry[i], name, number);
			if (name != ""  && name != "+" && name != "M")			kinetics.reaction[iReaction].nameInverse.push_back(name);
			if (name != "+" && iWarning == false && name != "M")	kinetics.reaction[iReaction].nuInverse.Append(number);

			if (name == "")		iWarning = true;
			else				iWarning = false;

			if (name == "M")	kinetics.reaction[iReaction].iThirdBody = true;
		}

		if (kinetics.reaction[iReaction].nameInverse.size()-1 != kinetics.reaction[iReaction].nuInverse.Size())
			ErrorMessage(iLine, "Error in reaction stoichiometry (product side)");

		kinetics.reaction[iReaction].nProducts = kinetics.reaction[iReaction].nameInverse.size()-1;

		//for(i=1;i<=kinetics.reaction[iReaction].nProducts;i++)
		//	cout << "P " << kinetics.reaction[iReaction].nameInverse[i] << " " << kinetics.reaction[iReaction].nuInverse[i] << endl;

		for(i=1;i<=kinetics.reaction[iReaction].nProducts;i++)
			kinetics.reaction[iReaction].indexInverse.Append(kinetics.RecognizeSpecies(kinetics.reaction[iReaction].nameInverse[i]));
	}

	kinetics.reaction[iReaction].CompactStoichiometry();
}

void OpenSMOKE_CHEMKINInterpreter::ParsingReactionsAdditional(const int i, const int iLine, vector<string> instructions)
{
	int number_of_instructions = instructions.size()-1;
	if (number_of_instructions >= 1)
	{
		string instruction = instructions[1];

		     if (caseInsCompare("CHEB", instruction))
			kinetics.reaction[i].Chebishev(instructions);
		else if (caseInsCompare("COLLEFF", instruction))
			kinetics.reaction[i].CollisionalEfficiency(instructions);
		else if (caseInsCompare("DUP", instruction))
			kinetics.reaction[i].Duplicate(instructions);
		else if (caseInsCompare("DUPLICATE", instruction))
			kinetics.reaction[i].Duplicate(instructions);
		else if (caseInsCompare("EXCI", instruction))
			kinetics.reaction[i].EnergyLossParameter(instructions);
		else if (caseInsCompare("FIT1", instruction))
			kinetics.reaction[i].Fit1(instructions);
		else if (caseInsCompare("FORD", instruction))
			kinetics.reaction[i].ForwardReactionKineticParameter(instructions);
		else if (caseInsCompare("HIGH", instruction))
			kinetics.reaction[i].HighPressureLimit(instructions);
		else if (caseInsCompare("JAN", instruction))
			kinetics.reaction[i].JanevLangerReactionRate(instructions);
		else if (caseInsCompare("LOW", instruction))
			kinetics.reaction[i].LowPressureLimit(instructions);
		else if (caseInsCompare("LT", instruction))
			kinetics.reaction[i].LandauTellerReaction(instructions);
		else if (caseInsCompare("MOME", instruction))
			kinetics.reaction[i].MomentumTransferCollisionFrequency(instructions);
		else if (caseInsCompare("PCHEB", instruction))
			kinetics.reaction[i].ChebishevPressureLimits(instructions);
		else if (caseInsCompare("PLOG", instruction))
			kinetics.reaction[i].PressureLogarithmicInterpolation(instructions);
		else if (caseInsCompare("REV", instruction))
			kinetics.reaction[i].ReverseRate(instructions);
		else if (caseInsCompare("RLT", instruction))
			kinetics.reaction[i].LandauTellerReverseReaction(instructions);
		else if (caseInsCompare("RORD", instruction))
			kinetics.reaction[i].BackwardReactionKineticParameter(instructions);
		else if (caseInsCompare("SRI", instruction))
			kinetics.reaction[i].SRIFallOffReaction(instructions);
		else if (caseInsCompare("TCHEB", instruction))
			kinetics.reaction[i].ChebishevTemperatureLimits(instructions);
		else if (caseInsCompare("TDEP", instruction))
			kinetics.reaction[i].SpeciesTemperatureDependence(instructions);
		else if (caseInsCompare("TROE", instruction))
			kinetics.reaction[i].TROEFallOffReaction(instructions);
		else if (caseInsCompare("UNITS", instruction))
			kinetics.reaction[i].Units(instructions);
		else if (caseInsCompare("USRPROG", instruction))
			kinetics.reaction[i].UserRateSubRoutine(instructions);
		else if (caseInsCompare("XSMI", instruction))
			kinetics.reaction[i].CollisionCrossSection(instructions);

		else if (caseInsCompare("CONSISTENCY", instruction))
			kinetics.reaction[i].ThermodynamicConsistency(instructions);

		else if (caseInsCompare("REACTAR", instruction))
			kinetics.reaction[i].TAR(instructions);
	
		else
			kinetics.reaction[i].ThirdBodyEfficiencies(instructions);
	}
}

void OpenSMOKE_CHEMKINInterpreter::SurfaceParsingReactionsAdditional(const int k, const int i, const int iLine, vector<string> instructions)
{
	int number_of_instructions = instructions.size()-1;
	if (number_of_instructions >= 1)
	{
		string instruction = instructions[1];

		if (caseInsCompare("STICK", instruction))
			material[k].surfaceKinetics.reaction[i].StickyReaction(instructions, thermo, material[k], indices_species_thermo_database, indices_site_thermo_database);
		else if (caseInsCompare("COV", instruction))
			material[k].surfaceKinetics.reaction[i].CoverageDependentReactions(instructions);
		else if (caseInsCompare("LANG", instruction))
			material[k].surfaceKinetics.reaction[i].LangmuirHinshelwoodReactions(instructions);
		else if (caseInsCompare("LHDE", instruction))
			material[k].surfaceKinetics.reaction[i].LangmuirHinshelwoodDenominatorExponentParameter(instructions);
//		else if (caseInsCompare("LHPR", instruction))
//			material[k].surfaceKinetics.reaction[i].LangmuirHinshelwoodEquilibriumPressure(instructions);
		else { cout << instruction << endl; ErrorMessage(1, "This keyword is not supported"); }
	}
}

OpenSMOKE_CHEMKINInterpreter::~OpenSMOKE_CHEMKINInterpreter()
{

}

void OpenSMOKE_CHEMKINInterpreter::AddSpaces(string &line, const char symbol)
{
	int length = line.size();

	int k=0;
	for(;;)
	{
		if (k>=length)	break;
		if (line.at(k) == symbol)
		{
			line.insert(k,1,' ');
			line.insert(k+2,1,' ');
			k+=2;
			length+=2;
		}
		k++;
	}
}

void OpenSMOKE_CHEMKINInterpreter::PrintIdealGasBinaryFile(const string file_name)
{
	OpenSMOKE_PreProcessorIdealGas preprocessor;
	preprocessor.SetupFromCHEMKINInterpreter(this, file_name);

	if (iVerboseMode == true)
	{/*
		{
			cout << " ** Writing New OutputFile..." << endl;

			ofstream fNewOutputFile;
			openOutputFileAndControl(fNewOutputFile, name_output_folder + "/NewOutputFile.out");

			preprocessor.WriteNewOutputFile(fNewOutputFile);

			fNewOutputFile.close();
		}
		*/
		{
			cout << " ** Writing ThermoCoefficients.out file..." << endl;

			ofstream fThermoCoefficients;
			openOutputFileAndControl(fThermoCoefficients, name_output_folder + "/ThermoCoefficients.out");
			
			preprocessor.WriteCpCoefficients(fThermoCoefficients);
			preprocessor.WriteEnthalpyCoefficients(fThermoCoefficients);
			preprocessor.WriteEntropyCoefficients(fThermoCoefficients);

			fThermoCoefficients.close();
		}

		if (transport.IsActivated() == true)
		{
			cout << " ** Writing FittingCoefficients.out file..." << endl;

			ofstream fFitting;
			openOutputFileAndControl(fFitting, name_output_folder + "/FittingCoefficients.out");
			
			preprocessor.WriteViscosityFittingCoefficients(fFitting);
			preprocessor.WriteConductivityFittingCoefficients(fFitting);
			preprocessor.WriteDiffusivityFittingCoefficients(fFitting);
			preprocessor.WriteThermalDiffusionRatiosFittingCoefficients(fFitting);

			fFitting.close();
		}
	}
}

void OpenSMOKE_CHEMKINInterpreter::WriteAdditionalFiles()
{
}

void OpenSMOKE_CHEMKINInterpreter::PrintSummaryFile(const string file_name)
{	
	ofstream fSummary;
	openOutputFileAndControl(fSummary, file_name);

	// Elements
	elements.SummaryOnFile(fSummary);
	
	// Gas thermodynamics
	thermo.SummaryOnFile(fSummary, indices_species_thermo_database, elements);
	
	// Surface thermodynamics
	if (iSurfaceKinetics == true)
		thermo.SummaryOnFile(fSummary, indices_site_thermo_database, indices_bulk_thermo_database, elements);
	
	// Gas transport properties
	if (transport.IsActivated() == true)
		transport.SummaryOnFile(fSummary, indices_species_transport_database);
	
	// Gas Kinetics
	kinetics.SummaryOnFile(fSummary);
	for(int i=1;i<=kinetics.number_of_reactions;i++)
			kinetics.reaction[i].SummaryOnFile(fSummary);

	// Surface Kinetics
	if (iSurfaceKinetics == true)
	{
		for (int k=1;k<=nSurfaceMaterials;k++)
		{
			material[k].surfaceKinetics.SummaryOnFile(fSummary);
			for(int i=1;i<=material[k].surfaceKinetics.number_of_reactions;i++)
					material[k].surfaceKinetics.reaction[i].SummaryOnFile(fSummary);
		}
	}
	fSummary.close();
}

void OpenSMOKE_CHEMKINInterpreter::PrintBuzziFiles()
{
	string MSDOS_command = "mkdir " + name_ascii_folder;
	system(MSDOS_command.c_str());
	
	kinetics.PrintNamesFile(name_ascii_folder + "/names.bzz");
	kinetics.PrintStoichiometryFile(name_ascii_folder + "/stoichiometry.bzz");
	kinetics.PrintReactionsFile(name_ascii_folder + "/reactions.bzz");
	thermo.PrintThermoFile(name_ascii_folder + "/termo.bzz", indices_species_thermo_database);
}

void OpenSMOKE_CHEMKINInterpreter::ParsingReactions()
{
	int i;

	for(i=1;i<=number_of_lines;i++)
		if (StringFindSubString(lines[i], "#REVERSIBLE#") || StringFindSubString(lines[i], "#DIRECT#"))
			kinetics.lineReactionIndex.Append(i);
	kinetics.number_of_reactions = kinetics.lineReactionIndex.Size();
	kinetics.reaction = new OpenSMOKE_CHEMKINInterpreter_ReactionData[2*kinetics.number_of_reactions+1];

	// Initialize Reactions
	for(i=1;i<=kinetics.number_of_reactions;i++)
		kinetics.reaction[i].Setup(i, kinetics.lineReactionIndex[i], &kinetics, &units);

	// Initialize Reactions: default units
	bool iEnergy;
	double energy_factor;
	double frequency_factor;
	units.GiveMeConversionFactor(units.current_energy, energy_factor, iEnergy);
	units.GiveMeConversionFactor(units.current_frequency_factor, frequency_factor, iEnergy);
	for(i=1;i<=kinetics.number_of_reactions;i++)
		kinetics.reaction[i].Setup(energy_factor, frequency_factor);
}

void OpenSMOKE_CHEMKINInterpreter::WriteReducedScheme(string fileName)
{
	int n;
	string checkstring;

	ifstream fInput;
	openInputFileAndControl(fInput, fileName);
	fInput >> n;
	BzzVectorInt iTarget(n);
	for (int i=1;i<=n;i++)
		fInput >> iTarget[i];
	fInput >> checkstring;
	if (checkstring != "END")
		ErrorMessage(0, "Wrong list of reactions to mantain");
	fInput.close();

	ofstream fOutput;
	openOutputFileAndControl(fOutput, "ReducedScheme.out");

	int counter = 0;
	for(int k=1;k<=iTarget.Size();k++)
	{
		cout << "Search " << iTarget[k] << endl;
		for(int i=1;i<=kinetics.number_of_reactions;i++)
		{
			if (i==iTarget[k])
			{
				cout << "Writing reaction " << i << endl;

				fOutput << "! " << iTarget[k] << endl;
				if (i < kinetics.number_of_reactions)
					for (int j=1;j<=kinetics.lineReactionIndex[i+1]-kinetics.lineReactionIndex[i];j++)
					{
						string line = lines[kinetics.lineReactionIndex[i]+j-1];
						StringSubstitutionAll(line, " #DIRECT# ", " => ");
						StringSubstitutionAll(line, " #REVERSIBLE# ", " = ");
						StringSubstitutionAll(line, " + M ", "(+M)");
						StringSubstitutionAll(line, " + M ", "(+m)");
						fOutput << line << endl;
					}
				else if (i == kinetics.number_of_reactions)
					for (int j=1;j<=lines.size()-kinetics.lineReactionIndex[i]-1;j++)
					{
						string line = lines[kinetics.lineReactionIndex[i]+j-1];
						StringSubstitutionAll(line, " #DIRECT# ", " => ");
						StringSubstitutionAll(line, " #REVERSIBLE# ", " = ");
						StringSubstitutionAll(line, " + M ", "(+M)");
						StringSubstitutionAll(line, " + M ", "(+m)");
						fOutput << line << endl;
					}
				else if (i > kinetics.number_of_reactions)
					ErrorMessage(i, "Wrong index in reduced reaction list");

				counter++;
				fOutput << endl;
				
				break;
			}
		}
	}

	fOutput.close();

	if (counter != iTarget.Size())
		ErrorMessage(0, "Something wrong in writing ReducedScheme.out");
}

void OpenSMOKE_CHEMKINInterpreter::ParsingReactionsAdditionalData()
{
	int i;

	for(i=1;i<=kinetics.number_of_reactions-1;i++)
		kinetics.additionalLineReactionIndex.Append(kinetics.lineReactionIndex[i+1]-kinetics.lineReactionIndex[i]-1);
	for(i=kinetics.lineReactionIndex[kinetics.number_of_reactions]+1;i<=number_of_lines;i++)
		if (StringFindSubString(lines[i], "END") == true)
			kinetics.additionalLineReactionIndex.Append(i-kinetics.lineReactionIndex[kinetics.number_of_reactions]-1);

	for(i=1;i<=kinetics.number_of_reactions;i++)
	{
		int j = kinetics.lineReactionIndex[i];

		vector<string> instructions;
		instructions.push_back("instructions");
	
		string dummy;
		stringstream parsed_string(lines[j]);
		for(;;)
		{
			parsed_string >> dummy;
			if (parsed_string.fail())
				break;
			instructions.push_back(dummy);
		}

		ParsingReactions(i, 1, instructions);
		kinetics.reaction[i].ReactionString();

		// Checking stoichiometry
		kinetics.reaction[i].CheckStoichiometry(thermo, indices_species_thermo_database, elements.elements_in_list, elements.isotope_in_list);
	}

	// Reactions parsing
	for(i=1;i<=kinetics.number_of_reactions;i++)
	{
		int j = kinetics.lineReactionIndex[i];

		for(int k=1;k<=kinetics.additionalLineReactionIndex[i];k++)
		{
			vector<string> instructions;
			instructions.push_back("instructions");
	
			string dummy;
			stringstream parsed_string(lines[j+k]);
			for(;;)
			{
				parsed_string >> dummy;
				if (parsed_string.fail())
					break;
				instructions.push_back(dummy);
			}
			ParsingReactionsAdditional(i, j, instructions);
		}
	}
}

void OpenSMOKE_CHEMKINInterpreter::ParsingSpecies()
{
	vector<string> instructions;
	instructions.push_back("instructions");
	
	for(int i=startSpecies;i<=startReactions-1;i++)
	{
		bool iSkip = false;

		for(int j=1;j<=finalUnits.Size();j++)
			if (i>=startUnits[j] && i<=finalUnits[j])	iSkip = true;

		if (iSkip == false)
		{
			string dummy;
			stringstream parsed_string(lines[i]);

			for(;;)
			{
				parsed_string >> dummy;
				if (parsed_string.fail())
					break;
				instructions.push_back(dummy);
			}
		}
	}
	
	ParsingSpecies(1, instructions);
	kinetics.species_list = species.species_in_list;
	species.Summary();
}

void OpenSMOKE_CHEMKINInterpreter::ParsingElements()
{
	vector<string> instructions;
	instructions.push_back("instructions");
	
	for(int i=startElements;i<=startSpecies-1;i++)
	{
		bool iSkip = false;
		for(int j=1;j<=finalUnits.Size();j++)
			if (i>=startUnits[j] && i<=finalUnits[j])	iSkip = true;

		if (iSkip == false)
		{
			string dummy;
			stringstream parsed_string(lines[i]);

			for(;;)
			{
				parsed_string >> dummy;
				if (parsed_string.fail())
					break;
				instructions.push_back(dummy);
			}
		}
	}

	ParsingElements(1, instructions);
	elements.Summary();
}

void OpenSMOKE_CHEMKINInterpreter::ParsingUnits()
{
	for(int i=1;i<=startUnits.Size();i++)
		if( startUnits[i] < startReactions)
		{
			int j;
			for(j=startUnits[i]+1;j<=number_of_lines;j++)
			{
				string dummy;
				stringstream parsed_string(lines[j]);

				parsed_string >> dummy;

				if (caseInsCompare("ELEM", dummy) || caseInsCompare("ELEMENTS", dummy))	break;
				if (caseInsCompare("UNITS", dummy) == true)								break;
				if (caseInsCompare("SPECIES", dummy) || caseInsCompare("SPEC", dummy))	break;
				if (caseInsCompare("REACTIONS", dummy) == true)							break;
			}
			
			finalUnits.Append(j-1);
			vector<string> instructions;
			instructions.push_back("instructions");
			string dummy;
			for(j=startUnits[i];j<=finalUnits[finalUnits.Size()];j++)
			{
				stringstream parsed_string(lines[j]);
				for(;;)
				{
					parsed_string >> dummy;
					if (parsed_string.fail())	break;
					instructions.push_back(dummy);
				}
			}

			if (instructions[1] != "UNITS")								ErrorMessage(startUnits[i], "Syntax error in UNITS declaration");
			if (instructions.size()-1 > 3)								ErrorMessage(startUnits[i], "Syntax error in UNITS declaration");
			if (instructions.size()-1 == 3 && instructions[3] != "END")	ErrorMessage(startUnits[i], "Syntax error in UNITS declaration");
			units.Parse_Units(instructions[2]);
		}
		

	vector<string> instructions;
	instructions.push_back("instructions");
	
	string dummy;
	stringstream parsed_string(lines[startReactions]);

	for(;;)
	{
		parsed_string >> dummy;
		if (parsed_string.fail())
			break;
		instructions.push_back(dummy);
	}

	for(int j=2;j<=int(instructions.size())-1;j++)
		units.Parse_Units(instructions[j]);
	
	units.Summary();
}

void OpenSMOKE_CHEMKINInterpreter::CleaningLines()
{
	for(int i=1;i<=number_of_lines;i++)
	{
		for(int j=0;j<int(lines[i].size());j++)
			if (lines[i].at(j) == '\t')	lines[i].at(j) = ' ';

		AddSpaces(lines[i], '/');
	}
}

void OpenSMOKE_CHEMKINInterpreter::SurfaceCleaningLines()
{
	for(int i=1;i<=surface_number_of_lines;i++)
	{
		for(int j=0;j<int(surface_lines[i].size());j++)
			if (surface_lines[i].at(j) == '\t')	surface_lines[i].at(j) = ' ';

		AddSpaces(surface_lines[i], '/');
	}
}

void OpenSMOKE_CHEMKINInterpreter::SyntaxCorrections()
{
	for(int i=1;i<=number_of_lines;i++)
	{
		int j;

		for(j=1;j<=int(units.energy_units.size())-1;j++)
			StringSubstitutionAll(lines[i], units.energy_units_original[j], units.energy_units[j]);

		for(j=1;j<=int(units.frequency_factor_units.size())-1;j++)
			StringSubstitutionAll(lines[i], units.frequency_factor_units_original[j], units.frequency_factor_units[j]);

		StringSubstitutionAll(lines[i], "<=>",  " #REVERSIBLE# ");
		StringSubstitutionAll(lines[i], "=>",   " #DIRECT# ");
		StringSubstitutionAll(lines[i], "=",    " #REVERSIBLE# ");
		StringSubstitutionAll(lines[i], "(+M)", " + M ");
		StringSubstitutionAll(lines[i], "(+m)", " + M ");

		bool found_thirdbodysingle = StringSubstitutionAll(lines[i], "(+",   " #THIRDBODYSINGLE# ");
		if (found_thirdbodysingle==true) StringSubstitutionAll(lines[i], ")",    " ");
	}
}

void OpenSMOKE_CHEMKINInterpreter::SurfaceSyntaxCorrections()
{
	for(int i=1;i<=surface_number_of_lines;i++)
	{
		StringSubstitutionAll(surface_lines[i], "<=>",      " #REVERSIBLE# ");
		StringSubstitutionAll(surface_lines[i], "=>",       " #DIRECT# ");
		StringSubstitutionAll(surface_lines[i], "=",        " #REVERSIBLE# ");
		StringSubstitutionAll(surface_lines[i], "SITESDEN", " SITE / SDEN ");
		StringSubstitutionAll(surface_lines[i], "SITE/",    " SITE / ");
		StringSubstitutionAll(surface_lines[i], "BULK/",    " BULK / ");
//		StringSubstitutionAll(surface_lines[i], "END",      " END ");
	}
}

void OpenSMOKE_CHEMKINInterpreter::ParsingSections()
{
	for(int i=1;i<=number_of_lines;i++)
	{
		string dummy;
		stringstream parsed_string(lines[i]);

		parsed_string >> dummy;
		
		if (startElements == 0)
			if (caseInsCompare("ELEM", dummy) || caseInsCompare("ELEMENTS", dummy))	
				startElements = i;

		if (caseInsCompare("UNITS", dummy) == true)	
			startUnits.Append(i);

		if (startSpecies == 0)
			if (caseInsCompare("SPECIES", dummy) || caseInsCompare("SPEC", dummy))	
				startSpecies = i;

		if (startReactions == 0)
			if (caseInsCompare("REACTIONS", dummy) == true)	
				startReactions = i;
	}
}

void OpenSMOKE_CHEMKINInterpreter::SurfaceParsingSections()
{
	BzzVectorInt surface_start_material;
	for(int i=1;i<=surface_number_of_lines;i++)
	{
		string dummy;
		stringstream parsed_string(surface_lines[i]);

		parsed_string >> dummy;
		
		if (caseInsCompare("MATERIAL", dummy))	
			surface_start_material.Append(i);
	}

	// Allocate memory for materials
	nSurfaceMaterials = surface_start_material.Size();
	material = new OpenSMOKE_PreProcessorSurfaceMaterial[nSurfaceMaterials+1];

	for(int k=1;k<=nSurfaceMaterials;k++)
	{
		material[k].start_material_section = surface_start_material[k];
		material[k].end_material_section = surface_number_of_lines;
		if (k<surface_start_material.Size()) material[k].end_material_section = surface_start_material[k+1]-1;

		material[k].nSurfaceSites = 0;
		material[k].nSurfaceBulks = 0;
		material[k].nSurfaceReactions = 0;

		for(int i=material[k].start_material_section;i<=material[k].end_material_section;i++)
		{
			string dummy;
			stringstream parsed_string(surface_lines[i]);

			parsed_string >> dummy;
		
			if (caseInsCompare("SITE", dummy))	
			{
				material[k].start_sites.Append(i);
				material[k].nSurfaceSites++;
			}

			if (caseInsCompare("BULK", dummy))	
			{
				material[k].start_bulks.Append(i);
				material[k].nSurfaceBulks++;
			}

			if (caseInsCompare("REACTIONS", dummy))	
			{
				material[k].start_reactions.Append(i);
				material[k].nSurfaceReactions++;
			}
		}
	}
}

void OpenSMOKE_CHEMKINInterpreter::ReadKineticScheme()
{
	char comment[Constants::COMMENT_SIZE];

	ifstream fInput;
	openInputFileAndControl(fInput, name_kinetics_file);

	lines.push_back("List of lines");
	while(!fInput.eof())
	{
		fInput.getline(comment, Constants::COMMENT_SIZE);
		lines.push_back(comment);
	}
	fInput.close();

	number_of_lines = lines.size()-1;

	for(int i=1;i<=number_of_lines;i++)
	{
		if (CheckForBlankLine(lines[i]) == true)		indexBlankLines.Append(i);
		else if (CheckForCommentLine(lines[i]) == true)	indexCommentLines.Append(i);
		else											indexLines.Append(i);
	}

	fLog << " ----------------------------------------------------------------"  << endl;
	fLog << "                     Kinetic Scheme File                         "  << endl;
	fLog << " ----------------------------------------------------------------"  << endl;
	fLog << "Total number of full lines:    " << indexLines.Size()				 << endl;
	fLog << "Total number of blank lines:   " << indexBlankLines.Size()			 << endl;
	fLog << "Total number of comment lines: " << indexCommentLines.Size()		 << endl;
	fLog << "Total number of lines:         " << number_of_lines << endl;
	fLog << " ----------------------------------------------------------------"  << endl;
}

void OpenSMOKE_CHEMKINInterpreter::ReadSurfaceKineticScheme()
{
	char comment[Constants::COMMENT_SIZE];

	ifstream fInput;
	openInputFileAndControl(fInput, name_surface_kinetics_file);

	surface_lines.push_back("List of lines");
	while(!fInput.eof())
	{
		fInput.getline(comment, Constants::COMMENT_SIZE);
		surface_lines.push_back(comment);
	}
	fInput.close();

	surface_number_of_lines = surface_lines.size()-1;

	

	for(int i=1;i<=surface_number_of_lines;i++)
	{
		if (CheckForBlankLine(surface_lines[i]) == true)		surfaceIndexBlankLines.Append(i);
		else if (CheckForCommentLine(surface_lines[i]) == true)	surfaceIndexCommentLines.Append(i);
		else													surfaceIndexLines.Append(i);
	}

	fLog << " ----------------------------------------------------------------"  << endl;
	fLog << "                  Surface Kinetic Scheme File                    "  << endl;
	fLog << " ----------------------------------------------------------------"  << endl;
	fLog << "Total number of full lines:    " << surfaceIndexLines.Size()				 << endl;
	fLog << "Total number of blank lines:   " << surfaceIndexBlankLines.Size()			 << endl;
	fLog << "Total number of comment lines: " << surfaceIndexCommentLines.Size()		 << endl;
	fLog << "Total number of lines:         " << surface_number_of_lines << endl;
	fLog << " ----------------------------------------------------------------"  << endl;
}

void OpenSMOKE_CHEMKINInterpreter::SurfaceParsingMaterials()
{
	for(int k=1;k<=nSurfaceMaterials;k++)
	{
		string dummy;
		stringstream parsed_string(surface_lines[material[k].start_material_section]);

		// Read Material Key word
		parsed_string >> dummy;
		
		// Optional: read name
		parsed_string >> dummy;
		if (parsed_string.fail())
		{
			stringstream index;
			index << k;
			material[k].name = "MATERIAL" + index.str();
		}
		else
			material[k].name = dummy;

		// Gas species (list)
		material[k].surfaceKinetics.gas_species_list = kinetics.species_list;

		material[k].surfaceKinetics.Setup(&fLog, &fWarning);
	}
}

void OpenSMOKE_CHEMKINInterpreter::SurfaceParsingSites()
{
	// Loop over materials
	for(int k=1;k<=nSurfaceMaterials;k++)
	{
		material[k].site_name.resize(material[k].nSurfaceSites+1);
		ChangeDimensions(material[k].nSurfaceSites+1, &material[k].site_density);
		
		material[k].site_occupancy = new BzzVector[material[k].nSurfaceSites+1];
		material[k].site_species_names = new vector<string>[material[k].nSurfaceSites+1];

		// Loop over sites
		for(int j=1;j<=material[k].nSurfaceSites;j++)
		{
			vector<string> instructions;
			instructions.push_back("instructions");
			material[k].site_species_names[j].push_back("species");
			
			// Recover instructions
			{
				int start_line = material[k].start_sites[j];
				int stop_line = material[k].start_reactions[1]-1;

				if (j<material[k].nSurfaceSites)
					stop_line = material[k].start_sites[j+1]-1;
				if (j==material[k].nSurfaceSites)
					if (material[k].start_bulks.Size() > 0)
						stop_line = material[k].start_bulks[1] -1;

				string site_lines;
				for(int i=start_line;i<=stop_line;i++)
					site_lines += surface_lines[i] + " ";

				string dummy;
				stringstream parsed_string(site_lines);

				// Read Material Key word
				for(;;)
				{
					parsed_string >> dummy;
					if (parsed_string.fail())
						break;
					else
						instructions.push_back(dummy);
				}

				// Read and interpret instructions
				int n = instructions.size()-1;
				if (instructions[2] == "/")
					if (instructions[4] != "/" || instructions[6] != "/")
						ErrorMessage(start_line, "Error in site definition...");
				
				int count=0;
				if (instructions[3] == "SDEN")
				{
					stringstream index;
					index << j;
					material[k].site_name[j] = "SITE" + index.str();
					material[k].site_density[j] = atof(instructions[5].c_str());
					count = 7;
				}
				else
				{
					if (instructions[5] != "SDEN" || instructions[8] != "/")
						ErrorMessage(start_line, "Error in site definition...");
					material[k].site_name[j] = instructions[3];
					material[k].site_density[j] = atof(instructions[7].c_str());
					count = 9;
				}

				for(;;)
				{
					if (count > n)
						break;

					if (instructions[count] == "END")
						break;
			
					material[k].site_species_names[j].push_back(instructions[count]);
					material[k].site_occupancy[j].Append(1.);
					count++;

					if (count <= n)
					{
						if (instructions[count] == "/")
						{
							material[k].site_occupancy[j][material[k].site_occupancy[j].Size()] = atof(instructions[count+1].c_str());
							if (instructions[count+2] != "/")
								ErrorMessage(start_line, "Error in site definition...");
							count+=3;
						}
					}

				}
			}

		}
	}

	// Build list of site species
	for(int k=1;k<=nSurfaceMaterials;k++)
	{
		material[k].surfaceKinetics.site_species_list.resize(0);
		material[k].surfaceKinetics.site_species_list.push_back("site_species_list");

		for(int j=1;j<=material[k].nSurfaceSites;j++)
		{
			for (int i=1;i<=material[k].site_occupancy[j].Size();i++)
			{
				bool iAlreadyAvailable = false;
				string focus = material[k].site_species_names[j][i];
				for(int n=1;n<=material[k].surfaceKinetics.site_species_list.size()-1;n++)
					if (material[k].surfaceKinetics.site_species_list[n] == focus)
					{
						iAlreadyAvailable=true;
						break;
					}
				if (iAlreadyAvailable == false)
					material[k].surfaceKinetics.site_species_list.push_back(focus);
			}
		}

		// Occupancy matrix
		ChangeDimensions(material[k].surfaceKinetics.site_species_list.size()-1, material[k].nSurfaceSites, &material[k].site_occupancy_matrix);
		for(int j=1;j<=material[k].nSurfaceSites;j++)
		{
			for (int n=1;n<=material[k].surfaceKinetics.site_species_list.size()-1;n++)
				for (int i=1;i<=material[k].site_occupancy[j].Size();i++)
				{
					if (material[k].surfaceKinetics.site_species_list[n] == material[k].site_species_names[j][i])
						material[k].site_occupancy_matrix[n][j] = material[k].site_occupancy[j][i];
				}
		}
	}


}

void OpenSMOKE_CHEMKINInterpreter::SurfaceParsingBulk()
{
	// Loop over materials
	for(int k=1;k<=nSurfaceMaterials;k++)
	{
		material[k].bulk_name.resize(material[k].nSurfaceBulks+1);
		
		material[k].bulk_density = new BzzVector[material[k].nSurfaceBulks+1];
		material[k].bulk_species_names = new vector<string>[material[k].nSurfaceBulks+1];

		// Loop over bulks
		for(int j=1;j<=material[k].nSurfaceBulks;j++)
		{
			vector<string> instructions;
			instructions.push_back("instructions");
			material[k].bulk_species_names[j].push_back("species");
			
			// Recover instructions
			{
				int start_line = material[k].start_bulks[j];
				int stop_line = material[k].start_reactions[1]-1;

				if (j<material[k].nSurfaceBulks)
					stop_line = material[k].start_bulks[j+1]-1;

				string bulk_lines;
				for(int i=start_line;i<=stop_line;i++)
					bulk_lines += surface_lines[i] + " ";

				string dummy;
				stringstream parsed_string(bulk_lines);

				// Read Material Key word
				for(;;)
				{
					parsed_string >> dummy;
					if (parsed_string.fail())
						break;
					else
						instructions.push_back(dummy);
				}

				// Read and interpret instructions
				int n = instructions.size()-1;

				int count=0;
				if (instructions[2] == "/")
				{
					if (instructions[4] != "/")
						ErrorMessage(start_line, "Error in bulk definition...");
					material[k].bulk_name[j] = instructions[3];
					count = 5;
				}
				else
				{
					stringstream index;
					index << j;
					material[k].bulk_name[j] = "BULK" + index.str();
					count = 2;
				}

				for(;;)
				{
					if (count > n)
						break;

					if (instructions[count] == "END")
						break;

					material[k].bulk_species_names[j].push_back(instructions[count]);
					material[k].bulk_density[j].Append(-1.);
					count++;

					if (count <= n)
					{
						if (instructions[count] == "/")
						{
							material[k].bulk_density[j][material[k].bulk_density[j].Size()] = atof(instructions[count+1].c_str());
							if (instructions[count+2] != "/")
								ErrorMessage(start_line, "Error in bulk definition...");
							count+=3;
						}
					}
				}
			}
		}
	}

	// Build list of bulk species
	for(int k=1;k<=nSurfaceMaterials;k++)
	{
		material[k].surfaceKinetics.bulk_species_list.resize(0);
		material[k].surfaceKinetics.bulk_species_list.push_back("bulk_species_list");

		for(int j=1;j<=material[k].nSurfaceBulks;j++)
		{
			for (int i=1;i<=material[k].bulk_density[j].Size();i++)
			{
				bool iAlreadyAvailable = false;
				string focus = material[k].bulk_species_names[j][i];
				for(int n=1;n<=material[k].surfaceKinetics.bulk_species_list.size()-1;n++)
					if (material[k].surfaceKinetics.bulk_species_list[n] == focus)
					{
						iAlreadyAvailable=true;
						break;
					}
				if (iAlreadyAvailable == false)
					material[k].surfaceKinetics.bulk_species_list.push_back(focus);
			}
		}

		// Occupancy matrix
		ChangeDimensions(material[k].surfaceKinetics.bulk_species_list.size()-1, material[k].nSurfaceBulks, &material[k].bulk_density_matrix);
		for(int j=1;j<=material[k].nSurfaceBulks;j++)
		{
			for (int n=1;n<=material[k].surfaceKinetics.bulk_species_list.size()-1;n++)
				for (int i=1;i<=material[k].bulk_density[j].Size();i++)
				{
					if (material[k].surfaceKinetics.bulk_species_list[n] == material[k].bulk_species_names[j][i])
						material[k].bulk_density_matrix[n][j] = material[k].bulk_density[j][i];
				}
		}

	}


	// Global list
	{
		bulk_species_global_list.resize(0);
		bulk_species_global_list.push_back("bulk_species_global_list");
		for(int k=1;k<=nSurfaceMaterials;k++)
		{
			for (int i=1;i<=material[k].surfaceKinetics.bulk_species_list.size()-1;i++)
			{
				bool iAlreadyAvailable = false;
				string focus = material[k].surfaceKinetics.bulk_species_list[i];
				for(int n=1;n<=bulk_species_global_list.size()-1;n++)
					if (bulk_species_global_list[n] == focus)
					{
						iAlreadyAvailable=true;
						break;
					}
				if (iAlreadyAvailable == false)
					bulk_species_global_list.push_back(focus);
				}
		}
	}

	// Site list
	{
		site_species_global_list.resize(0);
		site_species_global_list.push_back("site_species_global_list");
		for(int k=1;k<=nSurfaceMaterials;k++)
		{
			for (int i=1;i<=material[k].surfaceKinetics.site_species_list.size()-1;i++)
			{
				bool iAlreadyAvailable = false;
				string focus = material[k].surfaceKinetics.site_species_list[i];
				for(int n=1;n<=site_species_global_list.size()-1;n++)
					if (site_species_global_list[n] == focus)
					{
						iAlreadyAvailable=true;
						break;
					}
				if (iAlreadyAvailable == false)
					site_species_global_list.push_back(focus);
				}
		}
	}
}

void OpenSMOKE_CHEMKINInterpreter::SurfaceParsingReactions()
{
	int i;

	// Loop over materials
	for(int k=1;k<=nSurfaceMaterials;k++)
	{
		for(i=material[k].start_reactions[1];i<=material[k].end_material_section;i++)
			if (StringFindSubString(surface_lines[i], "#REVERSIBLE#") || StringFindSubString(surface_lines[i], "#DIRECT#"))
				material[k].surfaceKinetics.lineReactionIndex.Append(i);
		material[k].surfaceKinetics.number_of_reactions = material[k].surfaceKinetics.lineReactionIndex.Size();
		material[k].surfaceKinetics.reaction = new OpenSMOKE_CHEMKINInterpreter_SurfaceReactionData[2*material[k].surfaceKinetics.number_of_reactions+1];

		// Initialize Reactions
		for(i=1;i<=material[k].surfaceKinetics.number_of_reactions;i++)
			material[k].surfaceKinetics.reaction[i].Setup(i, material[k].surfaceKinetics.lineReactionIndex[i], &material[k].surfaceKinetics, &units);

		// Initialize Reactions: default units
	//	bool iEnergy;
		double energy_factor    = 1.;
		double frequency_factor = 1.;
	//	TODO
	//	units.GiveMeConversionFactor(units.current_energy, energy_factor, iEnergy);
	//	units.GiveMeConversionFactor(units.current_frequency_factor, frequency_factor, iEnergy);
		for(i=1;i<=material[k].surfaceKinetics.number_of_reactions;i++)
			material[k].surfaceKinetics.reaction[i].Setup(energy_factor, frequency_factor);
	}
}

void OpenSMOKE_CHEMKINInterpreter::SurfaceParsingReactionsAdditionalData()
{
	// Loop over materials

	for(int k=1;k<=nSurfaceMaterials;k++)
	{
		for(int i=1;i<=material[k].surfaceKinetics.number_of_reactions-1;i++)
			material[k].surfaceKinetics.additionalLineReactionIndex.Append(material[k].surfaceKinetics.lineReactionIndex[i+1]-material[k].surfaceKinetics.lineReactionIndex[i]-1);
		for(int i=material[k].surfaceKinetics.lineReactionIndex[material[k].surfaceKinetics.number_of_reactions]+1;i<=material[k].end_material_section;i++)
			if (StringFindSubString(surface_lines[i], "END") == true)
				material[k].surfaceKinetics.additionalLineReactionIndex.Append(i-material[k].surfaceKinetics.lineReactionIndex[material[k].surfaceKinetics.number_of_reactions]-1);

		for(int i=1;i<=material[k].surfaceKinetics.number_of_reactions;i++)
		{
			
			int j = material[k].surfaceKinetics.lineReactionIndex[i];

			vector<string> instructions;
			instructions.push_back("instructions");
	
			string dummy;
			stringstream parsed_string(surface_lines[j]);
			for(;;)
			{
				parsed_string >> dummy;
				if (parsed_string.fail())
					break;
				instructions.push_back(dummy);
			}

			SurfaceParsingReactions(k, i, 1, instructions);

			material[k].surfaceKinetics.reaction[i].ReactionString();

			// Checking stoichiometry
			material[k].surfaceKinetics.reaction[i].CheckStoichiometry(thermo, indices_species_thermo_database, indices_site_thermo_database, indices_bulk_thermo_database, elements.elements_in_list, elements.isotope_in_list);
		}

		// Reactions parsing
		for(int i=1;i<=material[k].surfaceKinetics.number_of_reactions;i++)
		{
			int j = material[k].surfaceKinetics.lineReactionIndex[i];

			for(int kk=1;kk<=material[k].surfaceKinetics.additionalLineReactionIndex[i];kk++)
			{
				vector<string> instructions;
				instructions.push_back("instructions");
	
				string dummy;
				stringstream parsed_string(surface_lines[j+kk]);
				for(;;)
				{
					parsed_string >> dummy;
					if (parsed_string.fail())
						break;
					instructions.push_back(dummy);
				}
				SurfaceParsingReactionsAdditional(k, i, j, instructions);
			}
		}
	}
}

void OpenSMOKE_CHEMKINInterpreter::SurfaceParsingReactions(const int k, const int iReaction, const int iLine, vector<string> instructions)
{
	int number_of_instructions = instructions.size()-1;
	{
		material[k].surfaceKinetics.reaction[iReaction].A		= atof(instructions[number_of_instructions-2].c_str());
		material[k].surfaceKinetics.reaction[iReaction].Beta	= atof(instructions[number_of_instructions-1].c_str());
		material[k].surfaceKinetics.reaction[iReaction].E		= atof(instructions[number_of_instructions].c_str());
	}

	int i;
	string whole_reaction;
	for(i=1;i<=number_of_instructions-3;i++)
		whole_reaction+= (instructions[i] + " ");

	AddSpaces(whole_reaction, '+');

	vector<string> stoichiometry;
	stoichiometry.push_back("stoichiometry");
		
	string dummy;
	string dummy_complete;
	stringstream parsed_string(whole_reaction);

	int index;
	for(;;)
	{
		parsed_string >> dummy;
		if (parsed_string.fail())
			break;
		stoichiometry.push_back(dummy);
		
		if (dummy == "#REVERSIBLE#")
		{
			material[k].surfaceKinetics.reaction[iReaction].iReversible = true;
			index = stoichiometry.size()-1;
		}	
		if (dummy == "#DIRECT#")
		{
			material[k].surfaceKinetics.reaction[iReaction].iReversible = false;
			index = stoichiometry.size()-1;
		}	
	}

	// Reactants
	{
		bool iWarning = false;
		for(i=1;i<=index-1;i++)
		{
			string name;
			double number;

			SeparateNumberFromString(stoichiometry[i], name, number);
			if (name != ""  && name != "+")			material[k].surfaceKinetics.reaction[iReaction].nameDirect.push_back(name);
			if (name != "+" && iWarning == false)	material[k].surfaceKinetics.reaction[iReaction].nuDirect.Append(number);

			if (name == "")		iWarning = true;
			else				iWarning = false;
		}

		if (material[k].surfaceKinetics.reaction[iReaction].nameDirect.size()-1 != material[k].surfaceKinetics.reaction[iReaction].nuDirect.Size())
			ErrorMessage(iLine, "Error in reaction stoichiometry (reactant side)");

		// Number of reactants
		material[k].surfaceKinetics.reaction[iReaction].nReactants = material[k].surfaceKinetics.reaction[iReaction].nameDirect.size()-1;

		// Recognize indices and phase
		material[k].surfaceKinetics.reaction[iReaction].phaseDirect.resize(0);
		material[k].surfaceKinetics.reaction[iReaction].phaseDirect.push_back('0');
		for(i=1;i<=material[k].surfaceKinetics.reaction[iReaction].nReactants;i++)
		{
			int index;
			char phase;
			material[k].surfaceKinetics.RecognizeSpecies(material[k].surfaceKinetics.reaction[iReaction].nameDirect[i], index, phase);
			material[k].surfaceKinetics.reaction[iReaction].indexDirect.Append(index);
			material[k].surfaceKinetics.reaction[iReaction].phaseDirect.push_back(phase);
		}
	}

	// Products
	{
		bool iWarning = false;
		for(i=index+1;i<=int(stoichiometry.size())-1;i++)
		{
			string name;
			double number;

			SeparateNumberFromString(stoichiometry[i], name, number);
			if (name != ""  && name != "+")			material[k].surfaceKinetics.reaction[iReaction].nameInverse.push_back(name);
			if (name != "+" && iWarning == false)	material[k].surfaceKinetics.reaction[iReaction].nuInverse.Append(number);

			if (name == "")		iWarning = true;
			else				iWarning = false;
		}

		// Additional checks
		if (material[k].surfaceKinetics.reaction[iReaction].nameInverse.size()-1 != material[k].surfaceKinetics.reaction[iReaction].nuInverse.Size())
			ErrorMessage(iLine, "Error in reaction stoichiometry (product side)");

		// Number of product species
		material[k].surfaceKinetics.reaction[iReaction].nProducts = material[k].surfaceKinetics.reaction[iReaction].nameInverse.size()-1;

		// Recognize indices and phase
		material[k].surfaceKinetics.reaction[iReaction].phaseInverse.resize(0);
		material[k].surfaceKinetics.reaction[iReaction].phaseInverse.push_back('0');
		for(i=1;i<=material[k].surfaceKinetics.reaction[iReaction].nProducts;i++)
		{
			int index;
			char phase;
			material[k].surfaceKinetics.RecognizeSpecies(material[k].surfaceKinetics.reaction[iReaction].nameInverse[i], index, phase);
			material[k].surfaceKinetics.reaction[iReaction].indexInverse.Append(index);
			material[k].surfaceKinetics.reaction[iReaction].phaseInverse.push_back(phase);
		}
	}

	// Stoichiometry
	material[k].surfaceKinetics.reaction[iReaction].CompactStoichiometry();
}

void OpenSMOKE_CHEMKINInterpreter::SurfaceSummaryOnFile()
{
	ofstream fOutput;
	openOutputFileAndControl(fOutput, "SurfaceChemistrySummary.out");
	fOutput.setf(ios::scientific);
	
	// Print on video
	for(int k=1;k<=nSurfaceMaterials;k++)
	{
		fOutput << "---------------------------------------------------------------------------------" << endl;
		fOutput << "                     " << material[k].name << "                                  " << endl;
		fOutput << "---------------------------------------------------------------------------------" << endl;
		fOutput << endl;

		fOutput << "  Site species: " << material[k].surfaceKinetics.site_species_list.size()-1 << endl;
		for(int j=0;j<=material[k].surfaceKinetics.site_species_list.size()-1;j++)
			fOutput << "                " << material[k].surfaceKinetics.site_species_list[j] << endl;
		fOutput << endl;

		fOutput << "  Bulk species: " << material[k].surfaceKinetics.bulk_species_list.size()-1 << endl;
		for(int j=0;j<=material[k].surfaceKinetics.bulk_species_list.size()-1;j++)
			fOutput << "                " << material[k].surfaceKinetics.bulk_species_list[j] << endl;
		fOutput << endl;

		for(int j=1;j<=material[k].nSurfaceSites;j++)
		{
			fOutput << "  SITE: " << material[k].site_name[j] << endl;
			fOutput << "        " << material[k].site_density[j] << " mol/cm2" << endl;
			for (int i=1;i<=material[k].site_occupancy[j].Size();i++)
			fOutput << "        " << material[k].site_species_names[j][i] << " ("
							      << material[k].site_occupancy[j][i] << ")" << endl;
		}
		fOutput << endl;

		for(int j=1;j<=material[k].nSurfaceBulks;j++)
		{
			fOutput << "  BULK: " << material[k].bulk_name[j] << endl;
			for (int i=1;i<=material[k].bulk_density[j].Size();i++)
			{
				fOutput << "        " << material[k].bulk_species_names[j][i] << " ";
				if (material[k].bulk_density[j][i] != -1.)
					fOutput << material[k].bulk_density[j][i] << " g/cm3";
				fOutput << endl;
			}
		}

		fOutput << endl;

		// Reactions
		for(int i=1;i<=material[k].surfaceKinetics.number_of_reactions;i++)
		{
			fOutput << "Reaction #" << i << " " << material[k].surfaceKinetics.reaction[i].A << " "
												<< material[k].surfaceKinetics.reaction[i].Beta << " "
												<< material[k].surfaceKinetics.reaction[i].E << endl;
			
			fOutput << "Reactants: " << endl;
			for(int j=1;j<=material[k].surfaceKinetics.reaction[i].nReactants;j++)
			{
				fOutput << material[k].surfaceKinetics.reaction[i].phaseDirect[j] << " ";
				if (material[k].surfaceKinetics.reaction[i].phaseDirect[j] == 'G')
					fOutput << material[k].surfaceKinetics.gas_species_list[material[k].surfaceKinetics.reaction[i].indexDirect[j]] << " ";
				else if (material[k].surfaceKinetics.reaction[i].phaseDirect[j] == 'S')
					fOutput << material[k].surfaceKinetics.site_species_list[material[k].surfaceKinetics.reaction[i].indexDirect[j]] << " ";
				else if (material[k].surfaceKinetics.reaction[i].phaseDirect[j] == 'B')
					fOutput << material[k].surfaceKinetics.bulk_species_list[material[k].surfaceKinetics.reaction[i].indexDirect[j]] << " ";
				fOutput << material[k].surfaceKinetics.reaction[i].nuDirect[j] << endl;
			}
			fOutput << "Products: " << endl;
			for(int j=1;j<=material[k].surfaceKinetics.reaction[i].nProducts;j++)
			{
				fOutput << material[k].surfaceKinetics.reaction[i].phaseInverse[j] << " ";
				if (material[k].surfaceKinetics.reaction[i].phaseInverse[j] == 'G')
					fOutput << material[k].surfaceKinetics.gas_species_list[material[k].surfaceKinetics.reaction[i].indexInverse[j]] << " ";
				else if (material[k].surfaceKinetics.reaction[i].phaseInverse[j] == 'S')
					fOutput << material[k].surfaceKinetics.site_species_list[material[k].surfaceKinetics.reaction[i].indexInverse[j]] << " ";
				else if (material[k].surfaceKinetics.reaction[i].phaseInverse[j] == 'B')
					fOutput << material[k].surfaceKinetics.bulk_species_list[material[k].surfaceKinetics.reaction[i].indexInverse[j]] << " ";
				fOutput << material[k].surfaceKinetics.reaction[i].nuInverse[j] << endl;
			}
			fOutput << endl;
		}
		

	}	

	fOutput.close();
}

void OpenSMOKE_CHEMKINInterpreter::PrintSurfaceBinaryFile(const string file_name)
{
	string dummy;
	char name[Constants::NAME_SIZE];

	// Write On Binary File
	BzzSave binaryFile;
	binaryFile('*', file_name);
	getchar();
	BzzSave asciiFile;
	string file_name_ascii = file_name+".ascii";
	asciiFile(file_name_ascii);
	

	// Version
	dummy = "Surface110514";
	strcpy(name, dummy.c_str());
	binaryFile.fileSave.write((char*) name, sizeof(name));
	asciiFile << name;

	// Number of species
	binaryFile << int(site_species_global_list.size())-1;
	binaryFile << int(bulk_species_global_list.size())-1;
	asciiFile << int(site_species_global_list.size())-1;
	asciiFile << int(bulk_species_global_list.size())-1;

	// List of site species
	for(int k=1;k<=site_species_global_list.size()-1;k++)
	{	
		strcpy(name, site_species_global_list[k].c_str());
		binaryFile.fileSave.write((char*) name, sizeof(name));
		asciiFile << site_species_global_list[k];
	}

	// List of Bulk Species
	for(int k=1;k<=bulk_species_global_list.size()-1;k++)
	{
		strcpy(name, bulk_species_global_list[k].c_str());
		binaryFile.fileSave.write((char*) name, sizeof(name));
		asciiFile << bulk_species_global_list[k];
	}

	// Materials and sites
	binaryFile << nSurfaceMaterials;
	asciiFile << nSurfaceMaterials;
	for(int k=1;k<=nSurfaceMaterials;k++)
	{	
		// Material name
		strcpy(name, material[k].name.c_str());
		binaryFile.fileSave.write((char*) name, sizeof(name));
		asciiFile << material[k].name;

		// Number of sites
		binaryFile << material[k].nSurfaceSites;
		asciiFile << material[k].nSurfaceSites;

		// Site data
		for(int j=1;j<=material[k].nSurfaceSites;j++)
		{
			strcpy(name, material[k].site_name[j].c_str());
			binaryFile.fileSave.write((char*) name, sizeof(name));
			asciiFile << material[k].site_name[j];
			
			binaryFile << material[k].site_density[j];	// [mol/cm2]
			asciiFile << material[k].site_density[j];
		}

		// Site occupancies
		binaryFile << material[k].site_occupancy_matrix;
		asciiFile << material[k].site_occupancy_matrix;
		
		// Number of sites
		binaryFile << material[k].nSurfaceBulks;
		asciiFile << material[k].nSurfaceBulks;
		if (material[k].nSurfaceBulks > 0)
		{
			// Bulk data
			for(int j=1;j<=material[k].nSurfaceBulks;j++)
			{
				strcpy(name, material[k].bulk_name[j].c_str());
				binaryFile.fileSave.write((char*) name, sizeof(name));
				asciiFile << material[k].bulk_name[j];
			}
			// Bulk densities
			binaryFile << material[k].bulk_density_matrix;
			asciiFile << material[k].bulk_density_matrix;
		}
	}

	// Write thermodynamic
	thermo.WriteSurfaceThermodynamicData(binaryFile, asciiFile, indices_site_thermo_database, indices_bulk_thermo_database);

	// Elements
	thermo.WriteSurfaceElementsData(binaryFile, asciiFile);

	// Write kinetics
	for (int k=1;k<=nSurfaceMaterials;k++)
		material[k].surfaceKinetics.PrintBinaryFile(binaryFile, asciiFile, thermo, *this, material[k]);

	binaryFile.End();
	asciiFile.End();
}