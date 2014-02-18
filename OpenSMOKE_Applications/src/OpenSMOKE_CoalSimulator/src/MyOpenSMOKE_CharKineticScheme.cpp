/***************************************************************************
 *   Copyright (C) 2009 by Tiziano Maffei and Alberto Cuoci  	           *
 *   tiziano.maffei@chem.polimi.it   				                       *
 *   alberto.cuoci@polimi.it   						                       *
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
#include <iomanip>
#include "basic/OpenSMOKE_Utilities.h"
#include "basic/OpenSMOKE_Constants.h"
#include "MyOpenSMOKE_CharKineticScheme.h"
#include "OpenSMOKE_SolidExperiment.h"

void CompactStoichiometries(BzzVectorInt &index, BzzVector &nu, int &n, vector<string> &names);

void OpenSMOKE_CharKineticScheme::ErrorMessage(const string message)
{
    cout << endl;
    cout << "Class:  OpenSMOKE_CharKineticScheme"		<< endl;
    cout << "Object: " << name_object				<< endl;
    cout << "Error:  " << message					<< endl;
    cout << "Press enter to continue... "			<< endl;
    getchar();
    exit(-1);
}

void OpenSMOKE_CharKineticScheme::WarningMessage(const string message)
{
    cout << endl;
    cout << "Class:	  OpenSMOKE_CharKineticScheme"	<< endl;
    cout << "Object:  " << name_object				<< endl;
    cout << "Warning: " << message					<< endl;
	cout << endl;
}

void OpenSMOKE_CharKineticScheme::ErrorMessage(const string message, const int iLine)
{
    cout << endl;
    cout << "Class:  OpenSMOKE_CharKineticScheme"		<< endl;
    cout << "Object: " << name_object				<< endl;
	cout << "Line:   " << iLine						<< endl;
    cout << "Error:  " << message					<< endl;
    cout << "Press enter to continue... "			<< endl;
    getchar();
    exit(-1);
}

void OpenSMOKE_CharKineticScheme::WarningMessage(const string message, const int iLine)
{
    cout << endl;
    cout << "Class:	  OpenSMOKE_CharKineticScheme"	<< endl;
    cout << "Object:  " << name_object				<< endl;
	cout << "Line:    " << iLine					<< endl;
    cout << "Warning: " << message					<< endl;
	cout << endl;
}

OpenSMOKE_CharKineticScheme::OpenSMOKE_CharKineticScheme()
{
	name_object = "[Name not assigned]";

	openOutputFileAndControl(fLog, "Kinetics.log");
	fLog.setf(ios::scientific);
}

void OpenSMOKE_CharKineticScheme::ReadKineticScheme(const string file_name)
{
	char comment[Constants::COMMENT_SIZE];
	ifstream fInput;

	name_kinetics_file = file_name;
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

	ParsingSections();
	SyntaxCorrections();
	CleaningLines();
	ParsingSpecies();
	ParsingReactions();
	MemoryAllocation();
}

void OpenSMOKE_CharKineticScheme::ParsingSections()
{
	startGas		= 0;	endGas			= 100000;
	startCarbon		= 0;	endCarbon		= 100000;
	startAdsorbed	= 0;	endBulk			= 100000;
	startBulk		= 0;	endAdsorbed		= 100000;
	startReactions	= 0;	endReactions	= 100000;
	startEnd		= 0;

	for(int i=1;i<=number_of_lines;i++)
	{
		string dummy;
		stringstream parsed_string(lines[i]);

		parsed_string >> dummy;
		
		if (caseInsCompare("GAS", dummy) == true)	
		{
			if (startGas != 0)	ErrorMessage("GAS section is declared more than once...");
			startGas = i;
		}

		if (caseInsCompare("CARBON", dummy) == true)	
		{
			if (startCarbon != 0)	ErrorMessage("CARBON section is declared more than once...");
			startCarbon = i;
		}

		if (caseInsCompare("ADSORBED", dummy) == true)	
		{
			if (startAdsorbed != 0)	ErrorMessage("ADSORBED section is declared more than once...");
			startAdsorbed = i;
		}

		if (caseInsCompare("BULK", dummy) == true)	
		{
			if (startBulk != 0)	ErrorMessage("BULK section is declared more than once...");
			startBulk = i;
		}

		if (caseInsCompare("REACTIONS", dummy) == true)	
		{
			if (startReactions != 0)	ErrorMessage("REACTION section is declared more than once...");
			startReactions = i;
		}

		if (caseInsCompare("END", dummy) == true)	
		{
			if (startEnd != 0)	ErrorMessage("END section is declared more than once...");
			startEnd = i;
		}
	}
	
	if (startCarbon > startGas			&& startCarbon < endGas)		endGas = startCarbon-1;
	if (startBulk > startGas			&& startBulk < endGas)			endGas = startBulk-1;
	if (startAdsorbed > startGas		&& startAdsorbed < endGas)		endGas = startAdsorbed-1;
	if (startReactions > startGas		&& startReactions < endGas)		endGas = startReactions-1;

	if (startGas > startCarbon			&& startGas < endCarbon)		endCarbon = startGas-1;
	if (startBulk > startCarbon			&& startBulk < endCarbon)		endCarbon = startBulk-1;
	if (startAdsorbed > startCarbon		&& startAdsorbed < endCarbon)	endCarbon = startAdsorbed-1;
	if (startReactions > startCarbon	&& startReactions < endCarbon)	endCarbon = startReactions-1;

	if (startGas > startBulk			&& startGas < endBulk)			endBulk = startGas-1;
	if (startCarbon > startBulk			&& startCarbon < endBulk)		endBulk = startCarbon-1;
	if (startAdsorbed > startBulk		&& startAdsorbed < endBulk)		endBulk = startAdsorbed-1;
	if (startReactions > startBulk		&& startReactions < endBulk)	endBulk = startReactions-1;

	if (startGas > startAdsorbed		&& startGas < endAdsorbed)			endAdsorbed = startGas-1;
	if (startCarbon > startAdsorbed		&& startCarbon < endAdsorbed)		endAdsorbed = startCarbon-1;
	if (startBulk > startAdsorbed		&& startBulk < endAdsorbed)			endAdsorbed = startBulk-1;
	if (startReactions > startAdsorbed	&& startReactions < endAdsorbed)	endAdsorbed = startReactions-1;

	endReactions = startEnd;

	fLog << startCarbon << " " << endCarbon << endl;
	fLog << startBulk << " " << endBulk << endl;
	fLog << startGas << " " << endGas << endl;
	fLog << startAdsorbed << " " << endAdsorbed << endl;
	fLog << startReactions << " " << endReactions << endl;
}

void OpenSMOKE_CharKineticScheme::SyntaxCorrections()
{
	for(int i=1;i<=number_of_lines;i++)
	{
		StringSubstitutionAll(lines[i], "<=>",  " #REVERSIBLE# ");
		StringSubstitutionAll(lines[i], "=>",   " #DIRECT# ");
		StringSubstitutionAll(lines[i], "=",    " #REVERSIBLE# ");
	}
}

void OpenSMOKE_CharKineticScheme::CleaningLines()
{
	for(int i=1;i<=number_of_lines;i++)
	{
		for(int j=0;j<int(lines[i].size());j++)
			if (lines[i].at(j) == '\t')	lines[i].at(j) = ' ';

		AddSpaces(lines[i], '/');
	}
}

void OpenSMOKE_CharKineticScheme::AddSpaces(string &line, const char symbol)
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

void OpenSMOKE_CharKineticScheme::ParsingSpecies()
{
	ParsingSpecies(startCarbon, endCarbon, speciesCarbon);
	ParsingSpecies(startBulk, endBulk, speciesBulk);
	ParsingSpecies(startGas, endGas, speciesGas);
	ParsingSpecies(startAdsorbed, endAdsorbed, speciesAdsorbed);

	nCarbon = speciesCarbon.size()-1;
	nBulk = speciesBulk.size()-1;
	nGas = speciesGas.size()-1;
	nAdsorbed = speciesAdsorbed.size()-1;

	int j;
	
	list_of_species.push_back("list of species");
	for(j=1;j<=nCarbon;j++)		list_of_species.push_back(speciesCarbon[j]);
	for(j=1;j<=nBulk;j++)		list_of_species.push_back(speciesBulk[j]);
	for(j=1;j<=nGas;j++)		list_of_species.push_back(speciesGas[j]);
	for(j=1;j<=nAdsorbed;j++)	list_of_species.push_back(speciesAdsorbed[j]);

	for(j=1;j<=int(list_of_species.size())-1;j++)
		for(int i=j+1;i<=int(list_of_species.size())-1;i++)
			if (list_of_species[j] == list_of_species[i])
				ErrorMessage("This species was declared more than once: " + list_of_species[i]);
}

void OpenSMOKE_CharKineticScheme::ParsingSpecies(const int start, const int end, vector<string> &species)
{
	vector<string> instructions;
	instructions.push_back("instructions");
	
	for(int i=start;i<=end;i++)
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
	
	ParsingSpecies(instructions, species);
}

void OpenSMOKE_CharKineticScheme::ParsingSpecies(vector<string> instructions, vector<string> &species)
{
	int number_of_instructions = instructions.size()-1;

	for(int i=1;i<=number_of_instructions;i++)
	{	
		species.push_back(instructions[i]);
		if (i>=number_of_instructions) break;
	}
}

void OpenSMOKE_CharKineticScheme::ParsingReactions()
{
	int i;

	for(i=1;i<=number_of_lines;i++)
		if (StringFindSubString(lines[i], "#REVERSIBLE#") || StringFindSubString(lines[i], "#DIRECT#"))
			lineReactionIndex.Append(i);
	number_of_reactions = lineReactionIndex.Size();
	reactions = new OpenSMOKE_CharKineticScheme_ReactionData[number_of_reactions+1];

	// Initialize Reactions
	for(i=1;i<=number_of_reactions;i++)
		reactions[i].Setup(i, lineReactionIndex[i], this);

	ParsingReactionsAdditionalData();
}

void OpenSMOKE_CharKineticScheme::ParsingReactionsAdditionalData()
{
	int i;

	for(i=1;i<=number_of_reactions-1;i++)
		additionalLineReactionIndex.Append(lineReactionIndex[i+1]-lineReactionIndex[i]-1);
	for(i=lineReactionIndex[number_of_reactions]+1;i<=number_of_lines;i++)
		if (StringFindSubString(lines[i], "END") == true)
			additionalLineReactionIndex.Append(i-lineReactionIndex[number_of_reactions]-1);

	for(i=1;i<=number_of_reactions;i++)
	{
		int j = lineReactionIndex[i];

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

		reactions[i].ReactionString();
		fLog << reactions[i].reaction_string_clean << endl;

		// Checking stoichiometry
		//kinetics.reaction[i].CheckStoichiometry(thermo, indices_species_thermo_database, elements.elements_in_list, elements.isotope_in_list);
	}
}


void OpenSMOKE_CharKineticScheme::ParsingReactions(const int iReaction, const int iLine, vector<string> instructions)
{
	int number_of_instructions = instructions.size()-1;
	{
		reactions[iReaction].A		= atof(instructions[number_of_instructions-2].c_str());
		reactions[iReaction].Beta	= atof(instructions[number_of_instructions-1].c_str());
		reactions[iReaction].E		= atof(instructions[number_of_instructions].c_str());
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
			reactions[iReaction].iReversible = true;
			index = stoichiometry.size()-1;
		}	
		if (dummy == "#DIRECT#")
		{
			reactions[iReaction].iReversible = false;
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
			if (name != ""  && name != "+" && name != "M")			reactions[iReaction].nameDirect.push_back(name);
			if (name != "+" && iWarning == false && name != "M")	reactions[iReaction].nuDirect.Append(number);

			if (name == "")		iWarning = true;
			else				iWarning = false;

			if (name == "M")	reactions[iReaction].iThirdBody = true;
		}

		if (reactions[iReaction].nameDirect.size()-1 != reactions[iReaction].nuDirect.Size())
			ErrorMessage("Error in reaction stoichiometry (reactant side)", iLine);
		reactions[iReaction].nReactants = reactions[iReaction].nameDirect.size()-1;	
	
		for(i=1;i<=reactions[iReaction].nReactants;i++)
			reactions[iReaction].indexDirect.Append(RecognizeSpecies(reactions[iReaction].nameDirect[i]));
	}

	// Products
	{
		bool iWarning = false;
		for(i=index+1;i<=int(stoichiometry.size())-1;i++)
		{
			string name;
			double number;

			SeparateNumberFromString(stoichiometry[i], name, number);
			if (name != ""  && name != "+" && name != "M")			reactions[iReaction].nameInverse.push_back(name);
			if (name != "+" && iWarning == false && name != "M")	reactions[iReaction].nuInverse.Append(number);

			if (name == "")		iWarning = true;
			else				iWarning = false;

			if (name == "M")	reactions[iReaction].iThirdBody = true;
		}

		if (reactions[iReaction].nameInverse.size()-1 != reactions[iReaction].nuInverse.Size())
			ErrorMessage("Error in reaction stoichiometry (product side)", iLine);

		reactions[iReaction].nProducts = reactions[iReaction].nameInverse.size()-1;

		for(i=1;i<=reactions[iReaction].nProducts;i++)
			reactions[iReaction].indexInverse.Append(RecognizeSpecies(reactions[iReaction].nameInverse[i]));
	}

	reactions[iReaction].CompactStoichiometry();
}

int OpenSMOKE_CharKineticScheme::RecognizeSpecies(const string name_species)
{
	for(int i=1;i<=int(list_of_species.size())-1;i++)
		if (name_species == list_of_species[i])	return i;

	ErrorMessage("The " + name_species + " species was not declared!");
	return -1;
}

int OpenSMOKE_CharKineticScheme::RecognizeCarbonSpecies(const string name_species)
{
	for(int i=1;i<=int(speciesCarbon.size())-1;i++)
		if (name_species == speciesCarbon[i])	return i;

	ErrorMessage("The " + name_species + " species was not declared!");
	return -1;
}

void OpenSMOKE_CharKineticScheme::MemoryAllocation()
{
	int i;

	ChangeDimensions(nAdsorbed, &eta_adsorbed);
	ChangeDimensions(nGas,		&alfa_gas);
	ChangeDimensions(number_of_reactions, &A);
	ChangeDimensions(number_of_reactions, &Beta);
	ChangeDimensions(number_of_reactions, &E);

	ChangeDimensions(nAdsorbed,	number_of_reactions,	&lambdaAdsorbed);
	ChangeDimensions(nCarbon,	number_of_reactions,	&lambdaCarbon);
	ChangeDimensions(nGas,		number_of_reactions,	&lambdaGas);
	ChangeDimensions(nBulk,		number_of_reactions,	&lambdaBulk);

	ChangeDimensions(nAdsorbed,	number_of_reactions,	&nuAdsorbed);
	ChangeDimensions(nCarbon,	number_of_reactions,	&nuCarbon);
	ChangeDimensions(nGas,		number_of_reactions,	&nuGas);
	ChangeDimensions(nBulk,		number_of_reactions,	&nuBulk);

	ChangeDimensions(number_of_reactions,	&kappa);
	ChangeDimensions(number_of_reactions,	&rr);
	ChangeDimensions(nAdsorbed,				&R);
	ChangeDimensions(nCarbon,				&Rcarbon);
	ChangeDimensions(nGas,					&Rgas);
	ChangeDimensions(nCarbon,				&eta_f);

	for(i=1;i<=nAdsorbed;i++)
		eta_adsorbed[i] = GetEta(speciesAdsorbed[i]);

	for(i=1;i<=nGas;i++)
		alfa_gas[i] = GetAlfa(speciesGas[i]);

	for(i=1;i<=nCarbon;i++)	// TODO
		eta_f[i] = 1.;

	for(i=1;i<=number_of_reactions;i++)
	{
		A[i] = reactions[i].A;
		Beta[i] = reactions[i].Beta;
		E[i] = reactions[i].E;

		int j;
		for(j=1;j<=reactions[i].indexDirect.Size();j++)
		{
			if (reactions[i].indexDirect[j] >= 1 && reactions[i].indexDirect[j]<1+nCarbon)
			{
				int index = reactions[i].indexDirect[j]; 
				nuCarbon[index][i]		= -reactions[i].nuDirect[j];
				lambdaCarbon[index][i]	= reactions[i].nuDirect[j];
			}
			else if (reactions[i].indexDirect[j] >= 1+nCarbon && reactions[i].indexDirect[j]<1+nCarbon+nBulk)
			{
				int index = reactions[i].indexDirect[j]-nCarbon;
				nuBulk[index][i]	 = -reactions[i].nuDirect[j];
				lambdaBulk[index][i] = reactions[i].nuDirect[j];
			}
			else if (reactions[i].indexDirect[j] >= 1+nCarbon+nBulk && reactions[i].indexDirect[j]<1+nCarbon+nBulk+nGas)
			{
				int index = reactions[i].indexDirect[j]-nCarbon-nBulk;
				nuGas[index][i]		= -reactions[i].nuDirect[j];
				lambdaGas[index][i] = reactions[i].nuDirect[j];
			}
			else if (reactions[i].indexDirect[j] >= 1+nCarbon+nBulk+nGas && reactions[i].indexDirect[j]<1+nCarbon+nBulk+nGas+nGas)
			{
				int index = reactions[i].indexDirect[j]-nCarbon-nBulk-nGas;
				nuAdsorbed[index][i]		= -reactions[i].nuDirect[j];
				lambdaAdsorbed[index][i]	= reactions[i].nuDirect[j];
			}
		}

		for(j=1;j<=reactions[i].indexInverse.Size();j++)
		{
			if (reactions[i].indexInverse[j] >= 1 && reactions[i].indexInverse[j]<1+nCarbon)
			{
				int index = reactions[i].indexInverse[j]; 
				nuCarbon[index][i] += reactions[i].nuInverse[j];
			}
			else if (reactions[i].indexInverse[j] >= 1+nCarbon && reactions[i].indexInverse[j]<1+nCarbon+nBulk)
			{
				int index = reactions[i].indexInverse[j]-nCarbon;
				nuBulk[index][i] += reactions[i].nuInverse[j];
			}
			else if (reactions[i].indexInverse[j] >= 1+nCarbon+nBulk && reactions[i].indexInverse[j]<1+nCarbon+nBulk+nGas)
			{
				int index = reactions[i].indexInverse[j]-nCarbon-nBulk;
				nuGas[index][i] += reactions[i].nuInverse[j];
			}
			else if (reactions[i].indexInverse[j] >= 1+nCarbon+nBulk+nGas && reactions[i].indexInverse[j]<1+nCarbon+nBulk+nGas+nGas)
			{
				int index = reactions[i].indexInverse[j]-nCarbon-nBulk-nGas;
				nuAdsorbed[index][i] += reactions[i].nuInverse[j];
			}
		}
	}

	BzzVector lambdaSurface(number_of_reactions);
	BzzVector lambdaVolumetric(number_of_reactions);
	for(i=1;i<=number_of_reactions;i++)
	{
		lambdaSurface[i]	= 0.;
		lambdaVolumetric[i]	= 0.;
		for(int j=1;j<=nCarbon;j++)		lambdaSurface[i]	+= lambdaCarbon[j][i];
		for(int j=1;j<=nAdsorbed;j++)	lambdaSurface[i]	+= lambdaAdsorbed[j][i];
		for(int j=1;j<=nGas;j++)		lambdaVolumetric[i] += lambdaGas[j][i];
	}

	ChangeDimensions(number_of_reactions, &moleDimension);
	ChangeDimensions(number_of_reactions, &lengthDimension);
	ChangeDimensions(number_of_reactions, &conversion_from_kmol_m_s__to__mol_m_s);
	ChangeDimensions(number_of_reactions, &conversion_from_kmol_m_s__to__mol_cm_s);

	for(i=1;i<=number_of_reactions;i++)
	{
		moleDimension[i]	= 1.-lambdaSurface[i]-lambdaVolumetric[i];
		lengthDimension[i]	= 2.-2.*lambdaSurface[i]-3.*lambdaVolumetric[i];
		conversion_from_kmol_m_s__to__mol_m_s[i]  = 1. /  pow(1.e-3, moleDimension[i]);
		conversion_from_kmol_m_s__to__mol_cm_s[i] = 1. / (pow(1.e-3, moleDimension[i])/pow(1.e-2, lengthDimension[i]));
	}

	Summary(fLog);


	//TODO
	//nuAdsorbed.BzzPrint("nuAds");
	//nuGas.BzzPrint("nuGas");
	//nuCarbon.BzzPrint("nuCarb");
	//nuBulk.BzzPrint("nuBulk");
	//lambdaAdsorbed.BzzPrint("lambdaAds");
	//lambdaGas.BzzPrint("lambdaGas");
	//lambdaCarbon.BzzPrint("lambdaCarb");
	//lambdaBulk.BzzPrint("lambdaBulk");
	//eta_adsorbed.BzzPrint("eta");
}

string OpenSMOKE_CharKineticScheme::DimensionString(const int i)
{
	string dimension;
	stringstream mole;
	stringstream length;
	
	mole << -moleDimension[i];
	length << -lengthDimension[i];
	
	if (moleDimension[i] == 0. && lengthDimension[i] == 0.)
		dimension = "1/s";
	else
		dimension = "m^" + length.str() + " / kmol^" + mole.str() + " / s"; 

	return dimension;
}

void OpenSMOKE_CharKineticScheme::ReadDatabase(const string file_name)
{
	string dummy;
	ifstream fDatabase;
	openInputFileAndControl(fDatabase, file_name);

	database_species.push_back("species");
	database_species_eta.push_back(0.);
	database_species_alfa.push_back(0.);

	while(!fDatabase.eof())
	{
		fDatabase >> dummy;
		database_species.push_back(dummy);

		fDatabase >> dummy;
		if (dummy != "Eta")	ErrorMessage("Expected Eta key word - Found: " + dummy);
		fDatabase >> dummy;
		database_species_eta.push_back(atof(dummy.c_str()));

		fDatabase >> dummy;
		if (dummy != "Alfa")	ErrorMessage("Expected Alfa key word - Found: " + dummy);
		fDatabase >> dummy;
		database_species_alfa.push_back(atof(dummy.c_str()));
	}

	fDatabase.close();
}

double OpenSMOKE_CharKineticScheme::GetEta(const string species_name)
{
	for(int j=1;j<=int(database_species.size())-1;j++)
		if (database_species[j] == species_name)
			return database_species_eta[j];
	ErrorMessage("This species was not declared in the database: " + species_name);
	return -1;
}

double OpenSMOKE_CharKineticScheme::GetAlfa(const string species_name)
{
	for(int j=1;j<=int(database_species.size())-1;j++)
		if (database_species[j] == species_name)
			return database_species_alfa[j];
	ErrorMessage("This species was not declared in the database: " + species_name);
	return -1;
}

void OpenSMOKE_CharKineticScheme::Summary(ofstream &fOutput)
{
	int i;

	fOutput << "--------------------------------------------------------------------------------------" << endl;
	fOutput << "                                 Kinetic Scheme Summary                               " << endl;
	fOutput << "--------------------------------------------------------------------------------------" << endl;
	fOutput << endl;
	fOutput << "Number of Reactions: " << number_of_reactions			<< endl;
	fOutput << "Number of Species:   " << nCarbon+nBulk+nGas+nAdsorbed	<< endl;
	fOutput << endl;
	fOutput << "Number of Gas Species:      " << nGas		<< endl;
	fOutput << "Number of Adsorbed Species: " << nAdsorbed	<< endl;
	fOutput << endl << endl;

	fOutput << "--------------------------------------------------------------------------------------" << endl;
	fOutput << "                         Carbon Stoichiometric coefficients                           " << endl;
	fOutput << "--------------------------------------------------------------------------------------" << endl;
	fOutput << setw(14) << "";
	for(i=1;i<=number_of_reactions;i++)
		fOutput << setw(8) << left << i;
	fOutput << endl;
	fOutput << "--------------------------------------------------------------------------------------" << endl;
	for(i=1;i<=nCarbon;i++)
	{
		fOutput << setw(14) << speciesCarbon[i];
		for(int j=1;j<=number_of_reactions;j++)
			fOutput << setw(8) << left << fixed << setprecision(3) << nuCarbon[i][j];
		fOutput << endl;
	}
	fOutput << endl << endl;

	fOutput << "--------------------------------------------------------------------------------------" << endl;
	fOutput << "                         Bulk Stoichiometric coefficients                             " << endl;
	fOutput << "--------------------------------------------------------------------------------------" << endl;
	fOutput << setw(14) << "";
	for(i=1;i<=number_of_reactions;i++)
		fOutput << setw(8) << left << i;
	fOutput << endl;
	fOutput << "--------------------------------------------------------------------------------------" << endl;
	for(i=1;i<=nBulk;i++)
	{
		fOutput << setw(14) << speciesBulk[i];	// TODO
		for(int j=1;j<=number_of_reactions;j++)
			fOutput << setw(8) << left << fixed << setprecision(3) << nuBulk[i][j];
		fOutput << endl;
	}
	fOutput << endl << endl;

	fOutput << "--------------------------------------------------------------------------------------" << endl;
	fOutput << "                         Gas Species Stoichiometric coefficients                      " << endl;
	fOutput << "--------------------------------------------------------------------------------------" << endl;
	fOutput << setw(14) << "";
	for(i=1;i<=number_of_reactions;i++)
		fOutput << setw(8) << left << i;
	fOutput << endl;
	fOutput << "--------------------------------------------------------------------------------------" << endl;
	for(i=1;i<=nGas;i++)
	{
		fOutput << setw(14) << speciesGas[i];
		for(int j=1;j<=number_of_reactions;j++)
			fOutput << setw(8) << left << fixed << setprecision(3) << nuGas[i][j];
		fOutput << endl;
	}
	fOutput << endl << endl;

	fOutput << "--------------------------------------------------------------------------------------" << endl;
	fOutput << "                       Adsorbed Species Stoichiometric coefficients                   " << endl;
	fOutput << "--------------------------------------------------------------------------------------" << endl;
	fOutput << setw(14) << "";
	for(i=1;i<=number_of_reactions;i++)
		fOutput << setw(8) << left << i;
	fOutput << endl;
	fOutput << "--------------------------------------------------------------------------------------" << endl;
	for(i=1;i<=nAdsorbed;i++)
	{
		fOutput << setw(14) << speciesAdsorbed[i];
		for(int j=1;j<=number_of_reactions;j++)
			fOutput << setw(8) << left << fixed << setprecision(3) << nuAdsorbed[i][j];
		fOutput << endl;
	}
	fOutput << endl << endl;

	fOutput << "--------------------------------------------------------------------------------------" << endl;
	fOutput << "                              Carbon Reaction Orders                                  " << endl;
	fOutput << "--------------------------------------------------------------------------------------" << endl;
	fOutput << setw(14) << "";
	for(i=1;i<=number_of_reactions;i++)
		fOutput << setw(8) << left << i;
	fOutput << endl;
	fOutput << "--------------------------------------------------------------------------------------" << endl;
	for(i=1;i<=nCarbon;i++)
	{
		fOutput << setw(14) << speciesCarbon[i];
		for(int j=1;j<=number_of_reactions;j++)
			fOutput << setw(8) << left << fixed << setprecision(3) << lambdaCarbon[i][j];
		fOutput << endl;
	}
	fOutput << endl << endl;

	fOutput << "--------------------------------------------------------------------------------------" << endl;
	fOutput << "                               Bulk Reaction Orders                                   " << endl;
	fOutput << "--------------------------------------------------------------------------------------" << endl;
	fOutput << setw(14) << "";
	for(i=1;i<=number_of_reactions;i++)
		fOutput << setw(8) << left << i;
	fOutput << endl;
	fOutput << "--------------------------------------------------------------------------------------" << endl;
	for(i=1;i<=nBulk;i++)
	{
		fOutput << setw(14) << speciesBulk[i];	// TODO
		for(int j=1;j<=number_of_reactions;j++)
			fOutput << setw(8) << left << fixed << setprecision(3) << lambdaBulk[i][j];
		fOutput << endl;
	}
	fOutput << endl << endl;

	fOutput << "--------------------------------------------------------------------------------------" << endl;
	fOutput << "                              Gas Species Reaction Orders                             " << endl;
	fOutput << "--------------------------------------------------------------------------------------" << endl;
	fOutput << setw(14) << "";
	for(i=1;i<=number_of_reactions;i++)
		fOutput << setw(8) << left << i;
	fOutput << endl;
	fOutput << "--------------------------------------------------------------------------------------" << endl;
	for(i=1;i<=nGas;i++)
	{
		fOutput << setw(14) << speciesGas[i];
		for(int j=1;j<=number_of_reactions;j++)
			fOutput << setw(8) << left << fixed << setprecision(3) << lambdaGas[i][j];
		fOutput << endl;
	}
	fOutput << endl << endl;

	fOutput << "--------------------------------------------------------------------------------------" << endl;
	fOutput << "                          Adsorbed Species Reaction Orders                            " << endl;
	fOutput << "--------------------------------------------------------------------------------------" << endl;
	fOutput << setw(14) << "";
	for(i=1;i<=number_of_reactions;i++)
		fOutput << setw(8) << left << i;
	fOutput << endl;
	fOutput << "--------------------------------------------------------------------------------------" << endl;
	for(i=1;i<=nAdsorbed;i++)
	{
		fOutput << setw(14) << speciesAdsorbed[i];
		for(int j=1;j<=number_of_reactions;j++)
			fOutput << setw(8) << left << fixed << setprecision(3) << lambdaAdsorbed[i][j];
		fOutput << endl;
	}
	fOutput << endl << endl;

	
	// Reactions
	fOutput << "--------------------------------------------------------------------------------------" << endl;
	fOutput << "                                   List of Reactions                                  " << endl;
	fOutput << endl;
	fOutput << "                            Units: [kmol, m, s] and [J/mol]                           " << endl;
	fOutput << "--------------------------------------------------------------------------------------" << endl;
	fOutput << endl << endl;
	for(i=1;i<=number_of_reactions;i++)	reactions[i].SummaryOnFile(fLog);

	fOutput << "--------------------------------------------------------------------------------------" << endl;
	fOutput << "                               Reaction kinetic parameters                            " << endl;
	fOutput << "--------------------------------------------------------------------------------------" << endl;
	fOutput << endl;
	
	for(i=1;i<=number_of_reactions;i++)
	{
		fOutput << "Reaction " << i << endl;
		fOutput << "   Frequency factor:  " << scientific << setprecision(5) << reactions[i].A << " [" << DimensionString(i) << "]" << endl;
		fOutput << "   Temp. exponent:    " << scientific << setprecision(5) << reactions[i].Beta << " " << "[-]"			<< endl;
		fOutput << "   Activation energy: " << scientific << setprecision(5) << reactions[i].E << " " << "[J/kmol]"			<< endl;
		fOutput << endl;
		fOutput << "   Frequency factor:  " << scientific << setprecision(5) << reactions[i].A*conversion_from_kmol_m_s__to__mol_m_s[i]  << " [mol, m, s]" << endl;
		fOutput << "                      " << scientific << setprecision(5) << reactions[i].A*conversion_from_kmol_m_s__to__mol_cm_s[i] << " [mol, cm, s]" << endl;
		fOutput << "   Activation energy: " << scientific << setprecision(5) << reactions[i].E/1000.  << " [J/mol]" << endl;
		fOutput << "                      " << scientific << setprecision(5) << reactions[i].E/1000./Constants::R_cal_mol << " [cal/mol]" << endl;
		fOutput << endl;
	}	
}

void OpenSMOKE_CharKineticScheme::UpdateKineticConstants(const double T)
{
	for(int j=1;j<=number_of_reactions;j++)
		kappa[j] = A[j]*pow(T, Beta[j])*exp(-E[j]/Constants::R_J_kmol/T);		// [kmol, m, s]
}

void OpenSMOKE_CharKineticScheme::UpdateFormationRates(BzzVector &Ccarbon, BzzVector &Cbulk, BzzVector &Cgas,
													   BzzVector &Cadsorbed)
{
	int i,j;

	// Reaction rates
	rr = kappa;	
	//for(j=1;j<=number_of_reactions;j++)	// Cf
	  // cout <<" j: " << kappa[j]<<endl;
   // getchar();

	// Carbon species
	for(j=1;j<=number_of_reactions;j++)	// Cf
		for(i=1;i<=nCarbon;i++)
			rr[j] *= pow(Ccarbon[i], lambdaCarbon[i][j]);				// [kmol/m2/s]

	//for(i=1;i<=nGas;i++)
	//cout <<" Ci: " << Cgas[i]<<endl;
	//for(i=1;i<=nCarbon;i++)
	//cout <<" Cf: " << Ccarbon[i]<<endl;

   //  for(j=1;j<=number_of_reactions;j++)
//	{		 cout <<" j: " << kappa[j]<<endl;
//	     cout <<" j: " << rr[j]<<endl;
//}
  //     getchar();
	// Bulk species 
	// for(j=1;j<=number_of_reactions;j++)
	//	for(i=1;i<=nBulk;i++)
	//		rr[j] *= pow(Cbulk[i], lambdaBulk[i][j]);					// [kmol/m2/s]
	

	// Gas species
	for(i=1;i<=nGas;i++)
	{
		if (Cgas[i]>0.)
			for(j=1;j<=number_of_reactions;j++)	// Gas
				if (lambdaGas[i][j]!=0.)	
				{
					if (lambdaGas[i][j] == 1.)	rr[j] *= Cgas[i];							// [kmol/m2/s]
					else						rr[j] *= pow(Cgas[i], lambdaGas[i][j]);		// [kmol/m2/s]
				}
	}

	// Adsorbed species
	for(i=1;i<=nAdsorbed;i++)
		for(j=1;j<=number_of_reactions;j++)
			if (lambdaAdsorbed[i][j]!=0.)	
			{
				if (lambdaAdsorbed[i][j] == 1.)	rr[j] *= Cadsorbed[i];		// [kmol/m2/s]
				else rr[j] *= pow(Cadsorbed[i], lambdaAdsorbed[i][j]);		// [kmol/m2/s]
			}

	// Formation rates											
	for(i=1;i<=nCarbon;i++)		
		Rcarbon[i] = Dot(rr, nuCarbon.GetRow(i));			// [kmol/m2/s]

	// Formation rates											
	for(i=1;i<=nGas;i++)		
		Rgas[i] = Dot(rr, nuGas.GetRow(i));					// [kmol/m2/s]

	// Formation rates											
	for(i=1;i<=nAdsorbed;i++)		
		R[i] = Dot(rr, nuAdsorbed.GetRow(i));				// [kmol/m2/s]

	// Char reaction rates	
	double Rgas_tot		= Dot(Rgas, alfa_gas);
	Rchar = Rgas_tot;

//	for(j=1;j<=alfa_gas.Size();j++)		
//		cout << "j " << j << " " << alfa_gas[j] << endl;
//	getchar();

//	for(j=1;j<=number_of_reactions;j++)
//	{	//	 cout <<" j: " << kappa[j]<<endl;
	 //    cout <<" j: " << rr[j]<<endl;
//}
     //  getchar();
}

void OpenSMOKE_CharKineticScheme::UpdateKineticParameters(const int iModel, solid_regression_parameters &parameters, BzzVector &b, BzzMatrixInt &user, const int nCases, OpenSMOKE_SolidExperiment *experiment)
{
	int j;

	if (parameters == SOLID_REGRESSION_ALL)
	{
		if (iModel == 0)
		{
			for(j=1;j<=number_of_reactions;j++)	A[j]	= b[j];
			for(j=1;j<=number_of_reactions;j++)	Beta[j] = b[number_of_reactions+j];
			for(j=1;j<=number_of_reactions;j++)	E[j]	= b[2*number_of_reactions+j];
		}
		else if (iModel == 1)
		{
			for(j=1;j<=number_of_reactions;j++)	A[j]	= exp(b[j]);
			for(j=1;j<=number_of_reactions;j++)	Beta[j] = b[number_of_reactions+j];
			for(j=1;j<=number_of_reactions;j++)	E[j]	= b[2*number_of_reactions+j]*Constants::R_J_kmol;
		}
		else if (iModel == 2)
		{
			ErrorMessage("Model 2 not yet implemented...");
		}
	}

	if (parameters == SOLID_REGRESSION_ALL_PHI)
	{
		if (iModel == 0)
		{
			for(j=1;j<=number_of_reactions;j++)	A[j]				= b[j];
			for(j=1;j<=number_of_reactions;j++)	Beta[j]				= b[number_of_reactions+j];
			for(j=1;j<=number_of_reactions;j++)	E[j]				= b[2*number_of_reactions+j];
			for(j=1;j<=nCases;j++)				experiment[j].PHI	= b[3*number_of_reactions+j];
		}
		else if (iModel == 1)
		{
			for(j=1;j<=number_of_reactions;j++)	A[j]				= exp(b[j]);
			for(j=1;j<=number_of_reactions;j++)	Beta[j]				= b[number_of_reactions+j];
			for(j=1;j<=number_of_reactions;j++)	E[j]				= b[2*number_of_reactions+j]*Constants::R_J_kmol;
			for(j=1;j<=nCases;j++)				experiment[j].PHI	= b[3*number_of_reactions+j];
		}
		else if (iModel == 2)
		{
			ErrorMessage("Model 2 not yet implemented...");
		}
	}
	
	else if (parameters == SOLID_REGRESSION_A)
	{
		if (iModel == 0)		for(j=1;j<=number_of_reactions;j++)	A[j]	= b[j];
		else if (iModel == 1)	for(j=1;j<=number_of_reactions;j++)	A[j]	= exp(b[j]);
		else if (iModel == 2)	ErrorMessage("Model 2 not yet implemented...");
	}
	
	else if (parameters == SOLID_REGRESSION_BETA)
	{
		if (iModel == 0)		for(j=1;j<=number_of_reactions;j++)	Beta[j] = b[j];
		else if (iModel == 1)	for(j=1;j<=number_of_reactions;j++)	Beta[j] = b[j];
		else if (iModel == 2)	ErrorMessage("Model 2 not yet implemented...");
	}
	
	else if (parameters == SOLID_REGRESSION_E)
	{
		if (iModel == 0)		for(j=1;j<=number_of_reactions;j++)	E[j]	= b[j];
		else if (iModel == 1)	for(j=1;j<=number_of_reactions;j++)	E[j]	= b[j]*Constants::R_J_kmol;
		else if (iModel == 2)	ErrorMessage("Model 2 not yet implemented...");
	}

	else if (parameters == SOLID_REGRESSION_A_E)
	{
		if (iModel == 0)
		{
			for(j=1;j<=number_of_reactions;j++)	A[j] = b[j];
			for(j=1;j<=number_of_reactions;j++)	E[j] = b[number_of_reactions+j];
		}
		else if (iModel == 1)
		{
			for(j=1;j<=number_of_reactions;j++)	A[j] = exp(b[j]);
			for(j=1;j<=number_of_reactions;j++)	E[j] = b[number_of_reactions+j]*Constants::R_J_kmol;
		}
		else if (iModel == 2)
		{
		//	for(j=1;j<=number_of_reactions;j++)	A[j] = exp(b[j]+C1*C2*b[j+number_of_reactions]);
		//	for(j=1;j<=number_of_reactions;j++)	E[j] = -C1*Constants::R_J_kmol*b[number_of_reactions+j];
		}
	}

	else if (parameters == SOLID_REGRESSION_A_E_PHI)
	{
		if (iModel == 0)
		{
			for(j=1;j<=number_of_reactions;j++)	A[j]				= b[j];
			for(j=1;j<=number_of_reactions;j++)	E[j]				= b[number_of_reactions+j];
			for(j=1;j<=nCases;j++)				experiment[j].PHI	= b[2*number_of_reactions+j];

		}
		else if (iModel == 1)
		{
			for(j=1;j<=number_of_reactions;j++)	A[j] = exp(b[j]);
			for(j=1;j<=number_of_reactions;j++)	E[j] = b[number_of_reactions+j]*Constants::R_J_kmol;
			for(j=1;j<=nCases;j++)				experiment[j].PHI	= b[2*number_of_reactions+j];
		}
		else if (iModel == 2)
		{
		//	for(j=1;j<=number_of_reactions;j++)	A[j] = exp(b[j]+C1*C2*b[j+number_of_reactions]);
		//	for(j=1;j<=number_of_reactions;j++)	E[j] = -C1*Constants::R_J_kmol*b[number_of_reactions+j];
		}
	}

	else if (parameters == SOLID_REGRESSION_A_BETA)
	{
		if (iModel == 0)
		{
			for(j=1;j<=number_of_reactions;j++)	A[j]	= b[j];
			for(j=1;j<=number_of_reactions;j++)	Beta[j] = b[number_of_reactions+j];
		}
		else if (iModel == 1)
		{
			for(j=1;j<=number_of_reactions;j++)	A[j]	= exp(b[j]);
			for(j=1;j<=number_of_reactions;j++)	Beta[j] = b[number_of_reactions+j];
		}
		else if (iModel == 2)
		{
			ErrorMessage("Model 2 not yet implemented...");
		}
	}
	
	else if (parameters == SOLID_REGRESSION_BETA_E)
	{
		if (iModel == 0)
		{
			for(j=1;j<=number_of_reactions;j++)	Beta[j] = b[j];
			for(j=1;j<=number_of_reactions;j++)	E[j]	= b[number_of_reactions+j];
		}
		else if (iModel == 1)
		{
			for(j=1;j<=number_of_reactions;j++)	Beta[j] = b[j];
			for(j=1;j<=number_of_reactions;j++)	E[j]	= b[number_of_reactions+j]*Constants::R_J_kmol;
		}
		else if (iModel == 2)
		{
			ErrorMessage("Model 2 not yet implemented...");
		}
	}

	else if (parameters == SOLID_REGRESSION_USERDEFINED)
	{
		for(j=1;j<=b.Size();j++)
		{
			int jReaction = user[j][1];

			if (iModel == 0)
			{
					 if (user[j][2] == 1)	A[jReaction]				= b[j];
				else if (user[j][2] == 2)	Beta[jReaction]				= b[j];
				else if (user[j][2] == 3)	E[jReaction]				= b[j];
				else if (user[j][2] == 4)	experiment[user[j][1]].PHI	= b[j];
				else ErrorMessage("Only A || E || Beta can be optimized...");
			}
			

			else if (iModel == 1)
			{
					 if (user[j][2] == 1)	A[jReaction]				= exp(b[j]);
				else if (user[j][2] == 2)	Beta[jReaction]				= b[j];
				else if (user[j][2] == 3)	E[jReaction]				= b[j]*Constants::R_J_kmol;
				else if (user[j][2] == 4)	experiment[user[j][1]].PHI	= b[j];
				else ErrorMessage("Only A || E || Beta can be optimized...");
			}
			
			else if (iModel == 2)
			{
				ErrorMessage("Model 2 not yet implemented...");
			}
		}
	}
}

void OpenSMOKE_CharKineticScheme::UpdateOptimizerParameters(const int iModel, solid_regression_parameters &parameters, BzzVector &b, BzzMatrixInt &user, const int nCases, OpenSMOKE_SolidExperiment *experiment)
{
	int j;

	// First guess values
	if (parameters == SOLID_REGRESSION_ALL)
	{
		if (iModel == 0)
		{
			for(j=1;j<=number_of_reactions;j++)	b[j]						= A[j];
			for(j=1;j<=number_of_reactions;j++)	b[number_of_reactions+j]	= Beta[j];		
			for(j=1;j<=number_of_reactions;j++)	b[2*number_of_reactions+j]	= E[j];	
		}
                                                
		else if (iModel == 1)
		{
			for(j=1;j<=number_of_reactions;j++)	b[j]						= log(A[j]);
			for(j=1;j<=number_of_reactions;j++)	b[number_of_reactions+j]	= Beta[j];
			for(j=1;j<=number_of_reactions;j++)	b[2*number_of_reactions+j]	= E[j]/Constants::R_J_kmol;
		}
		else if (iModel == 2)
		{
			ErrorMessage("Model 2 not yet implemented...");
		}
	}

	if (parameters == SOLID_REGRESSION_ALL_PHI)
	{
		if (iModel == 0)
		{
			for(j=1;j<=number_of_reactions;j++)	b[j]						= A[j];
			for(j=1;j<=number_of_reactions;j++)	b[number_of_reactions+j]	= Beta[j];		
			for(j=1;j<=number_of_reactions;j++)	b[2*number_of_reactions+j]	= E[j];	
			for(j=1;j<=nCases;j++)				b[3*number_of_reactions+j]	= experiment[j].PHI;	
		}                                   
		else if (iModel == 1)
		{
			for(j=1;j<=number_of_reactions;j++)	b[j]						= log(A[j]);
			for(j=1;j<=number_of_reactions;j++)	b[number_of_reactions+j]	= Beta[j];
			for(j=1;j<=number_of_reactions;j++)	b[2*number_of_reactions+j]	= E[j]/Constants::R_J_kmol;
			for(j=1;j<=nCases;j++)				b[3*number_of_reactions+j]	= experiment[j].PHI;
		}
		else if (iModel == 2)
		{
			ErrorMessage("Model 2 not yet implemented...");
		}
	}

	else if (parameters == SOLID_REGRESSION_A)
	{
		if (iModel == 0)		for(j=1;j<=number_of_reactions;j++)	b[j]	= A[j];
		else if (iModel == 1)	for(j=1;j<=number_of_reactions;j++)	b[j]	= log(A[j]);
		else if (iModel == 2)	ErrorMessage("Model 2 not yet implemented...");
	}

	else if (parameters == SOLID_REGRESSION_BETA)
	{
		if (iModel == 0)		for(j=1;j<=number_of_reactions;j++)	b[j]	= Beta[j];
		else if (iModel == 1)	for(j=1;j<=number_of_reactions;j++)	b[j]	= Beta[j];
		else if (iModel == 2)	ErrorMessage("Model 2 not yet implemented...");
	}

	else if (parameters == SOLID_REGRESSION_E)
	{
		if (iModel == 0)		for(j=1;j<=number_of_reactions;j++)	b[j]	= E[j];		
		else if (iModel == 1)	for(j=1;j<=number_of_reactions;j++)	b[j]	= E[j]/Constants::R_J_kmol;
		else if (iModel == 2)	ErrorMessage("Model 2 not yet implemented...");
	}

	else if (parameters == SOLID_REGRESSION_A_E)
	{
		if (iModel == 0)
		{
			for(j=1;j<=number_of_reactions;j++)	b[j]						= A[j];
			for(j=1;j<=number_of_reactions;j++)	b[number_of_reactions+j]	= E[j];		
		}
		else if (iModel == 1)
		{
			for(j=1;j<=number_of_reactions;j++)	b[j]						= log(A[j]);
			for(j=1;j<=number_of_reactions;j++)	b[number_of_reactions+j]	= E[j]/Constants::R_J_kmol;
		}
		else if (iModel == 2)
		{
		//	for(j=1;j<=number_of_reactions;j++)	b[j]						= log(A[j]) + C2*E[j]/Constants::R_J_kmol;
		//	for(j=1;j<=number_of_reactions;j++)	b[number_of_reactions+j]	= -E[j]/Constants::R_J_kmol/C1;
		}
	}

	else if (parameters == SOLID_REGRESSION_A_E_PHI)
	{
		if (iModel == 0)
		{
			for(j=1;j<=number_of_reactions;j++)	b[j]						= A[j];
			for(j=1;j<=number_of_reactions;j++)	b[number_of_reactions+j]	= E[j];		
			for(j=1;j<=nCases;j++)				b[2*number_of_reactions+j]	= experiment[j].PHI;
		}
		
		else if (iModel == 1)
		{
			for(j=1;j<=number_of_reactions;j++)	b[j]						= log(A[j]);
			for(j=1;j<=number_of_reactions;j++)	b[number_of_reactions+j]	= E[j]/Constants::R_J_kmol;
			for(j=1;j<=nCases;j++)				b[2*number_of_reactions+j]	= experiment[j].PHI;
		}

		else if (iModel == 2)
		{
			ErrorMessage("Model 2 not yet implemented...");
		}
	}

	else if (parameters == SOLID_REGRESSION_A_BETA)
	{
		if (iModel == 0)
		{
			for(j=1;j<=number_of_reactions;j++)	b[j]						= A[j];
			for(j=1;j<=number_of_reactions;j++)	b[number_of_reactions+j]	= Beta[j];		
		}
		else if (iModel == 1)
		{
			for(j=1;j<=number_of_reactions;j++)	b[j]						= log(A[j]);
			for(j=1;j<=number_of_reactions;j++)	b[number_of_reactions+j]	= Beta[j];
		}
		else if (iModel == 2)
		{
			ErrorMessage("Model 2 not yet implemented...");
		}
	}

	else if (parameters == SOLID_REGRESSION_BETA_E)
	{
		if (iModel == 0)
		{
			for(j=1;j<=number_of_reactions;j++)	b[j]						= Beta[j];		
			for(j=1;j<=number_of_reactions;j++)	b[number_of_reactions+j]	= E[j];		
		}
		else if (iModel == 1)
		{
			for(j=1;j<=number_of_reactions;j++)	b[j]						= Beta[j];
			for(j=1;j<=number_of_reactions;j++)	b[number_of_reactions+j]	= E[j]/Constants::R_J_kmol;
		}
		else if (iModel == 2)
		{
			ErrorMessage("Model 2 not yet implemented...");
		}
	}

	else if (parameters == SOLID_REGRESSION_USERDEFINED)
	{
		for(j=1;j<=b.Size();j++)
		{
			int jReaction = user[j][1];

			if (iModel == 0)
			{
					 if (user[j][2] == 1)	b[j]	= A[jReaction];
				else if (user[j][2] == 2)	b[j]	= Beta[jReaction];		
				else if (user[j][2] == 3)	b[j]	= E[jReaction];	
				else if (user[j][2] == 4)	b[j]	= experiment[user[j][1]].PHI;	

				else ErrorMessage("Only A || E || Beta can be optimized...");
			}


			else if (iModel == 1)
			{
					 if (user[j][2] == 1)	b[j]	= log(A[jReaction]);
				else if (user[j][2] == 2)	b[j]	= Beta[jReaction];		
				else if (user[j][2] == 3)	b[j]	= E[jReaction]/Constants::R_J_kmol;
				else if (user[j][2] == 4)	b[j]	= experiment[user[j][1]].PHI;	
				else ErrorMessage("Only A || E || Beta can be optimized...");
			}
			
			else if (iModel == 2)
			{
				ErrorMessage("Model 2 not yet implemented...");
			}
		}
	}
}



void OpenSMOKE_CharKineticScheme_ReactionData::ErrorMessage(const string message)
{
    cout << endl;
    cout << "Class:    OpenSMOKE_CharKineticScheme_ReactionData"	<< endl;
    cout << "Error:    " << message         << endl;
    cout << "Reaction: " << iReaction       << endl;
    cout << "Line:     " << iLine           << endl;    
	cout << "Press a key to continue... "   << endl;
    getchar();
    exit(-1);
}

void OpenSMOKE_CharKineticScheme_ReactionData::WarningMessage(const string message)
{
    cout << endl;
    cout << "Class:  OpenSMOKE_CharKineticScheme_ReactionData"	<< endl;
    cout << "Warning:  " << message         << endl;
    cout << "Reaction: " << iReaction       << endl;
    cout << "Line:     " << iLine           << endl;    
    cout << "Press a key to continue... "   << endl;
    getchar();
}

void OpenSMOKE_CharKineticScheme_ReactionData::Setup(const int _iReaction, const int _iLine, OpenSMOKE_CharKineticScheme *_ptKinetics)
{
	iReaction	= _iReaction;
	iLine		= _iLine;
	ptKinetics	= _ptKinetics;
}

OpenSMOKE_CharKineticScheme_ReactionData::OpenSMOKE_CharKineticScheme_ReactionData()
{
	A		= 0.;
	Beta	= 0.;
	E		= 0.;

	nReactants = 0;
	nameDirect.push_back("reactant names");

	nProducts = 0;
	nameInverse.push_back("product names");

	nGlobalReactants = 0;
	nameGlobalDirect.push_back("reactant names (global)");
	
	nGlobalProducts = 0;
	nameGlobalInverse.push_back("product names (global)");
	
	sumNu			= 0.;
	lambdaTotal		= 0.;

	iReversible		= false;
	iGlobal			= false;
	iThirdBody		= false;
	iReactionBis	= 0;
}

void OpenSMOKE_CharKineticScheme_ReactionData::CompactStoichiometry()
{
	CompactStoichiometries(indexDirect, nuDirect, nReactants, nameDirect);
	CompactStoichiometries(indexInverse, nuInverse, nProducts, nameInverse);
}

void CompactStoichiometries(BzzVectorInt &index, BzzVector &nu, int &n, vector<string> &names)
{
	int i;
	
	for(i=1;i<=index.Size();i++)
	{
		if (index[i] != 0)
		for(int j=1;j<=index.Size();j++)
			if (index[i] == index[j])
				if (i != j)
				{
					index[j]   = 0;
					nu[i]	  += nu[j];
					nu[j]	   = 0;
				}	
	}

	BzzVectorInt indices;
	for(i=1;i<=index.Size();i++)
		if (index[i] == 0)	indices.Append(i);

	if (indices.Size() != 0)
	{
		vector<string> aux_names; aux_names.push_back("names");
		for(i=1;i<=index.Size();i++)
			if (index[i] != 0)	aux_names.push_back(names[i]);

		names.clear();
		names = aux_names;

		index.DeleteElements(indices);
		nu.DeleteElements(indices);

		n = index.Size();
	}
}

void OpenSMOKE_CharKineticScheme_ReactionData::ReactionString()
{
	int i;
	reaction_string_clean = "";
	
	for(i=1;i<=nReactants;i++)
	{
		if (nuDirect[i]!=1.)	
		{
			stringstream number;
			number << nuDirect[i];
			reaction_string_clean += number.str();
		}
		reaction_string_clean += nameDirect[i];
		if (i<nReactants)		reaction_string_clean += "+";
	}
	if (iThirdBody==true)		reaction_string_clean += "+M";

	if (iReversible==true)		reaction_string_clean += "=";
	else						reaction_string_clean += "=>";

	for(i=1;i<=nProducts;i++)
	{
		if (nuInverse[i]!=1.)
		{
			stringstream number;
			number << nuInverse[i];
			reaction_string_clean += number.str();
		}
		reaction_string_clean += nameInverse[i];
		if (i<nProducts)			reaction_string_clean += "+";
	}
	if (iThirdBody==true)		reaction_string_clean += "+M";

	sumNu = 0.;
	for(i=1;i<=nProducts;i++)	sumNu += nuInverse[i];
	for(i=1;i<=nReactants;i++)	sumNu -= nuDirect[i];

	reaction_string_complete = reaction_string_clean;
	
	if(iReactionBis !=0)	
	{
		stringstream index;
		index << iReaction;
		reaction_string_complete  += "#";
		reaction_string_complete += index.str();
		reaction_string_complete += "bis";
		reaction_string_complete += "#";
	}
	{
		string dummy;
		stringstream index;
		if(iReactionBis != 0)	index << iReactionBis;
		else					index << iReaction;

		dummy  = "#";
		dummy += index.str();
		dummy += "#";

		reaction_string_complete.insert(0,dummy);
	}
}

void OpenSMOKE_CharKineticScheme_ReactionData::SummaryOnFile(ofstream &fOutput)
{
	fOutput << setw(7) << right << iReaction;
	fOutput << ". ";
	fOutput << reaction_string_clean << endl;

	// Conventional reactions
	fOutput << setw(9)  << " ";
	fOutput << setw(9)  << left << "k:";
	fOutput << scientific	<< setprecision(6) << right << A << "\t";
	fOutput << setw(8) << setprecision(2) << fixed << right << Beta;
	fOutput << setw(14) << setprecision(2) << fixed << right << E/1000. << endl;

	fOutput << endl << endl;
}

