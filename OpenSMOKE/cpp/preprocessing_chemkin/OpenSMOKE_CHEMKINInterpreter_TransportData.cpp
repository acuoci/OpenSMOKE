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
#include "basic/OpenSMOKE_Constants.h"
#include "preprocessing_chemkin/OpenSMOKE_CHEMKINInterpreter_TransportData.h"
#include "preprocessing_chemkin/OpenSMOKE_CHEMKINInterpreter_TransportSpecies.h"
#include <iomanip>

void OpenSMOKE_CHEMKINInterpreter_TransportData::ErrorMessage(const int iLine, const std::string message)
{
    cout << endl;
    cout << "Class:  OpenSMOKE_CHEMKINInterpreter_TransportData"	<< endl;
    cout << "Object: " << name_object		<< endl;
    cout << "Line:   " << iLine             << endl;
    cout << "Error:  " << message           << endl;
    cout << "Press a key to continue... "   << endl;
    getchar();
    exit(-1);
}

void OpenSMOKE_CHEMKINInterpreter_TransportData::ErrorMessage(const std::string message)
{
    cout << endl;
    cout << "Class:  OpenSMOKE_CHEMKINInterpreter_TransportData"	<< endl;
    cout << "Object: " << name_object		<< endl;
    cout << "Error:  " << message           << endl;
    cout << "Press a key to continue... "   << endl;
    getchar();
    exit(-1);
}

void OpenSMOKE_CHEMKINInterpreter_TransportData::WarningMessage(const std::string message)
{
    cout << endl;
    cout << "Class:  OpenSMOKE_CHEMKINInterpreter_ThermoData"	<< endl;
    cout << "Object: "		<< name_object	<< endl;
    cout << "Warning:  "	<< message      << endl;
    cout << "Press a key to continue... "   << endl;
    getchar();
}

OpenSMOKE_CHEMKINInterpreter_TransportData::OpenSMOKE_CHEMKINInterpreter_TransportData()
{
	isActivated = false;
}

OpenSMOKE_CHEMKINInterpreter_TransportData::~OpenSMOKE_CHEMKINInterpreter_TransportData()
{

}

bool OpenSMOKE_CHEMKINInterpreter_TransportData::IsActivated()
{
	return isActivated;
}


void OpenSMOKE_CHEMKINInterpreter_TransportData::ReadTransportData(const std::string file_name, ofstream *_fLog)
{
	const int SIZE = 2000;
	char comment[SIZE];

	isActivated = true;
	fLog = _fLog;

	ifstream fInput;
	openInputFileAndControl(fInput, file_name);
	
	// ---------------------------------------------------------------
	// Reading lines
	// ---------------------------------------------------------------
	lines.push_back("List of lines");
	
	int count = 1;
	while(!fInput.eof())
	{
		fInput.getline(comment, SIZE);
		lines.push_back(comment);
	}
	fInput.close();
	

	number_of_lines = lines.size()-1;
	

	// ---------------------------------------------------------------
	// Parsing lines
	// ---------------------------------------------------------------

	int i;
	for(i=1;i<=number_of_lines;i++)
	{
			 if (CheckForBlankLine(lines[i])			== true)	indexBlankLines.Append(i);
		else if (CheckForCommentLineFromStart(lines[i])	== true)	indexCommentLines.Append(i);
		else if (CheckForEndLine(lines[i])				== true)	indexCommentLines.Append(i);		
		else
		{
			CheckForCommentLine(lines[i]);
			indexLines.Append(i);
		}
	}
	

	total_number_of_species = indexLines.Size();

	*fLog << endl;
	*fLog << " ----------------------------------------------------------------" << endl;
	*fLog << "                 Transport Properties Database                   " << endl;
	*fLog << " ----------------------------------------------------------------" << endl;
	*fLog << "    Total number of full lines:    " << indexLines.Size()			<< endl;
	*fLog << "    Total number of blank lines:   " << indexBlankLines.Size()	<< endl;
	*fLog << "    Total number of comment lines: " << indexCommentLines.Size()	<< endl;
	*fLog << "    Total number of lines:         " << number_of_lines			<< endl;
	*fLog << "    Total number of species:       " << total_number_of_species	<< endl;
	*fLog << " ----------------------------------------------------------------" << endl;
	*fLog << endl;

	species = new OpenSMOKE_CHEMKINInterpreter_TransportSpecies[total_number_of_species+1];
	
	for(i=1;i<=total_number_of_species;i++)
	{
		species[i].ReadMainData(lines[indexLines[i]], indexLines[i]);
	}

	CheckForDuplicates();
}

BzzVectorInt OpenSMOKE_CHEMKINInterpreter_TransportData::GiveMeSpeciesIndices(vector<string> &list)
{
	BzzVectorInt indices(list.size()-1);
	for(int j=1;j<=int(list.size())-1;j++)
	{
		bool jFound = false;
		for(int i=1;i<=total_number_of_species;i++)
			if ( caseInsCompare(species[i].name, list[j]) == true)
			{
				jFound = true;
				indices[j] = i;
				break;
			}

		if (jFound == false)
		{
			//indices[j]=1;	//TODO
			//break;
			//ErrorMessage("The " + list[j] + " species is not available in the transport properties database");
			*fLog << list[j] << endl;
			WarningMessage("The " + list[j] + " species is not available in the transport properties database");
		}
	}

	return indices;
}

void OpenSMOKE_CHEMKINInterpreter_TransportData::CheckForDuplicates()
{
	int i;
	BzzVectorInt iDuplicateCouples;

	for(i=1;i<=total_number_of_species;i++)
		for(int j=i+1;j<=total_number_of_species;j++)
			if ( caseInsCompare(species[i].name, species[j].name)==true)
			{
				iDuplicateCouples.Append(i);
				iDuplicateCouples.Append(j);
			}

	if (iDuplicateCouples !=0)
	{
		*fLog << " ----------------------------------------------------------------" << endl;
		*fLog << "    WARNING: Duplicate species in transport properties database   " << endl;
		*fLog << " ----------------------------------------------------------------" << endl;
		
		int k=1;
		for(i=1;i<=int(iDuplicateCouples.Size()/2);i++)
		{
			int j;
			*fLog << endl;
			
			j = iDuplicateCouples[k];
			*fLog << "    (" << species[j].index_line << ")\t\t" << species[j].name << endl;
			k++;

			j = iDuplicateCouples[k];
			*fLog << "    (" << species[j].index_line << ")\t\t" << species[j].name << endl;
			k++;
		}

		*fLog << endl;
		*fLog << " ----------------------------------------------------------------" << endl;
	}
}


void OpenSMOKE_CHEMKINInterpreter_TransportData::SummaryOnFile(ofstream &fOutput, BzzVectorInt &indices)
{
	int i;

	fOutput << "------------------------------------------------------------------------------------------------------------------" << endl;
	fOutput << "  Species                  Shape       eps/kb          sigma             mu           alfa        zRot298         " << endl;
	fOutput << "------------------------------------------------------------------------------------------------------------------" << endl;
	for(i=1;i<=indices.Size();i++)
	{
		fOutput << right << setw(5) << i;
		fOutput << ". ";
		fOutput << setw(20) << left  << species[indices[i]].name;
		fOutput << setw(3)  << right << species[indices[i]].shape_factor;
		fOutput << setw(15) << fixed << right << setprecision(4)  << species[indices[i]].epsylon_over_kb;
		fOutput << setw(15) << fixed << right << setprecision(4)  << species[indices[i]].sigma;
		fOutput << setw(15) << fixed << right << setprecision(4)  << species[indices[i]].mu/1.e-6;
		fOutput << setw(15) << fixed << right << setprecision(4)  << species[indices[i]].alfa;
		fOutput << setw(15) << fixed << right << setprecision(4)  << species[indices[i]].zRot298;
		fOutput << endl;
	}
	fOutput << "------------------------------------------------------------------------";
	fOutput << endl << endl << endl;
}

void OpenSMOKE_CHEMKINInterpreter_TransportData::ProcessTransportData(BzzVectorInt &indices)
{
	int number_of_species = indices.Size();

	ChangeDimensions(number_of_species, &shape_factor);
	ChangeDimensions(number_of_species, &epsylon_over_kb);
	ChangeDimensions(number_of_species, &sigma);
	ChangeDimensions(number_of_species, &mu);
	ChangeDimensions(number_of_species, &alfa);
	ChangeDimensions(number_of_species, &zRot298);

	for(int k=1;k<=number_of_species;k++)
	{
		int j = indices[k];
		
		shape_factor[k]		=	species[j].shape_factor;
		epsylon_over_kb[k]	=	species[j].epsylon_over_kb;
		sigma[k]			=	species[j].sigma;
		mu[k]				=	species[j].mu;
		alfa[k]				=	species[j].alfa;
		zRot298[k]			=	species[j].zRot298;
	}
}

double OpenSMOKE_CHEMKINInterpreter_TransportData::ReducedDiameter(const int i1, const int i2)
{
	return 0.50*(sigma[i1]+sigma[i2]);	// [A]
}
