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

#include <sstream>
#include "basic/OpenSMOKE_Utilities.h"
#include "basic/OpenSMOKE_WarningFile.h"
#include "preprocessing_chemkin/OpenSMOKE_CHEMKINInterpreter.h"
#include "preprocessing_chemkin/OpenSMOKE_CHEMKINInterpreter_SurfaceKineticsData.h"
#include "preprocessing_chemkin/OpenSMOKE_CHEMKINInterpreter_TransportData.h"

void OpenSMOKE_CHEMKINInterpreter_SurfaceKineticsData::ErrorMessage(const std::string message)
{
    cout << endl;
    cout << "Class:  OpenSMOKE_CHEMKINInterpreter_SurfaceKineticsData"	<< endl;
    cout << "Object: " << name_object		<< endl;
    cout << "Error:  " << message           << endl;
    cout << "Press a key to continue... "   << endl;
    getchar();
    exit(-1);
}

void OpenSMOKE_CHEMKINInterpreter_SurfaceKineticsData::WarningMessage(const std::string message)
{
    cout << endl;
    cout << "Class:  OpenSMOKE_CHEMKINInterpreter_SurfaceKineticsData"	<< endl;
    cout << "Object: "		<< name_object	<< endl;
    cout << "Warning:  "	<< message      << endl;
    cout << "Press a key to continue... "   << endl;
    getchar();
}

OpenSMOKE_CHEMKINInterpreter_SurfaceKineticsData::~OpenSMOKE_CHEMKINInterpreter_SurfaceKineticsData()
{

}

OpenSMOKE_CHEMKINInterpreter_SurfaceKineticsData::OpenSMOKE_CHEMKINInterpreter_SurfaceKineticsData()
{

}

void OpenSMOKE_CHEMKINInterpreter_SurfaceKineticsData::Setup(ofstream *_fLog, OpenSMOKE_WarningFile *_fWarning)
{
	fLog		= _fLog;
	fWarning	= _fWarning; 
}

void OpenSMOKE_CHEMKINInterpreter_SurfaceKineticsData::RecognizeSpecies(const std::string name_species, int &index, char &phase)
{
	for(int i=1;i<=int(gas_species_list.size())-1;i++)
		if (name_species == gas_species_list[i])
		{
			index = i;
			phase = 'G';
			return;
		}

	for(int i=1;i<=int(site_species_list.size())-1;i++)
		if (name_species == site_species_list[i])
		{
			index = i;
			phase = 'S';
			return;
		}

	for(int i=1;i<=int(bulk_species_list.size())-1;i++)
		if (name_species == bulk_species_list[i])
		{
			index = i;
			phase = 'B';
			return;
		}

	ErrorMessage("The " + name_species + " species was not declared!");
}

void OpenSMOKE_CHEMKINInterpreter_SurfaceKineticsData::FinalChecks()
{
	for(int i=1;i<=number_of_reactions;i++)
	{
		reaction[i].Complete();
		reaction[i].Check();
	}
//	CheckDuplicateReactions();
	CheckReactionRatesReactions();
}

void OpenSMOKE_CHEMKINInterpreter_SurfaceKineticsData::CheckDuplicateReactions()
{
	ErrorMessage("CheckDuplicateReactions not yet available");
	/*
	int i;

	for(i=1;i<=number_of_reactions;i++)
	{	
		reaction[i].indexDirectOrdered	= reaction[i].indexDirect;
		reaction[i].indexInverseOrdered	= reaction[i].indexInverse;

		Sort(&reaction[i].indexDirectOrdered,  &reaction[i].indexDirectOrderedSequence);
		Sort(&reaction[i].indexInverseOrdered, &reaction[i].indexInverseOrderedSequence);
	}

	for(i=1;i<=number_of_reactions;i++)
		for(int j=i+1;j<=number_of_reactions;j++)
		{			
			if (CheckDuplicateReactions(i,j)==true)
			{
				aDuplicate.Append(i);
				bDuplicate.Append(j);
			}
		}

	*fLog << " ----------------------------------------------------------------" << endl;
	*fLog << "                 Duplicate Reactions (Couples)                   " << endl;
	*fLog << " ----------------------------------------------------------------" << endl;
	for(i=1;i<=aDuplicate.Size();i++)
		*fLog << i << "\t" << aDuplicate[i] << "\t" << bDuplicate[i] << endl;
	
	for(i=1;i<=aDuplicate.Size();i++)
	{
		if (reaction[aDuplicate[i]].iDuplicate == false && reaction[bDuplicate[i]].iDuplicate == false)
		{
			stringstream message;
			message << "Reaction " << aDuplicate[i] << " (line " << reaction[aDuplicate[i]].iLine << ") and Reaction " << bDuplicate[i] << " (line " << reaction[bDuplicate[i]].iLine << ") must be declared as DUPLICATE reactions"; 
			ErrorMessage(message.str());
		}
		else if (reaction[aDuplicate[i]].iDuplicate == false && reaction[bDuplicate[i]].iDuplicate == true)
		{
			stringstream message;
			message << "Reaction " << aDuplicate[i] << " (line " << reaction[aDuplicate[i]].iLine << ") is a DUPLICATE of Reaction " << bDuplicate[i] << " (line " << reaction[bDuplicate[i]].iLine << ")" ; 
			ErrorMessage(message.str());
		}
		else if (reaction[aDuplicate[i]].iDuplicate == true  && reaction[bDuplicate[i]].iDuplicate == false)
		{
			stringstream message;
			message << "Reaction " << bDuplicate[i] << " (line " << reaction[bDuplicate[i]].iLine << ") is a DUPLICATE of Reaction " << aDuplicate[i] << " (line " << reaction[aDuplicate[i]].iLine << ")" ; 
			ErrorMessage(message.str());
		}
	}

	for(i=1;i<=number_of_reactions;i++)
	{
		if (reaction[i].iDuplicate == true)
		{
			int k;
			bool kFound = false;
			for(k=1;k<=aDuplicate.Size();k++)
				if (aDuplicate[k] == i)	{kFound = true; break;}
			for(k=1;k<=bDuplicate.Size();k++)
				if (bDuplicate[k] == i)	{kFound = true; break;}
			if (kFound == false)
			{
				stringstream message;
				message << "Reaction " << i << " (line " << reaction[i].iLine << ") declared DUPLICATE, however, no duplicate for this reaction was found";
				ErrorMessage(message.str());
			}
		}
	}
	*/
}

void OpenSMOKE_CHEMKINInterpreter_SurfaceKineticsData::CheckReactionRatesReactions()
{
	int i;
	BzzVectorInt indices;
	for(i=1;i<=number_of_reactions;i++)
	{
		if(reaction[i].iReverseRate == true)
		{
			indices.Append(i);
			reaction[i].iReversible = false;
		}
	}

	for(i=1;i<=indices.Size();i++)
		reaction[number_of_reactions+i].ReverseRateFromReaction(number_of_reactions+i, reaction[indices[i]]);

	number_of_reactions += indices.Size();
}

void OpenSMOKE_CHEMKINInterpreter_SurfaceKineticsData::Statistics()
{
	int i;
	for(i=1;i<=number_of_reactions;i++)
	{
		if (reaction[i].iReversible	 == true)	iReversible.Append(i);
		if (reaction[i].iDuplicate == true)		iDuplicate.Append(i);
		if (reaction[i].iGlobal == true)		iGlobal.Append(i);
	}
}

void OpenSMOKE_CHEMKINInterpreter_SurfaceKineticsData::SummaryOnFile(ofstream &fOutput)
{
	fOutput << "-----------------------------------------------------------------"	<< endl;
	fOutput << "                    SURFACE CHEMICAL REACTIONS                   "	<< endl;
	fOutput << "-----------------------------------------------------------------"	<< endl;
	fOutput << endl;
	fOutput << "  Number of Reactions:                                  " << number_of_reactions << endl;
	fOutput << "  Number of Reversible Reactions:                       " << iReversible.Size() << endl;
	fOutput << "  Number of Irreversible Reactions:                     " << number_of_reactions-iReversible.Size() << endl;
	fOutput << "  Number of Duplicate Reactions:                        " << iDuplicate.Size() << endl;
	fOutput << "  Number of Global Reactions:                           " << iGlobal.Size() << endl;
	fOutput << endl;
	fOutput << "-----------------------------------------------------------------"	<< endl;
	fOutput << endl << endl;
}

void OpenSMOKE_CHEMKINInterpreter_SurfaceKineticsData::PrintBinaryFile(BzzSave &outputFile, BzzSave &asciiFile,OpenSMOKE_CHEMKINInterpreter_ThermoData &thermo, OpenSMOKE_CHEMKINInterpreter &interp, OpenSMOKE_PreProcessorSurfaceMaterial &material)
{
	// Writing number of species and number of reactions
	outputFile << number_of_reactions;
	asciiFile << number_of_reactions;

	// Writing reaction data
	for(int i=1;i<=number_of_reactions;i++)
		reaction[i].PrintOnBinaryFile(outputFile, asciiFile, thermo, material);

	if (iGlobal.Size() != 0)
	{
		ErrorMessage("Global not yet available");
		/*
		outputFile << iGlobal;

		for(i=1;i<=iGlobal.Size();i++)
		{
			outputFile << reaction[iGlobal[i]].indexGlobalDirect;
			outputFile << reaction[iGlobal[i]].lambdaGlobalDirect;
			if (reaction[iGlobal[i]].indexGlobalInverse.Size() !=0)
			{
				outputFile << 1;
				outputFile << reaction[iGlobal[i]].indexGlobalInverse;
				outputFile << reaction[iGlobal[i]].lambdaGlobalInverse;
			}
			else outputFile << 0;
		}
		*/
	}
}