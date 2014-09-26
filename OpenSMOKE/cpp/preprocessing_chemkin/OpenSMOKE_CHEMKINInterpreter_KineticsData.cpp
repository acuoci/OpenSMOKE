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
#include "preprocessing_chemkin/OpenSMOKE_CHEMKINInterpreter_KineticsData.h"
#include "preprocessing_chemkin/OpenSMOKE_CHEMKINInterpreter_TransportData.h"

void OpenSMOKE_CHEMKINInterpreter_KineticsData::ErrorMessage(const std::string message)
{
    cout << endl;
    cout << "Class:  OpenSMOKE_CHEMKINInterpreter_KineticsData"	<< endl;
    cout << "Object: " << name_object		<< endl;
    cout << "Error:  " << message           << endl;
    cout << "Press a key to continue... "   << endl;
    getchar();
    exit(-1);
}

void OpenSMOKE_CHEMKINInterpreter_KineticsData::WarningMessage(const std::string message)
{
    cout << endl;
    cout << "Class:  OpenSMOKE_CHEMKINInterpreter_KineticsData"	<< endl;
    cout << "Object: "		<< name_object	<< endl;
    cout << "Warning:  "	<< message      << endl;
    cout << "Press a key to continue... "   << endl;
    getchar();
}

OpenSMOKE_CHEMKINInterpreter_KineticsData::~OpenSMOKE_CHEMKINInterpreter_KineticsData()
{

}

OpenSMOKE_CHEMKINInterpreter_KineticsData::OpenSMOKE_CHEMKINInterpreter_KineticsData()
{

}

void OpenSMOKE_CHEMKINInterpreter_KineticsData::Setup(ofstream *_fLog, OpenSMOKE_WarningFile *_fWarning)
{
	fLog		= _fLog;
	fWarning	= _fWarning; 
}

int OpenSMOKE_CHEMKINInterpreter_KineticsData::RecognizeSpecies(const std::string name_species)
{
	for(int i=1;i<=int(species_list.size())-1;i++)
		if (name_species == species_list[i])	return i;

	ErrorMessage("The " + name_species + " species was not declared!");
	return -1;
}

void OpenSMOKE_CHEMKINInterpreter_KineticsData::PrintNamesFile(const std::string file_name)
{
	ofstream fOutput;
	openOutputFileAndControl(fOutput, file_name);
	fOutput.setf(ios::scientific);

	for(int i=1;i<=int(species_list.size())-1;i++)
		fOutput << species_list[i] << endl;
	fOutput << "END" << endl;

	fOutput.close();
}

void OpenSMOKE_CHEMKINInterpreter_KineticsData::CheckDuplicateReactions()
{
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
}

bool OpenSMOKE_CHEMKINInterpreter_KineticsData::CheckDuplicateReactions(const int i, const int j)
{
	// Third Body
	if (reaction[i].iThirdBody == true)
		if (reaction[j].iThirdBody == false && reaction[j].iHighPressure == false && reaction[j].iLowPressure == false && reaction[j].iChebishevPolynomials == false) return false;	
	if (reaction[j].iThirdBody == true)
		if (reaction[i].iThirdBody == false && reaction[i].iHighPressure == false && reaction[i].iLowPressure == false && reaction[i].iChebishevPolynomials == false) return false;
	if (reaction[j].iThirdBody == true && reaction[j].thirdbodysingle == "null" &&
		reaction[i].iThirdBody == true && reaction[i].thirdbodysingle != "null")	return false;
	if (reaction[j].iThirdBody == true && reaction[j].thirdbodysingle != "null" &&
		reaction[i].iThirdBody == true && reaction[i].thirdbodysingle == "null")	return false;


	// Fall-Off
	if (reaction[i].iLowPressure == true)
		if (reaction[j].iThirdBody == false && reaction[j].iHighPressure == false && reaction[j].iLowPressure == false && reaction[j].iChebishevPolynomials == false) return false;	
	if (reaction[j].iLowPressure == true)
		if (reaction[i].iThirdBody == false && reaction[i].iHighPressure == false && reaction[i].iLowPressure == false && reaction[i].iChebishevPolynomials == false) return false;

	// Chemically Activated Bimolecular Reactions
	if (reaction[i].iHighPressure == true)
		if (reaction[j].iThirdBody == false && reaction[j].iHighPressure == false && reaction[j].iLowPressure == false && reaction[j].iChebishevPolynomials == false) return false;	
	if (reaction[j].iHighPressure == true)
		if (reaction[i].iThirdBody == false && reaction[i].iHighPressure == false && reaction[i].iLowPressure == false && reaction[i].iChebishevPolynomials == false) return false;

	// Chebishev Reactions
	if (reaction[i].iChebishevPolynomials == true)
		if (reaction[j].iThirdBody == false && reaction[j].iHighPressure == false && reaction[j].iLowPressure == false && reaction[j].iChebishevPolynomials == false) return false;	
	if (reaction[j].iChebishevPolynomials == true)
		if (reaction[i].iThirdBody == false && reaction[i].iHighPressure == false && reaction[i].iLowPressure == false && reaction[i].iChebishevPolynomials == false) return false;

	if (reaction[i].iReversible == false && reaction[j].iReversible == false)
	{
		return CheckDuplicateReactionsGivenOrder(i,j);
	}
	else
	{
		bool iEqual = CheckDuplicateReactionsGivenOrder(i,j);
		if (iEqual == false)	iEqual = CheckDuplicateReactionsInverseOrder(i,j);
		return iEqual;
	}
}

bool OpenSMOKE_CHEMKINInterpreter_KineticsData::CheckDuplicateReactionsGivenOrder(const int i, const int j)
{
	if (reaction[i].indexDirect.Size()  != reaction[j].indexDirect.Size())	return false;
	if (reaction[i].indexInverse.Size() != reaction[j].indexInverse.Size())	return false;

	int k;
	for (k=1;k<=reaction[i].indexDirect.Size();k++)
		if (reaction[i].indexDirectOrdered[k] != reaction[j].indexDirectOrdered[k])		return false;
					
	for (k=1;k<=reaction[i].indexInverse.Size();k++)
		if (reaction[i].indexInverseOrdered[k] != reaction[j].indexInverseOrdered[k])	return false;

	BzzVector nuDirectOrdered_i  = reaction[i].nuDirect;
	Reorder(&nuDirectOrdered_i, reaction[i].indexDirectOrderedSequence);
	BzzVector nuDirectOrdered_j  = reaction[j].nuDirect;
	Reorder(&nuDirectOrdered_j, reaction[j].indexDirectOrderedSequence);

	for (k=1;k<=reaction[i].indexDirect.Size();k++)
		if (nuDirectOrdered_i[k] != nuDirectOrdered_j[k])		return false;

	BzzVector nuInverseOrdered_i  = reaction[i].nuInverse;
	Reorder(&nuInverseOrdered_i, reaction[i].indexInverseOrderedSequence);
	BzzVector nuInverseOrdered_j  = reaction[j].nuInverse;
	Reorder(&nuInverseOrdered_j, reaction[j].indexInverseOrderedSequence);
	
	for (k=1;k<=reaction[i].indexInverse.Size();k++)
		if (nuInverseOrdered_i[k] != nuInverseOrdered_j[k])		return false;

	return true;
}

bool OpenSMOKE_CHEMKINInterpreter_KineticsData::CheckDuplicateReactionsInverseOrder(const int i, const int j)
{
	if (reaction[i].indexDirect.Size()  != reaction[j].indexInverse.Size())	return false;
	if (reaction[i].indexInverse.Size() != reaction[j].indexDirect.Size())	return false;

	int k;
	for (k=1;k<=reaction[i].indexDirect.Size();k++)
		if (reaction[i].indexDirectOrdered[k] != reaction[j].indexInverseOrdered[k])		return false;
					
	for (k=1;k<=reaction[i].indexInverse.Size();k++)
		if (reaction[i].indexInverseOrdered[k] != reaction[j].indexDirectOrdered[k])	return false;

	BzzVector nuDirectOrdered_i  = reaction[i].nuDirect;
	Reorder(&nuDirectOrdered_i, reaction[i].indexDirectOrderedSequence);
	BzzVector nuInverseOrdered_j  = reaction[j].nuInverse;
	Reorder(&nuInverseOrdered_j, reaction[j].indexInverseOrderedSequence);
	for (k=1;k<=reaction[i].indexDirect.Size();k++)
		if (nuDirectOrdered_i[k] != nuInverseOrdered_j[k])		return false;

	BzzVector nuInverseOrdered_i  = reaction[i].nuInverse;
	Reorder(&nuInverseOrdered_i, reaction[i].indexInverseOrderedSequence);
	BzzVector nuDirectOrdered_j  = reaction[j].nuDirect;
	Reorder(&nuDirectOrdered_j, reaction[j].indexDirectOrderedSequence);
	
	for (k=1;k<=reaction[i].indexInverse.Size();k++)
		if (nuInverseOrdered_i[k] != nuDirectOrdered_j[k])		return false;

	return true;
}

void OpenSMOKE_CHEMKINInterpreter_KineticsData::CheckReactionRatesReactions()
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

void OpenSMOKE_CHEMKINInterpreter_KineticsData::PrintStoichiometryFile(const std::string file_name)
{
	ofstream fOutput;
	openOutputFileAndControl(fOutput, file_name);
	fOutput.setf(ios::scientific);

	fOutput << species_list.size()-1 << "\t\t"
			<< number_of_reactions << endl;

	for(int i=1;i<=int(species_list.size())-1;i++)
		for(int j=1;j<=species_reaction[i].Size();j++)
			fOutput << i						<< "\t\t"
					<< species_reaction[i][j]	<< "\t\t"
					<< species_nu[i][j]			<< endl;

	fOutput.close();
}

void OpenSMOKE_CHEMKINInterpreter_KineticsData::PrintReactionsFile(const std::string file_name)
{
	ofstream fOutput;
	openOutputFileAndControl(fOutput, file_name);
	fOutput.setf(ios::scientific);

	for(int i=1;i<=number_of_reactions;i++)
		reaction[i].PrintOnFile(fOutput);

	fOutput.close();
}

void OpenSMOKE_CHEMKINInterpreter_KineticsData::Statistics()
{
	int i;
	for(i=1;i<=number_of_reactions;i++)
	{
		if (reaction[i].iReversible	 == true)	iReversible.Append(i);
		if (reaction[i].iReverseRate == true)	iReverseRate.Append(i);
		if (reaction[i].iDuplicate == true)		iDuplicate.Append(i);
		if (reaction[i].iGlobal == true)		iGlobal.Append(i);
	}


	for(i=1;i<=number_of_reactions;i++)
	{
			 if (reaction[i].iLowPressure			== true) 
			 {
				 iFallOff.Append(i);
				 if (reaction[i].kindPressureDependence == PRESSURE_LINDEMANN) iFallOff_Lindemann.Append(i);
				 else if (reaction[i].kindPressureDependence == PRESSURE_TROE) iFallOff_Troe.Append(i);
				 else if (reaction[i].kindPressureDependence == PRESSURE_SRI)  iFallOff_SRI.Append(i);
			 }
		else if	(reaction[i].iHighPressure			== true)
			 {
				 iCABR.Append(i);
				 if (reaction[i].kindPressureDependence == PRESSURE_LINDEMANN) iCABR_Lindemann.Append(i);
				 else if (reaction[i].kindPressureDependence == PRESSURE_TROE) iCABR_Troe.Append(i);
				 else if (reaction[i].kindPressureDependence == PRESSURE_SRI)  iCABR_SRI.Append(i);
			 }	
		else if (reaction[i].iLandauTeller			== true) iLandauTeller.Append(i);
		else if (reaction[i].iPowerSeries			== true) iPowerSeries.Append(i);
		else if (reaction[i].iJanevLanger			== true) iJanevLanger.Append(i);
		else if (reaction[i].iChebishevPolynomials	== true) iChebishevPolynomials.Append(i);
		else if (reaction[i].iPressureLogarithmic	== true) iPressureLogarithmic.Append(i);
		else if (reaction[i].iCollisionEfficiency	== true) iCollisionEfficiency.Append(i);
	}
}

void OpenSMOKE_CHEMKINInterpreter_KineticsData::SummaryOnFile(ofstream &fOutput)
{
	fOutput << "-----------------------------------------------------------------"	<< endl;
	fOutput << "                        CHEMICAL REACTIONS                       "	<< endl;
	fOutput << "-----------------------------------------------------------------"	<< endl;
	fOutput << endl;
	fOutput << "  Number of Reactions:                                  " << number_of_reactions << endl;
	fOutput << "  Number of Reversible Reactions:                       " << iReversible.Size() << endl;
	fOutput << "  Number of Irreversible Reactions:                     " << number_of_reactions-iReversible.Size() << endl;
	fOutput << "  Number of Explicit Reverse Rate Reactions:            " << iReverseRate.Size() << endl;
	fOutput << "  Number of Fall-Off Reactions:                         " << iFallOff.Size() << endl;
	fOutput << "  Number of Chemically-Activated Bimolecular Reactions: " << iCABR.Size() << endl;
	fOutput << "  Number of Landau-Teller Reactions:                    " << iLandauTeller.Size() << endl;
	fOutput << "  Number of Janev-Langer Reactions:                     " << iJanevLanger.Size() << endl;
	fOutput << "  Number of Power-Series Reactions:                     " << iPowerSeries.Size() << endl;
	fOutput << "  Number of Chebichev Polynomial Expansion Reactions:   " << iChebishevPolynomials.Size() << endl;
	fOutput << "  Number of Duplicate Reactions:                        " << iDuplicate.Size() << endl;
	fOutput << "  Number of Global Reactions:                           " << iGlobal.Size() << endl;
	fOutput << "  Number of Bimolecular Coll. Reaction:                 " << iCollisionEfficiency.Size() << endl;
	fOutput << endl;

	if (iFallOff.Size() > 0)
	{
		fOutput << "  Fall-Off Reactions" << endl;
		fOutput << "    Number of Lindemann Reactions:      " << iFallOff_Lindemann.Size() << endl;
		fOutput << "    Number of Troe Reactions:           " << iFallOff_Troe.Size() << endl;
		fOutput << "    Number of SRI Reactions:            " << iFallOff_SRI.Size() << endl;
		fOutput << endl;
	}

	if (iCABR.Size() > 0)
	{
		fOutput << "  Chemically-Activated Bimolecular Reactions" << endl;
		fOutput << "    Number of Lindemann Reactions:      " << iCABR_Lindemann.Size() << endl;
		fOutput << "    Number of Troe Reactions:           " << iCABR_Troe.Size() << endl;
		fOutput << "    Number of SRI Reactions:            " << iCABR_SRI.Size() << endl;
		fOutput << endl;
	}

	fOutput << "-----------------------------------------------------------------"	<< endl;
	fOutput << endl << endl;
}

void OpenSMOKE_CHEMKINInterpreter_KineticsData::PrintBinaryFile(const std::string file_name, OpenSMOKE_CHEMKINInterpreter_ThermoData &thermo, OpenSMOKE_CHEMKINInterpreter_TransportData &transport, OpenSMOKE_CHEMKINInterpreter &interp)
{
	int number_of_species = species_list.size()-1;
	int i;

	BzzSave outputFile;
	BzzSave asciiFile;
	outputFile('*', file_name);
	std::string file_name_ascii = file_name+".ascii";
	asciiFile(file_name_ascii);
	
	// Writing number of species and number of reactions
	outputFile << number_of_species;
	outputFile << number_of_reactions;
	asciiFile << number_of_species;
	asciiFile << number_of_reactions;

	// Writing reaction data
	for(i=1;i<=number_of_reactions;i++)
		reaction[i].PrintOnBinaryFile(outputFile, thermo, transport);
	for(i=1;i<=number_of_reactions;i++)
		reaction[i].PrintOnBinaryFile(asciiFile, thermo, transport);

	// Writing stoichiometry
	{
		outputFile << numDir1;
		outputFile << numDir2;
		outputFile << numDir3;
		outputFile << numDir4;
		outputFile << numDir5;

		outputFile << numInvTot1;
		outputFile << numInvTot2;
		outputFile << numInvTot3;
		outputFile << numInvTot4;
		outputFile << numInvTot5;

		outputFile << numInvEq1;
		outputFile << numInvEq2;
		outputFile << numInvEq3;
		outputFile << numInvEq4;
		outputFile << numInvEq5;

		outputFile << jDir1;
		outputFile << jDir2;
		outputFile << jDir3;
		outputFile << jDir4;
		outputFile << jDir5;
		outputFile << valDir5;

		outputFile << jInvTot1;
		outputFile << jInvTot2;
		outputFile << jInvTot3;
		outputFile << jInvTot4;
		outputFile << jInvTot5;
		outputFile << valInvTot5;

		outputFile << jInvEq1;
		outputFile << jInvEq2;
		outputFile << jInvEq3;
		outputFile << jInvEq4;
		outputFile << jInvEq5;
		outputFile << valInvEq5;

		outputFile << sumNuij;
	}

	{
		asciiFile << numDir1;
		asciiFile << numDir2;
		asciiFile << numDir3;
		asciiFile << numDir4;
		asciiFile << numDir5;

		asciiFile << numInvTot1;
		asciiFile << numInvTot2;
		asciiFile << numInvTot3;
		asciiFile << numInvTot4;
		asciiFile << numInvTot5;

		asciiFile << numInvEq1;
		asciiFile << numInvEq2;
		asciiFile << numInvEq3;
		asciiFile << numInvEq4;
		asciiFile << numInvEq5;

		asciiFile << jDir1;
		asciiFile << jDir2;
		asciiFile << jDir3;
		asciiFile << jDir4;
		asciiFile << jDir5;
		asciiFile << valDir5;

		asciiFile << jInvTot1;
		asciiFile << jInvTot2;
		asciiFile << jInvTot3;
		asciiFile << jInvTot4;
		asciiFile << jInvTot5;
		asciiFile << valInvTot5;

		asciiFile << jInvEq1;
		asciiFile << jInvEq2;
		asciiFile << jInvEq3;
		asciiFile << jInvEq4;
		asciiFile << jInvEq5;
		asciiFile << valInvEq5;

		asciiFile << sumNuij;
	}

	// Writing reaction names on file
	for(i=1;i<=number_of_reactions;i++)
	{
		char name[Constants::REACTION_NAME_SIZE];
		strcpy(name, reaction[i].reaction_string_clean.c_str());
		outputFile.fileSave.write((char*) name, sizeof(name));
		asciiFile << name;
	}

	//if (iGlobal.Size() != 0)
	{
		cout << "    Preparing reaction orders..." << endl;
		PrepareReactionOrders();
	}

	// Writing additional info on stoichiometry
	outputFile << forwardOrders;
	outputFile << backwardOrders;
	outputFile << iGlobal.Size();
	if (interp.iSootMode == true)	outputFile << number_of_species;
	else							outputFile << 0;

	asciiFile << forwardOrders;
	asciiFile << backwardOrders;
	asciiFile << iGlobal.Size();
	if (interp.iSootMode == true)	asciiFile << number_of_species;
	else							asciiFile << 0;

	if (iGlobal.Size() != 0)
	{
/*		outputFile << lambda_numDir1;
		outputFile << lambda_numDir2;
		outputFile << lambda_numDir3;
		outputFile << lambda_numDir4;
		outputFile << lambda_numDir5;

		outputFile << lambda_numInvEq1;
		outputFile << lambda_numInvEq2;
		outputFile << lambda_numInvEq3;
		outputFile << lambda_numInvEq4;
		outputFile << lambda_numInvEq5;

		outputFile << lambda_jDir1;
		outputFile << lambda_jDir2;
		outputFile << lambda_jDir3;
		outputFile << lambda_jDir4;
		outputFile << lambda_jDir5;
		outputFile << lambda_valDir5;

		outputFile << lambda_jInvEq1;
		outputFile << lambda_jInvEq2;
		outputFile << lambda_jInvEq3;
		outputFile << lambda_jInvEq4;
		outputFile << lambda_jInvEq5;
		outputFile << lambda_valInvEq5;
*/
		outputFile << iGlobal;
		asciiFile  << iGlobal;

		for(i=1;i<=iGlobal.Size();i++)
		{
			outputFile << reaction[iGlobal[i]].indexGlobalDirect;
			outputFile << reaction[iGlobal[i]].lambdaGlobalDirect;
			asciiFile << reaction[iGlobal[i]].indexGlobalDirect;
			asciiFile << reaction[iGlobal[i]].lambdaGlobalDirect;
			if (reaction[iGlobal[i]].indexGlobalInverse.Size() !=0)
			{
				outputFile << 1;
				outputFile << reaction[iGlobal[i]].indexGlobalInverse;
				outputFile << reaction[iGlobal[i]].lambdaGlobalInverse;
				asciiFile << 1;
				asciiFile << reaction[iGlobal[i]].indexGlobalInverse;
				asciiFile << reaction[iGlobal[i]].lambdaGlobalInverse;
			}
			else
			{
				outputFile << 0;
				asciiFile  << 0;
			}
		}
	}

	outputFile.End();
	asciiFile.End();
}

void OpenSMOKE_CHEMKINInterpreter_KineticsData::PrepareStoichiometry()
{
	int number_of_species = species_list.size()-1;

	ChangeDimensions(number_of_species, &numDir1);
	ChangeDimensions(number_of_species, &numDir2);
	ChangeDimensions(number_of_species, &numDir3);
	ChangeDimensions(number_of_species, &numDir4);
	ChangeDimensions(number_of_species, &numDir5);

	ChangeDimensions(number_of_species, &numInvTot1);
	ChangeDimensions(number_of_species, &numInvTot2);
	ChangeDimensions(number_of_species, &numInvTot3);
	ChangeDimensions(number_of_species, &numInvTot4);
	ChangeDimensions(number_of_species, &numInvTot5);

	ChangeDimensions(number_of_species, &numInvEq1);
	ChangeDimensions(number_of_species, &numInvEq2);
	ChangeDimensions(number_of_species, &numInvEq3);
	ChangeDimensions(number_of_species, &numInvEq4);
	ChangeDimensions(number_of_species, &numInvEq5);

	ChangeDimensions(number_of_reactions, &sumNuij);
	ChangeDimensions(number_of_reactions, &forwardOrders);
	ChangeDimensions(number_of_reactions, &backwardOrders);

	species_reaction	= new BzzVectorInt[number_of_species+1];
	species_nu			= new BzzVector[number_of_species+1];

	int i;
	for(i=1;i<=number_of_reactions;i++)
	{
		int j;

		for(j=1;j<=reaction[i].nReactants;j++)
		{
			species_reaction[reaction[i].indexDirect[j]].Append(i);
			species_nu[reaction[i].indexDirect[j]].Append(-reaction[i].nuDirect[j]);
		}

		for(j=1;j<=reaction[i].nProducts;j++)
		{
			species_reaction[reaction[i].indexInverse[j]].Append(i);
			species_nu[reaction[i].indexInverse[j]].Append(reaction[i].nuInverse[j]);
		}
	}

	for(i=1;i<=number_of_species;i++)
		for(int k=1;k<=species_reaction[i].Size();k++)
		{
			int j = species_reaction[i][k];

			if(species_nu[i][k] == -1.)			numDir1[i]++;
			else if(species_nu[i][k] == -2.)	numDir2[i]++;
			else if(species_nu[i][k] == -3.)	numDir3[i]++;
			else if(species_nu[i][k] == -0.50)	numDir4[i]++;
			else if(species_nu[i][k] < 0.)		numDir5[i]++;
			
			else if(species_nu[i][k] == 1.)
			{
				numInvTot1[i]++;
				if(reaction[j].iReversible == true)	numInvEq1[i]++;
			}
			else if(species_nu[i][k] == 2.)
			{
				numInvTot2[i]++;
				if(reaction[j].iReversible == true)	numInvEq2[i]++;
			}
			else if(species_nu[i][k] == 3.)
			{
				numInvTot3[i]++;
				if(reaction[j].iReversible == true)	numInvEq3[i]++;
			}
			else if(species_nu[i][k] == 0.50)
			{
				numInvTot4[i]++;
				if(reaction[j].iReversible == true)	numInvEq4[i]++;
			}
			else if(species_nu[i][k] > 0.)
			{
				numInvTot5[i]++;
				if(reaction[j].iReversible == true)	numInvEq5[i]++;
			}			
		}

	ChangeDimensions(numDir1.Sum(),&jDir1);
	ChangeDimensions(numDir2.Sum(),&jDir2);
	ChangeDimensions(numDir3.Sum(),&jDir3);
	ChangeDimensions(numDir4.Sum(),&jDir4);
	ChangeDimensions(numDir5.Sum(),&jDir5);
	ChangeDimensions(numDir5.Sum(),&valDir5);

	ChangeDimensions(numInvTot1.Sum(),&jInvTot1);
	ChangeDimensions(numInvTot2.Sum(),&jInvTot2);
	ChangeDimensions(numInvTot3.Sum(),&jInvTot3);
	ChangeDimensions(numInvTot4.Sum(),&jInvTot4);
	ChangeDimensions(numInvTot5.Sum(),&jInvTot5);
	ChangeDimensions(numInvTot5.Sum(),&valInvTot5);

	ChangeDimensions(numInvEq1.Sum(),&jInvEq1);
	ChangeDimensions(numInvEq2.Sum(),&jInvEq2);
	ChangeDimensions(numInvEq3.Sum(),&jInvEq3);
	ChangeDimensions(numInvEq4.Sum(),&jInvEq4);
	ChangeDimensions(numInvEq5.Sum(),&jInvEq5);
	ChangeDimensions(numInvEq5.Sum(),&valInvEq5);

	int		*jD1 = jDir1.GetHandle();
	int		*jD2 = jDir2.GetHandle();
	int		*jD3 = jDir3.GetHandle();
	int		*jD4 = jDir4.GetHandle();
	int		*jD5 = jDir5.GetHandle();
	double	*vD5 = valDir5.GetHandle();

	int		*jIT1 = jInvTot1.GetHandle();
	int		*jIT2 = jInvTot2.GetHandle();
	int		*jIT3 = jInvTot3.GetHandle();
	int		*jIT4 = jInvTot4.GetHandle();
	int		*jIT5 = jInvTot5.GetHandle();
	double	*vIT5 = valInvTot5.GetHandle();
	
	int		*jIE1 = jInvEq1.GetHandle();
	int		*jIE2 = jInvEq2.GetHandle();
	int		*jIE3 = jInvEq3.GetHandle();
	int		*jIE4 = jInvEq4.GetHandle();
	int		*jIE5 = jInvEq5.GetHandle();
	double	*vIE5 = valInvEq5.GetHandle();

	for(i=1;i<=number_of_species;i++)
		for(int k=1;k<=species_reaction[i].Size();k++)
		{
			int j = species_reaction[i][k];

			if(species_nu[i][k] == -1.)			*jD1++ = j;
			else if(species_nu[i][k] == -2.)	*jD2++ = j;
			else if(species_nu[i][k] == -3.)	*jD3++ = j;
			else if(species_nu[i][k] == -0.50)	*jD4++ = j;
			else if(species_nu[i][k] < 0.)
			{
				*jD5++ = j;
				*vD5++ = fabs(species_nu[i][k]);
			}
		
			else if(species_nu[i][k] == 1.)
			{
				*jIT1++ = j;
				if(reaction[j].iReversible == true)	*jIE1++ = j;
			}
			else if(species_nu[i][k] == 2.)
			{
				*jIT2++ = j;
				if(reaction[j].iReversible == true)	*jIE2++ = j;
			}
			else if(species_nu[i][k] == 3.)
			{
				*jIT3++ = j;
				if(reaction[j].iReversible == true)	*jIE3++ = j;
			}
			else if(species_nu[i][k] == 0.50)
			{
				*jIT4++ = j;
				if(reaction[j].iReversible == true)	*jIE4++ = j;
			}
			else if(species_nu[i][k] > 0.)
			{
				*jIT5++ = j;
				*vIT5++ = species_nu[i][k];
				if(reaction[j].iReversible == true)
				{
					*jIE5++ = j;
					*vIE5++ = species_nu[i][k];
				}
			}	
		}
	{
		jD1 = jDir1.GetHandle();
		jD2 = jDir2.GetHandle();
		jD3 = jDir3.GetHandle();
		jD4 = jDir4.GetHandle();
		jD5 = jDir5.GetHandle();
		vD5 = valDir5.GetHandle();

		jIT1 = jInvTot1.GetHandle();
		jIT2 = jInvTot2.GetHandle();
		jIT3 = jInvTot3.GetHandle();
		jIT4 = jInvTot4.GetHandle();
		jIT5 = jInvTot5.GetHandle();
		vIT5 = valInvTot5.GetHandle();

		for(i=1;i<=number_of_species;i++)
		{
			int k;
			for(k=1;k<=numDir1[i];k++)		sumNuij[*jD1++] -= 1.;
			for(k=1;k<=numDir2[i];k++)		sumNuij[*jD2++] -= 2.;
			for(k=1;k<=numDir3[i];k++)		sumNuij[*jD3++] -= 3.;
			for(k=1;k<=numDir4[i];k++)		sumNuij[*jD4++] -= 0.5;
			for(k=1;k<=numDir5[i];k++)		sumNuij[*jD5++] -= *vD5++;

			for(k=1;k <= numInvTot1[i];k++)	sumNuij[*jIT1++] += 1.;
			for(k=1;k <= numInvTot2[i];k++)	sumNuij[*jIT2++] += 2.;
			for(k=1;k <= numInvTot3[i];k++)	sumNuij[*jIT3++] += 3.;
			for(k=1;k <= numInvTot4[i];k++)	sumNuij[*jIT4++] += 0.5;
			for(k=1;k <= numInvTot5[i];k++)	sumNuij[*jIT5++] += *vIT5++;
		}
	}
}

void OpenSMOKE_CHEMKINInterpreter_KineticsData::FinalChecks()
{
	for(int i=1;i<=number_of_reactions;i++)
	{
		reaction[i].Complete();
		reaction[i].Check();
	}
	CheckDuplicateReactions();
	CheckReactionRatesReactions();
}

void OpenSMOKE_CHEMKINInterpreter_KineticsData::PrepareReactionOrders()
{
	int number_of_species = species_list.size()-1;

	ChangeDimensions(number_of_species, &lambda_numDir1);
	ChangeDimensions(number_of_species, &lambda_numDir2);
	ChangeDimensions(number_of_species, &lambda_numDir3);
	ChangeDimensions(number_of_species, &lambda_numDir4);
	ChangeDimensions(number_of_species, &lambda_numDir5);

	ChangeDimensions(number_of_species, &lambda_numInvEq1);
	ChangeDimensions(number_of_species, &lambda_numInvEq2);
	ChangeDimensions(number_of_species, &lambda_numInvEq3);
	ChangeDimensions(number_of_species, &lambda_numInvEq4);
	ChangeDimensions(number_of_species, &lambda_numInvEq5);

	lambda_species_reaction	= new BzzVectorInt[number_of_species+1];
	lambda_species_nu		= new BzzVector[number_of_species+1];

	int i;
	for(i=1;i<=number_of_reactions;i++)
	{
		if (reaction[i].iGlobal == false)
		{
			int j;

			for(j=1;j<=reaction[i].nReactants;j++)
			{
				lambda_species_reaction[reaction[i].indexDirect[j]].Append(i);
				lambda_species_nu[reaction[i].indexDirect[j]].Append(-reaction[i].nuDirect[j]);
			}

			for(j=1;j<=reaction[i].nProducts;j++)
			{
				lambda_species_reaction[reaction[i].indexInverse[j]].Append(i);
				lambda_species_nu[reaction[i].indexInverse[j]].Append(reaction[i].nuInverse[j]);
			}
		}
	}

	for(i=1;i<=number_of_species;i++)
		for(int k=1;k<=lambda_species_reaction[i].Size();k++)
		{
			int j = lambda_species_reaction[i][k];

			if(lambda_species_nu[i][k] == -1.)			lambda_numDir1[i]++;
			else if(lambda_species_nu[i][k] == -2.)		lambda_numDir2[i]++;
			else if(lambda_species_nu[i][k] == -3.)		lambda_numDir3[i]++;
			else if(lambda_species_nu[i][k] == -0.50)	lambda_numDir4[i]++;
			else if(lambda_species_nu[i][k] < 0.)		lambda_numDir5[i]++;
		
			else if(lambda_species_nu[i][k] == 1.)
				if(reaction[j].iReversible == true)		lambda_numInvEq1[i]++;
			else if(lambda_species_nu[i][k] == 2.)
				if(reaction[j].iReversible == true)		lambda_numInvEq2[i]++;
			else if(lambda_species_nu[i][k] == 3.)
				if(reaction[j].iReversible == true)		lambda_numInvEq3[i]++;
			else if(lambda_species_nu[i][k] == 0.50)
				if(reaction[j].iReversible == true)		lambda_numInvEq4[i]++;
			else if(lambda_species_nu[i][k] > 0.)
				if(reaction[j].iReversible == true)		lambda_numInvEq5[i]++;		
	}


	ChangeDimensions(lambda_numDir1.Sum(),&lambda_jDir1);
	ChangeDimensions(lambda_numDir2.Sum(),&lambda_jDir2);
	ChangeDimensions(lambda_numDir3.Sum(),&lambda_jDir3);
	ChangeDimensions(lambda_numDir4.Sum(),&lambda_jDir4);
	ChangeDimensions(lambda_numDir5.Sum(),&lambda_jDir5);
	ChangeDimensions(lambda_numDir5.Sum(),&lambda_valDir5);

	ChangeDimensions(lambda_numInvEq1.Sum(),&lambda_jInvEq1);
	ChangeDimensions(lambda_numInvEq2.Sum(),&lambda_jInvEq2);
	ChangeDimensions(lambda_numInvEq3.Sum(),&lambda_jInvEq3);
	ChangeDimensions(lambda_numInvEq4.Sum(),&lambda_jInvEq4);
	ChangeDimensions(lambda_numInvEq5.Sum(),&lambda_jInvEq5);
	ChangeDimensions(lambda_numInvEq5.Sum(),&lambda_valInvEq5);

	int		*jD1 = lambda_jDir1.GetHandle();
	int		*jD2 = lambda_jDir2.GetHandle();
	int		*jD3 = lambda_jDir3.GetHandle();
	int		*jD4 = lambda_jDir4.GetHandle();
	int		*jD5 = lambda_jDir5.GetHandle();
	double	*vD5 = lambda_valDir5.GetHandle();
	
	int		*jIE1 = lambda_jInvEq1.GetHandle();
	int		*jIE2 = lambda_jInvEq2.GetHandle();
	int		*jIE3 = lambda_jInvEq3.GetHandle();
	int		*jIE4 = lambda_jInvEq4.GetHandle();
	int		*jIE5 = lambda_jInvEq5.GetHandle();
	double	*vIE5 = lambda_valInvEq5.GetHandle();

	for(i=1;i<=number_of_species;i++)
		for(int k=1;k<=lambda_species_reaction[i].Size();k++)
		{
			int j = lambda_species_reaction[i][k];

			if(lambda_species_nu[i][k] == -1.)			*jD1++ = j;
			else if(lambda_species_nu[i][k] == -2.)		*jD2++ = j;
			else if(lambda_species_nu[i][k] == -3.)		*jD3++ = j;
			else if(lambda_species_nu[i][k] == -0.50)	*jD4++ = j;
			else if(lambda_species_nu[i][k] < 0.)
			{
				*jD5++ = j;
				*vD5++ = fabs(lambda_species_nu[i][k]);
			}
			else if(lambda_species_nu[i][k] == 1.)
				if(reaction[j].iReversible == true)	*jIE1++ = j;
			else if(lambda_species_nu[i][k] == 2.)
				if(reaction[j].iReversible == true)	*jIE2++ = j;
			else if(lambda_species_nu[i][k] == 3.)
				if(reaction[j].iReversible == true)	*jIE3++ = j;
			else if(lambda_species_nu[i][k] == 0.50)
				if(reaction[j].iReversible == true)	*jIE4++ = j;
			else if(lambda_species_nu[i][k] > 0.)
				if(reaction[j].iReversible == true)
				{
					*jIE5++ = j;
					*vIE5++ = lambda_species_nu[i][k];
				}
		}

	// Reaction orders
	{
		if (iGlobal.Size() == 0) 
		{
			int		*jD1 = jDir1.GetHandle();
			int		*jD2 = jDir2.GetHandle();
			int		*jD3 = jDir3.GetHandle();
			int		*jD4 = jDir4.GetHandle();
			int		*jD5 = jDir5.GetHandle();
			double	*vD5 = valDir5.GetHandle();

			int		*jIT1 = jInvTot1.GetHandle();
			int		*jIT2 = jInvTot2.GetHandle();
			int		*jIT3 = jInvTot3.GetHandle();
			int		*jIT4 = jInvTot4.GetHandle();
			int		*jIT5 = jInvTot5.GetHandle();
			double  *vIT5 = valInvTot5.GetHandle();

			for(i=1;i<=number_of_species;i++)
			{
				int k;
				for(k=1;k<=numDir1[i];k++)		forwardOrders[*jD1++] += 1.;
				for(k=1;k<=numDir2[i];k++)		forwardOrders[*jD2++] += 2.;
				for(k=1;k<=numDir3[i];k++)		forwardOrders[*jD3++] += 3.;
				for(k=1;k<=numDir4[i];k++)		forwardOrders[*jD4++] += 0.5;
				for(k=1;k<=numDir5[i];k++)		forwardOrders[*jD5++] += *vD5++;

				for(k=1;k <= numInvTot1[i];k++)	backwardOrders[*jIT1++] += 1.;
				for(k=1;k <= numInvTot2[i];k++)	backwardOrders[*jIT2++] += 2.;
				for(k=1;k <= numInvTot3[i];k++)	backwardOrders[*jIT3++] += 3.;
				for(k=1;k <= numInvTot4[i];k++)	backwardOrders[*jIT4++] += 0.5;
				for(k=1;k <= numInvTot5[i];k++)	backwardOrders[*jIT5++] += *vIT5++;
			}
		}
		else
		{
			int		*jD1 = lambda_jDir1.GetHandle();
			int		*jD2 = lambda_jDir2.GetHandle();
			int		*jD3 = lambda_jDir3.GetHandle();
			int		*jD4 = lambda_jDir4.GetHandle();
			int		*jD5 = lambda_jDir5.GetHandle();
			double  *vD5 = lambda_valDir5.GetHandle();

			int		*jIE1 = lambda_jInvEq1.GetHandle();
			int		*jIE2 = lambda_jInvEq2.GetHandle();
			int		*jIE3 = lambda_jInvEq3.GetHandle();
			int		*jIE4 = lambda_jInvEq4.GetHandle();
			int		*jIE5 = lambda_jInvEq5.GetHandle();
			double  *vIE5 = lambda_valInvEq5.GetHandle();

			for(i=1;i<=number_of_species;i++)
			{
				int k;
				for(k=1;k<=lambda_numDir1[i];k++)		forwardOrders[*jD1++] += 1.;
				for(k=1;k<=lambda_numDir2[i];k++)		forwardOrders[*jD2++] += 2.;
				for(k=1;k<=lambda_numDir3[i];k++)		forwardOrders[*jD3++] += 3.;
				for(k=1;k<=lambda_numDir4[i];k++)		forwardOrders[*jD4++] += 0.5;
				for(k=1;k<=lambda_numDir5[i];k++)		forwardOrders[*jD5++] += *vD5++;

				for(k=1;k<=lambda_numInvEq1[i];k++)		backwardOrders[*jIE1++] += 1.;
				for(k=1;k<=lambda_numInvEq2[i];k++)		backwardOrders[*jIE2++] += 2.;
				for(k=1;k<=lambda_numInvEq3[i];k++)		backwardOrders[*jIE3++] += 3.;
				for(k=1;k<=lambda_numInvEq4[i];k++)		backwardOrders[*jIE4++] += 0.5;
				for(k=1;k<=lambda_numInvEq5[i];k++)		backwardOrders[*jIE5++] += *vIE5++;
			}

			for(i=1;i<=number_of_reactions;i++)
				if (reaction[i].iGlobal == true)
				{
					int j;
					
					forwardOrders[i] = 0.;
					backwardOrders[i] = 0.;
					for(j=1;j<=reaction[i].lambdaGlobalDirect.Size();j++)
						forwardOrders[i]  += reaction[i].lambdaGlobalDirect[j];
					for(j=1;j<=reaction[i].lambdaGlobalInverse.Size();j++)
						backwardOrders[i] += reaction[i].lambdaGlobalInverse[j];
				}

			for(i=1;i<=number_of_reactions;i++)
				if (reaction[i].iThirdBody == true)
				{	
					forwardOrders[i] +=1.;
					backwardOrders[i] +=1.;
				}
		}
	}
}


