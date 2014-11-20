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

#include <iomanip>
#include "sstream"
#include "idealreactors/flame1d/OpenSMOKE_Flame1D.h"
#include "idealreactors/flame1d/OpenSMOKE_Flame1D_FlameSpeedManager.h" 

void OpenSMOKE_Flame1D_FlameSpeedManager::ErrorMessage(const std::string message)
{
    cout << endl;
    cout << "Class:  OpenSMOKE_Flame1D_FlameSpeedManager"	<< endl;
    cout << "Object: " << name_object			<< endl;
    cout << "Error:  " << message				<< endl;
    cout << "Press a key to continue... "		<< endl;
    getchar();
    exit(-1);
}

void OpenSMOKE_Flame1D_FlameSpeedManager::WarningMessage(const std::string message)
{
    cout << endl;
    cout << "Class:  OpenSMOKE_Flame1D_FlameSpeedManager"	<< endl;
    cout << "Object: "		<< name_object		<< endl;
    cout << "Warning:  "	<< message			<< endl;
    cout << "Press a key to continue... "		<< endl;
    getchar();
}

OpenSMOKE_Flame1D_FlameSpeedManager::OpenSMOKE_Flame1D_FlameSpeedManager()
{
	name_object			= "[Name not assigned]";
	iEquivalenceRatios = false;
	iMoleFractions = false;
	iMassFractions = false;
	iMoles = false;
	iMasses = false;
}

void OpenSMOKE_Flame1D_FlameSpeedManager::SetName(const std::string name)
{
	name_object = name;
}

void OpenSMOKE_Flame1D_FlameSpeedManager::SetupFromFile(OpenSMOKE_Flame1D *_flame, OpenSMOKE_Flame1D_DataManager *_data)
{
	flame = _flame;
	data  = _data;

	std::string dummy;


	ifstream fInput;
	openInputFileAndControl(fInput, data->flameSpeedAnalysisFileName);

	fInput >> dummy;
	if		(dummy == "#EquivalenceRatios")		iEquivalenceRatios = true;
	else if	(dummy == "#MoleFractions")			iMoleFractions = true;
	else if	(dummy == "#MassFractions")			iMassFractions = true;
	else if	(dummy == "#Moles")					iMoles = true;
	else if	(dummy == "#Masses")				iMasses = true;
	else	ErrorMessage("#EquivalenceRatios || #MoleFractions || #MassFractions || #Moles || #Masses keywords are expected!");

	if (iEquivalenceRatios == true)
	{
		for(;;)
		{
			fInput >> dummy;
			if (dummy == "#End")	break;
			equivalence_ratios.Append(atof(dummy.c_str()));
		}
	}
	else
	{
		char line[Constants::COMMENT_SIZE];
		fInput.getline(line, Constants::COMMENT_SIZE);

		for(;;)
		{
			vector<string> list;
			fInput.getline(line, Constants::COMMENT_SIZE);

			if (CheckForBlankLine(line)==false)
			{
				TokenizeString(line, list);
				
				if (list[0] == "#End")	break;

				vector<string> names;
				vector<double> values;

				int jStart = 1;
				if (list.size()%2 != 0)
				{
					tag_flame.push_back(list[0]);
					jStart=2;
				}
				
				int j;
				double sum = 0.;
				for(j=jStart;j<=int(list.size());j+=2)
				{
					double number = atof(list[j].c_str());
					names.push_back(list[j-1]);
					values.push_back(number);
					sum+=number;
				}			
				
				if (iMoles == true || iMasses == true)
					for(int j=1;j<=int(values.size());j++)
						values[j-1] /= sum;
				
				list_names.push_back(names);
				list_values.push_back(values);
			}
		}


	}
	fInput.close();

	inlet_stream.AssignKineticScheme(*flame->mix);
	inlet_stream.AssignTemperature(data->TC, "K");
	inlet_stream.AssignPressure(data->P_Pascal, "Pa");
	inlet_stream.AssignMassFlowRate(data->MassFlowRate, "kg/s");

	if (iEquivalenceRatios == true)
	{
		if (data->iEquivalenceRatioForPremixedFlames == false)
			ErrorMessage("#EquivalenceRatios option is allowed only if the inlet composition is specified through the #EquivalenceRatio option...");
		
		N = equivalence_ratios.Size();

		inlet_stream.AssignEquivalenceRatio(equivalence_ratios[1], data->fuel_names, data->moles_fuel, data->oxidizer_names, data->moles_oxidizer);
		inlet_stream.lock();
	}
	else
	{
		if (data->iEquivalenceRatioForPremixedFlames == true)
			ErrorMessage("#MoleFractions option is allowed only if the inlet composition is not specified through the #EquivalenceRatio option...");
		
		N = list_names.size();

		if (iMoleFractions==true || iMoles==true)
			inlet_stream.AssignMoleFractions(list_names[1-1], list_values[1-1]);
		if (iMassFractions==true || iMasses==true)
			inlet_stream.AssignMassFractions(list_names[1-1], list_values[1-1]);
		inlet_stream.lock();
	}

	flame->BCW_C = inlet_stream.omega;

	openOutputFileAndControl(fFlameSpeed, flame->nameFolderAdditionalData + "/FlameSpeeds.out");
	openOutputFileAndControl(fLog, flame->nameFolderAdditionalData + "/FlameSpeeds.log");
	fFlameSpeed.setf(ios::scientific);
	fLog.setf(ios::scientific);
//	if (iEquivalenceRatios == true)	fFlameSpeed << "Phi" << "\t\t\t" << "Velocity[m/s]" << "\t\t" << "MaxT[K]" << endl;
//	else							fFlameSpeed << "Index" << "\t\t\t" << "Velocity[m/s]" << "\t\t" << "MaxT[K]" << endl;

//	if (equivalence_ratios[1] <= 0.99999 || equivalence_ratios[1] >= 1.0001)
//		ErrorMessage("First equivalence ratio must be equal to 1...");

	tag_unit_equivalence_ratio = 0;
	for (int j=2;j<=N;j++)
		if (equivalence_ratios[j] >= 0.99999*equivalence_ratios[1] && equivalence_ratios[j]*equivalence_ratios[1] <= 1.0001)
			tag_unit_equivalence_ratio = j;

	if (equivalence_ratios[2]>equivalence_ratios[1] && equivalence_ratios[N]<equivalence_ratios[1])
		if (tag_unit_equivalence_ratio == 0)
			ErrorMessage("First equivalence ratio must be repeated for starting the lean side...");

	if (equivalence_ratios[2]<equivalence_ratios[1] && equivalence_ratios[N]>equivalence_ratios[1])
		if (tag_unit_equivalence_ratio == 0)
			ErrorMessage("First equivalence ratio must be repeated for starting the rich side...");
}

void OpenSMOKE_Flame1D_FlameSpeedManager::Run()
{
	OpenSMOKE_Flame1D_Solution unit_solution;

	ChangeDimensions(N, &list_of_phi);
	ChangeDimensions(N, &list_of_flame_T);
	ChangeDimensions(N, &list_of_flame_speeds);

	double tag_number;
	for (int count=1;count<=N;count++)
	{
		if (tag_unit_equivalence_ratio == count)
		{
			// Paste solution
			cout << "Paste from Unit_Solution..." << endl;
			flame->PasteFromExternalSolution(unit_solution);
			continue;
		}

		if (iEquivalenceRatios == true)	tag_number = equivalence_ratios[count];
		else							tag_number = count;

		if (iEquivalenceRatios == true)
			inlet_stream.ChangeEquivalenceRatio(equivalence_ratios[count],data->fuel_names, data->moles_fuel, data->oxidizer_names, data->moles_oxidizer);
		if (iMoleFractions==true || iMoles==true)
			inlet_stream.AssignMoleFractions(list_names[count-1], list_values[count-1]);
		if (iMassFractions==true || iMasses==true)
			inlet_stream.AssignMassFractions(list_names[count-1], list_values[count-1]);
		inlet_stream.lock();
			
		
		flame->BCW_C = inlet_stream.omega;
		flame->data->AssignFuelMassFractions(inlet_stream.omega);

		PrintComposition(fLog, tag_number);

		// Solve new flame
		fLog << "0" << "\t" << flame->Np << "\t" << flame->H[1]/flame->rho[1] << endl;
		flame->solveDAE_Premixed(PREMIXED_FLAMESPEED, 1.e5);
		fLog << "1" << "\t" << flame->Np << "\t" << flame->H[1]/flame->rho[1] << endl;
		
		// Refine if needed
		for (int j=1;j<=3;j++)
		{	
			bool iCalculate = true;
			iCalculate = flame->newPoints("TEMPERATURE", 'D');
			if (iCalculate == true)	flame->solveDAE_Premixed(PREMIXED_FLAMESPEED, 1.e5);
			else break;
			fLog << j+2 << "\t" << flame->Np << "\t" << flame->H[1]/flame->rho[1] << endl;
		}
		for (int j=1;j<=3;j++)
		{	
			bool iCalculate = true;
			iCalculate = flame->newPoints("TEMPERATURE", 'C');
			if (iCalculate == true)	flame->solveDAE_Premixed(PREMIXED_FLAMESPEED, 1.e5);
			else break;
			fLog << j+2 << "\t" << flame->Np << "\t" << flame->H[1]/flame->rho[1] << endl;
		}
		for (int j=1;j<=3;j++)
		{	
			bool iCalculate = true;
			iCalculate = flame->newPoints("TEMPERATURE", 'D');
			if (iCalculate == true)	flame->solveDAE_Premixed(PREMIXED_FLAMESPEED, 1.e5);
			else break;
			fLog << j+2 << "\t" << flame->Np << "\t" << flame->H[1]/flame->rho[1] << endl;
		}

		stringstream number;
		number << tag_number;
		std::string nameFileSolutionComplete = flame->nameFolderSteadyData + "/Solution_"; + "Phi_" + number.str() + ".out";
		flame->printOnFile(nameFileSolutionComplete);

		fFlameSpeed << tag_number					<< "\t" 
					<< flame->H[1]/flame->rho[1]	<< "\t"
					<< flame->T.Max()				<< "\t";
		
		list_of_phi[count] = tag_number;
		list_of_flame_speeds[count] = flame->H[1]/flame->rho[1];
		list_of_flame_T[count] = flame->T.Max();

		if (tag_flame.size() > 0)	fFlameSpeed << tag_flame[count-1];
					
		fFlameSpeed << endl;

		Print(tag_number);

		if (count == 1)
		{
			cout << "Save Unit Solution..." << endl;
			unit_solution.PasteFromExternalSolution(*flame, *flame->data);
		}
	}

	BzzVectorInt iRemove;
	for (int i=1;i<=N;i++)
		if (list_of_phi[i] <= 0.1)	iRemove.Append(i);
	if (iRemove.Size()>0)
	{
		list_of_phi.DeleteElements(iRemove);
		list_of_flame_speeds.DeleteElements(iRemove);
		list_of_flame_T.DeleteElements(iRemove);
	}
	BzzVectorInt iSort;
	Sort(&list_of_phi, &iSort);
	Reorder(&list_of_flame_speeds, iSort);
	Reorder(&list_of_flame_T, iSort);

	fFlameSpeed.close();
}

void OpenSMOKE_Flame1D_FlameSpeedManager::PrintComposition(ofstream &fOut, const double tag_number)
{
	int j;

	std::string tag_string;
	if (iEquivalenceRatios == true)	tag_string = " Equivalence ratio: ";
	else							tag_string = " Flame Index: ";

	fOut << endl;
	fOut << "----------------------------------------------------------------" << endl;
	fOut << tag_string << tag_number			                   << endl;
	fOut << "----------------------------------------------------------------" << endl;
	fOut << " " << setw(16) << left << "Name" << "omega" << "\t\t" << "x" << endl;
	for(j=1;j<=flame->mix->NumberOfSpecies();j++)
		if (inlet_stream.omega[j] > 1.e-16)
			fOut << " " << setw(16) << left << flame->mix->names[j] << inlet_stream.omega[j] << "\t" << inlet_stream.x[j] << endl;
	fOut << endl;

	cout << endl;
	cout << "----------------------------------------------------------------" << endl;
	cout << tag_string << tag_number			                   << endl;
	cout << "----------------------------------------------------------------" << endl;
	cout << " " << setw(16) << left << "Name" << "omega" << "\t\t" << "x" << endl;
	for(j=1;j<=flame->mix->NumberOfSpecies();j++)
		if (inlet_stream.omega[j] > 1.e-16)
			cout << " " << setw(16) << left << flame->mix->names[j] << inlet_stream.omega[j] << "\t" << inlet_stream.x[j] << endl;
	cout << endl;
}

void OpenSMOKE_Flame1D_FlameSpeedManager::Print(const double tag_number)
{
	stringstream phi;
	std::string fileName_BackUp;
	std::string fileName_BackUpInputData;
	std::string fileNameOutput;
			
	if (iEquivalenceRatios == true)
	{
		phi << tag_number;
		fileNameOutput				= flame->nameFolderAdditionalData + "/Solution_Phi_" + phi.str() + ".out";
		fileName_BackUp				= flame->nameFolderAdditionalData + "/BackUp_Phi_" + phi.str();
		fileName_BackUpInputData	= flame->nameFolderAdditionalData + "/BackUp_Phi_" + phi.str();
	}
	else
	{
		int tag = int(tag_number);
		if (tag_flame.size() == 0)
		{
			phi << tag;
			fileNameOutput				= flame->nameFolderAdditionalData + "/Solution_Index_" + phi.str() + ".out";
			fileName_BackUp				= flame->nameFolderAdditionalData + "/BackUp_Index_" + phi.str();
			fileName_BackUpInputData	= flame->nameFolderAdditionalData + "/BackUp_Index_" + phi.str();
		}
		else
		{
			fileNameOutput				= flame->nameFolderAdditionalData + "/Solution_" + tag_flame[tag-1] + ".out";
			fileName_BackUp				= flame->nameFolderAdditionalData + "/BackUp_" + tag_flame[tag-1];
			fileName_BackUpInputData	= flame->nameFolderAdditionalData + "/BackUp_" + tag_flame[tag-1];
		}
	}

	cout << "Printing BackUp Data on file: " << fileName_BackUp << ".out" << endl;
	flame->printBackUpOnlyInputData(fileName_BackUpInputData);
	flame->printBackUpOnlyData(fileName_BackUp);
	flame->printOnFile(fileNameOutput);
}
