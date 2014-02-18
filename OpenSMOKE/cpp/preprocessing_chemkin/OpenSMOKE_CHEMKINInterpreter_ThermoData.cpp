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
#include "preprocessing_chemkin/OpenSMOKE_CHEMKINInterpreter_ThermoData.h"
#include "preprocessing_chemkin/OpenSMOKE_CHEMKINInterpreter_ThermoSpecies.h"
#include "preprocessing_chemkin/OpenSMOKE_CHEMKINInterpreter_ElementsData.h"
#include <iomanip>

void OpenSMOKE_CHEMKINInterpreter_ThermoData::ErrorMessage(const string message)
{
    cout << endl;
    cout << "Class:  OpenSMOKE_CHEMKINInterpreter_ThermoData"	<< endl;
    cout << "Object: " << name_object		<< endl;
    cout << "Error:  " << message           << endl;
    cout << "Press a key to continue... "   << endl;
    getchar();
    exit(-1);
}

void OpenSMOKE_CHEMKINInterpreter_ThermoData::WarningMessage(const string message)
{
    cout << endl;
    cout << "Class:  OpenSMOKE_CHEMKINInterpreter_ThermoData"	<< endl;
    cout << "Object: "		<< name_object	<< endl;
    cout << "Warning:  "	<< message      << endl;
    cout << "Press a key to continue... "   << endl;
    getchar();
}

OpenSMOKE_CHEMKINInterpreter_ThermoData::OpenSMOKE_CHEMKINInterpreter_ThermoData()
{

}

OpenSMOKE_CHEMKINInterpreter_ThermoData::~OpenSMOKE_CHEMKINInterpreter_ThermoData()
{

}

void OpenSMOKE_CHEMKINInterpreter_ThermoData::ReadThermoData(const string file_name, ofstream *_fLog)
{
	const int SIZE = 400;
	char comment[SIZE];

	fLog = _fLog;

	ifstream fInput;
	openInputFileAndControl(fInput, file_name);

	// ---------------------------------------------------------------
	// Reading lines
	// ---------------------------------------------------------------
	lines.push_back("List of lines");

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
		if (CheckForBlankLine(lines[i]) == true)		indexBlankLines.Append(i);
		else if (CheckForCommentLine(lines[i]) == true)	indexCommentLines.Append(i);
		else if (CheckForEndLine(lines[i]) == true)		indexCommentLines.Append(i);		
		else											indexLines.Append(i);
	}

	int count_additional=0;
	for(i=1;i<=indexLines.Size();i++)
		count_additional += StringFind(lines[indexLines[i]], "1&");
	total_number_of_species = (indexLines.Size()-count_additional)/4;

	*fLog << " ----------------------------------------------------------------" << endl;
	*fLog << "                     Thermodynamic Database                      " << endl;
	*fLog << " ----------------------------------------------------------------" << endl;
	*fLog << "    Total number of full lines:    " << indexLines.Size()			<< endl;
	*fLog << "    Total number of blank lines:   " << indexBlankLines.Size()		<< endl;
	*fLog << "    Total number of comment lines: " << indexCommentLines.Size()	<< endl;
	*fLog << "    Total number of lines:         " << number_of_lines			<< endl;
	*fLog << "    Total number of species:       " << total_number_of_species			<< endl;
	*fLog << " ----------------------------------------------------------------" << endl;

	species = new OpenSMOKE_CHEMKINInterpreter_ThermoSpecies[total_number_of_species+1];
	
	vector<string> instructions;
	
	// ---------------------------------------------------------------
	// First line
	// ---------------------------------------------------------------
	SeparateInstructions(lines[indexLines[1]], instructions);
	if (instructions.size()-1 > 2)		ErrorMessage("Too many arguments in the first line");
	if (instructions[1] != "THERMO")	ErrorMessage("Expected: THERMO - Found: " + instructions[1]);
	if (instructions.size()-1 == 2)
		if (instructions[2] != "ALL")	ErrorMessage("Expected: ALL - Found: " + instructions[2]);

	// ---------------------------------------------------------------
	// Second line
	// ---------------------------------------------------------------
	SeparateInstructions(lines[indexLines[2]], instructions);
	if (instructions.size()-1 != 3)		ErrorMessage("Too many arguments in the second line");
	tmin  = atof(instructions[1].c_str());
	tmean = atof(instructions[2].c_str());
	tmax  = atof(instructions[3].c_str());

	int j=3;
	for(i=1;i<=total_number_of_species;i++)
	{
		species[i].AssignTemperatures(tmin, tmean, tmax);
		species[i].ReadMainData(lines[indexLines[j]], indexLines[j]);	j++;
		if (species[i].iContinuation==true)
		{	
			species[i].ReadAdditionalLine(lines[indexLines[j]]);		j++;
		}
		species[i].ReadFirstLine(lines[indexLines[j]]);					j++;
		species[i].ReadSecondLine(lines[indexLines[j]]);				j++;
		species[i].ReadThirdLine(lines[indexLines[j]]);					j++;

		species[i].Analyze();
	}
}

BzzVectorInt OpenSMOKE_CHEMKINInterpreter_ThermoData::GiveMeSpeciesIndices(vector<string> &list)
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
		//	indices[j] = 1;//TODO
		//	break;
			ErrorMessage("The " + list[j] + " species is not available in the thermodynamics database");
		}
	}

	return indices;
}

void OpenSMOKE_CHEMKINInterpreter_ThermoData::GiveMeMolecularWeights(BzzVectorInt &indices, OpenSMOKE_CHEMKINInterpreter_ElementsData &elements)
{
	int i;

	BzzVector mw_elements(elements.elements_in_list_indices.size()-1);
	for(i=1;i<=int(elements.elements_in_list.size())-1;i++)
		mw_elements[i] = elements.elements_mw[elements.elements_in_list_indices[i]];

	BzzVector mw_isotope(elements.isotope_in_list.size()-1);
	for(i=1;i<=int(elements.isotope_in_list.size())-1;i++)
		mw_isotope[i] = elements.isotope_mw_in_list[i];
	
	for(int j=1;j<=indices.Size();j++)
	{
		int h = indices[j];
		ChangeDimensions(elements.elements_in_list.size()-1, &species[h].element_indices);
		ChangeDimensions(elements.isotope_in_list.size()-1,  &species[h].isotope_indices);

		for(int k=1;k<=int(species[h].element_names.size())-1;k++)
		{
			int  i;
			bool jFound = false;
			for(i=1;i<=int(elements.elements_in_list.size())-1;i++)
				if (caseInsCompare(species[h].element_names[k], elements.elements_in_list[i])==true)
				{
					jFound = true;
					species[h].element_indices[i] = species[h].element_numbers[k];
					break;
				}

			if (jFound == false)
			{

				for(i=1;i<=int(elements.isotope_in_list.size())-1;i++)
					if (caseInsCompare(species[h].element_names[k], elements.isotope_in_list[i])==true)
					{
						jFound = true;
						species[h].isotope_indices[i] = species[h].element_numbers[k];
						break;
					}
			}

			if (jFound == false)
			{
				for (int kk=0;kk<=species[h].element_names.size()-1;kk++)
					cout << kk << species[h].element_names[k] << endl;
				cout << species[h].name << endl;
				ErrorMessage("The " + species[h].element_names[k] + " element is not available in the database");

			}
		}

		species[h].mw  = Dot(species[h].element_indices, mw_elements);
		species[h].mw += Dot(species[h].isotope_indices, mw_isotope);
	}
}

void OpenSMOKE_CHEMKINInterpreter_ThermoData::PrintThermoFile(const string file_name, BzzVectorInt &indices)
{
	ofstream fOutput;
	openOutputFileAndControl(fOutput, file_name);
	fOutput.setf(ios::scientific);

	for(int j=1;j<=indices.Size();j++)
	{
		int i;

		fOutput << species[indices[j]].name	<< "\t";
		fOutput.precision(12);
		fOutput << species[indices[j]].mw	<< endl;
		fOutput.precision(8);
		for(i=1;i<=7;i++)	
			fOutput << species[indices[j]].upper[i] << "\t";
		fOutput << endl;
		for(i=1;i<=7;i++)	
			fOutput << species[indices[j]].lower[i] << "\t";
		fOutput << endl;
		fOutput << species[indices[j]].t_min	<< "\t";
		fOutput << species[indices[j]].t_mean	<< "\t";
		fOutput << species[indices[j]].t_max	<< endl;
	}

	fOutput.close();
}

void OpenSMOKE_CHEMKINInterpreter_ThermoData::SummaryOnFile(ofstream &fOutput, BzzVectorInt &indices, OpenSMOKE_CHEMKINInterpreter_ElementsData &elements)
{
	int i;

	fOutput << "------------------------------------------------------------------------";
	fOutput << setfill('-');
	fOutput << setw(8*(elements.elements_in_list.size()-1)+2) << "-" << endl;
	fOutput << setfill(' ');
	fOutput << "                                                Molecular        Temperature        Elements" << endl;
	fOutput << "  Species                  Phase  Charge         weight         Low      High";
	for(i=1;i<=int(elements.elements_in_list.size())-1;i++)
		fOutput << setw(8) << right << elements.elements_in_list[i];
	fOutput << endl;
	fOutput << "------------------------------------------------------------------------";
	fOutput << setfill('-');
	fOutput << setw(8*(elements.elements_in_list.size()-1)+2) << "-" << endl;
	fOutput << setfill(' ');
	for(i=1;i<=indices.Size();i++)
	{
		fOutput << right << setw(5) << i;
		fOutput << ". ";
		fOutput << setw(20) << left  << species[indices[i]].name;
		fOutput << setw(3)  << right << species[indices[i]].phase;
		fOutput << setw(7)  << right << 0;
		fOutput << setw(20) << fixed << right << setprecision(6)  << species[indices[i]].mw;
		fOutput << setw(10) << right << setprecision(2)  << species[indices[i]].t_min;
		fOutput << setw(10) << right << setprecision(2)  << species[indices[i]].t_max;
		for(int j=1;j<=int(elements.elements_in_list.size())-1;j++)
			fOutput << setw(8) << setprecision(0) << fixed << right << species[indices[i]].element_indices[j];
		fOutput << endl;
	}
	fOutput << "------------------------------------------------------------------------";
	fOutput << setfill('-');
	fOutput << setw(8*(elements.elements_in_list.size()-1)+2) << "-" << endl;
	fOutput << setfill(' ');
	fOutput << endl;


//	for(i=1;i<=indices.Size();i++)
//	{
//		fOutput << species[indices[i]].name << "\t  ";
//		for(int j=1;j<=elements.elements_in_list.size()-1;j++)
//		{
//			if (species[indices[i]].element_indices[j] != 0.)
//				fOutput << elements.elements_in_list[j] << "_Element  " << species[indices[i]].element_indices[j] << "  ";
//		}
//		fOutput << "/" << endl;
//	}
}

void OpenSMOKE_CHEMKINInterpreter_ThermoData::SummaryOnFile(ofstream &fOutput, BzzVectorInt &site_indices, BzzVectorInt &bulk_indices, OpenSMOKE_CHEMKINInterpreter_ElementsData &elements)
{
	int i;

	fOutput << "------------------------------------------------------------------------";
	fOutput << setfill('-');
	fOutput << setw(8*(elements.elements_in_list.size()-1)+2) << "-" << endl;
	fOutput << setfill(' ');
	fOutput << "                                                Molecular        Temperature        Elements" << endl;
	fOutput << "  Species                  Phase  Charge         weight         Low      High";
	for(i=1;i<=int(elements.elements_in_list.size())-1;i++)
		fOutput << setw(8) << right << elements.elements_in_list[i];
	fOutput << endl;
	fOutput << "------------------------------------------------------------------------";
	fOutput << setfill('-');
	fOutput << setw(8*(elements.elements_in_list.size()-1)+2) << "-" << endl;
	fOutput << setfill(' ');
	
	// Site
	for(i=1;i<=site_indices.Size();i++)
	{
		fOutput << right << setw(5) << i;
		fOutput << ". ";
		fOutput << setw(20) << left  << species[site_indices[i]].name;
		fOutput << setw(3)  << right << species[site_indices[i]].phase;
		fOutput << setw(7)  << right << 0;
		fOutput << setw(20) << fixed << right << setprecision(6)  << species[site_indices[i]].mw;
		fOutput << setw(10) << right << setprecision(2)  << species[site_indices[i]].t_min;
		fOutput << setw(10) << right << setprecision(2)  << species[site_indices[i]].t_max;
		for(int j=1;j<=int(elements.elements_in_list.size())-1;j++)
			fOutput << setw(8) << setprecision(0) << fixed << right << species[site_indices[i]].element_indices[j];
		fOutput << endl;
	}

	// Bulk
	for(i=1;i<=bulk_indices.Size();i++)
	{
		fOutput << right << setw(5) << i;
		fOutput << ". ";
		fOutput << setw(20) << left  << species[bulk_indices[i]].name;
		fOutput << setw(3)  << right << species[bulk_indices[i]].phase;
		fOutput << setw(7)  << right << 0;
		fOutput << setw(20) << fixed << right << setprecision(6)  << species[bulk_indices[i]].mw;
		fOutput << setw(10) << right << setprecision(2)  << species[bulk_indices[i]].t_min;
		fOutput << setw(10) << right << setprecision(2)  << species[bulk_indices[i]].t_max;
		for(int j=1;j<=int(elements.elements_in_list.size())-1;j++)
			fOutput << setw(8) << setprecision(0) << fixed << right << species[bulk_indices[i]].element_indices[j];
		fOutput << endl;
	}

	fOutput << "------------------------------------------------------------------------";
	fOutput << setfill('-');
	fOutput << setw(8*(elements.elements_in_list.size()-1)+2) << "-" << endl;
	fOutput << setfill(' ');
	fOutput << endl;


//	for(i=1;i<=indices.Size();i++)
//	{
//		fOutput << species[indices[i]].name << "\t  ";
//		for(int j=1;j<=elements.elements_in_list.size()-1;j++)
//		{
//			if (species[indices[i]].element_indices[j] != 0.)
//				fOutput << elements.elements_in_list[j] << "_Element  " << species[indices[i]].element_indices[j] << "  ";
//		}
//		fOutput << "/" << endl;
//	}
}

void OpenSMOKE_CHEMKINInterpreter_ThermoData::ProcessThermoData(BzzVectorInt &indices)
{
	int k;
	double RGAS_su_PM;

	int number_of_species = indices.Size();

	aDH		=	new BzzVector[number_of_species+1];
	bDH		=   new BzzVector[number_of_species+1];
	aDS		=   new BzzVector[number_of_species+1];
	bDS		=   new BzzVector[number_of_species+1];
	CpHT	=	new BzzVector[number_of_species+1];
	CpLT	=	new BzzVector[number_of_species+1];

	for(k=1;k<=number_of_species;k++)
	{
		ChangeDimensions(6,&aDH[k]);
		ChangeDimensions(6,&bDH[k]);
		ChangeDimensions(6,&aDS[k]);
		ChangeDimensions(6,&bDS[k]);
		ChangeDimensions(5,&CpHT[k]);
		ChangeDimensions(5,&CpLT[k]);
	}

	ChangeDimensions(number_of_species,&M);
	ChangeDimensions(number_of_species,&T1);
	ChangeDimensions(number_of_species,&T2);
	ChangeDimensions(number_of_species,&T3);

	for(k=1;k<=number_of_species;k++)
	{
		int j = indices[k];

		M[k]=species[j].mw;
		RGAS_su_PM = Constants::R_J_mol/(M[k]*1.e-3);

		CpHT[k][1] = species[j].upper[1]*RGAS_su_PM;		// [J/kgK]
		CpHT[k][2] = species[j].upper[2]*RGAS_su_PM;
		CpHT[k][3] = species[j].upper[3]*RGAS_su_PM;
		CpHT[k][4] = species[j].upper[4]*RGAS_su_PM;
		CpHT[k][5] = species[j].upper[5]*RGAS_su_PM;

		CpLT[k][1] = species[j].lower[1]*RGAS_su_PM;
		CpLT[k][2] = species[j].lower[2]*RGAS_su_PM;
		CpLT[k][3] = species[j].lower[3]*RGAS_su_PM;
		CpLT[k][4] = species[j].lower[4]*RGAS_su_PM;
		CpLT[k][5] = species[j].lower[5]*RGAS_su_PM;

		aDH[k][1] = species[j].upper[1];				// L'entalpia viene calcolata in forma
		aDH[k][2] = .5 * species[j].upper[2];			// adimensionale, cio?dividendo per Rgas e
		aDH[k][3] = 1./3. * species[j].upper[3];		// per la temperatura in modo tale da
		aDH[k][4] = .25 * species[j].upper[4];			// poter calcolare pi velocemente le costanti
		aDH[k][5] = .2 * species[j].upper[5];			// di equilibrio
		aDH[k][6] = species[j].upper[6];

		bDH[k][1] = species[j].lower[1];
		bDH[k][2] = .5 * species[j].lower[2];
		bDH[k][3] = 1./3. * species[j].lower[3];
		bDH[k][4] = .25 * species[j].lower[4];
		bDH[k][5] = .2 * species[j].lower[5];
		bDH[k][6] = species[j].lower[6];

		aDS[k][1] = species[j].upper[1];				// L' entropia viene calcolata in questo modo
		aDS[k][2] = species[j].upper[2];				// cio?adimensionale perch?ci?velocizza
		aDS[k][3] = .5 * species[j].upper[3];		// il calcolo delle costanti di equilibrio
		aDS[k][4] = 1./3. * species[j].upper[4];
		aDS[k][5] = .25 * species[j].upper[5];
		aDS[k][6] = species[j].upper[7];

		bDS[k][1] = species[j].lower[1];
		bDS[k][2] = species[j].lower[2];
		bDS[k][3] = .5 * species[j].lower[3];
		bDS[k][4] = 1./3. * species[j].lower[4];
		bDS[k][5] = .25 * species[j].lower[5];
		bDS[k][6] = species[j].lower[7];

		T1[k] = species[j].t_min;
		T2[k] = species[j].t_mean;
		T3[k] = species[j].t_max;
	}
}

void OpenSMOKE_CHEMKINInterpreter_ThermoData::ProcessThermoData(BzzVectorInt &site_indices, BzzVectorInt &bulk_indices)
{
	// Site
	{
		int number_of_species = site_indices.Size();

		site_aDH	=	new BzzVector[number_of_species+1];
		site_bDH	=   new BzzVector[number_of_species+1];
		site_aDS	=   new BzzVector[number_of_species+1];
		site_bDS	=   new BzzVector[number_of_species+1];
		site_CpHT	=	new BzzVector[number_of_species+1];
		site_CpLT	=	new BzzVector[number_of_species+1];

		for(int k=1;k<=number_of_species;k++)
		{
			ChangeDimensions(6,&site_aDH[k]);
			ChangeDimensions(6,&site_bDH[k]);
			ChangeDimensions(6,&site_aDS[k]);
			ChangeDimensions(6,&site_bDS[k]);
			ChangeDimensions(5,&site_CpHT[k]);
			ChangeDimensions(5,&site_CpLT[k]);
		}

		ChangeDimensions(number_of_species,&site_M);
		ChangeDimensions(number_of_species,&site_T1);
		ChangeDimensions(number_of_species,&site_T2);
		ChangeDimensions(number_of_species,&site_T3);

		for(int k=1;k<=number_of_species;k++)
		{
			int j = site_indices[k];

			site_M[k]=species[j].mw;
			double RGAS_su_PM = Constants::R_J_mol/(site_M[k]*1.e-3);

			site_CpHT[k][1] = species[j].upper[1]*RGAS_su_PM;		// [J/kgK]
			site_CpHT[k][2] = species[j].upper[2]*RGAS_su_PM;
			site_CpHT[k][3] = species[j].upper[3]*RGAS_su_PM;
			site_CpHT[k][4] = species[j].upper[4]*RGAS_su_PM;
			site_CpHT[k][5] = species[j].upper[5]*RGAS_su_PM;

			site_CpLT[k][1] = species[j].lower[1]*RGAS_su_PM;
			site_CpLT[k][2] = species[j].lower[2]*RGAS_su_PM;
			site_CpLT[k][3] = species[j].lower[3]*RGAS_su_PM;
			site_CpLT[k][4] = species[j].lower[4]*RGAS_su_PM;
			site_CpLT[k][5] = species[j].lower[5]*RGAS_su_PM;

			site_aDH[k][1] = species[j].upper[1];				// L'entalpia viene calcolata in forma
			site_aDH[k][2] = .5 * species[j].upper[2];			// adimensionale, cio?dividendo per Rgas e
			site_aDH[k][3] = 1./3. * species[j].upper[3];		// per la temperatura in modo tale da
			site_aDH[k][4] = .25 * species[j].upper[4];			// poter calcolare pi velocemente le costanti
			site_aDH[k][5] = .2 * species[j].upper[5];			// di equilibrio
			site_aDH[k][6] = species[j].upper[6];

			site_bDH[k][1] = species[j].lower[1];
			site_bDH[k][2] = .5 * species[j].lower[2];
			site_bDH[k][3] = 1./3. * species[j].lower[3];
			site_bDH[k][4] = .25 * species[j].lower[4];
			site_bDH[k][5] = .2 * species[j].lower[5];
			site_bDH[k][6] = species[j].lower[6];

			site_aDS[k][1] = species[j].upper[1];				// L' entropia viene calcolata in questo modo
			site_aDS[k][2] = species[j].upper[2];				// cio?adimensionale perch?ci?velocizza
			site_aDS[k][3] = .5 * species[j].upper[3];		// il calcolo delle costanti di equilibrio
			site_aDS[k][4] = 1./3. * species[j].upper[4];
			site_aDS[k][5] = .25 * species[j].upper[5];
			site_aDS[k][6] = species[j].upper[7];

			site_bDS[k][1] = species[j].lower[1];
			site_bDS[k][2] = species[j].lower[2];
			site_bDS[k][3] = .5 * species[j].lower[3];
			site_bDS[k][4] = 1./3. * species[j].lower[4];
			site_bDS[k][5] = .25 * species[j].lower[5];
			site_bDS[k][6] = species[j].lower[7];

			site_T1[k] = species[j].t_min;
			site_T2[k] = species[j].t_mean;
			site_T3[k] = species[j].t_max;
		}
	}

	// Bulk
	{
		int number_of_species = bulk_indices.Size();

		bulk_aDH	=	new BzzVector[number_of_species+1];
		bulk_bDH	=   new BzzVector[number_of_species+1];
		bulk_aDS	=   new BzzVector[number_of_species+1];
		bulk_bDS	=   new BzzVector[number_of_species+1];
		bulk_CpHT	=	new BzzVector[number_of_species+1];
		bulk_CpLT	=	new BzzVector[number_of_species+1];

		for(int k=1;k<=number_of_species;k++)
		{
			ChangeDimensions(6, &bulk_aDH[k]);
			ChangeDimensions(6, &bulk_bDH[k]);
			ChangeDimensions(6, &bulk_aDS[k]);
			ChangeDimensions(6, &bulk_bDS[k]);
			ChangeDimensions(5, &bulk_CpHT[k]);
			ChangeDimensions(5, &bulk_CpLT[k]);
		}

		ChangeDimensions(number_of_species, &bulk_M);
		ChangeDimensions(number_of_species, &bulk_T1);
		ChangeDimensions(number_of_species, &bulk_T2);
		ChangeDimensions(number_of_species, &bulk_T3);

		for(int k=1;k<=number_of_species;k++)
		{
			int j = site_indices[k];

			bulk_M[k]=species[j].mw;
			double RGAS_su_PM = Constants::R_J_mol/(bulk_M[k]*1.e-3);

			bulk_CpHT[k][1] = species[j].upper[1]*RGAS_su_PM;		// [J/kgK]
			bulk_CpHT[k][2] = species[j].upper[2]*RGAS_su_PM;
			bulk_CpHT[k][3] = species[j].upper[3]*RGAS_su_PM;
			bulk_CpHT[k][4] = species[j].upper[4]*RGAS_su_PM;
			bulk_CpHT[k][5] = species[j].upper[5]*RGAS_su_PM;

			bulk_CpLT[k][1] = species[j].lower[1]*RGAS_su_PM;
			bulk_CpLT[k][2] = species[j].lower[2]*RGAS_su_PM;
			bulk_CpLT[k][3] = species[j].lower[3]*RGAS_su_PM;
			bulk_CpLT[k][4] = species[j].lower[4]*RGAS_su_PM;
			bulk_CpLT[k][5] = species[j].lower[5]*RGAS_su_PM;

			bulk_aDH[k][1] = species[j].upper[1];				// L'entalpia viene calcolata in forma
			bulk_aDH[k][2] = .5 * species[j].upper[2];			// adimensionale, cio?dividendo per Rgas e
			bulk_aDH[k][3] = 1./3. * species[j].upper[3];		// per la temperatura in modo tale da
			bulk_aDH[k][4] = .25 * species[j].upper[4];			// poter calcolare pi velocemente le costanti
			bulk_aDH[k][5] = .2 * species[j].upper[5];			// di equilibrio
			bulk_aDH[k][6] = species[j].upper[6];

			bulk_bDH[k][1] = species[j].lower[1];
			bulk_bDH[k][2] = .5 * species[j].lower[2];
			bulk_bDH[k][3] = 1./3. * species[j].lower[3];
			bulk_bDH[k][4] = .25 * species[j].lower[4];
			bulk_bDH[k][5] = .2 * species[j].lower[5];
			bulk_bDH[k][6] = species[j].lower[6];

			bulk_aDS[k][1] = species[j].upper[1];				// L' entropia viene calcolata in questo modo
			bulk_aDS[k][2] = species[j].upper[2];				// cio?adimensionale perch?ci?velocizza
			bulk_aDS[k][3] = .5 * species[j].upper[3];		// il calcolo delle costanti di equilibrio
			bulk_aDS[k][4] = 1./3. * species[j].upper[4];
			bulk_aDS[k][5] = .25 * species[j].upper[5];
			bulk_aDS[k][6] = species[j].upper[7];

			bulk_bDS[k][1] = species[j].lower[1];
			bulk_bDS[k][2] = species[j].lower[2];
			bulk_bDS[k][3] = .5 * species[j].lower[3];
			bulk_bDS[k][4] = 1./3. * species[j].lower[4];
			bulk_bDS[k][5] = .25 * species[j].lower[5];
			bulk_bDS[k][6] = species[j].lower[7];

			bulk_T1[k] = species[j].t_min;
			bulk_T2[k] = species[j].t_mean;
			bulk_T3[k] = species[j].t_max;
		}
	}
}

void OpenSMOKE_CHEMKINInterpreter_ThermoData::WriteSurfaceThermodynamicData(BzzSave &binaryFile, BzzSave &asciiFile, BzzVectorInt &site_indices, BzzVectorInt &bulk_indices)
{
	{
		// Site species
		for(int k=1;k<=site_indices.Size();k++)
		{
			binaryFile << site_CpHT[k];
			binaryFile << site_CpLT[k];
			binaryFile << site_aDH[k];
			binaryFile << site_bDH[k];
			binaryFile << site_aDS[k];
			binaryFile << site_bDS[k];
		}
		binaryFile << site_M;
		binaryFile << site_T1;
		binaryFile << site_T2;
		binaryFile << site_T3;

		// Bulk species
		if (bulk_indices.Size() > 0)
		{
			for(int k=1;k<=bulk_indices.Size();k++)
			{
				binaryFile << bulk_CpHT[k];
				binaryFile << bulk_CpLT[k];
				binaryFile << bulk_aDH[k];
				binaryFile << bulk_bDH[k];
				binaryFile << bulk_aDS[k];
				binaryFile << bulk_bDS[k];
			}
			binaryFile << bulk_M;
			binaryFile << bulk_T1;
			binaryFile << bulk_T2;
			binaryFile << bulk_T3;
		}
	}

	{
		// Site species
		for(int k=1;k<=site_indices.Size();k++)
		{
			asciiFile << site_CpHT[k];
			asciiFile << site_CpLT[k];
			asciiFile << site_aDH[k];
			asciiFile << site_bDH[k];
			asciiFile << site_aDS[k];
			asciiFile << site_bDS[k];
		}
		asciiFile << site_M;
		asciiFile << site_T1;
		asciiFile << site_T2;
		asciiFile << site_T3;

		// Bulk species
		if (bulk_indices.Size() > 0)
		{
			for(int k=1;k<=bulk_indices.Size();k++)
			{
				asciiFile << bulk_CpHT[k];
				asciiFile << bulk_CpLT[k];
				asciiFile << bulk_aDH[k];
				asciiFile << bulk_bDH[k];
				asciiFile << bulk_aDS[k];
				asciiFile << bulk_bDS[k];
			}
			asciiFile << bulk_M;
			asciiFile << bulk_T1;
			asciiFile << bulk_T2;
			asciiFile << bulk_T3;
		}
	}
}

void OpenSMOKE_CHEMKINInterpreter_ThermoData::WriteSurfaceElementsData(BzzSave &binaryFile,BzzSave &asciiFile)
{
	{
		// Site
		binaryFile << site_elements_matrix;
		binaryFile << site_elements_mw;

		// Bulk
		if (bulk_elements_matrix.Columns() > 0)
		{
			binaryFile << bulk_elements_matrix;
			binaryFile << bulk_elements_mw;
		}
	}

	{
		// Site
		asciiFile << site_elements_matrix;
		asciiFile << site_elements_mw;

		// Bulk
		if (bulk_elements_matrix.Columns() > 0)
		{
			asciiFile << bulk_elements_matrix;
			asciiFile << bulk_elements_mw;
		}
	}
}

void OpenSMOKE_CHEMKINInterpreter_ThermoData::ProcessElementData(BzzVectorInt &indices, OpenSMOKE_CHEMKINInterpreter_ElementsData &elements)
{
	int i;
	int number_of_species  = indices.Size();
	int number_of_elements = elements.elements_in_list.size()-1;

	ChangeDimensions(number_of_elements, &elements_mw);
	ChangeDimensions(number_of_elements, number_of_species, &elements_matrix);

	for(i=1;i<=number_of_elements;i++)
	{
		list_of_elements.push_back(elements.elements_in_list[i]);
		elements_mw[i] = elements.elements_mw[elements.elements_in_list_indices[i]];
	}

	for(i=1;i<=number_of_species;i++)
	{
		int k = indices[i];	
		for(int j=1;j<=int(elements.elements_in_list.size())-1;j++)
			elements_matrix[j][i] = species[k].element_indices[j];
	}
}

void OpenSMOKE_CHEMKINInterpreter_ThermoData::ProcessElementData(BzzVectorInt &site_indices, BzzVectorInt &bulk_indices, OpenSMOKE_CHEMKINInterpreter_ElementsData &elements)
{
	// Sites
	{
		int number_of_species  = site_indices.Size();
		int number_of_elements = elements.elements_in_list.size()-1;

		ChangeDimensions(number_of_elements, &site_elements_mw);
		ChangeDimensions(number_of_elements, number_of_species, &site_elements_matrix);

		for(int i=1;i<=number_of_species;i++)
		{
			int k = site_indices[i];	
			for(int j=1;j<=int(elements.elements_in_list.size())-1;j++)
				site_elements_matrix[j][i] = species[k].element_indices[j];
		}
	}

	// Bulk
	{
		int number_of_species  = bulk_indices.Size();
		int number_of_elements = elements.elements_in_list.size()-1;

		ChangeDimensions(number_of_elements, &bulk_elements_mw);
		ChangeDimensions(number_of_elements, number_of_species, &bulk_elements_matrix);

		for(int i=1;i<=number_of_species;i++)
		{
			int k = bulk_indices[i];	
			for(int j=1;j<=int(elements.elements_in_list.size())-1;j++)
				bulk_elements_matrix[j][i] = species[k].element_indices[j];
		}
	}
}

double OpenSMOKE_CHEMKINInterpreter_ThermoData::ReducedMolecularWeight(const string name1, const string name2)
{
	vector<string> names_list;
	BzzVectorInt indices;
	names_list.push_back("names");
	names_list.push_back(name1);
	names_list.push_back(name2);
	indices = GiveMeSpeciesIndices(names_list);
	return   species[indices[1]].mw*species[indices[2]].mw /
			(species[indices[1]].mw+species[indices[2]].mw);
}

