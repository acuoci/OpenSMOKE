/***************************************************************************
 *   Copyright (C) 2010 by Alberto Cuoci         	                       *
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

#include <vector>
#include <iomanip>
#include <sstream>
#include "basic/OpenSMOKE_Utilities.h"
#include "basic/OpenSMOKE_Constants.h"
#include "addons/OpenSMOKE_PostProcessor.h"
#include "addons/OpenSMOKE_PostProcessor_Flame1D.h"
#include "addons/OpenSMOKE_PostProcessor_SensitivityAnalysis_Flame1D.h"

void OpenSMOKE_PostProcessor_SensitivityAnalysis_Flame1D::ErrorMessage(const std::string message)
{
    cout << endl;
    cout << "Class:  OpenSMOKE_PostProcessor_SensitivityAnalysis_Flame1D"	<< endl;
    cout << "Object: " << name_object								<< endl;
    cout << "Error:  " << message									<< endl;
    cout << "Press a key to continue... "							<< endl;
    getchar();
    exit(-1);
}

void OpenSMOKE_PostProcessor_SensitivityAnalysis_Flame1D::WarningMessage(const std::string message)
{
    cout << endl;
    cout << "Class:  OpenSMOKE_PostProcessor_SensitivityAnalysis_Flame1D"	<< endl;
    cout << "Object: "		<< name_object							<< endl;
    cout << "Warning:  "	<< message								<< endl;
    cout << "Press a key to continue... "							<< endl;
    getchar();
}

OpenSMOKE_PostProcessor_SensitivityAnalysis_Flame1D::OpenSMOKE_PostProcessor_SensitivityAnalysis_Flame1D(OpenSMOKE_PostProcessor *post_processor_) 
: OpenSMOKE_PostProcessor_SensitivityAnalysis_General(post_processor_)
{
}

void OpenSMOKE_PostProcessor_SensitivityAnalysis_Flame1D::ReadFromBinaryFile(BzzLoad &fLoad, const int index)
{
	int j;
	char dummy[Constants::NAME_SIZE];
	char name_reaction[Constants::REACTION_NAME_SIZE];

	fLoad.fileLoad.read((char*) dummy, sizeof(dummy));
	std::string version = dummy;
	if (version != "V20100417")
		ErrorMessage("This version post processing file is not supported: " + version);
	cout << "Version: " << version << endl;

	std::string tag = NextInBinaryFile(fLoad);
		 if (tag == "OPPOSED")		kind = OPPOSED;
	else if (tag == "FLAMESPEED")	kind = FLAMESPEED;
	else if (tag == "PREMIXED")		kind = PREMIXED;
	else if (tag == "TWIN")			kind = TWIN;
	else if (tag == "PFR-ISOT")		kind = PFR_ISOTHERMAL;
	else if (tag == "PFR-NONISOT")	kind = PFR_NONISOTHERMAL;
	else ErrorMessage("Post processing not available for this kind: " + tag);

	// Indices of species
	CheckInBinaryFile(fLoad, "INDICES");
	fLoad >> indexTransfer;
	S_NC = indexTransfer.Size();

	CheckInBinaryFile(fLoad, "NP");
	fLoad >> NP;
	cout << "Number of parameters: "  << NP << endl;

	CheckInBinaryFile(fLoad, "NV");
	fLoad >> NV;
	S_NV = NV-NC+S_NC;
	cout << "Number of variables per point: "  << NV << endl;

	CheckInBinaryFile(fLoad, "N");
	fLoad >> N;
	cout << "Number of points: "  << N << endl;

	CheckInBinaryFile(fLoad, "PARAMETERS");
	fLoad >> parameters;

	// Frequency factors
	if (index == 0)
	{
		CheckInBinaryFile(fLoad, "FALLOFF");
		reactions = new std::string[NP + 1];
		for(j=1;j<=NR;j++)		// Conventional (0)
			reactions[j] = post_processor->reactions[j];
		for(j=NR+1;j<=NP;j++)	// Fall-Off (inf)
		{
			fLoad.fileLoad.read((char*) name_reaction, sizeof(name_reaction));
			reactions[j] = name_reaction;
		}
	}
	// Diffusivity
	else if (index == 1)
	{
		nExtracted = Min(20,NP);
		reactions = new std::string[NP + 1];
		int NPtransport = NP/NC;
		for(int kk=1;kk<=NPtransport;kk++)
		{
			     if (kk==1) tag = "[e/k]";
			else if (kk==2) tag = "[sigma]";
			else if (kk==3) tag = "[mu]";
			else if (kk==4) tag = "[alfa]";
			else if (kk==5) tag = "[zRot298]";
			for(j=1;j<=NC;j++)		
				reactions[(kk-1)*NC+j] = post_processor->names[j]+tag;
		}

		for(j=1;j<=NP;j++)		
			cout << reactions[j] << endl;
	}


	CheckInBinaryFile(fLoad, "S_S");
	fLoad >> S;

	// Video information
	cout << endl;
	cout << "Number of species: "				<< NC			<< endl;
	cout << "Number of reactions: "				<< NR			<< endl;
	cout << "Number of species (sensitivity): " << S_NC			<< endl;
	cout << "Number of parameters: "			<< NP			<< endl;
	cout << "Number of points: "				<< N			<< endl;
	cout << "Block dimension: "					<< NV			<< endl;
	cout << "Block dimension (sensitivity): "	<< NV-NC+S_NC	<< endl;
	cout << endl;

	cout << " * Memory allocation... " << endl;
	ChangeDimensions(S_NV,		&vectorMax);
	ChangeDimensions(S_NV,		&vectorMin);
	ChangeDimensions(N, S_NV,	&matrixVariables);
	ChangeDimensions(N, NP,		&Slocal);
	ChangeDimensions(NP,		&SlocalVector);

	if (post_processor->kind == POST_PROCESSOR_FLAME1D)
	{
		if (kind == FLAMESPEED)
		{
			index_temperature	= 1;
			index_massflowrate	= 2;
			index_species		= 3;
			post_processor->startSensitivityAdditional=0;

			matrixVariables.SetColumn(index_temperature,  post_processor->get_flame().T);
			matrixVariables.SetColumn(index_massflowrate, post_processor->get_flame().H);
			for(int j=1;j<=S_NC;j++)
				matrixVariables.SetColumn(index_species+j-1, post_processor->get_flame().omega[indexTransfer[j]]);

			
		}

		else if (kind == PREMIXED)
		{
			index_temperature		= 1;
			index_species			= 2;
			post_processor->startSensitivityAdditional=0;

			matrixVariables.SetColumn(index_temperature, post_processor->get_flame().T);
			for(int j=1;j<=S_NC;j++)
				matrixVariables.SetColumn(index_species+j-1, post_processor->get_flame().omega[indexTransfer[j]]);
			
		}
		
		else if (kind == OPPOSED || kind == TWIN)
		{
			index_velocity			= 1;
			index_massflowrate		= 2;
			index_pressurecurvature	= 3;
			index_temperature		= 4;
			index_species			= 5;
			post_processor->startSensitivityAdditional=0;

			matrixVariables.SetColumn(index_velocity, post_processor->get_flame().u);
			matrixVariables.SetColumn(index_massflowrate, post_processor->get_flame().G);
			matrixVariables.SetColumn(index_pressurecurvature, post_processor->get_flame().H);
			matrixVariables.SetColumn(index_temperature, post_processor->get_flame().T);
			for(int j=1;j<=S_NC;j++)
				matrixVariables.SetColumn(index_species+j-1, post_processor->get_flame().omega[indexTransfer[j]]);
		}
	}

	if (post_processor->kind == POST_PROCESSOR_PFR)
	{
		if (kind == PFR_ISOTHERMAL)
		{
			index_species		= 1;
			post_processor->startSensitivityAdditional=0;

			for(int j=1;j<=S_NC;j++)
				matrixVariables.SetColumn(index_species+j-1, post_processor->get_omega_profile(indexTransfer[j]));
		}

		else if (kind == PFR_NONISOTHERMAL)
		{
			index_species			= 1;
			index_temperature		= S_NC+1;
			post_processor->startSensitivityAdditional=S_NC;
			
			cout << index_species << " " << index_temperature << endl;
			cout << matrixVariables.Rows() << " " << matrixVariables.Columns() << endl;
			cout << "omega" << endl;
			for(int j=1;j<=S_NC;j++)
				matrixVariables.SetColumn(index_species+j-1, post_processor->get_omega_profile(indexTransfer[j]));
			cout << "T" << endl;
			BzzVector cc;
			cc=post_processor->get_temperature_profile();
			cout << cc.Size() << endl;
			matrixVariables.SetColumn(index_temperature, post_processor->get_temperature_profile());
			cout << "done" << endl;
		}
	}

	// Populate Min and Max Vectors
	MinMaxCalculations();

	// Prepare Additional names
	PrepareAdditionalNames();

	// Prepare Additional names
	PrepareSensitivitySpeciesNames();
	
	// Prepare Additional data
	Prepare();
}

void OpenSMOKE_PostProcessor_SensitivityAnalysis_Flame1D::PrepareAdditionalNames()
{
	if (kind == OPPOSED || kind == TWIN)
	{
		additional_names.resize(4);
		additional_names[0] = "Velocity";
		additional_names[1] = "Mass flow rate";
		additional_names[2] = "Pressure curvature";
		additional_names[3] = "Temperature";
	}
	else if (kind == FLAMESPEED)
	{
		additional_names.resize(2);
		additional_names[0] = "Temperature";
		additional_names[1] = "Flame speed";
	}
	else if (kind == PREMIXED)
	{
		additional_names.resize(1);
		additional_names[0] = "Temperature";
	}
	else if (kind == PFR_ISOTHERMAL)
	{
	}
	else if (kind == PFR_NONISOTHERMAL)
	{
		additional_names.resize(1);
		additional_names[0] = "Temperature";
	}
}
