/***************************************************************************
 *   Copyright (C) 2006-2008 by Alberto Cuoci   	                       *
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
#include "basic/OpenSMOKE_Utilities.h"
#include "basic/OpenSMOKE_Constants.h"
#include "engine/OpenSMOKE_ReactingGas.h"
#include "addons/OpenSMOKE_SensitivityAnalysis_Diffusion_Flame1D.h"

void OpenSMOKE_SensitivityAnalysis_Diffusion_Flame1D::ErrorMessage(const std::string message)
{
    cout << endl;
    cout << "Class:  OpenSMOKE_SensitivityAnalysis_Diffusion_Flame1D"	<< endl;
    cout << "Object: " << name_object						<< endl;
    cout << "Error:  " << message							<< endl;
    cout << "Press a key to continue... "					<< endl;
    getchar();
    exit(-1);
}

void OpenSMOKE_SensitivityAnalysis_Diffusion_Flame1D::WarningMessage(const std::string message)
{
    cout << endl;
    cout << "Class:  OpenSMOKE_SensitivityAnalysis_Diffusion_Flame1D"	<< endl;
    cout << "Object: "		<< name_object					<< endl;
    cout << "Warning:  "	<< message						<< endl;
    cout << "Press a key to continue... "					<< endl;
    getchar();
}

void OpenSMOKE_SensitivityAnalysis_Diffusion_Flame1D::SetName(const std::string name)
{
	name_object = name;
}

OpenSMOKE_SensitivityAnalysis_Diffusion_Flame1D::OpenSMOKE_SensitivityAnalysis_Diffusion_Flame1D()
{
	name_object					= "name not assigned";
	indexVelocity				= 0;			// index of velocity
	indexMassFlowRate			= 0;			// index of mass flow rate
	indexPressureCurvature		= 0;			// index of pressure curvature		
	indexTemperature			= 0;			// index of temperature
	indexSpecies				= 0;			// index of species
}

void OpenSMOKE_SensitivityAnalysis_Diffusion_Flame1D::Initialize(const int kind_of_flame, OpenSMOKE_ReactingGas *_mix, BzzVectorInt &_indices, const int _N)
{
	mix				= _mix;

	NC = mix->NumberOfSpecies();				// Number of Species
	NP = NC;									// Number of parameters
	N  = _N;									// Number of points


	if (kind_of_flame == 1)						// FLAME_SPEED
	{
		reacting_system	= "FLAME_SPEED";
		NV = NC+2;								// Number of variables
		NE = N*NV;								// Number of equations

		indexTemperature		= 1;			// index of temperature
		indexMassFlowRate		= 2;			// index of mass flow rate
		indexSpecies			= 3;			// index of species
	}
	else if (kind_of_flame == 2)				// PREMIXED
	{
		reacting_system	= "PREMIXED";
		NV = NC+1;								// Number of variables
		NE = N*NV;								// Number of equations

		indexTemperature	= 1;				// index of temperature
		indexSpecies		= 2;				// index of species
	}
	else if (kind_of_flame == 3)				// OPPOSED
	{
		reacting_system	= "OPPOSED";
		NV = NC+4;								// Number of variables
		NE = N*NV;								// Number of equations

		indexVelocity				= 1;		// index of velocity
		indexMassFlowRate			= 2;		// index of mass flow rate
		indexPressureCurvature		= 3;		// index of pressure curvature		
		indexTemperature			= 4;		// index of temperature
		indexSpecies				= 5;		// index of species
	}
	else
		ErrorMessage("Not yet implemented");

	indices_print_species = _indices;
}


void OpenSMOKE_SensitivityAnalysis_Diffusion_Flame1D::PrintOnFile(	BzzVector &x, BzzVector &T, BzzVector &H, BzzMatrix &omega, 
																	BzzVector &U, BzzVector &G, BzzVector &rho, BzzVector &MWtot,
																	vector<string> list_of_names, BzzMatrix &S, BzzMatrix &Diffusivity)
{
	int j;
	std::string dummy;
	char name[Constants::NAME_SIZE];
	char name_reaction[Constants::REACTION_NAME_SIZE];

	BzzVector parameters(NP);
	for(int i=1;i<=Diffusivity.Columns();i++)
	{
		BzzVector aux = Diffusivity.GetColumn(i);
		parameters[i] = aux.GetSumElements()/double(Diffusivity.Columns());
		cout << i << " " << parameters[i] << endl;
	}

	BzzSave fOutput;
	fOutput('*', "sensitivity.bin");

	dummy = "V20090713";
	strcpy(name, dummy.c_str());
	fOutput.fileSave.write((char*) name, sizeof(name));
	cout << dummy << endl;
	
	std::string building_date = GiveMeTimeAndDate();
	dummy = building_date;
	strcpy(name, dummy.c_str());
	fOutput.fileSave.write((char*) name, sizeof(name));
	cout << dummy << endl;

	dummy = reacting_system;
	strcpy(name, dummy.c_str());
	fOutput.fileSave.write((char*) name, sizeof(name));
	cout << dummy << endl;

	dummy = "NP";
	strcpy(name, dummy.c_str());
	fOutput.fileSave.write((char*) name, sizeof(name));
	fOutput << NP;
	cout << dummy << endl;

	dummy = "NV";
	strcpy(name, dummy.c_str());
	fOutput.fileSave.write((char*) name, sizeof(name));
	fOutput << NV;
	cout << dummy << " " << NV << endl;

	dummy = "NC";
	strcpy(name, dummy.c_str());
	fOutput.fileSave.write((char*) name, sizeof(name));
	fOutput << NC;
	cout << dummy << " " << NC << endl;

	dummy = "NR";
	strcpy(name, dummy.c_str());
	fOutput.fileSave.write((char*) name, sizeof(name));
	fOutput << NP;
	cout << dummy << " " << NP << endl;

	dummy = "N";
	strcpy(name, dummy.c_str());
	fOutput.fileSave.write((char*) name, sizeof(name));
	fOutput << N;
	cout << dummy << " " << N << endl;

	dummy = "GRID";
	strcpy(name, dummy.c_str());
	fOutput.fileSave.write((char*) name, sizeof(name));
	fOutput << x;

	dummy = "T";
	strcpy(name, dummy.c_str());
	fOutput.fileSave.write((char*) name, sizeof(name));
	fOutput << T;

	dummy = "H";
	strcpy(name, dummy.c_str());
	fOutput.fileSave.write((char*) name, sizeof(name));
	fOutput << H;

	dummy = "OMEGA";
	strcpy(name, dummy.c_str());
	fOutput.fileSave.write((char*) name, sizeof(name));
	fOutput << omega;

	dummy = "U";
	strcpy(name, dummy.c_str());
	fOutput.fileSave.write((char*) name, sizeof(name));
	fOutput << U;

	dummy = "G";
	strcpy(name, dummy.c_str());
	fOutput.fileSave.write((char*) name, sizeof(name));
	fOutput << G;

	dummy = "RHO";
	strcpy(name, dummy.c_str());
	fOutput.fileSave.write((char*) name, sizeof(name));
	fOutput << rho;

	dummy = "Mtot";
	strcpy(name, dummy.c_str());
	fOutput.fileSave.write((char*) name, sizeof(name));
	fOutput << MWtot;

	dummy = "M";
	strcpy(name, dummy.c_str());
	fOutput.fileSave.write((char*) name, sizeof(name));
	fOutput << mix->M;

	dummy = "PARAMETERS";
	strcpy(name, dummy.c_str());
	fOutput.fileSave.write((char*) name, sizeof(name));
	fOutput << parameters;

	dummy = "SPECIES";
	strcpy(name, dummy.c_str());
	fOutput.fileSave.write((char*) name, sizeof(name));
	for(j=1;j<=NC;j++)
	{
		char name[Constants::NAME_SIZE];
		strcpy(name, mix->names[j].c_str());
		fOutput.fileSave.write((char*) name, sizeof(name));
	}

	dummy = "REACTIONS";
	strcpy(name, dummy.c_str());
	fOutput.fileSave.write((char*) name, sizeof(name));
	for(j=1;j<=NC;j++)
	{
		char name[Constants::REACTION_NAME_SIZE];
		strcpy(name, mix->names[j].c_str());
		fOutput.fileSave.write((char*) name, sizeof(name));
	}

	int S_NC = list_of_names.size();
	dummy = "S_NC";
	strcpy(name, dummy.c_str());
	fOutput.fileSave.write((char*) name, sizeof(name));
	fOutput << S_NC;

	dummy = "S_SPECIES";
	strcpy(name, dummy.c_str());
	fOutput.fileSave.write((char*) name, sizeof(name));
	for(j=1;j<=S_NC;j++)
	{
		char name[Constants::NAME_SIZE];
		strcpy(name, list_of_names[j-1].c_str());
		fOutput.fileSave.write((char*) name, sizeof(name));
	}

	int S_NV = NV-NC+S_NC;
	int S_NE = S_NV*N;
	BzzMatrix S_S(S_NE,NP);
	
	// First Rows 
	for(j=1;j<=indexSpecies-1;j++)
		for(int i=1;i<=N;i++)
			S_S.SetRow((i-1)*S_NV+j, S.GetRow((i-1)*NV+j));

	// Species Rows
	for(j=1;j<=S_NC;j++)
	{
		int index = mix->recognize_species(list_of_names[j-1]);
		for(int i=1;i<=N;i++)
			S_S.SetRow((i-1)*S_NV+(j+indexSpecies-1), S.GetRow( (i-1)*NV + (index+indexSpecies-1) ) );
	}

	dummy = "S_S";
	strcpy(name, dummy.c_str());
	fOutput.fileSave.write((char*) name, sizeof(name));
	fOutput << S_S;

	fOutput.End();

	S_S.BzzPrint("S_S");
}
