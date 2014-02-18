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
#include "addons/OpenSMOKE_SensitivityAnalysis_Fast_Flame1D.h"

void OpenSMOKE_SensitivityAnalysis_Fast_Flame1D::ErrorMessage(const string message)
{
    cout << endl;
    cout << "Class:  OpenSMOKE_SensitivityAnalysis_Fast_Flame1D"	<< endl;
    cout << "Object: " << name_object						<< endl;
    cout << "Error:  " << message							<< endl;
    cout << "Press a key to continue... "					<< endl;
    getchar();
    exit(-1);
}

void OpenSMOKE_SensitivityAnalysis_Fast_Flame1D::WarningMessage(const string message)
{
    cout << endl;
    cout << "Class:  OpenSMOKE_SensitivityAnalysis_Fast_Flame1D"	<< endl;
    cout << "Object: "		<< name_object					<< endl;
    cout << "Warning:  "	<< message						<< endl;
    cout << "Press a key to continue... "					<< endl;
    getchar();
}

void OpenSMOKE_SensitivityAnalysis_Fast_Flame1D::SetName(const string name)
{
	name_object = name;
}

OpenSMOKE_SensitivityAnalysis_Fast_Flame1D::OpenSMOKE_SensitivityAnalysis_Fast_Flame1D()
{
	name_object					= "name not assigned";
	indexVelocity				= 0;			// index of velocity
	indexMassFlowRate			= 0;			// index of mass flow rate
	indexPressureCurvature		= 0;			// index of pressure curvature		
	indexTemperature			= 0;			// index of temperature
	indexSpecies				= 0;			// index of species
}

void OpenSMOKE_SensitivityAnalysis_Fast_Flame1D::Initialize(const int kind_of_flame, OpenSMOKE_ReactingGas *_mix, BzzVectorInt &_indices, const int _N, const kindOfSensitivityParameter _kindOfParameter)
{
	int i,k;

	mix				= _mix;

	NC = mix->NumberOfSpecies();				// Number of Species
	NR = mix->NumberOfReactions();				// Number of reactions
	NP = NR + mix->kinetics.numFallOff;			// Number of parameters
	N  = _N;									// Number of points
	cout << NC << endl;
	cout << NR << endl;
	cout << NP << endl;
	cout << N << endl;

	double MB_MAX = 256.;							// max RAM available
	NP_BLOCK = int(MB_MAX*1.e6/(8.*N*NC));			// maximum number of parameters per block
	N_BLOCKS = NP/NP_BLOCK;							// number of blocks
	NP_RESIDUAL = NP - N_BLOCKS*NP_BLOCK;			// number of residual parameters

	cout << NP_BLOCK << endl;
	cout << N_BLOCKS << endl;
	cout << NP_RESIDUAL << endl;

	if (NP_RESIDUAL == 0)
	{
		NP_BLOCK += 1;									// maximum number of parameters per block
		N_BLOCKS = NP/NP_BLOCK;							// number of blocks
		NP_RESIDUAL = NP - N_BLOCKS*NP_BLOCK;			// number of residual parameters
	}
	cout << NP_BLOCK << endl;
	cout << N_BLOCKS << endl;
	cout << NP_RESIDUAL << endl;

	kindOfParameter	= _kindOfParameter;

	if (kind_of_flame == 1)					// FLAME_SPEED
	{
		cout << "FLAME_SPEED" << endl;
		reacting_system	= "FLAME_SPEED";
		NV = NC+2;									// Number of variables
		NE = N*NV;									// Number of equations

		indexTemperature	= 1;			// index of temperature
		indexMassFlowRate	= 2;			// index of mass flow rate
		indexSpecies		= 3;			// index of species
	}
	else if (kind_of_flame == 2)			// PREMIXED
	{
		cout << "PREMIXED" << endl;
		reacting_system	= "PREMIXED";
		NV = NC+1;								// Number of variables
		NE = N*NV;								// Number of equations

		indexTemperature	= 1;				// index of temperature
		indexSpecies		= 2;				// index of species
	}
	else if (kind_of_flame == 3)			// OPPOSED
	{
		cout << "OPPOSED" << endl;
		reacting_system	= "OPPOSED";
		NV = NC+4;									// Number of variables
		NE = N*NV;									// Number of equations

		indexVelocity				= 1;			// index of velocity
		indexMassFlowRate			= 2;			// index of mass flow rate
		indexPressureCurvature		= 3;			// index of pressure curvature		
		indexTemperature			= 4;			// index of temperature
		indexSpecies				= 5;			// index of species
	}
	else if (kind_of_flame == 4)			// OPPOSED
	{
		cout << "PREMIXED_NOENERGY" << endl;
		reacting_system	= "PREMIXED_NOENERGY";
		NV = NC;									// Number of variables
		NE = N*NV;									// Number of equations

		indexVelocity				= 0;			// index of velocity
		indexMassFlowRate			= 0;			// index of mass flow rate
		indexPressureCurvature		= 0;			// index of pressure curvature		
		indexTemperature			= 0;			// index of temperature
		indexSpecies				= 1;			// index of species
	}
	else
		ErrorMessage("Not yet implemented");


	ChangeDimensions(NP,			 &parameters);
	ChangeDimensions(NV,   NP,		 &JAlfaPoint);
	ChangeDimensions(NE,   NP_BLOCK, &S);
	ChangeDimensions(NE,   NP_BLOCK, &JAlfa);

	M		= mix->M;							// Molecular weigths [kg/kmol]
	uM		= mix->uM;							// uMolecular weigths [kmol/kg]
		
	nu = new OpenSMOKE_NuManager[NC+1];
	for(i =1;i<=NC;i++)	
		nu[i].Set(i, mix->names[i]);
	BuildNuMatrix(&mix->kinetics);
	for(k=1;k<=NC;k++)
		nu[k].Clean();
	
	indices_print_species = _indices;

	threshold_normalization = 1.e-10;
}

void OpenSMOKE_SensitivityAnalysis_Fast_Flame1D::BuildNuMatrix(OpenSMOKE_Kinetics *kinetics)
{
	int i,k;

	int*	jD1 = kinetics->jDir1.GetHandle();
	int*	jD2 = kinetics->jDir2.GetHandle();
	int*	jD3 = kinetics->jDir3.GetHandle();
	int*	jD4 = kinetics->jDir4.GetHandle();
	int*	jD5 = kinetics->jDir5.GetHandle();
	double* vD5 = kinetics->valDir5.GetHandle();

	int*	jIT1 = kinetics->jInvTot1.GetHandle();
	int*	jIT2 = kinetics->jInvTot2.GetHandle();
	int*	jIT3 = kinetics->jInvTot3.GetHandle();
	int*	jIT4 = kinetics->jInvTot4.GetHandle();
	int*	jIT5 = kinetics->jInvTot5.GetHandle();
	double* vIT5 = kinetics->valInvTot5.GetHandle();
	
	for(i = 1;i <= NC;i++)
	{
		for(k = 1;k <= kinetics->numDir1[i];k++)
			nu[i].Set(-1., *jD1++);
		for(k = 1;k <= kinetics->numDir2[i];k++)
			nu[i].Set(-2., *jD2++);
		for(k = 1;k <= kinetics->numDir3[i];k++)
			nu[i].Set(-4., *jD3++);
		for(k = 1;k <= kinetics->numDir4[i];k++)
			nu[i].Set(-0.50, *jD4++);
		for(k = 1;k <= kinetics->numDir5[i];k++)
			nu[i].Set(-(*vD5++), *jD5++);

		for(k = 1;k <= kinetics->numInvTot1[i];k++)
			nu[i].Set(1., *jIT1++);
		for(k = 1;k <= kinetics->numInvTot2[i];k++)
			nu[i].Set(2., *jIT2++);
		for(k = 1;k <= kinetics->numInvTot3[i];k++)
			nu[i].Set(3., *jIT3++);
		for(k = 1;k <= kinetics->numInvTot4[i];k++)
			nu[i].Set(0.50, *jIT4++);
		for(k = 1;k <= kinetics->numInvTot5[i];k++)
			nu[i].Set(*vIT5++, *jIT5++);
	}
}

void OpenSMOKE_SensitivityAnalysis_Fast_Flame1D::BuildJAlfaMatrix(BzzVector &T, BzzVector &rho, BzzVector &Cp, vector<string> list_of_names,
																  const double P_Pascal, BzzMatrix &X)
{
	int S_NC = list_of_names.size();
	int S_NV = NV-NC+S_NC;
	int S_NE = S_NV*N;
	ChangeDimensions(S_NE, NP, &S_S);
	BzzVector xVector(NC);
	BzzVector cVector(NC);
	BzzVector RVector(NC);

	for(int n=1;n<=N_BLOCKS+1;n++)
	{
		int iStart;
		int iEnd;
		if (n<N_BLOCKS+1)
		{
			cout << "Internal block..." << endl;
			iStart	= NP_BLOCK*(n-1)+1;
			iEnd	= NP_BLOCK*n;
		}
		else
		{
			cout << "Last block..." << endl;
			iStart	= NP_BLOCK*(n-1)+1;
			iEnd	= NP;
			ChangeDimensions(NE, NP_RESIDUAL, &JAlfa);
			ChangeDimensions(NE, NP_RESIDUAL, &S);
		}

		cout << "Sensitivity Block #" << n << "/" << N_BLOCKS+1 << endl;
		
		JAlfa=0;
		for(int iPoint=1;iPoint<=N;iPoint++)
		{
			int index = (iPoint-1)*NV;

			JAlfaPoint = 0.;
			
			if (iPoint == 1 || iPoint == N)
			{
			}
			else
			{

				double cTot = P_Pascal  / (Constants::R_J_kmol*T[iPoint]);
				X.GetRow(iPoint,&xVector);
				cVector = cTot*xVector;

				mix->ComputeKineticParameters(T[iPoint], log(T[iPoint]), 1./T[iPoint], P_Pascal);
				mix->ComputeFromConcentrations(T[iPoint], cVector, cTot, &RVector);

				if (kindOfParameter == FREQUENCY_FACTOR)						mix->GiveMe_Jalfa_A(JAlfaPoint, nu, T[iPoint], indexSpecies, indexTemperature, parameters);
				else if (kindOfParameter == FREQUENCY_FACTOR_AND_TRANSPORT_PROPERTIES)	mix->GiveMe_Jalfa_A(JAlfaPoint, nu, T[iPoint], indexSpecies, indexTemperature, parameters);
				else ErrorMessage("Wrong Sensitivity Parameter");
				//if (kindOfParameter == BETA)				mix->GiveMe_Jalfa_Beta(JAlfaPoint, nu, T[iPoint], indexSpecies, indexTemperature, parameters);
				//if (kindOfParameter == ACTIVATION_ENERGY)	mix->GiveMe_Jalfa_Eatt(JAlfaPoint, nu, T[iPoint], indexSpecies, indexTemperature, parameters);

				for(int k=1;k<=NP;k++)
					if (parameters[k] <= 0.)	ErrorMessage("Parameters must be larger than zero!");
			}

			// Temperature equations
			if (Cp[1]<0) // this means no energy equation
			{
				for(int j=iStart;j<=iEnd;j++)
					JAlfa[index+indexTemperature][j-iStart+1] = 0;
			}
			else
			{
				for(int j=iStart;j<=iEnd;j++)
					JAlfa[index+indexTemperature][j-iStart+1] = JAlfaPoint[indexTemperature][j]/(rho[iPoint]*Cp[iPoint]);
			}

			// Species equations
			for(int i=1;i<=NC;i++)
				for(int j=iStart;j<=iEnd;j++)
					JAlfa[index+indexSpecies+i-1][j-iStart+1] = JAlfaPoint[indexSpecies+i-1][j]/rho[iPoint];
		}

		JAlfa *= -1.;
		Solve(AFactorized, JAlfa, &S);

		for(int k=iStart;k<=iEnd;k++)
		{
			int j;

			// First Rows 
			for(j=1;j<=indexSpecies-1;j++)
				for(int i=1;i<=N;i++)
					S_S[(i-1)*S_NV+j][k] = S[(i-1)*NV+j][k-iStart+1];

			// Species Rows
			for(j=1;j<=S_NC;j++)
			{
				int index = mix->recognize_species(list_of_names[j-1]);
				for(int i=1;i<=N;i++)
					S_S[(i-1)*S_NV+(j+indexSpecies-1)][k] = S[(i-1)*NV+(index+indexSpecies-1)][k-iStart+1];
			}
		}
	}
}


void OpenSMOKE_SensitivityAnalysis_Fast_Flame1D::VideoSummary()
{
	cout << endl;
	cout << " - Number of points:                                   "	<< N			<< endl;
	cout << " - Number of species:                                  "	<< NC			<< endl;
	cout << " - Number of reactions:                                "	<< NR			<< endl;
	cout << " - Number of fall-off reactions:                       "	<< mix->kinetics.numFallOff	<< endl;
	cout << " - Number of parameters:                               "	<< NP			<< endl;
	cout << " - Temperature index:                                  "	<< indexTemperature		<< endl;
	cout << " - Species index:                                      "	<< indexSpecies			<< endl;
	cout << " - Block dimension:                                    "	<< NV			<< endl;
	cout << " - Total number of variables:                          "	<< NE			<< endl;
	cout << " - Total number of sensitivity coefficients per point: "	<< NP*NV		<< endl;
	cout << " - Total number of sensitivity coefficients:           "	<< NP*NE		<< endl;
	cout << " - Parameter block dimension:                          "   << NP_BLOCK	<< endl;
	cout << " - Number of parameter blocks:                         "   << N_BLOCKS	<< endl;
	cout << " - Jacobian matrix dimension:                          "   << (NV*NV*3*N)* (8./1024000.)		<< " MB" << endl;
	cout << " - Sensitivity matrix dimension:                       "   << (NE*NP_BLOCK)* (8./1024000.)		<< " MB" << endl;
	cout << " - Jalfa matrix dimension:                             "   << (NE*NP_BLOCK)* (8./1024000.)		<< " MB" << endl;

}


void OpenSMOKE_SensitivityAnalysis_Fast_Flame1D::PrintOnFile(	const string nameFolderAdditionalData, BzzVector &x, BzzVector &T, BzzVector &H, BzzMatrix &omega, 
																BzzVector &U, BzzVector &G, BzzVector &rho, BzzVector &MWtot,
																vector<string> list_of_names)
{
/*	int j;
	string dummy;
	char name[Constants::NAME_SIZE];
	char name_reaction[Constants::REACTION_NAME_SIZE];

	BzzSave fOutput;
	fOutput('*', nameFolderAdditionalData + "/sensitivity.sen");

	dummy = "V20090713";
	strcpy(name, dummy.c_str());
	fOutput.fileSave.write((char*) name, sizeof(name));
	
	string building_date = GiveMeTimeAndDate();
	dummy = building_date;
	strcpy(name, dummy.c_str());
	fOutput.fileSave.write((char*) name, sizeof(name));

	dummy = reacting_system;
	strcpy(name, dummy.c_str());
	fOutput.fileSave.write((char*) name, sizeof(name));

	dummy = "NP";
	strcpy(name, dummy.c_str());
	fOutput.fileSave.write((char*) name, sizeof(name));
	fOutput << NP;

	dummy = "NV";
	strcpy(name, dummy.c_str());
	fOutput.fileSave.write((char*) name, sizeof(name));
	fOutput << NV;

	dummy = "NC";
	strcpy(name, dummy.c_str());
	fOutput.fileSave.write((char*) name, sizeof(name));
	fOutput << NC;

	dummy = "NR";
	strcpy(name, dummy.c_str());
	fOutput.fileSave.write((char*) name, sizeof(name));
	fOutput << NP;

	dummy = "N";
	strcpy(name, dummy.c_str());
	fOutput.fileSave.write((char*) name, sizeof(name));
	fOutput << N;

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
	for(j=1;j<=NR;j++)
	{
		if (mix->kinetics.IsAFallOffReaction(j)!=0)
		{
			strcpy(name_reaction, "(0)");
			strcat(name_reaction, mix->kinetics.reactionRates->strReaction[j]);
		}
		else
			strcpy(name_reaction, mix->kinetics.reactionRates->strReaction[j]);
		fOutput.fileSave.write((char*) name_reaction, sizeof(name_reaction));
	}
	for(j=1;j<=mix->kinetics.numFallOff;j++)
	{
		int k=mix->kinetics.iFallOff[j];
		strcpy(name_reaction, "(inf)");
		strcat(name_reaction, mix->kinetics.reactionRates->strReaction[k]);
		fOutput.fileSave.write((char*) name_reaction, sizeof(name_reaction));
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
*/
/*	int S_NV = NV-NC+S_NC;
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
*/
/*	dummy = "S_S";
	strcpy(name, dummy.c_str());
	fOutput.fileSave.write((char*) name, sizeof(name));
	fOutput << S_S;

	fOutput.End();
*/
}

void OpenSMOKE_SensitivityAnalysis_Fast_Flame1D::SaveOnBinaryFile(BzzSave &fOutput)
{
	string dummy;
	char name[Constants::NAME_SIZE];
	char name_reaction[Constants::REACTION_NAME_SIZE];

	dummy = "SENSITIVITY";
	strcpy(name, dummy.c_str());
	fOutput.fileSave.write((char*) name, sizeof(name));

	dummy = "V20100417";
	strcpy(name, dummy.c_str());
	fOutput.fileSave.write((char*) name, sizeof(name));

	if (reacting_system	== "PREMIXED")		dummy = "PREMIXED";
	if (reacting_system	== "FLAME_SPEED")	dummy = "FLAMESPEED";
	if (reacting_system	== "OPPOSED")		dummy = "OPPOSED";
	if (reacting_system	== "TWIN")			dummy = "TWIN";
	strcpy(name, dummy.c_str());
	fOutput.fileSave.write((char*) name, sizeof(name));

	dummy = "INDICES";
	strcpy(name, dummy.c_str());
	fOutput.fileSave.write((char*) name, sizeof(name));
	fOutput << indices_print_species;

	dummy = "NP";
	strcpy(name, dummy.c_str());
	fOutput.fileSave.write((char*) name, sizeof(name));
	fOutput << NP;

	dummy = "NV";
	strcpy(name, dummy.c_str());
	fOutput.fileSave.write((char*) name, sizeof(name));
	fOutput << NV;

	dummy = "N";
	strcpy(name, dummy.c_str());
	fOutput.fileSave.write((char*) name, sizeof(name));
	fOutput << N;

	dummy = "PARAMETERS";
	strcpy(name, dummy.c_str());
	fOutput.fileSave.write((char*) name, sizeof(name));
	fOutput << parameters;

	dummy = "FALLOFF";
	strcpy(name, dummy.c_str());
	fOutput.fileSave.write((char*) name, sizeof(name));
	for(int j=1;j<=mix->kinetics.numFallOff;j++)
	{
		int k=mix->kinetics.iFallOff[j];
		strcpy(name_reaction, "(inf)");
		strcat(name_reaction, mix->kinetics.reactionRates->strReaction[k]);
		fOutput.fileSave.write((char*) name_reaction, sizeof(name_reaction));
	}

	dummy = "S_S";
	strcpy(name, dummy.c_str());
	fOutput.fileSave.write((char*) name, sizeof(name));
	fOutput << S_S;
}

void OpenSMOKE_SensitivityAnalysis_Fast_Flame1D::SaveOnXMLFile(const std::string folder_name)
{
	unsigned int number_of_variables = 0;
	
	if (indexTemperature!=0)		number_of_variables++;
	if (indexMassFlowRate!=0)		number_of_variables++;
	if (indexVelocity!=0)			number_of_variables++;
	if (indexPressureCurvature!=0)	number_of_variables++;
	if (indexSpecies!=0)			number_of_variables+=indices_print_species.Size();

	unsigned int number_of_points = S_S.Rows()/number_of_variables;

	// Main file
	{
		std::string file_name = folder_name + "\\Sensitivities.xml";
		ofstream fXML;
		fXML.open(file_name.c_str(), ios::out);
		fXML.setf(ios::scientific);

		fXML << "<?xml version=\"1.0\" encoding=\"utf-8\"?>" << endl;
		fXML << "<opensmoke version=\"0.1a\">" << endl;
		
		fXML << "<variables>" << endl;
		if (indexMassFlowRate==0)
			fXML << indices_print_species.Size()+1 << endl;
		else
			fXML << indices_print_species.Size()+2 << endl;
		for(int i=1;i<=indices_print_species.Size();i++)
			fXML << mix->names[indices_print_species[i]] << " " << i-1 << " " << indices_print_species[i] << std::endl;
		fXML << "temperature" << " " << indices_print_species.Size() << " " << mix->NumberOfSpecies()+1 << std::endl;
		if (indexMassFlowRate!=0)
			fXML << "mass-flow-rate" << " " << indices_print_species.Size()+1 << " " << mix->NumberOfSpecies()+2 << std::endl;
		fXML << "</variables>" << endl;

		fXML << "<n-parameters>" << endl;
		fXML << parameters.Size() << endl;
		fXML << "</n-parameters>" << endl;
		fXML << "<points>" << endl;
		fXML << number_of_points << std::endl;
		fXML << "</points>" << endl;
		fXML << "<constant-parameters>" << std::endl;
		for(int i=1;i<=parameters.Size();i++)
			fXML << parameters[i] << endl;
		fXML << "</constant-parameters>" << std::endl;
		fXML << "</opensmoke>" << endl;

		fXML.close();
	}

	// Temperature
	if (indexTemperature!=0)
	{
		std::string file_name = folder_name + "\\Sensitivities.temperature.xml";
		ofstream fXML;
		fXML.open(file_name.c_str(), ios::out);
		fXML.setf(ios::scientific);

		fXML << "<?xml version=\"1.0\" encoding=\"utf-8\"?>" << endl;
		fXML << "<opensmoke version=\"0.1a\">" << endl;
		fXML << "<coefficients>" << endl;
		for (int j=1;j<=number_of_points;j++)
		{
			int k = (j-1)*number_of_variables + indexTemperature;
			for (int i=1;i<=parameters.Size();i++)
				fXML << S_S[k][i] << " ";
			fXML << endl;
		}
		fXML << "</coefficients>" << endl;
		fXML << "</opensmoke>" << endl;

		fXML.close();
	}

	// Species
	if (indexSpecies!=0)
	{
		for(int z=1;z<=indices_print_species.Size();z++)
		{
			std::string file_name = folder_name + "\\Sensitivities." + mix->names[indices_print_species[z]] + ".xml";
			ofstream fXML;
			fXML.open(file_name.c_str(), ios::out);
			fXML.setf(ios::scientific);

			fXML << "<?xml version=\"1.0\" encoding=\"utf-8\"?>" << endl;
			fXML << "<opensmoke version=\"0.1a\">" << endl;
			fXML << "<coefficients>" << endl;
			for (int j=1;j<=number_of_points;j++) 
			{
				int k = (j-1)*number_of_variables + indexSpecies + (z-1);
				for (int i=1;i<=parameters.Size();i++)
					fXML << S_S[k][i] << " ";
				fXML << endl;
			}
			fXML << "</coefficients>" << endl;
			fXML << "</opensmoke>" << endl;

			fXML.close();
		}
	}

	// Axial-velocity
	if (indexMassFlowRate!=0)
	{
		std::string file_name = folder_name + "\\Sensitivities.mass-flow-rate.xml";
		ofstream fXML;
		fXML.open(file_name.c_str(), ios::out);
		fXML.setf(ios::scientific);

		fXML << "<?xml version=\"1.0\" encoding=\"utf-8\"?>" << endl;
		fXML << "<opensmoke version=\"0.1a\">" << endl;
		fXML << "<coefficients>" << endl;
		for (int j=1;j<=number_of_points;j++)
		{
			int k = (j-1)*number_of_variables + indexMassFlowRate;
			for (int i=1;i<=parameters.Size();i++)
				fXML << S_S[k][i] << " ";
			fXML << endl;
		}
		fXML << "</coefficients>" << endl;
		fXML << "</opensmoke>" << endl;

		fXML.close();
	}
}