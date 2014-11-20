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
#include <iomanip>
#include <sstream>
#include "basic/OpenSMOKE_Utilities.h"
#include "basic/OpenSMOKE_Constants.h"
#include "engine/OpenSMOKE_ReactingGas.h"
#include "addons/OpenSMOKE_SensitivityAnalysis_Fast_Unsteady0D.h"

void OpenSMOKE_SensitivityAnalysis_Fast_Unsteady0D::ErrorMessage(const std::string message)
{
    cout << endl;
    cout << "Class:  OpenSMOKE_SensitivityAnalysis_Fast_Unsteady0D"	<< endl;
    cout << "Object: " << name_object						<< endl;
    cout << "Error:  " << message							<< endl;
    cout << "Press a key to continue... "					<< endl;
    getchar();
    exit(-1);
}

void OpenSMOKE_SensitivityAnalysis_Fast_Unsteady0D::WarningMessage(const std::string message)
{
    cout << endl;
    cout << "Class:  OpenSMOKE_SensitivityAnalysis_Fast_Unsteady0D"	<< endl;
    cout << "Object: "		<< name_object					<< endl;
    cout << "Warning:  "	<< message						<< endl;
    cout << "Press a key to continue... "					<< endl;
    getchar();
}

void OpenSMOKE_SensitivityAnalysis_Fast_Unsteady0D::SetName(const std::string name)
{
	name_object = name;
}

OpenSMOKE_SensitivityAnalysis_Fast_Unsteady0D::OpenSMOKE_SensitivityAnalysis_Fast_Unsteady0D()
{
	name_object					= "name not assigned";
	indexTemperature			= 0;			// index of temperature
	indexSpecies				= 0;			// index of species
	indexPressure				= 0;			// index of pressure
	indexHeatRelease			= 0;			// index of heat release
	additional_variables		= 0;			// index of additional variables
}

void OpenSMOKE_SensitivityAnalysis_Fast_Unsteady0D::Initialize(const kind_of_sensitivity_analysis kind_, const kind_of_integration integration_, OpenSMOKE_ReactingGas *_mix, BzzVectorInt &_indices, const std::string path_name)
{
	int i,k;

	mix				= _mix;

	// Species to save
	_indices_print_species = _indices;

	// Dimensions
	NC     = mix->NumberOfSpecies();				// Number of Species
	NR     = mix->NumberOfReactions();				// Number of reactions
	NP     = NR + mix->kinetics.numFallOff;			// Number of parameters
	NC_RED = _indices_print_species.Size();

	// Integration
	integration = integration_;

	// ISOTHERMAL PFR
	if (kind_ == OPENSMOKE_SENSITIVITYANALYSIS_FAST_UNSTEADY0D_ONLY_SPECIES)					
	{
		reacting_system	= kind_;
		additional_variables = 2;
		indexSpecies		 = 1;			// index of species
	}
	
	// NON-ISOTHERMAL PFR
	else if (kind_ == OPENSMOKE_SENSITIVITYANALYSIS_FAST_UNSTEADY0D_SPECIES_PLUS_TEMPERATURE)					
	{
		reacting_system	= kind_;
		additional_variables = 3;
		indexSpecies		 = 1;			// index of species
		indexTemperature	 = NC+1;		// index of temperature
	}

	// ICEM MULTIZONE
	else if (kind_ == OPENSMOKE_SENSITIVITYANALYSIS_FAST_UNSTEADY0D_ICEM_MULTIZONE)					
	{
		reacting_system	= kind_;
		additional_variables = 4;
		indexSpecies		 = 1;			// index of species
		indexTemperature	 = NC+1;			// index of temperature
		indexPressure		 = NC+3;			// index of temperature
		indexHeatRelease	 = NC+4;			// index of temperature
	}
	
	else
		ErrorMessage("Not yet implemented");

	NV     = NC+additional_variables;							// Number of variables
	NV_RED = NC_RED+additional_variables;


	// Memory Allocation
	ChangeDimensions(NP,     &_parameters);
	ChangeDimensions(NV, NP, &_S);
	ChangeDimensions(NV, NP, &_JAlfa);
	ChangeDimensions(NV_RED, NP, &_S_S);
	
	
	if (integration == OPENSMOKE_SENSITIVITYANALYSIS_FAST_UNSTEADY0D_EULER_EXPLICIT)
	{
		ChangeDimensions(NV, NV, &_A_Euler_Explicit);
		ChangeDimensions(NV,     &_I);
		ChangeDimensions(NV, NP, &_Sold);
		_I = 1.;	
	}
	else if (integration == OPENSMOKE_SENSITIVITYANALYSIS_FAST_UNSTEADY0D_EULER_IMPLICIT)
	{
		ChangeDimensions(NV, NV, &_A_Euler_Implicit);
		ChangeDimensions(NV,     &_I);
		ChangeDimensions(NV, NP, &_Sold);
		_I = 1.;	
	}
	else if (integration == OPENSMOKE_SENSITIVITYANALYSIS_FAST_UNSTEADY0D_ACCURATE)
	{
		ChangeDimensions(NV, NP, &_dS_over_dt);
	}

	// Build nu matrix
	nu = new OpenSMOKE_NuManager[NC+1];
	for(i =1;i<=NC;i++)	
		nu[i].Set(i, mix->names[i]);
	BuildNuMatrix(&mix->kinetics);
	for(k=1;k<=NC;k++)
		nu[k].Clean();

	// Normalization threshold
	threshold_normalization = 1.e-10;

	// Time
	tOld = 0.;

	// Summary
	VideoSummary();
	std::string message = "mkdir " + path_name;
    system(message.c_str());
	PrepareOutputFile(path_name + "\\Sensitivity.log");
	fBinary('*', path_name + "\\Sensitivity.bin");
}

void OpenSMOKE_SensitivityAnalysis_Fast_Unsteady0D::BuildNuMatrix(OpenSMOKE_Kinetics *kinetics)
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


void OpenSMOKE_SensitivityAnalysis_Fast_Unsteady0D::BuildJAlfaMatrix(const double T, const double P_Pascal, 
																	 const double rho, const double denominator, const double v,
																	 BzzVector &X)
{
	BzzVector cVector(NC);
	BzzVector RVector(NC);

	// Build Jalfa matrix
	_JAlfa=0.;
	{
		double cTot = P_Pascal  / (Constants::R_J_kmol*T);
		cVector = cTot*X;

		mix->ComputeKineticParameters(T, log(T), 1./T, P_Pascal);
		mix->ComputeFromConcentrations(T, cVector, cTot, &RVector);

		mix->GiveMe_Jalfa_A(_JAlfa, nu, T, indexSpecies, indexTemperature, _parameters);
	}

	// Temperature equations
	if (indexTemperature > 0)
		for(int j=1;j<=NP;j++)
			_JAlfa[indexTemperature][j] /= (rho*denominator);	// (Cp+v2/T) if PFR, Cv if ICEM

	if (indexHeatRelease > 0)
		for(int j=1;j<=NP;j++)
			_JAlfa[indexTemperature][j] *= v;			// in this case v is the volume

	// Species equations
	for(int i=1;i<=NC;i++)
		for(int j=1;j<=NP;j++)
			_JAlfa[indexSpecies+i-1][j] /= rho;
}

void OpenSMOKE_SensitivityAnalysis_Fast_Unsteady0D::Update(const double tCurrent, BzzMatrix &J,
														   const double T, const double P_Pascal, 
														   const double rho, const double denominator, const double v, BzzVector &X)
{
	double deltat = tCurrent - tOld;
	tOld = tCurrent;

	if (_JOld.Rows() == 0)	_JOld = J;

	BuildJAlfaMatrix(T, P_Pascal, rho, denominator, v, X);

	if (integration == OPENSMOKE_SENSITIVITYANALYSIS_FAST_UNSTEADY0D_EULER_EXPLICIT)
		Euler_Explicit(deltat, J);

	else if (integration == OPENSMOKE_SENSITIVITYANALYSIS_FAST_UNSTEADY0D_EULER_IMPLICIT)
		Euler_Implicit(deltat, J);

	UpdateOutputFile(tCurrent);

	SaveOnBinaryFile(fBinary);

	// Save old values
	_JOld = J;
}

void OpenSMOKE_SensitivityAnalysis_Fast_Unsteady0D::Euler_Explicit(const double deltat, BzzMatrix &J)
{
	_A_Euler_Explicit = deltat*J;			//	A = J*deltat
	Sum(_I, &_A_Euler_Explicit);			//  A = I + J*deltat
	Product(_A_Euler_Explicit, &_Sold);		//  Sold = (I+J*deltat)*Sold
	_S = deltat*_JAlfa;						//  S = Jalfa*deltat
	Sum(_Sold,&_S);							//  S = (I+J*deltat)*Sold + deltat*JAlfa
	_Sold = _S;
}

void OpenSMOKE_SensitivityAnalysis_Fast_Unsteady0D::Euler_Implicit(const double deltat, BzzMatrix &J)
{
	_S = deltat*_JAlfa;							//  S = Jalfa*deltat
	Sum(_S,&_Sold);								//  Sold = Jalfa*deltat + Sold
	_A_Euler_Implicit = -deltat*J;				//  A = -J*deltat
	Sum(_I, &_A_Euler_Implicit);				//	A = I-J*deltat
	ReplaceBzzMatrixWithBzzFactorized(&_A_Euler_Implicit, &_A_Euler_Implicit_Factorized);
	Solve(_A_Euler_Implicit_Factorized, _Sold, &_S);		//  
	_Sold = _S;
}

void OpenSMOKE_SensitivityAnalysis_Fast_Unsteady0D::VideoSummary()
{
	cout << endl;
	cout << " - Number of species:                                  "	<< NC						<< endl;
	cout << " - Number of reactions:                                "	<< NR						<< endl;
	cout << " - Number of fall-off reactions:                       "	<< mix->kinetics.numFallOff	<< endl;
	cout << " - Number of parameters:                               "	<< NP						<< endl;
	cout << " - Temperature index:                                  "	<< indexTemperature			<< endl;
	cout << " - Pressure index:                                     "	<< indexPressure			<< endl;
	cout << " - Heat release index:                                 "	<< indexHeatRelease			<< endl;
	cout << " - Species index:                                      "	<< indexSpecies				<< endl;
	cout << " - Block dimension:                                    "	<< NV						<< endl;
	cout << " - Total number of variables:                          "	<< NV						<< endl;
	cout << " - Total number of sensitivity coefficients per point: "	<< NP*NV					<< endl;
	cout << " - Total number of sensitivity coefficients:           "	<< NP*NV					<< endl;
}

	
void OpenSMOKE_SensitivityAnalysis_Fast_Unsteady0D::PrepareOutputFile(const std::string file_name)
{
	openOutputFileAndControl(fOutput, file_name);
	fOutput.setf(ios::scientific);

	fOutput << setw(14) << left << "time[s]";

	int count = 2;

	if (indexTemperature > 0)
		for(int k=1;k<=NP;k++)
		{
			stringstream number_k;		number_k << k;
			stringstream number_count;	number_count << count++;
			std::string label = "S_T_"+ number_k.str() + "(" + number_count.str() + ")";
			fOutput << setw(16) << left << label;
		}

	if (indexPressure > 0)
		for(int k=1;k<=NP;k++)
		{
			stringstream number_k;		number_k << k;
			stringstream number_count;	number_count << count++;
			std::string label = "S_P_"+ number_k.str() + "(" + number_count.str() + ")";
			fOutput << setw(16) << left << label;
		}
	
	if (indexHeatRelease > 0)
		for(int k=1;k<=NP;k++)
		{
			stringstream number_k;		number_k << k;
			stringstream number_count;	number_count << count++;
			std::string label = "S_Q_"+ number_k.str() + "(" + number_count.str() + ")";
			fOutput << setw(16) << left << label;
		}

	for(int j=1;j<=_indices_print_species.Size();j++)
		for(int k=1;k<=NP;k++)
		{
			stringstream number_k;		number_k << k;
			stringstream number_count;	number_count << count++;
			std::string label = "S_" + mix->names[_indices_print_species[j]] + "_" + number_k.str() + "(" + number_count.str() + ")";
			fOutput << setw(16) << left << label;
		}

	fOutput << endl << endl;

}
void OpenSMOKE_SensitivityAnalysis_Fast_Unsteady0D::UpdateOutputFile(const double t)
{
	fOutput << setw(14) << left << t;

	if (indexTemperature > 0)
		for(int k=1;k<=NP;k++)
			fOutput << setw(16) << left << _S[indexTemperature][k];

	if (indexPressure > 0)
		for(int k=1;k<=NP;k++)
			fOutput << setw(16) << left << _S[indexPressure][k];

	if (indexHeatRelease > 0)
		for(int k=1;k<=NP;k++)
			fOutput << setw(16) << left << _S[indexHeatRelease][k];

	for(int j=1;j<=_indices_print_species.Size();j++)
		for(int k=1;k<=NP;k++)
			fOutput << setw(16) << left << _S[_indices_print_species[j]][k];
	
	fOutput << endl;
}

void OpenSMOKE_SensitivityAnalysis_Fast_Unsteady0D::UpdateOutputFile(const double t, BzzVector &x)
{
	// Reconstruct S from x
	int k=NV+1;
	for(int i=1;i<=NV;i++)
		for(int j=1;j<=NP;j++)
			_S[i][j] = x[k++];

	UpdateOutputFile(t);
}

void OpenSMOKE_SensitivityAnalysis_Fast_Unsteady0D::Update(BzzMatrix &J, const double T, const double P_Pascal, 
														   const double rho, const double denominator, const double v, BzzVector &X,
														   BzzVector &x, BzzVector &r)
{
	// Reconstruct S from x (x size is NV + NV*NP)
	int k=NV+1;
	for(int i=1;i<=NV;i++)
		for(int j=1;j<=NP;j++)
			_S[i][j] = x[k++];

	// Calculate
	BuildJAlfaMatrix(T, P_Pascal, rho, denominator, v, X);

	// dS_over_dt = J*S
	if (J.Rows()>0)
	{
		BzzMatrix JJ(NV, NV);
		for(int i=1;i<=NV;i++)
			for(int j=1;j<=NV;j++)
				JJ[i][j] = J[i][j];
		Product(JJ,_S,&_dS_over_dt);
	}
	else
		_dS_over_dt = 0.;

	// dS_over_dt = Jalfa + J*S
	Sum(_JAlfa,&_dS_over_dt);

	// Reconstruct r from dS_over_dt (r size is NV*NP)
	{
		r = 0.;
		int k=1;
		for(int i=1;i<=NV;i++)
			for(int j=1;j<=NP;j++)
				r[k++] = _dS_over_dt[i][j];
	}
}

void OpenSMOKE_SensitivityAnalysis_Fast_Unsteady0D::Close()
{
	cout << "closing binary file" << endl;
	CloseBinaryFile(fBinary);
}

void OpenSMOKE_SensitivityAnalysis_Fast_Unsteady0D::CloseBinaryFile(BzzSave &fOutput)
{
	std::string dummy;
	char name[Constants::NAME_SIZE];

	dummy = "END";
	strcpy(name, dummy.c_str());
	fOutput.fileSave.write((char*) name, sizeof(name));

	fOutput.End();
}

void OpenSMOKE_SensitivityAnalysis_Fast_Unsteady0D::SaveOnBinaryFile(BzzSave &fOutput)
{
	std::string dummy;
	char name[Constants::NAME_SIZE];

	Populate_S_S();

	dummy = "S_S";
	strcpy(name, dummy.c_str());
	fOutput.fileSave.write((char*) name, sizeof(name));
	fOutput << _S_S;
}

void OpenSMOKE_SensitivityAnalysis_Fast_Unsteady0D::PrepareBinaryFile(BzzSave &fOutput, const int N, BzzMatrix &S, const std::string tag)
{
	std::string dummy;
	char name[Constants::NAME_SIZE];
	char name_reaction[Constants::REACTION_NAME_SIZE];

	dummy = "SENSITIVITY";
	strcpy(name, dummy.c_str());
	fOutput.fileSave.write((char*) name, sizeof(name));

	dummy = "V20100417";
	strcpy(name, dummy.c_str());
	fOutput.fileSave.write((char*) name, sizeof(name));

	if ( reacting_system == OPENSMOKE_SENSITIVITYANALYSIS_FAST_UNSTEADY0D_ONLY_SPECIES || 
		 reacting_system == OPENSMOKE_SENSITIVITYANALYSIS_FAST_UNSTEADY0D_SPECIES_PLUS_TEMPERATURE )					
	{
		if (tag == "pfr-isothermal")			dummy = "PFR-ISOT";
		else if (tag == "pfr-non-isothermal")	dummy = "PFR-NONISOT";
		else ErrorMessage("Only pfr-isothermal or pfr-non-isothermal are currently allowed");

		strcpy(name, dummy.c_str());
		fOutput.fileSave.write((char*) name, sizeof(name));
	}

	else if (reacting_system == OPENSMOKE_SENSITIVITYANALYSIS_FAST_UNSTEADY0D_ICEM_MULTIZONE)					
	{
		dummy = "ICEM-MULTI";
		strcpy(name, dummy.c_str());
		fOutput.fileSave.write((char*) name, sizeof(name));
		fOutput << indexTemperature;
		fOutput << indexPressure;
		fOutput << indexHeatRelease;
	}

	dummy = "INDICES";
	strcpy(name, dummy.c_str());
	fOutput.fileSave.write((char*) name, sizeof(name));
	fOutput << _indices_print_species;

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
	fOutput << _parameters;

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
	fOutput << S;
}

void OpenSMOKE_SensitivityAnalysis_Fast_Unsteady0D::Populate_S_S()
{
	for(int i=1;i<=NC_RED;i++)
	{
		int k = _indices_print_species[i];
		for(int j=1;j<=NP;j++)
			_S_S[i][j] = _S[k][j];
	}

	if (indexTemperature > 0)
		for(int j=1;j<=NP;j++)
			_S_S[indexTemperature-NC+NC_RED][j] = _S[indexTemperature][j];

	if (indexHeatRelease > 0)
		for(int j=1;j<=NP;j++)
			_S_S[indexHeatRelease-NC+NC_RED][j] = _S[indexHeatRelease][j];

	if (indexPressure > 0)
		for(int j=1;j<=NP;j++)
			_S_S[indexPressure-NC+NC_RED][j] = _S[indexPressure][j];
}

void OpenSMOKE_SensitivityAnalysis_Fast_Unsteady0D::MoveFromFiles(const std::string nameFileLoadProvisional, const int N, BzzSave &fOutput, const std::string tag)
{	 
	Close();

	BzzMatrix S_S(NV_RED, NP);
	BzzMatrix S(NV_RED*N, NP);

	BzzLoad fLoad;
	fLoad('*', nameFileLoadProvisional);
	for(int i=1;i<=N;i++)
	{
		//cout << i << endl;
		CheckInBinaryFile(fLoad, "S_S");
		fLoad >> S_S;

		int index = (i-1)*NV_RED;
		for(int k=1;k<=NV_RED;k++)
			for(int j=1;j<=NP;j++)
				S[index+k][j] = S_S[k][j];
	}
	fLoad.End();

	PrepareBinaryFile(fOutput, N, S, tag);
}