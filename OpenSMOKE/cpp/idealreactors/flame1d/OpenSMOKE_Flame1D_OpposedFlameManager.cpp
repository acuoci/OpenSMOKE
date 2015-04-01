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
#include "idealreactors/flame1d/OpenSMOKE_Flame1D_OpposedFlameManager.h" 

void OpenSMOKE_Flame1D_OpposedFlameManager::ErrorMessage(const std::string message)
{
    cout << endl;
    cout << "Class:  OpenSMOKE_Flame1D_OpposedFlameManager"	<< endl;
    cout << "Object: " << name_object			<< endl;
    cout << "Error:  " << message				<< endl;
    cout << "Press a key to continue... "		<< endl;
    getchar();
    exit(-1);
}

void OpenSMOKE_Flame1D_OpposedFlameManager::WarningMessage(const std::string message)
{
    cout << endl;
    cout << "Class:  OpenSMOKE_Flame1D_OpposedFlameManager"	<< endl;
    cout << "Object: "		<< name_object		<< endl;
    cout << "Warning:  "	<< message			<< endl;
    cout << "Press a key to continue... "		<< endl;
    getchar();
}

OpenSMOKE_Flame1D_OpposedFlameManager::OpenSMOKE_Flame1D_OpposedFlameManager()
{
	name_object	= "[Name not assigned]";
	iAssignedFuelVelocity			= false;
	iAssignedOxidizerVelocity		= false;
	iAssignedFuelTemperature		= false;
	iAssignedOxidizerTemperature	= false;
	iAssignedPressure				= false;
	N = 0;
}

void OpenSMOKE_Flame1D_OpposedFlameManager::SetName(const std::string name)
{
	name_object = name;
}

void OpenSMOKE_Flame1D_OpposedFlameManager::SetupFromFile(OpenSMOKE_Flame1D *_flame, OpenSMOKE_Flame1D_DataManager *_data)
{
	flame = _flame;
	data  = _data;

    std::string  string_value;
	vector<string> string_vector;
	vector<double> double_vector;

	OpenSMOKE_Dictionary_Flame1D_OpposedFlameManager dictionary;
    dictionary.ParseFile(data->opposedFlameAnalysisFileName);

	if (dictionary.Return("#FuelVelocity", string_vector))
        SetFuelVelocity(string_vector);

	if (dictionary.Return("#FuelTemperature", string_vector))
        SetFuelTemperature(string_vector);

	if (dictionary.Return("#OxidizerVelocity", string_vector))
        SetOxidizerVelocity(string_vector);

	if (dictionary.Return("#OxidizerTemperature", string_vector))
        SetOxidizerTemperature(string_vector);

	if (dictionary.Return("#Pressure", string_vector))
        SetPressure(string_vector);

	Lock();

	openOutputFileAndControl(fOpposedFlameManager, flame->nameFolderAdditionalData + "/OpposedFlameManager.out");
	openOutputFileAndControl(fLog, flame->nameFolderAdditionalData + "/OpposedFlameManager.log");
	fOpposedFlameManager.setf(ios::scientific);
	fLog.setf(ios::scientific);

	fOpposedFlameManager << setw(16) << left << "#Index(1)" 
						 << setw(16) << left << "VF[cm/s](2)" 
						 << setw(16) << left << "VO[cm/s](3)" 
						 << setw(16) << left << "TF[K](4)" 
						 << setw(16) << left << "TO[K](5)" 
						 << setw(16) << left << "X[1/s](6)" 
						 << setw(16) << left << "TFlame[K](7)" 
						 << setw(16) << left << "X(P.F.)[1/s](8)" 
						 << endl << endl;
}

double OpenSMOKE_Flame1D_OpposedFlameManager::FuelDensity()
{
	double sum  = data->XC.GetSumElements();
	double diff = fabs(1.-sum);

	cout << "Sum of mass fractions: " << sum << endl;
	if (diff > 1.e-4)	ErrorMessage("Wrong composition on fuel side...");
	else				data->XC /= sum;

	OpenSMOKE_GasStream	fuel_inlet_stream;
	fuel_inlet_stream.AssignKineticScheme(*flame->mix);
	fuel_inlet_stream.AssignTemperature(data->TC, "K");
	fuel_inlet_stream.AssignPressure(data->P_Pascal, "Pa");
	fuel_inlet_stream.AssignMassFlowRate(1.0, "kg/s");
	fuel_inlet_stream.AssignMoleFractions(data->XC);
	fuel_inlet_stream.lock();
	return fuel_inlet_stream.rho;
}

double OpenSMOKE_Flame1D_OpposedFlameManager::OxidizerDensity()
{
	double sum  = data->XO.GetSumElements();
	double diff = fabs(1.-sum);

	cout << "Sum of mass fractions: " << sum << endl;
	if (diff > 1.e-4)	ErrorMessage("Wrong composition on oxidizer side...");
	else				data->XO /= sum;

	OpenSMOKE_GasStream	oxidizer_inlet_stream;
	oxidizer_inlet_stream.AssignKineticScheme(*flame->mix);
	oxidizer_inlet_stream.AssignTemperature(data->TO, "K");
	oxidizer_inlet_stream.AssignPressure(data->P_Pascal, "Pa");
	oxidizer_inlet_stream.AssignMassFlowRate(1.0, "kg/s");
	oxidizer_inlet_stream.AssignMoleFractions(data->XO);
	oxidizer_inlet_stream.lock();
	return oxidizer_inlet_stream.rho;
}

void OpenSMOKE_Flame1D_OpposedFlameManager::Update()
{
	double rhoFuel = FuelDensity();
	double rhoOxidizer = OxidizerDensity();

//	flame->UC =  rhoFuel*data->VC / (flame->nGeometry-1.);
//	flame->UO = -rhoOxidizer*data->VO / (flame->nGeometry-1.);

	flame->GC = -rhoFuel*data->radialGradientC;
	flame->GO = -rhoOxidizer*data->radialGradientO;

	for(int i=1;i<=flame->NC;i++)
	{
		flame->BCW_C[i] =  rhoFuel*data->VC*flame->WC[i];
		flame->BCW_O[i] = -rhoOxidizer*data->VO*flame->WO[i];
	}
}

void OpenSMOKE_Flame1D_OpposedFlameManager::Solution(const int count)
{
	// Solution
	flame->solveNLS_Opposed(1, OPPOSED_ONLY_MOMENTUM);
	flame->solveDAE_Opposed(1, OPPOSED_ALL, 1.e3);
	PrintOnFile(fLog, fOpposedFlameManager, count);

	stringstream number;
	number << count;
	std::string nameFileSolutionComplete = flame->nameFolderSteadyData + "/Solution_Number_" + number.str() + ".out";
	std::string nameFileBackUpCompleteInput	= flame->nameFolderBackupData + "/BackUp_Number_" + number.str() + ".inp";
	std::string nameFileBackUpCompleteData	= flame->nameFolderBackupData + "/BackUp_Number_" + number.str() + ".bin";
	flame->printOnFile(nameFileSolutionComplete);
	flame->printBackUpOnlyInputData(nameFileBackUpCompleteInput);
	flame->printBackUpOnlyData(nameFileBackUpCompleteData);
}

void OpenSMOKE_Flame1D_OpposedFlameManager::Run()
{
	double MinimumTemperatureDifference = 2.5;			// [K]
	double MinimumVelocityDifference	= 0.01;			// [m/s]
	double MinimumPressureDifference	= 0.05/101325;	// [Pa] = 0.05 atm

	int iCount = 0;
	if (iAssignedFuelTemperature == true)		iCount++;
	if (iAssignedOxidizerTemperature == true)	iCount++;
	if (iAssignedFuelVelocity == true)			iCount++;
	if (iAssignedOxidizerVelocity == true)		iCount++;
	if (iAssignedPressure == true)				iCount++;

//	if (iCount >  1)	ErrorMessage("Only one parameter per time can be changed...");
	if (iCount == 0)	ErrorMessage("At least one parameter per time must be changed...");

	OpenSMOKE_Flame1D_Solution flame_hot_solution;

	for (int count=1;count<=N;count++)
	{		
		if (iAssignedFuelTemperature == true)		data->TC = TFuel[count];
		if (iAssignedOxidizerTemperature == true)	data->TO = TOxidizer[count];
		if (iAssignedFuelVelocity == true)			data->VC = vFuel[count];
		if (iAssignedOxidizerVelocity == true)		data->VO = vOxidizer[count];
		if (iAssignedPressure == true)				
		{
			data->P_Pascal = Pressures[count];
			data->P_bar = Pressures[count]/1.e5;
			data->P_atm = Pressures[count]/101325.;
		}

		Update();
		Solution(count);

		if (flame->T.Max() >= Constants::TMinExtinction)
			flame_hot_solution.PasteFromExternalSolution(*flame, *data);
		else
		{
			int count_last_hot = count-1;
			
			double TC_Hot, TO_Hot, VC_Hot, VO_Hot, P_Hot;
			double TC_Cold, TO_Cold, VC_Cold, VO_Cold, P_Cold;
			
			if (iAssignedFuelTemperature == true)		TC_Hot = TFuel[count_last_hot];
			if (iAssignedOxidizerTemperature == true)	TO_Hot = TOxidizer[count_last_hot];
			if (iAssignedFuelVelocity == true)			VC_Hot = vFuel[count_last_hot];
			if (iAssignedOxidizerVelocity == true)		VO_Hot = vOxidizer[count_last_hot];
			if (iAssignedPressure == true)				P_Hot = Pressures[count_last_hot];
			
			if (iAssignedFuelTemperature == true)		TC_Cold = TFuel[count_last_hot+1];
			if (iAssignedOxidizerTemperature == true)	TO_Cold = TOxidizer[count_last_hot+1];
			if (iAssignedFuelVelocity == true)			VC_Cold = vFuel[count_last_hot+1];
			if (iAssignedOxidizerVelocity == true)		VO_Cold = vOxidizer[count_last_hot+1];
			if (iAssignedPressure == true)				P_Cold = Pressures[count_last_hot+1];
			
			flame->PasteFromExternalSolution(flame_hot_solution);

			for(;;)
			{
				count++;
				if (iAssignedFuelTemperature == true)		data->TC = 0.50*(TC_Hot+TC_Cold);
				if (iAssignedOxidizerTemperature == true)	data->TO = 0.50*(TO_Hot+TO_Cold);
				if (iAssignedFuelVelocity == true)			data->VC = 0.50*(VC_Hot+VC_Cold);
				if (iAssignedOxidizerVelocity == true)		data->VO = 0.50*(VO_Hot+VO_Cold);
				if (iAssignedPressure == true)				
				{
					data->P_Pascal = 0.50*(P_Hot+P_Cold);
					data->P_bar = data->P_Pascal/1.e5;
					data->P_atm = data->P_Pascal/101325.;
				}

				Update();
				Solution(count);

				if (flame->T.Max() >= Constants::TMinExtinction)
				{
					flame_hot_solution.PasteFromExternalSolution(*flame, *data);
					TC_Hot = data->TC;
					TO_Hot = data->TO;
					VC_Hot = data->VC;
					VO_Hot = data->VO;
				}
				else
				{
					flame->PasteFromExternalSolution(flame_hot_solution);
					TC_Cold = data->TC;
					TO_Cold = data->TO;
					VC_Cold = data->VC;
					VO_Cold = data->VO;
				}

				if (iAssignedFuelTemperature == true)		
					if (fabs(TC_Cold-TC_Hot) <= MinimumTemperatureDifference)	break;
				if (iAssignedOxidizerTemperature == true)	
					if (fabs(TO_Cold-TO_Hot) <= MinimumTemperatureDifference)	break;
				if (iAssignedFuelVelocity == true)
					if (fabs(VC_Cold-VC_Hot) <= MinimumVelocityDifference)		break;
				if (iAssignedOxidizerVelocity == true)
					if (fabs(VO_Cold-VO_Hot) <= MinimumVelocityDifference)		break;
				if (iAssignedPressure == true)	
					if (fabs(P_Cold-P_Hot) <= MinimumPressureDifference)		break;
			}

			cout << "--------------------------------------------" << endl;
			cout << "               Final Results                " << endl;
			cout << "--------------------------------------------" << endl;

			if (iAssignedFuelTemperature == true)
			{
				cout << " Minimum TC: " << TC_Cold << " K" << endl;
				cout << " Maximum TC: " << TC_Hot  << " K" << endl;
				cout << " Mean TC:    " << (TC_Cold+TC_Hot)/2.  << " K" << endl;
			}

			if (iAssignedOxidizerTemperature == true)
			{
				cout << " Minimum TO: " << TO_Cold << " K" << endl;
				cout << " Maximum TO: " << TO_Hot  << " K" << endl;
				cout << " Mean TO:    " << (TO_Cold+TO_Hot)/2.  << " K" << endl;
			}

			if (iAssignedFuelVelocity == true)
			{
				cout << " Minimum VC: " << VC_Cold*100. << " cm/s" << endl;
				cout << " Maximum VC: " << VC_Hot*100.  << " cm/s" << endl;
				cout << " Mean VC:    " << (VC_Cold+VC_Hot)/2.*100.  << " cm/s" << endl;
			}

			if (iAssignedOxidizerVelocity == true)
			{
				cout << " Minimum VO: " << VO_Cold*100. << " cm/s" << endl;
				cout << " Maximum VO: " << VO_Hot*100.  << " cm/s" << endl;
				cout << " Mean VO:    " << (VO_Cold+VO_Hot)/2.*100.  << " cm/s" << endl;
			}

			if (iAssignedPressure == true)
			{
				cout << " Minimum P: " << P_Cold << " Pa" << endl;
				cout << " Maximum P: " << P_Hot  << " Pa" << endl;
				cout << " Mean P:    " << (P_Cold+P_Hot)/2.*100.  << " Pa" << endl;
			}

			count = N+1;
		}
	}

	fLog.close();
	fOpposedFlameManager.close();
}

void OpenSMOKE_Flame1D_OpposedFlameManager::PrintOnFile(ofstream &fLog, ofstream &fFlame,const int index)
{
	double strainRate = 0;
	double strainRatePoolFire = 0;
		strainRate = 2.*data->VO/data->L*(1.+data->VC/data->VO*sqrt(flame->rho[1]/flame->rho[flame->Np]));

	if (data->iPoolFire != POOL_FIRE_NONE)
		strainRatePoolFire = 2.*data->VO/data->L;

	fLog << endl;
	fLog << "----------------------------------------------------------------" << endl;
	fLog << " Flame #" << index << endl;
	fLog << "----------------------------------------------------------------" << endl;
	fLog << " TFlame        " << flame->T.Max()			<< " K" << endl;
	fLog << " vFuel         " << data->VC*100			<< " cm/s" << endl;
	fLog << " vOxidizer     " << data->VO*100			<< " cm/s" << endl;
	fLog << " TFuel         " << data->TC				<< " K" << endl;
	fLog << " TOxidizer     " << data->TO				<< " K" << endl;
	fLog << " rhoFuel       " << flame->rho[1]			<< " kg/m3" << endl;
	fLog << " rhoOxidizer   " << flame->rho[flame->Np]	<< " kg/m3" << endl;
	fLog << " Chi           " << strainRate				<< " 1/s" << endl;
	fLog << " Chi_PoolFire  " << strainRatePoolFire		<< " 1/s" << endl;
	fLog << endl;

	fFlame << setw(16) << left << index 
		   << setw(16) << left << data->VC*100. 
		   << setw(16) << left << data->VO*100. 
		   << setw(16) << left << data->TC 
		   << setw(16) << left << data->TO 
		   << setw(16) << left << strainRate  
		   << setw(16) << left << flame->T.Max()
		   << setw(16) << left << strainRatePoolFire 
		   << endl;  
}


OpenSMOKE_Dictionary_Flame1D_OpposedFlameManager::OpenSMOKE_Dictionary_Flame1D_OpposedFlameManager()
{
    SetupBase();
	SetName("OpenSMOKE_Flame1D_OpposedFlameManager Dictionary");
	Add("#FuelVelocity",		'O', 'V', "List of fuel velocities");
	Add("#OxidizerVelocity",	'O', 'V', "List of oxidizer velocities");
	Add("#FuelTemperature",		'O', 'V', "List of fuel temperatures");
	Add("#OxidizerTemperature", 'O', 'V', "List of oxidizer temperatures");
	Add("#Pressure",			'O', 'V', "List of pressures");
    Lock();
}

void OpenSMOKE_Flame1D_OpposedFlameManager::Lock()
{
}

void OpenSMOKE_Flame1D_OpposedFlameManager::SetFuelVelocity(const vector<string> string_vector)
{
	int n = string_vector.size();
	std::string units = string_vector[n-1];

	if (N>0 && (n-1)!=N)	
		ErrorMessage("Wrong number of requested flames...");

	N=n-1;
	for(int i=0;i<N;i++)
	{
		double value = OpenSMOKE_Conversions::conversion_velocity(atof(string_vector[i].c_str()), units);
		vFuel.Append(value);
	}

	iAssignedFuelVelocity = true;
}

void OpenSMOKE_Flame1D_OpposedFlameManager::SetFuelTemperature(const vector<string> string_vector)
{
	int n = string_vector.size();
	std::string units = string_vector[n-1];

	if (N>0 && (n-1)!=N)	
		ErrorMessage("Wrong number of requested flames...");

	N=n-1;
	for(int i=0;i<N;i++)
	{
		double value = OpenSMOKE_Conversions::conversion_velocity(atof(string_vector[i].c_str()), units);
		TFuel.Append(value);
	}

	iAssignedFuelTemperature = true;
}

void OpenSMOKE_Flame1D_OpposedFlameManager::SetOxidizerVelocity(const vector<string> string_vector)
{
	int n = string_vector.size();
	std::string units = string_vector[n-1];

	if (N>0 && (n-1)!=N)	
		ErrorMessage("Wrong number of requested flames...");

	N=n-1;
	for(int i=0;i<N;i++)
	{
		double value = OpenSMOKE_Conversions::conversion_velocity(atof(string_vector[i].c_str()), units);
		vOxidizer.Append(value);
	}

	iAssignedOxidizerVelocity = true;
}

void OpenSMOKE_Flame1D_OpposedFlameManager::SetOxidizerTemperature(const vector<string> string_vector)
{
	int n = string_vector.size();
	std::string units = string_vector[n-1];

	if (N>0 && (n-1)!=N)	
		ErrorMessage("Wrong number of requested flames...");

	N=n-1;
	for(int i=0;i<N;i++)
	{
		double value = OpenSMOKE_Conversions::conversion_temperature(atof(string_vector[i].c_str()), units);
		TOxidizer.Append(value);
	}

	iAssignedOxidizerTemperature = true;
}

void OpenSMOKE_Flame1D_OpposedFlameManager::SetPressure(const vector<string> string_vector)
{
	int n = string_vector.size();
	std::string units = string_vector[n-1];

	if (N>0 && (n-1)!=N)	
		ErrorMessage("Wrong number of requested flames...");

	N=n-1;
	for(int i=0;i<N;i++)
	{
		double value = OpenSMOKE_Conversions::conversion_pressure(atof(string_vector[i].c_str()), units);
		Pressures.Append(value);
	}

	iAssignedPressure = true;
}