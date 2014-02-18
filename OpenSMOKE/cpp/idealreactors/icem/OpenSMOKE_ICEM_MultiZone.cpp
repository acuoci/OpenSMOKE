/***************************************************************************
 *   Copyright (C) 2006-2008 by Alberto Cuoci							   *
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
#include <string>
#include <sstream>
#include <iomanip>
#include "idealreactors/pfr/OpenSMOKE_PFR.h"
#include "idealreactors/icem/OpenSMOKE_ICEM_MultiZone.h"
#include "basic/OpenSMOKE_Conversions.h"
#include "addons/OpenSMOKE_PostProcessor.h"
#include "addons/OpenSMOKE_SensitivityAnalysis_Fast_Unsteady0D.h"

const double OpenSMOKE_ICEM_MultiZone::rad_to_degrees  = 360./(2.*Constants::pi);
const double OpenSMOKE_ICEM_MultiZone::degrees_to_rad  = (2.*Constants::pi)/360.;
const double OpenSMOKE_ICEM_MultiZone::min_temperature =  200.;	// [K]
const double OpenSMOKE_ICEM_MultiZone::max_temperature = 6000.;	// [K]
const double OpenSMOKE_ICEM_MultiZone::min_pressure    =  100.;	// [Pa]
const double OpenSMOKE_ICEM_MultiZone::max_pressure    =  1.e8;	// [Pa]

const double epsdtmax = 1e-8;	// [K]

OpenSMOKE_ICEM_MultiZone::OpenSMOKE_ICEM_MultiZone()
{
	class_name	= "OpenSMOKE_ICEM_MultiZone";
	out_name	= "ICEM_MultiZone.out";
	
	Setup();

	assignedEnd					= true;
	assignedTotalDisplacement   = false;
	assignedClearanceVolume		= false;
	assignedCompressionRatio	= false;
	assignedArmRatio			= false;
	assignedExhaustRatio		= false;
	assignedStartAngle			= false;
	assignedEndAngle			= false;
	assignedRotationRate		= false;
	assignedNumberOfZones		= false;
	assignedNumberOfCycles		= false;
	assignedCrevicesVolume		= false;

	timeOld					= 0.;	// [s]
	delta_step_reference	= 0.01;	// [deg]
	zone_count				= 0;	// [-]
	work					= 0.;	// [J]
	TauMixing				= 0.;	// [s]
	Cphi					= 2.;	// [-]

	// Heat transfer model
	icem_heat_transfer_model	= ICEM_HEAT_TRANSFER_NONE;
	icem_exchange_area_model	= ICEM_EXCHANGE_AREA_LATERAL;
	
	// Additional options
	iHot			  = true;
	iVerboseZones     = false;
	iVerboseGeometry  = false;
	iVerboseMasses    = false;
	iMonteCarloMixing = false;
	iRelaxation		  = false;
	iPlugFlowExhaustGases = false;
	iPlugFlowExhaustGasesInletTemperature = false;
	iResetSpecies = false;
}

void OpenSMOKE_ICEM_MultiZone::AssignEngineModel(const icem_multizone_models value)
{
	icem_multizone_model		= value;
	
	if (icem_multizone_model == ICEM_CHEMKIN_MULTIZONE)
	{
		jCrevices = 0;
		jOutmost  = 1;
	}

	if (icem_multizone_model == ICEM_KOMNINOS_SPLITTING)
	{
		jCrevices = 1;
		jOutmost  = 2;
	}
}

void OpenSMOKE_ICEM_MultiZone::AssignEnd(const string units, const double value)
{
}

void OpenSMOKE_ICEM_MultiZone::AssignRotationRate(const string units, const double value)
{
	rotation_rate = OpenSMOKE_Conversions::conversion_angular_velocity(value, units);
    assignedRotationRate	= true;
}

void OpenSMOKE_ICEM_MultiZone::AssignClearanceVolume(const string units, const double value)
{
	volume_clearance = OpenSMOKE_Conversions::conversion_volume(value, units);
    assignedClearanceVolume	= true;
}

void OpenSMOKE_ICEM_MultiZone::AssignTotalDisplacement(const string units, const double value)
{
	volume_displacement_total = OpenSMOKE_Conversions::conversion_volume(value, units);
    assignedTotalDisplacement = true;
}

void OpenSMOKE_ICEM_MultiZone::AssignCrevicesVolume(const string units, const double value)
{
	volume_crevices = OpenSMOKE_Conversions::conversion_volume(value, units);
    assignedCrevicesVolume	= true;
}

void OpenSMOKE_ICEM_MultiZone::AssignNumberOfZones(const int value)
{
	N = value;
	
	if (icem_multizone_model == ICEM_KOMNINOS_SPLITTING && N<3)
		ErrorMessage("Minimum number of zones is 3");

    assignedNumberOfZones	= true;
}

void OpenSMOKE_ICEM_MultiZone::AssignCompressionRatio(const double value)
{
	compression_ratio = value;
    assignedCompressionRatio = true;
}

void OpenSMOKE_ICEM_MultiZone::AssignArmRatio(const double value)
{
	arm_ratio = value;
    assignedArmRatio = true;
}

void OpenSMOKE_ICEM_MultiZone::AssignExhaustRatio(const double value)
{
	exhaust_ratio = value;
    assignedExhaustRatio = true;

	openOutputFileAndControl(fExhaust, outputFolderName + "/" + "Exhaust.out");
	fExhaust.setf(ios::scientific);
	openOutputFileAndControl(fRecycle, outputFolderName + "/" + "Recycle.out");
	fRecycle.setf(ios::scientific);
}

void OpenSMOKE_ICEM_MultiZone::AssignStartAngle(const string units, const double value)
{
	start_angle = OpenSMOKE_Conversions::conversion_angle(value, units);
    assignedStartAngle	= true;
}

void OpenSMOKE_ICEM_MultiZone::AssignEndAngle(const string units, const double value)
{
	end_angle = OpenSMOKE_Conversions::conversion_angle(value, units);
    assignedEndAngle	= true;
}

void OpenSMOKE_ICEM_MultiZone::AssignDiameter(const string units, const double value)
{
	diameter_cylinder  = OpenSMOKE_Conversions::conversion_length(value, units);
	area_cylinder_base = Constants::pi/4.*diameter_cylinder*diameter_cylinder;
    assignedDiameter   = true;
}

void OpenSMOKE_ICEM_MultiZone::SetPlugFlowTime(const string units, const double value)
{
	iPlugFlowExhaustGases = true;
	plugFlowExhaustGasesTime  = OpenSMOKE_Conversions::conversion_time(value, units);
}

void OpenSMOKE_ICEM_MultiZone::SetPlugFlowInletTemperature(const string units, const double value)
{
	iPlugFlowExhaustGasesInletTemperature = true;
	plugFlowExhaustGasesInletTemperature  = OpenSMOKE_Conversions::conversion_temperature(value, units);
}

void OpenSMOKE_ICEM_MultiZone::SetRelaxation(const int value)
{
	if (value>100 || value <1)
		ErrorMessage("Relaxation number must be between 1 and 100!");

	relaxation_index   = 0;
	relaxation_cycles  = value;
	relaxation_matrix  = new BzzMatrix[101];
	iRelaxation        = true;
}

void OpenSMOKE_ICEM_MultiZone::SetColdSimulation()
{
	iHot = false; 
}

void OpenSMOKE_ICEM_MultiZone::SetResetSpecies(const vector<string> names)
{
	iResetSpecies = true;
	ChangeDimensions(names.size(), &indicesResetSpecies);
	for(int i=1;i<=names.size();i++)
		indicesResetSpecies[i] = mix->recognize_species(names[i-1]);
}

void OpenSMOKE_ICEM_MultiZone::AssignNumberOfCycles(const int value)
{
	number_of_cycles		= value; 
	assignedNumberOfCycles	= true;
}

void OpenSMOKE_ICEM_MultiZone::SetConstantExchangeArea(const double value, const string units)
{
	iUserDefinedExchangeArea = CONSTANT;
	area_cylinder_total = OpenSMOKE_Conversions::conversion_area(value, units);
}

void OpenSMOKE_ICEM_MultiZone::SetExchangeAreaModel(const icem_exchange_area_models value)
{
	icem_exchange_area_model = value;
}

void OpenSMOKE_ICEM_MultiZone::SetUserDefinedExchangeArea(const string fileName)
{
    iUserDefinedExchangeArea = USERDEFINED;
    ud_A_profile.AssignFromFile(fileName, "AREA");
	ud_A_profile.SetName(name_object + " - Exchange Area Profile");
}

void OpenSMOKE_ICEM_MultiZone::UnsetUserDefinedExchangeArea()
{
    iUserDefinedExchangeArea = NONE;
}


void OpenSMOKE_ICEM_MultiZone::SetVerboseZones(const vector<string> _names)
{
	GiveMeIndicesAndNames(*mix, _names, index_species_zones, names_species_zones);
	iVerboseZones = true;

	openOutputFileAndControl(fZones, outputFolderName + "/" + "Zones.out");
	fZones.setf(ios::scientific);
}

void OpenSMOKE_ICEM_MultiZone::SetVerboseGeometry()
{
    iVerboseGeometry = true;

	openOutputFileAndControl(fGeometry, outputFolderName + "/" + "Geometry.out");
	fGeometry.setf(ios::scientific);
}

void OpenSMOKE_ICEM_MultiZone::SetVerboseMasses()
{
    iVerboseMasses = true;

	openOutputFileAndControl(fMass, outputFolderName + "/" + "Mass.out");
	fMass.setf(ios::scientific);
}

void OpenSMOKE_ICEM_MultiZone::SetMonteCarloMixing()
{
	iMonteCarloMixing = true;
}

void OpenSMOKE_ICEM_MultiZone::SetMicromixingTime(const double value, const string units)
{
	TauMixing = OpenSMOKE_Conversions::conversion_time(value, units);
}

void OpenSMOKE_ICEM_MultiZone::SetMicromixingConstant(const double value)
{
	Cphi = value;
}

void OpenSMOKE_ICEM_MultiZone::SetZoneInitialConditions(const vector<string> string_vector)
{
	int total_size = string_vector.size();
	int n = total_size / 5;

	for(int j=1;j<=n;j++)
	{
		zone_count++;

		int		index = atoi(string_vector[5*(j-1)].c_str());
		double  mass_fraction = atof(string_vector[5*(j-1)+1].c_str());
		double  temperature = atof(string_vector[5*(j-1)+2].c_str());
		string  temperature_units = string_vector[5*(j-1)+3];
		double  external_exchange = atof(string_vector[5*(j-1)+4].c_str());

		if (index != zone_count)	ErrorMessage("The zone indices must be ordered. Please check yor input for zones...");

		temperature = OpenSMOKE_Conversions::conversion_temperature(temperature, temperature_units);

		if (mass_fraction <=0.) ErrorMessage("Zone mass fraction outside the ranges. Please check yor input for zone mass fractions...");
		if (mass_fraction  >1.) ErrorMessage("Zone mass fraction outside the ranges. Please check yor input for zone mass fractions...");

		if (temperature   <= min_temperature)  ErrorMessage("Temperature outside the ranges. Please check yor input for zone temperatures...");
		if (temperature   >= max_temperature)  ErrorMessage("Temperature outside the ranges. Please check yor input for zone temperatures...");

		zone_list_temperature.Append(temperature);
		zone_list_mass_fraction.Append(mass_fraction);
		zone_list_external_exchange.Append(external_exchange);
	}

	if (zone_count != N)	ErrorMessage("The number of zones is not correct...");
	if (zone_list_mass_fraction.GetSumElements() > 1.00001 || zone_list_mass_fraction.GetSumElements() < 0.99999)
		ErrorMessage("The sum of zone mass fractions must be equal to 1..."); 
	if (zone_list_external_exchange.GetSumElements() > 1.00001 || zone_list_external_exchange.GetSumElements() < 0.99999)
		ErrorMessage("The sum of zone external exchange fractions must be equal to 1..."); 

	double tmean = Dot(zone_list_temperature, zone_list_mass_fraction);
	double sigmasquared = 0.;
	for(int j=1;j<=N;j++)
		sigmasquared += BzzPow2(zone_list_temperature[j]-tmean)*zone_list_mass_fraction[j];

	cout << endl;
	cout << "--------------------------------------------------------" << endl;
	cout << "               Temperature distribution                 " << endl;
	cout << "--------------------------------------------------------" << endl;
	cout << " Mean temperature: " << tmean				<< " K" << endl;
	cout << " Std deviation:    " << sqrt(sigmasquared)	<< " K" << endl;
	for(int j=1;j<=N;j++)
		cout << " " << j << "\t" << zone_list_temperature[j] << " " << zone_list_mass_fraction[j] << endl;
	cout << "--------------------------------------------------------" << endl;
	cout << "--------------------------------------------------------" << endl;
	cout << endl;

}

void OpenSMOKE_ICEM_MultiZone::SetGaussianInitialConditions(const vector<string> string_vector)
{
	if (string_vector.size() != 13)	ErrorMessage("Wrong number of parameters in Gauss distribution definition...");

	bool iSafe;
	double tmin;
	double tmax;

	int j=0;
	for(;;)
	{
		if (string_vector[j] == "mean")
		{
			
			gaussian_mean = OpenSMOKE_Conversions::conversion_temperature(atof(string_vector[j+1].c_str()), string_vector[j+2]);
			j+=3;
		}
		else if (string_vector[j] == "std")
		{
			gaussian_sigma =	OpenSMOKE_Conversions::conversion_temperature(atof(string_vector[j+1].c_str()), string_vector[j+2]) -
								OpenSMOKE_Conversions::conversion_temperature(0., string_vector[j+2]);
			j+=3;
		}
		else if (string_vector[j] == "min")
		{
			tmin = OpenSMOKE_Conversions::conversion_temperature(atof(string_vector[j+1].c_str()), string_vector[j+2]);
			j+=3;
		}
		else if (string_vector[j] == "max")
		{
			tmax = OpenSMOKE_Conversions::conversion_temperature(atof(string_vector[j+1].c_str()), string_vector[j+2]);
			j+=3;
		}
		else if (string_vector[j] == "safe")
		{
			iSafe = true;
			j+=1;
		}
		else if (string_vector[j] == "unsafe")
		{
			iSafe = false;
			j+=1;
		}

		if (j >= 13) break;
	}

	zone_count			= N;
	ChangeDimensions(N, &zone_list_temperature);
	ChangeDimensions(N, &zone_list_mass_fraction);
	ChangeDimensions(N, &zone_list_external_exchange);
	zone_list_external_exchange[1] = 1.;

	zone_list_temperature[1] = tmin + 0.50*(tmax-tmin)/double(N);
	for(j=2;j<=N;j++)
		zone_list_temperature[j] = zone_list_temperature[j-1] + (tmax-tmin)/double(N);

	double mu;
	double sigma;
	if (iSafe == true)
	{
		BzzVector first_guess(2);
		BzzVector minimum_vector(2);
		BzzVector solution(2);
		BzzVector residuals(2);
		first_guess[1] = gaussian_mean;
		first_guess[2] = gaussian_sigma;
		NLS_Gaussian_MultiZone_Object.assignICEM(this);
		BzzNonLinearSystemObject nls(first_guess, &NLS_Gaussian_MultiZone_Object);
		nls.SetMinimumConstraints(minimum_vector);
		nls();
		nls.GetSolution(&solution, &residuals);
		mu = solution[1];
		sigma = solution[2];
	}
	else
	{
		mu = gaussian_mean;
		sigma = gaussian_sigma;
	}

	for(int j=1;j<=N;j++)
		zone_list_mass_fraction[j] = 1./sqrt(2.*Constants::pi)/sigma*exp(-BzzPow2(zone_list_temperature[j]-mu)/2./sigma/sigma);

	double sum = zone_list_mass_fraction.GetSumElements();
	for(int j=1;j<=N;j++)
		zone_list_mass_fraction[j] /= sum;

	double tmean = Dot(zone_list_temperature, zone_list_mass_fraction);
	double sigmasquared = 0.;
	for(int j=1;j<=N;j++)
		sigmasquared += BzzPow2(zone_list_temperature[j]-mu)*zone_list_mass_fraction[j];

	cout << endl;
	cout << "--------------------------------------------------------" << endl;
	cout << "               Temperature distribution                 " << endl;
	cout << "--------------------------------------------------------" << endl;
	cout << " Mean temperature: " << tmean				<< " K" << endl;
	cout << " Std deviation:    " << sqrt(sigmasquared)	<< " K" << endl;
	for(int j=1;j<=N;j++)
		cout << " " << j << "\t" << zone_list_temperature[j] << " " << zone_list_mass_fraction[j] << endl;
	cout << "--------------------------------------------------------" << endl;
	cout << " " << j << "\t" << tmean << " " << zone_list_mass_fraction.GetSumElements() << endl;
	cout << "--------------------------------------------------------" << endl;
	cout << endl;
}

void OpenSMOKE_ICEM_MultiZone::NLSSystem_ICEM_Gaussian(BzzVector &x, BzzVector &f)
{
	int j;

	for(j=1;j<=N;j++)
		zone_list_mass_fraction[j] = 1./sqrt(2.*Constants::pi)/x[2]*exp(-BzzPow2(zone_list_temperature[j]-x[1])/2./x[2]/x[2]);
	double sum = zone_list_mass_fraction.GetSumElements();
	zone_list_mass_fraction /= sum;

	double mean = Dot(zone_list_mass_fraction, zone_list_temperature);

	double variance = 0.;
	for(j=1;j<=N;j++)
		variance += zone_list_mass_fraction[j]*BzzPow2(zone_list_temperature[j]-mean);

	f[1] = mean - gaussian_mean;
	f[2] = sqrt(variance) - gaussian_sigma;
}

void OpenSMOKE_ICEM_MultiZone::SetHeatTransferModel_WoschniParameters(const vector<string> string_vector)
{
	if (icem_heat_transfer_model != ICEM_HEAT_TRANSFER_WOSCHNI)
		ErrorMessage("The Woschni list can be used only if the WOSCHNI heat transfer model is enabled...");
	
	int n = string_vector.size();
	for(int j=0;j<n;j+=2)
	{
		string option_string = string_vector[j];
		double option_value  = atof(string_vector[j+1].c_str());

			 if (option_string == "alfa")	woschni_alfa		= option_value;	
		else if (option_string == "c11")	woschni_c1_phase1	= option_value;	
		else if (option_string == "c12")	woschni_c1_phase2	= option_value;	
		else if (option_string == "c21")	woschni_c2_phase1	= option_value;	
		else if (option_string == "c22")	woschni_c2_phase2	= option_value;	
		else ErrorMessage("Wrong parameter in Woschni list...");	
	}
}

void OpenSMOKE_ICEM_MultiZone::SetHeatTransferModel_AssanisParameters(const vector<string> string_vector)
{
	if (icem_heat_transfer_model != ICEM_HEAT_TRANSFER_ASSANIS)
		ErrorMessage("The Assanis list can be used only if the ASSANIS heat transfer model is enabled...");

	int n = string_vector.size();
	for(int j=0;j<n;j+=2)
	{
		string option_string = string_vector[j];
		double option_value  = atof(string_vector[j+1].c_str());

			 if (option_string == "alfa")	assanis_alfa		= option_value;	
		else if (option_string == "c11")	assanis_c1_phase1	= option_value;	
		else if (option_string == "c12")	assanis_c1_phase2	= option_value;	
		else if (option_string == "c21")	assanis_c2_phase1	= option_value;	
		else if (option_string == "c22")	assanis_c2_phase2	= option_value;	
		else ErrorMessage("Wrong parameter in Assanis list...");	
	}
}

void OpenSMOKE_ICEM_MultiZone::SetHeatTransferModel_HohenbergParameters(const vector<string> string_vector)
{
	if (icem_heat_transfer_model != ICEM_HEAT_TRANSFER_HOHENMBERG)
		ErrorMessage("The Hohenberg list can be used only if the HOHENBERG heat transfer model is enabled...");

	int n = string_vector.size();
	for(int j=0;j<n;j+=2)
	{
		string option_string = string_vector[j];
		double option_value  = atof(string_vector[j+1].c_str());

			 if (option_string == "alfa")	hohenberg_alfa	= option_value;	
		else if (option_string == "b")		hohenberg_b		= option_value;	
		else ErrorMessage("Wrong parameter in Hohenberg list...");	
	}
}

void OpenSMOKE_ICEM_MultiZone::SetHeatTransferModel_AnnandParameters(const vector<string> string_vector)
{
	if (icem_heat_transfer_model != ICEM_HEAT_TRANSFER_ANNAND)
		ErrorMessage("The Annand list can be used only if the ANNAND heat transfer model is enabled...");

	int n = string_vector.size();
	for(int j=0;j<n;j+=2)
	{
		string option_string = string_vector[j];
		double option_value  = atof(string_vector[j+1].c_str());

			 if (option_string == "alfa")	annand_alfa	= option_value;	
		else if (option_string == "Re")		annand_Re	= option_value;	
		else ErrorMessage("Wrong parameter in Annand list...");	
	}
}

void OpenSMOKE_ICEM_MultiZone::SetHeatTransferModel_Woschni()
{
	icem_heat_transfer_model = ICEM_HEAT_TRANSFER_WOSCHNI;

	woschni_p		  =  0.80;
	woschni_l		  = -0.20;
	woschni_t		  = -0.55;
	woschni_v		  =  0.80;

	woschni_alfa	  =  130;		// []

	woschni_c1_phase1 = 2.28;		// [-]
	woschni_c1_phase2 = 2.28;		// [-]

	woschni_c2_phase1 = 0.0;		// [m/s/K]
	woschni_c2_phase2 = 3.24e-3;	// [m/s/K]
}

void OpenSMOKE_ICEM_MultiZone::SetHeatTransferModel_Assanis()
{
	icem_heat_transfer_model = ICEM_HEAT_TRANSFER_ASSANIS;

	assanis_p		  =  0.80;
	assanis_l		  = -0.20;
	assanis_t		  = -0.73;
	assanis_v		  =  0.80;

	assanis_alfa	  =  130;

	assanis_c1_phase1 = 2.28;		// [-]
	assanis_c1_phase2 = 2.28;		// [-]

	assanis_c2_phase1 = 0.0;		// [m/s/K]
	assanis_c2_phase2 = 3.24e-3;	// [m/s/K]
}

void OpenSMOKE_ICEM_MultiZone::SetHeatTransferModel_Hohenberg()
{
	icem_heat_transfer_model = ICEM_HEAT_TRANSFER_HOHENMBERG;

	hohenberg_alfa	  =  130.;
	hohenberg_b		  =  1.40;

	hohenberg_p		  =  0.80;
	hohenberg_V		  = -0.06;
	hohenberg_t		  = -0.40;
	hohenberg_v		  =  0.80;
}

void OpenSMOKE_ICEM_MultiZone::SetHeatTransferModel_Annand()
{
	icem_heat_transfer_model = ICEM_HEAT_TRANSFER_ANNAND;
	
	annand_alfa	  =  0.55;
	annand_Re	  =  0.70;
}

void OpenSMOKE_ICEM_MultiZone::Lock()
{
    if (assignedKineticScheme == false)
        ErrorMessage("The kinetic scheme was not defined!");
    if (assignedGlobalKineticScheme == false)
        ErrorMessage("The global kinetic scheme was not defined!");
    if (assignedEnd == false)
		ErrorMessage("The reactor contact time was not defined!");
    if (assignedInletFlows == false)
        ErrorMessage("The inlet flow was not defined!");
	if (assignedSoot2EModel == false)
        ErrorMessage("The 2E Model was not defined!");
	if (assignedNumberOfZones == false)
        ErrorMessage("The number of zones was not defined!");
	if (assignedCrevicesVolume == false)
        ErrorMessage("The crevices volume was not defined!");

	Initialize();
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
//																								   //
//									GAS MIXTURE PROPERTIES										   //
//																								   //
/////////////////////////////////////////////////////////////////////////////////////////////////////

void OpenSMOKE_ICEM_MultiZone::UpdateProperties_isothermal(int memoIndex)
{
	ErrorMessage("Isothermal properties cannot be used...");
}

void OpenSMOKE_ICEM_MultiZone::UpdateProperties(const int jacobianIndex, const int indexT)
{
	ErrorMessage("Properties not available...");
}

void OpenSMOKE_ICEM_MultiZone::UpdatePropertiesMixing()
{
	int j;

	// Molecular weights
	mix->GetMWAndMoleFractionsFromMassFractions(zone_mw, zone_x, zone_omega);

	for(j=1;j<=N;j++)
	{
		BzzVector aux_omega = zone_omega.GetRow(j);
		BzzVector aux_x = zone_x.GetRow(j);;

		zone_mole[j] = zone_mass[j]/zone_mw[j];					// [kmol]
		zone_rho[j]  = zone_mass[j]/zone_volume[j];				// [kg/m3]
		zone_ctot[j] = zone_rho[j]/zone_mw[j];					// [kmol/m3]
		zone_c.SetRow(j, zone_ctot[j]*aux_x);		// [kmol/m3]

		// Specific heat
		mix->SpeciesCp(zone_temperature[j]);
		zone_Cp[j] = mix->MixCp_FromMassFractions(aux_omega);		// [J/kg/K]
		zone_Cv[j] = zone_Cp[j]-Constants::R_J_kmol/zone_mw[j];					// [J/kg/K]

		if (icem_heat_transfer_model != ICEM_HEAT_TRANSFER_NONE)
		{
			// Thermal conductivity 
			mix->SpeciesConductivityFromFitting(zone_temperature[j]);
			zone_lambda[j] = mix->MixConductivity_FromMolarFractions(aux_x);		// [W/m/K]

			// Dynamic viscosity
			mix->SpeciesViscosityFromFitting(zone_temperature[j]);
			zone_mu[j] = mix->MixViscosity_FromMolarFractions(aux_x);			// [Pa.s]
		}

		// Enthalpy
		mix->GetMixAveragedEnthalpy_Mass(dummy_vector_nc, zone_temperature[j]);		// [J/kg
		zone_hmass.SetRow(j, dummy_vector_nc);										// [J/kg]
		zone_enthalpy[j] = zone_mass[j]*Dot(dummy_vector_nc,aux_omega);	// [J]
		
		// Internal energy
		mix->GetMixAveragedInternalEnergy_Mass(dummy_vector_nc, zone_temperature[j]);		// [J/kg]
		zone_umass.SetRow(j, dummy_vector_nc);												// [J/kg]
		zone_internal_energy[j] = zone_mass[j]*Dot(dummy_vector_nc,aux_omega);	// [J]
	}
}

void OpenSMOKE_ICEM_MultiZone::UpdatePropertiesSingleReactor()
{
	BzzVector aux_omega = zone_omega.GetRow(jGlobal);
	BzzVector aux_x = zone_x.GetRow(jGlobal);;

	// Molecular weights and mole fractions
	mix->GetMWAndMoleFractionsFromMassFractions(zone_mw[jGlobal], dummy_vector_nc, aux_omega);
	zone_x.SetRow(jGlobal, dummy_vector_nc);

	// Density and concentrations
	zone_mole[jGlobal] = zone_mass[jGlobal]/zone_mw[jGlobal];				// [kmol]
	zone_rho[jGlobal]  = zone_mass[jGlobal]/zone_volume[jGlobal];			// [kg/m3]
	zone_ctot[jGlobal] = zone_rho[jGlobal]/zone_mw[jGlobal];				// [kmol/m3]
	zone_c.SetRow(jGlobal, zone_ctot[jGlobal]*aux_x);		// [kmol/m3]

	// Specific heats
	mix->SpeciesCp(zone_temperature[jGlobal]);
	zone_Cp[jGlobal] = mix->MixCp_FromMassFractions(aux_omega);		// [J/kg/K]
	zone_Cv[jGlobal] = zone_Cp[jGlobal]-Constants::R_J_kmol/zone_mw[jGlobal];			// [J/kg/K]

	// Kinetics
	if (iHot == true)
	{
		BzzVector aux_c = zone_c.GetRow(jGlobal);
		BzzVector aux_R = zone_R.GetRow(jGlobal);
		mix->ComputeKineticParameters(zone_temperature[jGlobal], log(zone_temperature[jGlobal]), 1./zone_temperature[jGlobal], zone_pressure[jGlobal]);	
		mix->ComputeFromConcentrations(zone_temperature[jGlobal], aux_c, zone_ctot[jGlobal], &dummy_vector_nc);			// [kmol/m3/s]
		ElementByElementProduct(dummy_vector_nc, mix->M, &dummy_vector_nc);					// [kg/m3/s]
		zone_R.SetRow(jGlobal, dummy_vector_nc);											// [kg/m3/s]
		zone_QReaction[jGlobal] = - mix->ComputeQReaction(zone_temperature[jGlobal]);		// [J/m3/s]

		mix->GetStandardInternalEnergy_Mass(dummy_vector_nc, zone_temperature[jGlobal]);	// [J/kg]
		zone_UReaction[jGlobal] = -Dot(dummy_vector_nc, aux_R);			// [J/m3/s]	
	}
}

void OpenSMOKE_ICEM_MultiZone::UpdateProperties(const int jacobianIndex, BzzVectorInt &jacobianVariables, int dimBlock, const int indexT)
{
	int j;
	int memoIndex;
	BzzVectorInt pointToUpdate;

	     if (jacobianIndex < 0)		memoIndex = -2;		// Calculate all without storing
	else if (jacobianIndex == 0)	memoIndex = -1;		// Calculate and storing
	else if (jacobianIndex > 0)							// 
	{
		memoIndex = 0;
		for(j=1;j<=jacobianVariables.Size();j++)
			if ( (jacobianVariables[j]%(dimBlock)) == indexT)
				pointToUpdate.Append(int(jacobianVariables[j]/dimBlock)+1);
	}

	// Molecular weights
	mix->GetMWAndMoleFractionsFromMassFractions(zone_mw, zone_x, zone_omega);

	// Moles and concentrations
	for(j=1;j<=N;j++)
	{
		zone_mole[j] = zone_mass[j]/zone_mw[j];
		zone_rho[j]  = zone_mass[j]/zone_volume[j];
		zone_ctot[j] = zone_rho[j]/zone_mw[j];
		zone_c.SetRow(j, zone_ctot[j]*zone_x.GetRow(j));
	}

	// --------------------------------------------------------------------------
	// Update every variable which does not depend on T
	// --------------------------------------------------------------------------
	if(memoIndex == -2)
	{
		for(j=1;j<=N;j++)
		{
			BzzVector aux_x = zone_x.GetRow(j);
			BzzVector aux_omega = zone_omega.GetRow(j);
			BzzVector aux_c = zone_c.GetRow(j);

			// Specific heat
			mix->SpeciesCp(zone_temperature[j]);
			zone_Cp[j] = mix->MixCp_FromMassFractions(aux_omega);		// [J/kg/K]
			zone_Cv[j] = zone_Cp[j]-Constants::R_J_kmol/zone_mw[j];					// [J/kg/K]

			if (icem_heat_transfer_model == ICEM_HEAT_TRANSFER_ANNAND)
			{
				// Thermal conductivity 
				mix->SpeciesConductivityFromFitting(zone_temperature[j]);
				zone_lambda[j] = mix->MixConductivity_FromMolarFractions(aux_x);		// [W/mK]

				// Dynamic viscosity
				mix->SpeciesViscosityFromFitting(zone_temperature[j]);
				zone_mu[j] = mix->MixViscosity_FromMolarFractions(aux_x);			// [Pa.s]
			}

			// Kinetics
			if (iHot == true)
			{
				mix->ComputeKineticParameters(zone_temperature[j], log(zone_temperature[j]), 1./zone_temperature[j], zone_pressure[j]);	
				mix->ComputeFromConcentrations(zone_temperature[j], aux_c, zone_ctot[j], &dummy_vector_nc);			// [kmol/m3/s]
				ElementByElementProduct(dummy_vector_nc, mix->M, &dummy_vector_nc);			// [kg/m3/s]
				zone_R.SetRow(j, dummy_vector_nc);											// [kg/m3/s]
				zone_QReaction[j] = - mix->ComputeQReaction(zone_temperature[j]);			// [J/m3/s]
				
				mix->GetStandardInternalEnergy_Mass(dummy_vector_nc, zone_temperature[j]);	// [J/kg]
				zone_UReaction[j] = -Dot(dummy_vector_nc, zone_R.GetRow(j));				// [J/m3/s]	
			}
		}
	}

	// --------------------------------------------------------------------------
	// Update only the variables depending on the temperature
	// --------------------------------------------------------------------------
	else if(memoIndex==-1)
	{	
		for(j=1;j<=N;j++)
		{
			// Specific Heat [J/kgK]
			// ----------------------------------------------------------------------
			mix->SpeciesCp(zone_temperature[j]);
			zone_CpMap.SetRow(j, mix->Cp);

			if (icem_heat_transfer_model == ICEM_HEAT_TRANSFER_ANNAND)
			{
				// Thermal conductivity [W/mK]
				mix->SpeciesConductivityFromFitting(zone_temperature[j]);
				zone_lambdaMap.SetRow(j, mix->lambda);

				// Dynamic viscosity [Pa.s]
				mix->SpeciesViscosityFromFitting(zone_temperature[j]);
				zone_muMap.SetRow(j, mix->eta);
			}

			if (iHot == true)
			{
				mix->ComputeKineticParameters(zone_temperature[j], log(zone_temperature[j]), 1./zone_temperature[j], zone_pressure[j]);	
				zone_k1Map.SetRow(j, mix->k1);
				zone_k2Map.SetRow(j, mix->k2);
				zone_uKeqMap.SetRow(j, mix->uKeq);
				zone_logFcentMap.SetRow(j, mix->logFcent);
				zone_reactionDHMap.SetRow(j, mix->reactionDH);
				zone_reactionDSMap.SetRow(j, mix->reactionDS);
			}
		}

		memoIndex = 0;
	}

	// --------------------------------------------------------------------------
	// Update every variable which does not depend on T
	// --------------------------------------------------------------------------
	if(memoIndex == 0)
	{	
		for(j=1;j<=N;j++)
		{
			BzzVector aux_omega = zone_omega.GetRow(j);
			BzzVector aux_x = zone_x.GetRow(j);;
			BzzVector aux_c = zone_c.GetRow(j);

			// Specific heat
			zone_CpMap.GetRow(j, &mix->Cp);
			zone_Cp[j] = mix->MixCp_FromMassFractions(aux_omega);					// [J/kg/K]
			zone_Cv[j] = zone_Cp[j]-Constants::R_J_kmol/zone_mw[j];					// [J/kg/K]

			if (icem_heat_transfer_model == ICEM_HEAT_TRANSFER_ANNAND)
			{
				// Thermal conductivity [W/m/K]
				zone_lambdaMap.GetRow(j, &mix->lambda);
				zone_lambda[j] = mix->MixConductivity_FromMolarFractions(aux_x);

				// Dynamic viscosity [Pa.s]
				zone_muMap.GetRow(j, &mix->eta);
				zone_mu[j] = mix->MixViscosity_FromMolarFractions(aux_x);
			}

			if (iHot == true)
			{
				zone_k1Map.GetRow(j, &mix->k1);
				zone_k2Map.GetRow(j, &mix->k2);
				zone_uKeqMap.GetRow(j, &mix->uKeq);
				zone_logFcentMap.GetRow(j, &mix->logFcent);
				zone_reactionDHMap.GetRow(j, &mix->reactionDH);
				zone_reactionDSMap.GetRow(j, &mix->reactionDS);
				mix->ComputeFromConcentrations( zone_temperature[j], aux_c, zone_ctot[j], &dummy_vector_nc);		// [kmol/m3/s]
				ElementByElementProduct(dummy_vector_nc, mix->M, &dummy_vector_nc);					// [kg/m3/s]
				zone_R.SetRow(j, dummy_vector_nc);
				zone_QReaction[j] = - mix->ComputeQReaction(zone_temperature[j]);					// [J/m3.s]

				mix->GetStandardInternalEnergy_Mass(dummy_vector_nc, zone_temperature[j]);	// [J/kg]
				zone_UReaction[j] = -Dot(dummy_vector_nc, zone_R.GetRow(j));				// [J/m3/s]		
			}
		}
	}

	// Points to update
	for(int i=1;i<=pointToUpdate.Size();i++)
	{	
		j = pointToUpdate[i];

		BzzVector aux_x = zone_x.GetRow(j);
		BzzVector aux_omega = zone_omega.GetRow(j);
		BzzVector aux_c = zone_c.GetRow(j);

		// Specific heat
		mix->SpeciesCp(zone_temperature[j]);
		zone_Cp[j] = mix->MixCp_FromMassFractions(aux_omega);		// [J/kg/K]
		zone_Cv[j] = zone_Cp[j]-Constants::R_J_kmol/zone_mw[j];					// [J/kg/K]

		if (icem_heat_transfer_model == ICEM_HEAT_TRANSFER_ANNAND)
		{
			// Thermal conductivity 
			mix->SpeciesConductivityFromFitting(zone_temperature[j]);
			zone_lambda[j] = mix->MixConductivity_FromMolarFractions(aux_x);					BzzVector aux_c = zone_c.GetRow(j);
			BzzVector aux_R = zone_R.GetRow(j);// [W/mK]

			// Dynamic viscosity
			mix->SpeciesViscosityFromFitting(zone_temperature[j]);
			zone_mu[j] = mix->MixViscosity_FromMolarFractions(aux_x);			// [Pa.s]
		}

		// Kinetics
		if (iHot == true)
		{
			mix->ComputeKineticParameters(zone_temperature[j], log(zone_temperature[j]), 1./zone_temperature[j], zone_pressure[j]);	
			mix->ComputeFromConcentrations(zone_temperature[j], aux_c, zone_ctot[j], &dummy_vector_nc);			// [kmol/m3/s]
			ElementByElementProduct(dummy_vector_nc, mix->M, &dummy_vector_nc);			// [kg/m3/s]
			zone_R.SetRow(j, dummy_vector_nc);											// [kg/m3/s]
			zone_QReaction[j] = - mix->ComputeQReaction(zone_temperature[j]);			// [J/m3/s]

			mix->GetStandardInternalEnergy_Mass(dummy_vector_nc, zone_temperature[j]);	// [J/kg]
			zone_UReaction[j] = -Dot(dummy_vector_nc, zone_R.GetRow(j));				// [J/m3/s]	
		}
	}
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
//																								   //
//									UPDATING FUNCTIONS											   //
//																								   //
/////////////////////////////////////////////////////////////////////////////////////////////////////

void OpenSMOKE_ICEM_MultiZone::UpdateTwoEquationModel(BzzVector &y, BzzVector &dy)
{
	ErrorMessage("The UpdateTwoEquationModel cannot be used");
}

void OpenSMOKE_ICEM_MultiZone::UpdateHeatFlux(const double tau, const double csi)
{
	ErrorMessage("The update heat flux cannot be used");
}

void OpenSMOKE_ICEM_MultiZone::UpdateExchangeArea(const double tau, const double csi)
{
	if (iUserDefinedExchangeArea == NONE)			area_cylinder_total = 2.*area_cylinder_base + area_cylinder_lateral;
	if (iUserDefinedExchangeArea == USERDEFINED)	area_cylinder_total = ud_A_profile.GiveMeValue(tau, MINUSONE);
}


void OpenSMOKE_ICEM_MultiZone::UpdateAngleAndVolume(const double teta)
{
	if (icem_multizone_model == ICEM_CHEMKIN_MULTIZONE)
	{
		// Volumes are already available
		double sqrt_coeff		= sqrt(arm_ratio*arm_ratio-BzzPow2(sin(teta)));
		volume_cylinder_total	= volume_clearance * (1.+0.50*(compression_ratio-1.)*(arm_ratio+1.-cos(teta)-sqrt_coeff));
		volume_total			= volume_cylinder_total + volume_crevices;
		height_cylinder			= volume_cylinder_total/area_cylinder_base; 
		area_cylinder_lateral	= 4.*volume_cylinder_total/diameter_cylinder;
		area_cylinder_total		= 2.*area_cylinder_base+area_cylinder_lateral;

		height_e = height_cylinder;
		height_i = height_cylinder;

		diameter_e[1]		=	diameter_cylinder; 
		for(int j=1;j<=N-1;j++)
		{
			diameter_i[j]		=	sqrt(diameter_e[j]*diameter_e[j] - 4.*zone_volume[j]/Constants::pi/height_cylinder); 
			diameter_e[j+1]     =   diameter_i[j];
		}
		diameter_i[N] = 0.;

		for(int j=1;j<=N;j++)
			area_crown[j] = Constants::pi/4.*(diameter_e[j]*diameter_e[j] - diameter_i[j]*diameter_i[j]);
	}

	else if (icem_multizone_model == ICEM_KOMNINOS_SPLITTING)
	{
		double sqrt_coeff		= sqrt(arm_ratio*arm_ratio-BzzPow2(sin(teta)));
		volume_cylinder_total	= volume_clearance * (1.+0.50*(compression_ratio-1.)*(arm_ratio+1.-cos(teta)-sqrt_coeff));
		height_cylinder			= volume_cylinder_total/area_cylinder_base; 
		area_cylinder_lateral	= 4.*volume_cylinder_total/diameter_cylinder;
		area_cylinder_total		= 2.*area_cylinder_base+area_cylinder_lateral;
		volume_total			= volume_cylinder_total + volume_crevices;

		// Update Exchange Areas and volumes
		zone_volume[1] = volume_crevices;
		for(int j=2;j<=N-1;j++)
		{
			height_e[j]		= height_cylinder - 2.*(j-2)*thickness;
			height_i[j]		= height_cylinder - 2.*(j-1)*thickness;
			area_e[j]		= Constants::pi/4.*BzzPow2(diameter_e[j])*2. + Constants::pi*diameter_e[j]*height_e[j];
			area_i[j]		= Constants::pi/4.*BzzPow2(diameter_i[j])*2. + Constants::pi*diameter_i[j]*height_i[j];
			zone_volume[j]	= Constants::pi/4.*( BzzPow2(diameter_e[j])*height_e[j]-BzzPow2(diameter_i[j])*height_i[j] );
		}
		height_i[N]		= 0.;
		area_i[N]		= 0.;
		height_e[N]		= height_cylinder - 2.*(N-2)*thickness;
		area_e[N]		= Constants::pi/4.*BzzPow2(diameter_e[N])*2. + Constants::pi*diameter_e[N]*height_e[N];
		zone_volume[N]	= Constants::pi/4.*BzzPow2(diameter_e[N])*height_e[N];

		// Check
		double sum_volume = zone_volume.GetSumElements();
		if (fabs(volume_total-sum_volume)>1.e-8)	
		{
			cout << "volume_total: " << volume_total << endl;
			cout << "sum_volume:   " << sum_volume << endl;
			cout << "difference:   " << volume_total - sum_volume << endl;
			ErrorMessage("Volume does not fit!");
		}
	}
}


void OpenSMOKE_ICEM_MultiZone::UpdateVolumeDerivatives(const double teta)
{
	if (icem_multizone_model == ICEM_CHEMKIN_MULTIZONE)
	{
		double t = (teta+start_angle)/rotation_rate;
		double deltat = t-timeOld;
		if (deltat < 1.e-16)	dvolume_over_dt = 0.;
		else
			for(int j=1;j<=N;j++)
				dvolume_over_dt[j] = (zone_volume[j]-zone_volume_old[j])/deltat;	
	}

	else if (icem_multizone_model == ICEM_KOMNINOS_SPLITTING)
	{
		double sqrt_coeff = sqrt(arm_ratio*arm_ratio-BzzPow2(sin(teta)));
		double fc = rotation_rate*0.50*(compression_ratio-1.)*sin(teta)*(1.+cos(teta)/sqrt_coeff);
	
		dvolume_over_dt[1] = 0.;
		for(int j=2;j<=N;j++)
			dvolume_over_dt[j]	= area_crown[j] * height_cylinder_minimum * fc;
	}
}

void OpenSMOKE_ICEM_MultiZone::UpdateHeatTransfer()
{
	if (icem_heat_transfer_model == ICEM_HEAT_TRANSFER_NONE)
	{	
		heat_transfer_velocity = 0.;
		U = 0.;
		return;
	}

	else if (icem_heat_transfer_model == ICEM_HEAT_TRANSFER_WOSCHNI)
	{
		double gamma = zone_Cp[jOutmost]/zone_Cv[jOutmost];										// [-]
		double volume_swept = volume_clearance*(compression_ratio-1.);							// [m3]
		double pressure_motored = pressure_reference*pow(volume_reference/volume_total, gamma);	// [Pa]
		heat_transfer_velocity = woschni_c1_phase1*velocity_piston + woschni_c2_phase1*volume_swept*temperature_reference/pressure_reference/volume_reference*(zone_pressure[jOutmost]-pressure_motored);
		
		U = woschni_alfa * pow(diameter_cylinder, woschni_l) * pow(zone_pressure[jOutmost]/1.e5, woschni_p) *
						   pow(zone_temperature[jOutmost], woschni_t) * pow(heat_transfer_velocity, woschni_v);	// [W/m2/K]
	}

	else if (icem_heat_transfer_model == ICEM_HEAT_TRANSFER_ASSANIS)
	{
		double gamma = zone_Cp[jOutmost]/zone_Cv[jOutmost];													// [-]
		double volume_swept = volume_clearance*(compression_ratio-1.);							// [m3]
		double pressure_motored = pressure_reference*pow(volume_reference/volume_total, gamma);	// [Pa]
		heat_transfer_velocity = assanis_c1_phase1*velocity_piston + assanis_c2_phase1/6.*volume_swept*temperature_reference/pressure_reference/volume_reference*(zone_pressure[jOutmost]-pressure_motored);
		
		U = assanis_alfa * pow(height_cylinder, assanis_l) * pow(zone_pressure[jOutmost]/1.e5, assanis_p) *
						   pow(zone_temperature[jOutmost], assanis_t) * pow(heat_transfer_velocity, assanis_v);	// [W/m2/K]
	}

	else if (icem_heat_transfer_model == ICEM_HEAT_TRANSFER_HOHENMBERG)
	{
		heat_transfer_velocity = velocity_piston + hohenberg_b;
		
		U = hohenberg_alfa * pow(volume_total, hohenberg_V) * pow(zone_pressure[jOutmost]/1.e5, hohenberg_p) *
						     pow(zone_temperature[jOutmost], hohenberg_t) * pow(heat_transfer_velocity, hohenberg_v);				// [W/m2/K]
	}

	else if (icem_heat_transfer_model == ICEM_HEAT_TRANSFER_ANNAND)
	{
		heat_transfer_velocity = velocity_piston;
		Re = zone_rho[jOutmost]*velocity_piston*diameter_cylinder / zone_mu[jOutmost];		// [-]
		U = annand_alfa * zone_lambda[jOutmost] * pow(Re, annand_Re)/diameter_cylinder;		// [W/m2/K]
	}

	if (icem_multizone_model == ICEM_CHEMKIN_MULTIZONE)
	{
		for(int j=1;j<=N;j++)
			zone_Qe[j] = U*(zone_temperature[j]-Tambient);

		double area_effective;
		if (icem_exchange_area_model == ICEM_EXCHANGE_AREA_LATERAL)
			area_effective = area_cylinder_lateral;
		else if (icem_exchange_area_model == ICEM_EXCHANGE_AREA_TOTAL)
			area_effective = area_cylinder_lateral + 2.*area_cylinder_base;

		if (zone_count > 0)
		{
			for(int j=1;j<=N;j++)
				zone_area_effective[j] = zone_list_external_exchange[j]*area_effective;
		}
		else
			zone_area_effective[jOutmost] = area_effective;
	}

	else if (icem_multizone_model == ICEM_KOMNINOS_SPLITTING)
	{
		zone_Qe[jCrevices] = U*(zone_temperature[jCrevices]-Tambient);
		zone_Qe[jOutmost]  = U*(zone_temperature[jOutmost]-Tambient);
	}
}

void OpenSMOKE_ICEM_MultiZone::InitializeGeometry()
{
	if (assignedClearanceVolume == true)
		volume_displacement_total = volume_clearance*compression_ratio;
	else if (assignedTotalDisplacement == true)
		volume_clearance = volume_displacement_total/compression_ratio;

	volume_displacement_swept = volume_clearance*(compression_ratio-1.);
	area_cylinder_base = Constants::pi/4.*BzzPow2(diameter_cylinder);			// [m2]
	height_cylinder_minimum = volume_clearance/area_cylinder_base;				// [m]
	arm_La = (compression_ratio-1.)/2.*volume_clearance/area_cylinder_base;		// [m]
	arm_Lc = arm_ratio*arm_La;													// [m]
	velocity_piston = 2.*arm_La*rotation_rate;									// [m/s]


	if (icem_multizone_model == ICEM_CHEMKIN_MULTIZONE)
	{
		thickness			= 0.;					// [m]
		height_core_minimum = 0.;					// [m]
		area_crevices		= 0.;					// [m2]
		
		diameter_e[1]	= diameter_cylinder;		// [m]
		diameter_i[N]	= 0.;						// [m]
	}

	else if (icem_multizone_model == ICEM_KOMNINOS_SPLITTING)
	{
		int j;

		thickness = height_cylinder_minimum/2./double(N-1);
		height_core_minimum = height_cylinder_minimum - 2*(N-2)*thickness;
		area_crevices = volume_crevices*(2.*diameter_cylinder-thickness)/thickness/(diameter_cylinder-thickness) + 
						Constants::pi*thickness*(diameter_cylinder-thickness);
		
		// External diameter
							diameter_e[1]	= 0.;
		for(j=2;j<=N-1;j++)	diameter_e[j]	= diameter_cylinder - 2.*(j-2)*thickness;
							diameter_e[N]	= diameter_cylinder - 2.*(N-2)*thickness;
		
		// Internal diameter	
							diameter_i[1]	= 0.;
		for(j=2;j<=N-1;j++)	diameter_i[j]	= diameter_cylinder - 2.*(j-1)*thickness;
							diameter_i[N]	= 0.;

							area_crown[1]	= 0.;
		for(j=2;j<=N-1;j++)	area_crown[j]	= Constants::pi/4.*(BzzPow2(diameter_e[j])-BzzPow2(diameter_i[j]));
							area_crown[N]	= Constants::pi/4.*BzzPow2(diameter_e[N]);
	}
}

void OpenSMOKE_ICEM_MultiZone::MemoryAllocation()
{
	// Geometry
	ChangeDimensions(N, &zone_volume);
	ChangeDimensions(N, &area_i);
	ChangeDimensions(N, &area_e);
	ChangeDimensions(N, &diameter_i);
	ChangeDimensions(N, &diameter_e);
	ChangeDimensions(N, &height_i);
	ChangeDimensions(N, &height_e);
	ChangeDimensions(N, &area_crown);


	// Zone variables
	ChangeDimensions(N, &zone_mass);
	ChangeDimensions(N, &zone_mole);
	ChangeDimensions(N, &zone_temperature);
	ChangeDimensions(N, &zone_pressure);
	ChangeDimensions(N,	&zone_g);
	ChangeDimensions(N, &zone_enthalpy);
	ChangeDimensions(N, &zone_internal_energy);
	ChangeDimensions(N,	&zone_heat_release);
	ChangeDimensions(N, NC, &zone_omega);
	ChangeDimensions(N, NC, &zone_x);
	ChangeDimensions(N, NC, &zone_c);
	ChangeDimensions(N, NC, &zone_hmass);
	ChangeDimensions(N, NC, &zone_umass);


	// Equations
	ChangeDimensions(N, NC, &domega_over_dt);
	ChangeDimensions(N, &dtemperature_over_dt);
	ChangeDimensions(N, &dmass_over_dt);
	ChangeDimensions(N,	&g_equation);
	ChangeDimensions(N, &p_equation);


	// Properties
	ChangeDimensions(N,  &zone_mw);
	ChangeDimensions(N,  &zone_ctot);
	ChangeDimensions(N,  &zone_rho);
	ChangeDimensions(N,  &zone_Cp);
	ChangeDimensions(N,  &zone_Cv);
	ChangeDimensions(N,  &zone_mu);
	ChangeDimensions(N,  &zone_lambda);
	ChangeDimensions(N,  &zone_QReaction);
	ChangeDimensions(N,  &zone_UReaction);
	ChangeDimensions(N,  NC, &zone_R);


	// Mean properties
	ChangeDimensions(NC, &mean_omega);
	ChangeDimensions(NC, &mean_x);


	// Heat exchange
	ChangeDimensions(N,  &zone_area_effective);
	ChangeDimensions(N,  &zone_Qe);


	// Additional vectors
	ChangeDimensions(N, &dvolume_over_dt);
	ChangeDimensions(N, &zone_volume_old);
	ChangeDimensions(NC, &dummy_vector_nc);
	
	
	if (icem_multizone_model == ICEM_CHEMKIN_MULTIZONE)
	{
		// Maps
		ChangeDimensions(N,	NC, &zone_CpMap);
		ChangeDimensions(N,	NC, &zone_lambdaMap);
		ChangeDimensions(N,	NC, &zone_muMap);
		ChangeDimensions(N, NR, &zone_k1Map);
		ChangeDimensions(N, NR, &zone_k2Map);
		ChangeDimensions(N, NR, &zone_uKeqMap);
		ChangeDimensions(N, NR, &zone_logFcentMap);
		ChangeDimensions(N, NR, &zone_reactionDSMap);
		ChangeDimensions(N, NR, &zone_reactionDHMap);
	}

	if (icem_multizone_model == ICEM_KOMNINOS_SPLITTING)
	{
		// Zone fluxes
		ChangeDimensions(N+1,		&zone_m_in);
		ChangeDimensions(N+1,		&zone_h_in);
		ChangeDimensions(N+1, NC,	&zone_hmass_in);
		ChangeDimensions(N+1, NC,	&zone_omega_in);

		// Thermal conduction
		ChangeDimensions(N+1,  &zone_Qdiffusion);
		ChangeDimensions(N+1,  &zone_lambda_in);
		ChangeDimensions(N+1,  &zone_temperature_gradient_in);
		ChangeDimensions(N+1,  &viscosity_ratio_in);
		ChangeDimensions(N+1,  &yPlus_in);
	}
}

void OpenSMOKE_ICEM_MultiZone::GeometrySummary()
{
	cout << endl;
	cout << "---------------------------------------------------------------------" << endl;
	cout << "           Internal Combustion Engine Model: Geometry Summary        " << endl;
    cout << "---------------------------------------------------------------------" << endl;
	cout << " Total Displacement [cm3]:  " << volume_displacement_total*1.e6 << endl;
	cout << " Swept Displacement [cm3]:  " << volume_displacement_swept*1.e6 << endl;
	cout << " Volume Clearance [cm3]:    " << volume_clearance*1.e6 << endl;
	cout << " Minimum Height [cm]:       " << height_cylinder_minimum*1.e2 << endl;
	cout << " Maximum Height [cm]:       " << height_cylinder*1.e2 << endl;
	cout << " Area Cyl base[cm2]:        " << area_cylinder_base*1.e4 << endl;
	cout << " Area Cyl lat [cm2]:        " << area_cylinder_lateral*1.e4 << endl;
	cout << " Area Cyl tot [cm2]:        " << area_cylinder_total*1.e4 << endl;
	cout << " La [cm]:                   " << arm_La*1e2 << endl;
	cout << " Lc [cm]:                   " << arm_Lc*1e2 << endl;
	cout << " v [m/s]:                   " << velocity_piston << endl;
	cout << "---------------------------------------------------------------------" << endl;
    cout << "---------------------------------------------------------------------" << endl;
	cout << endl;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
//																								   //
//									SOLVING FUNCTIONS											   //
//																								   //
/////////////////////////////////////////////////////////////////////////////////////////////////////

void OpenSMOKE_ICEM_MultiZone::Initialize()
{
	MemoryAllocation();
	
	InitializeGeometry();

	teta = rotation_rate*0 - start_angle;
	UpdateAngleAndVolume(teta);
	UpdateVolumeDerivatives(teta);
	volume_total_starting = volume_total;

	GeometrySummary();
	
    T		= 0.;
    omega	= 0.;
	P		= 0.;

	zone_temperature	= inletStream->T;
	zone_pressure		= inletStream->P;
	for(int j=1;j<=N;j++)	
		zone_omega.SetRow(j,inletStream->omega);

	if (icem_multizone_model == ICEM_CHEMKIN_MULTIZONE)
	{
		if (zone_count == 0)	
		{
			mass_total_initial = inletStream->P*inletStream->MW*volume_total/Constants::R_J_kmol/inletStream->T;
			mass_total = mass_total_initial;

			zone_volume = volume_total / double(N);

			int j;
			for(j=1;j<=N;j++)	zone_mass[j] = inletStream->rho*zone_volume[j];
		}
		else
		{
			mass_total_initial = inletStream->P*inletStream->MW*volume_total/Constants::R_J_kmol/Dot(zone_list_temperature, zone_list_mass_fraction);
			mass_total = mass_total_initial;

			int j;
			for(j=1;j<=N;j++)	zone_volume[j]		= (mass_total_initial*zone_list_mass_fraction[j])*Constants::R_J_kmol*zone_temperature[j]/inletStream->P/inletStream->MW;
			for(j=1;j<=N;j++)	zone_temperature[j] = zone_list_temperature[j];
			for(j=1;j<=N;j++)	zone_mass[j]		= zone_list_mass_fraction[j]*mass_total_initial;
		}
	}

	else if (icem_multizone_model == ICEM_KOMNINOS_SPLITTING)
	{
		int j;
		for(j=1;j<=N;j++)	zone_mass[j] = inletStream->rho*zone_volume[j];
		mass_total_initial	= zone_mass.GetSumElements();
	}
	
	// Old variables
	zone_mass_old		= zone_mass;
	zone_volume_old		= zone_volume;

	// Reference variables for heat transfer
	pressure_reference		= zone_pressure[jOutmost];		// [Pa]
	volume_reference		= volume_total;					// [m3]
	temperature_reference	= zone_temperature[jOutmost];	// [K]

	BzzVectorInt dummy_int;
	UpdateProperties(MINUSONE, dummy_int, 0, NC+1);


	if (icem_multizone_model == ICEM_CHEMKIN_MULTIZONE)
	{
		zone_g[1] = zone_mass[1]/zone_mw[1]*Constants::R_J_kmol*zone_temperature[1];
		for(int j=2;j<=N;j++)	
			zone_g[j] = zone_g[j-1] + zone_mass[j]/zone_mw[j]*Constants::R_J_kmol*zone_temperature[j];
	}

	if (icem_multizone_model == ICEM_CHEMKIN_MULTIZONE)
	{
		if (zone_count != 0)
		{
			int j;

			for(j=1;j<=N;j++)
			{
				BzzVector aux_omega = zone_omega.GetRow(j);

				// Enthalpy
				mix->GetMixAveragedEnthalpy_Mass(dummy_vector_nc, zone_temperature[j]);		// [J/kg
				zone_hmass.SetRow(j, dummy_vector_nc);										// [J/kg]
				zone_enthalpy[j] = zone_mass[j]*Dot(dummy_vector_nc,aux_omega);	// [J]
		
				// Internal energy
				mix->GetMixAveragedInternalEnergy_Mass(dummy_vector_nc, zone_temperature[j]);		// [J/kg]
				zone_umass.SetRow(j, dummy_vector_nc);												// [J/kg]
				zone_internal_energy[j] = zone_mass[j]*Dot(dummy_vector_nc,aux_omega);	// [J]
			}

			// Mean temperature
			double mean_temperature = Dot(zone_list_temperature, zone_list_mass_fraction);

			// Total enthalpy
			enthalpy_total = zone_enthalpy.GetSumElements();								// [J]

			// Mean composition
			mean_omega = 0.;
			for(j=1;j<=N;j++)
			{
				double mass_fraction = zone_mass[j]/mass_total;
				for(int i=1;i<=NC;i++)
					mean_omega[i] += mass_fraction*zone_omega[j][i];
			}

			mix->GetMWAndMoleFractionsFromMassFractions(mean_mw, mean_x, mean_omega);

			cout << " Mean molecular weight:          " << mean_mw << endl;
			cout << " Sum of mass fractions:          " << mean_omega.GetSumElements() << endl;
			cout << " Mean temperature (first guess): " << mean_temperature << endl;
			cout << " Total enthalpy:                 " << enthalpy_total << endl;
			cout << " Total mass:                     " << mass_total_initial << endl;

			double mean_enthalpy_temperature = mix->GetTemperatureFromMassEnthalpyAndMoleFractions(mean_temperature, enthalpy_total/mass_total_initial, mean_x, epsdtmax);

			cout << " Mean temperature (from enthalpy): " << mean_enthalpy_temperature << endl;
		}
	}
}

void OpenSMOKE_ICEM_MultiZone::SetMinimumValues(BzzVector &xMin)
{
	double minimum_volume		= 1.e-7;							// [m3]
	double minimum_temperature	= 200.;								// [K]
	double minimum_pressure		= 1.;								// [Pa]
	double minimum_g			= minimum_pressure*minimum_volume;	// [Pa.m3]
	double minimum_heat_release	= -1.e16;							// [J]
	
	if (icem_multizone_model == ICEM_CHEMKIN_MULTIZONE)
	{
		int k=1;
		for(int j=1;j<=N;j++)
		{
			for(int i=1;i<=NC;i++)
				xMin[k++] = ZERO;				// mass fractions [-]
			xMin[k++] = minimum_temperature;	// temperature [K]
			xMin[k++] = minimum_g;				// g [Pa.m3]
			xMin[k++] = minimum_pressure;		// pressure [Pa]
			xMin[k++] = minimum_heat_release;	// heat release [J]
		}
	}

	else if (icem_multizone_model == ICEM_KOMNINOS_SPLITTING)
	{
		int k=1;
		for(int i=1;i<=NC;i++)
			xMin[k++] = ZERO;				// mass fractions [-]
		xMin[k++] = minimum_temperature;	// temperature [K]
	}
}

void OpenSMOKE_ICEM_MultiZone::SetMaximumValues(BzzVector &xMax)
{
	double maximum_volume		= 1.e0;		// [m3]
	double maximum_temperature  = 6000.;	// [K]
	double maximum_pressure		= 101315e3;	// [Pa]
	double maximum_g			= maximum_pressure*maximum_volume;	// [Pa.m3]
	double maximum_heat_release	= 1.e16;							// [J]

	if (icem_multizone_model == ICEM_CHEMKIN_MULTIZONE)	
	{
		int k=1;
		for(int j=1;j<=N;j++)
		{
			for(int i=1;i<=NC;i++)
				xMax[k++] = ONE;				// mass fractions [-]
			xMax[k++] = maximum_temperature;	// temperature [K]
			xMax[k++] = maximum_g;				// g [Pa.m3]
			xMax[k++] = maximum_pressure;		// pressure [Pa]
			xMax[k++] = maximum_heat_release;	// heat release [J]
		}
	}

	else if (icem_multizone_model == ICEM_KOMNINOS_SPLITTING)	
	{
		int k=1;
		for(int i=1;i<=NC;i++)
			xMax[k++] = ONE;				// mass fractions [-]
		xMax[k++] = maximum_temperature;	// temperature [K]
	}

}

void OpenSMOKE_ICEM_MultiZone::SetInitialValues(BzzVector &xInitial)
{
	if (icem_multizone_model == ICEM_CHEMKIN_MULTIZONE)	
	{
		zone_omega_initial				= zone_omega;
		zone_temperature_initial		= zone_temperature;
		zone_enthalpy_initial			= zone_enthalpy;
		zone_g_initial					= zone_g;
		zone_pressure_initial			= zone_pressure;
		zone_heat_release_initial		= zone_heat_release;

		int k=1;
		for(int j=1;j<=N;j++)
		{
			for(int i=1;i<=NC;i++)
				xInitial[k++] = zone_omega[j][i];		// mass fractions [-]
			xInitial[k++] = zone_temperature[j];		// temperature [K]
			xInitial[k++] = zone_g[j];					// g [Pa.m3]
			xInitial[k++] = zone_pressure[j];			// pressure [Pa]
			xInitial[k++] = zone_heat_release[j];		// heat release [J]
		}

		work = 0.;
		volume_total_old = volume_total_starting;
	}

	else if (icem_multizone_model == ICEM_KOMNINOS_SPLITTING)	
	{
	}
}

void OpenSMOKE_ICEM_MultiZone::ReSetInitialValues(BzzVector &xInitial)
{
	if (icem_multizone_model == ICEM_CHEMKIN_MULTIZONE)	
	{
		int k=1;
		for(int j=1;j<=N;j++)
		{
			for(int i=1;i<=NC;i++)
				xInitial[k++] = zone_omega[j][i];			// mass fractions [-]
			xInitial[k++] = zone_temperature[j];			// temperature [K]
			xInitial[k++] = zone_g_initial[j];				// g [Pa.m3]
			xInitial[k++] = zone_pressure_initial[j];		// pressure [Pa]
			xInitial[k++] = zone_heat_release_initial[j];	// heat release [J]
		}

		work = 0.;
		volume_total_old = volume_total_starting;
	}

	else if (icem_multizone_model == ICEM_KOMNINOS_SPLITTING)	
	{
	}
}

void OpenSMOKE_ICEM_MultiZone::SetDifferentialEquations(BzzVectorInt &fDifferential)
{
	const int DIFFERENTIAL	= 1;
	const int ALGEBRAIC		= 0;

	if (icem_multizone_model == ICEM_CHEMKIN_MULTIZONE)
	{
		int k=1;
		for(int j=1;j<=N;j++)
		{
			for(int i=1;i<=NC;i++)
				fDifferential[k++] = DIFFERENTIAL;		// mass fractions [-]
			fDifferential[k++] = DIFFERENTIAL;			// temperature [K]
			fDifferential[k++] = ALGEBRAIC;				// g [Pa.m3]
			fDifferential[k++] = ALGEBRAIC;				// pressure [Pa]
			fDifferential[k++] = DIFFERENTIAL;			// heat release [J]
		}
	}

	else if (icem_multizone_model == ICEM_KOMNINOS_SPLITTING)
	{
		int k=1;
		for(int i=1;i<=NC;i++)
			fDifferential[k++] = DIFFERENTIAL;		// mass fractions [-]
		fDifferential[k++] = DIFFERENTIAL;			// temperature [K]
	}
}

void OpenSMOKE_ICEM_MultiZone::SetRelativeTolerances(BzzVector &fRelTol)
{
	const double omegaRelTol		= relativeTolerance*MachEpsFloat();
	const double temperatureRelTol	= relativeTolerance*MachEpsFloat();
	const double gRelTol			= relativeTolerance*MachEpsFloat();
	const double pressureRelTol		= relativeTolerance*MachEpsFloat();
	const double heatreleaseRelTol	= relativeTolerance*MachEpsFloat();

	if (icem_multizone_model == ICEM_CHEMKIN_MULTIZONE)
	{
		int k=1;
		for(int j=1;j<=N;j++)
		{
			for(int i=1;i<=NC;i++)
				fRelTol[k++] = omegaRelTol;		// mass fractions [-]
			fRelTol[k++] = temperatureRelTol;	// temperature [K]
			fRelTol[k++] = gRelTol;				// g [Pa.m3]
			fRelTol[k++] = pressureRelTol;		// pressure [Pa]
			fRelTol[k++] = heatreleaseRelTol;	// heat_release [J]
		}
	}

	else if (icem_multizone_model == ICEM_KOMNINOS_SPLITTING)
	{
		int k=1;
		for(int i=1;i<=NC;i++)
			fRelTol[k++] = omegaRelTol;		// mass fractions [-]
		fRelTol[k++] = temperatureRelTol;	// temperature [K]
	}
}
	
void OpenSMOKE_ICEM_MultiZone::SetEquations()
{
	if (icem_multizone_model == ICEM_CHEMKIN_MULTIZONE)
	{
		dimBlock = NC+4;
		nEquations = dimBlock*N;
	}

	else if (icem_multizone_model == ICEM_KOMNINOS_SPLITTING)
	{
		dimBlock = NC+1;
		nEquations = dimBlock;
	}
}

void OpenSMOKE_ICEM_MultiZone::Solve()
{
    double timeStart, timeEnd;

	BzzVector dummy;
	BzzVector xMin;
    BzzVector xMax;
    BzzVector xInitial;
    BzzVector fRelTol;
    BzzVectorInt fDifferential;

	PrepareFiles();
	PrepareAdditionalFiles();
    DAE_ICEM_MultiZone_Object.assignICEM_MultiZone(this);
    ODE_ICEM_MultiZone_Object.assignICEM_MultiZone(this);
    ODE_Reaction_ICEM_MultiZone_Object.assignICEM_MultiZone(this);
	SetEquations();
	ChangeDimensions(nEquations, &xMin);
	ChangeDimensions(nEquations, &xMax);
	ChangeDimensions(nEquations, &xInitial);
	ChangeDimensions(nEquations, &fDifferential);
	ChangeDimensions(nEquations, &fRelTol);

	SetInitialValues(xInitial);
	SetMinimumValues(xMin);
	SetMaximumValues(xMax);
	SetDifferentialEquations(fDifferential);
//	SetRelativeTolerances(fRelTol);
//	SetAbsoluteTolerances(fAbsTol);
	
	// ICEM_CHEMKIN_MULTIZONE
	if (icem_multizone_model == ICEM_CHEMKIN_MULTIZONE)
    {	
		cout << "StartCycles " << zone_temperature.Max() << endl;
		for(int n=1;n<=number_of_cycles;n++)
		{
			// History
			if (iHistory == true)
			{
				Tau_History = 0.;
				Csi_History = 0.;
				T_History = 0.;
				P_History = 0.;
				mass_History = 0.;
				mole_History = 0.;
			}
			cout << "Options" << endl;
			if ( iVerboseSensitivity == true )
			{
		//		BzzVector aux_x = zone_x.GetRow(1);
			//	sensitivity_fast->Initialize(0., zone_temperature[1], zone_pressure[1], zone_rho[1], zone_Cp[1], zone_volume[1], aux_x);
			}

			if (N==1)
			{
				dae_single(xInitial, 0., fDifferential, &DAE_ICEM_MultiZone_Object);
				dae_single.MyStepPrint();
				dae_single.SetMinimumConstraints(xMin);
				dae_single.SetMaximumConstraints(xMax);
				if (iAbsoluteTolerance == true)	dae_single.SetTollAbs(absoluteTolerance);
				if (iRelativeTolerance == true)	dae_single.SetTollRel(relativeTolerance*MachEpsFloat());
				
				{
					countGlobalIterations = -1;  // From -1 to avoid to store results from the first iteration
					nFileSteps		=   countFileSteps		= 1;
					nVideoSteps		=   countVideoSteps		= 1;
					countIterations = 0;

					timeStart = BzzGetCpuTime();
					TauTotal = 2.*Constants::pi/rotation_rate * (end_angle+start_angle)/Constants::pi/2.;
					dae_single(TauTotal,TauTotal);
					timeEnd = BzzGetCpuTime();
				 }

				if (iVerbose == true)
				{
					cout << endl;
					cout << "Number of function for Jacobian: " << dae_single.GetNumFunctionForJacobian()		<< endl;
					cout << "Numerical Jacobians: "				<< dae_single.GetNumNumericalJacobian()			<< endl;
					cout << "Time ODE solution: "				<< timeEnd - timeStart	<< " s"			<< endl << endl;
				}
			}
			else
			{
				dae_sparse(xInitial, 0., fDifferential, &DAE_ICEM_MultiZone_Object, dimBlock);
				dae_sparse.MyStepPrint();
				dae_sparse.SetMinimumConstraints(xMin);
				dae_sparse.SetMaximumConstraints(xMax);
				
				if (iAbsoluteTolerance == true)	
				{
					cout << "Absolute tolerance: " << absoluteTolerance << endl;
					dae_sparse.SetTollAbs(absoluteTolerance);
				}
				else
					cout << "Absolute tolerance (default): " << 1.e-10 << endl;

				if (iRelativeTolerance == true)	
				{
					cout << "Relative tolerance: " << relativeTolerance*MachEpsFloat() << endl;
					dae_sparse.SetTollRel(relativeTolerance*MachEpsFloat());
				}
				else 
					cout << "Relative tolerance (default): " << 100*MachEpsFloat() << endl;

				{
					countGlobalIterations = -1;  // From -1 to avoid to store results from the first iteration
					nFileSteps		=   countFileSteps		= 1;
					nVideoSteps		=   countVideoSteps		= 1;
					countIterations = 0;
					timeStart = BzzGetCpuTime();
					TauTotal = 2.*Constants::pi/rotation_rate * (end_angle+start_angle)/Constants::pi/2.;
					dae_sparse(TauTotal,TauTotal);
					timeEnd = BzzGetCpuTime();
				 }

				if (iVerbose == true)
				{
					cout << endl;
					cout << "Number of function for Jacobian: " << dae_sparse.GetNumFunctionForJacobian()		<< endl;
					cout << "Numerical Jacobians: "				<< dae_sparse.GetNumNumericalJacobian()			<< endl;
					cout << "Time ODE solution: "				<< timeEnd - timeStart	<< " s"			<< endl << endl;
				}
			}

			if (iHistory == true)
			{
				SaveOnBinaryFile(outputOSMName);
			}

			// PFR Exhaust gases
			BzzVector omega_exhaust_gases(NC);
			BzzMatrix omega_averaged(N,NC);

			if (iPlugFlowExhaustGases == true)
			{
				ErrorMessage("Plug-flow was disabled!");

				for (int j=1;j<=NC;j++)
					for (int i=1;i<=N;i++)
						omega_exhaust_gases[j] += zone_omega[i][j] * zone_mass[i];
				omega_exhaust_gases /= mass_total;
				double sum = omega_exhaust_gases.GetSumElements();
				omega_exhaust_gases /= sum; 

				OpenSMOKE_GasStream inlet_exhaust_gases;
				inlet_exhaust_gases.AssignKineticScheme(*mix);
				if (iPlugFlowExhaustGasesInletTemperature == false)
					inlet_exhaust_gases.AssignTemperature(zone_temperature[1], "K");
				else
					inlet_exhaust_gases.AssignTemperature(plugFlowExhaustGasesInletTemperature, "K");
				inlet_exhaust_gases.AssignPressure(zone_pressure[1], "Pa");
				inlet_exhaust_gases.AssignMassFlowRate(1.e-10, "kg/s");
				inlet_exhaust_gases.AssignMassFractions(omega_exhaust_gases);
				inlet_exhaust_gases.lock();

				OpenSMOKE_UD_Profile udf;
				if (iPlugFlowExhaustGasesInletTemperature == false)				
					udf.Setup("s", "K", 0., plugFlowExhaustGasesTime, zone_temperature[1], zone_temperature_initial[1]);
				else
					udf.Setup("s", "K", 0., plugFlowExhaustGasesTime, plugFlowExhaustGasesInletTemperature, zone_temperature_initial[1]);
				
				OpenSMOKE_PFR pfr;
				OpenSMOKE_GlobalKinetics global;
				OpenSMOKE_2EModel soot2EModel;

				pfr.SetOutputFolder(outputFolderName);
				pfr.AssignKineticScheme(*mix);
				pfr.AssignGlobalKineticScheme(global);
				pfr.AssignSoot2EModel(soot2EModel);
				
				pfr.AssignInletFlows(inlet_exhaust_gases);
				
				pfr.AssignEnergy(false);
				pfr.SetTimeIndependent();
				pfr.AssignEnd("s", plugFlowExhaustGasesTime);
				pfr.AssignDiameter("cm", 1.);
				pfr.SetUserDefinedTemperature(udf);
				pfr.Lock();

				pfr.VideoSummary();
				pfr.Solve();
			    pfr.VideoFinalResult();
				OpenSMOKE_GasStream outlet;
				pfr.OutletStream(outlet);
				omega_exhaust_gases = outlet.omega;
			}

			if (iRelaxation == false)
			{
				if (iPlugFlowExhaustGases == true)
				{
					ErrorMessage("Plug-flow was disabled!");
					for (int i=1;i<=N;i++)
						omega_averaged.SetRow(i,omega_exhaust_gases);
				}
				else
				{
					omega_averaged = zone_omega;
				}

				if (iResetSpecies == true)
					ResetSpecies(omega_averaged);
			}

			else if (iRelaxation == true)
			{
				relaxation_index++;
				
				if (iPlugFlowExhaustGases == true)
				{
					ErrorMessage("Plug-flow was disabled!");
					ChangeDimensions(N,NC, &relaxation_matrix[relaxation_index]);
					for (int i=1;i<=N;i++)
						relaxation_matrix[relaxation_index].SetRow(i,omega_exhaust_gases);
				}
				else
				{
					relaxation_matrix[relaxation_index] = zone_omega;
				}

				if (iResetSpecies == true)
					ResetSpecies(relaxation_matrix[relaxation_index]);
				
				if (relaxation_index<=relaxation_cycles)
				{
					for(int k=1;k<=relaxation_index;k++)
						for(int j=1;j<=N;j++)
							for(int i=1;i<=NC;i++)
								omega_averaged[j][i] += relaxation_matrix[k][j][i];
					omega_averaged /= double(relaxation_index);
				}
				else
				{
					for(int k=relaxation_index-relaxation_cycles+1;k<=relaxation_index;k++)
						for(int j=1;j<=N;j++)
							for(int i=1;i<=NC;i++)
								omega_averaged[j][i] += relaxation_matrix[k][j][i];
					omega_averaged /= double(relaxation_cycles);
				}
			}

			// Write on recycle file
			if (assignedExhaustRatio == true)
			{
				{
					double MWmix;
					BzzVector omega_averaged_mean(NC);
					for (int j=1;j<=NC;j++)
						for (int i=1;i<=N;i++)
							omega_averaged_mean[j] += omega_averaged[i][j] * zone_mass[i];
					omega_averaged_mean /= mass_total;

					BzzVector x_averaged_mean(NC);
					mix->GetMWAndMoleFractionsFromMassFractions(MWmix, x_averaged_mean, omega_averaged_mean);

					fRecycle << setw(5) << left << n;
					fRecycle << setw(16) << zone_temperature_initial[1];
					fRecycle << setw(16) << zone_pressure_initial[1]/101325.;
					for(int i=1;i<=NC;i++)
						fRecycle << setw(16) << left << x_averaged_mean[i];
					for(int i=1;i<=NC;i++)
						fRecycle << setw(16) << left << omega_averaged_mean[i];
					fRecycle << endl;
				}

				{
					double MWmix;
					BzzVector omega_averaged_mean(NC);
					for (int j=1;j<=NC;j++)
						for (int i=1;i<=N;i++)
							omega_averaged_mean[j] += zone_omega[i][j] * zone_mass[i];
					omega_averaged_mean /= mass_total;
					BzzVector x_averaged_mean(NC);
					mix->GetMWAndMoleFractionsFromMassFractions(MWmix, x_averaged_mean, omega_averaged_mean);
					
					fExhaust << setw(5)  << left << n;
					fExhaust << setw(16) << zone_temperature[1];
					fExhaust << setw(16) << zone_pressure[1]/101325.;
					for(int i=1;i<=NC;i++)
						fExhaust << setw(16) << left << x_averaged_mean[i];
					for(int i=1;i<=NC;i++)
						fExhaust << setw(16) << left << omega_averaged_mean[i];
					fExhaust << endl;
				}
			}

			for(int j=1;j<=N;j++)
			{
				cout << "Initial temperature: " << zone_temperature_initial[j] << "\t" 
					<<  "Zone temperature:    " << zone_temperature[j] << "\t";
	
				for(int i=1;i<=NC;i++)
					zone_omega[j][i] = (1.-exhaust_ratio)* zone_omega_initial[j][i] + 
										   exhaust_ratio * omega_averaged[j][i];

			// mass fractions [-]
				zone_enthalpy[j] = (1.-exhaust_ratio)* zone_enthalpy_initial[j] + 
									   exhaust_ratio * zone_enthalpy[j];						// enthalpy [J/kg]
				
				BzzVector aux_omega = zone_omega.GetRow(j);
				zone_temperature[j] = mix->GetTemperatureFromMassEnthalpyAndMassFractions(zone_temperature_initial[j], zone_enthalpy[j]/zone_mass[j], aux_omega, epsdtmax);		// temperature [K]
			//	zone_temperature[j] = zone_temperature_initial[j];
			}
			ReSetInitialValues(xInitial);
			fOutput << endl;
		}
	}

	// ICEM_KOMNINOS_SPLITTING
	if (icem_multizone_model == ICEM_KOMNINOS_SPLITTING)
    {	
		BzzOdeStiffObject *ode;
		ode = new BzzOdeStiffObject[N+1];

	//	ofstream fSplitting;
	//	openOutputFileAndControl(fSplitting, "Splitting.out");
	//	fSplitting.setf(ios::scientific);

		TauTotal = 2.*Constants::pi/rotation_rate;					// [s]
		double deltat = TauTotal/(360./delta_step_reference);		// [s]
		
		double t = -deltat;
		const int max_iterations = 50000;

		double timeReaction = 0.;
		timeStart = BzzGetCpuTime();
		for(int k=1;k<=max_iterations;k++)
		{
			int j;
			
			// Update time
			t += deltat;									// [s]
			teta = rotation_rate*t - start_angle;			// [rad]
			double teta_deg = teta*360./2./Constants::pi;	// [deg] 

			// Update Geometries
			UpdateAngleAndVolume(teta);
			UpdateVolumeDerivatives(teta);

			mix->GetMWAndMoleFractionsFromMassFractions(zone_mw, zone_x, zone_omega);
			double sum = 0.;
			for(j=1;j<=N;j++)
				sum += zone_volume[j]*zone_mw[j]/Constants::R_J_kmol/zone_temperature[j];
			zone_pressure = mass_total_initial/sum;

			for(j=1;j<=N;j++)	
				zone_mass[j] = zone_pressure[j]*zone_volume[j]*zone_mw[j]/Constants::R_J_kmol/zone_temperature[j];
			mass_total = zone_mass.GetSumElements();

			if (k%10 == 1)
				cout << k << "\t" << teta_deg << "\t" << zone_pressure[N]/101325. << "\t" << zone_temperature[N] << endl;

			// Update mass and pressure derivatives
			UpdateTimeDerivative(t);
			
			// Update Properties
			UpdatePropertiesMixing();
			
			// Update interface fluxes
			UpdateMassFluxes();
			UpdateInterfaces();
			UpdateEnthalpyFluxes();

			// Mixing
			for(j=1;j<=N;j++)
				for(int i=1;i<=NC;i++)
					zone_omega[j][i] += ( zone_m_in[j]/zone_mass[j] *(zone_omega_in[j][i]-zone_omega[j][i]) -
							              zone_m_in[j+1]/zone_mass[j]*(zone_omega_in[j+1][i]-zone_omega[j][i]) ) * deltat;

			// Temperature
			for(j=1;j<=N;j++)
				zone_temperature[j] += ( zone_m_in[j]*zone_h_in[j]-zone_m_in[j+1]*zone_h_in[j+1] -
										 zone_pressure[j]*dvolume_over_dt[j] ) / (zone_mass[j]*zone_Cv[j]) * deltat;
	
			if (icem_heat_transfer_model != ICEM_HEAT_TRANSFER_NONE)
			{
				// Update heat transfer coefficient
				UpdateHeatTransfer();

				// Internal thermal conduction
				zone_Qdiffusion = 0.;
		
				// Crevices
				if (zone_temperature[jCrevices] < Tambient)
					zone_temperature[jCrevices] += - zone_Qe[jCrevices]*area_crevices / 
												    (zone_mass[jCrevices]*zone_Cv[jCrevices]) * deltat;
				else
					zone_temperature[jCrevices]  = Tambient;

				// Outmost zone
				zone_temperature[jOutmost] +=  - zone_Qe[jOutmost]*area_cylinder_total / 
											    (zone_mass[jOutmost]*zone_Cv[jOutmost]) * deltat;
			}

			// Reaction
			double timeReactionStart = BzzGetCpuTime();
			for(j=1;j<=N;j++)
			{
				jGlobal = j;

				int k=1;
				for(int i=1;i<=NC;i++)
					xInitial[k++] = zone_omega[jGlobal][i];		// mass fractions [-]
				xInitial[k++] = zone_temperature[jGlobal];		// temperature [K]

				
				ode[jGlobal].Deinitialize();
				ode[jGlobal].SetInitialConditions(xInitial, t, &ODE_Reaction_ICEM_MultiZone_Object);

				ode[jGlobal].SetMinimumConstraints(xMin);
				ode[jGlobal].SetMaximumConstraints(xMax);

		//		if (iAbsoluteTolerance == true)	ode[jGlobal].SetTollAbs(absoluteTolerance);
		//		if (iRelativeTolerance == true)	ode[jGlobal].SetTollRel(fRelTol);

				ode[jGlobal](t+deltat,t+deltat);

				zone_heat_release[jGlobal] += zone_QReaction[jGlobal]*zone_volume[jGlobal]*deltat;
			}

			double timeReactionEnd = BzzGetCpuTime();
			timeReaction += (timeReactionEnd-timeReactionStart);

			if (k%5 == 1)
			{
			//	fSplitting << teta_deg << "\t" << zone_pressure[N]/101325. << "\t" << volume_total << "\t" << mass_total << "\t";
			//	for(j=1;j<=N;j++)
			//		fSplitting << zone_temperature[j] << "\t";
			//	fSplitting << endl;

				EquationSystemPrint(dummy, t);
			}


			// Update old values
			//zone_mass_old		= zone_mass;
			//timeOld				= t;

			     if (teta_deg <  -35.)						deltat = TauTotal/(360./delta_step_reference);
			else if (teta_deg >= -35. && teta_deg < -20.)	deltat = TauTotal/(360./(delta_step_reference/2.));
			else if (teta_deg >= -20. && teta_deg < -10.)	deltat = TauTotal/(360./(delta_step_reference/4.));
			else if (teta_deg >= -10. && teta_deg <  10.)	deltat = TauTotal/(360./(delta_step_reference/10.));
			else if (teta_deg >=  10. && teta_deg <  20.)	deltat = TauTotal/(360./(delta_step_reference/4.));
			else if (teta_deg >=  20. && teta_deg <  35.)	deltat = TauTotal/(360./(delta_step_reference/2.));
			else if (teta_deg >=  35. )						deltat = TauTotal/(360./delta_step_reference);

			if (t>TauTotal)	break;
		}
		timeEnd = BzzGetCpuTime();

		cout << " Time Total:    " << timeEnd - timeStart << "s" << endl;
		cout << " Time Reaction: " << timeReaction        << "s (" << timeReaction/(timeEnd-timeStart)*100. << "%)" << endl;
	}

	CloseFiles();
	if (iVerboseZones == true)		fZones.close();
	if (iVerboseGeometry == true)	fGeometry.close();
	if (iVerboseMasses == true)		fMass.close();
}

void OpenSMOKE_ICEM_MultiZone::ResetSpecies(BzzMatrix &species)
{
	if (species.Rows()!=N)		ErrorMessage("Wrong dimensions in Reset species (Rows)");
	if (species.Columns()!=NC)	ErrorMessage("Wrong dimensions in Reset species (Columns)");

	cout << "Resetting species: " << endl;
	for (int j=1;j<=indicesResetSpecies.Size();j++)
	{
		int k = indicesResetSpecies[j];
		cout << " * " << mix->names[k] << ": " << species[1][k] << endl; 
	}

	for (int i=1;i<=species.Rows();i++)
	{
		for (int j=1;j<=indicesResetSpecies.Size();j++)
		{
			int k= indicesResetSpecies[j];
			species[i][k] = 0.;
		}

		BzzVector aux = species.GetRow(i);
		double sum = aux.GetSumElements();
		for (int j=1;j<=NC;j++)
			species[i][j] /= sum;
	}

	cout << "Done!" << endl;
}

void OpenSMOKE_ICEM_MultiZone::DefineFromFile(const string inputFile)
{
    double			double_value;
    string			string_value;
    int				int_value;
	vector<string>  string_vector;

    OpenSMOKE_Dictionary_ICEM_MultiZone dictionary;
    dictionary.ParseFile(inputFile);

	// COMPULSORY: Model
	if (dictionary.Return("#Model", string_value))
	{
			  if (string_value == "CHEMKIN_MULTIZONE")		AssignEngineModel(ICEM_CHEMKIN_MULTIZONE);
		 else if (string_value == "KOMNINOS_MULTIZONE")		AssignEngineModel(ICEM_KOMNINOS_SPLITTING);
		 else ErrorMessage("Wrong engine model...");
	}

    // COMPULSORY: Rotation Rate
	if (dictionary.Return("#RotationRate", double_value, string_value))
		AssignRotationRate(string_value, double_value);

    // COMPULSORY: Compression Ratio
	if (dictionary.Return("#CompressionRatio", double_value))
		AssignCompressionRatio(double_value);

    // COMPULSORY: Reactor Clearance Volume
    if (dictionary.Return("#ClearanceVolume", double_value, string_value))
        AssignClearanceVolume(string_value, double_value);

    if (dictionary.Return("#TotalDisplacement", double_value, string_value))
        AssignTotalDisplacement(string_value, double_value);

	// COMPULSORY: Crevices volume
    if (dictionary.Return("#CrevicesVolume", double_value, string_value))
        AssignCrevicesVolume(string_value, double_value);

		// COMPULSORY: Crevices volume
    if (dictionary.Return("#NumberOfZones", int_value))
        AssignNumberOfZones(int_value);

    // COMPULSORY: Arm Ratio
    if (dictionary.Return("#ArmRatio", double_value))
        AssignArmRatio(double_value);

    // COMPULSORY: Start angle
    if (dictionary.Return("#StartAngle", double_value, string_value))
        AssignStartAngle(string_value, double_value);

	// COMPULSORY: End angle
    if (dictionary.Return("#EndAngle", double_value, string_value))
        AssignEndAngle(string_value, double_value);
	
	// COMPULSORY: Number of cycles
    if (dictionary.Return("#NumberOfCycles", int_value))
        AssignNumberOfCycles(int_value);

    // COMPULSORY: Exhaust ratio
	if (dictionary.Return("#ExhaustRatio", double_value))
		AssignExhaustRatio(double_value);

	// COMPULSORY: Start angle
    if (dictionary.Return("#Diameter", double_value, string_value))
        AssignDiameter(string_value, double_value);

	// OPTIONAL: KeySpecies
    if (dictionary.Return("#Key", string_vector))
        SetKeySpecies(string_vector);

	// OPTIONAL: Relaxation
	if (dictionary.Return("#Relaxation", int_value))
		SetRelaxation(int_value);

	// OPTIONAL: Plug Flow
	if (dictionary.Return("#PlugFlow", double_value, string_value))
        SetPlugFlowTime(string_value, double_value);
	if (dictionary.Return("#PlugFlowInletTemperature", double_value, string_value))
        SetPlugFlowInletTemperature(string_value, double_value);

	// OPTIONAL: KeySpecies
    if (dictionary.Return("#ResetSpecies", string_vector))
        SetResetSpecies(string_vector);

    // OPTIONAL: Output folder
    if (dictionary.Return("#Output", string_value))
        SetOutputFolder(string_value);

    // OPTIONAL: Print Options
    if (dictionary.Return("#Nvideo", int_value))
        SetVideoOptions(int_value);

    // OPTIONAL: Print Options
    if (dictionary.Return("#Nfile", int_value))
        SetFileOptions(int_value);

    // OPTIONAL: Print Options
    if (dictionary.Return("#NoVerbose"))
        UnsetVerbose();

	// OPTIONAL: Verbose Energy
    if (dictionary.Return("#VerboseEnergy"))
        SetVerboseEnergy();

    // OPTIONAL: Global kinetics
    if (dictionary.Return("#GlobalKinetics"))
        SetGlobalKinetics();

    // OPTIONAL: Relative tolerance
    if (dictionary.Return("#RelTol", double_value))
        SetRelativeTolerance(double_value);

    // OPTIONAL: Absolute tolerance
    if (dictionary.Return("#AbsTol", double_value))
        SetAbsoluteTolerance(double_value);

	// OPTIONAL: Heat flux
    if (dictionary.Return("#Qe", double_value, string_value))
        SetConstantHeatFlux(double_value, string_value);

	// OPTIONAL: Exchange area
    if (dictionary.Return("#A", double_value, string_value))
        SetConstantExchangeArea(double_value, string_value);

	// OPTIONAL: Exchange coefficient
    if (dictionary.Return("#U", double_value, string_value))
        SetConstantHeatExchangeCoefficient(double_value, string_value);

	// OPTIONAL: Ambient temperature
    if (dictionary.Return("#Twall", double_value, string_value))
        SetConstantAmbientTemperature(double_value, string_value);


	// OPTIONAL: UserDefined_Qe
    if (dictionary.Return("#UserDefined_Qe", string_value))
        SetUserDefinedHeatFlux(string_value);

	// OPTIONAL: UserDefined_Ae
    if (dictionary.Return("#UserDefined_A", string_value))
        SetUserDefinedExchangeArea(string_value);

	// OPTIONAL: UserDefined_U
    if (dictionary.Return("#UserDefined_U", string_value))
        SetUserDefinedHeatExchangeCoefficient(string_value);

	// OPTIONAL: UserDefined_Twall
    if (dictionary.Return("#UserDefined_Twall", string_value))
        SetUserDefinedAmbientTemperature(string_value);
 

	// OPTIONAL: TwoEquationModel
    if (dictionary.Return("#2E_Model"))
        SetTwoEquationModel();

	// OPTIONAL: Reaction Rates
    if (dictionary.Return("#ReactionRates", string_vector))
		SetReactionRatesOnFile(string_vector);

	// OPTIONAL: Formation Rates
    if (dictionary.Return("#FormationRates", string_vector))
		SetFormationRatesOnFile(string_vector);

	// OPTIONAL: Rate of Production Analysis
    if (dictionary.Return("#ROPA"))
		SetROPAOnFile();

	// OPTIONAL: Rate of Production Analysis
    if (dictionary.Return("#VerboseROPA", string_vector))
		SetVerboseROPAOnFile(string_vector);

	// OPTIONAL: Sensitivity Analysis
    if (dictionary.Return("#SensitivityOptions", string_vector))
		SetSensitivityOptions(string_vector);

	// OPTIONAL: Sensitivity Analysis
    if (dictionary.Return("#Sensitivity", string_vector))
		SetSensitivityOnFile(string_vector);

	// OPTIONAL: Verbose Zones
    if (dictionary.Return("#VerboseZones", string_vector))
		SetVerboseZones(string_vector);

	// OPTIONAL: Verbose geometry
    if (dictionary.Return("#VerboseGeometry"))
		SetVerboseGeometry();

	// OPTIONAL: Verbose masses
    if (dictionary.Return("#VerboseMasses"))
		SetVerboseMasses();

	// OPTIONAL: Heat flux
    if (dictionary.Return("#HeatTransfer", string_value))
	{
		     if (string_value == "NONE")	    { }
		else if (string_value == "WOSCHNI")		SetHeatTransferModel_Woschni();
		else if (string_value == "ASSANIS")		SetHeatTransferModel_Assanis();
		else if (string_value == "HOHENBERG")	SetHeatTransferModel_Hohenberg();
		else if (string_value == "ANNAND")		SetHeatTransferModel_Annand();
		else ErrorMessage("Wrong heat transfer model...");
	}

	// OPTIONAL: Initial conditions
	if (dictionary.Return("#Zones", string_vector))
		SetZoneInitialConditions(string_vector);

	// OPTIONAL: Initial conditions
	if (dictionary.Return("#GaussianZones", string_vector))
		SetGaussianInitialConditions(string_vector);

	// OPTIONAL: Cold simulation
	if (dictionary.Return("#Cold"))	SetColdSimulation();

	// OPTIONAL: Exchange area model
	if (dictionary.Return("#ExchangeArea", string_value))
	{
			  if (string_value == "LATERAL")	SetExchangeAreaModel(ICEM_EXCHANGE_AREA_LATERAL);
		 else if (string_value == "TOTAL")		SetExchangeAreaModel(ICEM_EXCHANGE_AREA_TOTAL);
		 else ErrorMessage("Wrong exchange area model...");
	}

	// OPTIONAL: Heat transfer model parameters
	if (dictionary.Return("#WoschniParameters", string_vector))		SetHeatTransferModel_WoschniParameters(string_vector);
	if (dictionary.Return("#AssanisParameters", string_vector))		SetHeatTransferModel_AssanisParameters(string_vector);
	if (dictionary.Return("#HohenbergParameters", string_vector))	SetHeatTransferModel_HohenbergParameters(string_vector);
	if (dictionary.Return("#AnnandParameters", string_vector))		SetHeatTransferModel_AnnandParameters(string_vector);
		
	if (dictionary.Return("#MonteCarloMixing"))
		SetMonteCarloMixing();

	if (dictionary.Return("#TauMixing", double_value, string_value))
        SetMicromixingTime(double_value, string_value);

	if (dictionary.Return("#Cphi", double_value))
        SetMicromixingConstant(double_value);

	if (dictionary.Return("#PostProcessing"))
        SetPostProcessing();

	Lock();
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
//																								   //
//									OUTPUT FUNCTIONS											   //
//																								   //
/////////////////////////////////////////////////////////////////////////////////////////////////////

void OpenSMOKE_ICEM_MultiZone::LabelODEFile()
{
	LabelEquationSystemFile();
}

void OpenSMOKE_ICEM_MultiZone::LabelEquationSystemFile()
{
    int i;

    if (iVerbose == true)
    {
        {
			fOutput << "Tau[s](1)     "		<< "\t";
            fOutput << "Angle[deg](2) "		<< "\t";
			fOutput << "T[K](3)       "		<< "\t";
			fOutput << "P[atm](4)     "		<< "\t";
			fOutput << "Utot[J](5)    "		<< "\t";
			fOutput << "Htot[J](6)    "		<< "\t";
			fOutput << "HRR[J/s](7)   "		<< "\t";
			fOutput << "HR[J](8)      "		<< "\t";
            fOutput << "QE[?](9)      "		<< "\t";
            fOutput << "rh[kg/m3](10) "	    << "\t";
            fOutput << "MW(11)        "	    << "\t";
			fOutput << "H[mm](12)     "		<< "\t";
            fOutput << "V[cm3](13)    "		<< "\t";
            fOutput << "L[J](14)      "		<< "\t";

			if(key_species_names.size() == 1)
				fOutput << "Eta_" << key_species_names[0] << "(14)" << "\t";
			else
	            fOutput << "dummy(15)    "	    << "\t";

            int fOutputCount = 16;
            for (i=1;i<=NC;i++)
                fOutput << mix->names[i] << " (x" << fOutputCount++ << ")  \t";
            for (i=1;i<=NC;i++)
                fOutput << mix->names[i] << " (w" << fOutputCount++ << ")  \t";

            fOutput << endl;
            fOutput << endl;
        }

        if (iTwoEquationModel == true)
		{
            fSoot2E	<< setw(20) << left << "Tau[s](1)";
            soot2EModel->GnuPlotInterface(fSoot2E, 4);

			fSoot2E << endl << endl;
		}
/*
        if (mix->pah_manager.iPAH == 1)
		{
			fPAH << "Tau[s](1)      "		<< "\t";
            fPAH << "Angle[dg](2)   "		<< "\t";
            fPAH << "T[K](3)        "		<< "\t";
			mix->pah_manager.GnuPlotInterface(fPAH, 3);
			fPAH << endl << endl;
		}
		*/
 /*       if (mix->soot_manager.iSoot == 1)
        {
			fSoot << "Tau[s](1)      "		<< "\t";
            fSoot << "Angle[dg](2)   "		<< "\t";
            fSoot << "T[K](3)        "		<< "\t";
            mix->soot_manager.GnuPlotInterface(fSoot, 4);
            fSoot << endl << endl;
        }
*/
		if (iKeySpecies == true && key_species_names.size()>1)
        {
			fConversions << "Tau[s](1)      "		<< "\t";
            fConversions << "Angle[dg](2)   "		<< "\t";
            fConversions << "T[K](3)        "		<< "\t";
			
			int fOutputCount = 4;
            for(i=0;i<int(key_species_names.size());i++)
	            fConversions << "Eta_" << key_species_names[i] << "(" << fOutputCount++ << ")  \t";
			fConversions << endl << endl;
        }

		if (iVerboseReactionRates == true)
		{
			fReactionRates << "Tau[s](1)      "		<< "\t";
            fReactionRates << "Angle[dg](2)   "		<< "\t";
            fReactionRates << "T[K](3)        "		<< "\t";
			
			int countFileOutput = 4;
            for (i=1;i<=index_reaction_rates.Size();i++)
                fReactionRates << "R_" << index_reaction_rates[i] << "[" << names_reaction_rates[i-1] << "](" << countFileOutput++ << ")  \t";
			fReactionRates << endl;
			fReactionRates << endl;
		}

		if (iVerboseFormationRates == true)
		{
			fFormationRates << "Tau[s](1)      "		<< "\t";
            fFormationRates << "Angle[dg](2)   "		<< "\t";
            fFormationRates << "T[K](3)        "		<< "\t";
			
			int countFileOutput = 4;
            for (i=1;i<=index_formation_rates.Size();i++)
                fFormationRates << names_formation_rates[i-1] << "(" << countFileOutput++ << ")  \t";
			fFormationRates << endl;
			fFormationRates << endl;
		}
    }
}

void OpenSMOKE_ICEM_MultiZone::VideoGeneralInfo()
{
    cout.setf(ios::scientific);
    cout << endl;

    cout << "xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx" << endl;
    cout << "xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx" << endl;
    cout << endl;

    cout << "---------------------------------------------------------------------" << endl;
    cout << "Internal Combustion Engine Model Summary" << endl;
    cout << "---------------------------------------------------------------------" << endl;
	cout << "Angle:\t\t"				<< teta*rad_to_degrees				<< " [deg]"		<< endl;
    cout << "Volume:\t\t\t"				<< volume_total*1.e6				<< " [cm3]"		<< endl;
    cout << "Starting volume:\t\t"		<< volume_total_starting*1.e6		<< " [cm3]"		<< endl;
	cout << "Mass:\t\t\t"				<< mass_total*1.e3					<< " [g]"		<< endl;
    cout << "Moles:\t\t\t"				<< mass_total/zone_mw[N]			<< " [kmol]"	<< endl;
	cout << "Height:\t\t\t"				<< height_cylinder*100.				<< " [cm]"		<< endl;
	cout << "Temperature:\t\t"			<< zone_temperature[N]				<< " [K]"		<< endl;
	cout << "Pressure:\t\t"				<< zone_pressure[N]					<< " [Pa]"		<< endl;
    cout << "Density:\t\t"				<< zone_rho[N]						<< " [kg/m3]"	<< endl;
    cout << "Molecular weight:\t"		<< zone_mw[N]						<< " [kg/kmol]" << endl;
    cout << "Mean velocity:\t\t"		<< velocity_piston					<< " [m/s]"		<< endl;
    cout << "Arm length (La):\t"		<< arm_La*100.						<< " [cm]"		<< endl;
    cout << "Arm length (Lc):\t"		<< arm_Lc*100.						<< " [cm]"		<< endl;
	cout <<  endl;
}

void OpenSMOKE_ICEM_MultiZone::VideoSummary()
{
    VideoGeneralInfo();

    cout << "#\tName\tx\t\tomega" << endl;
    cout << "---------------------------------------------------------------------" << endl;
    for (int i=1;i<=NC;i++)
        if (inletStream->omega[i]!=0.)
            cout << i << "\t\t" << mix->names[i] << "\t" << inletStream->x[i] << "\t" << inletStream->omega[i] << endl;
    cout << endl;
    cout << "xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx" << endl;
    cout << "xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx" << endl;
    cout << endl;
}

void OpenSMOKE_ICEM_MultiZone::VideoFinalResult()
{
    int i;

    // Vectors
    BzzVector conversion(NC);

    // Print General Information about Batch
    VideoGeneralInfo();

    // Calculate Conversion of Inlet Species
    for (i=1;i<=NC;i++)
        if (inletStream->omega[i]!=0.)
            conversion[i] = (1. - mean_omega[i]/inletStream->omega[i])*100.;

    // Final Composition
    cout << "#\tName\tx\t\tomega" << endl;
    cout << "---------------------------------------------------------------------" << endl;
    for (i=1;i<=NC;i++)
        if (mean_x[i]!=0.)
            cout << i << "\t" << mix->names[i] << "\t" << mean_x[i] << "\t" << mean_omega[i] << endl;
    cout << endl;

    // Conversion
    cout << "#\tName\tomegaInitial\tomegaFinal\tConversion(%)" << endl;
    cout << "---------------------------------------------------------------------" << endl;
    for (i=1;i<=NC;i++)
        if (conversion[i]>0.)
            cout << i << "\t" << mix->names[i] << "\t"
            << inletStream->omega[i] << "\t" << mean_omega[i] << "\t" << conversion[i] << endl;
    cout << endl;

    cout << "xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx" << endl;
    cout << "xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx" << endl;
    cout << endl;
}

void OpenSMOKE_ICEM_MultiZone::SummaryOnFile()
{
}

void OpenSMOKE_ICEM_MultiZone::ODEPrint(BzzVector &y, double t)
{
}

void OpenSMOKE_ICEM_MultiZone::EquationSystemPrint(BzzVector &y, double t)
{
	double dummy = ZERO;
	
	// Work
	work += zone_pressure[N]*(volume_total-volume_total_old);

	double mean_temperature;
    if (iVerbose == true)
    {
		double teta_deg = teta*rad_to_degrees;
		double mean_enthalpy_temperature;

		// Energy additional information
		//if ( countFileSteps == nFileSteps )
		{
			int j;
			for(j=1;j<=N;j++)
			{
				BzzVector aux_omega = zone_omega.GetRow(j);

				// Enthalpy
				mix->GetMixAveragedEnthalpy_Mass(dummy_vector_nc, zone_temperature[j]);				// [J/kg]
				zone_hmass.SetRow(j, dummy_vector_nc);												// [J/kg]
				zone_enthalpy[j] = zone_mass[j]*Dot(dummy_vector_nc,aux_omega);			// [J]
			
				// Internal energy
				mix->GetMixAveragedInternalEnergy_Mass(dummy_vector_nc, zone_temperature[j]);		// [J/kg]
				zone_umass.SetRow(j, dummy_vector_nc);												// [J/kg]
				zone_internal_energy[j] = zone_mass[j]*Dot(dummy_vector_nc,aux_omega);	// [J]
			}

			// Total enthalpy and internal energy
			enthalpy_total = zone_enthalpy.GetSumElements();									// [J]
			internal_energy_total = enthalpy_total;												// [J]
			for(j=1;j<=N;j++)
				internal_energy_total -= zone_pressure[j]/zone_rho[j]*zone_mass[j];				// [J]
		}

		// Mean properties
		if ( countFileSteps == nFileSteps || countVideoSteps == nVideoSteps || iHistory == true)
		{
			int j;
			BzzVector mass_fraction(N);
			mean_omega = 0.;
			for(j=1;j<=N;j++)
			{
				mass_fraction[j] = zone_mass[j]/mass_total;
				for(int i=1;i<=NC;i++)
					mean_omega[i] += mass_fraction[j]*zone_omega[j][i];
			}

			mix->GetMWAndMoleFractionsFromMassFractions(mean_mw, mean_x, mean_omega);
			mean_rho = mass_total/volume_total;
			QReaction_total = Dot(zone_QReaction, zone_volume);			// Total heat release rate [J/s]
			HeatRelease_total = zone_heat_release.GetSumElements();		// Total heat release [J]

			mean_temperature = Dot(mass_fraction, zone_temperature);

			//cout << " Sum of mass fractions:            " << mean_omega.GetSumElements() << endl;
			//cout << " Mean temperature (first guess):   " << mean_temperature << endl;
			//cout << " Total enthalpy:                   " << enthalpy_total << endl;
			//cout << " Total mass:                       " << mass_total << endl;

			mean_enthalpy_temperature = mix->GetTemperatureFromMassEnthalpyAndMoleFractions(mean_temperature, enthalpy_total/mass_total, mean_x, epsdtmax);
			//cout << " Mean temperature (from enthalpy): " << mean_enthalpy_temperature << endl;
		}

        if ( countIterations%(300*nVideoSteps) == 0)
        {
            cout    << endl
                    << setw(8)  << left << "#"
                    << setw(16) << left << "Angle[deg]"
                    << setw(16) << left << "TM[K]"
                    << setw(16) << left << "T1[K]"
                    << setw(16) << left << "TN[K]"
					<< setw(16) << left << "P[atm]";
                    
			cout	<< endl;
        }

        if (countVideoSteps == nVideoSteps)
        {
            cout	<< setw(8)  << left << countIterations
                    << setw(16) << left << teta_deg
                    << setw(16) << left << mean_enthalpy_temperature
                    << setw(16) << left << zone_temperature[1]
                    << setw(16) << left << zone_temperature[N]
                    << setw(16) << left << zone_pressure[N]/101325.;


            // Soot Two Equation Model
            if (iTwoEquationModel == true)
            {
                cout 	<< setw(16)  << left << soot2EModel->m0
                        << setw(16)  << left << soot2EModel->fv;
            }

            cout	<< endl;

            countVideoSteps = 0;
        }

        if (countFileSteps == nFileSteps)
        {
            int i;

            fOutput << t										<< "\t"			// [s]
                    << teta_deg									<< "\t"			// [deg]
                    << mean_enthalpy_temperature				<< "\t"			// [K]
                    << zone_pressure[N]/101325.					<< "\t"			// [atm]
                    << internal_energy_total					<< "\t"			// [J]
					<< enthalpy_total							<< "\t"			// [J]
					<< QReaction_total							<< "\t"			// [J/s]
					<< HeatRelease_total						<< "\t"			// [J]
					<< Qe*area_cylinder_total/volume_total     	<< "\t"			// [?]
                    << mean_rho									<< "\t"			// [kg/m3]
                    << mean_mw									<< "\t"			// [kg/kmol]
					<< height_cylinder*1000.					<< "\t"			// [mm]
                    << volume_total*1.e6						<< "\t"			// [cm3]
                    << work										<< "\t";		// [J]
			
			if(key_species_names.size() == 1)
				fOutput << 1.-mean_omega[key_species_index[0]]/inletStream->omega[key_species_index[0]]<< "\t";
			else
	            fOutput << dummy					<< "\t";

            for (i=1;i<=NC;i++)
                fOutput << mean_x[i]				<< "\t";
            for (i=1;i<=NC;i++)
                fOutput << mean_omega[i]			<< "\t";
            fOutput << endl;

			// Soot Two Equation Model
            if (iTwoEquationModel == true)
			{
				fSoot2E	<< setw(20) << left << t;
				soot2EModel->write_on_file(fSoot2E, y[indexTwoEquations], y[indexTwoEquations+1]);
			}

    /*        if (mix->soot_manager.iSoot == 1)
            {
                // Soot Distribution data
                mix->soot_manager.large_calculateAll(omega, x, c, rho);
                mix->soot_manager.small_calculateAll(omega, x, c, rho);

                fSoot << t			<< "\t";
                fSoot << teta_deg	<< "\t";
                fSoot << T			<< "\t";
                mix->soot_manager.print_main_data_on_file(fSoot);
                fSoot << endl;
            }
			*/
 /*           if (mix->pah_manager.iPAH == 1)
            {
                // PAH Distribution data
                mix->pah_manager.pah_calculateAll(omega, x, c, rho);

                fPAH << t			<< "\t";
                fPAH << teta_deg	<< "\t";
                fPAH << T			<< "\t";
                mix->pah_manager.print_main_data_on_file(fPAH);
                fPAH << endl;
            }
*/
			if (iKeySpecies == true && key_species_names.size()>1)
            {
                fConversions << t			<< "\t";
                fConversions << teta_deg	<< "\t";
                fConversions << T			<< "\t";

				for(i=0;i<int(key_species_names.size());i++)
	                fConversions << 1.-omega[key_species_index[i]]/inletStream->omega[key_species_index[i]]	<< "\t";
	        
                fConversions << endl;
            }

			if (iVerboseReactionRates == true)
			{
				fReactionRates	<< t		<< "\t"
								<< teta_deg	<< "\t"
								<< T		<< "\t";
				
				for (i=1;i<=index_reaction_rates.Size();i++)
					fReactionRates	<< mix->r[index_reaction_rates[i]]	<< "\t";
				
				fReactionRates	<< endl;
			}

			if (iVerboseFormationRates == true)
			{
				fFormationRates	<< t		<< "\t"
								<< teta_deg	<< "\t"
								<< T		<< "\t";
				
				for (i=1;i<=index_formation_rates.Size();i++)
					fFormationRates	<< R[index_formation_rates[i]]	<< "\t";
				
				fFormationRates	<< endl;
			}

			PrintAdditionalFiles(t, teta_deg);

			if (iVerboseSensitivity == true)
			{
				if ( sensitivity_fast->Integration() == OPENSMOKE_SENSITIVITYANALYSIS_FAST_UNSTEADY0D_EULER_IMPLICIT || 
					 sensitivity_fast->Integration() == OPENSMOKE_SENSITIVITYANALYSIS_FAST_UNSTEADY0D_EULER_EXPLICIT   )
				{
					double tStart,tEnd;
					BzzMatrix J;
					BzzVector aux_x = zone_x.GetRow(1);
					dae_single.GetLastJacobian(&J);
					sensitivity_fast->Update(t, J, zone_temperature[1], zone_pressure[1], zone_rho[1], zone_Cp[1], zone_volume[1], aux_x);
				}
				else if (sensitivity_fast->Integration() == OPENSMOKE_SENSITIVITYANALYSIS_FAST_UNSTEADY0D_ACCURATE)
				{
					ErrorMessage("Accurate Sensitivity Analysis not yet implemented...");
					//sensitivity_fast->UpdateOutputFile(Tau,y);
				}
			}

            countFileSteps = 0;
        }
    }

	// Update old values
	volume_total_old	= volume_total;
	zone_volume_old		= zone_volume;
	zone_mass_old		= zone_mass;
	timeOld				= t;
	
	// Update counters
    countVideoSteps++;
    countFileSteps++;
    countIterations++;

	// History
	if (iHistory == true)
    {
        countGlobalIterations++;

        if(countGlobalIterations>0)
        {
			int index = countGlobalIterations;

			if (countGlobalIterations>MAX_TIME_STEPS)
				index = MAX_TIME_STEPS;    
			Tau_History[index] = t;
			Csi_History[index] = teta*rad_to_degrees;
			P_History[index] = zone_pressure[N];
			T_History[index] = mean_temperature;
			mass_History.SetRow(index, mean_omega);
			mole_History.SetRow(index, mean_x);			
        }
    }
}

void OpenSMOKE_ICEM_MultiZone::EnergyAnalysis(OpenSMOKE_GasStream &outletStream)
{
	double H_mass_initial  = inletStream->massSpecificEnthalpy;			// [J/kg]
	double P_mass_initial  = inletStream->massSpecificPressureEnergy;	// [J/kg]
	double U_mass_initial  = inletStream->massSpecificInternalEnergy;	// [J/kg]
	
	double H_mass_final  = outletStream.massSpecificEnthalpy;			// [J/kg]
	double P_mass_final  = outletStream.massSpecificPressureEnergy;		// [J/kg]
	double U_mass_final  = outletStream.massSpecificInternalEnergy;		// [J/kg]


	cout << endl;
	cout << "---------------------------------------------------------------------" << endl;
	cout << "                         ENERGY ANALYSIS                             " << endl;
	cout << "---------------------------------------------------------------------" << endl;
	cout << "Initial - Enthalpy:        " << H_mass_initial		<< " [J/kg]" << endl;
	cout << "Initial - Pressure Energy: " << P_mass_initial		<< " [J/kg]" << endl;
	cout << "Initial - Internal energy: " << U_mass_initial		<< " [J/kg]" << endl;
	cout << "Final   - Enthalpy:        " << H_mass_final		<< " [J/kg]" << endl;
	cout << "Final   - Pressure Energy: " << P_mass_final		<< " [J/kg]" << endl;
	cout << "Final   - Internal energy: " << U_mass_final		<< " [J/kg]" << endl;
	cout << endl;

//	cout << "Initial energy:     "	<< U_mass_initial*mass								<< " [J]" << endl;
//	cout << "Final energy:       "	<< U_mass_final*mass								<< " [J]" << endl;
//	cout << "Balance (IN-OUT):   "  << U_mass_initial*mass - U_mass_final*mass			<< " [J]" << endl;
}

void OpenSMOKE_ICEM_MultiZone::UpdateTimeDerivative(const double t)
{
	double deltat = t-timeOld;

	if (deltat < 1.e-16)	
		dmass_over_dt = 0.;
	else
		for(int j=1;j<=N;j++)
			dmass_over_dt[j] = (zone_mass[j] - zone_mass_old[j]) / deltat;
}

void OpenSMOKE_ICEM_MultiZone::UpdateMassFluxes()
{
	zone_m_in[N]	 = dmass_over_dt[N];

	for(int j=N-1;j>=2;j--)
		zone_m_in[j] = zone_m_in[j+1]+dmass_over_dt[j];

	zone_m_in[1]	 = 0.;
}

void OpenSMOKE_ICEM_MultiZone::UpdateInterfaces()
{
	zone_omega_in  = 0.;
	zone_hmass_in  = 0.;

	for(int j=2;j<=N;j++)
	{
		if (zone_m_in[j]>=0.)
			for(int i=1;i<=NC;i++)
			{
				zone_omega_in[j][i] = zone_omega[j-1][i];
				zone_hmass_in[j][i] = zone_hmass[j-1][i];
			}
		else
			for(int i=1;i<=NC;i++)
			{
				zone_omega_in[j][i] = zone_omega[j][i];
				zone_hmass_in[j][i] = zone_hmass[j][i];
			}
	}

	{
	int j;

	// Thermal conductivity at interfaces
	zone_lambda_in[1] = 0.;
	for(j=2;j<=N-1;j++)
		zone_lambda_in[j] = 0.50*(zone_lambda[j-1]+zone_lambda[j]);
	zone_lambda_in[N]   = zone_lambda[N-1] + (zone_lambda[N]-zone_lambda[N-1])*thickness/(thickness+height_e[N]);
	zone_lambda_in[N+1] = zone_lambda[N-1] + (zone_lambda[N]-zone_lambda[N-1])*thickness/(thickness+diameter_e[N]);

	// Temperature gradient (interfaces)
	zone_temperature_gradient_in[1] = 0.;
	for(j=2;j<=N-1;j++)
		zone_temperature_gradient_in[j] = (zone_temperature[j-1]-zone_temperature[j])/thickness;
	zone_temperature_gradient_in[N]   = (zone_temperature[N-1]-zone_temperature[N])/0.50/(thickness+height_e[N]);
	zone_temperature_gradient_in[N+1] = (zone_temperature[N-1]-zone_temperature[N])/0.50/(thickness+diameter_e[N]);

	// Turbulent Distance from wall
	double uStar = 0.10*velocity_piston;
	double muWall = zone_mu[jCrevices];
	yPlus_in[jCrevices]  = 0.;
	yPlus_in[jOutmost]   = zone_rho[jOutmost] * (0.50*thickness);
	yPlus_in[jOutmost+1] = zone_rho[jOutmost] * thickness;
	for(j=jOutmost+2;j<=N;j++)
		yPlus_in[j] = yPlus_in[j-1] + zone_rho[j-1] * thickness;
	yPlus_in[N+1] = yPlus_in[N];
	for(j=1;j<=N+1;j++)
		yPlus_in[j] *= (uStar/muWall);

	// Viscosity ratio
	for(j=1;j<=N+1;j++)
		viscosity_ratio_in[j] = 0.41*yPlus_in[j]*(1.-exp(-2.*0.06*0.41*yPlus_in[j]));

	// Prandtl ratio
	BzzVector Prandtl_ratio(N+1);
	for(j=1;j<=N;j++)
		Prandtl_ratio[j] = zone_mu[j]*zone_Cp[j]/zone_lambda_in[j];
	Prandtl_ratio[N+1] = Prandtl_ratio[N];

//	for(j=1;j<=N+1;j++)
//		zone_lambda_in *= (1.+viscosity_ratio_in[j]);	// TODO

//	cout << muWall << "\t" << uStar << endl;
//	cout << yPlus_in[jOutmost] << " " << yPlus_in[N/2] << endl;
//	cout << Prandtl_ratio[jOutmost] << " " << Prandtl_ratio[N/2] << endl;
	//cout << viscosity_ratio_in[jOutmost] << " " << viscosity_ratio_in[N/2] << " " << viscosity_ratio_in[N] << " " << viscosity_ratio_in[N] << endl;
	}
}

void OpenSMOKE_ICEM_MultiZone::UpdateEnthalpyFluxes()
{
	zone_h_in	= 0.;

	for(int j=2;j<=N;j++)
		for(int i=1;i<=NC;i++)
//			zone_h_in[j] += zone_omega_in[j][i]*(zone_hmass_in[j][i]-zone_hmass[j][i]);
			zone_h_in[j] += zone_omega_in[j][i]*(zone_hmass_in[j][i]-zone_umass[j][i]);
}

void OpenSMOKE_ICEM_MultiZone::UpdateMeanMassFractions()
{
	mean_omega = 0.;
	for(int j=1;j<=N;j++)
	{
		double mass_fraction = zone_mass[j]/mass_total;
		for(int i=1;i<=NC;i++)
			mean_omega[i] += mass_fraction*zone_omega[j][i];
	}
}

void OpenSMOKE_ICEM_MultiZone::Equations_Species()
{
	if (icem_multizone_model == ICEM_CHEMKIN_MULTIZONE)
	{
		if (iMonteCarloMixing == true)
		{
			UpdateMeanMassFractions();
			for(int j=1;j<=N;j++)
				for(int i=1;i<=NC;i++)
					domega_over_dt[j][i] = zone_R[j][i]/zone_rho[j] - Cphi/TauMixing*(zone_omega[j][i]-mean_omega[i]);
		}
		else
		{
			for(int j=1;j<=N;j++)
				for(int i=1;i<=NC;i++)
					domega_over_dt[j][i] = zone_R[j][i]/zone_rho[j];
		}
	}

	else if (icem_multizone_model == ICEM_KOMNINOS_SPLITTING)
	{
		for(int j=1;j<=N;j++)
			for(int i=1;i<=NC;i++)
				domega_over_dt[j][i] = zone_m_in[j]/zone_mass[j] *(zone_omega_in[j][i]-zone_omega[j][i]) -
									   zone_m_in[j+1]/zone_mass[j]*(zone_omega_in[j+1][i]-zone_omega[j][i]) +
									   zone_R[j][i]/zone_rho[j];
	}
}

void OpenSMOKE_ICEM_MultiZone::Equations_Energy()
{
	if (icem_multizone_model == ICEM_CHEMKIN_MULTIZONE)
	{
		dtemperature_over_dt[1] = (zone_UReaction[1]*zone_volume[1] - zone_pressure[1]*dvolume_over_dt[1] - zone_Qe[1]*zone_area_effective[1])/(zone_mass[1]*zone_Cv[1]);

		for(int j=2;j<=N;j++)
			dtemperature_over_dt[j] =	(zone_UReaction[j]*zone_volume[j] - zone_pressure[j]*dvolume_over_dt[j]- zone_Qe[j]*zone_area_effective[j])/(zone_mass[j]*zone_Cv[j]);
	}

	else if (icem_multizone_model == ICEM_KOMNINOS_SPLITTING)
	{
		ErrorMessage("Not available for Komninos splitting");
	}
}

void OpenSMOKE_ICEM_MultiZone::Equations_Pressure()
{
	for(int j=1;j<=N-1;j++)
		p_equation[j]  = zone_pressure[j] - zone_pressure[j+1];
	p_equation[N]  = zone_pressure[N] - zone_g[N]/volume_total;

	g_equation[1] = zone_g[1] - zone_mass[1]/zone_mw[1]*Constants::R_J_kmol*zone_temperature[1];
	for(int j=2;j<=N;j++)
		g_equation[j]  = zone_g[j] - zone_g[j-1] - zone_mass[j]/zone_mw[j]*Constants::R_J_kmol*zone_temperature[j];
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
//																								   //
//					ODE-DAE SYSTEM EQUATIONS - NON ISOTHERMAL REACTOR							   //
//																								   //
/////////////////////////////////////////////////////////////////////////////////////////////////////
void OpenSMOKE_ICEM_MultiZone::ODESystem_ICEM_MultiZone(BzzVector &y, double t, BzzVector &dy)
{
	ErrorMessage("Not available");
}

void OpenSMOKE_ICEM_MultiZone::DAESystem_ICEM_MultiZone(BzzVector &y, double t, BzzVector &dy)
{
	int k, j;

	// 1. Recovering unknowns
	k=1;
	for(j=1;j<=N;j++)
	{
		for(int i=1;i<=NC;i++)
			zone_omega[j][i]	= y[k++];
		zone_temperature[j]		= y[k++];
		zone_g[j]				= y[k++];
		zone_pressure[j]		= y[k++];
		zone_heat_release[j]	= y[k++];
	}
	

	// 2. Recover volumes
	zone_volume[1] = zone_g[1]/zone_pressure[1];
	for(j=2;j<=N;j++)
		zone_volume[j] = (zone_g[j]-zone_g[j-1])/zone_pressure[j];

	// 3. Update Properties
	UpdateProperties(DAE_ICEM_MultiZone_Object.jacobianIndex, DAE_ICEM_MultiZone_Object.jacobianVariables, dimBlock, NC+1);
//	UpdateProperties(MINUSONE, DAE_ICEM_MultiZone_Object.jacobianVariables, dimBlock, NC+1);

	// 4. Update Geometries
	teta = rotation_rate*t - start_angle;
	UpdateAngleAndVolume(teta);
	UpdateVolumeDerivatives(teta);

	// 5. Update heat transfer
	UpdateHeatTransfer();

	// 6. Conservation equations
	Equations_Energy();

	Equations_Pressure();

	Equations_Species();

	// 7. Recover differential equations //ALBERTO
	k=1;
	for(j=1;j<=N;j++)
	{
		// Species equations
		for(int i=1;i<=NC;i++)
			dy[k++] = domega_over_dt[j][i];
		
		// Energy equation
		dy[k++] = dtemperature_over_dt[j];

		// Pressure equations
		dy[k++] = g_equation[j];
		dy[k++] = p_equation[j];

		// Heat release equation
		dy[k++] = zone_QReaction[j]*zone_volume[j];
	}

}


void OpenSMOKE_ICEM_MultiZone::ODESystem_Reaction_ICEM_MultiZone(BzzVector &y, double t, BzzVector &dy)
{
	int k;

	// Recovering unknowns
	k=1;
	for(int i=1;i<=NC;i++)
		zone_omega[jGlobal][i]	= y[k++];
	zone_temperature[jGlobal]	= y[k++];

	// Update Properties
	UpdatePropertiesSingleReactor();

	// Species
	for(int i=1;i<=NC;i++)
		domega_over_dt[jGlobal][i] = zone_R[jGlobal][i]/zone_rho[jGlobal];

	// Energy
	dtemperature_over_dt[jGlobal]  =  zone_UReaction[jGlobal]*zone_volume[jGlobal] / (zone_mass[jGlobal]*zone_Cv[jGlobal]);

	// Recover differential equations
	k=1;
	for(int i=1;i<=NC;i++)
		dy[k++] = domega_over_dt[jGlobal][i];	
	dy[k++] = dtemperature_over_dt[jGlobal];
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
//																								   //
//											DICTIONARY											   //
//																								   //
/////////////////////////////////////////////////////////////////////////////////////////////////////

OpenSMOKE_Dictionary_ICEM_MultiZone::OpenSMOKE_Dictionary_ICEM_MultiZone()
{

    SetupBase();
    
	// Compulsory
	Add("#Model",				'C', 'S', "Engine model: CHEMKIN_MULTIZONE || KOMNINOS_MULTIZONE");
    Add("#RotationRate",		'C', 'M', "Rotation rate");
    Add("#CompressionRatio",	'C', 'D', "Compression ratio");
    Add("#ClearanceVolume",		'O', 'M', "Clearance volume");
    Add("#TotalDisplacement",	'O', 'M', "Total displacement");
    Add("#ArmRatio",			'C', 'D', "Ratio between arms");
    Add("#StartAngle",			'C', 'M', "Starting angle");
	Add("#EndAngle",			'C', 'M', "End angle");
	Add("#Diameter",			'C', 'M', "Cylinder diameter");
    Add("#CrevicesVolume",		'C', 'M', "Crevices volume");
	Add("#NumberOfZones",		'C', 'I', "Number of zones");
	Add("#NumberOfCycles",		'C', 'I', "Number of cycles");
	Add("#ExhaustRatio",		'C', 'D', "Ratio between exhaust gases and fresh gases");
	Add("#Relaxation",			'O', 'I', "Relaxation: number of cycles to be averaged");
	Add("#PlugFlow",		    'O', 'M', "Plug flow residence time");
	Add("#PlugFlowInletTemperature", 'O', 'M', "Plug flow inlet temperature");
	Add("#ResetSpecies",        'O', 'V', "List of species to reset");

	// Optional
    Add("#Key",					'O', 'V', "Key species");
	Add("#Cold",				'O', 'N', "Cold simulation (without chemical reactions)");
	Add("#RelTol",				'O', 'D', "Relative tolerance (default 1e2)");
	Add("#AbsTol",				'O', 'D', "Absolute tolerance (default 1e-10)");

	// Optional
	Add("#ReactionRates",		'O', 'V', "Reaction rates");
	Add("#FormationRates",		'O', 'V', "Formation rates");

	// Heat transfer
	Add("#HeatTransfer",			'O', 'S', "Heat transfer model: NONE, WOSCHNI, ASSANIS, HOHENBERG, ANNAND");
	Add("#Twall",					'O', 'M', "Wall temperature");
	Add("#ExchangeArea",			'O', 'S', "Exchange area: LATERAL || TOTAL");
	Add("#WoschniParameters",		'O', 'V', "Woschni correlation parameters: alfa, c11, c12, c21, c22");
	Add("#AssanisParameters",		'O', 'V', "Assanis correlation parameters: alfa, c11, c12, c21, c22");
	Add("#HohenbergParameters",		'O', 'V', "Hohenberg correlation parameters: alfa, b");
	Add("#AnnandParameters",		'O', 'V', "Annand correlation parameters: alfa, Re");

	// Initial conditions
	Add("#Zones",					'O', 'V', "Zone initial conditions");
	Add("#GaussianZones",			'O', 'V', "Gaussian Temperature profile: mean, sigma, tmin, tmax");

	// Output
	Add("#VerboseEnergy",			'O', 'N', "Report on energy");
	Add("#VerboseGeometry",			'O', 'N', "Geometry output");
	Add("#VerboseMasses",			'O', 'N', "Masses output");
	Add("#VerboseZones",			'O', 'V', "Zones output (species are required)");
	
	
	// Optional (not yet available)
	Add("#ROPA",				'O', 'N', "Rate of Production Analysis");
	Add("#VerboseROPA",			'O', 'V', "Rate of Production Analysis (list of species names)");
	Add("#Sensitivity",			'O', 'V', "Sensitivity Analysis");
	Add("#SensitivityOptions",	'O', 'V', "Sensitivity Options");

    Add("#2E_Model",			'O', 'N', "The Soot 2E Model is solved");

	// IEM model
	Add("#MonteCarloMixing",	'O', 'N', "MonteCarlo Mixing is enabled");
	Add("#TauMixing",			'O', 'M', "Micromixing time");
	Add("#Cphi",				'O', 'D', "Micromixing constant (default 2.0)");

	// Optional
	Add("#Qe",						'O', 'M', "Heat flux");
	Add("#A",						'O', 'M', "Exchange area");
	Add("#U",						'O', 'M', "Heat thermal exchange coefficient");
	Add("#PostProcessing",			'O', 'N', "Post-Processing");

	// User defined options
	Add("#UserDefined_T",			'O', 'S', "The temperature is assigned from external file");
	Add("#UserDefined_Qe",			'O', 'S', "The heat flux is assigned from external file");
    Add("#UserDefined_A",			'O', 'S', "The exchange area is assigned from external file");
	Add("#UserDefined_U",			'O', 'S', "The heat exchange coefficient is assigned from external file");
    Add("#UserDefined_Twall",		'O', 'S', "The wall temperature is assigned from external file");

	// Compulsory
	Compulsory("#ClearanceVolume",	"#TotalDisplacement");

	// Conflicts
	Conflict("#UserDefined_Qe",			"#Qe");
	Conflict("#UserDefined_A",			"#A");
	Conflict("#UserDefined_U",			"#U");
	Conflict("#UserDefined_Twall",		"#Twall");

	// Conflicts
	Conflict("#Qe",						"#U");
	Conflict("#Qe",						"#Twall");
	Conflict("#Qe",						"#UserDefined_U");
	Conflict("#Qe",						"#UserDefined_Twall");
	Conflict("#UserDefined_Qe",			"#U");
	Conflict("#UserDefined_Qe",			"#Twall");
	Conflict("#UserDefined_Qe",			"#UserDefined_U");
	Conflict("#UserDefined_Qe",			"#UserDefined_Twall");
	Conflict("#Zones",					"#GaussianZones");

	Conflict("#WoschniParameters",		"#AssanisParameters");
	Conflict("#WoschniParameters",		"#HohenbergParameters");
	Conflict("#WoschniParameters",		"#AnnandParameters");
	Conflict("#AssanisParameters",		"#HohenbergParameters");
	Conflict("#AssanisParameters",		"#AnnandParameters");
	Conflict("#HohenbergParameters",	"#AnnandParameters");

    Lock();
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
//																								   //
//								ODE-DAE SYSTEM - CLASS DEFINITION								   //
//																								   //
/////////////////////////////////////////////////////////////////////////////////////////////////////

void MyOdeSystem_ICEM_MultiZone::ObjectBzzPrint(void)
{
}

void MyOdeSystem_ICEM_MultiZone::GetSystemFunctions(BzzVector &x, double eta, BzzVector &f)
{
	ptICEM_MultiZone->ODESystem_ICEM_MultiZone(x, eta, f);
}

void MyOdeSystem_ICEM_MultiZone::assignICEM_MultiZone(OpenSMOKE_ICEM_MultiZone *icem_multizone)
{
    ptICEM_MultiZone = icem_multizone;
}

void MyOdeSystem_ICEM_MultiZone::MyODEPrint(BzzVector &y, double eta)
{
	ptICEM_MultiZone->EquationSystemPrint(y, eta);
}


/////////////////////////////////////////////////////////////////////////////////////////////////////
//																								   //
//								ODE-DAE SYSTEM - CLASS DEFINITION								   //
//																								   //
/////////////////////////////////////////////////////////////////////////////////////////////////////

void MyOdeSystem_Reaction_ICEM_MultiZone::ObjectBzzPrint(void)
{
}

void MyOdeSystem_Reaction_ICEM_MultiZone::GetSystemFunctions(BzzVector &x, double eta, BzzVector &f)
{
	ptICEM_MultiZone->ODESystem_Reaction_ICEM_MultiZone(x, eta, f);
}

void MyOdeSystem_Reaction_ICEM_MultiZone::assignICEM_MultiZone(OpenSMOKE_ICEM_MultiZone *icem_multizone)
{
    ptICEM_MultiZone = icem_multizone;
}

void MyOdeSystem_Reaction_ICEM_MultiZone::MyODEPrint(BzzVector &y, double eta)
{
	ptICEM_MultiZone->EquationSystemPrint(y, eta);
}


/////////////////////////////////////////////////////////////////////////////////////////////////////
//																								   //
//									DAE SYSTEM - CLASS DEFINITION								   //
//																								   //
/////////////////////////////////////////////////////////////////////////////////////////////////////

void MyDaeSystem_ICEM_MultiZone::ObjectBzzPrint(void)
{
}


void MyDaeSystem_ICEM_MultiZone::GetSystemFunctions(BzzVector &x, double eta, BzzVector &f)
{
	ptICEM_MultiZone->DAESystem_ICEM_MultiZone(x, eta, f);
}


void MyDaeSystem_ICEM_MultiZone::assignICEM_MultiZone(OpenSMOKE_ICEM_MultiZone *icem_multizone)
{
    ptICEM_MultiZone = icem_multizone;
}


void MyDaeSystem_ICEM_MultiZone::MyDAEPrint(BzzVector &y, double eta)
{
	ptICEM_MultiZone->EquationSystemPrint(y, eta);
}


void OpenSMOKE_ICEM_MultiZone::PrepareAdditionalFiles()
{
    LabelAdditionalFiles();
}


void OpenSMOKE_ICEM_MultiZone::LabelAdditionalFiles()
{
	if (iVerboseGeometry == true)
	{
		fGeometry << "Tau[s](1)      "		<< "\t";
        fGeometry << "Angle[dg](2)   "		<< "\t";
        fGeometry << "Vtot[cm3](3)   "		<< "\t";
        fGeometry << "Vcyl[cm3](4)   "		<< "\t";
        fGeometry << "Height[cm](5)  "		<< "\t";
        fGeometry << "At_cyl[cm2](6) "		<< "\t";
        fGeometry << "Al_cyl[cm2](7) "		<< "\t";

        fGeometry << "Vcore[cm3](8)    "	<< "\t";
        fGeometry << "dVcore[cm3/s](9) "	<< "\t";
		fGeometry << "Hecore[cm](10)   "	<< "\t";
		fGeometry << "Aecore[cm2](11)  "	<< "\t";
		
		int count=12;
		for(int j=1;j<=N-1;j++)
		{
			fGeometry << "V_"  << j << "[cm3]("		<< count++ << ")"	<< "\t";
			fGeometry << "dV_" << j << "[cm3/s]("	<< count++ << ")"	<< "\t";
			fGeometry << "He_" << j << "[cm]("		<< count++ << ")"	<< "\t";
			fGeometry << "Ae_" << j << "[cm2]("		<< count++ << ")"	<< "\t";
			fGeometry << "De_" << j << "[cm]("		<< count++ << ")"	<< "\t";
		}

        fGeometry << endl;
        fGeometry << endl;
    }

	if (iVerboseMasses == true)
	{
		fMass << "Tau[s](1)          "		<< "\t";
        fMass << "Angle[dg](2)       "		<< "\t";
        fMass << "mass_total[g](3)   "		<< "\t";
        fMass << "mass_error[%](4)   "		<< "\t";
       
        fMass << "mass_core[g](5)     "	<< "\t";
        fMass << "dmass_core[g/s](6)  "	<< "\t";
        fMass << "m_in_core[g/s](7)   "	<< "\t";
		fMass << "m_out_core[g/s](8)  "	<< "\t";
		
		int count=9;
		for(int j=1;j<=N-1;j++)
		{
			fMass << "mass_"  << j << "[g]("	<< count++ << ")"	<< "\t";
			fMass << "dmass_" << j << "[g/s]("	<< count++ << ")"	<< "\t";
			fMass << "m_in_"  << j << "[g/s]("	<< count++ << ")"	<< "\t";
			fMass << "m_out_" << j << "[g/s]("	<< count++ << ")"	<< "\t";
		}

        fMass << endl;
        fMass << endl;
    }

	if (iVerboseZones == true)
	{
		fZones << "Tau[s](1)     "		<< "\t";
        fZones << "Angle[dg](2)  "		<< "\t";
		fZones << "P[atm](3)     "		<< "\t";

		int fOutputCount = 4;
		int i;
        for (i=1;i<=N;i++)
            fZones << "T_" << i << "[K]("	<< fOutputCount++ << ")\t";
        for (i=1;i<=N;i++)
            fZones << "n_" << i << "[mol](" << fOutputCount++ << ")\t";
        for (i=1;i<=N;i++)
            fZones << "m_" << i << "[g]("	<< fOutputCount++ << ")\t";

		int j;
		for (j=1;j<=index_species_zones.Size();j++)
			for (i=1;i<=N;i++)
				fZones << "n_" << names_species_zones[j-1] << "_" << i << "[mol](" << fOutputCount++ << ")\t";

        for (j=1;j<=index_species_zones.Size();j++)
			for (i=1;i<=N;i++)
				fZones << "m_" << names_species_zones[j-1] << "_" << i << "[g](" << fOutputCount++ << ")\t";

		fZones << endl << endl;
	}

	if (iVerboseEnergy == true)
	{
		fEnergy << "Tau[s](1)      "	<< "\t";
        fEnergy << "Angle[dg](2)   "	<< "\t";
        fEnergy << "U[J](3)        "	<< "\t"; 
        fEnergy << "H[J](4)        "	<< "\t";
		fEnergy << "Ek[J](5)       "	<< "\t";
        fEnergy << "PV[J](6)       "	<< "\t";
        fEnergy << "QR[J/s](7)     "	<< "\t";
        fEnergy << "QE[???](8)     "	<< "\t";
        fEnergy << "Cp[J/kg/K](9)  "	<< "\t";
        fEnergy << "dummy          "	<< "\t";
        fEnergy << "dummy          "	<< "\t";
        fEnergy << "h[W/m2/K](12)  "	<< "\t";
        fEnergy << "Re[-](13)      "	<< "\t";
        fEnergy << "v[-](14)       "	<< "\t";
        fEnergy << "dummy[-](15)   "	<< "\t";
        fEnergy << "A[cm2](16)     "	<< "\t";
        fEnergy << "Ab[cm2](17)    "	<< "\t";
        fEnergy << "Al[cm2](18)    "	<< "\t";
        fEnergy << "H[mm](19)      "	<< "\t";
		
		fEnergy << endl;
		fEnergy << endl;
	}

	if(assignedExhaustRatio == true)
	{
		int fOutputCount;

		fRecycle << setw(5)  << left << "#";
		fRecycle << setw(16) << left << "T[K](2)";
		fRecycle << setw(16) << left << "P[atm](3)";
		fOutputCount = 4;
		for(int i=1;i<=NC;i++)
			PrintTagOnGnuplotLabel(16, fRecycle, mix->names[i] + "_x",	fOutputCount);
		for(int i=1;i<=NC;i++)
			PrintTagOnGnuplotLabel(16, fRecycle, mix->names[i] + "_w",	fOutputCount);
		fRecycle << endl;

		fExhaust << setw(5) << left << "#";
		fExhaust << setw(16) << left << "T[K](2)";
		fExhaust << setw(16) << left << "P[atm](3)";
		fOutputCount = 4;
		for(int i=1;i<=NC;i++)
			PrintTagOnGnuplotLabel(16, fExhaust, mix->names[i] + "_x",	fOutputCount);
		for(int i=1;i<=NC;i++)
			PrintTagOnGnuplotLabel(16, fExhaust, mix->names[i] + "_w",	fOutputCount);
		fExhaust << endl;
	}
}

void OpenSMOKE_ICEM_MultiZone::PrintAdditionalFiles(const double t, const double teta)
{
	if (iVerboseGeometry == true)
	{
		fGeometry << t								<< "\t";
        fGeometry << teta							<< "\t";
        fGeometry << volume_total*1.e6				<< "\t";
        fGeometry << volume_cylinder_total*1.e6		<< "\t";
		fGeometry << height_cylinder*1.e2			<< "\t";
		fGeometry << area_cylinder_total*1.e4		<< "\t";
		fGeometry << area_cylinder_lateral*1.e4		<< "\t";

		fGeometry << zone_volume[N]*1.e6			<< "\t";
		fGeometry << dvolume_over_dt[N]*1.e6		<< "\t";
		fGeometry << height_e[N]*1.e2				<< "\t";
		fGeometry << area_e[N]*1.e4					<< "\t";

		for(int j=1;j<=N-1;j++)
		{
			fGeometry << zone_volume[j]*1.e6			<< "\t";
			fGeometry << dvolume_over_dt[j]*1.e6		<< "\t";
			fGeometry << height_e[j]*1.e2				<< "\t";
			fGeometry << area_e[j]*1.e4					<< "\t";
			fGeometry << diameter_e[j]*1.e2				<< "\t";
		}

		fGeometry << endl;
    }

	if (iVerboseMasses == true)
	{
		fMass << t								<< "\t";
        fMass << teta							<< "\t";
		fMass << mass_total*1.e3				<< "\t";
		fMass << (mass_total_initial-mass_total)/mass_total_initial*100.				<< "\t";

		fMass << zone_mass[N]*1.e3				<< "\t";
		fMass << dmass_over_dt[N]*1.e3			<< "\t";
		fMass << zone_m_in[N]*1.e3				<< "\t";
		fMass << zone_m_in[N+1]*1.e3			<< "\t";

		for(int j=1;j<=N-1;j++)
		{
			fMass << zone_mass[j]*1.e3			<< "\t";
			fMass << dmass_over_dt[j]*1.e3		<< "\t";
			fMass << zone_m_in[j]*1.e3			<< "\t";
			fMass << zone_m_in[j+1]*1.e3		<< "\t";
		}

		fMass << endl;
    }

	if (iVerboseZones == true)
	{
		fZones << t								<< "\t";
		fZones << teta							<< "\t";
		fZones << zone_pressure[N]/101325.		<< "\t";

		int i;
		for (i=1;i<=N;i++)
			fZones << zone_temperature[i] << "\t";
		for (i=1;i<=N;i++)
			fZones << zone_mole[i]*1.e3 << "\t";
		for (i=1;i<=N;i++)
			fZones << zone_mass[i]*1.e3 << "\t";

		int j;
		for (j=1;j<=index_species_zones.Size();j++)
			for (i=1;i<=N;i++)
				fZones << zone_mole[i]*1.e3*zone_x[i][index_species_zones[j]] << "\t";

		for (j=1;j<=index_species_zones.Size();j++)
			for (i=1;i<=N;i++)
				fZones << zone_mass[i]*1.e3*zone_omega[i][index_species_zones[j]] << "\t";

		fZones << endl;
	}

	if (iVerboseEnergy == true)
	{
		double Ek = 0.;		// [J]
		double dummy = 0.;

		fEnergy << t									<< "\t"		// [s]
				<< teta									<< "\t"		// [deg]
				<< internal_energy_total				<< "\t"		// [J]
				<< enthalpy_total						<< "\t"		// [J]
				<< Ek									<< "\t"		// [J]
				<< enthalpy_total-internal_energy_total	<< "\t"		// [J]
				<< QReaction_total						<< "\t"		// [J/s]
				<< Qe*area_cylinder_total		       	<< "\t"		// [???]
				<< Cp		       						<< "\t"		// [J/kg/K]
				<< dummy	       						<< "\t"		// [???]
				<< dummy	      						<< "\t"		// [???]
				<< U									<< "\t"		// [W/m2/K]
				<< Re									<< "\t"		// [-]
				<< heat_transfer_velocity				<< "\t"		// [m/s]		
				<< 0.									<< "\t"		// [-]
				<< area_cylinder_total*1.e4				<< "\t"		// [cm2]
				<< area_cylinder_base*1.e4				<< "\t"		// [cm2]
				<< area_cylinder_lateral*1.e4			<< "\t"		// [cm2]
				<< height_cylinder*1.e3					<< "\t"		// [mm]
				<< endl;
	}
}

void MyNLSSystem_ICEM_MultiZone_Gaussian::ObjectBzzPrint(void)
{
}

void MyNLSSystem_ICEM_MultiZone_Gaussian::GetResiduals(BzzVector &x, BzzVector &f)
{
	ptICEM->NLSSystem_ICEM_Gaussian(x, f);
}

void MyNLSSystem_ICEM_MultiZone_Gaussian::assignICEM(OpenSMOKE_ICEM_MultiZone *icem)
{
	ptICEM = icem;
}

void OpenSMOKE_ICEM_MultiZone::SaveOnBinaryFile(BzzSave &fOutput)
{
	string dummy;
	char name[Constants::NAME_SIZE];

	int nSteps = countGlobalIterations;			
	BzzVector taugrid	 = GetBzzVector(nSteps, 1, Tau_History);
	BzzVector tetagrid   = GetBzzVector(nSteps, 1, Csi_History);
	BzzVector tgrid		 = GetBzzVector(nSteps, 1, T_History);
	BzzVector pgrid		 = GetBzzVector(nSteps, 1, P_History);

	dummy = "ICEM-MULTI";
	strcpy(name, dummy.c_str());
	fOutput.fileSave.write((char*) name, sizeof(name));

	dummy = "V20100417";
	strcpy(name, dummy.c_str());
	fOutput.fileSave.write((char*) name, sizeof(name));

	// Time
	dummy = "TAU";
	strcpy(name, dummy.c_str());
	fOutput.fileSave.write((char*) name, sizeof(name));
	fOutput << taugrid;

	// Space coordinate
	dummy = "TETA";
	strcpy(name, dummy.c_str());
	fOutput.fileSave.write((char*) name, sizeof(name));
	fOutput << tetagrid;

	// Temperature
	dummy = "TEMPERATURE";
	strcpy(name, dummy.c_str());
	fOutput.fileSave.write((char*) name, sizeof(name));
	fOutput << tgrid;

	// Pressure
	dummy = "PRESSURE";
	strcpy(name, dummy.c_str());
	fOutput.fileSave.write((char*) name, sizeof(name));
	fOutput << pgrid;

	// Mixture molecular weight
	BzzVector MWmix(nSteps);
	for (int i=1;i<=nSteps;i++)
	{
		BzzVector aux = mole_History.GetRow(i);
		MWmix[i] = Dot(aux, mix->M);
	}
	dummy = "MW";
	strcpy(name, dummy.c_str());
	fOutput.fileSave.write((char*) name, sizeof(name));
	fOutput << MWmix;

	// Indices of species
	dummy = "INDICES";
	strcpy(name, dummy.c_str());
	fOutput.fileSave.write((char*) name, sizeof(name));
	fOutput << iOutputSpecies;

	// Indices of species
	dummy = "MOLEFRACTIONS";
	strcpy(name, dummy.c_str());
	fOutput.fileSave.write((char*) name, sizeof(name));
	for(int j=1;j<=NC;j++)
	{
		BzzVector aux = GetBzzVector(nSteps, 1, mole_History.GetColumn(j));
		fOutput << aux;
	}
}

void OpenSMOKE_ICEM_MultiZone::SaveOnBinaryFile(const string filename)
{
	cout << "Save on binary file..." << endl;
	BzzSave fOutput;
	fOutput('*', filename);
	PrintHeaderOnBinaryFile(fOutput);

	cout << "  -- mixture details..." << endl;
	mix->SaveOnBinaryFile(fOutput);

	cout << "  -- flame details..." << endl;
	this->SaveOnBinaryFile(fOutput);
	
	// ROPA
	if (iAssignedROPA == true)
	{
		cout << "  -- writing rate of production analysis..." << endl;

		string dummy;
		char name[Constants::NAME_SIZE];

		BzzVector RRvector(mix->NumberOfReactions());
		BzzMatrix RRmatrix(countGlobalIterations, mix->NumberOfReactions());
		for(int j=1;j<=countGlobalIterations;j++)
		{
			BzzVector aux = mole_History.GetRow(j);
			UpdateReactionRates(T_History[j], P_History[j], aux, RRvector);
			RRmatrix.SetRow(j, RRvector);
		}
		
		dummy = "ROPA";
		strcpy(name, dummy.c_str());
		fOutput.fileSave.write((char*) name, sizeof(name));

		dummy = "V20100417";
		strcpy(name, dummy.c_str());
		fOutput.fileSave.write((char*) name, sizeof(name));

		// Indices of species for local analysis
		BzzVectorInt aux_int(NC);
		for (int j=1;j<=NC;j++)	aux_int[j] = j;
		dummy = "INDICES";
		strcpy(name, dummy.c_str());
		fOutput.fileSave.write((char*) name, sizeof(name));
		fOutput << aux_int;

		// Reaction rates
		dummy = "REACTIONRATES";
		strcpy(name, dummy.c_str());
		fOutput.fileSave.write((char*) name, sizeof(name));
		fOutput << RRmatrix;

		// Verbose ROPA
		if (iVerboseAssignedROPA == true)
		{
			BzzVector xgrid = GetBzzVector(countGlobalIterations,1,Tau_History);
			ropa.SetNumberOfPoints(countGlobalIterations);
			ropa.Run(RRmatrix, xgrid);
			ropa.PrintIntegralRateOfProductionAnalyses(outputFolderName + "/ROPA_Integral.out");
		}
	}

	// Sensitivity
	if (iVerboseSensitivity == true)
	{
		cout << "  -- sensitivity analysis..." << endl; // Ciao

		if ( sensitivity_fast->Integration() != OPENSMOKE_SENSITIVITYANALYSIS_FAST_UNSTEADY0D_ACCURATE )
		{
			ErrorMessage("Sensitivity Analysis not yet available");
			//sensitivity_fast->MoveFromFiles(outputSensitivityName, countGlobalIterations, fOutput);
		}
	//	else if (sensitivity_fast->Integration() == OPENSMOKE_SENSITIVITYANALYSIS_FAST_UNSTEADY0D_ACCURATE)
	//	{
	//	}
	}
	
	PrintEndOnBinaryFile(fOutput);
	fOutput.End();

	// Post processor test
//	OpenSMOKE_PostProcessor post_processor;
//	post_processor.ReadFromBinaryFile(filename);
}

void OpenSMOKE_ICEM_MultiZone::UpdateReactionRates(const double TT, const double PP, BzzVector &xx, BzzVector &rr)
{
	BzzVector RR(mix->NumberOfSpecies());
	BzzVector cc(mix->NumberOfSpecies());
	double MWtot = mix->GetMWFromMoleFractions(xx);
	double cTot  = PP  / (Constants::R_J_kmol*TT);
    double rho   = cTot * MWtot;
            cc   = cTot*xx;
  
    mix->ComputeKineticParameters( TT, log(TT), 1./TT, PP);
	mix->ComputeFromConcentrations( TT, cc, cTot, &RR);			// [kmol/m3/s]
	rr = mix->r;
}
