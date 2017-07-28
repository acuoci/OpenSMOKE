/***************************************************************************
 *   Copyright (C) 2007 by Alberto Cuoci								   *
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

#include "idealreactors/flame1d/OpenSMOKE_Flame1D_OscillatingBoundary.h"


void OpenSMOKE_Flame1D_OscillatingBoundary::ErrorMessage(const std::string message)
{
    cout << endl;
    cout << "Class:  OpenSMOKE_Flame1D_OscillatingBoundary"	<< endl;
    cout << "Object: " << name_object		<< endl;
    cout << "Error:  " << message           << endl;
    cout << "Press a key to continue... "   << endl;
    getchar();
    exit(-1);
}

void OpenSMOKE_Flame1D_OscillatingBoundary::WarningMessage(const std::string message)
{
    cout << endl;
    cout << "Class:  OpenSMOKE_Flame1D_OscillatingBoundary"	<< endl;
    cout << "Object: "		<< name_object	<< endl;
    cout << "Warning:  "	<< message      << endl;
    cout << "Press a key to continue... "   << endl;
    getchar();
}

void OpenSMOKE_Flame1D_OscillatingBoundary::setup(double _USteadyFuel, double _USteadyAir, 
								double _rhoSteadyFuel, double _rhoSteadyAir, 
								double _MWSteadyFuel, double _MWSteadyAir,
								BzzVector& _YSteadyFuel, BzzVector& _YSteadyAir,
								double _TSteadyFuel, double _TSteadyAir,
								double _L, double _P, OpenSMOKE_ReactingGas *_mix, const std::string fileName)
{
	const int SIZE = 120;
	char commento[SIZE];
	double semiAmplitude;
	std::string dummy;
	next_target = 0;
	onPrint = 0;

	mix = _mix;

	ifstream fInput;
	openInputFileAndControl(fInput, fileName);
	
	fInput.getline(commento, SIZE);
	fInput.getline(commento, SIZE);
	fInput.getline(commento, SIZE);
	fInput >> dummy;				fInput.getline(commento, SIZE);

             if (dummy == "SIN")				unsteady_boundary_kind = OSCILLATING_BOUNDARY_SIN;
	else if (dummy == "SIN_EXP")				unsteady_boundary_kind = OSCILLATING_BOUNDARY_SIN_EXP;
	else if (dummy == "AIR_VELOCITY")			unsteady_boundary_kind = OSCILLATING_BOUNDARY_AIR_VELOCITY;
	else if (dummy == "SIN_TEMPERATURE")			unsteady_boundary_kind = OSCILLATING_BOUNDARY_AIR_TEMPERATURE;
	else if (dummy == "FUEL_AIR_VELOCITIES")		unsteady_boundary_kind = OSCILLATING_BOUNDARY_FUEL_AIR_VELOCITIES;
	else if (dummy == "FUEL_COMPOSITION")			unsteady_boundary_kind = OSCILLATING_BOUNDARY_FUEL_COMPOSITION;
	else if (dummy == "RELAXATION_FUEL_VELOCITY")		unsteady_boundary_kind = RELAXATION_BOUNDARY_FUEL_VELOCITY;
	else if (dummy == "RELAXATION_FUEL_TEMPERATURE")	unsteady_boundary_kind = RELAXATION_BOUNDARY_FUEL_TEMPERATURE;
	//else if (dummy == "RELAXATION_FUEL_EQRATIO")		unsteady_boundary_kind = RELAXATION_BOUNDARY_FUEL_EQRATIO;	// TODO

	else if (dummy == "DYNAMIC_TEMPERATURE")		unsteady_boundary_kind = DYNAMIC_TEMPERATURE;
	else if (dummy == "DYNAMIC_VELOCITY")			unsteady_boundary_kind = DYNAMIC_VELOCITY;

	else ErrorMessage("Wrong oscillation kind: " + dummy);
	
	if (	unsteady_boundary_kind == DYNAMIC_TEMPERATURE || unsteady_boundary_kind == DYNAMIC_VELOCITY       )
	{
		fInput >> slope_temperature_fuel_;		fInput.getline(commento, SIZE);		// temperature slope in K/s
		fInput >> slope_temperature_oxidizer_;		fInput.getline(commento, SIZE);		// temperature slope in K/s
		fInput >> slope_velocity_fuel_;			fInput.getline(commento, SIZE);		// velocity slope in cm/s/s
		fInput >> slope_velocity_oxidizer_;		fInput.getline(commento, SIZE);		// velocity slope in cm/s/s

		// Setting steady-state values
		P			= _P;						// [Pa]
		L			= _L;						// [m]
		USteadyFuel		= _USteadyFuel;					// [kg/m2/s]
		USteadyAir		= _USteadyAir;					// [kg/m2/s]
		rhoSteadyFuel		= _rhoSteadyFuel;				// [kg/m3]
		rhoSteadyAir		= _rhoSteadyAir;				// [kg/m3]
		MWSteadyFuel		= _MWSteadyFuel;				// [kg/kmol]
		MWSteadyAir		= _MWSteadyAir;					// [kg/kmol]
		YSteadyFuel = _YSteadyFuel;						// mass fractions
		YSteadyAir = _YSteadyAir;						// mass fractions
		
		ChangeDimensions(YSteadyFuel.Size(), &XSteadyFuel);
		ChangeDimensions(YSteadyAir.Size(), &XSteadyAir);
		mix->GetMoleFractionsFromMassFractions(XSteadyFuel, YSteadyFuel);
		mix->GetMoleFractionsFromMassFractions(XSteadyAir, YSteadyAir);
		
		TSteadyFuel		= _TSteadyFuel;					// [K]
		TSteadyAir		= _TSteadyAir;					// [k]
		vSteadyFuel    		=  200. * USteadyFuel / rhoSteadyFuel;		// [cm/s]
		vSteadyAir      	= -200. * USteadyAir  / rhoSteadyAir;		// [cm/s]
		KSteady         	= 2.e-2*vSteadyAir/L*(1.+vSteadyFuel/vSteadyAir*sqrt(rhoSteadyFuel/rhoSteadyAir));
		xSt			= L / ( 1. + _rhoSteadyAir/_rhoSteadyFuel*BzzPow2(vSteadyAir/vSteadyFuel) );
		betaCoefficient 	= _rhoSteadyFuel/_rhoSteadyAir * (L/xSt-1.);

		tDelay = 0.;
		onPrint = 0;
		lastPrint =-1;

		{
			std::cout << std::endl;
			std::cout << "-----------------------------------------------------------" << std::endl;
			std::cout << "                       DYNAMIC FLAME                       " << std::endl;
			std::cout << "-----------------------------------------------------------" << std::endl;
			std::cout << " * Slope of temperature (fuel side) [K/s]: " << slope_temperature_fuel_		<< std::endl; 
			std::cout << " * Slope of temperature (ox. side) [K/s]:  " << slope_temperature_oxidizer_	<< std::endl; 
			std::cout << " * Slope of velocity (fuel side) [cm/s/s]: " << slope_velocity_fuel_		<< std::endl; 
			std::cout << " * Slope of velocity (ox. side) [cm/s/s]:  " << slope_velocity_oxidizer_		<< std::endl; 
			std::cout << "-----------------------------------------------------------" << std::endl;
			std::cout << std::endl;
		}

		// Dummy variables
		frequency = 1.;
		semiAmplitude = 1.;		
		InPhase = 1;		
		slope = 1.;		
		relaxationTimeFuel = 1.;
		T = 1./frequency;

		// Equivalent Strain rate
		Ntimes  	= 1;
		K		= KSteady;
		Kvector.Append(KSteady);
		timeVector.Append(0.);
	}

	else
	{
		fInput >> frequency;				fInput.getline(commento, SIZE);	
		fInput >> semiAmplitude;			fInput.getline(commento, SIZE);
		fInput >> InPhase;				fInput.getline(commento, SIZE);
		fInput >> slope;				fInput.getline(commento, SIZE);	
		fInput >> relaxationTimeFuel;			fInput.getline(commento, SIZE);	

		// Read possible targets
		fInput >> vTargetFuel;				fInput.getline(commento, SIZE);	
		fInput >> TTargetFuel;				fInput.getline(commento, SIZE);
		fInput >> phiTargetFuel;			fInput.getline(commento, SIZE);

		// Read target velocities
		int  n_target_velocities;	
		fInput >> n_target_velocities;
		list_target_fuel_velocities.resize(n_target_velocities+1);
		for(unsigned int j=0;j<n_target_velocities;j++)
			fInput >> list_target_fuel_velocities[j];
		fInput.getline(commento, SIZE);	

		// Read target temperatures
		int  n_target_temperatures;	
		fInput >> n_target_temperatures;
		list_target_fuel_temperatures.resize(n_target_temperatures+1);
		for(unsigned int j=0;j<n_target_temperatures;j++)
			fInput >> list_target_fuel_temperatures[j];
		fInput.getline(commento, SIZE);	

		// Read target equivalence ratios
		int  n_target_equivalence_ratios;	
		fInput >> n_target_equivalence_ratios;
		list_target_fuel_equivalence_ratios.resize(n_target_equivalence_ratios+1);
		for(unsigned int j=0;j<n_target_equivalence_ratios;j++)
			fInput >> list_target_fuel_equivalence_ratios[j];
		fInput.getline(commento, SIZE);	

		
		fInput.close();

	

		// The sinusoidal shape is assumed by default to be 3 (the best and reliable choice)
		iShape = 3;

		// -------------------------------------------------------------------------------------------------
		// Setting steady state values	
		// -------------------------------------------------------------------------------------------------
		P			= _P;									// [Pa]
		L			= _L;									// [m]
		USteadyFuel		= _USteadyFuel;								// [kg/m2/s]
		USteadyAir		= _USteadyAir;								// [kg/m2/s]
		rhoSteadyFuel		= _rhoSteadyFuel;							// [kg/m3]
		rhoSteadyAir		= _rhoSteadyAir;							// [kg/m3]
		MWSteadyFuel		= _MWSteadyFuel;							// [kg/kmol]
		MWSteadyAir		= _MWSteadyAir;								// [kg/kmol]
		YSteadyFuel = _YSteadyFuel;						// mass fractions
		YSteadyAir = _YSteadyAir;						// mass fractions

		ChangeDimensions(YSteadyFuel.Size(), &XSteadyFuel);
		ChangeDimensions(YSteadyAir.Size(), &XSteadyAir);
		mix->GetMoleFractionsFromMassFractions(XSteadyFuel, YSteadyFuel);
		mix->GetMoleFractionsFromMassFractions(XSteadyAir, YSteadyAir);

		TSteadyFuel		= _TSteadyFuel;								// [K]
		TSteadyAir		= _TSteadyAir;								// [k]
		vSteadyFuel    		=  200. * USteadyFuel / rhoSteadyFuel;		// [cm/s]
		vSteadyAir      	= -200. * USteadyAir  / rhoSteadyAir;		// [cm/s]
		KSteady         	= 2.e-2*vSteadyAir/L*(1.+vSteadyFuel/vSteadyAir*sqrt(rhoSteadyFuel/rhoSteadyAir));
		xSt			= L / ( 1. + _rhoSteadyAir/_rhoSteadyFuel*BzzPow2(vSteadyAir/vSteadyFuel) );
		betaCoefficient 	= _rhoSteadyFuel/_rhoSteadyAir * (L/xSt-1.);
		// -------------------------------------------------------------------------------------------------

		if (unsteady_boundary_kind == RELAXATION_BOUNDARY_FUEL_VELOCITY)
		{
			if (vTargetFuel > vSteadyFuel)
				relaxationFuelCoefficient = fabs(relaxationTimeFuel / log(0.01*vTargetFuel/(vTargetFuel-vSteadyFuel)));
			else
				relaxationFuelCoefficient = fabs(relaxationTimeFuel / log(-0.01*vTargetFuel/(vTargetFuel-vSteadyFuel)));
		}

		if (unsteady_boundary_kind == RELAXATION_BOUNDARY_FUEL_TEMPERATURE)
		{
			if (TTargetFuel > TSteadyFuel)
				relaxationFuelCoefficient = fabs(relaxationTimeFuel / log(0.01*TTargetFuel/(TTargetFuel-TSteadyFuel)));
			else
				relaxationFuelCoefficient = fabs(relaxationTimeFuel / log(-0.01*TTargetFuel/(TTargetFuel-TSteadyFuel)));
		}

		if (unsteady_boundary_kind == RELAXATION_BOUNDARY_FUEL_EQRATIO)
		{
			// TODO
			/*
			if (phiTargetFuel > phiSteadyFuel)
				relaxationFuelCoefficient = fabs(relaxationTimeFuel / log(0.01*phiTargetFuel/(phiTargetFuel-phiSteadyFuel)));
			else
				relaxationFuelCoefficient = fabs(relaxationTimeFuel / log(-0.01*phiTargetFuel/(phiTargetFuel-phiSteadyFuel)));
			*/
		}

		// The semiamplitude for air and fuel is not the same
		// The location of stagnation plane is the same
		semiAmplitudeFuel = semiAmplitude;
		semiAmplitudeAir  = sqrt(betaCoefficient*BzzPow2(vSteadyFuel/vSteadyAir)) * 
							(1.+semiAmplitudeFuel) - 1.00;

		T = 1./frequency;

		_2pi = 2.*acos(-1.);

		tDelay = 0.;
		onPrint = 0;
		lastPrint =-1;

		{
			cout << "------------------------------------------------" << endl;
			cout << "                UNSTEADY FLAME                  " << endl;
			cout << "------------------------------------------------" << endl;

			if (unsteady_boundary_kind == OSCILLATING_BOUNDARY_SIN)
				cout << " SINUSOIDAL OSCILLATIONS"	<< endl;
			if (unsteady_boundary_kind == OSCILLATING_BOUNDARY_SIN_EXP)
				cout << " SINUSOIDAL OSCILLATIONS WITH EXP INCREASING"	<< endl;
	
			if (unsteady_boundary_kind == OSCILLATING_BOUNDARY_SIN || unsteady_boundary_kind == OSCILLATING_BOUNDARY_SIN_EXP)
			{
				cout << " Shape:              "	<< iShape				<< endl; 
				cout << " Frequency:          "	<< frequency			<< endl; 
				cout << " SemiAmplitude Fuel: "	<< semiAmplitudeFuel	<< endl; 
				cout << " SemiAmplitude Air:  "	<< semiAmplitudeAir		<< endl; 
				cout << " Period:             "	<< T					<< endl;
				cout << " In Phase?           " << InPhase               << endl;
			}

			else if (unsteady_boundary_kind == OSCILLATING_BOUNDARY_AIR_VELOCITY)
			{
				cout << " AIR VELOCITY INCREASING (XST not constant)"	<< endl;
				cout << " Slope: " << slope << " kg/m2/s" << endl;
			}

			else if (unsteady_boundary_kind == OSCILLATING_BOUNDARY_FUEL_AIR_VELOCITIES)
			{
				cout << " FUEL-AIR VELOCITIES INCREASING (XST constant)"	<< endl;
				cout << " Slope: " << slope << " kg/m2/s" <<  endl;
			}

			else if (unsteady_boundary_kind == OSCILLATING_BOUNDARY_AIR_TEMPERATURE)
			{
				cout << " AIR TEMPERATURE INCREASING (Xst Constant)"	<< endl;
				cout << " Slope: " << slope << " K/s" <<  endl;
			}

			else if (unsteady_boundary_kind == OSCILLATING_BOUNDARY_FUEL_COMPOSITION)
			{
				cout << " FUEL COMPOSITION INCREASING (Xst Constant)"	<< endl;
				cout << " phi:  " << frequency		<< endl;
				cout << " alfa: " << semiAmplitudeFuel	<<  endl;
			}

			else if (unsteady_boundary_kind == RELAXATION_BOUNDARY_FUEL_VELOCITY)
			{
				cout << " FUEL VELOCITY RELAXATION"<< endl;
				cout << "  * target velocity:        " << vSteadyFuel			<< " cm/s"	<< endl;
				cout << "  * target velocity:        " << vTargetFuel			<< " cm/s"	<< endl;
				cout << "  * relaxation time:        " << relaxationTimeFuel 		<< " s"		<<  endl;
				cout << "  * relaxation coefficient: " << relaxationFuelCoefficient 	<< " s"		<<  endl;
				cout << "  * output velocities:      ";
				for (unsigned int j=0; j<list_target_fuel_velocities.size()-1; j++)
					cout << list_target_fuel_velocities[j] << " ";
				cout << "cm/s" << endl;
			}

			else if (unsteady_boundary_kind == RELAXATION_BOUNDARY_FUEL_TEMPERATURE)
			{
				cout << " FUEL TEMPERATURE RELAXATION"	<< endl;
				cout << "  * target temperature:     "  << TSteadyFuel			<< " K"	<< endl;
				cout << "  * target temperature:     "  << TTargetFuel			<< " K"	<< endl;
				cout << "  * relaxation time:        "  << relaxationTimeFuel 		<< " s"	<<  endl;
				cout << "  * relaxation coefficient: "  << relaxationFuelCoefficient 	<< " s"	<<  endl;
				cout << "  * output temperatures:    ";
				for (unsigned int j=0; j<list_target_fuel_temperatures.size()-1; j++)
					cout << list_target_fuel_temperatures[j] << " ";
				cout << "K" << endl;
	 
			}

			else if (unsteady_boundary_kind == RELAXATION_BOUNDARY_FUEL_EQRATIO)
			{
				// TODO
				/*
				cout << " FUEL EQ. RATIO RELAXATION"	<< endl;
				cout << "  * target eq. ratio:       " << phiSteadyFuel			<< "  "	<< endl;
				cout << "  * target eq. ratio:       " << phiTargetFuel			<< "  "	<< endl;
				cout << "  * relaxation time:        " << relaxationTimeFuel 		<< " s"	<<  endl;
				cout << "  * relaxation coefficient: " << relaxationFuelCoefficient 	<< " s"	<<  endl;
				for (unsigned int j=0; j<list_target_fuel_equivalence_ratios.size()-1; j++)
					cout << list_target_fuel_equivalence_ratios[j] << " ";
				cout << endl;
				*/
			}

			cout << endl;
			cout << " Temperatures: " << TSteadyFuel	<< " - " << TSteadyAir		<< " K"		<< endl;	
			cout << " Velocities:   " << vSteadyFuel	<< " - " << vSteadyAir		<< " cm/s"	<< endl;	
			cout << " Densities:    " << rhoSteadyFuel	<< " - " << rhoSteadyAir	<< " kg/m3"	<< endl;	
			cout << " Mol. weights: " << MWSteadyFuel	<< " - " << MWSteadyAir		<< " kg/mol"<< endl;	
			cout << " Distance:     " << L*100.			<< " cm"  << endl;
			cout << " Pressure:     " << P/101325.		<< " atm" << endl;
			cout << " Strain rate:  " << KSteady		<< " Hz" << endl;
			cout << " Stag. plane:  " << xSt*100.		<< " cm" << endl;
			cout << " Au Fuel:      " << semiAmplitudeFuel*100.	<< " %" << endl;
			cout << " Au Air:       " << semiAmplitudeAir*100.	<< " %" << endl;
			cout << endl;
			cout << "------------------------------------------------" << endl;
			cout << endl;
		}

		// Equivalent Strain rate
		Ntimes  	= 1;
		K		= KSteady;
		Kvector.Append(KSteady);
		timeVector.Append(0.);
	}
}

void OpenSMOKE_Flame1D_OscillatingBoundary::update_boundary_conditions(double _time, 	double &UC,    double &UO, 
											double &TFuel, double &TAir,
											double &rhoC,  double &rhoO, 
											double _Tmax,
											BzzVector &WC, BzzVector &WO,
											OpenSMOKE_Flame1D_DataManager& data)
{
	const double pi = acos(-1.);

	time = _time+tDelay;
	double vAir, vFuel;

	// --------------------------------------------------------
	// 1. Oscillazioni sinusoidali
	//    Il piano di ristagno viene mantenuto nella stessa 
	//    posizione; se le correnti erano in partenza bilanciate,
	//    esse rimarranno tali durante tutto il corso delle
	//    oscillazioni. Anche per la (11) si ottiene lo stesso
	//    tipo di risultato
	// --------------------------------------------------------
	if (unsteady_boundary_kind == OSCILLATING_BOUNDARY_SIN)		// SINUSOIDAL CONSTANT
	{
		nPeriods = int(time/T);
		phase    = (time - nPeriods*T)/T * 360.;

		UC = give_oscillating_profile_Fuel(time);
		UO = give_oscillating_profile_Air(time);

		vAir  = -UO*200. / rhoO;		// [cm/s]
		vFuel =  UC*200. / rhoC;		// [cm/s]
	}

	if (unsteady_boundary_kind == OSCILLATING_BOUNDARY_SIN_EXP)	// SINUSOIDAL EXPONENTIAL
	{
		nPeriods = int(time/T);
		phase    = (time - nPeriods*T)/T * 360.;

		double increasing_factor = 1.1;
		UC = give_oscillating_profile_Fuel(time, increasing_factor);
		UO = give_oscillating_profile_Air(time, increasing_factor);

		vAir  = -UO*200. / rhoO;		// [cm/s]
		vFuel =  UC*200. / rhoC;		// [cm/s]
	}

	// --------------------------------------------------------
	// 2. Aumento la velocita' di ingresso solo dell'aria
	//    Le correnti non possono essere mantenute bilanciate
	// --------------------------------------------------------
	if (unsteady_boundary_kind == OSCILLATING_BOUNDARY_AIR_VELOCITY)
	{
		double q = 1.-1.50*pi;

		if (_time<=pi)
			UO = -USteadyAir;
		else if (_time>pi && _time<=1.50*pi)
			UO = -USteadyAir + slope * ( cos(_time-2.0*pi) + 1.00);
		else if (_time>1.50*pi)
			UO = -USteadyAir + slope * (_time + q);
		UO = UO*-1.0;

		vFuel =  vSteadyFuel;			// [cm/s]
		UC    =  vFuel*rhoC/200.;		// [kg/m2/s]
		vAir  = -200.*UO / rhoO;		// [cm/s]
	}

	// -----------------------------------------------------------
	// 3. Aumento lo strain rate mantenendo costante la posizione
	//    del piano di ristagno; se le correnti erano inizialmente
	//    bilanciate, rimarranno tali durante tutta la dinamica
	// -----------------------------------------------------------
	if (unsteady_boundary_kind == OSCILLATING_BOUNDARY_FUEL_AIR_VELOCITIES)
	{
		double q = 1.-1.50*pi;

		double slopeFuel = slope;
		double slopeAir  = slopeFuel*semiAmplitudeAir/semiAmplitudeFuel;

		if (_time>pi)
		{
			if (_time>pi && _time<=1.50*pi)
				UC = USteadyFuel + slope * ( cos(_time-2.0*pi) + 1.00);
			else if (_time>1.50*pi)
				UC = USteadyFuel + slope * (_time + q);
			
			vFuel = 200.*UC / rhoC;					// [cm/s]
			vAir  = sqrt(rhoC/rhoO*(L/xSt-1.)) * vFuel;		// [cm/s]
			UO = -rhoO * vAir / 200.;
		}
		else
		{
			UO = USteadyAir;
			UC = USteadyFuel;
			vFuel = 200.*UC / rhoC;					// [cm/s]
			vAir  = sqrt(rhoC/rhoO*(L/xSt-1.)) * vFuel;		// [cm/s]
		}

	//	double checkXst = L/(1.+_rhoO/_rhoC*vAir*vAir/vFuel/vFuel);
	//	cout << "Check Stagnation: " <<  checkXst  << " " << _rhoO << " " << _rhoC << " " << vAir << " " << vFuel << endl;
	}

	// --------------------------------------------------------
	// 4. Aumento la temperatura di ingresso dell'aria mantenendo
	//    costante il piano di ristagno (TODO)
	// --------------------------------------------------------
	if (unsteady_boundary_kind == OSCILLATING_BOUNDARY_AIR_TEMPERATURE)
	{
		double correction = 1.;

		double q = 1.-1.50*pi;
		
		if (_Tmax <= 1000.)
		{
			if (_time<=pi)
				TAir = TSteadyAir;
			else if (_time>pi && _time<=1.50*pi)
				TAir = TSteadyAir + slope * ( cos(_time-2.0*pi) + 1.00);
			else if (_time>1.50*pi)
				TAir = TSteadyAir + slope * (_time + q);
		}
		
		rhoO = MWSteadyAir*P/8314./TAir;
		
		// Correction of density because of the equation of state
		if (data.eos.type == EosModel::EOS_PR)
		{
			const double Z = data.eos.Z(TAir, P, XSteadyAir);
			rhoO *= Z;
		}

		betaCoefficient = rhoC/rhoO*(L/xSt-1.0);

		vAir  = vSteadyAir;
		vFuel = sqrt(1./betaCoefficient) * vAir;
		
		UC =  rhoC*vFuel/200.;
		UO = -rhoO*vAir/200.;

	}

	// -----------------------------------------------------------
	// 5. Jackson: modificazione della composizione
	// -----------------------------------------------------------
	if (unsteady_boundary_kind == OSCILLATING_BOUNDARY_FUEL_COMPOSITION)
	{
		double q = 1.-1.50*pi;

		double slopeFuel		= slope;
		double phiSteadyFuel	= frequency;
		double alfa				= semiAmplitudeFuel;
		double phi;

		if (_time>pi)
		{
			if (_time>pi && _time<=1.50*pi)
				phi = phiSteadyFuel - slopeFuel * ( cos(_time-2.0*pi) + 1.00);
			else if (_time>1.50*pi)
				phi = phiSteadyFuel - slopeFuel * (_time + q);
			
			double nTot = (1.-alfa)+4.*alfa+2./phi*(1.+79./21.);
			double xCH4 = (1.-alfa)/nTot;
			double xH2 = 4.*alfa/nTot;
			double xO2 = 2./phi/nTot;
			double xN2 = 2./phi*79./21./nTot;

			int iCH4	= mix->recognize_species("CH4");
			int iH2		= mix->recognize_species("H2");
			int iO2		= mix->recognize_species("O2");
			int iN2		= mix->recognize_species("N2");
			double MWCH4 = mix->M[iCH4];
			double MWH2  = mix->M[iH2];
			double MWO2  = mix->M[iO2];
			double MWN2  = mix->M[iN2];
			
			double MWmix = xCH4*MWCH4 + xH2*MWH2 + xO2*MWO2 +xN2*MWN2;
			
			WC = 0.;
			WC[iCH4]	= xCH4/MWmix*MWCH4;
			WC[iH2]		= xH2/MWmix*MWH2;
			WC[iO2]		= xO2/MWmix*MWO2;
			WC[iN2]		= xN2/MWmix*MWN2;
			WO = WC;

			// TODO: add the cubic equation of state
			rhoC		= P/Constants::R_J_kmol/TAir*MWmix;
			rhoO		= P/Constants::R_J_kmol/TAir*MWmix;

			cout << "phi: " << phi	<< "  wH2:"		<< WC[iH2]	<< " SR: " << 
				2.*vAir/L	*(1.+vFuel/vAir*sqrt(rhoC/rhoO)) / 100. << endl;
			
			UC =  rhoC*vSteadyFuel/200.;
			UO = -rhoO*vSteadyAir/200.;
		}
		else
		{
			rhoC	= rhoSteadyFuel;
			rhoO	= rhoSteadyAir;

			UC =  rhoC*vSteadyFuel/200.;
			UO = -rhoO*vSteadyAir/200.;
		}
	}

	if (unsteady_boundary_kind == RELAXATION_BOUNDARY_FUEL_VELOCITY)
	{
		if (vTargetFuel == vSteadyFuel)
			vFuel = vSteadyFuel;
		else
			vFuel = vTargetFuel + (vSteadyFuel - vTargetFuel)*exp(-time/relaxationFuelCoefficient);	// [cm/s]	
		vAir  = -200.*UO / rhoO;									// [cm/s]
		UC    =  rhoC*vFuel/200.;									// [kg/m2/s]
	}

	if (unsteady_boundary_kind == RELAXATION_BOUNDARY_FUEL_TEMPERATURE)
	{
		if (TTargetFuel == TSteadyFuel)
			TFuel = TSteadyFuel;
		else
			TFuel = TTargetFuel + (TSteadyFuel - TTargetFuel)*exp(-time/relaxationFuelCoefficient);	// [K]	

		rhoC = P/Constants::R_J_kmol/TFuel*MWSteadyFuel;					// [kg/m3]
			
		// Correction of density because of the equation of state
		if (data.eos.type == EosModel::EOS_PR)
		{
			const double Z = data.eos.Z(TFuel, P, XSteadyFuel);
			rhoC *= Z;
		}

		vAir  = -200.*UO / rhoO;								// [cm/s]
		UC    =  rhoC*vSteadyFuel/200.;								// [kg/m2/s]
	}

	if (unsteady_boundary_kind == RELAXATION_BOUNDARY_FUEL_EQRATIO)
	{
		// TODO
		/*
		if (phiTargetFuel == phiSteadyFuel)
			phiFuel = phiSteadyFuel;
		else
			phiFuel = phiTargetFuel + (phiSteadyFuel - phiTargetFuel)*exp(-time/relaxationFuelCoefficient);	
	
		vAir  = -200.*UO / rhoO;								// [cm/s]
		UC    =  rhoC*vSteadyFuel/200.;								// [kg/m2/s]
		*/
	}

	if (unsteady_boundary_kind == DYNAMIC_TEMPERATURE)
	{
		if (slope_temperature_fuel_ != 0.)
		{
			TFuel = TSteadyFuel + slope_temperature_fuel_*time;	// [K]

			rhoC  = P/Constants::R_J_kmol/TFuel*MWSteadyFuel;	// [kg/m3]
			// Correction of density because of the equation of state
			if (data.eos.type == EosModel::EOS_PR)
			{
				const double Z = data.eos.Z(TFuel, P, XSteadyFuel);
				rhoC *= Z;
			}
			UC    = rhoC*vSteadyFuel/200.;				// [kg/m2/s]
		}

		if (slope_temperature_oxidizer_ != 0.)
		{
			TAir  = TSteadyAir + slope_temperature_oxidizer_*time;	// [K]
			
			rhoO  = P/Constants::R_J_kmol/TAir*MWSteadyAir;		// [kg/m3]
			// Correction of density because of the equation of state
			if (data.eos.type == EosModel::EOS_PR)
			{
				const double Z = data.eos.Z(TAir, P, XSteadyAir);
				rhoO *= Z;
			}

			UO    = -rhoO*vSteadyAir/200.;				// [kg/m2/s]
		}

		vFuel =  200.*UC / rhoC;
		vAir  = -200.*UO / rhoO;
	}

	if (unsteady_boundary_kind == DYNAMIC_VELOCITY)
	{
		if (slope_velocity_fuel_ != 0.)
		{
			vFuel = vSteadyFuel + slope_velocity_fuel_*time;	// [cm/s]
			UC    = rhoC*vFuel/200.;				// [kg/m2/s]
			vAir  = -200.*UO / rhoO;
		}

		if (slope_velocity_oxidizer_ != 0.)
		{
			vAir  = vSteadyAir + slope_velocity_oxidizer_*time;	// [K]
			UO    = -rhoO*vAir/200.;				// [kg/m2/s]
			vFuel =  200.*UC / rhoC;
		}
	}

	// Strain rate (instantaneous)
	// --------------------------------------------------------
	K		= 2.*vAir/L	*(1.+vFuel/vAir*sqrt(rhoC/rhoO)) / 100.;
	KSeshadri	= 2.*vFuel/L*(1.-vAir/vFuel*sqrt(rhoO/rhoC)) / 100.;

}

void OpenSMOKE_Flame1D_OscillatingBoundary::update_time_target(const double time)
{
	onPrint = 0;

	if (unsteady_boundary_kind == RELAXATION_BOUNDARY_FUEL_VELOCITY)
	{
		double vFuel;

		if (vTargetFuel == vSteadyFuel)
			vFuel = vSteadyFuel;
		else
			vFuel = vTargetFuel + (vSteadyFuel - vTargetFuel)*exp(-time/relaxationFuelCoefficient);	// [cm/s]	
		
		if (vTargetFuel > vSteadyFuel)
		{
			list_target_fuel_velocities[list_target_fuel_velocities.size()-1] = 1.e32;
			if (vFuel >= list_target_fuel_velocities[next_target])
			{
				onPrint = 1;
				next_target++;
			}
		}
		else
		{
			list_target_fuel_velocities[list_target_fuel_velocities.size()-1] = -1.e32;
			if (vFuel <= list_target_fuel_velocities[next_target])
			{
				onPrint = 1;
				next_target++;
			}
		}						
	}

	if (unsteady_boundary_kind == RELAXATION_BOUNDARY_FUEL_TEMPERATURE)
	{
		double TFuel;

		if (TTargetFuel == TSteadyFuel)
			TFuel = TSteadyFuel;
		else
			TFuel = TTargetFuel + (TSteadyFuel - TTargetFuel)*exp(-time/relaxationFuelCoefficient);	// [K]	

		if (TTargetFuel > TSteadyFuel)
		{
			list_target_fuel_temperatures[list_target_fuel_temperatures.size()-1] = 1.e32;
			if (TFuel >= list_target_fuel_temperatures[next_target])
			{
				onPrint = 1;
				next_target++;
			}
		}
		else
		{
			list_target_fuel_temperatures[list_target_fuel_temperatures.size()-1] = -1.e32;
			if (TFuel <= list_target_fuel_temperatures[next_target])
			{
				onPrint = 1;
				next_target++;
			}
		}
	}
}

double OpenSMOKE_Flame1D_OscillatingBoundary::give_oscillating_profile_Fuel(double time)
{
	//if (iShape==3)
	{
		double alfa = -log(1.e-5) / T;
		return USteadyFuel*(1.+ (1.-exp(-alfa*time))*semiAmplitudeFuel*sin(_2pi*frequency*(time+tDelay)));
	}
}

double OpenSMOKE_Flame1D_OscillatingBoundary::give_oscillating_profile_Air(double time)
{
	//if (iShape==3)
	{
		double alfa = -log(1.e-5) / T;
		return USteadyAir*(1. + InPhase*(1.-exp(-alfa*time))*semiAmplitudeAir*sin(_2pi*frequency*(time+tDelay)));
	}
}

double OpenSMOKE_Flame1D_OscillatingBoundary::give_oscillating_profile_Fuel(double time, double incresing_factor)
{
	{
		double alfa = -log(1.e-5) / T;
		return USteadyFuel*(1.+ (1.-exp(-alfa*time))*semiAmplitudeFuel*sin(_2pi*frequency*(time+tDelay)) * 
			                    incresing_factor*exp(log(incresing_factor)/T*time));
	}
}

double OpenSMOKE_Flame1D_OscillatingBoundary::give_oscillating_profile_Air(double time, double incresing_factor)
{
	{
		double alfa = -log(1.e-5) / T;
		return USteadyAir*(1. + InPhase*(1.-exp(-alfa*time))*semiAmplitudeAir*sin(_2pi*frequency*(time+tDelay)) * 
			                    incresing_factor*exp(log(incresing_factor)/T*time));
	}
}


// ----------------------------------------------------------------------------------
// TODO
// ----------------------------------------------------------------------------------
void OpenSMOKE_Flame1D_OscillatingBoundary::update_K_equivalent(double _time)
{
	if (_time>0.)
	{
		Ntimes++;
		Kvector.Append(K);
		timeVector.Append(_time);
		Kequivalent = integral_b();
		cout << Ntimes << endl;;
	}
}

double OpenSMOKE_Flame1D_OscillatingBoundary::integral_a(int j)
{
	double sum = 0.;
	for(int i=j;i<=Ntimes-1;i++)
		sum += (timeVector[i+1]-timeVector[i]) * 0.50 * (Kvector[i+1]+Kvector[i]);
	return sum;
}

double OpenSMOKE_Flame1D_OscillatingBoundary::integral_b()
{
	double sum= 0.;
	for(int i=1;i<=Ntimes-1;i++)
	{
		double a = exp(-2.*integral_a(i));
		double b = exp(-2.*integral_a(i+1));
		sum += (timeVector[i+1]-timeVector[i]) * 0.50 * (a+b);
	}
	return 1./(2.*sum);
}
