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
#include "idealreactors/shocktube/OpenSMOKE_ShockTube_InitialConditions.h"
#include "idealreactors/shocktube/OpenSMOKE_ShockTube.h"
#include "basic/OpenSMOKE_Conversions.h"
//#include "addons/OpenSMOKE_PostProcessor.h"

OpenSMOKE_ShockTube	*ptShockTube;

void ODEPrintExternal_ShockTube(BzzVector &y, double csi)
{
    ptShockTube->ODEPrint(y, csi);
}

const double    OpenSMOKE_ShockTube::EPSILON			=  1.e-6;
const double    OpenSMOKE_ShockTube::_1_PLUS_EPSILON	=  1.+EPSILON;
const double    OpenSMOKE_ShockTube::_1_MINUS_EPSILON   =  1.-EPSILON;

OpenSMOKE_ShockTube::OpenSMOKE_ShockTube()
{
    ptShockTube		=  this;
	class_name		= "OpenSMOKE_ShockTube";
	out_name		= "ShockTube.out";
	UReflectedShock	= 0.;
	cout.setf(ios::scientific);
	
	Setup();

	assignedEnergy				= true;		// 
	assignedInletFlows			= false;	//
	assignedKindOfReactor		= false;	// 
	assignedShockVelocity		= false;	//
	assignedShockTemperature	= false;	//	
	assignedShockPressure		= false;	//	
	assignedShockDensity		= false;	//	
	assignedShockComposition	= false;	//	
	assignedGeometry			= false;	// Diameter or section
	
	iShockEquivalenceRatio		= false;
	iBoundaryLayerCorrection	= false;		// Viscosity correction
	iIncidentShock				= true;			// Incident or Reflected Shock
	iInitialCondition			= SHOCK_NONE;	// Kind of initial condition
}

void OpenSMOKE_ShockTube::AssignIncidentShock()
{
	iIncidentShock = true;	
    assignedKindOfReactor = true;
}

void OpenSMOKE_ShockTube::AssignReflectedShock()
{
	iIncidentShock = false;	
    assignedKindOfReactor = true;
}

void OpenSMOKE_ShockTube::AssignShockVelocity(const std::string units, const double value)
{
	UShock = OpenSMOKE_Conversions::conversion_velocity(value, units);	
    assignedShockVelocity = true;
}

void OpenSMOKE_ShockTube::SetReflectedShockVelocity(const std::string units, const double value)
{
	UReflectedShock = OpenSMOKE_Conversions::conversion_velocity(value, units);	
}

void OpenSMOKE_ShockTube::AssignShockTemperature(const std::string units, const double value, const kindInitialCondition _iInitialCondition)
{
	inletStream->AssignTemperature(value, units);	
    assignedShockTemperature	= true;
	iInitialCondition			= _iInitialCondition;
}

void OpenSMOKE_ShockTube::AssignShockPressure(const std::string units, const double value, const kindInitialCondition _iInitialCondition)
{
	inletStream->AssignPressure(value, units);	
    assignedShockPressure	= true;
	iInitialCondition		= _iInitialCondition;
}

void OpenSMOKE_ShockTube::AssignShockDensity(const std::string units, const double value, const kindInitialCondition _iInitialCondition)
{
	inletStream->AssignDensity(value, units);	
    assignedShockDensity	= true;
	iInitialCondition		= _iInitialCondition;
}

void OpenSMOKE_ShockTube::AssignShockMoleFractions(const vector<string> _names, const vector<double> _values)
{
	inletStream->AssignMoleFractions(_names, _values);
	assignedShockComposition = true;
}

void OpenSMOKE_ShockTube::AssignShockMassFractions(const vector<string> _names, const vector<double> _values)
{
	inletStream->AssignMassFractions(_names, _values);
	assignedShockComposition = true;
}

void OpenSMOKE_ShockTube::SetBoundaryLayerCorrection()
{
	iBoundaryLayerCorrection = true;
}

void OpenSMOKE_ShockTube::UnsetBoundaryLayerCorrection()
{
	iBoundaryLayerCorrection = false;
}

void OpenSMOKE_ShockTube::AssignEnd(const std::string units, const double value)
{
	TauTotal = OpenSMOKE_Conversions::conversion_time(value, units);	
    assignedEnd = true;
}

void OpenSMOKE_ShockTube::AssignDiameter(const std::string units, const double value)
{
    D       = OpenSMOKE_Conversions::conversion_length(value, units);
    Area    = Constants::pi * D*D / 4.;
	Area0	= Area;
	
    assignedGeometry = true;
}

void OpenSMOKE_ShockTube::AssignArea(const std::string units, const double value)
{
    Area    = OpenSMOKE_Conversions::conversion_area(value, units);
    D       = sqrt(4./Constants::pi*Area);
	Area0	= Area;

    assignedGeometry = true; 
}

void OpenSMOKE_ShockTube::Lock()
{
    if (assignedKineticScheme == false)
        ErrorMessage("The kinetic scheme was not defined!");
    if (assignedGlobalKineticScheme == false)
        ErrorMessage("The global kinetic scheme was not defined!");
    if (assignedEnd == false)
    	ErrorMessage("The reactor length was not defined!");
	if (assignedGeometry == false)
        ErrorMessage("The reactor geometry was not defined!");
    if (assignedInletFlows == false)
        ErrorMessage("The inlet flow was not defined!");
	if (assignedKindOfReactor == false)
        ErrorMessage("The kind of reactor was not defined!");
	if (assignedSoot2EModel == false)
        ErrorMessage("The Soot 2E Model was not defined!");

	if (iIncidentShock == true && assignedShockVelocity == false)
        ErrorMessage("The shock velocity must be assigned");

	if (iInitialCondition == SHOCK_BEFORE && assignedShockVelocity == false)
        ErrorMessage("The shock velocity must be assigned");

	if (iInitialCondition == SHOCK_AFTER && assignedShockVelocity == false)
        ErrorMessage("The shock velocity must be assigned");

	if (UReflectedShock > 0. && assignedShockVelocity == false)
        ErrorMessage("The shock velocity must be assigned when the reflected shock velocity is given");

	if (iTwoEquationModel == true && mix->iSootMode == false)
		ErrorMessage("This kinetic scheme was interpreted without the SootMode option and cannot be used together the 2E model.\nMore details on the OpenSMOKE's User Guide.");

	inletStream->AssignMassFlowRate(1.0, "kg/s");
	inletStream->lock();

	Initialize();
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
//																								   //
//									OUTPUT FUNCTIONS											   //
//																								   //
/////////////////////////////////////////////////////////////////////////////////////////////////////

void OpenSMOKE_ShockTube::VideoFinalResult()
{
    int i;

    // Vectors
    BzzVector conversion(NC);

    // Print General Information about ShockTube
    VideoGeneralInfo();

    // Calculate Conversion of Inlet Species
    for (i=1;i<=NC;i++)
        if (shockStream.omega[i]!=0.)
            conversion[i] = (1. - omega[i]/shockStream.omega[i])*100.;

    // Final Composition
    cout << "#\tName\tx\t\tomega" << endl;
    cout << "---------------------------------------------------------------------" << endl;
    for (i=1;i<=NC;i++)
        if (x[i]!=0.)
            cout << i << "\t" << mix->names[i] << "\t" << x[i] << "\t" << omega[i] << endl;
    cout << endl;

    // Conversion
    cout << "#\tName\tomegaInlet\tomegaOutlet\tConversion(%)" << endl;
    cout << "---------------------------------------------------------------------" << endl;
    for (i=1;i<=NC;i++)
        if (conversion[i]>0.)
            cout << i << "\t" << mix->names[i] << "\t"
            << inletStream->omega[i] << "\t" << omega[i] << "\t" << conversion[i] << endl;
    cout << endl;

    cout << "xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx" << endl;
    cout << "xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx" << endl;
    cout << endl;
}

void OpenSMOKE_ShockTube::VideoGeneralInfo()
{
    cout.setf(ios::scientific);
    cout << endl;

    cout << "xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx" << endl;
    cout << "xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx" << endl;
    cout << endl;

    cout << "---------------------------------------------------------------------" << endl;
    cout << "Shock Tube Summary" << endl;
    cout << "---------------------------------------------------------------------" << endl;
    cout << "Length:            "		<< Csi								<< " [m]"		<< endl;
    cout << "Residence time:    "		<< Tau								<< " [s]"		<< endl;
	cout << "Mass flow rate:    "		<< shockStream.massFlowRate			<< " [kg/s]"	<< endl;
    cout << "Mole flow rate:    "		<< shockStream.moleFlowRate			<< " [kmol/s]"	<< endl;
    cout << "Volume flow rate:  "		<< shockStream.volumetricFlowRate	<< " [m3/s]"	<< endl;
	cout << "Pressure:          "		<< P								<< " [Pa]"		<< endl;
    cout << "Density:           "		<< rho								<< " [kg/m3]"	<< endl;
    cout << "Molecular weight:  "		<< MWtot							<< " [kg/kmol]" << endl;
    cout << "Velocity:          "		<< v								<< " [m/s]"		<< endl;
    cout << "Area:              "		<< Area								<< " [m2]"		<< endl;
    cout << "Temperature:       "		<< T								<< " [K]"		<< endl;
    cout <<  endl;
}

void OpenSMOKE_ShockTube::VideoSummary()
{
	VideoGeneralInfo();

    cout << setw(6) << left << "#" << setw(20) << left << "name" << setw(16) << left << "x" << setw(16) << left << "omega" << endl;
    cout << "---------------------------------------------------------------------" << endl;
    for (int i=1;i<=NC;i++)
        if (inletStream->omega[i]!=0.)
            cout << setw(6)  << left << i 
				 << setw(20) << left << mix->names[i] 
				 << setw(16) << left << shockStream.x[i] 
				 << setw(16) << left << shockStream.omega[i]
				 << endl;
    cout << "---------------------------------------------------------------------" << endl;
    cout << endl << endl;
}

void OpenSMOKE_ShockTube::SummaryOnFile()
{
	int i;

	ofstream fSummary;
	std::string file_name = outputFolderName + "/Summary.out";
	openOutputFileAndControl(fSummary, file_name);
    fSummary.setf(ios::scientific);

    fSummary << "                  	Inlet        	  Outlet"										<< endl;
	fSummary << "----------------------------------------------------------------------------"		<< endl;
    fSummary << setw(20) << left << "Time[s]"			<< setw(16) << left << ZERO							<< setw(16) << left << Tau								<< endl;
    fSummary << setw(20) << left << "Length[m]"			<< setw(16) << left << ZERO							<< setw(16) << left << Csi								<< endl;
    fSummary << setw(20) << left << "Temperature[K]"	<< setw(16) << left << shockStream.T				<< setw(16) << left << T								<< endl;
    fSummary << setw(20) << left << "Pressure[Pa]"		<< setw(16) << left << shockStream.P				<< setw(16) << left << P								<< endl;
    fSummary << setw(20) << left << "Mass_flow[kg/s]"	<< setw(16) << left << shockStream.massFlowRate		<< setw(16) << left << shockStream.massFlowRate			<< endl;
    fSummary << setw(20) << left << "Mole_flow[kmol/s]"	<< setw(16) << left << shockStream.moleFlowRate		<< setw(16) << left << shockStream.massFlowRate/MWtot	<< endl;
    fSummary << setw(20) << left << "Density[kg/m3]"	<< setw(16) << left << shockStream.rho				<< setw(16) << left << rho								<< endl;
    fSummary << setw(20) << left << "MW[kg/kmol]"		<< setw(16) << left << shockStream.MW				<< setw(16) << left << MWtot							<< endl;
	fSummary << endl;
	
	
	fSummary << "Mole_fractions[-]" << endl;
	for(i=1;i<=NC;i++)
		fSummary << setw(20) << left << mix->names[i] << setw(16) << left << shockStream.x[i] << setw(16) << left << x[i]	<< endl;
	fSummary << endl;

	fSummary << endl;
	fSummary << "Mass_fractions[-]" << endl;
	for(i=1;i<=NC;i++)
		fSummary << setw(20) << left << mix->names[i] << setw(16) << left << shockStream.omega[i] << setw(16) << left << omega[i]	<< endl;
	fSummary << endl;
	
	fSummary << endl;
	fSummary << "Concentrations[kmol/m3]" << endl;
	for(i=1;i<=NC;i++)
		fSummary << setw(20) << left << mix->names[i] << setw(16) << left << shockStream.c[i] << setw(16) << left << c[i]	<< endl;
	fSummary << endl;

	fSummary << endl;
	fSummary << "Conversions[-]" << endl;
	for(i=1;i<=NC;i++)
		if (inletStream->omega[i] > 0.)
			fSummary << setw(20) << left << mix->names[i] << setw(16) << left << 1.-omega[i]/shockStream.omega[i]	<< endl;
	fSummary << endl;

    fSummary.close();
}


void OpenSMOKE_ShockTube::EnergyAnalysis(OpenSMOKE_GasStream &outletStream)
{
	double H_mass_inlet  = shockStream.massSpecificEnthalpy;			// [J/kg]
	double P_mass_inlet  = shockStream.massSpecificPressureEnergy;		// [J/kg]
	double U_mass_inlet  = shockStream.massSpecificInternalEnergy;		// [J/kg]
	double K_mass_inlet  = 0.50*BzzPow2(shockStream.velocity);			// [J/kg]
	double E_mass_inlet  = H_mass_inlet + K_mass_inlet;					// [J/kg]
	
	double H_mass_outlet  = outletStream.massSpecificEnthalpy;			// [J/kg]
	double P_mass_outlet  = outletStream.massSpecificPressureEnergy;	// [J/kg]
	double U_mass_outlet  = outletStream.massSpecificInternalEnergy;	// [J/kg]
	double K_mass_outlet  = 0.50*v*v;									// [J/kg]
	double E_mass_outlet  = H_mass_outlet + K_mass_outlet;				// [J/kg]

	cout << endl;
	cout << "---------------------------------------------------------------------" << endl;
	cout << "                         ENERGY ANALYSIS                             " << endl;
	cout << "---------------------------------------------------------------------" << endl;
	cout << "Inlet  - Enthalpy:       "	<< H_mass_inlet	  << " [J/kg]" << endl;
	cout << "Inlet  - Kinetic energy: " << K_mass_inlet	  << " [J/kg]" << endl;
	cout << "Inlet  - Total energy:   " << E_mass_inlet   << " [J/kg]" << endl;
	cout << "Outlet - Enthalpy:       "	<< H_mass_outlet  << " [J/kg]" << endl;
	cout << "Outlet - Kinetic energy: " << K_mass_outlet  << " [J/kg]" << endl;
	cout << "Outlet - Total energy:   " << E_mass_outlet  << " [J/kg]" << endl;
	cout << endl;

	cout << "Inlet power:     "			<< E_mass_inlet						<< " [J/kg]" << endl;
	cout << "Outlet power:    "			<< -E_mass_outlet					<< " [J/kg]" << endl;
	cout << "Exchanged power: "			<< -(E_mass_outlet-E_mass_inlet)	<< " [J/kg]" << endl;
}

void OpenSMOKE_ShockTube::ODEPrint(BzzVector &y, double eta)
{
	double dummy = ZERO;
	Tau = eta;

    if (iVerbose == true)
    {
        if ( countIterations%(300*nVideoSteps) == 0)
        {
            cout    << endl
                    << "#"			<< "\t"
                    << "Tau[s]"		<< "\t\t"
                    << "x[m]"		<< "\t\t"
                    << "T[K]"       << "\t\t"
                    << "v[m/s]"     << "\t\t"
					<< "P[atm]"		<< "\t\t";
                    
			cout	<< endl;
        }

        if (countVideoSteps == nVideoSteps)
        {
            cout	<< countIterations	<< "\t"
                    << Tau  			<< "\t"
                    << Csi				<< "\t"
                    << T				<< "\t"
                    << v				<< "\t"
					<< P/101325.		<< "\t";

            // Soot Two Equation Model
            if (iTwoEquationModel == true)
            {
                cout 	<< soot2EModel->m0	<< "\t"
                        << soot2EModel->fv	<< "\t";
            }

            cout	<< endl;

            countVideoSteps = 0;
        }

        if (countFileSteps == nFileSteps)
        {
            int i;

            fOutput << setw(20) << left << Tau
                    << setw(20) << left << Csi
                    << setw(20) << left << T
                    << setw(20) << left << tL
                    << setw(20) << left << P/101325.
                    << setw(20) << left << v
                    << setw(20) << left << QReaction
                    << setw(20) << left << dummy
                    << setw(20) << left << rho
                    << setw(20) << left << MWtot
                    << setw(20) << left << D
                    << setw(20) << left << Area
					<< setw(20) << left << dummy
					<< setw(20) << left << dummy;

			if(key_species_names.size() == 1)
				fOutput << setw(20) << left << 1.-omega[key_species_index[0]]/inletStream->omega[key_species_index[0]];
			else
	            fOutput << setw(20) << left << dummy;

			// Only user defined species
            for (i=1;i<=iOutputSpecies.Size();i++)
                fOutput << setw(20) << left << x[iOutputSpecies[i]];
            for (i=1;i<=iOutputSpecies.Size();i++)
                fOutput << setw(20) << left << omega[iOutputSpecies[i]];
            fOutput << endl;

            // Soot Two Equation Model
            if (iTwoEquationModel == true)
			{
				fSoot2E	<< setw(20) << left << Tau
						<< setw(20) << left << Csi
						<< setw(20) << left << T;

				soot2EModel->write_on_file(fSoot2E, y[indexTwoEquations], y[indexTwoEquations+1]);
			}

/*            if (mix->soot_manager.iSoot == true)
            {
                // Soot Distribution data
                mix->soot_manager.large_calculateAll(omega, x, c, rho);
                mix->soot_manager.small_calculateAll(omega, x, c, rho);

                fSoot << setw(20) << left << Tau
					  << setw(20) << left << Csi
					  << setw(20) << left << T;
                mix->soot_manager.print_main_data_on_file(fSoot);
                fSoot << endl;
            }
*/
/*            if (mix->pah_manager.iPAH == true)
            {
                // PAH Distribution data
                mix->pah_manager.pah_calculateAll(omega, x, c, rho);

                fPAH << setw(20) << left << Tau
					 << setw(20) << left << Csi
					 << setw(20) << left << T;
                mix->pah_manager.print_main_data_on_file(fPAH);
                fPAH << endl;
            }
*/
			if (iKeySpecies == true && key_species_names.size()>1)
            {
                fConversions << setw(20) << left << Tau
							 << setw(20) << left << Csi
							 << setw(20) << left << T;

				for(i=0;i<int(key_species_names.size());i++)
	                fConversions << setw(20) << left << 1.-omega[key_species_index[i]]/inletStream->omega[key_species_index[i]];
	        
                fConversions << endl;
            }

			if (iVerboseEnergy == true)
			{
				BzzVector h_mass(NC);
				mix->GetMixAveragedEnthalpy_Mass(h_mass, T);
	
				double H_mass = Dot(omega, h_mass);						// [J/kg]
				double Ek_mass = 0.50*v*v;
				double Ep_mass = P/rho;

				fEnergy << setw(20) << left << Tau
						<< setw(20) << left << Csi
						<< setw(20) << left << H_mass+Ep_mass
						<< setw(20) << left << H_mass
						<< setw(20) << left << Ek_mass
						<< setw(20) << left << Ep_mass
						<< setw(20) << left << QReaction
						<< setw(20) << left << dummy
						<< setw(20) << left << Cp
						<< endl;
			}

			if (iVerboseReactionRates == true)
			{
				fReactionRates	<< setw(20) << left << Tau
								<< setw(20) << left << Csi
								<< setw(20) << left << T;
				
				if (iGlobalKinetics == false)
					for (i=1;i<=index_reaction_rates.Size();i++)
						fReactionRates	<< setw(20) << left << mix->r[index_reaction_rates[i]];
				else
					for (i=1;i<=global->nReactions;i++)
						fReactionRates	<< setw(20) << left << global->r[i];

				fReactionRates	<< endl;
			}

			if (iVerboseFormationRates == true)
			{
				fFormationRates	<< setw(20) << left << Tau
								<< setw(20) << left << Csi
								<< setw(20) << left << T;
			
				for (i=1;i<=index_formation_rates.Size();i++)
					fFormationRates	<< setw(20) << left << R[index_formation_rates[i]]/mix->M[index_formation_rates[i]];
				
				fFormationRates	<< endl;
			}

			if (iAssignedROPA == true)
			{
			}

			if (iVerboseSensitivity == true)
			{
			}

            countFileSteps = 0;
        }

    }

    countVideoSteps++;
    countFileSteps++;
    countIterations++;

	if (iHistory == true)
    {
        countGlobalIterations++;

        if(countGlobalIterations>0)
        {
			int index = countGlobalIterations;

			if (countGlobalIterations>MAX_TIME_STEPS)
				index = MAX_TIME_STEPS;    
			
			Tau_History[index] = Tau;
			Csi_History[index] = Csi;
			T_History[index]   = T;
			mass_History.SetRow(index, omega);
			mole_History.SetRow(index, x);			
        }
    }
}

void OpenSMOKE_ShockTube::LabelODEFile()
{
    int i;

    if (iVerbose == true)
    {
        {
			int fOutputCount = 1;
			PrintTagOnGnuplotLabel(20, fOutput, "Tau[s]",		fOutputCount);
			PrintTagOnGnuplotLabel(20, fOutput, "x[m]",			fOutputCount);
			PrintTagOnGnuplotLabel(20, fOutput, "T[K]",			fOutputCount);
			PrintTagOnGnuplotLabel(20, fOutput, "tL[s]",		fOutputCount);
			PrintTagOnGnuplotLabel(20, fOutput, "P[atm]",		fOutputCount);
			PrintTagOnGnuplotLabel(20, fOutput, "v[m/s]",		fOutputCount);
			PrintTagOnGnuplotLabel(20, fOutput, "QR[W/m3]",		fOutputCount);
			PrintTagOnGnuplotLabel(20, fOutput, "dummy",		fOutputCount);
			PrintTagOnGnuplotLabel(20, fOutput, "rho[kg/m3]",	fOutputCount);
			PrintTagOnGnuplotLabel(20, fOutput, "MW",			fOutputCount);
			PrintTagOnGnuplotLabel(20, fOutput, "D[m]",			fOutputCount);
			PrintTagOnGnuplotLabel(20, fOutput, "A[m2]",		fOutputCount);
			PrintTagOnGnuplotLabel(20, fOutput, "dummy",		fOutputCount);
			PrintTagOnGnuplotLabel(20, fOutput, "dummy",		fOutputCount);

			if(key_species_names.size() == 1)
				PrintTagOnGnuplotLabel(20, fOutput, "Eta_" + key_species_names[0],	fOutputCount);
			else
				PrintTagOnGnuplotLabel(20, fOutput, "dummy",	fOutputCount);

			for(i=1;i<=iOutputSpecies.Size();i++)
				PrintTagOnGnuplotLabel(20, fOutput,  mix->names[iOutputSpecies[i]]+"_x",	fOutputCount);
			for(i=1;i<=iOutputSpecies.Size();i++)
				PrintTagOnGnuplotLabel(20, fOutput,  mix->names[iOutputSpecies[i]]+"_w",	fOutputCount);
	
            fOutput << endl;
            fOutput << endl;
        }

        if (iTwoEquationModel == true)
		{
            fSoot2E	<< setw(20) << left << "Tau[s](1)"
					<< setw(20) << left << "x[m](2)"
					<< setw(20) << left << "T[K](3)";

            soot2EModel->GnuPlotInterface(fSoot2E, 4);
			fSoot2E << endl << endl;
		}
/*
        if (mix->pah_manager.iPAH == true)
		{
            fPAH	<< setw(20) << left << "Tau[s](1)"
					<< setw(20) << left << "x[m](2)"
					<< setw(20) << left << "T[K](3)";

			mix->pah_manager.GnuPlotInterface(fPAH, 4);
			fPAH << endl << endl;
		}
		*/
  /*      if (mix->soot_manager.iSoot == true)
        {
            fSoot	<< setw(20) << left << "Tau[s](1)"
					<< setw(20) << left << "x[m](2)"
					<< setw(20) << left << "T[K](3)";

            mix->soot_manager.GnuPlotInterface(fSoot, 4);
            fSoot << endl << endl;
        }
		*/
		if (iKeySpecies == true && key_species_names.size()>1)
        {
            fConversions	<< setw(20) << left << "Tau[s](1)"
							<< setw(20) << left << "x[m](2)"
							<< setw(20) << left << "T[K](3)";
			
			int fOutputCount = 4;
			for(i=0;i<int(key_species_names.size());i++)
				PrintTagOnGnuplotLabel(20, fConversions, "Eta_" + key_species_names[i],	fOutputCount);
			fConversions << endl;
			fConversions << endl;
        }
		
		if (iVerboseEnergy == true)
		{
			int fOutputCount=1;
			PrintTagOnGnuplotLabel(20, fEnergy, "Tau[s]",		fOutputCount);
			PrintTagOnGnuplotLabel(20, fEnergy, "x[m]",			fOutputCount);
			PrintTagOnGnuplotLabel(20, fEnergy, "U[J/s]",		fOutputCount);
			PrintTagOnGnuplotLabel(20, fEnergy, "H[J/kg]" ,		fOutputCount);
			PrintTagOnGnuplotLabel(20, fEnergy, "Ek[J/kg]",		fOutputCount);
			PrintTagOnGnuplotLabel(20, fEnergy, "Ep[J/kg]",		fOutputCount);
			PrintTagOnGnuplotLabel(20, fEnergy, "QR[W/m3]",		fOutputCount);
			PrintTagOnGnuplotLabel(20, fEnergy, "dummy",		fOutputCount);
			PrintTagOnGnuplotLabel(20, fEnergy, "Cp[J/kg/K]",	fOutputCount);
			fEnergy << endl;
			fEnergy << endl;
		}


		if (iVerboseReactionRates == true)
		{
			if (iGlobalKinetics == false)
			{
				for (i=1;i<=index_reaction_rates.Size();i++)
					fReactionRates << "[" << index_reaction_rates[i] << "]  " <<  names_reaction_rates[i-1] << endl;
				fReactionRates << endl << endl;
			}

            fReactionRates	<< setw(20) << left << "Tau[s](1)"
							<< setw(20) << left << "x[m](2)"
							<< setw(20) << left << "T[K](3)";
			
			int fOutputCount = 4;
			if (iGlobalKinetics == false)
				for (i=1;i<=index_reaction_rates.Size();i++)
				{
					stringstream number_reaction;
					number_reaction << index_reaction_rates[i];
					PrintTagOnGnuplotLabel(20, fReactionRates, "r_" + number_reaction.str(),	fOutputCount);
				}
			else
				for (i=1;i<=global->nReactions;i++)
				{
					stringstream number_reaction;
					number_reaction << i;
					PrintTagOnGnuplotLabel(20, fReactionRates, "r_" + number_reaction.str(),	fOutputCount);
				}
			
			fReactionRates << endl;
			fReactionRates << endl;
		}

		if (iVerboseFormationRates == true)
		{
            fFormationRates	<< setw(20) << left << "Tau[s](1)"
							<< setw(20) << left << "x[m](2)"
							<< setw(20) << left << "T[K](3)";
			
			int fOutputCount = 4;
			for (i=1;i<=index_formation_rates.Size();i++)
				PrintTagOnGnuplotLabel(20, fFormationRates, names_formation_rates[i-1],	fOutputCount);

			fFormationRates << endl;
			fFormationRates << endl;
		}
    }
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
//																								   //
//									SOLVING FUNCTIONS											   //
//																								   //
/////////////////////////////////////////////////////////////////////////////////////////////////////

void OpenSMOKE_ShockTube::Initialize()
{
	shockStream.AssignKineticScheme(*mix);
	ics.SetPointer();
	ics.SetGasMixture(mix);

	if (iIncidentShock == true)
	{
		if (iInitialCondition == SHOCK_BEFORE)
		{
			ics.Set_IncidentShock_CaseA(UShock, inletStream->T, inletStream->P, inletStream->omega);
			ics.Solve_IncidentShock_CaseA(shockStream);
			ics.VideoSummary();
	
			if (iBoundaryLayerCorrection == true)
				lm = ics.LengthBoundaryLayerCorrection(D);
		}
		
		else if (iInitialCondition == SHOCK_AFTER)
		{
			ics.Set_IncidentShock_CaseB(UShock, inletStream->T, inletStream->P, inletStream->omega);
			ics.Solve_IncidentShock_CaseB(shockStream);
			ics.VideoSummary();
	
			if (iBoundaryLayerCorrection == true)
				lm = ics.LengthBoundaryLayerCorrection(D);
		}
	}

	if (iIncidentShock == false)
	{
		if (iInitialCondition == SHOCK_REFLECTED)
		{
			ics.Set_ReflectedShock_CaseA(inletStream->T, inletStream->P, inletStream->omega);
			ics.Solve_ReflectedShock_CaseA(shockStream);
			ics.VideoSummary();
		}
		
		else if (iInitialCondition == SHOCK_BEFORE)
		{
			ics.Set_ReflectedShock_CaseB(shockStream, UShock, UReflectedShock, inletStream->T, inletStream->P, inletStream->omega);
			ics.Solve_ReflectedShock_CaseB(shockStream);
			ics.VideoSummary();
		}
		
		else if (iInitialCondition == SHOCK_AFTER)
		{
			ics.Set_ReflectedShock_CaseC(shockStream, UShock, UReflectedShock, inletStream->T, inletStream->P, inletStream->omega);
			ics.Solve_ReflectedShock_CaseC(shockStream);
			ics.VideoSummary();
		}
	}

	shockStream.lock();

	rho			 = shockStream.rho;
	v			 = shockStream.velocity;
	omega        = shockStream.omega;
	T            = shockStream.T;
	Csi			 = 0.;
	tL			 = 0.;
	P            = shockStream.P;

	UpdateProperties(MINUSONE, NC+3);

	// History
	BzzVector		ZeroVector;
	BzzMatrix		ZeroMatrix;
	Tau_History		= ZeroVector;
	Csi_History		= ZeroVector;
	T_History		= ZeroVector;
	mass_History	= ZeroMatrix;
	mole_History	= ZeroMatrix;
	ChangeDimensions(MAX_TIME_STEPS, &Tau_History);
    ChangeDimensions(MAX_TIME_STEPS, &Csi_History);
    ChangeDimensions(MAX_TIME_STEPS, &T_History);
    ChangeDimensions(MAX_TIME_STEPS, mix->NumberOfSpecies(), &mass_History);
    ChangeDimensions(MAX_TIME_STEPS, mix->NumberOfSpecies(), &mole_History);

	countGlobalIterations = -1;
}

void OpenSMOKE_ShockTube::ReSolve()
{
	Initialize();
	Solve();
}

void OpenSMOKE_ShockTube::DefineFromFile(const std::string inputFile)
{
    double  double_value;
    std::string  string_value;
    int     int_value;
	vector<string> string_vector;
	vector<double>  double_vector;

    OpenSMOKE_Dictionary_ShockTube dictionary;
    dictionary.ParseFile(inputFile);

   // SEMI-COMPULSORY: Kind Of Simulation
	if (dictionary.Return("#IncidentShock"))
		AssignIncidentShock();

   // SEMI-COMPULSORY: Kind Of Simulation
	if (dictionary.Return("#ReflectedShock"))
		AssignReflectedShock();

    // COMPULSORY: Reactor Residence Time
	if (dictionary.Return("#ResidenceTime", double_value, string_value))
		AssignEnd(string_value, double_value);

    // COMPULSORY: Shock wave velocity
	if (dictionary.Return("#ShockVelocity", double_value, string_value))
		AssignShockVelocity(string_value, double_value);

    // COMPULSORY: Reactor Diameter or Area
    if (dictionary.Return("#Diameter", double_value, string_value))
        AssignDiameter(string_value, double_value);
    else if (dictionary.Return("#Area", double_value, string_value))
        AssignArea(string_value, double_value);

    // COMPULSORY: Shock Temperature
	if (dictionary.Return("#TemperatureBefore", double_value, string_value))
		AssignShockTemperature(string_value, double_value, SHOCK_BEFORE);
	if (dictionary.Return("#TemperatureAfter", double_value, string_value))
		AssignShockTemperature(string_value, double_value, SHOCK_AFTER);
	if (dictionary.Return("#TemperatureReflected", double_value, string_value))
		AssignShockTemperature(string_value, double_value, SHOCK_REFLECTED);

    // COMPULSORY: Shock Pressure
	if (dictionary.Return("#PressureBefore", double_value, string_value))
		AssignShockPressure(string_value, double_value, SHOCK_BEFORE);
	if (dictionary.Return("#PressureAfter", double_value, string_value))
		AssignShockPressure(string_value, double_value, SHOCK_AFTER);
	if (dictionary.Return("#PressureReflected", double_value, string_value))
		AssignShockPressure(string_value, double_value, SHOCK_REFLECTED);

    // COMPULSORY: Shock Density
	if (dictionary.Return("#DensityBefore", double_value, string_value))
		AssignShockDensity(string_value, double_value, SHOCK_BEFORE);
	if (dictionary.Return("#DensityAfter", double_value, string_value))
		AssignShockDensity(string_value, double_value, SHOCK_AFTER);
	if (dictionary.Return("#DensityReflected", double_value, string_value))
		AssignShockDensity(string_value, double_value, SHOCK_REFLECTED);

    // SEMI-COMPULSORY: Mass fractions
    if (dictionary.Return("#MassFractions", double_vector, string_vector))
		AssignShockMassFractions(string_vector, double_vector);

    // SEMI-COMPULSORY: Mole fractions
    else if (dictionary.Return("#MoleFractions", double_vector, string_vector))
		AssignShockMoleFractions(string_vector, double_vector);

    // SEMI-COMPULSORY: EquivalenceRatio
   if (dictionary.Return("#EquivalenceRatio", double_value))
		AssignShockEquivalenceRatio(double_value);
    if (dictionary.Return("#FuelMassFractions", double_vector, string_vector))
		AssignShockFuelMassFractions(string_vector, double_vector);
    if (dictionary.Return("#FuelMoleFractions", double_vector, string_vector))
		AssignShockFuelMoleFractions(string_vector, double_vector);
    if (dictionary.Return("#OxidizerMassFractions", double_vector, string_vector))
		AssignShockOxidizerMassFractions(string_vector, double_vector);
    if (dictionary.Return("#OxidizerMoleFractions", double_vector, string_vector))
		AssignShockOxidizerMoleFractions(string_vector, double_vector);
	AssignShockStreamComposition();

	// OPTIONAL: KeySpecies
    if (dictionary.Return("#BoundaryLayerCorrection", string_vector))
        SetBoundaryLayerCorrection();

	// OPTIONAL: Reflected Shock Velocity
	if (dictionary.Return("#ReflectedShockVelocity", double_value, string_value))
		SetReflectedShockVelocity(string_value, double_value);

	// OPTIONAL: KeySpecies
    if (dictionary.Return("#Key", string_vector))
        SetKeySpecies(string_vector);

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

    // OPTIONAL: Global kinetics
    if (dictionary.Return("#GlobalKinetics"))
        SetGlobalKinetics();

    // OPTIONAL: Relative tolerance
    if (dictionary.Return("#RelativeTolerance", double_value))
        SetRelativeTolerance(double_value);

    // OPTIONAL: Absolute tolerance
    if (dictionary.Return("#AbsoluteTolerance", double_value))
        SetAbsoluteTolerance(double_value);

	// OPTIONAL: Verbose Energy
    if (dictionary.Return("#VerboseEnergy"))
        SetVerboseEnergy();

	// OPTIONAL: Inlet Viscosity
    if (dictionary.Return("#Viscosity", double_value, string_value))
        SetViscosity(double_value, string_value);

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
    if (dictionary.Return("#Sensitivity", string_vector))
		SetSensitivityOnFile(string_vector);

	// OPTIONAL: OutputSpecies
	if (dictionary.Return("#OutputSpecies", string_vector))
		SetOutputSpecies(string_vector);

	Lock();
}

void OpenSMOKE_ShockTube::Solve()
{
    ptShockTube = this;

    double timeStart, timeEnd;
    BzzVector xMin;
    BzzVector xMax;
    BzzVector xInitial;

	PrepareFiles();

    ODE_ShockTube_Object.assignShockTube(this, iBoundaryLayerCorrection);

    // 1.A Incident Shock Tube
    if (iTwoEquationModel == false)
    {
        ChangeDimensions(NC+5, &xMin);  xMin=ZERO;					
        ChangeDimensions(NC+5, &xMax);  xMax=ONE;				// Mass fractions
                                        xMax[NC+1] = 1e3;		// Density
                                        xMax[NC+2] = 1e5;		// Velocity
                                        xMax[NC+3] = 7000.;		// Temperature
                                        xMax[NC+4] = 1e16;		// Axial coordinate
                                        xMax[NC+5] = 1e16;		// Laboratory time

		
        xInitial = omega;
        xInitial.Append(rho);
        xInitial.Append(v);
        xInitial.Append(T);
        xInitial.Append(0.);
        xInitial.Append(0.);
    }
	    // 1.A Incident Shock Tube
    if (iTwoEquationModel == true)
    {
		soot2EModel->initial_values(rho);
        indexTwoEquations = NC+6;

        ChangeDimensions(NC+7, &xMin);  xMin=ZERO;					
        ChangeDimensions(NC+7, &xMax);  xMax=ONE;				// Mass fractions
                                        xMax[NC+1] = 1e3;		// Density
                                        xMax[NC+2] = 1e5;		// Velocity
                                        xMax[NC+3] = 7000.;		// Temperature
                                        xMax[NC+4] = 1e16;		// Axial coordinate
                                        xMax[NC+5] = 1e16;		// Laboratory time
                                        xMax[NC+6] = ONE;		// Soot phiN
                                        xMax[NC+7] = ONE;		// Soot phiM

		
        xInitial = omega;
        xInitial.Append(rho);
        xInitial.Append(v);
        xInitial.Append(T);
        xInitial.Append(0.);
        xInitial.Append(0.);
        xInitial.Append(soot2EModel->phiNStart);
        xInitial.Append(soot2EModel->phiMStart);
    }

    BzzOdeStiffObject o(xInitial, 0., &ODE_ShockTube_Object);

    o.StepPrint(ODEPrintExternal_ShockTube);
    o.SetMinimumConstraints(xMin);
    o.SetMaximumConstraints(xMax);
	//o.SetMaxStep(MAX_TIME_STEPS);

    if (iRelativeTolerance == true)	o.SetTolRel(relativeTolerance);
    if (iAbsoluteTolerance == true)	o.SetTolAbs(absoluteTolerance);

    {
        countGlobalIterations = -1;  // From -1 to avoid to store results from the first iteration

        timeStart = BzzGetCpuTime();
		
		o(TauTotal,TauTotal);

		status = o.GetCalculationState();
        
		timeEnd = BzzGetCpuTime();
    }


    if (iVerbose == true)
    {
        cout << endl;
        cout << "Number of function for Jacobian: " << o.GetNumFunctionForJacobian()		<< endl;
        cout << "Numerical Jacobians: "				<< o.GetNumNumericalJacobian()			<< endl;
        cout << "Time DAE solution: "				<< timeEnd - timeStart	<< " s"						<< endl << endl;
    }

	if (iHistory == true)
	{
		SaveOnBinaryFile(outputOSMName);
	}

	CloseFiles();
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
//																								   //
//									UPDATING FUNCTIONS											   //
//																								   //
/////////////////////////////////////////////////////////////////////////////////////////////////////

void OpenSMOKE_ShockTube::UpdateProperties(int jacobianIndex, int indexT)
{
	int memoIndex;

	if		(jacobianIndex<0)			memoIndex = -2;		// Calculate all without storing
	else if (jacobianIndex==0)			memoIndex = -1;		// Calculate and storing
	else if (jacobianIndex >= indexT)	memoIndex = -2;		// Calculate all without storing
	else								memoIndex = 0;		// Recover from maps

	double cTot;

	// Mole fractions
	// ----------------------------------------------------------------------
	mix->GetMWAndMoleFractionsFromMassFractions(MWtot, x, omega);

	// --------------------------------------------------------------------------
	// Update every variable which does not depend on T
	// --------------------------------------------------------------------------
	if(memoIndex == -2)
	{
		// Calcolo della pressione e della concentrazione [Pa] [kmol/m3]
		// ----------------------------------------------------------------------
		P		= rho*Constants::R_J_kmol*T/MWtot;  
		cTot	= P  / (Constants::R_J_kmol*T);
		c		= cTot*x;

		// Specific Heat [J/kgK]
		// ----------------------------------------------------------------------
		mix->SpeciesCp(T);
		Cp = mix->MixCp_FromMassFractions(omega);

		// Table Creation (Memorization)
		// ----------------------------------------------------------------------
		mix->ComputeKineticParameters( T, log(T), 1./T, P);
		mix->ComputeFromConcentrations( T, c, cTot, &R);				// [kmol/m3/s]
		ElementByElementProduct(R, mix->M, &R);							// [kg/m3/s]
		QReaction = - mix->ComputeQReaction(T);							// [J/m3.s]
	}

	// --------------------------------------------------------------------------
	// Update only the variables depending on the temperature
	// --------------------------------------------------------------------------
	else if(memoIndex==-1)
	{
		// Specific Heat [J/kgK]
		// ----------------------------------------------------------------------
		mix->SpeciesCp(T);
		CpMap = mix->Cp;

		// Table Creation (Memorization)
		// ----------------------------------------------------------------------
		mix->ComputeKineticParameters( T, log(T), 1./T, P);
		k1Map			= mix->k1;
		k2Map			= mix->k2;
		uKeqMap			= mix->uKeq;
		logFcentMap		= mix->logFcent;
		reactionDHMap	= mix->reactionDH;
		reactionDSMap	= mix->reactionDS;

		memoIndex = 0;
	}

	// --------------------------------------------------------------------------
	// Update every variable which does not depend on T
	// --------------------------------------------------------------------------
	if(memoIndex == 0)
	{
		// Calcolo della pressione e della concentrazione [Pa] [kmol/m3]
		// ----------------------------------------------------------------------
		P		= rho*Constants::R_J_kmol*T/MWtot;  
		cTot	= P  / (Constants::R_J_kmol*T);
		c		= cTot*x;

		// Specific Heat [J/kgK]
		// ----------------------------------------------------------------------
		mix->Cp = CpMap;
		Cp = mix->MixCp_FromMassFractions(omega);


		// Table Creation (Memorization)
		// ----------------------------------------------------------------------
		mix->k1			= k1Map;
		mix->k2			= k2Map;
		mix->uKeq		= uKeqMap;
		mix->logFcent	= logFcentMap;
		mix->reactionDH	= reactionDHMap;
		mix->reactionDS	= reactionDSMap;
		mix->ComputeFromConcentrations( T, c, cTot, &R);		// [kmol/m3/s]
		ElementByElementProduct(R, mix->M, &R);					// [kg/m3/s]
		QReaction = - mix->ComputeQReaction(T);					// [J/m3.s]
	}

	// Global Kinetics
	if (iGlobalKinetics == true)
	{
		global->GiveMeFormationRates(T,c,R);
		global->GiveMeReactionHeat(T, R, QReaction);
	}
}

void OpenSMOKE_ShockTube::UpdateTwoEquationModel(BzzVector &y, BzzVector &dy)
{
    soot2EModel->update(T, P/101325., rho, 0., x, y[indexTwoEquations], y[indexTwoEquations+1]);
    soot2EModel->formation_rates();

    for (int i=1;i<=NC-1;i++)
		domega[i] += soot2EModel->SGas[i] / rho;	// Gas species
	domega[NC] += soot2EModel->S / rho;				// Soot

    dy[indexTwoEquations]	= soot2EModel->s / rho;
    dy[indexTwoEquations+1] = soot2EModel->S / rho;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
//																								   //
//						ODE SYSTEM EQUATIONS - ISOTHERMAL REACTOR								   //
//																								   //
/////////////////////////////////////////////////////////////////////////////////////////////////////

void OpenSMOKE_ShockTube::ODESystem_Isothermal_ShockTube(BzzVector &y, double t, BzzVector &dy)
{
    int i;
	double AreaCorrection = 0.;

    for (i=1;i<=NC;i++)	omega[i]    = y[i];
    rho								= y[NC+1];
    v								= y[NC+2];
    T								= y[NC+3];
    Csi								= y[NC+4];
    tL								= y[NC+5];

    // Updating area
	if (iBoundaryLayerCorrection == true)
		AreaCorrection=AreaBoundaryLayerCorrection();

    // Evaluation of properties
    UpdateProperties(ODE_ShockTube_Object.jacobianIndex, NC+3);

    // Conservation equations for all the species
    domega	= R / rho;

    // Continuity equation
	double coefficient_rho = P+rho*v*v/Cp/T-rho*v*v;
	double sumReactions = Dot(R,mix->uM);
    drho	= ( -Constants::R_J_kmol*rho*(QReaction/Cp/MWtot+T*sumReactions)      + 
				 rho*rho*v*v*v*(1.-Constants::R_J_kmol/Cp/MWtot)*AreaCorrection ) / 
				coefficient_rho ;

	// Momentum equation
	dv		= -v/rho*drho - v*v*AreaCorrection;

	// Energy equation
	dT		= (v*v*drho+QReaction)/(rho*Cp) + v*v*v/Cp*AreaCorrection;

	// Axial coordinate
	dCsi	= v;

	// Axial coordinate
	dtL		= ics.rho1/rho;

    // Soot Two Equation Model
    if (iTwoEquationModel == true)
        UpdateTwoEquationModel(y, dy);

    // Recovering residuals (Soot 2E are already updated)
    for (i=1;i<=NC;i++)	dy[i]   = domega[i];
    dy[NC+1]                    = drho;
    dy[NC+2]                    = dv;
    dy[NC+3]                    = dT;
    dy[NC+4]                    = dCsi;
    dy[NC+5]                    = dtL;
}

double OpenSMOKE_ShockTube::AreaBoundaryLayerCorrection()
{
	double CsiX	= lm/10000.; 
	double tg	= 0.50*(1.+tanh(10.*(Csi-0.50*CsiX)/CsiX)); 
	
	double coeff = sqrt(Csi/lm);
	double AreaCorrection = Csi>0. ? 1./(2.*coeff*lm*(1.-coeff))*tg : 0.;
	Area = Area0/(1.-coeff);
	D	 = sqrt(Area/(Constants::pi/4.));

	return AreaCorrection;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
//																								   //
//									DICTIONARY													   //
//																								   //
/////////////////////////////////////////////////////////////////////////////////////////////////////

OpenSMOKE_Dictionary_ShockTube::OpenSMOKE_Dictionary_ShockTube()
{
    SetupBase();
	SetName("OpenSMOKE_ShockTube Dictionary");
    
	Add("#IncidentShock",			'O', 'N', "Incident Shock Simulation");
	Add("#ReflectedShock",			'O', 'N', "Reflected Shock Simulation");

    Add("#ResidenceTime",			'C', 'M', "Reactor residence time");
    Add("#ShockVelocity",			'O', 'M', "Shock wave velocity");

	Add("#TemperatureBefore",		'O', 'M', "Temperature before shock");
	Add("#TemperatureAfter",		'O', 'M', "Temperature after shock");
	Add("#TemperatureReflected",	'O', 'M', "Temperature of reflected shock");

	Add("#PressureBefore",			'O', 'M', "Pressure before shock");
	Add("#PressureAfter",			'O', 'M', "Pressure after shock");
	Add("#PressureReflected",		'O', 'M', "Pressure of reflected shock");

	Add("#DensityBefore",			'O', 'M', "Density before shock");
	Add("#DensityAfter",			'O', 'M', "Density after shock");
	Add("#DensityReflected",		'O', 'M', "Density of reflected shock");

    Add("#Diameter",				'O', 'M', "Reactor diameter");
    Add("#Area",					'O', 'M', "Reactor cross section area");
     
	Add("#MassFractions",			'O', 'L', "Mass fractions");
	Add("#MoleFractions",			'O', 'L', "Mole fractions");

	Add("#EquivalenceRatio",		'O', 'D', "Equivalence ratio");
	Add("#FuelMassFractions",		'O', 'L', "Fuel stream mass fractions");
    Add("#FuelMoleFractions",		'O', 'L', "Fuel stream mole fractions");
	Add("#OxidizerMassFractions",   'O', 'L', "Oxidizer stream mass fractions");
    Add("#OxidizerMoleFractions",	'O', 'L', "Oxidizer stream mole fractions");

	// Optional
	Add("#BoundaryLayer",			'O', 'N', "Boundary Layer Correction");
	Add("#ReflectedShockVelocity",	'O', 'M', "Reflected shock wave velocity");
	Add("#VerboseEnergy",			'O', 'N', "Report on energy of ShockTube");
    Add("#Key",						'O', 'V', "Key species");
	Add("#Viscosity",				'O', 'M', "Inlet viscosity");
	Add("#ReactionRates",			'O', 'V', "Reaction rates on file (list of reaction indices)");
	Add("#FormationRates",			'O', 'V', "Formation rates on file (list of species names) ");
	Add("#ROPA",				'O', 'N', "Rate of Production Analysis");
	Add("#VerboseROPA",			'O', 'V', "Rate of Production Analysis (list of species names)");
	Add("#Sensitivity",				'O', 'V', "Sensitivity Analysis");
	Add("#OutputSpecies",			'O', 'V', "Output species");

	Add("#RelativeTolerance",		'O', 'D', "Relative tolerance (default: 1.2e-5)");
	Add("#AbsoluteTolerance",		'O', 'D', "Absolute tolerance (default: 1.0e-10)");

	// Compulsory
	Compulsory("#Diameter", "#Area");
	Compulsory("#IncidentShock", "#ReflectedShock");
	Compulsory("#MassFractions", "#MoleFractions", "#EquivalenceRatio");

	// Conflicts
	Conflict("#TemperatureBefore",		"#TemperatureAfter");
	Conflict("#TemperatureBefore",		"#TemperatureReflected");
	Conflict("#TemperatureReflected",	"#TemperatureAfter");
	
	Conflict("#PressureBefore",			"#PressureAfter");
	Conflict("#PressureBefore",			"#PressureReflected");
	Conflict("#PressureReflected",		"#PressureAfter");

	Conflict("#DensityBefore",			"#DensityAfter");
	Conflict("#DensityBefore",			"#DensityReflected");
	Conflict("#DensityReflected",		"#DensityAfter");

	Conflict("#IncidentShock",			"#TemperatureReflected");
	Conflict("#IncidentShock",			"#PressureReflected");
	Conflict("#IncidentShock",			"#DensityReflected");

	Conflict("#Diameter",		"#Area");
	Conflict("#IncidentShock",	"#ReflectedShockVelocity");
	Conflict("#ReflectedShock",	"#BoundaryLayer");

	Conflict("#EquivalenceRatio", "#MoleFractions");
	Conflict("#EquivalenceRatio", "#MassFractions");
	Conflict("#MassFractions", "#MoleFractions");

	
	Conflict("#EquivalenceRatio", "#MassFractions");
	Conflict("#EquivalenceRatio", "#MoleFractions");
	Conflict("#FuelMassFractions", "#MoleFractions");
	Conflict("#FuelMassFractions", "#MassFractions");
	Conflict("#OxidizerMassFractions", "#MoleFractions");
	Conflict("#OxidizerMassFractions", "#MassFractions");
	Conflict("#FuelMassFractions", "#FuelMoleFractions");
	Conflict("#OxidizerMoleFractions", "#OxidizerMassFractions");

	
    Lock();
}

void OpenSMOKE_ShockTube::SaveOnBinaryFile(BzzSave &fOutput)
{
	std::string dummy;
	char name[Constants::NAME_SIZE];

	int nSteps = countGlobalIterations;			
	BzzVector timegrid = GetBzzVector(nSteps, 1, Tau_History);
	BzzVector csigrid = GetBzzVector(nSteps, 1, Csi_History);
	BzzVector tgrid    = GetBzzVector(nSteps, 1, T_History);

	dummy = "SHOCKTUBE";
	strcpy(name, dummy.c_str());
	fOutput.fileSave.write((char*) name, sizeof(name));

	dummy = "V20100417";
	strcpy(name, dummy.c_str());
	fOutput.fileSave.write((char*) name, sizeof(name));

	// Time
	dummy = "TIME";
	strcpy(name, dummy.c_str());
	fOutput.fileSave.write((char*) name, sizeof(name));
	fOutput << timegrid;

	// Space coordinate
	dummy = "CSI";
	strcpy(name, dummy.c_str());
	fOutput.fileSave.write((char*) name, sizeof(name));
	fOutput << csigrid;

	// Temperature
	dummy = "TEMPERATURE";
	strcpy(name, dummy.c_str());
	fOutput.fileSave.write((char*) name, sizeof(name));
	fOutput << tgrid;

	// Pressure
	BzzVector P_Pa(nSteps); P_Pa = P;
	dummy = "PRESSURE";
	strcpy(name, dummy.c_str());
	fOutput.fileSave.write((char*) name, sizeof(name));
	fOutput << P_Pa;

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

void OpenSMOKE_ShockTube::SaveOnBinaryFile(const std::string filename)
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

		std::string dummy;
		char name[Constants::NAME_SIZE];

		BzzVector RRvector(mix->NumberOfReactions());
		BzzMatrix RRmatrix(countGlobalIterations, mix->NumberOfReactions());
		for(int j=1;j<=countGlobalIterations;j++)
		{
			BzzVector aux = mole_History.GetRow(j);
			UpdateReactionRates(T_History[j], P, aux, RRvector);
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
		cout << "  -- sensitivity analysis..." << endl;
	}

	PrintEndOnBinaryFile(fOutput);
	fOutput.End();
			
	// Post processor test
	//OpenSMOKE_PostProcessor post_processor;
	//post_processor.ReadFromBinaryFile(filename);
}

void OpenSMOKE_ShockTube::UpdateReactionRates(const double TT, const double PP, BzzVector &xx, BzzVector &rr)
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

/////////////////////////////////////////////////////////////////////////////////////////////////////
//																								   //
//									ODE SYSTEM - CLASS DEFINITION								   //
//																								   //
/////////////////////////////////////////////////////////////////////////////////////////////////////

void MyOdeSystem_ShockTube::ObjectBzzPrint(void)
{
    ::BzzPrint("\n OpenSMOKE_ShockTube Class\n");
}

void MyOdeSystem_ShockTube::GetSystemFunctions(BzzVector &x, double t, BzzVector &f)
{
	ptShockTube->ODESystem_Isothermal_ShockTube(x, t, f);
}

void MyOdeSystem_ShockTube::assignShockTube(OpenSMOKE_ShockTube *shocktube, bool _iBoundaryLayerCorrection)
{
    ptShockTube = shocktube;
	iBoundaryLayerCorrection = _iBoundaryLayerCorrection;
}

void OpenSMOKE_ShockTube::AssignShockEquivalenceRatio(const double value)
{
	iShockEquivalenceRatio = true;
	equivalence_ratio = value;
}

void OpenSMOKE_ShockTube::AssignShockFuelMassFractions(const vector<string> names, const vector<double> values)
{
	int i;
	double MWmix;
	BzzVector y(mix->NumberOfSpecies());
	BzzVector x(mix->NumberOfSpecies());
	
	for(i=1;i<=int(names.size());i++)
		y[mix->recognize_species(names[i-1])] = values[i-1];
	mix->GetMWAndMoleFractionsFromMassFractions(MWmix, x, y);

	vector<string> _names;
	vector<double> _values;
	for(i=1;i<=mix->NumberOfSpecies();i++)
		if (x[i]>0.)
		{
			_names.push_back(mix->names[i]);
			_values.push_back(x[i]);
		}
	AssignShockFuelMoleFractions(_names, _values);
}
 
void OpenSMOKE_ShockTube::AssignShockFuelMoleFractions(const vector<string> names, const vector<double> values)
{
	fuel_names.resize(names.size()+1);
	ChangeDimensions(names.size(), &moles_fuel);
	for(int i=1;i<=int(names.size());i++)
	{
		fuel_names[i] = names[i-1];
		moles_fuel[i] = values[i-1];
	}
}

void OpenSMOKE_ShockTube::AssignShockOxidizerMassFractions(const vector<string> names, const vector<double> values)
{
	int i;
	double MWmix;
	BzzVector y(mix->NumberOfSpecies());
	BzzVector x(mix->NumberOfSpecies());
	
	for(i=1;i<=int(names.size());i++)
		y[mix->recognize_species(names[i-1])] = values[i-1];
	mix->GetMWAndMoleFractionsFromMassFractions(MWmix, x, y);

	vector<string> _names;
	vector<double> _values;
	for(i=1;i<=mix->NumberOfSpecies();i++)
		if (x[i]>0.)
		{
			_names.push_back(mix->names[i]);
			_values.push_back(x[i]);
		}
	AssignShockOxidizerMoleFractions(_names, _values);
}
 
void OpenSMOKE_ShockTube::AssignShockOxidizerMoleFractions(const vector<string> names, const vector<double> values)
{
	oxidizer_names.resize(names.size()+1);
	ChangeDimensions(names.size(), &moles_oxidizer);
	for(int i=1;i<=int(names.size());i++)
	{
		oxidizer_names[i] = names[i-1];
		moles_oxidizer[i] = values[i-1];
	}
}

void OpenSMOKE_ShockTube::AssignShockStreamComposition()
{
	if (iShockEquivalenceRatio == true)
	{
		OpenSMOKE_GasStream gas_stream;
		gas_stream.AssignKineticScheme(*mix);
		gas_stream.AssignTemperature(Constants::T_Reference, "K");
		gas_stream.AssignPressure(Constants::P_Reference, "Pa");
		gas_stream.AssignMassFlowRate(1., "kg/s");
		gas_stream.AssignEquivalenceRatio( equivalence_ratio, fuel_names, moles_fuel, oxidizer_names, moles_oxidizer);
		gas_stream.lock();

		vector<string> names;
		vector<double> values;
		for(int i=1;i<=mix->NumberOfSpecies();i++)
			if (gas_stream.x[i] > 0.)
			{
				names.push_back(mix->names[i]);
				values.push_back(gas_stream.x[i]);
			}

		AssignShockMoleFractions(names, values);
	}
}