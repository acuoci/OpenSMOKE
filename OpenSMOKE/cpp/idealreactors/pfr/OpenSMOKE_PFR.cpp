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
#include "basic/OpenSMOKE_Conversions.h"
#include "addons/OpenSMOKE_PostProcessor.h"
#include "addons/OpenSMOKE_SensitivityAnalysis_Fast_Unsteady0D.h"
#include "addons/OpenSMOKE_PolimiSoot.h"

OpenSMOKE_PFR::OpenSMOKE_PFR()
{
	class_name	= "OpenSMOKE_PFR";
	out_name	= "PFR.out";
	
	Setup();

	assignedGeometry			= false;	//
	iTimeIndipendent            = false;    // false = Space as indipendent variable true = Time as indipendent variable
	iMomentum			        = false;    // Momentum equation enabled
	iLogExperiment				= false;

	initialTime					= 0.;
	initialLenght				= 0.;

	Beta						= 1.0;		//
	iODE_Solver					= ODE_SOLVER_BZZ;
}

void OpenSMOKE_PFR::AssignEnd(const string units, const double value)
{
	if (iTimeIndipendent == true)
		TauTotal = OpenSMOKE_Conversions::conversion_time(value, units);
    else
		L = OpenSMOKE_Conversions::conversion_length(value, units);
	
    assignedEnd = true;
}

void OpenSMOKE_PFR::AssignDiameter(const string units, const double value)
{
    D       = OpenSMOKE_Conversions::conversion_length(value, units);
    Area    = Constants::pi * D*D / 4.;
	
	if (iUserDefinedExchangeArea == NONE)		Ae = 4./D;

    assignedGeometry = true;
}

void OpenSMOKE_PFR::AssignArea(const string units, const double value)
{
    Area    = OpenSMOKE_Conversions::conversion_area(value, units);
    D       = sqrt(4./Constants::pi*Area);

	if (iUserDefinedExchangeArea == NONE)		Ae = 4./D;

    assignedGeometry = true; 
}

void OpenSMOKE_PFR::SetUserDefinedDiameter(const string fileName)
{
    //geometry.Setup(fileName, "LENGTH");
    geometry.Setup(fileName);
	geometry.Update(0., D, Area);

	if (iUserDefinedExchangeArea == NONE)		Ae = 4./D;

	assignedGeometry = true;
}

void OpenSMOKE_PFR::SetUserDefinedArea(const string fileName)
{
    //geometry.Setup(fileName, "AREA");
    geometry.Setup(fileName);
	geometry.Update(0., D, Area);

	if (iUserDefinedExchangeArea == NONE)		Ae = 4./D;

	assignedGeometry = true;
}

void OpenSMOKE_PFR::SetMomentumEquation()
{
    iMomentum = true;
}

void OpenSMOKE_PFR::UnsetMomentumEquation()
{
    iMomentum = false;
}

void OpenSMOKE_PFR::SetUserDefinedSpecificExchangeArea(const string fileName)
{
    iUserDefinedExchangeArea = USERDEFINED;
    ud_Ae_profile.AssignFromFile(fileName, "SPECIFIC_AREA");
	ud_Ae_profile.SetName(name_object + " - Exchange Area Profile");
}

void OpenSMOKE_PFR::SetConstantSpecificExchangeArea(const double value, const string units)
{
	iUserDefinedExchangeArea = CONSTANT;
	Ae = OpenSMOKE_Conversions::conversion_u_length(value, units);
}

void OpenSMOKE_PFR::UnsetUserDefinedSpecificExchangeArea()
{
    iUserDefinedExchangeArea = NONE;
}

void OpenSMOKE_PFR::SetLogExperiment(const int int_value)
{
	iLogExperiment = true;
	nExperiments = int_value;
}

void OpenSMOKE_PFR::SetODESolver(const string solver)
{
	if (solver == "BZZ")			iODE_Solver = ODE_SOLVER_BZZ;
	else if (solver == "RADAU")		iODE_Solver = ODE_SOLVER_RADAU;
	else if (solver == "DLSODE")	iODE_Solver = ODE_SOLVER_DLSODE;
	else ErrorMessage("ODE Solvers: BZZ || RADAU || DLSODE");
}

void OpenSMOKE_PFR::SetInitialTime(const string units, const double value)
{
    initialTime       = OpenSMOKE_Conversions::conversion_time(value, units);
}

void OpenSMOKE_PFR::SetInitialLenght(const string units, const double value)
{
    initialLenght       = OpenSMOKE_Conversions::conversion_length(value, units);
}

void GnuPlotSootDistributionInterface(ofstream &fSoot);

void OpenSMOKE_PFR::Lock()
{
    if (assignedKineticScheme == false)
        ErrorMessage("The kinetic scheme was not defined!");
    if (assignedGlobalKineticScheme == false)
        ErrorMessage("The global kinetic scheme was not defined!");
    if (assignedEnd == false)
		if (iTimeIndipendent == true)	ErrorMessage("The reactor residence time was not defined!");
    	else							ErrorMessage("The reactor lenght was not defined!");
	if (assignedGeometry == false)
        ErrorMessage("The reactor geometry was not defined!");
    if (assignedInletFlows == false)
        ErrorMessage("The inlet flow was not defined!");
    if (assignedEnergy == false)
        ErrorMessage("The kind of reactor was not defined!");
	if (assignedSoot2EModel == false)
        ErrorMessage("The Soot 2E Model was not defined!");
	if (iTwoEquationModel == true && mix->iSootMode == false)
		ErrorMessage("This kinetic scheme was interpreted without the SootMode option and cannot be used together the 2E model.\nMore details on the OpenSMOKE's User Guide.");

	Initialize();

	if (mix->polimiSoot->IsSoot() == true)
	{
		const string fileName = outputFolderName + "/SootDistribution.out";
		openOutputFileAndControl(fSootDistribution, fileName);
		fSootDistribution.setf(ios::scientific);
		GnuPlotSootDistributionInterface(fSootDistribution);
	}
}

void GnuPlotSootDistributionInterface(ofstream &fSoot)
{
	fSoot	<< setw(20) << left << "t[s](1)"
			<< setw(20) << left << "bin(2)";

    fSoot	<< setw(20) << left << "d[m](3)"
			<< setw(20) << left << "m[?](4)"
			<< setw(20) << left << "V[?](5)"
			<< setw(20) << left << "log10d(6)"
			<< setw(20) << left << "log10m(7)"
			<< setw(20) << left << "log10V(8)"
			<< setw(20) << left << "dlog10d(9)"
			<< setw(20) << left << "dlog10m(10)"
			<< setw(20) << left << "dlog10V(11)";

    fSoot	<< setw(20) << left << "fv(12)"
			<< setw(20) << left << "x(13)"
			<< setw(20) << left << "y(14)"
			<< setw(20) << left << "rho[kg/m3](15)"
			<< setw(20) << left << "N[#/m3](16)";

    fSoot	<< setw(20) << left << "fvN(17)"
			<< setw(20) << left << "xN(18)"
			<< setw(20) << left << "yN(19)"
			<< setw(20) << left << "rhoN[](20)"
			<< setw(20) << left << "NN[](21)";

	fSoot << endl << endl;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
//																								   //
//									OUTPUT FUNCTIONS											   //
//																								   //
/////////////////////////////////////////////////////////////////////////////////////////////////////

void OpenSMOKE_PFR::VideoFinalResult()
{
    int i;

    // Vectors
    BzzVector conversion(NC);

    // Print General Information about PFR
    VideoGeneralInfo();

    // Calculate Conversion of Inlet Species
    for (i=1;i<=NC;i++)
        if (inletStream->omega[i]!=0.)
            conversion[i] = (1. - omega[i]/inletStream->omega[i])*100.;

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

void OpenSMOKE_PFR::VideoGeneralInfo()
{
    cout.setf(ios::scientific);
    cout << endl;

    cout << "xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx" << endl;
    cout << "xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx" << endl;
    cout << endl;

    cout << "---------------------------------------------------------------------" << endl;
    cout << " PFR Summary " << endl;
    cout << "---------------------------------------------------------------------" << endl;
    cout << "Length:              "		<< L								<< " [m]"		<< endl;
    cout << "Residence time:      "		<< Tau								<< " [s]"		<< endl;
	cout << "Mass flow rate:      "		<< inletStream->massFlowRate		<< " [kg/s]"	<< endl;
    cout << "Mole flow rate:      "		<< inletStream->moleFlowRate		<< " [kmol/s]"	<< endl;
    cout << "Volume flow rate:    "		<< inletStream->volumetricFlowRate	<< " [m3/s]"	<< endl;
	cout << "Pressure:            "		<< P								<< " [Pa]"		<< endl;
    cout << "Density:             "		<< rho								<< " [kg/m3]"	<< endl;
    cout << "Molecular weight:    "		<< MWtot							<< " [kg/kmol]" << endl;
    cout << "Velocity:            "		<< v								<< " [m/s]"		<< endl;
    cout << "Area:                "		<< Area								<< " [m2]"		<< endl;
    cout << "Temperature:         "		<< T								<< " [K]"		<< endl;
    cout <<  endl;
}

void OpenSMOKE_PFR::VideoSummary()
{
    VideoGeneralInfo();

    cout << setw(6) << left << "#" << setw(20) << left << "name" << setw(16) << left << "x" << setw(16) << left << "omega" << endl;
    cout << "---------------------------------------------------------------------" << endl;
    for (int i=1;i<=NC;i++)
        if (inletStream->omega[i]!=0.)
            cout << setw(6)  << left << i 
				 << setw(20) << left << mix->names[i] 
				 << setw(16) << left << inletStream->x[i] 
				 << setw(16) << left << inletStream->omega[i] 
				 << endl;
    cout << "---------------------------------------------------------------------" << endl;
    cout << endl << endl;
}

void OpenSMOKE_PFR::SummaryOnFile()
{
	int i;
	string file_name = outputFolderName + "/Summary.out";
	ofstream fSummary;
	openOutputFileAndControl(fSummary, file_name);
    fSummary.setf(ios::scientific);


    fSummary << "                  	Inlet        	  Outlet"											<< endl;
	fSummary << "----------------------------------------------------------------------------"		<< endl;
    fSummary << setw(20) << left << "Time[s]"			<< setw(16) << left << ZERO							<< setw(16) << left << Tau								<< endl;
	fSummary << setw(20) << left << "Length[m]"			<< setw(16) << left << ZERO							<< setw(16) << left << Csi								<< endl;	
	fSummary << setw(20) << left << "Temperature[K]"	<< setw(16) << left << inletStream->T				<< setw(16) << left << T								<< endl;
	fSummary << setw(20) << left << "Pressure[Pa]"		<< setw(16) << left << inletStream->P				<< setw(16) << left << P								<< endl;
	fSummary << setw(20) << left << "Mass_flow[kg/s]"	<< setw(16) << left << inletStream->massFlowRate	<< setw(16) << left << inletStream->massFlowRate		<< endl;
	fSummary << setw(20) << left << "Mole_flow[kmol/s]"	<< setw(16) << left << inletStream->moleFlowRate	<< setw(16) << left << inletStream->massFlowRate/MWtot	<< endl;	
	fSummary << setw(20) << left << "Density[kg/m3]"	<< setw(16) << left << inletStream->rho				<< setw(16) << left << rho								<< endl;
	fSummary << setw(20) << left << "MW[kg/kmol]"		<< setw(16) << left << inletStream->MW				<< setw(16) << left << MWtot							<< endl;
	fSummary << endl;
	
	
	fSummary << "Mole_fractions[-]" << endl;
	for(i=1;i<=NC;i++)
		fSummary << setw(20) << left << mix->names[i] << setw(16) << left << inletStream->x[i] << setw(16) << left << x[i]	<< endl;
	fSummary << endl;

	fSummary << endl;
	fSummary << "Mass_fractions[-]" << endl;
	for(i=1;i<=NC;i++)
		fSummary << setw(20) << left << mix->names[i] << setw(16) << left << inletStream->omega[i] << setw(16) << left << omega[i]	<< endl;
	fSummary << endl;
	
	fSummary << endl;
	fSummary << "Concentrations[kmol/m3]" << endl;
	for(i=1;i<=NC;i++)
		fSummary << setw(20) << left << mix->names[i] << setw(16) << left << inletStream->c[i] << setw(16) << left << c[i]	<< endl;
	fSummary << endl;

	fSummary << endl;
	fSummary << "Conversions[-]" << endl;
	for(i=1;i<=NC;i++)
		if (inletStream->omega[i] > 0.)
			fSummary << setw(20) << left << mix->names[i] << setw(16) << left << 1.-omega[i]/inletStream->omega[i]	<< endl;
	fSummary << endl;

    fSummary.close();

	if (iVerboseExperiment == true)
		ExperimentOnFile();
}

void OpenSMOKE_PFR::ExperimentOnFile()
{
	int		iStart = 20;
	int		nSteps = 10;
	double	threshold = 1.e-6;


	ofstream fExperiment;
	openOutputFileAndControl(fExperiment, outputExperimentName);
	fExperiment.setf(ios::scientific);

	fExperiment << "REACTOR PFR" << endl;
		
	fExperiment << "TIME   s" << endl;
	fExperiment << "DATA    "   << 1+names_Experiment.size() << endl;
	fExperiment << endl;

	if (iLogExperiment == false)
	{
		if (T_History[iStart]!=T_History[countGlobalIterations])
		{
			fExperiment << "T TEMPERATURE 1.0 CROSS" << endl;
			for(int i=iStart;i<=countGlobalIterations;i+=nSteps)
				fExperiment << Tau_History[i] << "\t" << T_History[i] << endl;
			fExperiment << "//" << endl << endl;
		}
		else
		{
			fExperiment << "T TEMPERATURE 1.0 CROSS" << endl;
			fExperiment << Tau_History[iStart] << "\t" << T_History[iStart] << endl;
			fExperiment << Tau_History[countGlobalIterations] << "\t" << T_History[countGlobalIterations] << endl;
			fExperiment << "//" << endl << endl;
		}

		for(int k=0;k<int(names_Experiment.size());k++)
		{
			fExperiment << names_Experiment[k] << " MOLE_FRACTION 1.0 CROSS" << endl;
			for(int i=iStart;i<=countGlobalIterations;i+=nSteps)
				if (mole_History[i][index_Experiment[k+1]] > threshold)
					fExperiment << Tau_History[i] << "\t" << mole_History[i][index_Experiment[k+1]] << endl;
			fExperiment << "//" << endl << endl;
		}
	}
	else
	{
		double alfa;

		BzzVector support_t(nExperiments);
		BzzVector interpolation_x(countGlobalIterations);
		BzzVector interpolation_y(countGlobalIterations);

		alfa = pow(Tau_History[countGlobalIterations-1]/Tau_History[iStart], 1./double(nExperiments));

		{
			int j;
			support_t[1] = Tau_History[iStart];
			for(j=2;j<=nExperiments+1;j++)
				support_t[j] = alfa*support_t[j-1];
			for(j=1;j<=countGlobalIterations;j++)
				interpolation_x[j] = Tau_History[j];
		}

		if (T_History[iStart]!=T_History[countGlobalIterations])
		{
			int i;

			for(i=1;i<=countGlobalIterations;i++)
				interpolation_y[i] = T_History[i];

			LinearInterpolation interpolation;
			interpolation(interpolation_x,interpolation_y);

			fExperiment << "T TEMPERATURE 1.0 CROSS" << endl;
			for(i=1;i<=nExperiments+1;i++)
				fExperiment << support_t[i] << "\t" << interpolation(support_t[i]) << endl;
			fExperiment << "//" << endl << endl;
		}
		else
		{
			fExperiment << "T TEMPERATURE 1.0 CROSS" << endl;
			fExperiment << Tau_History[iStart] << "\t" << T_History[iStart] << endl;
			fExperiment << Tau_History[countGlobalIterations] << "\t" << T_History[countGlobalIterations] << endl;
			fExperiment << "//" << endl << endl;
		}

		for(int k=0;k<int(names_Experiment.size());k++)
		{
			int i;

			for(i=1;i<=countGlobalIterations;i++)
				interpolation_y[i] = mole_History[i][index_Experiment[k+1]];

			LinearInterpolation interpolation;
			interpolation(interpolation_x,interpolation_y);

			fExperiment << names_Experiment[k] << " MOLE_FRACTION 1.0 CROSS" << endl;
			for(i=1;i<=nExperiments+1;i++)
				if (interpolation(support_t[i]) > threshold)
					fExperiment << support_t[i] << "\t" << interpolation(support_t[i]) << endl;
			fExperiment << "//" << endl << endl;
		}
	}

	fExperiment.close();
}

void OpenSMOKE_PFR::EnergyAnalysis(OpenSMOKE_GasStream &outletStream)
{
	double H_mass_inlet  = inletStream->massSpecificEnthalpy;			// [J/kg]
	double P_mass_inlet  = inletStream->massSpecificPressureEnergy;		// [J/kg]
	double U_mass_inlet  = inletStream->massSpecificInternalEnergy;		// [J/kg]
	double K_mass_inlet  = 0.50*BzzPow2(inletStream->velocity);			// [J/kg]
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

	cout << "Inlet power:     "			<< E_mass_inlet*massFlowRate					<< " [W]" << endl;
	cout << "Outlet power:    "			<< -E_mass_outlet*massFlowRate					<< " [W]" << endl;
	cout << "Exchanged power: "			<< -(E_mass_outlet-E_mass_inlet)*massFlowRate	<< " [W]" << endl;
}

void OpenSMOKE_PFR::ODEPrint(BzzVector &y, double eta)
{
	double dummy = ZERO;

    if (iVerbose == true)
    {
        if ( countIterations%(300*nVideoSteps) == 0)
        {
			if (inletStream->iUndefinedFlowRate == false)
				cout    << endl
						<< "#"			<< "\t"
						<< "Tau[s]"		<< "\t\t"
						<< "x[m]"		<< "\t\t"
						<< "T[K]"       << "\t\t"
						<< "v[m/s]"     << "\t\t";
			else
				cout    << endl
						<< "#"			<< "\t"
						<< "Tau[s]"		<< "\t\t"
						<< "T[K]"       << "\t\t";
			
					if (iMomentum == true)	cout << "P[atm]"	<< "\t\t";
                    
			cout	<< endl;
        }

        if (countVideoSteps == nVideoSteps)
        {
			if (inletStream->iUndefinedFlowRate == false)
				cout	<< countIterations	<< "\t"
						<< Tau  			<< "\t"
						<< Csi				<< "\t"
						<< T				<< "\t"
						<< v				<< "\t";
			else
				cout	<< countIterations	<< "\t"
						<< Tau  			<< "\t"
						<< T				<< "\t";

			if (iMomentum == true)	
				cout << P/101325.		<< "\t";

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
                    << setw(20) << left << dummy
                    << setw(20) << left << P/101325.
                    << setw(20) << left << v
                    << setw(20) << left << QReaction
                    << setw(20) << left << Qe*Ae
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

            if (mix->polimiSoot->IsBin() == true)
            {
                // Soot Distribution data
               /* mix->soot_manager.large_calculateAll(omega, x, c, rho);
                mix->soot_manager.small_calculateAll(omega, x, c, rho);

                fSoot << setw(20) << left << Tau
					  << setw(20) << left << Csi
					  << setw(20) << left << T;
                mix->soot_manager.print_main_data_on_file(fSoot);
                fSoot << endl;
				*/

				mix->polimiSoot->Analysis(*mix, P, T, omega);
				mix->polimiSoot->Distribution();
				mix->polimiSoot->ProcessDistribution();

                fSoot << setw(20) << left << Tau
					  << setw(20) << left << Csi
					  << setw(20) << left << T;

				fSoot << setw(20) << left << mix->polimiSoot->fv_large();
				fSoot << setw(20) << left << mix->polimiSoot->x_large();
				fSoot << setw(20) << left << mix->polimiSoot->omega_large();
				fSoot << setw(20) << left << mix->polimiSoot->rho_large();
				fSoot << setw(20) << left << mix->polimiSoot->N_large();
				fSoot << setw(20) << left << mix->polimiSoot->h_over_c_large();
				fSoot << setw(20) << left << mix->polimiSoot->dmean_N_large();
				fSoot << setw(20) << left << mix->polimiSoot->d32_N_large();
				fSoot << setw(20) << left << mix->polimiSoot->dstd_N();

				fSoot << setw(20) << left << mix->polimiSoot->fv_small();
				fSoot << setw(20) << left << mix->polimiSoot->x_small();
				fSoot << setw(20) << left << mix->polimiSoot->omega_small();
				fSoot << setw(20) << left << mix->polimiSoot->rho_small();
				fSoot << setw(20) << left << mix->polimiSoot->N_small();

				fSoot << endl;

				// Distribution
				{
					for (int k=1;k<=mix->polimiSoot->baskets_d().Size();k++)
					{
						fSootDistribution << setw(20) << left << Tau;
						fSootDistribution << setw(20) << left << k;
					
						fSootDistribution << setw(20) << left << mix->polimiSoot->baskets_d()[k];
						fSootDistribution << setw(20) << left << mix->polimiSoot->baskets_m()[k];
						fSootDistribution << setw(20) << left << mix->polimiSoot->baskets_V()[k];

						fSootDistribution << setw(20) << left << log10(mix->polimiSoot->baskets_d()[k]);
						fSootDistribution << setw(20) << left << log10(mix->polimiSoot->baskets_m()[k]);
						fSootDistribution << setw(20) << left << log10(mix->polimiSoot->baskets_V()[k]);

						fSootDistribution << setw(20) << left << log10(mix->polimiSoot->baskets_dlog10d()[k]);
						fSootDistribution << setw(20) << left << log10(mix->polimiSoot->baskets_dlog10m()[k]);
						fSootDistribution << setw(20) << left << log10(mix->polimiSoot->baskets_dlog10V()[k]);

						fSootDistribution << setw(20) << left << mix->polimiSoot->baskets_fv()[k];
						fSootDistribution << setw(20) << left << mix->polimiSoot->baskets_x()[k];
						fSootDistribution << setw(20) << left << mix->polimiSoot->baskets_omega()[k];
						fSootDistribution << setw(20) << left << mix->polimiSoot->baskets_rho()[k];
						fSootDistribution << setw(20) << left << mix->polimiSoot->baskets_N()[k];
					
						fSootDistribution << setw(20) << left << mix->polimiSoot->baskets_fv()[k]/(mix->polimiSoot->fv_large()+mix->polimiSoot->fv_small());
						fSootDistribution << setw(20) << left << mix->polimiSoot->baskets_x()[k]/(mix->polimiSoot->x_large()+mix->polimiSoot->x_small());
						fSootDistribution << setw(20) << left << mix->polimiSoot->baskets_omega()[k]/(mix->polimiSoot->omega_large()+mix->polimiSoot->omega_small());
						fSootDistribution << setw(20) << left << mix->polimiSoot->baskets_rho()[k]/(mix->polimiSoot->rho_large()+mix->polimiSoot->rho_small());
						fSootDistribution << setw(20) << left << mix->polimiSoot->baskets_N()[k]/(mix->polimiSoot->N_large()+mix->polimiSoot->N_small());

						fSootDistribution << endl;
					}
				}
				fSootDistribution << endl;
            }

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
						<< setw(20) << left << Qe*Ae
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

			if (iVerboseSelectivity == true)
			{
				int index_c = mix->recognize_element("c");

			//	cout << "index key species: " << key_species_index[0] << endl;
			//	cout << "inlet mass: " << inletStream->omega[key_species_index[0]] << endl;
			//	cout << "n carbon in key species: " << mix->elements[index_c][key_species_index[0]] << endl;
			//	getchar();

				double conversion = 1.-omega[key_species_index[0]]/inletStream->omega[key_species_index[0]];
				fSelectivity	<< setw(20) << left << Tau
								<< setw(20) << left << Csi
								<< setw(20) << left << T
								<< setw(20) << left << conversion;
			
				for (i=1;i<=NC;i++)
					selectivity_mass[i] = omega[i]/inletStream->omega[key_species_index[0]]/conversion;

				for (i=1;i<=NC;i++)
					fSelectivity	<< setw(20) << left << selectivity_mass[i];
				for (i=1;i<=NC;i++)
					fSelectivity	<< setw(20) << left << selectivity_mass[i]*mix->M[key_species_index[0]]/mix->M[i];
				for (i=1;i<=NC;i++)
					fSelectivity	<< setw(20) << left << selectivity_mass[i]*mix->M[key_species_index[0]]/mix->M[i]*mix->elements[index_c][i]/mix->elements[index_c][key_species_index[0]];

				fSelectivity	<< endl;
			}

			if (iAssignedROPA == true)
			{
			}

			if (iVerboseSensitivity == true)
			{
				if ( sensitivity_fast->Integration() == OPENSMOKE_SENSITIVITYANALYSIS_FAST_UNSTEADY0D_EULER_IMPLICIT || 
					 sensitivity_fast->Integration() == OPENSMOKE_SENSITIVITYANALYSIS_FAST_UNSTEADY0D_EULER_EXPLICIT   )
				{
					double tStart,tEnd;
					BzzMatrix J;
					o.GetLastJacobian(&J);
					sensitivity_fast->Update(Tau, J, T, P, rho, Cp+v*v/T, v, x);
				}
				else if (sensitivity_fast->Integration() == OPENSMOKE_SENSITIVITYANALYSIS_FAST_UNSTEADY0D_ACCURATE)
				{
					sensitivity_fast->UpdateOutputFile(Tau,y);
				}
			}

			if (iVerboseElementFluxAnalysis == true)
			{
			}

            countFileSteps = 0;
        }

		if (iVerboseLocalElementFluxAnalysis == true)
		{
			if (index_local_ElementFluxAnalysis.Size()>0)
				if (Tau>=index_local_ElementFluxAnalysis[1])
				{
					BzzVector rForward(mix->NumberOfReactions());
					BzzVector rBackward(mix->NumberOfReactions());
					mix->ComputeFromConcentrations(T, c, c.GetSumElements(), rForward, rBackward);

					stringstream tag;
					tag << index_local_ElementFluxAnalysis[1];
					string name_file = outputFolderName + "/FluxAnalysis_" + tag.str();
					ElementFluxAnalysis.Run(name_file, rForward, rBackward);
					index_local_ElementFluxAnalysis.DeleteElement(1);
				}
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

void OpenSMOKE_PFR::LabelODEFile()
{
    int i;

    if (iVerbose == true)
    {
        {
			int fOutputCount = 1;
			PrintTagOnGnuplotLabel(20, fOutput, "Tau[s]",		fOutputCount);
			PrintTagOnGnuplotLabel(20, fOutput, "x[m]",			fOutputCount);
			PrintTagOnGnuplotLabel(20, fOutput, "T[K]",			fOutputCount);
			PrintTagOnGnuplotLabel(20, fOutput, "dummy",		fOutputCount);
			PrintTagOnGnuplotLabel(20, fOutput, "P[atm]",		fOutputCount);
			PrintTagOnGnuplotLabel(20, fOutput, "v[m/s]",		fOutputCount);
			PrintTagOnGnuplotLabel(20, fOutput, "QR[W/m3]",		fOutputCount);
			PrintTagOnGnuplotLabel(20, fOutput, "QE[W/m3]",		fOutputCount);
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
		if (mix->polimiSoot->IsSoot() == true)
        {
            fSoot	<< setw(20) << left << "Tau[s](1)"
					<< setw(20) << left << "x[m](2)"
					<< setw(20) << left << "T[K](3)";

//          mix->soot_manager.GnuPlotInterface(fSoot, 4);

            fSoot	<< setw(20) << left << "L_fv(4)"
					<< setw(20) << left << "L_x(5)"
					<< setw(20) << left << "L_y(6)"
					<< setw(20) << left << "L_rho[kg/m3](7)"
					<< setw(20) << left << "L_N[#/m3](8)"
					<< setw(20) << left << "L_H/C[-](9)"
					<< setw(20) << left << "L_d10[m](10)"
					<< setw(20) << left << "L_d32[m](11)"
					<< setw(20) << left << "L_dstd[m](12)";

            fSoot	<< setw(20) << left << "S_fv(13)"
					<< setw(20) << left << "S_x(14)"
					<< setw(20) << left << "S_y(15)"
					<< setw(20) << left << "S_rho[kg/m3](16)"
					<< setw(20) << left << "S_N[#/m3](17)";
            fSoot << endl << endl;
        }

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
			PrintTagOnGnuplotLabel(20, fEnergy, "QE[W/m3]",		fOutputCount);
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

		if (iVerboseSelectivity == true)
		{
            fSelectivity	<< setw(20) << left << "Tau[s](1)"
							<< setw(20) << left << "x[m](2)"
							<< setw(20) << left << "T[K](3)"
							<< setw(20) << left << "eta(4)";
			
			int fOutputCount = 5;
			for (i=1;i<=NC;i++)
				PrintTagOnGnuplotLabel(20, fSelectivity, mix->names[i]+"_w",	fOutputCount);
			for (i=1;i<=NC;i++)
				PrintTagOnGnuplotLabel(20, fSelectivity, mix->names[i]+"_x",	fOutputCount);
			for (i=1;i<=NC;i++)
				PrintTagOnGnuplotLabel(20, fSelectivity, mix->names[i]+"_xC",	fOutputCount);

			fSelectivity << endl;
			fSelectivity << endl;
		}
    }
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
//																								   //
//									SOLVING FUNCTIONS											   //
//																								   //
/////////////////////////////////////////////////////////////////////////////////////////////////////

void OpenSMOKE_PFR::Initialize()
{
	if (inletStream->iUndefinedFlowRate == true)
	{
		if (iTimeIndipendent == false)
			ErrorMessage("If reactor length is specified, mass flow rate must be specified for the inlet stream...");

		if (iMomentum == true)
			ErrorMessage("If momentum equation is enabled, mass flow rate must be specified for the inlet stream...");

		if (iUserDefinedHeatFlux != NONE)
			ErrorMessage("If heat flux option is enabled, mass flow rate must be specified for the inlet stream...");

		if (iUserDefinedExchangeArea != NONE)
			ErrorMessage("If exchange area option is enabled, mass flow rate must be specified for the inlet stream...");

		if (iUserDefinedHeatExchangeCoefficient != NONE)
			ErrorMessage("If heat exchange coefficient option is enabled, mass flow rate must be specified for the inlet stream...");

		if (iUserDefinedAmbientTemperature != NONE)
			ErrorMessage("If ambient temperature option is enabled, mass flow rate must be specified for the inlet stream...");
	}

    T            = inletStream->T;
    P            = inletStream->P;
    omega        = inletStream->omega;
    massFlowRate = inletStream->massFlowRate;

    if (iUserDefinedTemperature == true)
        ud_temperature_profile.Check(0., T);

    if (iEnergy == false)
    {
        UpdateProperties_isothermal(0);
        indexProperties = 1;
    }
    else if (iEnergy == true)
        UpdateProperties(MINUSONE, NC+1);

    UpdateVelocity();

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

void OpenSMOKE_PFR::ReSolve()
{
	ODE_PFR_Object.assignPFR(this, iEnergy, iMomentum);
	Initialize();
	Solve();
}

void OpenSMOKE_PFR::DefineFromFile(const string inputFile)
{
    double  double_value;
    string  string_value;
    int     int_value;
	vector<string> string_vector;

    OpenSMOKE_Dictionary_PFR dictionary;
    dictionary.ParseFile(inputFile);


    // COMPULSORY: Energy Equation solution
    AssignEnergy(dictionary.Return("#Energy"));

    // SEMI-COMPULSORY: Reactor Length || Reactor Residence Time
    if (dictionary.Return("#Length", double_value, string_value))
    {
		iTimeIndipendent = false;
		AssignEnd(string_value, double_value);
	}	
	else if (dictionary.Return("#ResidenceTime", double_value, string_value))
    {
		iTimeIndipendent = true;
		AssignEnd(string_value, double_value);
	}

    // COMPULSORY: Reactor Diameter or Area
    if (dictionary.Return("#Diameter", double_value, string_value))
        AssignDiameter(string_value, double_value);
    else if (dictionary.Return("#Area", double_value, string_value))
        AssignArea(string_value, double_value);

    // OPTIONAL: Momentum Equation
    if (dictionary.Return("#Momentum"))
        SetMomentumEquation();

	// OPTIONAL: KeySpecies
    if (dictionary.Return("#Key", string_vector))
        SetKeySpecies(string_vector);

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

    // OPTIONAL: UserDefinedTemperature
    if (dictionary.Return("#UserDefined_T", string_value))
        SetUserDefinedTemperature(string_value);

	// OPTIONAL: Verbose Energy
    if (dictionary.Return("#VerboseEnergy"))
        SetVerboseEnergy();

	// OPTIONAL: Specific heat flux
    if (dictionary.Return("#Qe", double_value, string_value))
        SetConstantHeatFlux(double_value, string_value);

	// OPTIONAL: Specific exchange area
    if (dictionary.Return("#Ae", double_value, string_value))
        SetConstantSpecificExchangeArea(double_value, string_value);

	// OPTIONAL: Exchange coefficient
    if (dictionary.Return("#U", double_value, string_value))
        SetConstantHeatExchangeCoefficient(double_value, string_value);

	// OPTIONAL: Ambient temperature
    if (dictionary.Return("#Tambient", double_value, string_value))
        SetConstantAmbientTemperature(double_value, string_value);


	// OPTIONAL: UserDefined_Qe
    if (dictionary.Return("#UserDefined_Qe", string_value))
        SetUserDefinedHeatFlux(string_value);

	// OPTIONAL: UserDefined_Ae
    if (dictionary.Return("#UserDefined_Ae", string_value))
        SetUserDefinedSpecificExchangeArea(string_value);

	// OPTIONAL: UserDefined_U
    if (dictionary.Return("#UserDefined_U", string_value))
        SetUserDefinedHeatExchangeCoefficient(string_value);

	// OPTIONAL: UserDefined_Tambient
    if (dictionary.Return("#UserDefined_Tambient", string_value))
        SetUserDefinedAmbientTemperature(string_value);
 
	// OPTIONAL: Inlet Viscosity
    if (dictionary.Return("#Viscosity", double_value, string_value))
        SetViscosity(double_value, string_value);


	// OPTIONAL: UserDefined_Diameter
    if (dictionary.Return("#UserDefined_Diameter", string_value))
        SetUserDefinedDiameter(string_value);

	// OPTIONAL: UserDefined_Area
    if (dictionary.Return("#UserDefined_Area", string_value))
        SetUserDefinedArea(string_value);

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
    if (dictionary.Return("#Selectivity"))
		SetSelectivitiesOnFile();

	// OPTIONAL: Rate of Production Analysis
    if (dictionary.Return("#VerboseROPA", string_vector))
		SetVerboseROPAOnFile(string_vector);

	// OPTIONAL: Sensitivity Analysis
    if (dictionary.Return("#SensitivityOptions", string_vector))
		SetSensitivityOptions(string_vector);

	// OPTIONAL: Sensitivity Analysis
    if (dictionary.Return("#Sensitivity", string_vector))
		SetSensitivityOnFile(string_vector);

	// OPTIONAL: Rate of Production Analysis
    if (dictionary.Return("#ElementFluxAnalysis", string_vector))
		SetElementFluxAnalysisOnFile(string_vector);

	// OPTIONAL: Local Rate of Production Analysis
    if (dictionary.Return("#LocalElementFluxAnalysis", string_vector))
		SetLocalElementFluxAnalysisOnFile(string_vector);

	// OPTIONAL: Experiment
    if (dictionary.Return("#Experiment", string_vector))
	{
		SetHistory();
		SetExperimentOnFile(string_vector);
	}

	// OPTIONAL: LogExperiment
	if (dictionary.Return("#LogExperiment", int_value))
		SetLogExperiment(int_value);

	// OPTIONAL: OutputSpecies
	if (dictionary.Return("#OutputSpecies", string_vector))
		SetOutputSpecies(string_vector);

	// OPTIONAL: ODESolver
    if (dictionary.Return("#ODESolver", string_value))
        SetODESolver(string_value);

	Lock();
}

/*
#include "C:\Politecnico\Temp\RadauInterface\RadauInterface\OpenSMOKE_Interface_RADAU.h"
#include "C:\Politecnico\Temp\RadauInterface\RadauInterface\OpenSMOKE_DLSODE.h"

OpenSMOKE_PFR *ptPFRForExternalODESolvers;
BzzVector yBzz;
BzzVector dyBzz;

void odesytemisothermal_radau(int *n, double *x, double *y, double *fy, double *rpar, int *ipar)
{
	double etaBzz = x[0];

	for(int i=0;i<n[0];i++)
		yBzz[i+1] = y[i];

	ptPFRForExternalODESolvers->ODESystem_Isothermal_PFR(yBzz, etaBzz, dyBzz);

	for(int i=0;i<n[0];i++)
		fy[i]=dyBzz[i+1];
}

void dummy_solout(int *nr,double *xold,double *x,double *y,double *cont,int *lrc,int *n,double *rpar,int *ipar,int *irtrn)
{
	double etaBzz = x[0];
	for(int i=0;i<n[0];i++)
		yBzz[i+1] = y[i];
	ptPFRForExternalODESolvers->ODEPrint(yBzz, etaBzz);
}

void jacobian_radau(int *n, double *x, double *y, double *dfy, int *ldfy, double *rpar, double *ipar) { }
void mass_dummy(int *n,double *am, int *lmas,int *rpar, int *ipar) { }

void odesytemisothermal_dlsode(int *n, double *x, double *y, double *fy)
{
	double etaBzz = x[0];

	for(int i=0;i<n[0];i++)
		yBzz[i+1] = y[i];

	ptPFRForExternalODESolvers->ODESystem_Isothermal_PFR(yBzz, etaBzz, dyBzz);

	for(int i=0;i<n[0];i++)
		fy[i]=dyBzz[i+1];
}

void jacobian_dlsode(int *n, double *x, double *y, int *ml, int *mu, double *pd, int *nrowpd) {}
*/


void OpenSMOKE_PFR::Solve()
{
    double timeStart, timeEnd;
    BzzVector xMin;
    BzzVector xMax;
    BzzVector xInitial;

	ODE_PFR_Object.assignPFR(this, iEnergy, iMomentum);
	PrepareFiles();

    // 1.A Isothermal PFR
    if (iEnergy == false && iMomentum == false && iTwoEquationModel == false)
    {
		if ( ( iVerboseSensitivity == false ) || 
			 ( iVerboseSensitivity == true && sensitivity_fast->Integration() != OPENSMOKE_SENSITIVITYANALYSIS_FAST_UNSTEADY0D_ACCURATE) )
		{
			ChangeDimensions(NC+2, &xMin);  xMin=ZERO;					
			ChangeDimensions(NC+2, &xMax);  xMax=ONE;					// Mole fractions
				                            xMax[NC+1] = MAXNUMBER;		// Mass Flow Rate
					                        xMax[NC+2] = MAXNUMBER;		// Eta
	        xInitial = omega;
		    xInitial.Append(massFlowRate);
			xInitial.Append(0.);
		}
		else
		{
			int NSENSITIVITY = sensitivity_fast->NumberOfCoefficients();
			int NEQS = NC+2+NSENSITIVITY;
			ChangeDimensions(NEQS, &xMin);  xMin=-MAXNUMBER;					
			ChangeDimensions(NEQS, &xMax);  xMax= MAXNUMBER;					
	        
			for(int i=1;i<=NC;i++)	xMin[i] = ZERO;
			for(int i=1;i<=NC;i++)	xMax[i] = ONE;
			xMin[NC+1] = ZERO;
			xMin[NC+2] = ZERO;

			ChangeDimensions(NSENSITIVITY, &S_right);
			ChangeDimensions(NEQS, &xInitial);
			for(int i=1;i<=NC;i++)	xInitial[i] = omega[i];
			xInitial[NC+1] = massFlowRate;
			xInitial[NC+2] = 0.;
		}
	}

    // 1.B Isothermal PFR + Momentum
    if (iEnergy == false && iMomentum == true && iTwoEquationModel == false)
    {
        ChangeDimensions(NC+3, &xMin);
		xMin=ZERO;
		xMin[NC+3] = MINPRESSURE;			// Pressure

        ChangeDimensions(NC+3, &xMax);
        xMax=ONE;							// Mass fractions
        xMax[NC+1] = MAXNUMBER;				// Mass Flow Rate
        xMax[NC+2] = MAXNUMBER;				// Eta
        xMax[NC+3] = MAXNUMBER;				// Pressure

        xInitial = omega;
        xInitial.Append(massFlowRate);
        xInitial.Append(0.);
        xInitial.Append(P);
    }


    // 1.C Isothermal PFR	+ Two Equation Model
    if (iEnergy == false && iMomentum == false && iTwoEquationModel == true)
    {
        soot2EModel->initial_values(rho);
        indexTwoEquations = NC+3;

        ChangeDimensions(NC+4, &xMin);
        xMin=ZERO;
        
		ChangeDimensions(NC+4, &xMax);
        xMax=ONE;
        xMax[NC+1] = MAXNUMBER;
        xMax[NC+2] = MAXNUMBER;
        xMax[NC+3] = ONE;
        xMax[NC+4] = ONE;

        xInitial = omega;
        xInitial.Append(massFlowRate);
        xInitial.Append(0.);
        xInitial.Append(soot2EModel->phiNStart);
        xInitial.Append(soot2EModel->phiMStart);
    }

    // 1.D Isothermal PFR	+ Momentum + Two Equation Model
    if (iEnergy == false && iMomentum == true && iTwoEquationModel == true)
    {
        soot2EModel->initial_values(rho);
        indexTwoEquations = NC+4;

        ChangeDimensions(NC+5, &xMin);
        xMin=ZERO;
		xMin[NC+3] = MINPRESSURE;			// Pressure
        
		ChangeDimensions(NC+5, &xMax);
        xMax=ONE;
        xMax[NC+1] = MAXNUMBER;
        xMax[NC+2] = MAXNUMBER;
        xMax[NC+3] = MAXNUMBER;
        xMax[NC+4] = ONE;
        xMax[NC+5] = ONE;

        xInitial = omega;           				// Species
        xInitial.Append(massFlowRate);        		// Mass flow rate
        xInitial.Append(0.);						// Residence time
        xInitial.Append(P);				          	// Pressure
        xInitial.Append(soot2EModel->phiNStart);	// Soot particle number density
        xInitial.Append(soot2EModel->phiMStart);	// Soot mass fraction
    }


    // 2.A NonIsothermal PFR
    if (iEnergy == true && iMomentum == false && iTwoEquationModel == false)
    {
        ChangeDimensions(NC+3, &xMin);
        xMin=ZERO;
        xMin[NC+1] = 200.;
        ChangeDimensions(NC+3, &xMax);
        xMax=ONE;
        xMax[NC+1] = 6000.;
        xMax[NC+2] = MAXNUMBER;
        xMax[NC+3] = MAXNUMBER;

        xInitial = omega;
        xInitial.Append(T);
        xInitial.Append(massFlowRate);
        xInitial.Append(0.);
    }

    // 2.B NonIsothermal PFR + Momentum
    if (iEnergy == true && iMomentum == true && iTwoEquationModel == false)
    {
        ChangeDimensions(NC+4, &xMin);
        xMin=ZERO;
        xMin[NC+1] = 200.;
        xMin[NC+4] = MINPRESSURE;			// Pressure

		ChangeDimensions(NC+4, &xMax);
        xMax=ONE;
        xMax[NC+1] = 6000.;
        xMax[NC+2] = MAXNUMBER;
        xMax[NC+3] = MAXNUMBER;
        xMax[NC+4] = MAXNUMBER;

        xInitial = omega;
        xInitial.Append(T);
        xInitial.Append(massFlowRate);
        xInitial.Append(0.);
        xInitial.Append(P);
    }

    // 2.C NonIsothermal PFR + TwoEquations
    if (iEnergy == true && iMomentum == false && iTwoEquationModel == true)
    {
        soot2EModel->initial_values(rho);
        indexTwoEquations = NC+4;

        ChangeDimensions(NC+5, &xMin);
        xMin=ZERO;
        xMin[NC+1] = 200.;
        ChangeDimensions(NC+5, &xMax);
        xMax=ONE;
        xMax[NC+1] = 6000.;
        xMax[NC+2] = MAXNUMBER;
        xMax[NC+3] = MAXNUMBER;
        xMax[NC+4] = ONE;
        xMax[NC+5] = ONE;

        xInitial = omega;
        xInitial.Append(T);
        xInitial.Append(massFlowRate);
        xInitial.Append(0.);
        xInitial.Append(soot2EModel->phiNStart);
        xInitial.Append(soot2EModel->phiMStart);
    }

    // 2.D NonIsothermal PFR + Momentum + TwoEquations
    if (iEnergy == true && iMomentum == true && iTwoEquationModel == true)
    {
        soot2EModel->initial_values(rho);
        indexTwoEquations = NC+5;

        ChangeDimensions(NC+6, &xMin);
        xMin=ZERO;
        xMin[NC+1] = 200.;
        xMin[NC+4] = MINPRESSURE;			// Pressure

		ChangeDimensions(NC+6, &xMax);
        xMax=ONE;
        xMax[NC+1] = 6000.;
        xMax[NC+2] = MAXNUMBER;
        xMax[NC+3] = MAXNUMBER;
        xMax[NC+4] = MAXNUMBER;
        xMax[NC+5] = ONE;
        xMax[NC+6] = ONE;

        xInitial = omega;           				    // Species
        xInitial.Append(T);			                    // Temperature
        xInitial.Append(massFlowRate);	                // Mass flow rate

		if (iTimeIndipendent == true)
			xInitial.Append(initialLenght);		    		// Initial Length
		else
			xInitial.Append(initialTime);		    		// Initial Time

        xInitial.Append(P);					            // Pressure
        xInitial.Append(soot2EModel->phiNStart);		// Soot particle number density
        xInitial.Append(soot2EModel->phiMStart);		// Soot mass fraction
    }

	if (iTimeIndipendent == true)	
		o(xInitial, initialTime, &ODE_PFR_Object);
	else
		o(xInitial, initialLenght, &ODE_PFR_Object);

//    o.StepPrint(ODEPrintExternal_PFR);
    o.MyStepPrint();
    o.SetMinimumConstraints(xMin);
    o.SetMaximumConstraints(xMax);
//	o.SetMaxStep(MAX_TIME_STEPS);

    if (iRelativeTolerance == true)	o.SetTollRel(relativeTolerance);
    if (iAbsoluteTolerance == true)	o.SetTollAbs(absoluteTolerance);

//	o.SetHMax(2.e-6);

	if ( iVerboseSensitivity == true )
	{
		if (iTimeIndipendent == false)
			ErrorMessage("The sensitivity analysis can be used only if the contact time is the independent variable!");
	}

	// Additional instructions
	inletStream->SetVelocity(v, "m/s");
	
	// BzzODE
	if (iODE_Solver == ODE_SOLVER_BZZ)
    {
        countGlobalIterations = -1;  // From -1 to avoid to store results from the first iteration

        timeStart = BzzGetCpuTime();

		if (iTimeIndipendent == true)	o(TauTotal+initialTime,TauTotal+initialTime);
		else							o(L+initialLenght, L+initialLenght);

		status = o.GetCalculationState();
		timeEnd = BzzGetCpuTime();

		if (iVerbose == true)
		{
			cout << endl;
			cout << "Number of function for Jacobian: " << o.GetNumFunctionForJacobian()		<< endl;
			cout << "Numerical Jacobians: "				<< o.GetNumNumericalJacobian()			<< endl;
			cout << "Time DAE solution: "				<< timeEnd - timeStart	<< " s"						<< endl << endl;
		}

		cout << "BzzODE CPU time: " << timeEnd - timeStart << endl;
		cout << "Press enter to continue..." << endl;
		getchar();

    }
	
	
	// RADAU with Bzz interface
/*	else if (iODE_Solver == ODE_SOLVER_RADAU)
	{
		countGlobalIterations = -1;  // From -1 to avoid to store results from the first iteration

		ptPFRForExternalODESolvers = this;
		ChangeDimensions(xInitial.Size(), &yBzz);
		ChangeDimensions(xInitial.Size(), &dyBzz);

		OpenSMOKE_Interface_RADAU osmradau;
		osmradau.SetSystem(odesytemisothermal_radau);
		osmradau.SetJacobian(jacobian_radau);
		osmradau.SetMassMatrix(mass_dummy);
		osmradau.SetOutputFunction(dummy_solout);
	    
		if (iRelativeTolerance == true)	osmradau.SetRelativeTolerance(relativeTolerance);
		if (iAbsoluteTolerance == true)	osmradau.SetAbsoluteTolerance(absoluteTolerance);

		osmradau.SetAnalyticalJacobian(false);
		osmradau.SetInitialValues(0., xInitial);
		
		double startTime = BzzGetCpuTime();
		if (iTimeIndipendent == true)	osmradau.Solve(TauTotal);
		else							osmradau.Solve(L);
		double endTime = BzzGetCpuTime();
		osmradau.Status();
		cout << "RADAU CPU time: " << endTime - startTime << endl;
		cout << "Press enter to continue..." << endl;
		getchar();
	}
	
	else if (iODE_Solver == ODE_SOLVER_DLSODE)
	{
		countGlobalIterations = -1;  // From -1 to avoid to store results from the first iteration

		ptPFRForExternalODESolvers = this;
		ChangeDimensions(xInitial.Size(), &yBzz);
		ChangeDimensions(xInitial.Size(), &dyBzz);

		OpenSMOKE_DLSODE osm;
		osm.SetSystem(odesytemisothermal_dlsode);
		osm.SetJacobian(jacobian_dlsode);
	    
		if (iRelativeTolerance == true)	osm.SetRelativeTolerance(relativeTolerance);
		if (iAbsoluteTolerance == true)	osm.SetAbsoluteTolerance(absoluteTolerance);

		osm.SetAnalyticalJacobian(false);
		osm.SetInitialValues(0., xInitial);

		double alfa = 5.;
		int nOutIntervals = 15;
		double sum = 1.;
		for (int i=1;i<=nOutIntervals-1;i++)
			sum += pow(alfa,i);
		
		double startTime = BzzGetCpuTime();
		if (iTimeIndipendent == true)
		{
			double dx = TauTotal/sum;
			double xend = dx;
			for (int i=1;i<=nOutIntervals;i++)
			{
				osm.Solve(xend);
				dx*=alfa;
				xend += dx;

			}
		}
		else
		{
			double dx = L/sum;
			double xend = dx;
			for (int i=1;i<=nOutIntervals;i++)
			{
				osm.Solve(xend);
				dx*=alfa;
				xend += dx;

			}
		}

		double endTime = BzzGetCpuTime();
		osm.Status();
		cout << "DLSODE CPU time: " << endTime - startTime << endl;
		cout << "Press enter to continue..." << endl;
		getchar();
	}
	*/
	if (iHistory == true)
	{
		SaveOnBinaryFile(outputOSMName);
	}

	CloseFiles();

	if (iVerboseElementFluxAnalysis == true)
	{
	}
}


/////////////////////////////////////////////////////////////////////////////////////////////////////
//																								   //
//									MOMENTUM EQUATIONS											   //
//																								   //
/////////////////////////////////////////////////////////////////////////////////////////////////////

double OpenSMOKE_PFR::FrictionFactor()
{
    double mu = SimplifiedViscosity();
    double Re = rho*v*D / mu;

    if (Re<=2100.)	return 16./Re;
    else			return 0.0791/pow(Re, 0.25);
}

double OpenSMOKE_PFR::FrictionForce()
{
    return  4./D * 0.50*rho*v*v * FrictionFactor();
}

double OpenSMOKE_PFR::Delta()
{
    double delta = 0.;
    for (int i=1; i<=NC; i++)
        delta += mix->uM[i]*domega[i];
    delta *= MWtot;

	return delta;
}

void OpenSMOKE_PFR::AdditionalExtraTerms()
{
    delta	= Delta();
	
	if (iMomentum == true)
	{
		FForce	= FrictionForce();
		Beta = (1.-rho*v*v/P);

		T_coefficient   = Cp + v*v/Beta/T;
		T_extra_terms   = -v*v/Beta*(delta + FForce/P);
	}
	else
	{
		Beta = 1.0;

		T_coefficient   = Cp + v*v/Beta/T;
		T_extra_terms   = -v*v/Beta*delta;
	}
}

double OpenSMOKE_PFR::dv()
{
	double derivative;

	if (iMomentum == true)
		derivative = v/Beta*delta + FForce/P;
	else
		derivative = v/Beta*delta;
	
	if (iEnergy == true) derivative   += v/Beta * 1./T*dT;

	return derivative;
}

void OpenSMOKE_PFR::UpdateVelocity()
{
    v = massFlowRate / (Area*rho);
}

void OpenSMOKE_PFR::CheckPressure()
{
	if (P <= MINPRESSURE)
		ErrorMessage("The pressure is too low. Please check the reactor length.");
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
//																								   //
//									UPDATING FUNCTIONS											   //
//																								   //
/////////////////////////////////////////////////////////////////////////////////////////////////////

void OpenSMOKE_PFR::UpdateHeatFlux(const double tau, const double csi)
{
	if (iUserDefinedHeatExchangeCoefficient == USERDEFINED)		U			= ud_U_profile.GiveMeValue(tau, csi);
	if (iUserDefinedAmbientTemperature		== USERDEFINED)		Tambient	= ud_Tambient_profile.GiveMeValue(tau, csi);
	if (iUserDefinedHeatFlux				== USERDEFINED)		Qe			= ud_Qe_profile.GiveMeValue(tau, csi);
	
	if (iUserDefinedHeatExchangeCoefficient != NONE)			Qe = U*(T-Tambient);
}

void OpenSMOKE_PFR::UpdateExchangeArea(const double tau, const double csi)
{
	if (iUserDefinedExchangeArea == NONE)			Ae = 4./D;
	if (iUserDefinedExchangeArea == USERDEFINED)	Ae = ud_Ae_profile.GiveMeValue(tau, csi);
}

void OpenSMOKE_PFR::UpdateProperties_isothermal(int memoIndex)
{
    double cTot;

    mix->GetMWAndMoleFractionsFromMassFractions(MWtot, x, omega);

    // --------------------------------------------------------------------------
    // PROPERTIES FOR CONSTANT T + MEMORIZATION
    // --------------------------------------------------------------------------
    if (memoIndex==0)
    {
        // Specific Heat [J/kgK]
        // ----------------------------------------------------------------------
        //mix->SpeciesCp(T);

        // Kinetic parameters (the enthalpies are calculated in this step)
        // ----------------------------------------------------------------------
        mix->ComputeKineticParameters( T, log(T), 1./T, P);

        memoIndex = 1;
    }

    // --------------------------------------------------------------------------
    // PROPERTIES FOR CONSTANT T (Recover from tables)
    // --------------------------------------------------------------------------
    if (memoIndex == 1)
    {
        // a. Calcolo della concentrazione e della densita [kmol/m3] [kg/m3]
        cTot    = P  / (Constants::R_J_kmol*T);
        rho     = cTot * MWtot;
        c       = cTot*x;

        // c. Reaction Rates
        mix->ComputeFromConcentrations( T, c, cTot, &R);		// [kmol/m3/s]
        ElementByElementProduct(R, mix->M, &R);					// [kg/m3/s]
    }

    // Global Kinetics
    if (iGlobalKinetics==true)
        global->GiveMeFormationRates(T,c,R);
}

void OpenSMOKE_PFR::UpdateProperties(int jacobianIndex, int indexT)
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
		// Calcolo della concentrazione e della densita [kmol/m3] [kg/m3]
		// ----------------------------------------------------------------------
		cTot = P  / (Constants::R_J_kmol*T);
		rho = cTot * MWtot;
		c = cTot*x;

		// Specific Heat [J/kgK]
		// ----------------------------------------------------------------------
		mix->SpeciesCp(T);
		Cp = mix->MixCp_FromMassFractions(omega);

		// Table Creation (Memorization)
		// ----------------------------------------------------------------------
		mix->ComputeKineticParameters( T, log(T), 1./T, P);
		mix->ComputeFromConcentrations( T, c, cTot, &R);		// [kmol/m3/s]
		ElementByElementProduct(R, mix->M, &R);					// [kg/m3/s]
		QReaction = - mix->ComputeQReaction(T);					// [J/m3.s]
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
		// Calcolo della concentrazione e della densita [kmol/m3] [kg/m3]
		// ----------------------------------------------------------------------
		cTot = P  / (Constants::R_J_kmol*T);
		rho = cTot * MWtot;
		c = cTot*x;

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

void OpenSMOKE_PFR::UpdateTwoEquationModel(BzzVector &y, BzzVector &dy)
{
    soot2EModel->update(T, P/101325., rho, 0., x, y[indexTwoEquations], y[indexTwoEquations+1]);
    soot2EModel->formation_rates();

    for (int i=1;i<=NC-1;i++)
		domega[i] += soot2EModel->SGas[i] / (rho*v);	// Gas species
	domega[NC] += soot2EModel->S / (rho*v);				// Soot

    dy[indexTwoEquations]	= soot2EModel->s / (rho*v);
    dy[indexTwoEquations+1] = soot2EModel->S / (rho*v);
}

void OpenSMOKE_PFR::UpdateReactionRates(const double TT, const double PP, BzzVector &xx, BzzVector &rr)
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
//						ODE SYSTEM EQUATIONS - ISOTHERMAL REACTOR								   //
//																								   //
/////////////////////////////////////////////////////////////////////////////////////////////////////

void OpenSMOKE_PFR::ODESystem_Isothermal_PFR(BzzVector &y, double eta, BzzVector &dy)
{
    int i;

    for (i=1;i<=NC;i++)	omega[i]    = y[i];
    massFlowRate	                = y[NC+1];
	/*
	double sum;
	for(int i=1;i<=NC;i++)
	{
		if (y[i]<0.)	omega[i]=1e-16;
		sum+=omega[i];
	}
	omega /= sum;*/

	if		(iTimeIndipendent == true)	{	Csi = y[NC+2]; Tau = eta;	}
	else								{	Tau = y[NC+2]; Csi = eta;	}

    // User Defined Temperature profile
    if (iUserDefinedTemperature == true)
    {
        T = ud_temperature_profile.GiveMeValue(Tau, Csi);
        indexProperties = 0;
    }

    // Updating geometry
    geometry.Update(Csi, D, Area);

    // Properties + Reaction rates (constant temperature)
    UpdateProperties_isothermal(indexProperties);

    // Velocity in PFR
    UpdateVelocity();

    // Conservation equations for all the species (PFR)
    domega = R / (rho*v);

    // Mass conservation
    dmassFlowRate = 0.;

    // Residence time
	if		(iTimeIndipendent == true)	dEta = 1.;
	else								dEta = 1./v;
    
    // Soot Two Equation Model
    if (iTwoEquationModel == true)
        UpdateTwoEquationModel(y, dy);

    // Recovering residuals (Soot 2E are already updated)
    for (i=1;i<=NC;i++)	dy[i]   = domega[i];
    dy[NC+1]                    = dmassFlowRate;
    dy[NC+2]                    = dEta;

	if	(iTimeIndipendent == true)	dy *= v;

	// TODO Sensitivity
	if ( iVerboseSensitivity == true)
		if (sensitivity_fast->Integration() == OPENSMOKE_SENSITIVITYANALYSIS_FAST_UNSTEADY0D_ACCURATE)
		{
			BzzMatrix J;
			o.GetLastJacobian(&J);
			sensitivity_fast->Update(J,T,P,rho,Cp,v,x, y, S_right);
			for (i=1;i<=S_right.Size();i++)
				dy[NC+2+i] = S_right[i];
		}
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
//																								   //
//						ODE SYSTEM EQUATIONS - ISOTHERMAL REACTOR - MOMENTUM					   //
//																								   //
/////////////////////////////////////////////////////////////////////////////////////////////////////

void OpenSMOKE_PFR::ODESystem_IsothermalMomentum_PFR(BzzVector &y, double eta, BzzVector &dy)
{
    int i;

    // Recovering variables
    for (i=1;i<=NC;i++)	omega[i]    = y[i];
    massFlowRate	                = y[NC+1];
	if		(iTimeIndipendent == true)	{	Csi = y[NC+2]; Tau = eta;	}
	else								{	Tau = y[NC+2]; Csi = eta;	}
    P					            = y[NC+3];

	// CheckPressure
	CheckPressure();

	// User Defined Temperature profile
    if (iUserDefinedTemperature == true)
    {
        T = ud_temperature_profile.GiveMeValue(Tau, Csi);
        indexProperties = 0;
    }

    // Updating geometry
    geometry.Update(Csi, D, Area);

    // Properties + Reaction rates (constant temperature)
    UpdateProperties_isothermal(indexProperties);

    // Velocity in PFR
    UpdateVelocity();

    // Mass conservation
    dmassFlowRate = 0.;

    // Conservation equations for all the species (PFR)
    domega = R / (rho*v);

    // Residence time
    if		(iTimeIndipendent == true)	dEta = 1.;
	else								dEta = 1./v;

    // Momentum equation
	AdditionalExtraTerms();
    dP = - FForce - rho*v*dv();

    // Soot Two Equation Model
    if (iTwoEquationModel == true)
        UpdateTwoEquationModel(y, dy);

    // Recovering residuals (Soot 2E are already updated)
    for (i=1;i<=NC;i++)	dy[i]   = domega[i];
    dy[NC+1]                    = dmassFlowRate;
    dy[NC+2]                    = dEta;
    dy[NC+3]                    = dP;

	if		(iTimeIndipendent == true)	dy *= v;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
//																								   //
//						ODE SYSTEM EQUATIONS - NON ISOTHERMAL REACTOR								   //
//																								   //
/////////////////////////////////////////////////////////////////////////////////////////////////////

void OpenSMOKE_PFR::ODESystem_NonIsothermal_PFR(BzzVector &y, double eta, BzzVector &dy)
{
    int i;

    // Recovering variables;
    for (i=1;i<=NC;i++)	omega[i]    = y[i];
    T				                = y[NC+1];
    massFlowRate	                = y[NC+2];
    if		(iTimeIndipendent == true)	{	Csi = y[NC+3]; Tau = eta;	}
	else								{	Tau = y[NC+3]; Csi = eta;	}

    // Updating geometry
	geometry.Update(Csi, D, Area);

	// Updating Exchange area
	UpdateExchangeArea(Tau, Csi);
    UpdateHeatFlux(Tau, Csi);

    // Evaluation of properties
    UpdateProperties(ODE_PFR_Object.jacobianIndex, NC+1);

    // Velocity in PFR
    UpdateVelocity();

    // Mass flow rate (velocity)
    dmassFlowRate = 0.;

    // Conservation equations for all the species (PFR)
    domega	= R / (rho*v);

	// Residence time
    if		(iTimeIndipendent == true)	dEta = 1.;
	else								dEta = 1./v;

    // Conservation equations for the energy (PFR)
    AdditionalExtraTerms();
    dT  = QReaction / (rho*v) -  Ae * Qe / (rho*v)  +  T_extra_terms;
    dT /= T_coefficient;
	// Attenzione: se la velocita' e' molto alta il termine v2 conta un casino
 
    // Soot Two Equation Model
    if (iTwoEquationModel == true)	
		UpdateTwoEquationModel(y, dy);

    // Recovering residuals (Soot 2E are already updated)
    for (i=1;i<=NC;i++)	dy[i] = domega[i];
    dy[NC+1] = dT;
    dy[NC+2] = dmassFlowRate;
    dy[NC+3] = dEta;

	if		(iTimeIndipendent == true)	dy *= v;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
//																								   //
//						ODE SYSTEM EQUATIONS - ADIABATIC REACTOR								   //
//																								   //
/////////////////////////////////////////////////////////////////////////////////////////////////////

void OpenSMOKE_PFR::ODESystem_NonIsothermalMomentum_PFR(BzzVector &y, double eta, BzzVector &dy)
{
    int i;

    // Recovering variables;
    for (i=1;i<=NC;i++)	omega[i]	= y[i];
    T								= y[NC+1];
    massFlowRate					= y[NC+2];
    if		(iTimeIndipendent == true)	{	Csi = y[NC+3]; Tau = eta;	}
	else								{	Tau = y[NC+3]; Csi = eta;	}
    P								= y[NC+4];

	// CheckPressure
	CheckPressure();

    // Updating geometry
    geometry.Update(Csi, D, Area);

    // Updating External source term
	UpdateExchangeArea(Tau, Csi);
    UpdateHeatFlux(Tau, Csi);	
	
    // Evaluation of properties
    UpdateProperties(ODE_PFR_Object.jacobianIndex, NC+1);

    // Velocity in PFR
    UpdateVelocity();

	// Mass flow rate (velocity)
    dmassFlowRate = 0.;

    // Residence time
    if		(iTimeIndipendent == true)	dEta = 1.;
	else								dEta = 1./v;

    // Conservation equations for all the species (PFR)
    domega	= R / (rho*v);

    // Conservation equations for the energy (PFR)
	AdditionalExtraTerms();
    dT  = QReaction / (rho*v)  -  Ae*Qe / (rho*v)  +  T_extra_terms;
    dT /= T_coefficient;

	// Momentum Equation
    dP = - FForce - rho*v*dv();

    // Soot Two Equation Model
    if (iTwoEquationModel == true)
        UpdateTwoEquationModel(y, dy);

    // Recovering residuals (Soot 2E are already updated)
    for (i=1;i<=NC;i++)	dy[i]   = domega[i];
    dy[NC+1]                    = dT;
    dy[NC+2]                    = dmassFlowRate;
    dy[NC+3]                    = dEta;
    dy[NC+4]                    = dP;

	if		(iTimeIndipendent == true)	dy *= v;
}


/////////////////////////////////////////////////////////////////////////////////////////////////////
//																								   //
//									DICTIONARY													   //
//																								   //
/////////////////////////////////////////////////////////////////////////////////////////////////////

OpenSMOKE_Dictionary_PFR::OpenSMOKE_Dictionary_PFR()
{
    SetupBase();
	SetName("OpenSMOKE_PFR Dictionary");

    Add("#Energy",        'O', 'N', "The energy equation is solved");
    Add("#NoEnergy",      'O', 'N', "The energy equation is not solved");
    
	Add("#Length",        'O', 'M', "Reactor length");
    Add("#ResidenceTime",   'O', 'M', "Reactor residence time");

    Add("#Diameter",		'O', 'M', "Reactor diameter");
    Add("#Area",			'O', 'M', "Reactor cross section area");
     
	Add("#Momentum",		'O', 'N', "The momentum equation is solved");
	Add("#VerboseEnergy",	'O', 'N', "Report on energy of PFR");

    Add("#Key",					'O', 'V', "Key species");

	Add("#Qe",				'O', 'M', "Specific heat flux");
	Add("#Ae",				'O', 'M', "Specific exchange area");
	Add("#U",				'O', 'M', "Heat transfer exchange coefficient");
	Add("#Tambient",		'O', 'M', "Ambient temperature");
	Add("#Viscosity",		'O', 'M', "Inlet viscosity");
	Add("#OutputSpecies",   'O', 'V', "Output species");

	// Optional
	Add("#ReactionRates",		'O', 'V', "Reaction rates on file (list of reaction indices)");
	Add("#FormationRates",		'O', 'V', "Formation rates on file (list of species names) ");
	Add("#ROPA",				'O', 'N', "Rate of Production Analysis");
	Add("#Selectivity",			'O', 'N', "Selectivity (towards the kkey species)");
	Add("#VerboseROPA",			'O', 'V', "Rate of Production Analysis (list of species names)");
	Add("#Sensitivity",			'O', 'V', "Sensitivity Analysis");
	Add("#SensitivityOptions",	'O', 'V', "Sensitivity Options");
	Add("#ElementFluxAnalysis",		 'O', 'V', "Element Flux Analysis (list of element names, lower case)");
	Add("#LocalElementFluxAnalysis", 'O', 'V', "Local Element Flux Analysis (list of times)");
	Add("#Experiment",			'O', 'V', "Experiment Analysis");
	Add("#LogExperiment",		'O', 'I', "Experiment Analysis: number of points");
	Add("#ODESolver",			'O', 'S', "ODE Solver: BZZ || RADAU || DLSODE");
	
	Add("#UserDefined_T",			'O', 'S', "The temperature is assigned from external file");
	Add("#UserDefined_Qe",			'O', 'S', "The specific heat flux is assigned from external file");
    Add("#UserDefined_Ae",			'O', 'S', "The specific exchange area is assigned from external file");
	Add("#UserDefined_U",			'O', 'S', "The heat exchange coefficient is assigned from external file");
    Add("#UserDefined_Tambient",	'O', 'S', "The ambient temperature is assigned from external file");
    Add("#UserDefined_Diameter",	'O', 'S', "The diameter is assigned from file");
    Add("#UserDefined_Area",		'O', 'S', "The cross section area is assigned from file");

	Add("#RelativeTolerance",		'O', 'D', "Relative tolerance (default: 1.2e-5)");
	Add("#AbsoluteTolerance",		'O', 'D', "Absolute tolerance (default: 1.0e-10)");
	
    Conflict("#Length",					"#ResidenceTime");
    Conflict("#Energy",					"#NoEnergy");
    Conflict("#Energy",					"#UserDefined_T");

	Conflict("#Diameter",				"#Area");
	Conflict("#UserDefined_Area",		"#Area");
	Conflict("#UserDefined_Diameter",	"#Diameter");
	Conflict("#UserDefined_Area",		"#Diameter");
	Conflict("#UserDefined_Diameter",	"#Area");
	Conflict("#UserDefined_Area",		"#UserDefined_Diameter");

 	Conflict("#UserDefined_Qe",			"#Qe");
	Conflict("#UserDefined_Ae",			"#Ae");
	Conflict("#UserDefined_U",			"#U");
	Conflict("#UserDefined_Tambient",	"#Tambient");

	Conflict("#NoEnergy",				"#Qe");
	Conflict("#NoEnergy",				"#Ae");
	Conflict("#NoEnergy",				"#U");
	Conflict("#NoEnergy",				"#Tambient");
	Conflict("#NoEnergy",				"#UserDefined_Qe");
	Conflict("#NoEnergy",				"#UserDefined_Ae");
	Conflict("#NoEnergy",				"#UserDefined_U");
	Conflict("#NoEnergy",				"#UserDefined_Tambient");

	Conflict("#Qe",						"#U");
	Conflict("#Qe",						"#Tambient");
	Conflict("#Qe",						"#UserDefined_U");
	Conflict("#Qe",						"#UserDefined_Tambient");
	Conflict("#UserDefined_Qe",			"#U");
	Conflict("#UserDefined_Qe",			"#Tambient");
	Conflict("#UserDefined_Qe",			"#UserDefined_U");
	Conflict("#UserDefined_Qe",			"#UserDefined_Tambient");


	Compulsory("#Energy",   "#NoEnergy");
	Compulsory("#Length",   "#ResidenceTime");
	Compulsory("#Diameter", "#Area", "#UserDefined_Area", "#UserDefined_Diameter");
	

    Lock();
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
//																								   //
//									ODE SYSTEM - CLASS DEFINITION								   //
//																								   //
/////////////////////////////////////////////////////////////////////////////////////////////////////

void MyOdeSystem_PFR::ObjectBzzPrint(void)
{
}

void MyOdeSystem_PFR::GetSystemFunctions(BzzVector &x, double eta, BzzVector &f)
{
    if (iEnergy == false && iMomentum == false )
        ptPFR->ODESystem_Isothermal_PFR(x, eta, f);

    else if (iEnergy == false && iMomentum == true )
        ptPFR->ODESystem_IsothermalMomentum_PFR(x, eta, f);

    else if (iEnergy == true && iMomentum == false )
        ptPFR->ODESystem_NonIsothermal_PFR(x, eta, f);

    else if (iEnergy == true && iMomentum == true )
        ptPFR->ODESystem_NonIsothermalMomentum_PFR(x, eta, f);
}

void MyOdeSystem_PFR::MyODEPrint(BzzVector &y, double eta)
{
	ptPFR->ODEPrint(y, eta);
}

void MyOdeSystem_PFR::assignPFR(OpenSMOKE_PFR *pfr, bool _iEnergy, bool _iMomentum)
{
    ptPFR       = pfr;
    iEnergy		= _iEnergy;
    iMomentum   = _iMomentum;
}

void OpenSMOKE_PFR::SaveOnBinaryFile(BzzSave &fOutput)
{
	string dummy;
	char name[Constants::NAME_SIZE];

	int nSteps = countGlobalIterations;			
	BzzVector taugrid = GetBzzVector(nSteps, 1, Tau_History);
	BzzVector csigrid = GetBzzVector(nSteps, 1, Csi_History);
	BzzVector tgrid   = GetBzzVector(nSteps, 1, T_History);
		
	dummy = "PFR";
	strcpy(name, dummy.c_str());
	fOutput.fileSave.write((char*) name, sizeof(name));

	dummy = "V20100417";
	strcpy(name, dummy.c_str());
	fOutput.fileSave.write((char*) name, sizeof(name));

	dummy = "NO-MOMENTUM";
	strcpy(name, dummy.c_str());
	fOutput.fileSave.write((char*) name, sizeof(name));

	// Contact time
	dummy = "TAU";
	strcpy(name, dummy.c_str());
	fOutput.fileSave.write((char*) name, sizeof(name));
	fOutput << taugrid;

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

void OpenSMOKE_PFR::SaveOnBinaryFile(const string filename)
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

		if ( sensitivity_fast->Integration() != OPENSMOKE_SENSITIVITYANALYSIS_FAST_UNSTEADY0D_ACCURATE )
		{
			if (iEnergy == false)
				sensitivity_fast->MoveFromFiles(outputFolderName + "\\" + "Sensitivity.bin", countGlobalIterations, fOutput, "pfr-isothermal");
			else
				sensitivity_fast->MoveFromFiles(outputFolderName + "\\" + "Sensitivity.bin", countGlobalIterations, fOutput, "pfr-non-isothermal");
		}
		else if (sensitivity_fast->Integration() == OPENSMOKE_SENSITIVITYANALYSIS_FAST_UNSTEADY0D_ACCURATE)
		{
		}
	}
	
	PrintEndOnBinaryFile(fOutput);
	fOutput.End();

	// Post processor test
	OpenSMOKE_PostProcessor post_processor;
	post_processor.ReadFromBinaryFile(filename);
}
