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
#include "idealreactors/icem/OpenSMOKE_ICEM.h"
#include "basic/OpenSMOKE_Conversions.h"

const double OpenSMOKE_ICEM::rad_to_degrees = 360./(2.*Constants::pi);
const double OpenSMOKE_ICEM::degrees_to_rad = (2.*Constants::pi)/360.;

OpenSMOKE_ICEM::OpenSMOKE_ICEM()
{
	class_name	= "OpenSMOKE_ICEM";
	out_name	= "ICEM.out";
	
	Setup();

	assignedEnd = true;
	assignedClearanceVolume	= false;
	assignedCompressionRatio = false;
	assignedArmRatio = false;
	assignedStartAngle	= false;
	assignedRotationRate = false;

	number_of_cycles = 1.;
	
	heat_transfer_model = ICEM_BASE;
	a_Nu = 0.035;
	b_Nu = 0.80;
	c_Nu = 0.;
}


void OpenSMOKE_ICEM::AssignEnd(const std::string units, const double value)
{
}

void OpenSMOKE_ICEM::AssignRotationRate(const std::string units, const double value)
{
	rotation_rate = OpenSMOKE_Conversions::conversion_angular_velocity(value, units);
    assignedRotationRate	= true;
}

void OpenSMOKE_ICEM::AssignClearanceVolume(const std::string units, const double value)
{
	volume_clearance = OpenSMOKE_Conversions::conversion_volume(value, units);
    assignedClearanceVolume	= true;
}

void OpenSMOKE_ICEM::AssignCompressionRatio(const double value)
{
	compression_ratio = value;
    assignedCompressionRatio = true;
}

void OpenSMOKE_ICEM::AssignArmRatio(const double value)
{
	arm_ratio = value;
    assignedArmRatio = true;
}

void OpenSMOKE_ICEM::AssignStartAngle(const std::string units, const double value)
{
	start_angle = OpenSMOKE_Conversions::conversion_angle(value, units);
    assignedStartAngle	= true;
}

void OpenSMOKE_ICEM::AssignDiameter(const std::string units, const double value)
{
	diameter = OpenSMOKE_Conversions::conversion_length(value, units);
	area_base = Constants::pi/4.*diameter*diameter;
    assignedDiameter	= true;
}

void OpenSMOKE_ICEM::SetNumberOfCycles(const double value)
{
	number_of_cycles = value; 
}

void OpenSMOKE_ICEM::SetConstantExchangeArea(const double value, const std::string units)
{
	iUserDefinedExchangeArea = CONSTANT;
	A = OpenSMOKE_Conversions::conversion_area(value, units);
}

void OpenSMOKE_ICEM::SetUserDefinedExchangeArea(const std::string fileName)
{
    iUserDefinedExchangeArea = USERDEFINED;
    ud_A_profile.AssignFromFile(fileName, "AREA");
	ud_A_profile.SetName(name_object + " - Exchange Area Profile");
}

void OpenSMOKE_ICEM::UnsetUserDefinedExchangeArea()
{
    iUserDefinedExchangeArea = NONE;
}

void OpenSMOKE_ICEM::SetAdiabatic()
{
	heat_transfer_model = ICEM_NONE;
	U = 0.;
	Nu = 0.;
	Re = 0.;
	Pr = 0.;
}

void OpenSMOKE_ICEM::Lock()
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

	Initialize();
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
//																								   //
//									GAS MIXTURE PROPERTIES										   //
//																								   //
/////////////////////////////////////////////////////////////////////////////////////////////////////

void OpenSMOKE_ICEM::UpdateProperties_isothermal(int memoIndex)
{
	ErrorMessage("Isothermal properties cannot be used...");
}

void OpenSMOKE_ICEM::UpdateProperties(const int jacobianIndex, const int indexT)
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
        // a. Calcolo della concentrazione e della densita [kmol/m3] [kg/m3]
		moles	= mass / MWtot;						// [kmol]
		rho		= mass / volume;					// [kg/m3]
		cTot	= rho  / MWtot;						// [kmol/m3]
		P		= cTot * (Constants::R_J_kmol*T);	// [Pa]
		c       = cTot * x;							// [kmol/m3]
			
		// Specific Heat [J/kgK]
		// ----------------------------------------------------------------------
		mix->SpeciesCp(T);
		Cp = mix->MixCp_FromMassFractions(omega);
		Cv = Cp-Constants::R_J_kmol/MWtot;

		if (heat_transfer_model != ICEM_NONE)
		{
			// Thermal conductivity [W/mK]
			mix->SpeciesConductivityFromFitting(T);
			lambda = mix->MixConductivity_FromMolarFractions(x);

			// Dynamic viscosity [Pa.s]
			mix->SpeciesViscosityFromFitting(T);
			mu = mix->MixViscosity_FromMolarFractions(x);
		}

		// Table Creation (Memorization)
		// ----------------------------------------------------------------------
		mix->ComputeKineticParameters( T, log(T), 1./T, P);
		mix->ComputeFromConcentrations( T, c, cTot, &R);		// [kmol/m3/s]
		ElementByElementProduct(R, mix->M, &R);					// [kg/m3/s]
		QReaction = - mix->ComputeQReaction(T);					// [J/m3/s]

		mix->GetStandardInternalEnergy_Mass(u, T);					// [J/kg]
		UReaction = -Dot(u,R);									// [J/m3/s]	
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

		if (heat_transfer_model != ICEM_NONE)
		{
			// Thermal conductivity [W/mK]
			mix->SpeciesConductivityFromFitting(T);
			lambdaMap = mix->lambda;

			// Dynamic viscosity [Pa.s]
			mix->SpeciesViscosityFromFitting(T);
			muMap = mix->eta;
		}

		// Table Creation (Memorization)
		// ----------------------------------------------------------------------
		mix->ComputeKineticParameters(T, log(T), 1./T, P);
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
        // a. Calcolo della concentrazione e della densita [kmol/m3] [kg/m3]
		moles	= mass / MWtot;						// [kmol]
		rho		= mass / volume;					// [kg/m3]
		cTot	= rho  / MWtot;						// [kmol/m3]
		P		= cTot * (Constants::R_J_kmol*T);	// [Pa]
		c       = cTot * x;							// [kmol/m3]

		// Specific Heat [J/kgK]
		// ----------------------------------------------------------------------
		mix->Cp = CpMap;
		Cp = mix->MixCp_FromMassFractions(omega);
		Cv = Cp-Constants::R_J_kmol/MWtot;

		if (heat_transfer_model != ICEM_NONE)
		{
			// Thermal conductivity [W/mK]
			mix->lambda = lambdaMap;
			lambda = mix->MixConductivity_FromMolarFractions(x);

			// Dynamic viscosity [Pa.s]
			mix->eta = muMap;
			mu = mix->MixViscosity_FromMolarFractions(x);
		}


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

		mix->GetStandardInternalEnergy_Mass(u, T);				// [J/kg]
		UReaction = -Dot(u,R);									// [J/m3/s]	
	}

	// Global Kinetics
	if (iGlobalKinetics == true)
	{
		global->GiveMeFormationRates(T,c,R);
		global->GiveMeReactionHeat(T, R, QReaction);
	}
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
//																								   //
//									UPDATING FUNCTIONS											   //
//																								   //
/////////////////////////////////////////////////////////////////////////////////////////////////////

void OpenSMOKE_ICEM::UpdateTwoEquationModel(BzzVector &y, BzzVector &dy)
{
    soot2EModel->update(T, P/101325., rho, 0., x, y[indexTwoEquations], y[indexTwoEquations+1]);
    soot2EModel->formation_rates();

    for (int i=1;i<=NC;i++)
        domega[i] += soot2EModel->SGas[i] / rho;

    dy[indexTwoEquations]	= soot2EModel->s / rho;
    dy[indexTwoEquations+1] = soot2EModel->S / rho;
}

void OpenSMOKE_ICEM::UpdateHeatFlux(const double tau, const double csi)
{
	if (iUserDefinedHeatFlux == USERDEFINED)
		Qe = ud_Qe_profile.GiveMeValue(tau, MINUSONE);
	
	else if (heat_transfer_model != NONE)
	{
		if (iUserDefinedAmbientTemperature		== USERDEFINED)		Tambient	= ud_Tambient_profile.GiveMeValue(tau, MINUSONE);
		
		if (iUserDefinedHeatExchangeCoefficient == USERDEFINED)		U			= ud_U_profile.GiveMeValue(tau, MINUSONE);
		else UpdateHeatTransfer();

		Qe = U*(T-Tambient);
	}
}

void OpenSMOKE_ICEM::UpdateExchangeArea(const double tau, const double csi)
{
	if (iUserDefinedExchangeArea == NONE)			A = 2.*area_base + area_lateral;
	if (iUserDefinedExchangeArea == USERDEFINED)	A = ud_A_profile.GiveMeValue(tau, MINUSONE);
}

void OpenSMOKE_ICEM::UpdateAngleAndVolume(const double t)
{
	double sqrt_coeff;
	teta	= rotation_rate*t - start_angle; 

	sqrt_coeff = sqrt(arm_ratio*arm_ratio-BzzPow2(sin(teta)));
	volume	= volume_clearance * (1.+0.50*(compression_ratio-1.)*(arm_ratio+1.-cos(teta)-sqrt_coeff));
	dvolume_over_dt = volume_clearance * rotation_rate*0.50*(compression_ratio-1.)*sin(teta)*(1.+cos(teta)/sqrt_coeff);

	height	= volume/area_base; 
	area_lateral = 4.*volume/diameter;
}

void OpenSMOKE_ICEM::UpdateHeatTransfer()
{
	if (heat_transfer_model == ICEM_NONE)
	{
		U = 0.;
	}

	else if (heat_transfer_model == ICEM_BASE)
	{
	//	mu = 5.42e-5;
	//	mu = 1.07e-1;
		Re = rho*velocity_piston*diameter / mu;		// [-]
		Nu = a_Nu*pow(Re, b_Nu);//*pow(Pr, c_Nu);		// [-]
		U  = Nu*lambda/diameter;					// [W/m2/K]
	//	cout << mu << endl;
	//	cout << lambda << endl;
	//	getchar();
	}
}


/////////////////////////////////////////////////////////////////////////////////////////////////////
//																								   //
//									SOLVING FUNCTIONS											   //
//																								   //
/////////////////////////////////////////////////////////////////////////////////////////////////////

void OpenSMOKE_ICEM::Initialize()
{
	ChangeDimensions(NC, &u);
	ChangeDimensions(NC, &lambdaMap);
	ChangeDimensions(NC, &muMap);
	
	UpdateAngleAndVolume(0.);
	arm_La = (compression_ratio-1.)/2.*volume_clearance/area_base;
	arm_Lc = arm_ratio*arm_La;
	velocity_piston = 2.*arm_La*rotation_rate;	// [m/s]

    T            = inletStream->T;
	P            = inletStream->P;
    omega        = inletStream->omega;
	mass         = inletStream->rho*volume;

	UpdateProperties(MINUSONE, NC+1);
}

void OpenSMOKE_ICEM::Solve()
{
    double timeStart, timeEnd;
    BzzVector xMin;
    BzzVector xMax;
    BzzVector xInitial;

	PrepareFiles();


    ODE_ICEM_Object.assignICEM(this);

    // 1.A 
    if (iTwoEquationModel == false)
    {
        ChangeDimensions(NC+1, &xMin);
        xMin=ZERO;
        xMin[NC+1] = 250.;
        
		ChangeDimensions(NC+1, &xMax);
        xMax=ONE;
        xMax[NC+1] = 6000.;
        
        xInitial = omega;
        xInitial.Append(T);
    }

    // 1.B
    if (iTwoEquationModel == true)
    {
        soot2EModel->initial_values(rho);
        indexTwoEquations = NC+2;

        ChangeDimensions(NC+3, &xMin);
        xMin=ZERO;
        xMin[NC+1] = 250.;
        
		ChangeDimensions(NC+3, &xMax);
        xMax=ONE;
        xMax[NC+1] = 6000.;

        xInitial = omega;
        xInitial.Append(T);
        xInitial.Append(soot2EModel->phiNStart);
        xInitial.Append(soot2EModel->phiMStart);
    }



    o(xInitial, 0., &ODE_ICEM_Object);

    o.MyStepPrint();
    o.SetMinimumConstraints(xMin);
    o.SetMaximumConstraints(xMax);

    if (iRelativeTolerance == true)	o.SetTolRel(MachEpsFloat()*relativeTolerance);
    if (iAbsoluteTolerance == true)	o.SetTolAbs(absoluteTolerance);

    {
        countGlobalIterations = -1;  // From -1 to avoid to store results from the first iteration

        timeStart = BzzGetCpuTime();
		
		TauTotal = 2.*Constants::pi/rotation_rate * number_of_cycles;
		o(TauTotal,TauTotal);
		
		timeEnd = BzzGetCpuTime();
    }


    if (iVerbose == true)
    {
        cout << endl;
        cout << "Number of function for Jacobian: " << o.GetNumFunctionForJacobian()		<< endl;
        cout << "Numerical Jacobians: "				<< o.GetNumNumericalJacobian()			<< endl;
        cout << "Time DAE solution: "				<< timeEnd - timeStart	<< " s"			<< endl << endl;
    }

	if (iVerboseSensitivity == true)
	{
	}

	CloseFiles();
}

void OpenSMOKE_ICEM::DefineFromFile(const std::string inputFile)
{
    double			double_value;
    std::string			string_value;
    int				int_value;
	vector<string>  string_vector;

    OpenSMOKE_Dictionary_ICEM dictionary;
    dictionary.ParseFile(inputFile);


    // COMPULSORY: Rotation Rate
	if (dictionary.Return("#RotationRate", double_value, string_value))
		AssignRotationRate(string_value, double_value);

    // COMPULSORY: Compression Ratio
	if (dictionary.Return("#CompressionRatio", double_value))
		AssignCompressionRatio(double_value);

    // COMPULSORY: Reactor Clearance Volume
    if (dictionary.Return("#ClearanceVolume", double_value, string_value))
        AssignClearanceVolume(string_value, double_value);

    // COMPULSORY: Arm Ratio
    if (dictionary.Return("#ArmRatio", double_value))
        AssignArmRatio(double_value);

    // COMPULSORY: Start angle
    if (dictionary.Return("#StartAngle", double_value, string_value))
        AssignStartAngle(string_value, double_value);
	
	// COMPULSORY: Start angle
    if (dictionary.Return("#Diameter", double_value, string_value))
        AssignDiameter(string_value, double_value);

	// OPTIONAL: Number of cycles
    if (dictionary.Return("#NumberOfCycles", double_value))
        SetNumberOfCycles(double_value);



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

	// OPTIONAL: Adiabatic
    if (dictionary.Return("#Adiabatic"))
        SetAdiabatic();


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

	// OPTIONAL: Sensitivity Analysis
    if (dictionary.Return("#Sensitivity", string_vector))
		SetSensitivityOnFile(string_vector);

	Lock();
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
//																								   //
//									OUTPUT FUNCTIONS											   //
//																								   //
/////////////////////////////////////////////////////////////////////////////////////////////////////

void OpenSMOKE_ICEM::LabelODEFile()
{
    int i;

    if (iVerbose == true)
    {
        {
			fOutput << "Tau[s](1)    "		<< "\t";
            fOutput << "Angle[dg](2) "		<< "\t";
			fOutput << "T[K](3)      "		<< "\t";
            fOutput << "V[l](4)      "		<< "\t";
			fOutput << "P[atm](5)    "		<< "\t";
            fOutput << "H[m](6)      "		<< "\t";
            fOutput << "QR[W/m3](7)  "		<< "\t";
            fOutput << "QE[W/m3](8)  "		<< "\t";
            fOutput << "rho[kg/m3](9)"	    << "\t";
            fOutput << "MW(10)       "	    << "\t";
            fOutput << "UR[W/m3](11) "	    << "\t";
            fOutput << "dummy(12)    "    << "\t";
			if(key_species_names.size() == 1)
				fOutput << "Eta_" << key_species_names[0] << "(13)" << "\t";
			else
	            fOutput << "dummy(13)    "	    << "\t";

            int fOutputCount = 14;
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
/*        if (mix->soot_manager.iSoot == 1)
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

		if (iVerboseEnergy == true)
		{
			fEnergy << "Tau[s](1)      "	<< "\t";
            fEnergy << "Angle[dg](2)   "	<< "\t";
            fEnergy << "U[J/s](3)      "	<< "\t"; 
            fEnergy << "H[J/kg](4)     "	<< "\t";
			fEnergy << "Ek[J/kg](5)    "	<< "\t";
            fEnergy << "Ep[J/kg](6)    "	<< "\t";
            fEnergy << "QR[W/m3](7)    "	<< "\t";
            fEnergy << "QE[W/m3](8)    "	<< "\t";
            fEnergy << "Cp[J/kg/K](9)  "	<< "\t";
            fEnergy << "Cv[J/kg/K](10) "	<< "\t";
            fEnergy << "UR[W/m3](11)   "	<< "\t";
            fEnergy << "h[W/m2/K](12)  "	<< "\t";
            fEnergy << "Nu[-](13)      "	<< "\t";
            fEnergy << "Re[-](14)      "	<< "\t";
            fEnergy << "Pr[-](15)      "	<< "\t";
            fEnergy << "A[cm2](16)     "	<< "\t";
            fEnergy << "Ab[cm2](17)    "	<< "\t";
            fEnergy << "Al[cm2](18)    "	<< "\t";
            fEnergy << "H[cm](19)      "	<< "\t";
			fEnergy << endl;
			fEnergy << endl;
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

		if (iAssignedROPA == true)
		{
		}

		if (iVerboseSensitivity == true)
		{
		}
    }
}

void OpenSMOKE_ICEM::VideoGeneralInfo()
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
    cout << "Volume:\t\t\t"				<< volume*1.e6						<< " [cm3]"		<< endl;
    cout << "Height:\t\t\t"				<< height*100.						<< " [cm]"		<< endl;
	cout << "Temperature:\t\t"			<< T								<< " [K]"		<< endl;
	cout << "Pressure:\t\t"				<< P								<< " [Pa]"		<< endl;
    cout << "Density:\t\t"				<< rho								<< " [kg/m3]"	<< endl;
    cout << "Molecular weight:\t"		<< MWtot							<< " [kg/kmol]" << endl;
    cout << "Mass:\t\t\t"				<< mass								<< " [kg]"		<< endl;
    cout << "Moles:\t\t\t"				<< moles							<< " [kmol]"	<< endl;
    cout << "Mean velocity:\t\t"		<< velocity_piston					<< " [m/s]"		<< endl;
    cout << "Arm length (La):\t"		<< arm_La*100.						<< " [cm]"		<< endl;
    cout << "Arm length (Lc):\t"		<< arm_Lc*100.						<< " [cm]"		<< endl;
	cout <<  endl;
}

void OpenSMOKE_ICEM::VideoSummary()
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

void OpenSMOKE_ICEM::VideoFinalResult()
{
    int i;

    // Vectors
    BzzVector conversion(NC);

    // Print General Information about Batch
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
    cout << "#\tName\tomegaInitial\tomegaFinal\tConversion(%)" << endl;
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

void OpenSMOKE_ICEM::SummaryOnFile()
{
	int i;

	ofstream fSummary;
	std::string file_name = outputFolderName + "/Summary.out";
	openOutputFileAndControl(fSummary, file_name);
    fSummary.setf(ios::scientific);

	double initial_volume	= mass / inletStream->rho;	// [m3]
	double initial_moles	= mass / inletStream->MW ;	// [kmol]


    fSummary << "               \tInitial      \tEnd"								<< endl;
	fSummary << "----------------------------------------------------------------------------"	<< endl;
	fSummary << "Angle[deg]     \t"	<< start_angle*rad_to_degrees << "\t" << teta*rad_to_degrees	<< endl;
	fSummary << "Volume[m3]     \t"	<< initial_volume		<< "\t" << volume			<< endl;	
	fSummary << "Temperature[K] \t"	<< inletStream->T		<< "\t" << T				<< endl;
	fSummary << "Pressure[Pa]   \t"	<< inletStream->P		<< "\t" << P				<< endl;
	fSummary << "Mass[kg]       \t"	<< mass					<< "\t" << mass				<< endl;
	fSummary << "Moles[kmol]    \t"	<< initial_moles		<< "\t" << moles			<< endl;	
	fSummary << "Density[kg/m3] \t"	<< inletStream->rho		<< "\t" << rho				<< endl;
	fSummary << "MW[kg/kmol]    \t"	<< inletStream->MW		<< "\t" << MWtot			<< endl;
	fSummary << endl;
	
	
	fSummary << "Mole fractions[-]" << endl;
	for(i=1;i<=NC;i++)
	{
		fSummary.width(12);
		fSummary << mix->names[i] << "    \t"	<<	inletStream->x[i]		<< "\t" << x[i]			<< endl;
	}
	fSummary << endl;

	fSummary << endl;
	fSummary << "Mass fractions[-]" << endl;
	for(i=1;i<=NC;i++)
	{
		fSummary.width(12);
		fSummary << mix->names[i] << "    \t"	<<	inletStream->omega[i]	<< "\t" << omega[i]		<< endl;
	}
	fSummary << endl;
	
	fSummary << endl;
	fSummary << "Concentrations[kmol/m3]" << endl;
	for(i=1;i<=NC;i++)
	{
		fSummary.width(12);
		fSummary << mix->names[i] << "    \t"	<<	inletStream->c[i]		<< "\t" << c[i]			<< endl;
	}
	fSummary << endl;

	fSummary << endl;
	fSummary << "Conversions[-]" << endl;
	for(i=1;i<=NC;i++)
	{
		fSummary.width(12);
		fSummary << mix->names[i] << "    \t"	<<	1.-omega[i]/inletStream->omega[i]	<< endl;
	}
	fSummary << endl;

    fSummary.close();
}

void OpenSMOKE_ICEM::ODEPrint(BzzVector &y, double t)
{
	double dummy = ZERO;

    if (iVerbose == true)
    {
		double teta_deg = teta*rad_to_degrees;

        if ( countIterations%(300*nVideoSteps) == 0)
        {
            cout    << endl
                    << setw(8)  << left << "#"
                    << setw(16) << left << "Tau[s]"
                    << setw(16) << left << "Angle[deg]"
                    << setw(16) << left << "T[K]"
					<< setw(16) << left << "P[atm]";
                    
			cout	<< endl;
        }

        if (countVideoSteps == nVideoSteps)
        {
            cout	<< setw(8)  << left << countIterations
                    << setw(16) << left << t
                    << setw(16) << left << teta_deg
                    << setw(16) << left << T
                    << setw(16) << left << P/101325.;

            // Soot Two Equation Model
            if (iTwoEquationModel == true)
            {
                cout 	<< soot2EModel->m0   << "\t"
                        << soot2EModel->fv	<< "\t";
            }

            cout	<< endl;

            countVideoSteps = 0;
        }

        if (countFileSteps == nFileSteps)
        {
            int i;

            fOutput << t					<< "\t"
                    << teta_deg				<< "\t"
                    << T					<< "\t"
                    << volume*1000.			<< "\t"
                    << P/101325.			<< "\t"
                    << height				<< "\t"
                    << QReaction			<< "\t"
                    << Qe*A/volume     		<< "\t"
                    << rho					<< "\t"
                    << MWtot				<< "\t"
                    << UReaction			<< "\t"
                    << dummy				<< "\t";

			if(key_species_names.size() == 1)
				fOutput << 1.-omega[key_species_index[0]]/inletStream->omega[key_species_index[0]]<< "\t";
			else
	            fOutput << dummy			<< "\t";

            for (i=1;i<=NC;i++)
                fOutput << x[i]				<< "\t";
            for (i=1;i<=NC;i++)
                fOutput << omega[i]			<< "\t";
            fOutput << endl;

			// Soot Two Equation Model
            if (iTwoEquationModel == true)
			{
				fSoot2E	<< setw(20) << left << t;
				soot2EModel->write_on_file(fSoot2E, y[indexTwoEquations], y[indexTwoEquations+1]);
			}

   /*         if (mix->soot_manager.iSoot == 1)
            {
                // Soot Distribution data
                mix->soot_manager.large_calculateAll(omega, x, c, rho);
                mix->soot_manager.small_calculateAll(omega, x, c, rho);

                fSoot << t			<< "\t";
                fSoot << teta_deg	<< "\t";
                fSoot << T			<< "\t";
                mix->soot_manager.print_main_data_on_file(fSoot);
                fSoot << endl;
            }*/
/*
            if (mix->pah_manager.iPAH == 1)
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

			if (iVerboseEnergy == true)
			{
				BzzVector h_mass(NC);
				mix->GetMixAveragedEnthalpy_Mass(h_mass, T);
	
				double H_mass  = Dot(omega, h_mass);			// [J/kg]
				double Ep_mass = P/rho;							// [J/kg]
				double Ek_mass = 0.;							// [J/kg]

				fEnergy << t									<< "\t"
						<< teta_deg								<< "\t"
						<< H_mass+Ep_mass						<< "\t"
						<< H_mass								<< "\t"
						<< Ek_mass								<< "\t"
						<< Ep_mass								<< "\t"
						<< QReaction							<< "\t"
						<< Qe*A		       						<< "\t"
						<< Cp		       						<< "\t"
						<< Cv		       						<< "\t"
						<< UReaction      						<< "\t"
						<< U									<< "\t"
						<< Nu									<< "\t"
						<< Re									<< "\t"
						<< Pr									<< "\t"
						<< A*1.e4								<< "\t"
						<< area_base*1.e4						<< "\t"
						<< area_lateral*1.e4					<< "\t"
						<< height*1.e2							<< "\t"
						<< endl;
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

		/*	if (iVerboseROPA == true)
			{
				ROPA.Run(mix->r);

				fROPA	<< t		<< "\t"
						<< teta_deg	<< "\t"
						<< T		<< "\t";

				ROPA.PrintRateOfProductionAnalyses(fROPA);
			}
*/
			if (iVerboseSensitivity == true)
			{
				if ( (t-timeOld)>-1 && t>0. )
				{
					fSensitivity	<< t		<< "\t"
									<< teta_deg	<< "\t"
									<< T		<< "\t";
					Sensitivity.PrintOnFile_SensitivityCoefficients(fSensitivity);
				}
			}

            countFileSteps = 0;
        }

    }

    countVideoSteps++;
    countFileSteps++;
    countIterations++;
	timeOld = t;

    if (iHistory == true)
    {
        countGlobalIterations++;

        if(countGlobalIterations>0)
        {
            Tau_History[countGlobalIterations] = t;
            T_History[countGlobalIterations]   = T;
            mass_History.SetRow(countGlobalIterations, omega);
            mole_History.SetRow(countGlobalIterations, x);
        }
    }
}

void OpenSMOKE_ICEM::EnergyAnalysis(OpenSMOKE_GasStream &outletStream)
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

	cout << "Initial energy:     "	<< U_mass_initial*mass								<< " [J]" << endl;
	cout << "Final energy:       "	<< U_mass_final*mass								<< " [J]" << endl;
	cout << "Balance (IN-OUT):   "  << U_mass_initial*mass - U_mass_final*mass			<< " [J]" << endl;
}


/////////////////////////////////////////////////////////////////////////////////////////////////////
//																								   //
//						ODE SYSTEM EQUATIONS - NON ISOTHERMAL REACTOR - CONST V					   //
//																								   //
/////////////////////////////////////////////////////////////////////////////////////////////////////

void OpenSMOKE_ICEM::ODESystem_ICEM(BzzVector &y, double t, BzzVector &dy)
{
    int i;

	UpdateAngleAndVolume(t);

    // Recovering variables;
    for (i=1;i<=NC;i++)	omega[i]    = y[i];
    T				                = y[NC+1];

	// Updating Exchange area
	UpdateExchangeArea(t, MINUSONE);
    UpdateHeatFlux(t, MINUSONE);


    // Evaluation of properties
    UpdateProperties(ODE_ICEM_Object.jacobianIndex, NC+1);

    // Conservation equations for all the species (Batch)
    domega	= R / rho ;

    // Conservation equations for energy (Batch)
    dT   = volume*UReaction - P*dvolume_over_dt - A*Qe;
    dT  /= (mass*Cv);
 
    // Soot Two Equation Model
    if (iTwoEquationModel == true)	UpdateTwoEquationModel(y, dy);

    // Recovering residuals (Soot 2E are already updated)
    for (i=1;i<=NC;i++)	dy[i]	= domega[i];
    dy[NC+1]					= dT;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
//																								   //
//											DICTIONARY											   //
//																								   //
/////////////////////////////////////////////////////////////////////////////////////////////////////

OpenSMOKE_Dictionary_ICEM::OpenSMOKE_Dictionary_ICEM()
{
    SetupBase();
    
	// Compulsory
    Add("#RotationRate",		'C', 'M', "Rotation rate");
    Add("#CompressionRatio",	'C', 'D', "Compression ratio");
    Add("#ClearanceVolume",		'C', 'M', "Clearance volume");
    Add("#ArmRatio",			'C', 'D', "Ratio between arms");
    Add("#StartAngle",			'C', 'M', "Starting angle");
    Add("#Diameter",			'C', 'M', "Cylinder diameter");

	// Optional
    Add("#Key",					'O', 'V', "Key species");
	Add("#NumberOfCycles",		'O', 'D', "Number of cycles");

	// Optional
	Add("#ReactionRates",		'O', 'V', "Reaction rates");
	Add("#FormationRates",		'O', 'V', "Formation rates");
	Add("#ROPA",				'O', 'V', "Rate of Production Analysis");
	Add("#Sensitivity",			'O', 'V', "Sensitivity Analysis");

     
	// Optional
	Add("#Adiabatic",				'O', 'N', "No thermal exchange");
	Add("#VerboseEnergy",			'O', 'N', "Report on energy");
	Add("#Qe",						'O', 'M', "Heat flux");
	Add("#A",						'O', 'M', "Exchange area");
	Add("#U",						'O', 'M', "Heat thermal exchange coefficient");
	Add("#Twall",					'O', 'M', "Wall temperature");

	Add("#UserDefined_T",			'O', 'S', "The temperature is assigned from external file");
	Add("#UserDefined_Qe",			'O', 'S', "The heat flux is assigned from external file");
    Add("#UserDefined_A",			'O', 'S', "The exchange area is assigned from external file");
	Add("#UserDefined_U",			'O', 'S', "The heat exchange coefficient is assigned from external file");
    Add("#UserDefined_Twall",		'O', 'S', "The wall temperature is assigned from external file");

	Add("#2E_Model",				'O', 'N', "The Soot 2E Model is solved");

	Conflict("#UserDefined_Qe",			"#Qe");
	Conflict("#UserDefined_A",			"#A");
	Conflict("#UserDefined_U",			"#U");
	Conflict("#UserDefined_Twall",		"#Twall");

	Conflict("#Qe",						"#U");
	Conflict("#Qe",						"#Twall");
	Conflict("#Qe",						"#UserDefined_U");
	Conflict("#Qe",						"#UserDefined_Twall");
	Conflict("#UserDefined_Qe",			"#U");
	Conflict("#UserDefined_Qe",			"#Twall");
	Conflict("#UserDefined_Qe",			"#UserDefined_U");
	Conflict("#UserDefined_Qe",			"#UserDefined_Twall");

	Compulsory("#Adiabatic", "#Twall", "#UserDefined_Twall", "#UserDefined_Qe");


    Lock();
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
//																								   //
//									ODE SYSTEM - CLASS DEFINITION								   //
//																								   //
/////////////////////////////////////////////////////////////////////////////////////////////////////

void MyOdeSystem_ICEM::ObjectBzzPrint(void)
{
}

void MyOdeSystem_ICEM::GetSystemFunctions(BzzVector &x, double eta, BzzVector &f)
{
	ptICEM->ODESystem_ICEM(x, eta, f);
}

void MyOdeSystem_ICEM::assignICEM(OpenSMOKE_ICEM *icem)
{
    ptICEM				= icem;
}

void MyOdeSystem_ICEM::MyODEPrint(BzzVector &y, double eta)
{
	ptICEM->ODEPrint(y, eta);
}
