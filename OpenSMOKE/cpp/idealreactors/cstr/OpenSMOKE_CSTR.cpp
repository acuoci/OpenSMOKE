/***************************************************************************
 *   Copyright (C) 2006-2009 by Alberto Cuoci							   *
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
#include "idealreactors/cstr/OpenSMOKE_CSTR.h"
#include "basic/OpenSMOKE_Conversions.h"
#include "surfaceChemistry/OpenSMOKE_ReactingSurface.h"
#include "surfaceChemistry/OpenSMOKE_SurfaceMaterial.h"
#include "surfaceChemistry/OpenSMOKE_SurfaceKinetics.h"
#include "surfaceChemistry/OpenSMOKE_SurfaceReaction.h"
//#include "addons/OpenSMOKE_PostProcessor.h"
#include "addons/OpenSMOKE_SensitivityAnalysis_Fast_Unsteady0D.h"

OpenSMOKE_CSTR	*ptCSTR;

void ODEPrintExternal_CSTR(BzzVector &y, double t)
{
    ptCSTR->ODEPrint(y, t);
}

OpenSMOKE_CSTR::OpenSMOKE_CSTR()
{
    ptCSTR		=  this;
	class_name	= "OpenSMOKE_CSTR";
	out_name	= "CSTR.out";

	tEnd					= 1.e9;		// integration time [s]
	iSetInitialTemperature	= false;
	iSetInitialComposition	= false;
	iSetFluctuations		= false;
	iSetGasDischargeLaw		= false;
	iHistoryPartial			= false;
	
	Setup();
}

void OpenSMOKE_CSTR::AssignEnd(const std::string units, const double value)
{
	tEnd = OpenSMOKE_Conversions::conversion_time(value, units);
}

void OpenSMOKE_CSTR::SetValveCoefficient(const double value, const std::string units)
{
	KValve = OpenSMOKE_Conversions::conversion_valve_flow_coefficient(value, units);
	PressureOutlet = inletStream->P;
	iSetGasDischargeLaw = true;
}

void OpenSMOKE_CSTR::SetPressureOutlet(const double value, const std::string units)
{
	PressureOutlet = OpenSMOKE_Conversions::conversion_pressure(value, units);
}

void OpenSMOKE_CSTR::AssignResidenceTime(const std::string units, const double value)
{
	iResidenceTime = true;

	Tau = OpenSMOKE_Conversions::conversion_time(value, units);
	assignedGeometry	= false;
	assignedEnd			= true;
}

void OpenSMOKE_CSTR::AssignVolume(const std::string units, const double value)
{
	iResidenceTime = false;

	Volume  = OpenSMOKE_Conversions::conversion_volume(value, units);
	D		= pow(Volume/(Constants::pi/6.), 1./3.);
	Area	= Constants::pi * (D*D);

	if (iUserDefinedExchangeArea == NONE)		Ae = Area;
	assignedGeometry	= false;
	assignedEnd			= true;
}

void OpenSMOKE_CSTR::AssignDiameter(const std::string units, const double value)
{
	iResidenceTime = false;

    D       = OpenSMOKE_Conversions::conversion_length(value, units);
    Area    = Constants::pi * (D*D);
    Volume  = Constants::pi/6. * (D*D*D);
	
	if (iUserDefinedExchangeArea == NONE)		Ae = Area;
	assignedGeometry	= false;
	assignedEnd			= true;
}

void OpenSMOKE_CSTR::SetUserDefinedDiameter(const std::string fileName)
{
	iResidenceTime = false;

	geometry.Setup(fileName);
	geometry.Update(0., D, Volume, Area);

	if (iUserDefinedExchangeArea == NONE)		Ae = Area;

	assignedGeometry	= true;
	assignedEnd			= true;
}

void OpenSMOKE_CSTR::SetUserDefinedResidenceTime(const std::string fileName)
{
	iResidenceTime = true;

	geometry.Setup(fileName);
	geometry.Update(0., Tau);

	assignedGeometry	= true;
	assignedEnd			= true;
}

void OpenSMOKE_CSTR::SetUserDefinedVolume(const std::string fileName)
{
	iResidenceTime = false;

	geometry.Setup(fileName);
	geometry.Update(0., D, Volume, Area);

	if (iUserDefinedExchangeArea == NONE)		Ae = Area;

	assignedGeometry	= true;
	assignedEnd			= true;
}

void OpenSMOKE_CSTR::SetInitialTemperature(const double value, const std::string units)
{
	Tinitial = OpenSMOKE_Conversions::conversion_temperature(value, units);
	iSetInitialTemperature = true;
}

void OpenSMOKE_CSTR::SetInitialMassFractions(const vector<string> _names, const vector<double> _values)
{
	ChangeDimensions(NC, &omega_initial);
	
	for(int i=1;i<=int(_names.size());i++)
		omega_initial[mix->recognize_species(_names[i-1])] = _values[i-1];
	iSetInitialComposition = true;
}

void OpenSMOKE_CSTR::SetInitialMoleFractions(const vector<string> _names, const vector<double> _values)
{
	double MWmix;

	BzzVector x_initial(NC);
	ChangeDimensions(NC, &omega_initial);
	
	for(int i=1;i<=int(_names.size());i++)
		x_initial[mix->recognize_species(_names[i-1])] = _values[i-1];

	mix->GetMWAndMassFractionsFromMoleFractions(MWmix, omega_initial, x_initial);
	iSetInitialComposition = true;
}

void OpenSMOKE_CSTR::Update(const double t)
{
	if (iSetGasDischargeLaw == false)
	{
		mass_flow_out = inletStream->massFlowRate;	// [kg/s]

		if (assignedGeometry == false)
		{
			// Case A: volume and mass flow rate are assigned; the residence time is calculated
			if (iResidenceTime == false)					
				Tau = rho*Volume / inletStream->massFlowRate;
			// Case B: residence time and mass flow rate are assigned; the volume is calculated
			else											 
			{
				Volume = inletStream->massFlowRate * Tau / rho;
				D		= pow(Volume/(Constants::pi/6.), 1./3.);
				Area	= Constants::pi * (D*D);
				if (iUserDefinedExchangeArea == NONE)		Ae = Area;
			}
		}
		else
		{
			// Case A: volume and mass flow rate are assigned; the residence time is calculated
			if (iResidenceTime == false)					
				geometry.Update(t, Tau);
			// Case B: residence time and mass flow rate are assigned; the volume is calculated
			else											 
			{
				geometry.Update(t, D, Volume, Area);
				if (iUserDefinedExchangeArea == NONE)		Ae = Area;
			}
		}
	}
	else
	{
		if  (P<PressureOutlet)	mass_flow_out = 0.;
		else					mass_flow_out = KValve*sqrt(P-PressureOutlet);

		if (assignedGeometry == false)
		{
			// Case A: volume and inlet mass flow rate are assigned
			//         the residence time is calculated accordingly
			if (iResidenceTime == false)					
				Tau = mass_total/inletStream->massFlowRate;
			// Case B: residence time and mass flow rate are assigned; the volume is calculated
			else											 
				ErrorMessage("If gas discharge law is enabled, the residence time cannot be specified...");
		}
		else
			ErrorMessage("If gas discharge law is enabled, the volume must be constant...");
	}
}

void OpenSMOKE_CSTR::SetUserDefinedExchangeArea(const std::string fileName)
{
    iUserDefinedExchangeArea = USERDEFINED;
    ud_Ae_profile.AssignFromFile(fileName, "AREA");
	ud_Ae_profile.SetName(name_object + " - Exchange Area Profile");
}

void OpenSMOKE_CSTR::SetConstantExchangeArea(const double value, const std::string units)
{
	iUserDefinedExchangeArea = CONSTANT;
	Ae = OpenSMOKE_Conversions::conversion_area(value, units);
}

void OpenSMOKE_CSTR::UnsetUserDefinedExchangeArea()
{
    iUserDefinedExchangeArea = NONE;
}

void OpenSMOKE_CSTR::SetFluctuations(const vector<string> string_vector)
{
	int total_size = string_vector.size();
	int n = total_size / 3;
	if (total_size != 12)
		ErrorMessage("Wrong list of arguments in #Fluctuations option...");

	for(int j=1;j<=n;j++)
	{
		double mean;
		double frequency;
		double semiAmplitude;
		double tWall;

		std::string s1 = string_vector[3*(j-1)+0];
		std::string s2 = string_vector[3*(j-1)+1];
		std::string s3 = string_vector[3*(j-1)+2];
		
		if (s1 == "T")
			mean = OpenSMOKE_Conversions::conversion_temperature(atof(s2.c_str()), s3);
		else if (s1 == "frequency")
			frequency = OpenSMOKE_Conversions::conversion_frequency(atof(s2.c_str()), s3);
		else if (s1 == "semiAmplitude")
			semiAmplitude = OpenSMOKE_Conversions::conversion_temperature(atof(s2.c_str()), s3);
		else if (s1 == "tWall")
			tWall = OpenSMOKE_Conversions::conversion_time(atof(s2.c_str()), s3);
		else
			ErrorMessage("Wrong argument in #Fluctuations option: " + s1);
		
		iSetFluctuations = true;
		fluctuations.Setup(mean, frequency, semiAmplitude, tWall);
		fluctuations.SetPointer(T);
	}
}


void OpenSMOKE_CSTR::Lock()
{
    if (assignedKineticScheme == false)
        ErrorMessage("The kinetic scheme was not defined!");
    if (assignedGlobalKineticScheme == false)
        ErrorMessage("The global kinetic scheme was not defined!");
	if (assignedEnd == false)
		ErrorMessage("The reactor residence time or volume were not defined!");
    if (assignedInletFlows == false)
        ErrorMessage("The inlet flow was not defined!");
    if (assignedEnergy == false)
        ErrorMessage("The kind of reactor was not defined!");
	if (assignedSoot2EModel == false)
        ErrorMessage("The Soot 2E Model was not defined!");

	if (iTwoEquationModel == true && mix->iSootMode == false)
		ErrorMessage("This kinetic scheme was interpreted without the SootMode option and cannot be used together the 2E model.\nMore details on the OpenSMOKE's User Guide.");

	Initialize();
}


/////////////////////////////////////////////////////////////////////////////////////////////////////
//																								   //
//									SOLVING FUNCTIONS											   //
//																								   //
/////////////////////////////////////////////////////////////////////////////////////////////////////

void OpenSMOKE_CSTR::Initialize()
{
	ChangeDimensions(NC, &h);
	ChangeDimensions(NC, &hMap);

    T            = inletStream->T;
    omega        = inletStream->omega;
    P            = inletStream->P;
    massFlowRate = inletStream->massFlowRate;

	if (iSetInitialTemperature == true)	T = Tinitial;
	if (iSetInitialComposition == true)	omega = omega_initial;

	if (iResidenceTime == false)
	{
		double MWinlet = mix->GetMWFromMassFractions(omega);
		mass_total = P*MWinlet / (Constants::R_J_kmol*T) * Volume;
	}
	else
		mass_total = inletStream->massFlowRate*Tau;
	mass_total_initial = mass_total;

    if (iUserDefinedTemperature == true)
        ud_temperature_profile.Check(0., T);

    if (iEnergy == false)
    {
        UpdateProperties_isothermal(0);
        indexProperties = 1;
    }
    else if (iEnergy == true)
        UpdateProperties(MINUSONE, NC+1);

	// History
	if (iHistory == true)
	{
		if (iHistoryPartial == true)
		{
			ChangeDimensions(1, &Tau_History);
			ChangeDimensions(1, &T_History);
			ChangeDimensions(1, mix->NumberOfSpecies(), &mass_History);
			ChangeDimensions(1, mix->NumberOfSpecies(), &mole_History);
		}
		else
		{
			BzzVector		ZeroVector;
			BzzMatrix		ZeroMatrix;
			Tau_History		= ZeroVector;
			T_History		= ZeroVector;
			mass_History	= ZeroMatrix;
			mole_History	= ZeroMatrix;
			ChangeDimensions(MAX_TIME_STEPS, &Tau_History);
			ChangeDimensions(MAX_TIME_STEPS, &T_History);
			ChangeDimensions(MAX_TIME_STEPS, mix->NumberOfSpecies(), &mass_History);
			ChangeDimensions(MAX_TIME_STEPS, mix->NumberOfSpecies(), &mole_History);
		}
	}

	countGlobalIterations = -1;
}

void OpenSMOKE_CSTR::ReSolve()
{
	ODE_CSTR_Object.assignCSTR(this, iEnergy, iSetGasDischargeLaw);
	Initialize();
	Solve();
}

void OpenSMOKE_CSTR::Solve()
{
    ptCSTR = this;

    double timeStart, timeEnd;
    BzzVector xMin;
    BzzVector xMax;
    BzzVector xInitial;

	PrepareFiles();

    ODE_CSTR_Object.assignCSTR(this, iEnergy, iSetGasDischargeLaw);

    // 1.A Isothermal CSTR
    if (iEnergy == false && iTwoEquationModel == false && iSetGasDischargeLaw == false)
    {
        ChangeDimensions(NC, &xMin);  xMin=ZERO;					
        ChangeDimensions(NC, &xMax);  xMax=ONE;					// Mole fractions
                                     
        xInitial = omega;
    }

    // 1.B Isothermal CSTR	+ Two Equation Model
    if (iEnergy == false && iTwoEquationModel == true && iSetGasDischargeLaw == false)
    {
        soot2EModel->initial_values(rho);
        indexTwoEquations = NC+1;

        ChangeDimensions(NC+2, &xMin);
        xMin=ZERO;
        
		ChangeDimensions(NC+2, &xMax);
        xMax=ONE;
        xMax[NC+1] = ONE;
        xMax[NC+2] = ONE;

        xInitial = omega;
        xInitial.Append(soot2EModel->phiNStart);
        xInitial.Append(soot2EModel->phiMStart);
    }

 
    // 2.A NonIsothermal CSTR
    if (iEnergy == true && iTwoEquationModel == false && iSetGasDischargeLaw == false)
    {
        ChangeDimensions(NC+1, &xMin);
        xMin=ZERO;
        xMin[NC+1] = 200.;
        ChangeDimensions(NC+1, &xMax);
        xMax=ONE;
        xMax[NC+1] = 6000.;

        xInitial = omega;
        xInitial.Append(T);
    }

  
    // 2.C NonIsothermal CSTR + TwoEquations
    if (iEnergy == true && iTwoEquationModel == true && iSetGasDischargeLaw == false)
    {
        soot2EModel->initial_values(rho);
        indexTwoEquations = NC+2;

        ChangeDimensions(NC+3, &xMin);
        xMin=ZERO;
        xMin[NC+1] = 200.;
        ChangeDimensions(NC+3, &xMax);
        xMax=ONE;
        xMax[NC+1] = 6000.;
        xMax[NC+2] = ONE;
        xMax[NC+3] = ONE;

        xInitial = omega;
        xInitial.Append(T);
        xInitial.Append(soot2EModel->phiNStart);
        xInitial.Append(soot2EModel->phiMStart);
    }

    // 1.D Isothermal CSTR + Gas Discharge Law
    if (iEnergy == false && iTwoEquationModel == false && iSetGasDischargeLaw == true)
    {
        ChangeDimensions(NC, &xMin);  xMin=ZERO;	// Mass fractions	
		xMin.Append(ZERO);							// Total mass
        ChangeDimensions(NC, &xMax);  xMax=ONE;		// Mass fractions
		xMax.Append(1.e16);							// Total mass

        xInitial = omega;							// Mass fractions
		xInitial.Append(mass_total);				// Total mass
    }

    // 1.D Isothermal CSTR + Gas Discharge Law
    if (iEnergy == true && iTwoEquationModel == false && iSetGasDischargeLaw == true)
    {
        ChangeDimensions(NC, &xMin);  xMin=ZERO;	// Mass fractions	
		xMin.Append(200);							// Temperature
		xMin.Append(ZERO);							// Total mass
        ChangeDimensions(NC, &xMax);  xMax=ONE;		// Mass fractions
		xMax.Append(6000.);							// Temperature
		xMax.Append(1.e16);							// Total mass

        xInitial = omega;							// Mass fractions
		xInitial.Append(T);							// Temperature
		xInitial.Append(mass_total);				// Total mass
    }

    BzzOdeStiffObject o(xInitial, 0., &ODE_CSTR_Object);

    o.StepPrint(ODEPrintExternal_CSTR);
    o.SetMinimumConstraints(xMin);
    o.SetMaximumConstraints(xMax);
	//o.SetMaxStep(MAX_TIME_STEPS);

    if (iRelativeTolerance == true)	o.SetTolRel(relativeTolerance);
    if (iAbsoluteTolerance == true)	o.SetTolAbs(absoluteTolerance);

    {
        countGlobalIterations = -1;  // From -1 to avoid to store results from the first iteration

        timeStart = BzzGetCpuTime();
		
		o(tEnd);

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

	CloseFiles();

	if (iHistory == true && iVerbose == true)
	{
		SaveOnBinaryFile(outputOSMName);
	}
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
//																								   //
//									UPDATING FUNCTIONS											   //
//																								   //
/////////////////////////////////////////////////////////////////////////////////////////////////////

void OpenSMOKE_CSTR::UpdateHeatFlux(const double tau, const double csi)
{
	if (iUserDefinedHeatExchangeCoefficient == USERDEFINED)		U			= ud_U_profile.GiveMeValue(tau, csi);
	if (iUserDefinedAmbientTemperature		== USERDEFINED)		Tambient	= ud_Tambient_profile.GiveMeValue(tau, csi);
	if (iUserDefinedHeatFlux				== USERDEFINED)		Qe			= ud_Qe_profile.GiveMeValue(tau, csi);
	
	if (iUserDefinedHeatExchangeCoefficient != NONE)			Qe = U*(T-Tambient);
}

void OpenSMOKE_CSTR::UpdateExchangeArea(const double tau, const double csi)
{
	if (iUserDefinedExchangeArea == NONE)			Ae = Area;
	if (iUserDefinedExchangeArea == USERDEFINED)	Ae = ud_Ae_profile.GiveMeValue(tau, csi);
}

void OpenSMOKE_CSTR::UpdateProperties_isothermal(int memoIndex)
{
    double cTot;

    mix->GetMWAndMoleFractionsFromMassFractions(MWtot, x, omega);

    // --------------------------------------------------------------------------
    // PROPERTIES FOR CONSTANT T + MEMORIZATION
    // --------------------------------------------------------------------------
    if (memoIndex==0)
    {
        // Update pressure
        // ----------------------------------------------------------------------
		if (iSetGasDischargeLaw == true)
			P = mass_total*Constants::R_J_kmol*T/(MWtot*Volume);

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
        // Update pressure
        // ----------------------------------------------------------------------
		if (iSetGasDischargeLaw == true)
			P = mass_total*Constants::R_J_kmol*T/(MWtot*Volume);

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

	// TODO
	/*
	{
		OpenSMOKE_ReactingSurface surface;
		surface.ReadFromBinaryFile("C:\\OpenSMOKE_Suite\\tutorials\\OpenSMOKE_CSTR\\OpenSMOKE\\Test12-Surface\\KineticScheme_Prova01\\Surface01\\surface.bin");

		BzzVector Z(surface.NumberSiteSpecies());
		BzzVector aBulk(surface.NumberBulkSpecies());

		BzzVector RGas(mix->NumberOfSpecies());
		BzzVector RSurface(surface.NumberSiteSpecies());
		BzzVector RBulk(surface.NumberBulkSpecies());

		Z[1] = 0.50;
		Z[2] = 0.50;
	//	Z[3] = 0.20;
	//	Z[4] = 0.20;
	//	Z[5] = 0.20;
		double QReactionSurface;

		surface.material()[1].ReactionRates(T, c, Z, aBulk, RGas, RSurface, RBulk, QReactionSurface);

		getchar();
	}
	*/
}

void OpenSMOKE_CSTR::UpdateProperties(int jacobianIndex, int indexT)
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
        // Update pressure
        // ----------------------------------------------------------------------
		if (iSetGasDischargeLaw == true)
			P = mass_total*Constants::R_J_kmol*T/(MWtot*Volume);

		// Calcolo della concentrazione e della densita [kmol/m3] [kg/m3]
		// ----------------------------------------------------------------------
		cTot = P  / (Constants::R_J_kmol*T);
		rho = cTot * MWtot;
		c = cTot*x;

		// Specific Heat [J/kgK]
		// ----------------------------------------------------------------------
		mix->SpeciesCp(T);
		Cp = mix->MixCp_FromMassFractions(omega);
		Cv = Cp - Constants::R_J_kmol/MWtot;

		// Table Creation (Memorization)
		// ----------------------------------------------------------------------
		mix->ComputeKineticParameters( T, log(T), 1./T, P);
		mix->ComputeFromConcentrations( T, c, cTot, &Rtilde);	// [kmol/m3/s]
		ElementByElementProduct(Rtilde, mix->M, &R);			// [kg/m3/s]
		QReaction = - mix->ComputeQReaction(T);					// [J/m3/s]

		Product((Constants::R_J_kmol*T), mix->h, &h);			// [J/kmol]
		ElementByElementProduct(h, mix->uM, &h);				// [J/kg]
		specificEnthalpy = Dot(omega,h);						// [J/kg]
		specificInternalEnergy = specificEnthalpy - Constants::R_J_kmol*T/MWtot;
	}

	// --------------------------------------------------------------------------
	// Update only the variables depending on the temperature
	// --------------------------------------------------------------------------
	else if(memoIndex==-1)
	{
        // Update pressure
        // ----------------------------------------------------------------------
		if (iSetGasDischargeLaw == true)
			P = mass_total*Constants::R_J_kmol*T/(MWtot*Volume);

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

		// Specific Enthalpies [J/kg/K]
		// ----------------------------------------------------------------------
		Product((Constants::R_J_kmol*T), mix->h, &h);			// [J/kmol]
		ElementByElementProduct(h, mix->uM, &h);				// [J/kg]
		hMap = h;

		memoIndex = 0;
	}

	// --------------------------------------------------------------------------
	// Update every variable which does not depend on T
	// --------------------------------------------------------------------------
	if(memoIndex == 0)
	{
        // Update pressure
        // ----------------------------------------------------------------------
		if (iSetGasDischargeLaw == true)
			P = mass_total*Constants::R_J_kmol*T/(MWtot*Volume);

		// Calcolo della concentrazione e della densita [kmol/m3] [kg/m3]
		// ----------------------------------------------------------------------
		cTot = P  / (Constants::R_J_kmol*T);
		rho = cTot * MWtot;
		c = cTot*x;

		// Specific Heat [J/kgK]
		// ----------------------------------------------------------------------
		mix->Cp = CpMap;
		Cp = mix->MixCp_FromMassFractions(omega);
		Cv = Cp - Constants::R_J_kmol/MWtot;


		// Table Creation (Memorization)
		// ----------------------------------------------------------------------
		mix->k1			= k1Map;
		mix->k2			= k2Map;
		mix->uKeq		= uKeqMap;
		mix->logFcent	= logFcentMap;
		mix->reactionDH	= reactionDHMap;
		mix->reactionDS	= reactionDSMap;
		mix->ComputeFromConcentrations( T, c, cTot, &Rtilde);		// [kmol/m3/s]
		ElementByElementProduct(Rtilde, mix->M, &R);				// [kg/m3/s]
		QReaction = - mix->ComputeQReaction(T);						// [J/m3/s]

		// Specific Enthalpy [J/kg]
		h = hMap;
		specificEnthalpy = Dot(omega,h);	// [J/kg]
		specificInternalEnergy = specificEnthalpy - Constants::R_J_kmol*T/MWtot;
	}

	// Global Kinetics
	if (iGlobalKinetics == true)
	{
		global->GiveMeFormationRates(T,c,R);
		global->GiveMeReactionHeat(T, R, QReaction);
		global->GiveMeSpecificEnthalpy(T, omega, specificEnthalpy);
		specificInternalEnergy = specificEnthalpy - Constants::R_J_kmol*T/MWtot;
	}
}

void OpenSMOKE_CSTR::UpdateTwoEquationModel(BzzVector &y, BzzVector &dy)
{
	soot2EModel->update(T, P/101325., rho, 0., x, y[indexTwoEquations], y[indexTwoEquations+1]);
    soot2EModel->formation_rates();

    for (int i=1;i<=NC-1;i++)
		domega[i] += soot2EModel->SGas[i] / rho;	// Gas species
	domega[NC] += soot2EModel->S / rho;				// Soot

    dy[indexTwoEquations]	= (soot2EModel->phiNStart-y[indexTwoEquations])   / Tau + soot2EModel->s/rho;
    dy[indexTwoEquations+1] = (soot2EModel->phiMStart-y[indexTwoEquations+1]) / Tau + soot2EModel->S/rho;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
//																								   //
//						ODE SYSTEM EQUATIONS - ISOTHERMAL REACTOR								   //
//																								   //
/////////////////////////////////////////////////////////////////////////////////////////////////////

void OpenSMOKE_CSTR::ODESystem_Isothermal_CSTR(BzzVector &y, double t, BzzVector &dy)
{
    int i;

    for (i=1;i<=NC;i++)	omega[i]    = y[i];

    // User Defined Temperature profile
    if (iUserDefinedTemperature == true)
    {
        T = ud_temperature_profile.GiveMeValue(t,t);
        indexProperties = 0;
    }

	// Temperature fluctuations
    if (iSetFluctuations == true)
    {
		fluctuations.Update(t);
        indexProperties = 0;
    }

    // Properties + Reaction rates (constant temperature)
    UpdateProperties_isothermal(indexProperties);

	// Update residence time and volume
	Update(t);

	// Update total mass [kg]
	mass_total = rho*Volume;

	// Conservation equations for all the species
	for(i=1;i<=NC;i++)
		domega[i] = (inletStream->omega[i] - omega[i]) / Tau + R[i]/rho;

    // Soot Two Equation Model
    if (iTwoEquationModel == true)
        UpdateTwoEquationModel(y, dy);

    // Recovering residuals (Soot 2E are already updated)
    for (i=1;i<=NC;i++)	dy[i] = domega[i];
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
//																								   //
//						ODE SYSTEM EQUATIONS - NON ISOTHERMAL REACTOR								   //
//																								   //
/////////////////////////////////////////////////////////////////////////////////////////////////////

void OpenSMOKE_CSTR::ODESystem_NonIsothermal_CSTR(BzzVector &y, double t, BzzVector &dy)
{
    int i;

    // Recovering variables;
    for (i=1;i<=NC;i++)	omega[i]    = y[i];
    T				                = y[NC+1];

	// Updating Exchange area
	UpdateExchangeArea(t,t);
    UpdateHeatFlux(t,t);

    // Evaluation of properties
    UpdateProperties(ODE_CSTR_Object.jacobianIndex, NC+1);

	// Update residence time and volume
	Update(t);

	// Update total mass [kg]
	mass_total = rho*Volume;
	double dmass_total = 0.;

	// Conservation equations for all the species
	for(i=1;i<=NC;i++)
		domega[i] = (inletStream->omega[i] - omega[i]) / Tau + R[i]/rho;

	// Conservation equations for energy (CSTR)	
	double hout = mass_flow_out*Dot(omega,h);		
	double reaction_term = mass_total*(Dot(h,domega) - Constants::R_J_kmol*T*Dot(mix->uM,domega));
	dT  = inletStream->enthalpy - hout - reaction_term - specificInternalEnergy*dmass_total;
	dT -= Qe*Ae;
	dT /= (Cv*mass_total);

    // Soot Two Equation Model
    if (iTwoEquationModel == true)	UpdateTwoEquationModel(y, dy);

    // Recovering residuals (Soot 2E are already updated)
    for (i=1;i<=NC;i++)	dy[i] = domega[i];
    dy[NC+1] = dT;
}

void OpenSMOKE_CSTR::ODESystem_Isothermal_GasDischargeLaw_CSTR(BzzVector &y, double t, BzzVector &dy)
{
	int i;

    for (i=1;i<=NC;i++)	omega[i]    = y[i];
						mass_total	= y[NC+1];

    // User Defined Temperature profile
    if (iUserDefinedTemperature == true)
    {
        T = ud_temperature_profile.GiveMeValue(t,t);
        indexProperties = 0;
    }

    // Properties + Reaction rates (constant temperature)
    UpdateProperties_isothermal(indexProperties);

	// Update residence time, volume and mass flow rates
	Update(t);
	
	// Total mass equation
	double dmass_total = inletStream->massFlowRate - mass_flow_out;

	// Conservation equations for all the species
	for(i=1;i<=NC;i++)
		domega[i] = (inletStream->massFlowRate*inletStream->omega[i] - mass_flow_out*omega[i])/mass_total +
					R[i]/rho - omega[i]/mass_total*dmass_total;

    // Recovering residuals
    for (i=1;i<=NC;i++)	
		dy[i] = domega[i];
	dy[NC+1] = dmass_total;
}

void OpenSMOKE_CSTR::ODESystem_NonIsothermal_GasDischargeLaw_CSTR(BzzVector &y, double t, BzzVector &dy)
{
	int i;

    for (i=1;i<=NC;i++)	omega[i]    = y[i];
						T			= y[NC+1];
						mass_total	= y[NC+2];

	// Updating Exchange area
	UpdateExchangeArea(t,t);
    UpdateHeatFlux(t,t);

    // Evaluation of properties
    UpdateProperties(ODE_CSTR_Object.jacobianIndex, NC+1);

	// Update residence time, volume and mass flow rates
	Update(t);
	
	// Total mass equation
	double dmass_total = inletStream->massFlowRate - mass_flow_out;

	// Conservation equations for all the species
	for(i=1;i<=NC;i++)
		domega[i] = (inletStream->massFlowRate*inletStream->omega[i] - mass_flow_out*omega[i])/mass_total +
					R[i]/rho - omega[i]/mass_total*dmass_total;

	double hout = mass_flow_out*Dot(omega,h);		
	double reaction_term = mass_total*(Dot(h,domega) - Constants::R_J_kmol*T*Dot(mix->uM,domega));
	dT  = inletStream->enthalpy - hout - reaction_term - specificInternalEnergy*dmass_total;
	dT -= Qe*Ae;
	dT /= (Cv*mass_total);


    // Recovering residuals
    for (i=1;i<=NC;i++)	
		dy[i] = domega[i];
	dy[NC+1] = dT;
	dy[NC+2] = dmass_total;
}

void OpenSMOKE_CSTR::DefineFromFile(const std::string inputFile)
{
    double  double_value;
    std::string  string_value;
    int     int_value;
	vector<string> string_vector;
	vector<double> double_vector;

    OpenSMOKE_Dictionary_CSTR dictionary;
    dictionary.ParseFile(inputFile);


    // COMPULSORY: Energy Equation solution
    AssignEnergy(dictionary.Return("#Energy"));

    // SEMI-COMPULSORY: Reactor Volume || Reactor Residence Time
    if (dictionary.Return("#Volume", double_value, string_value))
    {
		AssignVolume(string_value, double_value);
	}	
	else if (dictionary.Return("#ResidenceTime", double_value, string_value))
    {
		AssignResidenceTime(string_value, double_value);
	}
	else if (dictionary.Return("#Diameter", double_value, string_value))
    {
		AssignDiameter(string_value, double_value);
	}

    // SEMI-COMPULSORY: Reactor Volume || Reactor Residence Time
    if (dictionary.Return("#UserDefined_ResidenceTime", string_value))
    {
		SetUserDefinedResidenceTime(string_value);
	}
    else if (dictionary.Return("#UserDefined_Volume", string_value))
    {
		SetUserDefinedVolume(string_value);
	}
    else if (dictionary.Return("#UserDefined_Diameter", string_value))
    {
		SetUserDefinedDiameter(string_value);
	}

	// OPTIONAL: Integration time
    if (dictionary.Return("#IntegrationTime", double_value, string_value))
		AssignEnd(string_value, double_value);

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
        SetConstantExchangeArea(double_value, string_value);

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
        SetUserDefinedExchangeArea(string_value);

	// OPTIONAL: UserDefined_U
    if (dictionary.Return("#UserDefined_U", string_value))
        SetUserDefinedHeatExchangeCoefficient(string_value);

	// OPTIONAL: UserDefined_Tambient
    if (dictionary.Return("#UserDefined_Tambient", string_value))
        SetUserDefinedAmbientTemperature(string_value);

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

	// OPTIONAL: Rate of Production Analysis
    if (dictionary.Return("#ElementFluxAnalysis", string_vector))
		SetElementFluxAnalysisOnFile(string_vector);

	// OPTIONAL: Local Rate of Production Analysis
    if (dictionary.Return("#LocalElementFluxAnalysis", string_vector))
		SetLocalElementFluxAnalysisOnFile(string_vector);

    if (dictionary.Return("#InitialTemperature", double_value, string_value))
        SetInitialTemperature(double_value, string_value);

    if (dictionary.Return("#InitialMassFractions", double_vector, string_vector))
		SetInitialMassFractions(string_vector, double_vector);

    if (dictionary.Return("#InitialMoleFractions", double_vector, string_vector))
		SetInitialMoleFractions(string_vector, double_vector);

	 if (dictionary.Return("#Fluctuations", string_vector))
		 SetFluctuations(string_vector);

	 if (dictionary.Return("#ValveCoefficient",  double_value, string_value))
		 SetValveCoefficient(double_value, string_value);

	 if (dictionary.Return("#PressureOutlet",  double_value, string_value))
		 SetPressureOutlet(double_value, string_value);

    if (dictionary.Return("#Experiment", string_vector))
	{
		SetHistoryPartial();
		SetExperimentOnFile(string_vector);
	}

	// OPTIONAL: OutputSpecies
	if (dictionary.Return("#OutputSpecies", string_vector))
		SetOutputSpecies(string_vector);

	Lock();
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
//																								   //
//									DICTIONARY													   //
//																								   //
/////////////////////////////////////////////////////////////////////////////////////////////////////

OpenSMOKE_Dictionary_CSTR::OpenSMOKE_Dictionary_CSTR()
{
    SetupBase();
	SetName("OpenSMOKE_CSTR Dictionary");

    Add("#Energy",        'O', 'N', "The energy equation is solved");
    Add("#NoEnergy",      'O', 'N', "The energy equation is not solved");
    
    Add("#ResidenceTime",	'O', 'M', "Reactor residence time");
	Add("#Volume",			'O', 'M', "Reactor volume");
    Add("#Diameter",		'O', 'M', "Reactor diameter");
	Add("#IntegrationTime",	'O', 'M', "Integration time");
    
	Add("#UserDefined_Diameter",	'O', 'S', "The diameter is assigned from file");
    Add("#UserDefined_Volume",		'O', 'S', "The reactor volume is assigned from file");
    Add("#UserDefined_ResidenceTime",	'O', 'S', "The reactor residence time is assigned from file");

    Add("#Key",				'O', 'V', "Key species");

	Add("#Qe",				'O', 'M', "Specific heat flux");
	Add("#Ae",				'O', 'M', "Exchange area");
	Add("#U",				'O', 'M', "Heat transfer exchange coefficient");
	Add("#Tambient",		'O', 'M', "Ambient temperature");
	Add("#OutputSpecies",   'O', 'V', "Output species");

	Add("#InitialTemperature",		'O', 'M', "Initial reactor temperature");
	Add("#InitialMassFractions",	'O', 'L', "Initial reactor mass fractions");
    Add("#InitialMoleFractions",	'O', 'L', "Initial reactor mole fractions");

	// Optional
	Add("#VerboseEnergy",	'O', 'N', "Report on energy of CSTR");

	Add("#ReactionRates",		'O', 'V', "Reaction rates on file (list of reaction indices)");
	Add("#FormationRates",		'O', 'V', "Formation rates on file (list of species names) ");
	Add("#ROPA",				'O', 'N', "Rate of Production Analysis");
	Add("#VerboseROPA",			'O', 'V', "Rate of Production Analysis (list of species names)");
	Add("#Sensitivity",			'O', 'V', "Sensitivity Analysis");
	Add("#ElementFluxAnalysis",		 'O', 'V', "Element Flux Analysis (list of element names, lower case)");
	Add("#LocalElementFluxAnalysis", 'O', 'V', "Local Element Flux Analysis (list of times)");

	Add("#UserDefined_T",			'O', 'S', "The temperature is assigned from external file");
    
	Add("#UserDefined_Qe",			'O', 'S', "The specific heat flux is assigned from external file");
    Add("#UserDefined_Ae",			'O', 'S', "The specific exchange area is assigned from external file");
	Add("#UserDefined_U",			'O', 'S', "The heat exchange coefficient is assigned from external file");
    Add("#UserDefined_Tambient",	'O', 'S', "The ambient temperature is assigned from external file");

	Add("#Fluctuations",			'O', 'V', "Fluctuations (example T 1500 K f 100 Hz semiAmplitude 200 K tWall 1 s");

	Add("#ValveCoefficient",		'O', 'M', "Valve flow coefficient");
	Add("#PressureOutlet",			'O', 'M', "PressureOutlet");
	
	Add("#RelativeTolerance",		'O', 'D', "Relative tolerance (default: 1.2e-5)");
	Add("#AbsoluteTolerance",		'O', 'D', "Absolute tolerance (default: 1.0e-10)");

	Add("#Experiment",				'O', 'V', "Experiment Analysis");

	// Conflicts
    Conflict("#Energy",					"#NoEnergy");
    Conflict("#Energy",					"#UserDefined_T");
    Conflict("#Energy",					"#Fluctuations");
    Conflict("#UserDefined_T",			"#Fluctuations");
	
    Conflict("#Volume",					"#ResidenceTime");
    Conflict("#Diameter",		    	"#ResidenceTime");
    Conflict("#Volume",					"#Diameter");
    
	Conflict("#ResidenceTime",			"#ValveCoefficient");
    Conflict("#ResidenceTime",			"#PressureOutlet");

	Conflict("#UserDefined_Diameter",		"#ValveCoefficient");
	Conflict("#UserDefined_Volume",			"#ValveCoefficient");
	Conflict("#UserDefined_ResidenceTime",	"#ValveCoefficient");

	Conflict("#UserDefined_Diameter",		"#Diameter");
	Conflict("#UserDefined_Volume",			"#Volume");
	Conflict("#UserDefined_ResidenceTime",	"#ResidenceTime");

    Conflict("#UserDefined_Diameter",	"#ResidenceTime");
    Conflict("#UserDefined_Diameter",	"#Volume");
    Conflict("#UserDefined_Diameter",	"#UserDefined_Volume");
    Conflict("#UserDefined_Diameter",	"#UserDefined_ResidenceTime");

    Conflict("#UserDefined_Volume",	"#ResidenceTime");
    Conflict("#UserDefined_Volume",	"#Diameter");
    Conflict("#UserDefined_Volume",	"#UserDefined_Diameter");

    Conflict("#UserDefined_ResidenceTime",	"#Volume");
    Conflict("#UserDefined_ResidenceTime",	"#Diameter");
    
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

	Conflict("#NoEnergy",				"#InitialTemperature");
	Conflict("#InitialMassFractions",	"#InitialMoleFractions");

	// Compulsory
	Compulsory("#Energy",   "#NoEnergy");
	
    Lock();
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
//																								   //
//									ODE SYSTEM - CLASS DEFINITION								   //
//																								   //
/////////////////////////////////////////////////////////////////////////////////////////////////////

void MyOdeSystem_CSTR::ObjectBzzPrint(void)
{
    ::BzzPrint("\nOpenSMOKE_CSTR\n");
}

void MyOdeSystem_CSTR::GetSystemFunctions(BzzVector &x, double t, BzzVector &f)
{
    if (iEnergy == false && iSetGasDischargeLaw == false)
        ptCSTR->ODESystem_Isothermal_CSTR(x, t, f);

    else if (iEnergy == true && iSetGasDischargeLaw == false)
        ptCSTR->ODESystem_NonIsothermal_CSTR(x, t, f);

    else if (iEnergy == false && iSetGasDischargeLaw == true)
        ptCSTR->ODESystem_Isothermal_GasDischargeLaw_CSTR(x, t, f);

    else if (iEnergy == true && iSetGasDischargeLaw == true)
        ptCSTR->ODESystem_NonIsothermal_GasDischargeLaw_CSTR(x, t, f);
}

void MyOdeSystem_CSTR::assignCSTR(OpenSMOKE_CSTR *_cstr, const bool _iEnergy, const bool _iSetGasDischargeLaw)
{
    ptCSTR				= _cstr;
    iEnergy				= _iEnergy;
	iSetGasDischargeLaw	= _iSetGasDischargeLaw;
}

void OpenSMOKE_CSTR::VideoFinalResult()
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

void OpenSMOKE_CSTR::VideoGeneralInfo()
{
    cout.setf(ios::scientific);
    cout << endl;

    cout << "xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx" << endl;
    cout << "xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx" << endl;
    cout << endl;

    cout << "---------------------------------------------------------------------" << endl;
    cout << " CSTR Summary" << endl;
    cout << "---------------------------------------------------------------------" << endl;
    cout << " Volume:            "		<< Volume							<< " [m3]"		<< endl;
    cout << " Residence time:    "		<< Tau								<< " [s]"		<< endl;
	cout << " Mass flow rate:    "		<< inletStream->massFlowRate		<< " [kg/s]"	<< endl;
    cout << " Mole flow rate:    "		<< inletStream->moleFlowRate		<< " [kmol/s]"	<< endl;
    cout << " Volume flow rate:  "		<< inletStream->volumetricFlowRate	<< " [m3/s]"	<< endl;
	cout << " Pressure:          "		<< P								<< " [Pa]"		<< endl;
    cout << " Density:           "		<< rho								<< " [kg/m3]"	<< endl;
    cout << " Molecular weight:  "		<< MWtot							<< " [kg/kmol]" << endl;
    cout << " Area:              "		<< Area								<< " [m2]"		<< endl;
    cout << " Temperature:       "		<< T								<< " [K]"		<< endl;
    cout <<  endl;
}

void OpenSMOKE_CSTR::VideoSummary()
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

void OpenSMOKE_CSTR::SummaryOnFile()
{
	int i;
	std::string file_name = outputFolderName + "/Summary.out";
	ofstream fSummary;
	openOutputFileAndControl(fSummary, file_name);
    fSummary.setf(ios::scientific);

    fSummary << "                  	Inlet        	  Outlet"										<< endl;
	fSummary << "----------------------------------------------------------------------------"		<< endl;
    fSummary << setw(20) << left << "Residence_time[s]"	<< setw(16) << left << ZERO							<< setw(16) << left << Tau								<< endl;
	fSummary << setw(20) << left << "Volume[m3]"		<< setw(16) << left << ZERO							<< setw(16) << left << Volume							<< endl;	
	fSummary << setw(20) << left << "Mass[kg]"			<< setw(16) << left << mass_total_initial			<< setw(16) << left << mass_total						<< endl;	
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


void OpenSMOKE_CSTR::EnergyAnalysis(OpenSMOKE_GasStream &outletStream)
{
	double H_mass_inlet  = inletStream->massSpecificEnthalpy;			// [J/kg]
	double P_mass_inlet  = inletStream->massSpecificPressureEnergy;		// [J/kg]
	double U_mass_inlet  = inletStream->massSpecificInternalEnergy;		// [J/kg]
	
	double H_mass_outlet  = outletStream.massSpecificEnthalpy;			// [J/kg]
	double P_mass_outlet  = outletStream.massSpecificPressureEnergy;	// [J/kg]
	double U_mass_outlet  = outletStream.massSpecificInternalEnergy;	// [J/kg]

	cout << endl;
	cout << "---------------------------------------------------------------------" << endl;
	cout << "                         ENERGY ANALYSIS                             " << endl;
	cout << "---------------------------------------------------------------------" << endl;
	cout << "Inlet  - Enthalpy:       "	<< H_mass_inlet	  << " [J/kg]" << endl;
	cout << "Outlet - Enthalpy:       "	<< H_mass_outlet  << " [J/kg]" << endl;
	cout << endl;

	cout << "Inlet power:     "			<< H_mass_inlet*massFlowRate					<< " [W]" << endl;
	cout << "Outlet power:    "			<< -H_mass_outlet*massFlowRate					<< " [W]" << endl;
	cout << "Exchanged power: "			<< -(H_mass_outlet-H_mass_inlet)*massFlowRate	<< " [W]" << endl;
}

void OpenSMOKE_CSTR::ODEPrint(BzzVector &y, double t)
{
	double dummy = ZERO;

    if (iVerbose == true)
    {
        if ( countIterations%(300*nVideoSteps) == 0)
        {
            cout    << endl
                    << "#"			<< "\t"
                    << "Tau[s]"		<< "\t\t"
                    << "T[K]"       << "\t\t"
                    << "P[atm]"     << "\t\t";
                    
			cout	<< endl;
        }

        if (countVideoSteps == nVideoSteps)
        {
            cout	<< countIterations	<< "\t"
                    << t  				<< "\t"
                    << T				<< "\t"
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

            fOutput << setw(20) << left << t
                    << setw(20) << left << Volume
					<< setw(20) << left << T
					<< setw(20) << left << Tau
                    << setw(20) << left << P/101325.
					<< setw(20) << left << Qe
                    << setw(20) << left << QReaction
                    << setw(20) << left << Qe*Ae
                    << setw(20) << left << rho
                    << setw(20) << left << MWtot
                    << setw(20) << left << D 
                    << setw(20) << left << Ae
					<< setw(20) << left << mass_total
					<< setw(20) << left << mass_flow_out;

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
				fSoot2E	<< setw(20) << left << t
						<< setw(20) << left << Volume
						<< setw(20) << left << T;

				soot2EModel->write_on_file(fSoot2E, y[indexTwoEquations], y[indexTwoEquations+1]);
			}

 /*           if (mix->soot_manager.iSoot == true)
            {
                // Soot Distribution data
                mix->soot_manager.large_calculateAll(omega, x, c, rho);
                mix->soot_manager.small_calculateAll(omega, x, c, rho);

				fSoot	<< setw(20) << left << t
						<< setw(20) << left << Volume
						<< setw(20) << left << T;
                
				mix->soot_manager.print_main_data_on_file(fSoot);
                fSoot << endl;
            }
			*/
/*            if (mix->pah_manager.iPAH == true)
            {
                // PAH Distribution data
                mix->pah_manager.pah_calculateAll(omega, x, c, rho);

				fPAH	<< setw(20) << left << t
						<< setw(20) << left << Volume
						<< setw(20) << left << T;
                
				mix->pah_manager.print_main_data_on_file(fPAH);
                fPAH << endl;
            }
*/
			if (iKeySpecies == true && key_species_names.size()>1)
            {
				fConversions	<< setw(20) << left << t
								<< setw(20) << left << Volume
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
				double Ep_mass = P/rho;
				double Ek_mass = 0.;

				fEnergy << setw(20) << left << t
						<< setw(20) << left << Volume
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
				fReactionRates	<< setw(20) << left << t
								<< setw(20) << left << Volume
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
				fFormationRates	<< setw(20) << left << t
								<< setw(20) << left << Volume
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

		if (iVerboseLocalElementFluxAnalysis == true)
		{
			if (index_local_ElementFluxAnalysis.Size()>0)
				if (t >= index_local_ElementFluxAnalysis[1])
				{
					BzzVector rForward(mix->NumberOfReactions());
					BzzVector rBackward(mix->NumberOfReactions());
					mix->ComputeFromConcentrations(T, c, c.GetSumElements(), rForward, rBackward);

					stringstream tag;
					tag << index_local_ElementFluxAnalysis[1];
					std::string name_file = outputFolderName + "/FluxAnalysis_" + tag.str();
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
		
		if (iHistoryPartial == true)
		{
			Tau_History[1] = t;
			T_History[1]   = T;
			mass_History.SetRow(1, omega);
			mole_History.SetRow(1, x);
		}
		else
		{
			if(countGlobalIterations>0)
			{
				int index = countGlobalIterations;

				if (countGlobalIterations>MAX_TIME_STEPS)
					index = MAX_TIME_STEPS;    
				
				Tau_History[index] = t;
				T_History[index]   = T;
				mass_History.SetRow(index, omega);
				mole_History.SetRow(index, x);			
			}
		}
    }
}

void OpenSMOKE_CSTR::LabelODEFile()
{
    int i;

    if (iVerbose == true)
    {
        {
			int fOutputCount = 1;
			PrintTagOnGnuplotLabel(20, fOutput, "t[s]",			fOutputCount);
			PrintTagOnGnuplotLabel(20, fOutput, "V[m3]",		fOutputCount);
			PrintTagOnGnuplotLabel(20, fOutput, "T[K]",			fOutputCount);
			PrintTagOnGnuplotLabel(20, fOutput, "Tau[s]",		fOutputCount);
			PrintTagOnGnuplotLabel(20, fOutput, "P[atm]",		fOutputCount);
			PrintTagOnGnuplotLabel(20, fOutput, "Qe[W/m2]",		fOutputCount);
			PrintTagOnGnuplotLabel(20, fOutput, "QR[W/m3]",		fOutputCount);
			PrintTagOnGnuplotLabel(20, fOutput, "QE[W]",		fOutputCount);
			PrintTagOnGnuplotLabel(20, fOutput, "rho[kg/m3]",	fOutputCount);
			PrintTagOnGnuplotLabel(20, fOutput, "MW[g/mol]",	fOutputCount);
			PrintTagOnGnuplotLabel(20, fOutput, "D[m]",			fOutputCount);
			PrintTagOnGnuplotLabel(20, fOutput, "A[m2]",		fOutputCount);
			PrintTagOnGnuplotLabel(20, fOutput, "m[kg]",		fOutputCount);
			PrintTagOnGnuplotLabel(20, fOutput, "mout[kg/s]",		fOutputCount);

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
            fSoot2E	<< setw(20) << left << "t[s](1)"
					<< setw(20) << left << "V[m3](2)"
					<< setw(20) << left << "T[K](3)";

            soot2EModel->GnuPlotInterface(fSoot2E, 4);
			fSoot2E << endl << endl;
		}
/*
        if (mix->pah_manager.iPAH == true)
		{
            fPAH	<< setw(20) << left << "t[s](1)"
					<< setw(20) << left << "V[m3](2)"
					<< setw(20) << left << "T[K](3)";

			mix->pah_manager.GnuPlotInterface(fPAH, 4);
			fPAH << endl << endl;
		}
		*/
/*        if (mix->soot_manager.iSoot == true)
        {
            fSoot	<< setw(20) << left << "t[s](1)"
					<< setw(20) << left << "V[m3](2)"
					<< setw(20) << left << "T[K](3)";

            mix->soot_manager.GnuPlotInterface(fSoot, 4);
            fSoot << endl << endl;
        }
*/
		if (iKeySpecies == true && key_species_names.size()>1)
        {
            fConversions	<< setw(20) << left << "t[s](1)"
							<< setw(20) << left << "V[m3](2)"
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
			PrintTagOnGnuplotLabel(20, fEnergy, "t[s]",			fOutputCount);
			PrintTagOnGnuplotLabel(20, fEnergy, "V[m3]",		fOutputCount);
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

			fReactionRates	<< setw(20) << left << "t[s](1)"
							<< setw(20) << left << "V[m3](2)"
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
			fFormationRates	<< setw(20) << left << "t[s](1)"
							<< setw(20) << left << "V[m3](2)"
							<< setw(20) << left << "T[K](3)";
			
			int fOutputCount = 4;
			for (i=1;i<=index_formation_rates.Size();i++)
				PrintTagOnGnuplotLabel(20, fFormationRates, names_formation_rates[i-1],	fOutputCount);

			fFormationRates << endl;
			fFormationRates << endl;
		}
    }
}

void OpenSMOKE_CSTR::SetHistoryPartial()
{
	SetHistory();
	iHistoryPartial = true;
}

void OpenSMOKE_CSTR::ExperimentOnFile()
{
	ofstream fExperiment;
	openOutputFileAndControl(fExperiment, outputExperimentName);
	fExperiment.setf(ios::scientific);

	fExperiment << "REACTOR CSTR" << endl;
		
	fExperiment << "TIME   s" << endl;
	fExperiment << "DATA    "   << 1+names_Experiment.size() << endl;
	fExperiment << endl;

	fExperiment << "T TEMPERATURE 1.0 CROSS" << endl;
	fExperiment << 0. << "\t" << T_History[1] << endl;
	fExperiment << 1. << "\t" << T_History[1] << endl;
	fExperiment << "//" << endl << endl;

	for(int k=0;k<int(names_Experiment.size());k++)
	{
		fExperiment << names_Experiment[k] << " MOLE_FRACTION 1.0 CROSS" << endl;
		fExperiment << 0. << "\t" << mole_History[1][index_Experiment[k+1]] << endl;
		fExperiment << 1. << "\t" << mole_History[1][index_Experiment[k+1]] << endl;
		fExperiment << "//" << endl << endl;
	}

	fExperiment.close();
}

// ------------------------------------------------------------------------------ //
// 											UNSTEADY VARIABLES																				//
// ------------------------------------------------------------------------------ //

OpenSMOKE_Fluctuations::OpenSMOKE_Fluctuations() : pointer(NULL)
{
	iKind = 1;	// Sin fluctuations
	pointer = NULL;
}

void OpenSMOKE_Fluctuations::Setup(const double _mean, const double _frequency, const double _semiAmplitude, const double _tWall)
{	
	mean = _mean;
	frequency = _frequency;
	semiAmplitude = _semiAmplitude;
	tWall = _tWall;						// wall time [s]	
	
	alfa = -log(0.01) * frequency;		// damping factor
	period = Constants::_2pi/frequency;	// period [s]
}

void OpenSMOKE_Fluctuations::SetPointer(double &_pointer)
{
	pointer = &_pointer;
}

void OpenSMOKE_Fluctuations::Update(const double t)
{
	double value = mean;

	if (iKind == 1)
	{
		if (t>tWall)
		{
			double damping =  1.-exp(-alfa*(t-tWall));
			double sin_function = sin( Constants::_2pi*frequency*(t-tWall) );
			value = mean + damping*semiAmplitude*sin_function;
		}
	}

	*pointer = value;
}

void OpenSMOKE_CSTR::SaveOnBinaryFile(BzzSave &fOutput)
{
	std::string dummy;
	char name[Constants::NAME_SIZE];

	int nSteps = countGlobalIterations;			
	BzzVector timegrid = GetBzzVector(nSteps, 1, Tau_History);
	BzzVector tgrid    = GetBzzVector(nSteps, 1, T_History);

	dummy = "CSTR";
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

void OpenSMOKE_CSTR::SaveOnBinaryFile(const std::string filename)
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

void OpenSMOKE_CSTR::UpdateReactionRates(const double TT, const double PP, BzzVector &xx, BzzVector &rr)
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
