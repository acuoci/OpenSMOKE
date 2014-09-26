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
#include "idealreactors/batch/OpenSMOKE_Batch.h"
#include "basic/OpenSMOKE_Conversions.h"
//#include "addons/OpenSMOKE_PostProcessor.h"

OpenSMOKE_Batch	*ptBatch;

void ODEPrintExternal_Batch(BzzVector &y, double csi)
{
    ptBatch->ODEPrint(y, csi);
}

OpenSMOKE_Batch::OpenSMOKE_Batch()
{
    ptBatch		=  this;
	class_name	= "OpenSMOKE_Batch";
	out_name	= "Batch.out";
	
	Setup();
	
	assignedVolume				= false;
	iConstantPressure			= false;	// 0 = Constant volume (default) 1=Constant Pressure
	iVolumeLaw					= false;
}

void OpenSMOKE_Batch::AssignConstantPressure()
{
	iConstantPressure = true;
}

void OpenSMOKE_Batch::AssignEnd(const std::string units, const double value)
{
	TauTotal	= OpenSMOKE_Conversions::conversion_time(value, units);
    assignedEnd = true;
}

void OpenSMOKE_Batch::AssignVolume(const std::string units, const double value)
{
	volume			= OpenSMOKE_Conversions::conversion_volume(value, units);
	volumeInitial   = volume;
    assignedVolume	= true;
}

void OpenSMOKE_Batch::SetConstantExchangeArea(const double value, const std::string units)
{
	iUserDefinedExchangeArea = CONSTANT;
	A = OpenSMOKE_Conversions::conversion_area(value, units);
}

void OpenSMOKE_Batch::SetVolumeLaw(const double value, const std::string units)
{
	iVolumeLaw = true;
	alfaVolumeLaw = OpenSMOKE_Conversions::conversion_frequency(value, units);
}

void OpenSMOKE_Batch::SetUserDefinedExchangeArea(const std::string fileName)
{
    iUserDefinedExchangeArea = USERDEFINED;
    ud_A_profile.AssignFromFile(fileName, "AREA");
	ud_A_profile.SetName(name_object + " - Exchange Area Profile");
}

void OpenSMOKE_Batch::UnsetUserDefinedExchangeArea()
{
    iUserDefinedExchangeArea = NONE;
}

void OpenSMOKE_Batch::Lock()
{
    if (assignedKineticScheme == false)
        ErrorMessage("The kinetic scheme was not defined!");
    if (assignedGlobalKineticScheme == false)
        ErrorMessage("The global kinetic scheme was not defined!");
    if (assignedEnd == false)
		ErrorMessage("The reactor residence time was not defined!");
    if (assignedInletFlows == false)
        ErrorMessage("The inlet flow was not defined!");
    if (assignedEnergy == false)
        ErrorMessage("The kind of reactor was not defined!");
	if (assignedSoot2EModel == false)
        ErrorMessage("The 2E Model was not defined!");
	if (iTwoEquationModel == true && mix->iSootMode == false)
		ErrorMessage("This kinetic scheme was interpreted without the SootMode option and cannot be used together the 2E model.\nMore details on the OpenSMOKE's User Guide.");


	Initialize();
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
//																								   //
//									GAS MIXTURE PROPERTIES										   //
//																								   //
/////////////////////////////////////////////////////////////////////////////////////////////////////

void OpenSMOKE_Batch::UpdateProperties_isothermal(int memoIndex)
{
    double cTot;

    mix->GetMWAndMoleFractionsFromMassFractions(MWtot, x, omega);

    // --------------------------------------------------------------------------
    // PROPERTIES FOR CONSTANT T + MEMORIZATION
    // --------------------------------------------------------------------------
    if (memoIndex==0)
    {
        // Kinetic parameters (the enthalpies are calculated in this step)
        // ----------------------------------------------------------------------
		rho		= mass / volume;							// [kg/m3]
		cTot	= rho  / MWtot;								// [kmol/m3]
		P		= cTot * (Constants::R_J_kmol*T);			// [Pa]
		mix->ComputeKineticParameters(T, log(T), 1./T, P);	

        memoIndex = 1;
    }

    // --------------------------------------------------------------------------
    // PROPERTIES FOR CONSTANT T (Recover from tables)
    // --------------------------------------------------------------------------
    if (memoIndex == 1)
    {
        // a. Calcolo della concentrazione e della densita [kmol/m3] [kg/m3]
		if (iConstantPressure == false)	// Constant volume
		{
			moles	= mass / MWtot;						// [kmol]
			rho		= mass / volume;					// [kg/m3]
			cTot	= rho  / MWtot;						// [kmol/m3]
			P		= cTot * (Constants::R_J_kmol*T);	// [Pa]
			c       = cTot * x;							// [kmol/m3]	
		}
		else
		{
			moles	= mass  / MWtot;					// [kmol]
			volume	= moles * Constants::R_J_kmol*T/P;	// [m3]
			rho		= mass  / volume;					// [kg/m3]
			cTot	= rho   / MWtot;					// [kmol/m3]
			c       = cTot  * x;						// [kmol/m3]
		}			

        // c. Reaction Rates
        mix->ComputeFromConcentrations( T, c, cTot, &R);		// [kmol/m3/s]
        ElementByElementProduct(R, mix->M, &R);					// [kg/m3/s]
    }

    // Global Kinetics
    if (iGlobalKinetics==true)
        global->GiveMeFormationRates(T,c,R);					// [kg/m3/s]
}

void OpenSMOKE_Batch::UpdateProperties(const int jacobianIndex, const int indexT)
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
		if (iConstantPressure == false)	// Constant volume
		{
			moles	= mass / MWtot;						// [kmol]
			rho		= mass / volume;					// [kg/m3]
			cTot	= rho  / MWtot;						// [kmol/m3]
			P		= cTot * (Constants::R_J_kmol*T);	// [Pa]
			c       = cTot * x;							// [kmol/m3]
		}
		else
		{
			moles	= mass  / MWtot;					// [kmol]
			volume	= moles * Constants::R_J_kmol*T/P;	// [m3]
			rho		= mass  / volume;					// [kg/m3]
			cTot	= rho   / MWtot;					// [kmol/m3]
			c       = cTot  * x;						// [kmol/m3]
		}

		// Specific Heat [J/kgK]
		// ----------------------------------------------------------------------
		mix->SpeciesCp(T);
		Cp = mix->MixCp_FromMassFractions(omega);
		Cv = Cp-Constants::R_J_kmol/MWtot;

		// Table Creation (Memorization)
		// ----------------------------------------------------------------------
		mix->ComputeKineticParameters( T, log(T), 1./T, P);
		mix->ComputeFromConcentrations( T, c, cTot, &R);		// [kmol/m3/s]
		ElementByElementProduct(R, mix->M, &R);					// [kg/m3/s]
		QReaction = - mix->ComputeQReaction(T);					// [J/m3.s]

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
		// Specific Heat [J/kgK]
		// ----------------------------------------------------------------------
		mix->SpeciesCp(T);
		CpMap = mix->Cp;

		// Table Creation (Memorization)
		// ----------------------------------------------------------------------
		mix->ComputeKineticParameters(T, log(T), 1./T, P);
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
        // a. Calcolo della concentrazione e della densita [kmol/m3] [kg/m3]
		if (iConstantPressure == false)	// Constant volume
		{
			moles	= mass / MWtot;						// [kmol]
			rho		= mass / volume;					// [kg/m3]
			cTot	= rho  / MWtot;						// [kmol/m3]
			P		= cTot * (Constants::R_J_kmol*T);	// [Pa]
			c       = cTot * x;							// [kmol/m3]
			
		}
		else							// Constant Pressure
		{
			moles	= mass  / MWtot;					// [kmol]
			volume	= moles * Constants::R_J_kmol*T/P;	// [m3]
			rho		= mass  / volume;					// [kg/m3]
			cTot	= rho   / MWtot;					// [kmol/m3]
			c       = cTot  * x;						// [kmol/m3]
		}

		// Specific Heat [J/kgK]
		// ----------------------------------------------------------------------
		mix->Cp = CpMap;
		Cp = mix->MixCp_FromMassFractions(omega);
		Cv = Cp-Constants::R_J_kmol/MWtot;

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

/////////////////////////////////////////////////////////////////////////////////////////////////////
//																								   //
//									UPDATING FUNCTIONS											   //
//																								   //
/////////////////////////////////////////////////////////////////////////////////////////////////////

void OpenSMOKE_Batch::UpdateTwoEquationModel(BzzVector &y, BzzVector &dy)
{
    soot2EModel->update(T, P/101325., rho, 0., x, y[indexTwoEquations], y[indexTwoEquations+1]);
    soot2EModel->formation_rates();

    for (int i=1;i<=NC-1;i++)
        domega[i] += soot2EModel->SGas[i] / rho;		// Gas species
	domega[NC] += soot2EModel->S / rho;					// Soot

    dy[indexTwoEquations]	= soot2EModel->s / rho;
    dy[indexTwoEquations+1] = soot2EModel->S / rho;
}


void OpenSMOKE_Batch::UpdateHeatFlux(const double tau, const double csi)
{
	if (iUserDefinedHeatExchangeCoefficient == USERDEFINED)		U			= ud_U_profile.GiveMeValue(tau, MINUSONE);
	if (iUserDefinedAmbientTemperature		== USERDEFINED)		Tambient	= ud_Tambient_profile.GiveMeValue(tau, MINUSONE);
	if (iUserDefinedHeatFlux				== USERDEFINED)		Qe			= ud_Qe_profile.GiveMeValue(tau, MINUSONE);
	
	if (iUserDefinedHeatExchangeCoefficient != NONE)			Qe = U*(T-Tambient);
}

void OpenSMOKE_Batch::UpdateExchangeArea(const double tau, const double csi)
{
	if (iUserDefinedExchangeArea == NONE)			A = pow(36.*Constants::pi*volume*volume, 1./3.);
	if (iUserDefinedExchangeArea == USERDEFINED)	A = ud_A_profile.GiveMeValue(tau, MINUSONE);
}


/////////////////////////////////////////////////////////////////////////////////////////////////////
//																								   //
//									SOLVING FUNCTIONS											   //
//																								   //
/////////////////////////////////////////////////////////////////////////////////////////////////////

void OpenSMOKE_Batch::Initialize()
{
	ChangeDimensions(NC, &h);
	ChangeDimensions(NC, &hMap);

    T            = inletStream->T;
	P            = inletStream->P;
    omega        = inletStream->omega;
	mass         = inletStream->rho*volume;

    if (iUserDefinedTemperature == true)
        ud_temperature_profile.Check(0., T);

    if (iEnergy == false)					// Isothermal reactor
    {
        UpdateProperties_isothermal(0);
        indexProperties = 1;
    }
    else if (iEnergy == true)				// Non-isothermal reactor
        UpdateProperties(MINUSONE, NC+1);
}

void OpenSMOKE_Batch::Solve()
{
    ptBatch = this;

    double timeStart, timeEnd;
    BzzVector xMin;
    BzzVector xMax;
    BzzVector xInitial;

	PrepareFiles();


    ODE_Batch_Object.assignBatch(this, iEnergy, iConstantPressure);

    // 1.A Isothermal Batch
    if (iEnergy == false && iTwoEquationModel == false)
    {
        ChangeDimensions(NC, &xMin);  xMin=ZERO;					
        ChangeDimensions(NC, &xMax);  xMax=ONE;					// Mass fractions

        xInitial = omega;
    }

    // 1.C Isothermal Batch	+ Two Equation Model
    if (iEnergy == false && iTwoEquationModel == true)
    {
        soot2EModel->initial_values(rho);
        indexTwoEquations = NC+1;

        ChangeDimensions(NC+2, &xMin);
        xMin=ZERO;
        
		ChangeDimensions(NC+2, &xMax);
        xMax=ONE;

        xInitial = omega;
        xInitial.Append(soot2EModel->phiNStart);
        xInitial.Append(soot2EModel->phiMStart);
    }

    // 2.A NonIsothermal Batch
    if (iEnergy == true && iTwoEquationModel == false)
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

    // 2.C NonIsothermal Batch + TwoEquations
    if (iEnergy == true && iTwoEquationModel == true)
    {
        soot2EModel->initial_values(rho);
        indexTwoEquations = NC+2;

        ChangeDimensions(NC+3, &xMin);
        xMin=ZERO;
        xMin[NC+1] = 200.;
        
		ChangeDimensions(NC+3, &xMax);
        xMax=ONE;
        xMax[NC+1] = 6000.;

        xInitial = omega;
        xInitial.Append(T);
        xInitial.Append(soot2EModel->phiNStart);
        xInitial.Append(soot2EModel->phiMStart);
    }


    o(xInitial, 0., &ODE_Batch_Object);

    o.StepPrint(ODEPrintExternal_Batch);
    o.SetMinimumConstraints(xMin);
    o.SetMaximumConstraints(xMax);

    if (iRelativeTolerance == true)	o.SetTolRel(relativeTolerance);
    if (iAbsoluteTolerance == true)	o.SetTolAbs(absoluteTolerance);

    {
        countGlobalIterations = -1;  // From -1 to avoid to store results from the first iteration

        timeStart = BzzGetCpuTime();
		
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

	if (iHistory == true)
	{
		SaveOnBinaryFile(outputOSMName);
	}

	CloseFiles();
}

void OpenSMOKE_Batch::DefineFromFile(const std::string inputFile)
{
    double			double_value;
    std::string			string_value;
    int				int_value;
	vector<string>  string_vector;

    OpenSMOKE_Dictionary_Batch dictionary;
    dictionary.ParseFile(inputFile);


    // COMPULSORY: Energy Equation solution
	AssignEnergy(dictionary.Return("#Energy"));

    // COMPULSORY: Reactor Residence Time
	if (dictionary.Return("#ResidenceTime", double_value, string_value))
		AssignEnd(string_value, double_value);

    // COMPULSORY: Reactor Volume
    if (dictionary.Return("#Volume", double_value, string_value))
        AssignVolume(string_value, double_value);

	// COMPULSORY: Constant Volume
	if (dictionary.Return("#ConstantVolume"))
	{}

	// COMPULSORY: Constant Volume
	if (dictionary.Return("#ConstantPressure"))
		AssignConstantPressure();

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

	// OPTIONAL: Heat flux
    if (dictionary.Return("#Qe", double_value, string_value))
        SetConstantHeatFlux(double_value, string_value);

	// OPTIONAL: Exchange area
    if (dictionary.Return("#Ae", double_value, string_value))
        SetConstantExchangeArea(double_value, string_value);

	// OPTIONAL: Exchange coefficient
    if (dictionary.Return("#U", double_value, string_value))
        SetConstantHeatExchangeCoefficient(double_value, string_value);

	// OPTIONAL: Ambient temperature
    if (dictionary.Return("#Tambient", double_value, string_value))
        SetConstantAmbientTemperature(double_value, string_value);

	// OPTIONAL: VolumeLaw
    if (dictionary.Return("#VolumeLaw", double_value, string_value))
		SetVolumeLaw(double_value, string_value);

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

/////////////////////////////////////////////////////////////////////////////////////////////////////
//																								   //
//									OUTPUT FUNCTIONS											   //
//																								   //
/////////////////////////////////////////////////////////////////////////////////////////////////////

void OpenSMOKE_Batch::LabelODEFile()
{
    int i;

    if (iVerbose == true)
    {
        {
			int fOutputCount = 1;
			PrintTagOnGnuplotLabel(20, fOutput, "Tau[s]",		fOutputCount);
			PrintTagOnGnuplotLabel(20, fOutput, "dummy",		fOutputCount);
			PrintTagOnGnuplotLabel(20, fOutput, "T[K]",			fOutputCount);
			PrintTagOnGnuplotLabel(20, fOutput, "V[l]",			fOutputCount);
			PrintTagOnGnuplotLabel(20, fOutput, "P[atm]",		fOutputCount);
			PrintTagOnGnuplotLabel(20, fOutput, "Tambient[K]",	fOutputCount);
			PrintTagOnGnuplotLabel(20, fOutput, "QR[W/m3]",		fOutputCount);
			PrintTagOnGnuplotLabel(20, fOutput, "QE[W/m3]",		fOutputCount);
			PrintTagOnGnuplotLabel(20, fOutput, "rho[kg/m3]",	fOutputCount);
			PrintTagOnGnuplotLabel(20, fOutput, "MW",			fOutputCount);
			PrintTagOnGnuplotLabel(20, fOutput, "h[W/m2/K]",	fOutputCount);
			PrintTagOnGnuplotLabel(20, fOutput, "A[m2]",		fOutputCount);
			PrintTagOnGnuplotLabel(20, fOutput, "H[J/kg]",		fOutputCount);
			PrintTagOnGnuplotLabel(20, fOutput, "U[J/kg]",		fOutputCount);

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
					<< setw(20) << left << "dummy(2)"
					<< setw(20) << left << "T[K](3)";

            soot2EModel->GnuPlotInterface(fSoot2E, 4);
			fSoot2E << endl << endl;
		}
/*
        if (mix->pah_manager.iPAH == true)
		{
            fPAH	<< setw(20) << left << "Tau[s](1)"
					<< setw(20) << left << "dummy(2)"
					<< setw(20) << left << "T[K](3)";

			mix->pah_manager.GnuPlotInterface(fPAH, 4);
			fPAH << endl << endl;
		}
*/
       /* if (mix->soot_manager.iSoot == true)
        {
            fSoot	<< setw(20) << left << "Tau[s](1)"
					<< setw(20) << left << "dummy(2)"
					<< setw(20) << left << "T[K](3)";

            mix->soot_manager.GnuPlotInterface(fSoot, 4);
            fSoot << endl << endl;
        }
		*/
		if (iKeySpecies == true && key_species_names.size()>1)
        {
            fConversions	<< setw(20) << left << "Tau[s](1)"
							<< setw(20) << left << "dummy(2)"
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
			PrintTagOnGnuplotLabel(20, fEnergy, "dummy",			fOutputCount);
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
							<< setw(20) << left << "dummy(2)"
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
							<< setw(20) << left << "dummy(2)"
							<< setw(20) << left << "T[K](3)";
			
			int fOutputCount = 4;
			for (i=1;i<=index_formation_rates.Size();i++)
				PrintTagOnGnuplotLabel(20, fFormationRates, names_formation_rates[i-1],	fOutputCount);

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

void OpenSMOKE_Batch::VideoGeneralInfo()
{
    cout.setf(ios::scientific);
    cout << endl;

    cout << "xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx" << endl;
    cout << "xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx" << endl;
    cout << endl;

    cout << "---------------------------------------------------------------------" << endl;
    cout << "Batch Summary" << endl;
    cout << "---------------------------------------------------------------------" << endl;
    cout << "Residence time:    "	<< TauTotal		<< " [s]"		<< endl;
	cout << "Pressure:          "	<< P			<< " [Pa]"		<< endl;
    cout << "Density:           "	<< rho			<< " [kg/m3]"	<< endl;
    cout << "Molecular weight:  "	<< MWtot		<< " [kg/kmol]" << endl;
    cout << "Volume:            "	<< volume		<< " [m3]"		<< endl;
    cout << "Mass:              "	<< mass			<< " [kg]"		<< endl;
    cout << "Moles:             "	<< moles		<< " [kmol]"	<< endl;
	cout << "Temperature:       "	<< T			<< " [K]"		<< endl;
    cout <<  endl;
}

void OpenSMOKE_Batch::VideoSummary()
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

void OpenSMOKE_Batch::VideoFinalResult()
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

void OpenSMOKE_Batch::SummaryOnFile()
{
	int i;

	ofstream fSummary;
	std::string file_name = outputFolderName + "/Summary.out";
	openOutputFileAndControl(fSummary, file_name);
    fSummary.setf(ios::scientific);

	double initial_volume	= mass / inletStream->rho;	// [m3]
	double initial_moles	= mass / inletStream->MW ;	// [kmol]

    fSummary << "                  	Inlet        	  Outlet"										<< endl;
	fSummary << "----------------------------------------------------------------------------"		<< endl;
    fSummary << setw(20) << left << "Time[s]"			<< setw(16) << left << ZERO					<< setw(16) << left << TauTotal								<< endl;
	fSummary << setw(20) << left << "Volume[m3]"		<< setw(16) << left << initial_volume		<< setw(16) << left << volume								<< endl;	
	fSummary << setw(20) << left << "Temperature[K]"	<< setw(16) << left << inletStream->T		<< setw(16) << left << T								<< endl;
	fSummary << setw(20) << left << "Pressure[Pa]"		<< setw(16) << left << inletStream->P		<< setw(16) << left << P								<< endl;
	fSummary << setw(20) << left << "Mass[kg]"			<< setw(16) << left << mass					<< setw(16) << left << mass		<< endl;
	fSummary << setw(20) << left << "Moles[kmol]"		<< setw(16) << left << initial_moles		<< setw(16) << left << moles	<< endl;	
	fSummary << setw(20) << left << "Density[kg/m3]"	<< setw(16) << left << inletStream->rho		<< setw(16) << left << rho								<< endl;
	fSummary << setw(20) << left << "MW[kg/kmol]"		<< setw(16) << left << inletStream->MW		<< setw(16) << left << MWtot							<< endl;
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
}


void OpenSMOKE_Batch::ODEPrint(BzzVector &y, double t)
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
					<< "P[atm]"		<< "\t\t"
					<< "V[l]"		<< "\t\t"
					<< "n[kmol]"	<< "\t\t";
                    
			cout	<< endl;
        }

        if (countVideoSteps == nVideoSteps)
        {
            cout	<< countIterations	<< "\t"
                    << t	  			<< "\t"
                    << T				<< "\t"
                    << P/101325.		<< "\t"
                    << volume/1000.		<< "\t"
                    << moles			<< "\t";

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

            fOutput << setw(20) << left << t
                    << setw(20) << left << dummy
                    << setw(20) << left << T
                    << setw(20) << left << volume/1000.
                    << setw(20) << left << P/101325.
                    << setw(20) << left << Tambient
                    << setw(20) << left << QReaction
                    << setw(20) << left << Qe*A/volume
                    << setw(20) << left << rho
                    << setw(20) << left << MWtot
                    << setw(20) << left << U
                    << setw(20) << left << A
					<< setw(20) << left << specificEnthalpy
					<< setw(20) << left << specificInternalEnergy;

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
						<< setw(20) << left << dummy
						<< setw(20) << left << T;

				soot2EModel->write_on_file(fSoot2E, y[indexTwoEquations], y[indexTwoEquations+1]);
			}

       /*     if (mix->soot_manager.iSoot == true)
            {
                // Soot Distribution data
                mix->soot_manager.large_calculateAll(omega, x, c, rho);
                mix->soot_manager.small_calculateAll(omega, x, c, rho);

                fSoot << setw(20) << left << t
					  << setw(20) << left << dummy
					  << setw(20) << left << T;
                mix->soot_manager.print_main_data_on_file(fSoot);
                fSoot << endl;
            }*/
/*
            if (mix->pah_manager.iPAH == true)
            {
                // PAH Distribution data
                mix->pah_manager.pah_calculateAll(omega, x, c, rho);

                fPAH << setw(20) << left << t
					 << setw(20) << left << dummy
					 << setw(20) << left << T;
                mix->pah_manager.print_main_data_on_file(fPAH);
                fPAH << endl;
            }
			*/
			if (iKeySpecies == true && key_species_names.size()>1)
            {
                fConversions << setw(20) << left << t
							 << setw(20) << left << dummy
							 << setw(20) << left << T;

				for(i=0;i<int(key_species_names.size());i++)
	                fConversions << setw(20) << left << 1.-omega[key_species_index[i]]/inletStream->omega[key_species_index[i]];
	        
                fConversions << endl;
            }

			if (iVerboseEnergy == true)
			{
				BzzVector h_mass(NC);
				mix->GetMixAveragedEnthalpy_Mass(h_mass, T);
	
				double H_mass  = Dot(omega, h_mass);			// [J/kg]
				double Ek_mass = 0;
				double Ep_mass = P/rho;

				fEnergy << setw(20) << left << t
						<< setw(20) << left << dummy
						<< setw(20) << left << H_mass+Ep_mass
						<< setw(20) << left << H_mass
						<< setw(20) << left << Ek_mass
						<< setw(20) << left << Ep_mass
						<< setw(20) << left << QReaction
						<< setw(20) << left << Qe*A
						<< setw(20) << left << Cp
						<< endl;
			}

			if (iVerboseReactionRates == true)
			{
				fReactionRates	<< setw(20) << left << t
								<< setw(20) << left << dummy
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
								<< setw(20) << left << dummy
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
			
			Tau_History[index] = t;
			T_History[index]   = T;
			mass_History.SetRow(index, omega);
			mole_History.SetRow(index, x);			
        }
    }
}

void OpenSMOKE_Batch::EnergyAnalysis(OpenSMOKE_GasStream &outletStream)
{
	double H_mass_initial  = inletStream->massSpecificEnthalpy;			// [J/kg]
	double P_mass_initial  = inletStream->massSpecificPressureEnergy;	// [J/kg]
	double U_mass_initial  = inletStream->massSpecificInternalEnergy;	// [J/kg]
	
	double H_mass_final  = outletStream.massSpecificEnthalpy;			// [J/kg]
	double P_mass_final  = outletStream.massSpecificPressureEnergy;		// [J/kg]
	double U_mass_final  = outletStream.massSpecificInternalEnergy;		// [J/kg]

	double Work;	// Work done on the system (IN)
	if (iConstantPressure == true)	Work = -P*(volume-mass/inletStream->rho);
	else							Work = 0.;

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
	cout << "Work on the system: "	<< Work												<< " [J]" << endl;
	cout << "Balance (IN-OUT):   "  << Work + U_mass_initial*mass - U_mass_final*mass	<< " [J]" << endl;
}


/////////////////////////////////////////////////////////////////////////////////////////////////////
//																								   //
//						ODE SYSTEM EQUATIONS - ISOTHERMAL REACTOR								   //
//																								   //
/////////////////////////////////////////////////////////////////////////////////////////////////////

void OpenSMOKE_Batch::ODESystem_Isothermal_Batch(BzzVector &y, double t, BzzVector &dy)
{
	// Recover variables
    omega    = y;

    // User Defined Temperature profile
    if (iUserDefinedTemperature == true)
    {
        T = ud_temperature_profile.GiveMeValue(t, MINUSONE);
        indexProperties = 0;
    }

    // Properties + Reaction rates (constant temperature)
    UpdateProperties_isothermal(indexProperties);

    // Conservation equations for all the species (Batch)
    domega = R / rho;

    // Soot Two Equation Model
    if (iTwoEquationModel == true)
        UpdateTwoEquationModel(y, dy);

    // Recovering residuals (Soot 2E are already updated)
    dy   = domega;
}



/////////////////////////////////////////////////////////////////////////////////////////////////////
//																								   //
//						ODE SYSTEM EQUATIONS - NON ISOTHERMAL REACTOR - CONST P					   //
//																								   //
/////////////////////////////////////////////////////////////////////////////////////////////////////

void OpenSMOKE_Batch::ODESystem_NonIsothermal_ConstantPressure_Batch(BzzVector &y, double t, BzzVector &dy)
{
    int i;

    // Recovering variables;
    for (i=1;i<=NC;i++)	omega[i]    = y[i];
    T				                = y[NC+1];

	// Updating Exchange area
	UpdateExchangeArea(t, MINUSONE);
    UpdateHeatFlux(t, MINUSONE);

    // Evaluation of properties
    UpdateProperties(ODE_Batch_Object.jacobianIndex, NC+1);

    // Conservation equations for all the species (Batch)
    domega	= R / rho ;											// [1/s]

    // Conservation equations for the energy (Batch)
	dT  = (volume*QReaction - A*Qe) / (mass*Cp);				// [K/s]
 
    // Soot Two Equation Model
    if (iTwoEquationModel == true)	UpdateTwoEquationModel(y, dy);

    // Recovering residuals (Soot 2E are already updated)
    for (i=1;i<=NC;i++)	dy[i]	= domega[i];
    dy[NC+1]					= dT;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
//																								   //
//						ODE SYSTEM EQUATIONS - NON ISOTHERMAL REACTOR - CONST V					   //
//																								   //
/////////////////////////////////////////////////////////////////////////////////////////////////////

void OpenSMOKE_Batch::ODESystem_NonIsothermal_ConstantVolume_Batch(BzzVector &y, double t, BzzVector &dy)
{
    int i;

	if (iVolumeLaw == true)
		volume = volumeInitial/(1.+alfaVolumeLaw*t);

    // Recovering variables;
    for (i=1;i<=NC;i++)	omega[i]    = y[i];
    T				                = y[NC+1];

	// Updating Exchange area
	UpdateExchangeArea(t, MINUSONE);
    UpdateHeatFlux(t, MINUSONE);

    // Evaluation of properties
    UpdateProperties(ODE_Batch_Object.jacobianIndex, NC+1);

    // Conservation equations for all the species (Batch)
    domega	= R / rho ;

	// Number of moles
	double sumRmole = 0.;
	for (i=1;i<=NC;i++)
		sumRmole += R[i]/mix->M[i];

    // Conservation equations for the energy (Batch)
    dT  = (volume*(QReaction + Constants::R_J_kmol*T*sumRmole) - A*Qe) / (mass*Cv);
 
    // Soot Two Equation Model
    if (iTwoEquationModel == true)	UpdateTwoEquationModel(y, dy);

    // Recovering residuals (Soot2E are already updated)
    for (i=1;i<=NC;i++)	dy[i]	= domega[i];
    dy[NC+1]					= dT;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
//																								   //
//											DICTIONARY											   //
//																								   //
/////////////////////////////////////////////////////////////////////////////////////////////////////

OpenSMOKE_Dictionary_Batch::OpenSMOKE_Dictionary_Batch()
{
    SetupBase();

	// Compulsory
    Add("#Energy",			'O', 'N', "The energy equation is solved");
    Add("#NoEnergy",		'O', 'N', "The energy equation is not solved");
    
	// Compulsory
    Add("#ResidenceTime",		'C', 'M', "Reactor residence time");

	// Compulsory
    Add("#ConstantVolume",		'O', 'N', "The volume is kept constant");
    Add("#ConstantPressure",	'O', 'N', "The pressure is kept constant");

	// Compulsory
    Add("#Volume",				'O', 'M', "Reactor volume");

	// Optional
    Add("#Key",					'O', 'V', "Key species");
	Add("#OutputSpecies",		'O', 'V', "Output species");

	// Optional
	Add("#ReactionRates",		'O', 'V', "Reaction rates");
	Add("#FormationRates",		'O', 'V', "Formation rates");
	Add("#ROPA",				'O', 'N', "Rate of Production Analysis");
	Add("#VerboseROPA",			'O', 'V', "Rate of Production Analysis (list of species names)");
	Add("#Sensitivity",			'O', 'V', "Sensitivity Analysis");

     
	// Optional
	Add("#VerboseEnergy",			'O', 'N', "Report on energy");
	Add("#Qe",						'O', 'M', "Heat flux");
	Add("#Ae",						'O', 'M', "Exchange area");
	Add("#U",						'O', 'M', "Heat thermal exchange coefficient");
	Add("#Tambient",				'O', 'M', "Ambient temperature");
	Add("#UserDefined_T",			'O', 'S', "The temperature is assigned from external file");
	Add("#UserDefined_Qe",			'O', 'S', "The heat flux is assigned from external file");
    Add("#UserDefined_Ae",			'O', 'S', "The exchange area is assigned from external file");
	Add("#UserDefined_U",			'O', 'S', "The heat exchange coefficient is assigned from external file");
    Add("#UserDefined_Tambient",	'O', 'S', "The ambient temperature is assigned from external file");

	Add("#RelativeTolerance",		'O', 'D', "Relative tolerance (default: 1.2e-5)");
	Add("#AbsoluteTolerance",		'O', 'D', "Absolute tolerance (default: 1.0e-10)");

	Add("#VolumeLaw",				'O', 'M', "Alfa Coefficient in volume law: V(t)=V0/(1+alfa*t)");
	
    Conflict("#Energy",					"#NoEnergy");
    Conflict("#Energy",					"#UserDefined_T");
	Conflict("#NoEnergy",				"#Qe");
	Conflict("#NoEnergy",				"#Ae");
	Conflict("#UserDefined_Qe",			"#Qe");
	Conflict("#UserDefined_Ae",			"#Ae");
	Conflict("#UserDefined_U",			"#U");
	Conflict("#UserDefined_Tambient",	"#Tambient");

	Conflict("#NoEnergy",				"#U");
	Conflict("#NoEnergy",				"#Tambient");
	Conflict("#Qe",						"#U");
	Conflict("#Qe",						"#Tambient");
	Conflict("#Qe",						"#UserDefined_U");
	Conflict("#Qe",						"#UserDefined_Tambient");
	Conflict("#UserDefined_Qe",			"#U");
	Conflict("#UserDefined_Qe",			"#Tambient");
	Conflict("#UserDefined_Qe",			"#UserDefined_U");
	Conflict("#UserDefined_Qe",			"#UserDefined_Tambient");

	Compulsory("#Energy",			"#NoEnergy");
	Compulsory("#ConstantVolume",   "#ConstantPressure");
	Conflict("#ConstantPressure",	"#VolumeLaw");
	Conflict("#NoEnergy",			"#VolumeLaw");

    Lock();
}

void OpenSMOKE_Batch::SaveOnBinaryFile(BzzSave &fOutput)
{
	std::string dummy;
	char name[Constants::NAME_SIZE];

	int nSteps = countGlobalIterations;			
	BzzVector timegrid = GetBzzVector(nSteps, 1, Tau_History);
	BzzVector tgrid    = GetBzzVector(nSteps, 1, T_History);

	dummy = "BATCH";
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

void OpenSMOKE_Batch::SaveOnBinaryFile(const std::string filename)
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

void OpenSMOKE_Batch::UpdateReactionRates(const double TT, const double PP, BzzVector &xx, BzzVector &rr)
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

void MyOdeSystem_Batch::ObjectBzzPrint(void)
{
    ::BzzPrint("\n OpenSMOKE_Batch Class\n");
}

void MyOdeSystem_Batch::GetSystemFunctions(BzzVector &x, double eta, BzzVector &f)
{
    if (iEnergy == false)
        ptBatch->ODESystem_Isothermal_Batch(x, eta, f);

    else if (iEnergy == true && iConstantPressure == false)
        ptBatch->ODESystem_NonIsothermal_ConstantVolume_Batch(x, eta, f);

    else if (iEnergy == true && iConstantPressure == true)
        ptBatch->ODESystem_NonIsothermal_ConstantPressure_Batch(x, eta, f);
}

void MyOdeSystem_Batch::assignBatch(OpenSMOKE_Batch *batch, const bool _iEnergy, const bool _iConstantPressure)
{
    ptBatch				= batch;
    iEnergy				= _iEnergy;
    iConstantPressure	= _iConstantPressure;
}