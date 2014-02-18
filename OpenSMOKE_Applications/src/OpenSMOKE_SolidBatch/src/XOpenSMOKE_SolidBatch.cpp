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
#include "XOpenSMOKE_SolidBatch.h"
#include "solid/OpenSMOKE_SolidMixture.h"
#include "basic/OpenSMOKE_Conversions.h"

OpenSMOKE_SolidBatch	*ptSolidBatch;

void ODEPrintExternal_SolidBatch(BzzVectorDouble &y, double csi)
{
    ptSolidBatch->ODEPrint(y, csi);
}

void OpenSMOKE_SolidBatch::UpdateProperties(int indexJacobian, int indexT) {};
void OpenSMOKE_SolidBatch::UpdateProperties_isothermal(int indexMemo) {};

OpenSMOKE_SolidBatch::OpenSMOKE_SolidBatch()
{
    ptSolidBatch	=  this;
	class_name		= "OpenSMOKE_SolidBatch";
	out_name		= "SolidBatch.out";
	
	Setup();
	
	// Compulsory
	assignedVolume				= false;
	assignedSolidKineticScheme	= false;			// solid kinetic scheme
	reactor_type				= NOSOLIDREACTOR;	// TG, CONSTP, CONSTV
	particle_model				= NOSOLIDMODEL;		// NORESISTANCE, SHRINKING, MODEL1D
	
	// Optional
	iVerboseSolid				= true;				// Solid output
}

void OpenSMOKE_SolidBatch::AssignSolidKineticScheme(OpenSMOKE_SolidMixture &_solid)
{
    solid     = &_solid;

    NC_solid = solid->NC_solid;

	ChangeDimensions(NC_solid, &solid_species_mass);
	ChangeDimensions(NC_solid, &dsolid_species_mass);
	ChangeDimensions(NC_solid, &solid_omega);
	ChangeDimensions(NC_solid, &solid_Rsolid);
	ChangeDimensions(NC,       &solid_Rgas);
	ChangeDimensions(NC_solid, &solid_C);

	ChangeDimensions(NC, &gas_out_h);		// TODO
	ChangeDimensions(NC, &gas_out_mass);	// TODO
	ChangeDimensions(NC, &Dmix_gas);		// TODO
	ChangeDimensions(NC, &h);				// TODO

	ChangeDimensions(NC, &omega_surface);	// TODO
	ChangeDimensions(NC, &C_surface);		// TODO
	ChangeDimensions(NC, &x_surface);		// TODO
	ChangeDimensions(NC, &domega_surface);	// TODO
	ChangeDimensions(NC, &km);				// TODO
	ChangeDimensions(NC, &gas_species_mass);// TODO

    // Control Variables
    assignedSolidKineticScheme = true;
}

void OpenSMOKE_SolidBatch::AssignConstantPressure()
{
	reactor_type = CONSTP;
}

void OpenSMOKE_SolidBatch::AssignConstantVolume()
{
	reactor_type = CONSTV;
}

void OpenSMOKE_SolidBatch::AssignThermogravimetry()
{
	reactor_type = TG_NONISO;
}

void OpenSMOKE_SolidBatch::AssignThermogravimetryIsothermal()
{
	if (reactor_type == TG_NONISO)	reactor_type = TG_ISO;
	else ErrorMessage("The #Isothermal option can be applied only to thermogravimetric analysis");
}

void OpenSMOKE_SolidBatch::AssignSolidModel_NoResistance()
{
	particle_model = NORESISTANCE;
}

void OpenSMOKE_SolidBatch::AssignSolidModel_ShrinkingCore()
{
	particle_model = SHRINKING;
}

void OpenSMOKE_SolidBatch::AssignSolidModel_Model1D()
{
	particle_model = MODEL1D;
}

void OpenSMOKE_SolidBatch::AssignSolidTemperature(const string units, const double value)
{
	solid_temperature_initial = OpenSMOKE_Conversions::conversion_temperature(value, units);
}

void OpenSMOKE_SolidBatch::AssignSolidDiameter(const string units, const double value)
{
	solid_diameter_initial = OpenSMOKE_Conversions::conversion_length(value, units);
}

void OpenSMOKE_SolidBatch::AssignSolidMassFractions(const vector<string> _names, const vector<double> _values)
{
	for(int k=0; k<_names.size(); k++)
	{
		int index = solid->recognize_species_from_biomass(_names[k]);
		solid_omega[index] = _values[k];
	}
}

void OpenSMOKE_SolidBatch::AssignSolidTotalMass(const string units, const double value)
{
	solid_total_mass_initial = OpenSMOKE_Conversions::conversion_mass(value, units);
}

void OpenSMOKE_SolidBatch::AssignVolume(const string units, const double value)
{
	volume			= OpenSMOKE_Conversions::conversion_volume(value, units);
    assignedVolume	= true;
}

void OpenSMOKE_SolidBatch::AssignEnd(const string units, const double value)
{
	TauTotal	= OpenSMOKE_Conversions::conversion_time(value, units);
    assignedEnd = true;
}

void OpenSMOKE_SolidBatch::SetConstantExchangeArea(const double value, const string units)
{
	iUserDefinedExchangeArea = CONSTANT;
	A = OpenSMOKE_Conversions::conversion_area(value, units);
}

void OpenSMOKE_SolidBatch::SetUserDefinedExchangeArea(const string fileName)
{
    iUserDefinedExchangeArea = USERDEFINED;
    ud_A_profile.AssignFromFile(fileName, "AREA");
	ud_A_profile.SetName(name_object + " - Exchange Area Profile");
}

void OpenSMOKE_SolidBatch::UnsetUserDefinedExchangeArea()
{
    iUserDefinedExchangeArea = NONE;
}

void OpenSMOKE_SolidBatch::Lock()
{
    if (assignedKineticScheme == false)
        ErrorMessage("The kinetic scheme was not defined!");
    if (assignedGlobalKineticScheme == false)
        ErrorMessage("The global kinetic scheme was not defined!");
    if (assignedEnd == false)
		ErrorMessage("The reactor contact time was not defined!");
    if (assignedInletFlows == false)
        ErrorMessage("The inlet flow was not defined!");
    if (assignedEnergy == false)
        ErrorMessage("The kind of reactor was not defined!");
	if (assignedSoot2EModel == false)
        ErrorMessage("The 2E Model was not defined!");
	if (assignedSolidKineticScheme == false)
        ErrorMessage("The solid kinetic scheme was not defined!");

	Initialize();
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
//																								   //
//									GAS MIXTURE PROPERTIES										   //
//																								   //
/////////////////////////////////////////////////////////////////////////////////////////////////////

void OpenSMOKE_SolidBatch::UpdateProperties_Thermogravimetry_NoResistance()
{
	// Bulk phase
	{
		// Gas mole fractions
		mix->GetMWAndMoleFractionsFromMassFractions(MWtot, x, omega);

		// Calcolo dei calori specifici [J/kgK]
		mix->SpeciesCp(T);
		Cp_gas = mix->MixCp_FromMassFractions(omega);

		// Enthalpy of gases
		mix->GetStandardEnthalpy_Mass(h, T);
	
		// *** Calcolo della conducibilitita termica [W/mK]
		mix->SpeciesConductivityFromFitting(T);										// TODO
		lambda_gas = mix->MixConductivity_FromMolarFractions(x);					// TODO

		// *** Calcolo dei coefficienti di diffusione [m2/s]
        mix->SpeciesDiffusivityFromFitting(T, P/1.e5);								// TODO
	    mix->MixDiffusionCoefficient_FromMolarFractionsAndPM(Dmix_gas, x);			// TODO
	}

	// Gas phase in the particle
	{
		// Solid total mass and solid mass fractions
		solid_mass	= solid_species_mass.GetSumElements();
		solid_omega	= solid_species_mass/solid_mass;

		// Solid volume and diameter
		solid_density	= solid->GetDensityApparent(solid_omega);
		solid_volume	= solid_mass / solid_density;
		solid_diameter	= pow(6./Constants::pi*solid_volume, 1./3.);
		solid_surface	= Constants::pi*solid_diameter*solid_diameter;
		solid_epsilon	= solid->GetEpsilon(solid_omega, solid_density);

		// Specific heat
		solid_Cp = solid->GetSpecificHeat(solid_omega, solid_temperature);

		// Solid Concentrations [kmol/m3]
		solid_C = solid_species_mass / solid_volume;
		ElementByElementProduct(solid_C, solid->uMW_solid, &solid_C);

		// Gasification and volatilization reactions
		solid->GetReactionRates(solid_C, x, solid_temperature);				// [kmol, m3, s]
		solid->GetFormationRates(solid_Rsolid, solid_Rgas);					// [kg/m3/s]
		solid->GetReactionEnthalpies(solid_temperature, gas_out_h);			// [kg/m3/s]

		// Total production
		solid_sumRgas	= solid_Rgas.GetSumElements()   * solid_volume;	// [kg/s]		
		solid_sumRsolid	= solid_Rsolid.GetSumElements() * solid_volume;	// [kg/s]
		gas_out_mass	= solid_volume*solid_Rgas;						// [kg/s] 

		// Thermal exchange
		solid_hexchange = 2.*lambda_gas/solid_diameter;					// [W/m2/K]
		solid_Q = solid_hexchange*solid_surface*(T-solid_temperature);	// [W]
	}
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
//																								   //
//									UPDATING FUNCTIONS											   //
//																								   //
/////////////////////////////////////////////////////////////////////////////////////////////////////

void OpenSMOKE_SolidBatch::UpdateTwoEquationModel(BzzVectorDouble &y, BzzVectorDouble &dy)
{
    soot2EModel->update(T, P/101325., rho, 0., x, y[indexTwoEquations], y[indexTwoEquations+1]);
    soot2EModel->formation_rates();

    for (int i=1;i<=NC;i++)
        domega[i] += soot2EModel->SGas[i] / rho;

    dy[indexTwoEquations]	= soot2EModel->s / rho;
    dy[indexTwoEquations+1] = soot2EModel->S / rho;
}

void OpenSMOKE_SolidBatch::UpdateHeatFlux(const double tau, const double csi)
{
	if (iUserDefinedHeatExchangeCoefficient == USERDEFINED)		U			= ud_U_profile.GiveMeValue(tau, MINUSONE);
	if (iUserDefinedAmbientTemperature		== USERDEFINED)		Tambient	= ud_Tambient_profile.GiveMeValue(tau, MINUSONE);
	if (iUserDefinedHeatFlux				== USERDEFINED)		Qe			= ud_Qe_profile.GiveMeValue(tau, MINUSONE);
	
	if (iUserDefinedHeatExchangeCoefficient != NONE)			Qe = U*(T-Tambient);
}

void OpenSMOKE_SolidBatch::UpdateExchangeArea(const double tau, const double csi)
{
	if (iUserDefinedExchangeArea == NONE)			A = pow(36.*Constants::pi*volume*volume, 1./3.);
	if (iUserDefinedExchangeArea == USERDEFINED)	A = ud_A_profile.GiveMeValue(tau, MINUSONE);
}


/////////////////////////////////////////////////////////////////////////////////////////////////////
//																								   //
//									SOLVING FUNCTIONS											   //
//																								   //
/////////////////////////////////////////////////////////////////////////////////////////////////////

void OpenSMOKE_SolidBatch::Initialize()
{
	// Gas phase
	{
		T            = inletStream->T;
		P            = inletStream->P;
		omega        = inletStream->omega;
		mass         = inletStream->rho*volume;

		omega_surface = omega;			// TODO
	}

	// Solid phase
	{
		solid_temperature	= solid_temperature_initial;
		solid_diameter		= solid_diameter_initial;
		solid_volume		= Constants::pi/6.*pow(solid_diameter, 3.);
		solid_density		= solid->GetDensityApparent(solid_omega);
		solid_mass			= solid_density*solid_volume;
		solid_mass_initial	= solid_mass;
		solid_species_mass	= solid_omega; solid_species_mass *= solid_mass;
		solid_ash			= solid_species_mass[solid->recognize_species("Ash")];

		double MWtot_surface;
		double rho_surface;
		mix->GetMWAndMoleFractionsFromMassFractions(MWtot_surface, x_surface, omega_surface);
		rho_surface = P/Constants::R_J_kmol/solid_temperature;
		gas_species_mass = (solid_epsilon*rho_surface*solid_volume)*omega_surface;
	}

    if (iUserDefinedTemperature == true)
        ud_temperature_profile.Check(0., T);

    if (iEnergy == false && (reactor_type==TG_ISO || reactor_type==TG_NONISO) && particle_model==NORESISTANCE)
        UpdateProperties_Thermogravimetry_NoResistance();
	else
		ErrorMessage("Properties not yet available!");

//    else if (iEnergy == true)				// Non-isothermal reactor
//      UpdateProperties(MINUSONE, NC+1);
}

void OpenSMOKE_SolidBatch::Solve()
{
    ptSolidBatch = this;

    double timeStart, timeEnd;
    BzzVectorDouble xMin;
    BzzVectorDouble xMax;
    BzzVectorDouble xInitial;

	PrepareFiles();


    ODE_SolidBatch_Object.assignSolidBatch(this, iEnergy, reactor_type, particle_model);

    // 1.A Isothermal Thermogravimetry
    if (iEnergy == false && reactor_type == TG_ISO && particle_model==NORESISTANCE)
    {
		int j,k;

        ChangeDimensions(NC+NC_solid, &xMin);  					
        ChangeDimensions(NC+NC_solid, &xMax);
        ChangeDimensions(NC+NC_solid, &xInitial);

		k=1;
		for(j=1;j<=NC_solid;j++)	xMin[k++] = ZERO;
		for(j=1;j<=NC;j++)			xMin[k++] = ZERO;


		k=1;
		for(j=1;j<=NC_solid;j++)	xMax[k++] = 1e16;
		for(j=1;j<=NC;j++)			xMax[k++] = ONE;

		k=1;
		for(j=1;j<=NC_solid;j++)	xInitial[k++] = solid_species_mass[j];
		for(j=1;j<=NC;j++)			xInitial[k++] = omega[j];
    }


    // 1.B NonIsothermal Thermogravimetry
/*    if (iEnergy == false && reactor_type == TG_NONISO && particle_model==NORESISTANCE)
    {
		int j,k;

        ChangeDimensions(NC+NC_solid+1, &xMin);  					
        ChangeDimensions(NC+NC_solid+1, &xMax);
        ChangeDimensions(NC+NC_solid+1, &xInitial);

		k=1;
		for(j=1;j<=NC_solid;j++)	xMin[k++] = ZERO;
									xMin[k++] = mix->TMIN;
		for(j=1;j<=NC;j++)			xMin[k++] = ZERO;


		k=1;
		for(j=1;j<=NC_solid;j++)	xMax[k++] = 1e16;
									xMax[k++] = mix->TMAX+5000;
		for(j=1;j<=NC;j++)			xMax[k++] = ONE;

		k=1;
		for(j=1;j<=NC_solid;j++)	xInitial[k++] = solid_species_mass[j];
									xInitial[k++] = solid_temperature;
		for(j=1;j<=NC;j++)			xInitial[k++] = omega[j];
    }
*/
	if (iEnergy == false && reactor_type == TG_NONISO && particle_model==NORESISTANCE)
    {
		int DIFFERENTIAL = 1;
		int ALGEBRAIC = 0;
		int j,k;

		DAE_SolidBatch_Object.assignSolidBatch(this, iEnergy, reactor_type, particle_model);

		
		BzzVectorInt DifAlg;
        ChangeDimensions(4, &xMin);  					
        ChangeDimensions(4, &xMax);
        ChangeDimensions(4, &xInitial);
        ChangeDimensions(4, &DifAlg);

		k=1;
		xMin[k++] = -1e16;
		xMin[k++] = -1e16;
		xMin[k++] = -1e16;
		xMin[k++] = mix->TMIN;


		k=1;
		xMax[k++] = 1e16;
		xMax[k++] = 1e16;
		xMax[k++] = 1e16;
		xMax[k++] = mix->TMAX;

		k=1;
		xInitial[k++] = solid_species_mass[solid->recognize_species("Char")];
		xInitial[k++] = ZERO;
		xInitial[k++] = omega_surface[mix->recognize_species("O2")];
		xInitial[k++] = solid_temperature;

		k=1;
		DifAlg[k++] = DIFFERENTIAL;
		DifAlg[k++] = ALGEBRAIC;
		DifAlg[k++] = ALGEBRAIC;
		DifAlg[k++] = ALGEBRAIC;


		dae(xInitial, 0., DifAlg, &DAE_SolidBatch_Object);

	    dae.StepPrint(ODEPrintExternal_SolidBatch);
		dae.SetMinimumConstraints(xMin);
		dae.SetMaximumConstraints(xMax);

		if (iRelativeTolerance == true)	dae.SetTollRel(MachEps()*relativeTolerance);
		if (iAbsoluteTolerance == true)	dae.SetTollAbs(absoluteTolerance);

		{
			countGlobalIterations = -1;  // From -1 to avoid to store results from the first iteration

			timeStart = BzzGetCpuTime();
		
			dae(TauTotal,TauTotal);
		
			timeEnd = BzzGetCpuTime();
		}

		if (iVerbose == true)
		{
			cout << endl;
			cout << "Number of function for Jacobian: " << dae.GetNumFunctionForJacobian()		<< endl;
			cout << "Numerical Jacobians: "				<< dae.GetNumNumericalJacobian()		<< endl;
			cout << "Time DAE solution: "				<< timeEnd - timeStart	<< " s"			<< endl << endl;
		}
	}

/*
    o(xInitial, 0., &ODE_SolidBatch_Object);

    o.StepPrint(ODEPrintExternal_SolidBatch);
    o.SetMinimumConstraints(xMin);
    o.SetMaximumConstraints(xMax);

    if (iRelativeTolerance == true)	o.SetTollRel(MachEps()*relativeTolerance);
    if (iAbsoluteTolerance == true)	o.SetTollAbs(absoluteTolerance);

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
*/
	CloseFiles();
}

void OpenSMOKE_SolidBatch::DefineFromFile(const string inputFile)
{
    double			double_value;
    string			string_value;
    int				int_value;
	vector<string>  string_vector;
	vector<double>  double_vector;

    OpenSMOKE_Dictionary_SolidBatch dictionary;
    dictionary.ParseFile(inputFile);


    // 1. COMPULSORY: Energy Equation solution
	AssignEnergy(dictionary.Return("#Energy"));


	// 2. COMPULSORY: Kind of reactor
	if (dictionary.Return("#ConstantPressure"))	AssignConstantPressure();
	if (dictionary.Return("#ConstantVolume"))	AssignConstantVolume();
	if (dictionary.Return("#Thermogravimetry"))	AssignThermogravimetry();
	if (dictionary.Return("#Isothermal"))		AssignThermogravimetryIsothermal();


	// 3. COMPULSORY: Particle model
    if (dictionary.Return("#NoResistance"))		AssignSolidModel_NoResistance();
    if (dictionary.Return("#ShrinkingCore"))	AssignSolidModel_ShrinkingCore();
    if (dictionary.Return("#Model1D"))			AssignSolidModel_Model1D();


	// 4. COMPULSORY: Particle properties
    if (dictionary.Return("#SolidTemperature", double_value, string_value))		AssignSolidTemperature(string_value, double_value);
    if (dictionary.Return("#SolidDiameter", double_value, string_value))		AssignSolidDiameter(string_value, double_value);
    if (dictionary.Return("#SolidMassFractions", double_vector, string_vector))	AssignSolidMassFractions(string_vector, double_vector);


	// SEMI-COMPULSORY: Solid Total Mass
    if (dictionary.Return("#SolidTotalMass", double_value, string_value))	AssignSolidTotalMass(string_value, double_value);
    if (dictionary.Return("#Volume", double_value, string_value))			AssignVolume(string_value, double_value);


    // COMPULSORY: Reactor Contact Time
	if (dictionary.Return("#ContactTime", double_value, string_value))	AssignEnd(string_value, double_value);


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
    if (dictionary.Return("#RelTol", double_value))
        SetRelativeTolerance(double_value);

    // OPTIONAL: Absolute tolerance
    if (dictionary.Return("#AbsTol", double_value))
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
    if (dictionary.Return("#A", double_value, string_value))
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
    if (dictionary.Return("#UserDefined_A", string_value))
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
    if (dictionary.Return("#ROPA", string_vector))
		SetROPAOnFile(string_vector);

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

void OpenSMOKE_SolidBatch::LabelODEFile()
{
    int i;

    if (iVerbose == true)
    {

        {
			fOutput << "Tau[s](1)    "		<< "\t";
            fOutput << "dummy(2)     "		<< "\t";
			fOutput << "T[K](3)      "		<< "\t";
            fOutput << "V[l](4)      "		<< "\t";
			fOutput << "P[atm](5)    "		<< "\t";
            fOutput << "dummy(6)     "		<< "\t";
            fOutput << "QR[W/m3](7)  "		<< "\t";
            fOutput << "QE[W/m3](8)  "		<< "\t";
            fOutput << "rho[kg/m3](9)"	    << "\t";
            fOutput << "MW(10)       "	    << "\t";
            fOutput << "dummy(11)    "	    << "\t";
            fOutput << "dummy(12)    "	    << "\t";
			if(key_species_names.size() == 1)
				fOutput << "Eta_" << key_species_names[0] << "(13)" << "\t";
			else
	            fOutput << "dummy(13)    "	    << "\t";

            fOutputCount = 14;
            for (i=1;i<=NC;i++)
                fOutput << mix->names[i] << " (x" << fOutputCount++ << ")  \t";
            for (i=1;i<=NC;i++)
                fOutput << mix->names[i] << " (w" << fOutputCount++ << ")  \t";

            fOutput << endl;
            fOutput << endl;
        }

        if (iTwoEquationModel == true)
            soot2EModel->write_label_on_file(fSoot2E, 1, "Tau[s]");

        if (mix->pah_manager.iPAH == 1)
		{
			fPAH << "Tau[s](1)      "		<< "\t";
            fPAH << "dummy(2)       "		<< "\t";
            fPAH << "T[K](3)        "		<< "\t";
			mix->pah_manager.GnuPlotInterface(fPAH, 3);
			fPAH << endl << endl;
		}

        if (mix->soot_manager.iSoot == 1)
        {
			fSoot << "Tau[s](1)      "		<< "\t";
            fSoot << "dummy(2)       "		<< "\t";
            fSoot << "T[K](3)        "		<< "\t";
            mix->soot_manager.GnuPlotInterface(fSoot, 4);
            fSoot << endl << endl;
        }

		if (iKeySpecies == true && key_species_names.size()>1)
        {
			fConversions << "Tau[s](1)      "		<< "\t";
            fConversions << "dummy(2)       "		<< "\t";
            fConversions << "T[K](3)        "		<< "\t";
			
			fOutputCount = 4;
            for(i=0;i<key_species_names.size();i++)
	            fConversions << "Eta_" << key_species_names[i] << "(" << fOutputCount++ << ")  \t";
			fConversions << endl << endl;
        }

		if (iVerboseEnergy == true)
		{
			fEnergy << "Tau[s](1)      "	<< "\t";
            fEnergy << "dummy(2)       "	<< "\t";
            fEnergy << "U[J/s](3)      "	<< "\t"; 
            fEnergy << "H[J/kg](4)     "	<< "\t";
			fEnergy << "Ek[J/kg](5)    "	<< "\t";
            fEnergy << "Ep[J/kg](6)    "	<< "\t";
            fEnergy << "QR[W/m3](7)    "	<< "\t";
            fEnergy << "QE[W/m3](8)    "	<< "\t";
            fEnergy << "Cp[J/kg/K](9)  "	<< "\t";
			fEnergy << endl;
			fEnergy << endl;
		}

		if (iVerboseReactionRates == true)
		{
			fReactionRates << "Tau[s](1)      "		<< "\t";
            fReactionRates << "dummy(2)       "		<< "\t";
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
            fFormationRates << "dummy(2)       "		<< "\t";
            fFormationRates << "T[K](3)        "		<< "\t";
			
			int countFileOutput = 4;
            for (i=1;i<=index_formation_rates.Size();i++)
                fFormationRates << names_formation_rates[i-1] << "(" << countFileOutput++ << ")  \t";
			fFormationRates << endl;
			fFormationRates << endl;
		}

		if (iVerboseROPA == true)
		{
			fROPA << "Tau[s](1)      "		<< "\t";
            fROPA << "dummy(2)       "		<< "\t";
            fROPA << "T[K](3)        "		<< "\t";
			
			int countFileOutput = 4;
			ROPA.PrintRateOfProductionAnalyses_Label(fROPA, countFileOutput);
		}

		if (iVerboseSensitivity == true)
		{
		}

		// Solid file
		{
			fSolid << "Tau[s](1)    "  << "\t";
            fSolid << "dummy(2)     "  << "\t";
            fSolid << "dummy(3)     "  << "\t";
            fSolid << "Tsolid[K](4) "  << "\t";
            fSolid << "Tgas[K](5)   "  << "\t";
            fSolid << "Dp[cm](6)    "  << "\t";
			fSolid << "TG[%](7)     "  << "\t";
			fSolid << "mSol[g](8)   "  << "\t";
			fSolid << "VSol[cm3](9) "  << "\t";
			fSolid << "epsilon(10)  "  << "\t";
			fSolid << "rho(11)      "  << "\t";
			fSolid << "Eta[%](12)   "  << "\t";
			fSolid << "EtaFree[%](13)" << "\t";
			fSolid << "dm_sol[kg/s](14)"	<< "\t";
			fSolid << "dm_gas[kg/s](15)"	<< "\t";
			fSolid << "h[W/m2/K](16)"		<< "\t";
			fSolid << "Q[W](17)     "		<< "\t";
			fSolid << "QR[W/m3](18) "		<< "\t";
			
			int count = 19;
			for(i=1;i<=NC_solid;i++)
				fSolid << solid->solid_names[i] << "(S-" << count++ << ")  " << "\t";
			for(i=1;i<=NC;i++)
				fSolid << mix->names[i] << "(G-" << count++ << ")  " << "\t";
			for(i=1;i<=NC_solid;i++)
				fSolid << solid->solid_names[i] << "(HS-" << count++ << ")  " << "\t";
			for(i=1;i<=NC;i++)
				fSolid << mix->names[i] << "(HG-" << count++ << ")  " << "\t";
			for(i=1;i<=solid->NR;i++)
				fSolid << "DH_" << i << "(" << count++ << ")  " << "\t";
			for(i=1;i<=solid->NR;i++)
				fSolid << "rSol_" << i << "(" << count++ << ")  " << "\t";
			for(i=1;i<=NC_solid;i++)
				fSolid << "CSol_" << solid->solid_names[i] << "(" << count++ << ")  " << "\t";
			
			fSolid << endl << endl;
		}
    }
}

void OpenSMOKE_SolidBatch::VideoGeneralInfo()
{
    cout.setf(ios::scientific);
    cout << endl;

    cout << "xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx" << endl;
    cout << "xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx" << endl;
    cout << endl;

    cout << "---------------------------------------------------------------------" << endl;
    cout << "Solid Batch Summary" << endl;
    cout << "---------------------------------------------------------------------" << endl;
    cout << "Contact time:\t\t"			<< TauTotal							<< " [s]"		<< endl;
	cout << "Pressure:\t\t"				<< P								<< " [Pa]"		<< endl;
    cout << "Density:\t\t"				<< rho								<< " [kg/m3]"	<< endl;
    cout << "Molecular weight:\t"		<< MWtot							<< " [kg/kmol]" << endl;
    cout << "Volume:\t\t\t"				<< volume							<< " [m3]"		<< endl;
    cout << "Mass:\t\t\t"				<< mass								<< " [kg]"		<< endl;
    cout << "Moles:\t\t\t"				<< moles							<< " [kmol]"	<< endl;
	cout << "Temperature:\t\t"			<< T								<< " [K]"		<< endl;
    cout <<  endl;
}

void OpenSMOKE_SolidBatch::VideoSummary()
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

void OpenSMOKE_SolidBatch::VideoFinalResult()
{
    int i;

    // Vectors
    BzzVectorDouble conversion(NC);

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

void OpenSMOKE_SolidBatch::SummaryOnFile(const string file_name)
{
	int i;

	ofstream fSummary;
	openOutputFileAndControl(fSummary, file_name);
    fSummary.setf(ios::scientific);

	double initial_volume	= mass / inletStream->rho;	// [m3]
	double initial_moles	= mass / inletStream->MW ;	// [kmol]


    fSummary << "               \tInitial      \tEnd"								<< endl;
	fSummary << "----------------------------------------------------------------------------"	<< endl;
    fSummary << "Time[s]        \t"	<< ZERO					<< "\t" << TauTotal			<< endl;
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

void OpenSMOKE_SolidBatch::ODEPrint(BzzVectorDouble &y, double t)
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
                    << moles			<< "\t"
                    << solid_mass		<< "\t"
                    << solid_diameter	<< "\t";

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
                    << dummy				<< "\t"
                    << T					<< "\t"
                    << volume/1000.			<< "\t"
                    << P/101325.			<< "\t"
                    << dummy				<< "\t"
                    << QReaction			<< "\t"
                    << Qe*A/volume     		<< "\t"
                    << rho					<< "\t"
                    << MWtot				<< "\t"
                    << dummy				<< "\t"
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
                soot2EModel->write_on_file(fSoot2E, t, y[indexTwoEquations], y[indexTwoEquations+1]);

            if (mix->soot_manager.iSoot == 1)
            {
                // Soot Distribution data
                mix->soot_manager.large_calculateAll(omega, x, c, rho);
                mix->soot_manager.small_calculateAll(omega, x, c, rho);

                fSoot << t		<< "\t";
                fSoot << dummy	<< "\t";
                fSoot << T		<< "\t";
                mix->soot_manager.print_main_data_on_file(fSoot);
                fSoot << endl;
            }

            if (mix->pah_manager.iPAH == 1)
            {
                // PAH Distribution data
                mix->pah_manager.pah_calculateAll(omega, x, c, rho);

                fPAH << t		<< "\t";
                fPAH << dummy	<< "\t";
                fPAH << T		<< "\t";
                mix->pah_manager.print_main_data_on_file(fPAH);
                fPAH << endl;
            }

			if (iKeySpecies == true && key_species_names.size()>1)
            {
                fConversions << t		<< "\t";
                fConversions << dummy	<< "\t";
                fConversions << T		<< "\t";

				for(i=0;i<key_species_names.size();i++)
	                fConversions << 1.-omega[key_species_index[i]]/inletStream->omega[key_species_index[i]]	<< "\t";
	        
                fConversions << endl;
            }

			if (iVerboseEnergy == true)
			{
				BzzVectorDouble h_mass(NC);
				mix->GetMixAveragedEnthalpy_Mass(h_mass, T);
	
				double H_mass  = Dot(omega, h_mass);			// [J/kg]
				double Ep_mass = P/rho;							// [J/kg]
				double Ek_mass = 0.;							// [J/kg]

				fEnergy << t									<< "\t"
						<< dummy								<< "\t"
						<< H_mass+Ep_mass						<< "\t"
						<< H_mass								<< "\t"
						<< Ek_mass								<< "\t"
						<< Ep_mass								<< "\t"
						<< QReaction							<< "\t"
						<< Qe*A		       						<< "\t"
						<< Cp		       						<< "\t"
						<< endl;
			}

			if (iVerboseReactionRates == true)
			{
				fReactionRates	<< t		<< "\t"
								<< dummy	<< "\t"
								<< T		<< "\t";
				
				for (i=1;i<=index_reaction_rates.Size();i++)
					fReactionRates	<< mix->r[index_reaction_rates[i]]	<< "\t";
				
				fReactionRates	<< endl;
			}

			if (iVerboseFormationRates == true)
			{
				fFormationRates	<< t		<< "\t"
								<< dummy	<< "\t"
								<< T		<< "\t";
				
				for (i=1;i<=index_formation_rates.Size();i++)
					fFormationRates	<< R[index_formation_rates[i]]	<< "\t";
				
				fFormationRates	<< endl;
			}

			if (iVerboseROPA == true)
			{
				ROPA.Run(mix->r);

				fROPA	<< t		<< "\t"
						<< dummy	<< "\t"
						<< T		<< "\t";

				ROPA.PrintRateOfProductionAnalyses(fROPA);
			}

			if (iVerboseSensitivity == true)
			{
				if ( (t-timeOld)>-1 && t>0. )
				{
					Sensitivity.Jacobian = 0.;
					o.GetLastJacobian(&Sensitivity.Jacobian);
					Sensitivity.BuildJAlfaMatrix(mix->r, rho);
					Sensitivity.GiveMe_SensitivityCoefficients(t-timeOld);
					Sensitivity.Normalize_SensitivityCoefficients(omega,x,MWtot);
				}

				fSensitivity << t << "\t" << Sensitivity.Sx[3][1] << "\t" << Sensitivity.Sx[3][3] << endl;
			}

			// Solid file
			{
				double conversion		= (solid_mass_initial-solid_mass)/solid_mass_initial;
				double conversion_free	= (solid_mass_initial-solid_mass)/(solid_mass_initial-solid_ash);

				fSolid << t										<< "\t";
				fSolid << dummy									<< "\t";
				fSolid << dummy									<< "\t";
				fSolid << T										<< "\t";
				fSolid << solid_temperature						<< "\t";
				fSolid << solid_diameter						<< "\t";
				fSolid << solid_mass/solid_mass_initial*100.	<< "\t";
				fSolid << solid_mass							<< "\t";
				fSolid << solid_volume							<< "\t";
				fSolid << solid_epsilon							<< "\t";
				fSolid << solid_density							<< "\t";
				fSolid << conversion*100.						<< "\t";
				fSolid << conversion_free*100.					<< "\t";
				fSolid << -solid_sumRsolid						<< "\t";
				fSolid << solid_sumRgas							<< "\t";
				fSolid << solid_hexchange						<< "\t";
				fSolid << solid_Q								<< "\t";
				fSolid << solid_QReaction						<< "\t";
				
				for(i=1;i<=NC_solid;i++)
					fSolid << solid_species_mass[i]		<< "\t";
				
				for(i=1;i<=NC;i++)
					fSolid << gas_out_mass[i]			<< "\t";

				for(i=1;i<=NC_solid;i++)
					fSolid << solid->h_solid[i]			<< "\t";
				
				for(i=1;i<=NC;i++)
					fSolid << gas_out_h[i]*mix->M[i]	<< "\t";

				for(i=1;i<=solid->NR;i++)
					fSolid << solid->DH[i]				<< "\t";

				for(i=1;i<=solid->NR;i++)
					fSolid << solid->r[i]				<< "\t";

				for(i=1;i<=NC_solid;i++)
					fSolid << solid_C[i]				<< "\t";

				fSolid << endl;
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

void OpenSMOKE_SolidBatch::EnergyAnalysis(OpenSMOKE_GasStream &outletStream)
{
	double H_mass_initial  = inletStream->massSpecificEnthalpy;			// [J/kg]
	double P_mass_initial  = inletStream->massSpecificPressureEnergy;	// [J/kg]
	double U_mass_initial  = inletStream->massSpecificInternalEnergy;	// [J/kg]
	
	double H_mass_final  = outletStream.massSpecificEnthalpy;			// [J/kg]
	double P_mass_final  = outletStream.massSpecificPressureEnergy;		// [J/kg]
	double U_mass_final  = outletStream.massSpecificInternalEnergy;		// [J/kg]

	double Work;				// Work done on the system (IN)
	if (reactor_type = CONSTP)	Work = -P*(volume-mass/inletStream->rho);
	else						Work = 0.;

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

void OpenSMOKE_SolidBatch::ODESystem_ThermogravimetryIsothermal_NoResistance(BzzVectorDouble &y, double t, BzzVectorDouble &dy)
{
	int k, j;

	// Recover variables
	k=1;
	for(j=1;j<=NC_solid;j++)	solid_species_mass[j]   = y[k++];
	for(j=1;j<=NC;j++)			omega[j]				= y[k++];


    // User Defined Temperature profile
    if (iUserDefinedTemperature == true)
        T = ud_temperature_profile.GiveMeValue(t, MINUSONE);
	solid_temperature = T;

    // Properties + Reaction rates (constant temperature)
    UpdateProperties_Thermogravimetry_NoResistance();


    // 1. Conservation equations for solid species
    dsolid_species_mass = solid_volume*solid_Rsolid;

    // 2. Conservation equations for gas species
    domega = ZERO;


    // Recovering residuals
    k=1;
	for(j=1;j<=NC_solid;j++)	dy[k++] = dsolid_species_mass[j];
	for(j=1;j<=NC;j++)			dy[k++] = domega[j];
}

void OpenSMOKE_SolidBatch::ODESystem_ThermogravimetryNonIsothermal_NoResistance(BzzVectorDouble &y, double t, BzzVectorDouble &dy)
{
	int k, j;

	// Recover variables
	k=1;
	for(j=1;j<=NC_solid;j++)	solid_species_mass[j]   = y[k++];
								solid_temperature		= y[k++];
	for(j=1;j<=NC;j++)			omega[j]				= y[k++];


    // User Defined Temperature profile
    if (iUserDefinedTemperature == true)
        T = ud_temperature_profile.GiveMeValue(t, MINUSONE);

    // Properties + Reaction rates (constant temperature)
    UpdateProperties_Thermogravimetry_NoResistance();


    // 1. Conservation equations for solid species
    dsolid_species_mass = solid_volume*solid_Rsolid;

	// 2. Energy equation for the solid phase
	dsolid_temperature  = -solid_volume * ( Dot(solid->h_solid, solid_Rsolid)+Dot(gas_out_h, solid_Rgas) ) + solid_Q;
	dsolid_temperature /= (solid_mass*solid_Cp);


    // 2. Conservation equations for gas species
    domega = ZERO;


    // Recovering residuals
    k=1;
	for(j=1;j<=NC_solid;j++)	dy[k++] = dsolid_species_mass[j];
								dy[k++]	= dsolid_temperature;
	for(j=1;j<=NC;j++)			dy[k++] = domega[j];
}

void OpenSMOKE_SolidBatch::DAESystem_ThermogravimetryNonIsothermal_NoResistance(BzzVectorDouble &y, double t, BzzVectorDouble &dy)
{
	int k;
	double eq_1, eq_2, eq_3;
	double mc;
	double nuI = 1.;
	int    iO2_gas		= mix->recognize_species("O2");		
	int    iChar_solid	= solid->recognize_species("Char");	
	double MC			= solid->MW_solid[iChar_solid];
	double MO2			= mix->M[iO2_gas];
	double DHC			= 3.2765e7;	//[J/kg]

	// Recover variables
	k=1;
	solid_species_mass[iChar_solid]	= y[k++];
	mc								= y[k++];
	omega_surface[iO2_gas]			= y[k++];
	solid_temperature				= y[k++];


    // User Defined Temperature profile
    if (iUserDefinedTemperature == true)
        T = ud_temperature_profile.GiveMeValue(t, MINUSONE);

    // Properties + Reaction rates (constant temperature)
    UpdateProperties_Thermogravimetry_NoResistance();

	// Kinetic constant
	double Ctot_bulk	= P/Constants::R_J_kmol/T;						// [kmol/m3]
	double Ctot_surface = P/Constants::R_J_kmol/T;
	double surface		= (4.*Constants::pi)*BzzPow2(solid_diameter*0.50);
	double surface_red	= surface/(solid_diameter*0.50);
	//double Kc			= solid->kR[2]/MC * solid_C[iChar_solid]/Ctot_bulk * solid_diameter/6.;
	double Kc			= 1.*exp(-20000./solid_temperature);
	double Kkin			= surface*Kc * MC *MWtot/ MO2 * Ctot_surface;

	// 
	dsolid_species_mass[iChar_solid]	= -mc;

	// Char mass flow rate
	eq_1		=	mc -	surface_red * rho*Dmix_gas[iO2_gas] * 
							log((1.+omega[iO2_gas]/nuI)/(1.+omega_surface[iO2_gas]/nuI));
	
	// Omega surface
	eq_2		=	mc - Kkin*omega_surface[iO2_gas];

	// Energy
	eq_3		=	mc*DHC - mc*Cp_gas*(solid_temperature-T) * 
					exp(-mc*Cp_gas/surface_red/lambda_gas) / (1.-exp(-mc*Cp_gas/surface_red/lambda_gas)); 
	
	//eq_3		=	solid_temperature - T;

	// Energy equation for the solid phase
	//solid_QReaction = -Dot(solid->r, solid->DH);	// [W/m3]
	//dsolid_temperature = ( solid_QReaction*solid_volume + solid_Q ) / (solid_mass*solid_Cp) ;


    // Recovering residuals
    k=1;
	dy[k++] = dsolid_species_mass[iChar_solid];
	dy[k++] = eq_1;
	dy[k++] = eq_2;
	dy[k++] = eq_3;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
//																								   //
//											DICTIONARY											   //
//																								   //
/////////////////////////////////////////////////////////////////////////////////////////////////////

OpenSMOKE_Dictionary_SolidBatch::OpenSMOKE_Dictionary_SolidBatch()
{
    SetupBase();

	// 1. Compulsory: if energy balances are solved
    {
		Add("#Energy",		'O', 'N', "The energy equation is solved");
		Add("#NoEnergy",	'O', 'N', "The energy equation is not solved");
	
		Compulsory("#Energy",	"#NoEnergy");
	    Conflict("#Energy",		"#NoEnergy");
	}

	// 2. Compulsory: kind of simulation
    {
		Add("#ConstantVolume",		'O', 'N', "The volume is kept constant");
		Add("#ConstantPressure",	'O', 'N', "The pressure is kept constant");
		Add("#Thermogravimetry",	'O', 'N', "Thermogravimetry");
		Add("#Isothermal",			'O', 'N', "Thermogravimetry");


   		Compulsory("#ConstantVolume", "#ConstantPressure", "#Thermogravimetry");
		Conflict("#ConstantVolume",		"#ConstantPressure");
		Conflict("#ConstantVolume",		"#Thermogravimetry");
		Conflict("#ConstantPressure",	"#Thermogravimetry");
		Conflict("#ConstantVolume",		"#Isothermal");
		Conflict("#ConstantPressure",	"#Isothermal");
	}

	// 3. Compulsory: kind of mass particle
    {
		Add("#NoResistance",		'O', 'N', "Solid particle: no diffusive resistance");
		Add("#ShrinkingCore",		'O', 'N', "Solid particle: shrinking core");
		Add("#Model1D",				'O', 'N', "Solid particle: 1D model");
	
   		Compulsory("#NoResistance", "#ShrinkingCore", "#Model1D");
		Conflict("#NoResistance",	"#ShrinkingCore");
		Conflict("#NoResistance",	"#Model1D");
		Conflict("#Model1D",		"#ShrinkingCore");
	}

	// 4. Compulsory: particle data
	{
		Add("#SolidTemperature",	'C', 'M', "Solid particle initial temperature");
		Add("#SolidDiameter",		'C', 'M', "Solid particle initial diameter");
		Add("#SolidMassFractions",	'C', 'L', "Solid mass fractions");
	}

	// 5. Semi-Compulsory (only if we are not permorfing a thermogravimetry)
    {
		Add("#Volume",				'O', 'M', "Reactor volume");
		Add("#SolidTotalMass",		'O', 'M', "Solid total mass (initial)");

		Compulsory("#Volume",			"#Thermogravimetry");
		Compulsory("#SolidTotalMass",	"#Thermogravimetry");
		Conflict("#Volume",				"#Thermogravimetry");
		Conflict("#SolidTotalMass",		"#Thermogravimetry");
	}

	// 6. Compulsory
    Add("#ContactTime",			'C', 'M', "Reactor contact time");


	// Optional
    Add("#Key",					'O', 'V', "Key species");

	// Optional
	Add("#ReactionRates",		'O', 'V', "Reaction rates");
	Add("#FormationRates",		'O', 'V', "Formation rates");
	Add("#ROPA",				'O', 'V', "Rate of Production Analysis");
	Add("#Sensitivity",			'O', 'V', "Sensitivity Analysis");

     
	// Optional
	Add("#VerboseEnergy",			'O', 'N', "Report on energy");
	Add("#Qe",						'O', 'M', "Heat flux");
	Add("#A",						'O', 'M', "Exchange area");
	Add("#U",						'O', 'M', "Heat thermal exchange coefficient");
	Add("#Tambient",				'O', 'M', "Ambient temperature");
	Add("#Viscosity",				'O', 'M', "Inlet viscosity");
	Add("#UserDefined_T",			'O', 'S', "The temperature is assigned from external file");
	Add("#UserDefined_Qe",			'O', 'S', "The heat flux is assigned from external file");
    Add("#UserDefined_A",			'O', 'S', "The exchange area is assigned from external file");
	Add("#UserDefined_U",			'O', 'S', "The heat exchange coefficient is assigned from external file");
    Add("#UserDefined_Tambient",	'O', 'S', "The ambient temperature is assigned from external file");
	Add("#2E_Model",				'O', 'N', "The Soot 2E Model is solved");

    Conflict("#Energy",					"#UserDefined_T");
	Conflict("#NoEnergy",				"#Qe");
	Conflict("#NoEnergy",				"#A");
	Conflict("#UserDefined_Qe",			"#Qe");
	Conflict("#UserDefined_A",			"#A");
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

    Lock();
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
//																								   //
//									ODE SYSTEM - CLASS DEFINITION								   //
//																								   //
/////////////////////////////////////////////////////////////////////////////////////////////////////

void MyOdeSystem_SolidBatch::ObjectBzzPrint(void)
{
    ::BzzPrint("\n OpenSMOKE_SolidBatch Class\n");
}

void MyOdeSystem_SolidBatch::GetSystemFunctions(BzzVectorDouble &x, double eta, BzzVectorDouble &f)
{
    if (iEnergy == false && reactor_type == TG_ISO		&&	particle_model == NORESISTANCE)
        ptSolidBatch->ODESystem_ThermogravimetryIsothermal_NoResistance(x, eta, f);
    else if (iEnergy == false && reactor_type == TG_NONISO	&&	particle_model == NORESISTANCE)
        ptSolidBatch->ODESystem_ThermogravimetryNonIsothermal_NoResistance(x, eta, f);
	else
	{
		cout << "MyOdeSystem_SolidBatch: Error AA01" << endl;
	}
}

void MyOdeSystem_SolidBatch::assignSolidBatch(OpenSMOKE_SolidBatch *solid_batch, const bool _iEnergy, 
											  const SolidReactorKind _reactor_type, const SolidParticleModel _particle_model)
{
    ptSolidBatch		= solid_batch;
    iEnergy				= _iEnergy;
    reactor_type		= _reactor_type;
    particle_model		= _particle_model;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
//																								   //
//									DAE SYSTEM - CLASS DEFINITION								   //
//																								   //
/////////////////////////////////////////////////////////////////////////////////////////////////////

void MyDaeSystem_SolidBatch::ObjectBzzPrint(void)
{

}

void MyDaeSystem_SolidBatch::GetSystemFunctions(BzzVectorDouble &x, double eta, BzzVectorDouble &f)
{
	ptSolidBatch->DAESystem_ThermogravimetryNonIsothermal_NoResistance(x, eta, f);
}

void MyDaeSystem_SolidBatch::assignSolidBatch(OpenSMOKE_SolidBatch *solid_batch, const bool _iEnergy, 
											  const SolidReactorKind _reactor_type, const SolidParticleModel _particle_model)
{
    ptSolidBatch		= solid_batch;
    iEnergy				= _iEnergy;
    reactor_type		= _reactor_type;
    particle_model		= _particle_model;
}