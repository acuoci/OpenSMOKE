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
#include "idealreactors/flame1d/OpenSMOKE_Flame1D_DataManager.h" 
#include "addons/OpenSMOKE_ElementFluxAnalysis.h"
#include "liquid/OpenSMOKE_LiquidSpecies.h"
#include "liquid/OpenSMOKE_LiquidProperties_Database.h"

void OpenSMOKE_Flame1D_DataManager::ErrorMessage(const std::string message)
{
    cout << endl;
    cout << "Class:  OpenSMOKE_Flame1D_DataManager"	<< endl;
    cout << "Object: " << name_object			<< endl;
    cout << "Error:  " << message				<< endl;
    cout << "Press a key to continue... "		<< endl;
    getchar();
    exit(-1);
}

void OpenSMOKE_Flame1D_DataManager::WarningMessage(const std::string message)
{
    cout << endl;
    cout << "Class:  OpenSMOKE_Flame1D_DataManager"	<< endl;
    cout << "Object: "		<< name_object		<< endl;
    cout << "Warning:  "	<< message			<< endl;
    cout << "Press a key to continue... "		<< endl;
    getchar();
}

OpenSMOKE_Flame1D_DataManager::OpenSMOKE_Flame1D_DataManager()
{
	name_object			= "[Name not assigned]";

	iGlobalKinetics						= false;
	iUnsteady							= false;
	iBackUp								= false;
	iEquivalenceRatioForPremixedFlames	= false;
	iAssignedAdaptiveGridCoefficients	= false;
	iUserDefinedFolderName				= false;

	SetDefaultValues();

	kind_of_subphysics					= FLAME1D_SUBPHYSICS_NONE;
	kind_of_flame						= FLAME1D_PHYSICS_OPPOSED;
	kind_of_adaptive_grid				= ADAPTIVE_CHEMKIN_POW;
	sensitivityLennardJonesMode			= OPENSMOKE_FITTING_ALL;

	// Grid refinement
	nDiff = 4;
	nGrad = 2;	
	deltaDiff = 0.05;
	deltaGrad = 0.05;

	// Derivatives
	iDerG = 'U';
	iDerT = 'U';
	iDerW = 'U';

	// Radial gradients
	radialGradientC = 0.;
	radialGradientO = 0.;

	// Output
	nStepsVideo  = 10;		
	nStepsFile	 = 1000;
	nStepsBackUp = 1000;

	// Radiation
	iGasRadiation = false;
	environmentTemperature = 298.15;	
	
	// Tolerances
	rel_nls_Tolerances = sqrt(MachEpsFloat());				
	abs_nls_Tolerances = 1.e-10;	

	rel_dae_Tolerances = 100.*MachEpsFloat();				
	abs_dae_Tolerances = 1.e-10;	

	// Premixed flames
	alfa				= 1.100;
	iAreaProfile		= 0;
	MassFlowRate		= 0.;
	iFixedTemperature	= 0;

	// Corrections
	correctionReactionRates = 1.;
	
	iAssignedFlameThickness				= false;
	iAssignedFlamePosition				= false;
	iAssignedGrid						= false;
	iAssignedPeaks						= false;
	iAssignedInitialTemperatureProfile	= false;
	iAssignedFixedTemperatureProfile	= false;
	iAssignedVideoSteps					= false;
	iAssignedFileSteps					= false;
	iAssignedBackupSteps				= false;
	iAssignedNlsRelativeTolerance		= false;
	iAssignedNlsAbsoluteTolerance		= false;
	iAssignedDaeRelativeTolerance		= false;
	iAssignedDaeAbsoluteTolerance		= false;
	iAssignedFuelRadialGradient			= false;
	iAssignedOxidizerRadialGradient		= false;
	iAssignedEnvironmentTemperature		= false;
	iAssignedGasRadiation				= false;
	iAssignedSootRadiation				= false;
	iAssignedGridRefineGradient			= false;
	iAssignedGridRefineCurvature		= false;
	iAssignedGridRefineGradientStep		= false;
	iAssignedGridRefineCurvatureStep	= false;
	iAssignedDerivativeG				= false;
	iAssignedDerivativeT				= false;
	iAssignedDerivativeW				= false;
	iAssignedSensitivity				= false;
	iAssignedFormationRates				= false;
	iAssignedROPA						= false;
	iVerboseAssignedROPA				= false;
	iAssignedReactionRates				= false;
	iAssignedExperiment					= false;
	iAssignedGlobalElementFluxAnalysis  = false;
	iFocusReactionRates					= false;
	iFocusElemetFluxAnalysis			= false;

	iAssignedOutletTemperature			= false;
	iAssignedOutletMassFractions		= false;
	iAssignedOutletMoleFractions		= false;
	iAssignedFlameSpeedIndex			= false;
	iAssignedFlameSpeedTemperature		= false;
	iAssignedStretchingFactor			= false;
	iAssignedCrossSectionProfile		= false;
	iAssignedCrossSectionMIT			= false;
	iUserDefinedLewisNumbers            = false;

	initial_time_step = -1;
	max_integration_order = -1;

	// TODO
	iAssignedFixedTemperatureProfileProvisional	= false;

	counterUnchanged = 0;
	TMaxOld = 300.;
	iRobustTemperature = false;

	poolfire_grid_alfa_fuel = 1.10; 
	poolfire_grid_alfa_oxidizer = 1.10; 
	poolfire_grid_point_fraction = 0.80; 
	poolfire_grid_distance_fraction = 0.60; 
	correctionFactorVaporizationHeat = 1.0;
	correctionFactorVaporPressure = 1.0;
	correctionFactorSpecificHeat = 1.0;
	correctionFactorThermalConductivity = 1.0;

	// Soot
	iRadiativeSootModel = RADIATIVE_SOOT_MODEL_NONE;
	iThermophoreticEffect = false;
	iDepositionWall = false;

	bin_index_zero = 0;
	bin_density_A = 0.;
	bin_index_final = 0;
	bin_density_B = 0.;

	iCorrectDiffusionFormulation = true;
	iPhysicalSootDiffusionCoefficients = 0;
}

void OpenSMOKE_Flame1D_DataManager::SetName(const std::string name)
{
	name_object = name;
}

void OpenSMOKE_Flame1D_DataManager::Assign(OpenSMOKE_ReactingGas *_mix)
{
	mix = _mix;
}

void OpenSMOKE_Flame1D_DataManager::Assign(OpenSMOKE_Flame1D *_flame)
{
	flame = _flame;
}

void OpenSMOKE_Flame1D_DataManager::SetDefaultValues()
{
	i2E									= false;
	iQMOM								= false;
	iSoretEffect						= false;
	iThermophoreticEffect				= false;
	iLewisNumbers						= false;
	iTurbulentDiffusivity				= false;
	iSingleContributions				= false;
	iCorrectionDiffusivity				= false;
	iCorrectionFormationEnthalpy		= false;
	iFlameSpeedAnalysis					= false;
	iOpposedFlameAnalysis				= false;
	iCorrectionReactionRates			= false;
	iUnityLewisNumbers					= false;
	iVerboseMixtureProperties			= false;
	iVerboseFluxes						= false;
	iFakeTemperatureThermalConductivity = false;
	iPhysicalSootDiffusionCoefficients  = 0;
	Df = 1.80;
}

void OpenSMOKE_Flame1D_DataManager::Setup(const std::string kind)
{
	if (kind == "PREMIXED")
	{
		SetDefaultValues();
		kind_of_subphysics			= FLAME1D_SUBPHYSICS_NONE;
		kind_of_flame				= FLAME1D_PHYSICS_PREMIXED;
	}
    else if (kind == "OPPOSED")
	{
		SetDefaultValues();
		kind_of_subphysics			= FLAME1D_SUBPHYSICS_NONE;
		kind_of_flame				= FLAME1D_PHYSICS_OPPOSED;
	}
    else if (kind == "TWIN")
	{
		SetDefaultValues();
		kind_of_subphysics			= FLAME1D_SUBPHYSICS_NONE;
		kind_of_flame				= FLAME1D_PHYSICS_TWIN;
	}
	else if (kind == "PREMIXED_SOOT")
	{
		SetDefaultValues();
		i2E							= true;
		kind_of_subphysics			= FLAME1D_SUBPHYSICS_SOOT;
		kind_of_flame				= FLAME1D_PHYSICS_PREMIXED;
	}
    else if (kind == "OPPOSED_SOOT")
	{
		SetDefaultValues();
		i2E							= true;
		kind_of_subphysics			= FLAME1D_SUBPHYSICS_SOOT;
		kind_of_flame				= FLAME1D_PHYSICS_OPPOSED;
	}
    else if (kind == "TWIN_SOOT")
	{
		SetDefaultValues();
		i2E							= true;
		kind_of_subphysics			= FLAME1D_SUBPHYSICS_SOOT;
		kind_of_flame				= FLAME1D_PHYSICS_TWIN;
	}
	else if (kind == "PREMIXED_QMOM")
	{
		SetDefaultValues();
		iQMOM						= true;
		kind_of_subphysics			= FLAME1D_SUBPHYSICS_QMOM;
		kind_of_flame				= FLAME1D_PHYSICS_PREMIXED;
	}
    else if (kind == "OPPOSED_QMOM")
	{
		SetDefaultValues();
		iQMOM						= true;
		kind_of_subphysics			= FLAME1D_SUBPHYSICS_QMOM;
		kind_of_flame				= FLAME1D_PHYSICS_OPPOSED;
	}
    else if (kind == "TWIN_QMOM")
	{
		SetDefaultValues();
		iQMOM						= true;
		kind_of_subphysics			= FLAME1D_SUBPHYSICS_QMOM;
		kind_of_flame				= FLAME1D_PHYSICS_TWIN;
	}
}

void OpenSMOKE_Flame1D_DataManager::readFromFileForOpposed(const std::string fileName)
{
	DefineFromFileOpposed(fileName);

	iteration = 1;
	iterationVideoCounter	= nStepsVideo;
	iterationFileCounter	= nStepsFile;
	iterationBackUpCounter	= nStepsBackUp;

	nCold = iX.Size();
	
	// ASSIGNED PROFILES
	if (iTemperatureProfile==1 || iTemperatureProfile==2) 
	{
		ud_temperature_profile.AssignFromFile(ud_temperature_profile_file_name, "TEMPERATURE");
		ud_temperature_profile.SetName(name_object + " - User Defined Temperature Profile");

		ud_temperature_profile.Check(0., TC);
		ud_temperature_profile.Check(L,  TO);
	}
}

void OpenSMOKE_Flame1D_DataManager::PasteInputFile(const std::string fileName)
{
	ofstream fOutput;
	openOutputFileAndControl(fOutput, fileName);
	fOutput.setf(ios::scientific);

	for(unsigned int i=1;i<=input_file_lines.size()-1;i++)
	{
		     if (i == indexLinePoints)				fOutput << "	#Points               " << flame->Np << endl;
		else if (i == indexLineFuelVelocity)		fOutput << "	#FuelVelocity         " << VC*100. << "   cm/s" << endl;
		else if (i == indexLineOxidizerVelocity)	fOutput << "	#OxidizerVelocity     " << VO*100. << "   cm/s" << endl;
		else if (i == indexLineFlameSpeedIndex)		fOutput << "	#FlameSpeedIndex      " << iFixedTemperature << endl;
		else if (i == indexLineEquivalenceRatio)	fOutput << "	#EquivalenceRatio     " << "TODO" << endl;

		else										fOutput << input_file_lines[i] << endl;
	}

	fOutput.close();
}

void OpenSMOKE_Flame1D_DataManager::printFileForOpposed(const std::string fileName, double vFuel, double vAir)
{
	PasteInputFile(fileName);
}

void OpenSMOKE_Flame1D_DataManager::printFileForPremixed(const std::string fileName, double MassFlowRate_Updated)
{
	PasteInputFile(fileName);
}

void OpenSMOKE_Flame1D_DataManager::readFromFileForPremixed(const std::string fileName)
{
	DefineFromFilePremixed(fileName);

	iteration = 1;
	iterationVideoCounter = nStepsVideo;
	iterationFileCounter = nStepsFile;
	iterationBackUpCounter = nStepsBackUp;

	if (iTemperatureProfile==1 || iTemperatureProfile==2 || iAssignedFixedTemperatureProfileProvisional==true) 
	{
		ud_temperature_profile.AssignFromFile(ud_temperature_profile_file_name, "TEMPERATURE");
		ud_temperature_profile.SetName(name_object + " - User Defined Temperature Profile");
		ud_temperature_profile.Check(0., TC);
	}

	if (iAreaProfile==1) 
	{
		ud_cross_section_profile.AssignFromFile(ud_cross_section_profile_file_name, "AREA");
		ud_cross_section_profile.SetName(name_object + " - User Defined Cross Section Profile");
		CrossSectionalArea = ud_cross_section_profile.GiveMeValue(0., 0.);
	}
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
//									DICTIONARY													   //
/////////////////////////////////////////////////////////////////////////////////////////////////////

OpenSMOKE_Dictionary_Flame1D::OpenSMOKE_Dictionary_Flame1D()
{
    SetupBase();
	SetName("OpenSMOKE_Flame1D_Dictionary");

	Add("#Fuel",          'C', 'S', "Fuel name");
	Add("#Oxidizer",      'C', 'S', "Oxidizer name");
	Add("#Inert",         'C', 'S', "Inert name");

	Add("#Pressure",		  'C', 'M', "Pressure");
	Add("#Distance",		  'C', 'M', "Distance between nozzles");

	Add("#Points",		      'C', 'I', "Number of grid points");
	Add("#OutputSpecies",     'C', 'V', "Output species");

   	Add("#FlamePosition",	  'O', 'M', "Flame position from fuel side");
	Add("#FlameThickness",	  'O', 'M', "Flame thickness");
	Add("#Grid",			  'O', 'S', "Grid kind: equispaced || centered || user");
	Add("#PeakMassFractions", 'O', 'L', "Peak species (mass fractions)");

	Add("#InitialTemperatureProfile", 'O', 'S', "Initial temperature profile from file");
	Add("#FixedTemperatureProfile",   'O', 'S', "Fixed temperature profile from file");
	Add("#FixedTemperatureProfileProvisional",   'O', 'S', "Fixed temperature profile from file");
 
	Add("#NlsRelativeTolerance",	'O', 'D', "Nls relative tolerance (default sqrt(MachEpsFloat()))");
	Add("#NlsAbsoluteTolerance",	'O', 'D', "Nls absolute tolerance (default 1.e-10)");
	Add("#DaeRelativeTolerance",	'O', 'D', "Dae relative tolerance (default 1e2*MachEpsFloat())");
	Add("#DaeAbsoluteTolerance",	'O', 'D', "Dae absolute tolerance (default 1.e-10)");

	Add("#InitialTimeStep",			'O', 'M', "Initial time step");
	Add("#MaximumIntegrationOrder",	'O', 'I', "Maximum integration order");

	Add("#SoretEffect",				'O', 'N', "Soret effect");
	Add("#ThermophoreticEffect",	'O', 'N', "Thermophoretic effect");
	Add("#VerboseMixtureProperties",'O', 'N', "Verbose mixture properties");

	Add("#PhysicalSootDiffusionCoefficients", 'O', 'I', "Physical soot diffusion coefficients (BIN where to cut)");

	Add("#UncorrectDiffusionFormulation", 'O', 'N', "Uncorrect diffusion formulation: omega*V = -D*grad(omega)");

	Add("#LewisNumbers",			'O', 'N', "Lewis numbers on file");
	Add("#UnityLewisNumbers",		'O', 'N', "Lewis numbers are assumed equal to 1 for every species");
	Add("#UserDefinedLewisNumbers",	'O', 'L', "User defined Lewis numbers");
	Add("#CorrectionReactionRates", 'O', 'D', "Correction of reaction rates");
	Add("#TurbulentDiffusivity",	'O', 'D', "Turbulent Diffusivity Enhancing Factor (with respect to the inert specie)");

	Add("#GasRadiation",			'O', 'N', "Radiative heat transfer from gas phase");
	Add("#SootRadiation",			'O', 'S', "Radiative heat transfer from soot: none, fluent, bilger, widmann");
	Add("#EnvironmentTemperature",	'O', 'M', "Environment temperature");

	Add("#ReactionRates",			'O', 'V', "Reaction rates on file (list of reaction indices)");
	Add("#FormationRates",			'O', 'V', "Formation rates on file (list of species names)");
	Add("#ROPA",					'O', 'N', "Rate of Production Analysis");
	Add("#VerboseROPA",				'O', 'V', "Rate of Production Analysis (list of species names)");
	Add("#Sensitivity",				'O', 'V', "Sensitivity Analysis (list of species names)");
	Add("#Experiment",				'O', 'V', "Experiment Analysis (list of species names)");
	Add("#ElementFluxAnalysis",		'O', 'V', "Element Flux Analysis (list of elements, lower case)");
	Add("#GlobalElementFluxAnalysis", 'O', 'V', "Global Element Flux Analysis (tag xA xB units)");

	Add("#GridRefineGradient",      'O', 'I', "Maximum number of new grid points in grid refinement (gradient)");
	Add("#GridRefineCurvature",     'O', 'I', "Maximum number of new grid points in grid refinement (curvature)");

	Add("#GridRefineGradientStep",  'O', 'D', "Grid refinement control parameter (gradient)");
	Add("#GridRefineCurvatureStep", 'O', 'D', "Grid refinement control parameter (curvature)");

	Add("#DerivativeT",             'O', 'H', "Temperature derivative: U || C || B || F");
	Add("#DerivativeW",             'O', 'H', "Mass fractions derivative: U || C || B || F");

	Add("#SingleContributions",		'O', 'V', "Single contributions");
	
	Add("#AdaptiveGrid",				'O', 'S', "Adaptive Grid kind: CHEMKIN_POW, CHEMKIN, CHEMKIN_SQRT, GRADIENT, EISEMAN");
	Add("#AdaptiveGridCoefficients",	'O', 'V', "Adaptive Grid Coefficients: Alfa, Beta, Gamma");

	Add("#ChangeFrequencyFactor",		'O', 'V', "Change frequency factor (reaction, variation)");
	Add("#ChangeActivationEnergy",		'O', 'V', "Change activation energy (int reaction, double new value, std::string units)");
	Add("#ChangeDiffusivity",			'O', 'V', "Change diffusion coefficient (species, variation)");
	Add("#LennardJonesMode",			'O', 'S', "Lennard-Jones Mode: ALL, ONLY_DIFFUSIVITIES, ONLY_CONDUCTIVITY, ONLY_VISCOSITY");
	Add("#ChangeFormationEnthalpy",		'O', 'V', "Change formation enthalpy (species, variation, units)");

	Add("#FakeTemperatureIncrementThermalConductivity",  'O', 'M', "Fake temperature increment for thermal conductivity");

	Add("#DepositionWall",			'O', 'N', "Deposition wall");

	Add("#BINDensities",			'O', 'V', "BIN densities in kg/m3 (index0 rho0 indexF rhoF) [default 10 1500 20 1700]");

	Add("#VerboseFluxes",			'O', 'N', "Verbose fluxes (Fick + Soret)");
}

void OpenSMOKE_Dictionary_Flame1D::PrepareForOpposedFlames()
{
	Add("#Geometry",				'C', 'S', "Burner geometry: axis || planar");
	
	Add("#FuelVelocity",			'O', 'M', "Fuel stream velocity");
	Add("#OxidizerVelocity",		'C', 'M', "Oxidizer stream velocity");

	Add("#FuelTemperature",			'O', 'M', "Fuel stream temperature");
	Add("#OxidizerTemperature",		'C', 'M', "Oxidizer stream temperature");

	Add("#FuelMassFractions",		'O', 'L', "Fuel stream mass fractions");
    Add("#FuelMoleFractions",		'O', 'L', "Fuel stream mole fractions");
	Add("#OxidizerMassFractions",   'O', 'L', "Oxidizer stream mass fractions");
    Add("#OxidizerMoleFractions",	'O', 'L', "Oxidizer stream mole fractions");

	Add("#FuelRadialGradient",		'O', 'D', "Fuel stream radial gradient");
	Add("#OxidizerRadialGradient",	'O', 'D', "Oxidizer stream radial gradient");
	Add("#TemperaturePeak",			'C', 'M', "Peak temperature");

	Add("#DerivativeG",             'O', 'H', "Mass flow rate derivative: U || C || B || F");

	Add("#PoolFire",				'O', 'V', "Pool fire: database folder, ...");
	Add("#PoolFireGridOptions",		'O', 'V', "Pool fire grid options");

	Add("#StretchingFactor",		'O', 'D', "Grid stretching factor");

	Add("#RobustTemperature",		'O', 'N', "Robust temperature checks");
	Add("#PoolFireCorrectionVaporizationHeat",	  'O', 'D', "Correction factor for pool vaporization heat: H = alfa*H");
	Add("#PoolFireCorrectionVaporPressure",	      'O', 'D', "Correction factor for pool vapor pressure: Pv = alfa*Pv");
	Add("#PoolFireCorrectionSpecificHeat",	      'O', 'D', "Correction factor for pool specific heat: Cp = alfa*Cp");
	Add("#PoolFireCorrectionThermalConductivity", 'O', 'D', "Correction factor for pool thermal conductivity: k = alfa*k");

	Compulsory("#FuelMassFractions",     "#FuelMoleFractions",		"#PoolFire");
	Compulsory("#OxidizerMassFractions", "#OxidizerMoleFractions");

	Compulsory("#FuelTemperature",  "#PoolFire");
	Compulsory("#FuelVelocity",     "#PoolFire");

	Conflict("#FuelTemperature",    "#PoolFire");
	Conflict("#FuelVelocity",		"#PoolFire");
	Conflict("#FuelMassFractions",  "#PoolFire");
	Conflict("#FuelMoleFractions",  "#PoolFire");

    Conflict("#FuelMassFractions",		"#FuelMoleFractions");
    Conflict("#OxidizerMassFractions",	"#OxidizerMoleFractions");
	Conflict("#InitialTemperatureProfile", "#FixedTemperatureProfile");

	Conflict("#UnityLewisNumbers", "#UserDefinedLewisNumbers");

    Lock();
}

void OpenSMOKE_Dictionary_Flame1D::PrepareForPremixedFlames()
{
	Add("#MassFlowRate",			'O', 'M', "Mass flow rate");
	Add("#InletVelocity",			'O', 'M', "Inlet velocity");
	Add("#InletTemperature",		'C', 'M', "Inlet temperature");
	Add("#InletMassFractions",		'O', 'L', "Inlet mass fractions");
    Add("#InletMoleFractions",		'O', 'L', "Inlet mole fractions");

	Add("#EquivalenceRatio",		'O', 'D', "Equivalence ratio");
	Add("#FuelMassFractions",		'O', 'L', "Fuel stream mass fractions");
    Add("#FuelMoleFractions",		'O', 'L', "Fuel stream mole fractions");
	Add("#OxidizerMassFractions",   'O', 'L', "Oxidizer stream mass fractions");
    Add("#OxidizerMoleFractions",	'O', 'L', "Oxidizer stream mole fractions");

	Add("#CrossSection",			'O', 'M', "Cross section");
	Add("#CrossSectionProfile",		'O', 'S', "Cross section profile from file");
	Add("#CrossSectionMIT",			'O', 'N', "Cross section profile from MIT formula");

	Add("#OutletTemperature",		'C', 'M', "Outlet temperature");
    Add("#OutletMoleFractions",		'O', 'L', "Outlet mole fractions");
    Add("#OutletMassFractions",		'O', 'L', "Outlet mass fractions");

    Add("#StretchingFactor",		'O', 'D', "Grid stretching factor");

    Add("#FlameSpeedIndex",			'O', 'I', "Index point for fixed temperature");
    Add("#FlameSpeedTemperature",	'O', 'M', "Fixed temperature");

	Compulsory("#MassFlowRate",			"#InletVelocity");
	Compulsory("#InletMassFractions",	"#InletMoleFractions", "#EquivalenceRatio");
	Compulsory("#OutletMassFractions",	"#OutletMoleFractions");
	Compulsory("#CrossSection",			"#CrossSectionProfile");

    Conflict("#MassFlowRate",			"#InletVelocity");
    Conflict("#InletMassFractions",		"#InletMoleFractions");
    Conflict("#OutletMoleFractions",	"#OutletMassFractions");
	Conflict("#CrossSection",			"#CrossSectionProfile");
	Conflict("#CrossSectionProfile",	"#CrossSectionMIT");
	Conflict("#InitialTemperatureProfile", "#FixedTemperatureProfile");

	Conflict("#EquivalenceRatio", "#InletMassFractions");
	Conflict("#EquivalenceRatio", "#InletMoleFractions");
	Conflict("#FuelMassFractions", "#InletMoleFractions");
	Conflict("#FuelMassFractions", "#InletMassFractions");
	Conflict("#OxidizerMassFractions", "#InletMoleFractions");
	Conflict("#OxidizerMassFractions", "#InletMassFractions");
	Conflict("#FuelMassFractions", "#FuelMoleFractions");
	Conflict("#OxidizerMoleFractions", "#OxidizerMassFractions");

    Lock();
}

void OpenSMOKE_Flame1D_DataManager::CheckDictionaryForOpposedFlames(OpenSMOKE_Dictionary_Flame1D &dictionary)
{
    char    char_value;
	double  double_value;
    std::string  string_value;
	vector<string> string_vector;
	vector<double> double_vector;

	if (dictionary.Return("#FuelVelocity", double_value, string_value))
        AssignFuelVelocity(string_value, double_value);

	if (dictionary.Return("#OxidizerVelocity", double_value, string_value))
        AssignOxidizerVelocity(string_value, double_value);

	if (dictionary.Return("#FuelTemperature", double_value, string_value))
        AssignFuelTemperature(string_value, double_value);

	if (dictionary.Return("#OxidizerTemperature", double_value, string_value))
        AssignOxidizerTemperature(string_value, double_value);

	if (dictionary.Return("#TemperaturePeak", double_value, string_value))
        AssignTemperaturePeak(string_value, double_value);

	if (dictionary.Return("#Geometry", string_value))
        AssignGeometry(string_value);

    if (dictionary.Return("#FuelMassFractions", double_vector, string_vector))
		AssignFuelMassFractions(string_vector, double_vector);
    if (dictionary.Return("#FuelMoleFractions", double_vector, string_vector))
		AssignFuelMoleFractions(string_vector, double_vector);

    if (dictionary.Return("#OxidizerMassFractions", double_vector, string_vector))
		AssignOxidizerMassFractions(string_vector, double_vector);
    if (dictionary.Return("#OxidizerMoleFractions", double_vector, string_vector))
		AssignOxidizerMoleFractions(string_vector, double_vector);

	if (dictionary.Return("#PoolFire", string_vector))
        SetPoolFire(string_vector);

	if (dictionary.Return("#PoolFireGridOptions", string_vector))
        SetPoolFireGridOptions(string_vector);
	
	if (dictionary.Return("#PoolFireCorrectionVaporizationHeat", double_value))
		SetPoolFireCorrectionFactorVaporizationHeat(double_value);
	if (dictionary.Return("#PoolFireCorrectionVaporPressure", double_value))
		SetPoolFireCorrectionFactorVaporPressure(double_value);
	if (dictionary.Return("#PoolFireCorrectionThermalConductivity", double_value))
		SetPoolFireCorrectionFactorThermalConductivity(double_value);
	if (dictionary.Return("#PoolFireCorrectionSpecificHeat", double_value))
		SetPoolFireCorrectionFactorSpecificHeat(double_value);
	

	if (dictionary.Return("#FuelRadialGradient", double_value, string_value))
		SetFuelRadialGradient(string_value, double_value);
	
	if (dictionary.Return("#OxidizerRadialGradient", double_value, string_value))
		SetOxidizerRadialGradient(string_value, double_value);

	if (dictionary.Return("#DerivativeG", char_value))
		SetDerivativeG(char_value);

	if (dictionary.Return("#RobustTemperature"))
		SetRobustTemperature();

	if (dictionary.Return("#StretchingFactor", double_value))
        SetStretchingFactor(double_value);
}

void OpenSMOKE_Flame1D_DataManager::CheckDictionaryForPremixedFlames(OpenSMOKE_Dictionary_Flame1D &dictionary)
{
    int     int_value;
	double  double_value;
    std::string  string_value;
	vector<string> string_vector;
	vector<double> double_vector;

	if (dictionary.Return("#MassFlowRate", double_value, string_value))
        AssignMassFlowRate(string_value, double_value);

	if (dictionary.Return("#InletVelocity", double_value, string_value))
        AssignInletVelocity(string_value, double_value);
	
	if (dictionary.Return("#InletTemperature", double_value, string_value))
        AssignInletTemperature(string_value, double_value);
	
	if (dictionary.Return("#CrossSection", double_value, string_value))
        AssignCrossSection(string_value, double_value);

    if (dictionary.Return("#InletMassFractions", double_vector, string_vector))
		AssignInletMassFractions(string_vector, double_vector);
    if (dictionary.Return("#InletMoleFractions", double_vector, string_vector))
		AssignInletMoleFractions(string_vector, double_vector);

    if (dictionary.Return("#EquivalenceRatio", double_value))
		AssignEquivalenceRatioForPremixedFlames(double_value);
    if (dictionary.Return("#FuelMassFractions", double_vector, string_vector))
		AssignFuelMassFractionsForPremixedFlames(string_vector, double_vector);
    if (dictionary.Return("#FuelMoleFractions", double_vector, string_vector))
		AssignFuelMoleFractionsForPremixedFlames(string_vector, double_vector);
    if (dictionary.Return("#OxidizerMassFractions", double_vector, string_vector))
		AssignOxidizerMassFractionsForPremixedFlames(string_vector, double_vector);
    if (dictionary.Return("#OxidizerMoleFractions", double_vector, string_vector))
		AssignOxidizerMoleFractionsForPremixedFlames(string_vector, double_vector);
	AssignInletStreamCompositionForPremixedFlames();


	if (dictionary.Return("#OutletTemperature", double_value, string_value))
        SetOutletTemperature(string_value, double_value);

    if (dictionary.Return("#OutletMassFractions", double_vector, string_vector))
		SetOutletMassFractions(string_vector, double_vector);
    if (dictionary.Return("#OutletMoleFractions", double_vector, string_vector))
		SetOutletMoleFractions(string_vector, double_vector);

	if (dictionary.Return("#FlameSpeedIndex", int_value))
        SetFlameSpeedIndex(int_value);
	
	if (dictionary.Return("#FlameSpeedTemperature", double_value, string_value))
        SetFlameSpeedTemperature(string_value, double_value);

	if (dictionary.Return("#CrossSectionProfile", string_value))
        SetCrossSectionProfile(string_value);

	if (dictionary.Return("#CrossSectionMIT"))
        SetCrossSectionMIT();

	if (dictionary.Return("#StretchingFactor", double_value))
        SetStretchingFactor(double_value);
}

void OpenSMOKE_Flame1D_DataManager::CheckDictionary(OpenSMOKE_Dictionary_Flame1D &dictionary)
{
    int     int_value;
    char    char_value;
	double  double_value;
    std::string  string_value;
	vector<string> string_vector;
	vector<double> double_vector;

	if (dictionary.Return("#Pressure", double_value, string_value))
        AssignPressure(string_value, double_value);

	if (dictionary.Return("#Distance", double_value, string_value))
        AssignDistance(string_value, double_value);

	if (dictionary.Return("#Fuel", string_value))
        AssignFuel(string_value);

	if (dictionary.Return("#Oxidizer", string_value))
        AssignOxidizer(string_value);

	if (dictionary.Return("#Inert", string_value))
        AssignInert(string_value);

	if (dictionary.Return("#OutputSpecies", string_vector))
		AssignOutputSpecies(string_vector);

	if (dictionary.Return("#Points", int_value))
        AssignGridPoints(int_value);
	
	if (dictionary.Return("#Grid", string_value))
        SetGrid(string_value);

	if (dictionary.Return("#FlamePosition", double_value, string_value))
        SetFlamePosition(string_value, double_value);

	if (dictionary.Return("#FlameThickness", double_value, string_value))
        SetFlameThickness(string_value, double_value);

    if (dictionary.Return("#PeakMassFractions", double_vector, string_vector))
		SetPeaks(string_vector, double_vector);

	if (dictionary.Return("#InitialTemperatureProfile", string_value))
        SetInitialTemperatureProfile(string_value);

	if (dictionary.Return("#FixedTemperatureProfile", string_value))
        SetFixedTemperatureProfile(string_value);

	if (dictionary.Return("#FixedTemperatureProfileProvisional", string_value))
        SetFixedTemperatureProfileProvisional(string_value);

	if (dictionary.Return("#Nvideo", int_value))
		SetVideoSteps(int_value);

	if (dictionary.Return("#Nfile", int_value))
		SetFileSteps(int_value);

	if (dictionary.Return("#Nbackup", int_value))
		SetBackupSteps(int_value);

	if (dictionary.Return("#NlsRelativeTolerance", double_value))
		SetNlsRelativeTolerance(double_value);

	if (dictionary.Return("#NlsAbsoluteTolerance", double_value))
		SetNlsAbsoluteTolerance(double_value);

	if (dictionary.Return("#DaeRelativeTolerance", double_value))
		SetDaeRelativeTolerance(double_value);

	if (dictionary.Return("#DaeAbsoluteTolerance", double_value))
		SetDaeAbsoluteTolerance(double_value);

	if (dictionary.Return("#InitialTimeStep", double_value, string_value))
		SetInitialTimeStep(string_value, double_value);

	if (dictionary.Return("#MaximumIntegrationOrder", int_value))
		SetMaximumIntegrationOrder(int_value);

	if (dictionary.Return("#EnvironmentTemperature", double_value, string_value))
		SetEnvironmentTemperature(string_value, double_value);

	if (dictionary.Return("#GasRadiation"))	
		SetGasRadiation();

	if (dictionary.Return("#DepositionWall"))	
		SetDepositionWall();

	if (dictionary.Return("#SootRadiation", string_value))	
		SetSootRadiation(string_value);

	if (dictionary.Return("#SoretEffect"))	
		SetSoretEffect();

	if (dictionary.Return("#PhysicalSootDiffusionCoefficients", int_value))
		SetPhysicalSootDiffusionCoefficients(int_value);

	if (dictionary.Return("#ThermophoreticEffect"))	
		SetThermophoreticEffect();

	if (dictionary.Return("#UncorrectDiffusionFormulation"))
		SetUncorrectDiffusionFormulation();

	if (dictionary.Return("#VerboseMixtureProperties"))	
		SetVerboseMixtureProperties();

	if (dictionary.Return("#VerboseFluxes"))
		SetVerboseFluxes();

	if (dictionary.Return("#LewisNumbers"))	
		SetLewisNumbers();

	if (dictionary.Return("#UnityLewisNumbers"))	
		SetUnityLewisNumbers();

	if (dictionary.Return("#UserDefinedLewisNumbers", double_vector, string_vector))	
		SetUserDefinedLewisNumbers(string_vector, double_vector);

	if (dictionary.Return("#TurbulentDiffusivity", double_value))
		SetTurbulentDiffusivity(double_value);

	if (dictionary.Return("#GridRefineGradient", int_value))
		SetGridRefineGradient(int_value);

	if (dictionary.Return("#GridRefineCurvature", int_value))
		SetGridRefineCurvature(int_value);

	if (dictionary.Return("#GridRefineGradientStep", double_value))
		SetGridRefineGradientStep(double_value);

	if (dictionary.Return("#GridRefineCurvatureStep", double_value))
		SetGridRefineCurvatureStep(double_value);

	if (dictionary.Return("#DerivativeT", char_value))
		SetDerivativeT(char_value);

	if (dictionary.Return("#DerivativeW", char_value))
		SetDerivativeW( char_value);

    if (dictionary.Return("#ReactionRates", string_vector))
		SetReactionRatesOnFile(string_vector);

    if (dictionary.Return("#BINDensities", string_vector))
		SetBINDensities(string_vector);

    if (dictionary.Return("#FormationRates", string_vector))
		SetFormationRatesOnFile(string_vector);

    if (dictionary.Return("#VerboseROPA", string_vector))
		SetVerboseROPAOnFile(string_vector);

	if (dictionary.Return("#ROPA"))
		SetROPAOnFile();

    if (dictionary.Return("#ElementFluxAnalysis", string_vector))
		SetElementFluxAnalysisOnFile(string_vector);

    if (dictionary.Return("#GlobalElementFluxAnalysis", string_vector))
		SetGlobalElementFluxAnalysisOnFile(string_vector);

    if (dictionary.Return("#Experiment", string_vector))
		SetExperimentOnFile(string_vector);

    if (dictionary.Return("#Sensitivity", string_vector))
		SetSensitivityOnFile(string_vector); 

    if (dictionary.Return("#SingleContributions", string_vector))
		SetSingleContributions(string_vector); 

	if (dictionary.Return("#AdaptiveGrid", string_value))
		SetAdaptiveGrid(string_value); 

    if (dictionary.Return("#AdaptiveGridCoefficients", string_vector))
		SetAdaptiveGridCoefficients(string_vector); 

	if (dictionary.Return("#CorrectionReactionRates", double_value))
		SetCorrectionReactionRates(double_value);

	if (dictionary.Return("#ChangeFrequencyFactor", string_vector))
		SetChangeFrequencyFactor(string_vector);

	if (dictionary.Return("#ChangeActivationEnergy", string_vector))
		SetChangeActivationEnergy(string_vector);

	if (dictionary.Return("#ChangeDiffusivity", string_vector))
		SetChangeDiffusivity(string_vector);

	if (dictionary.Return("#LennardJonesMode", string_value))
        SetLennardJonesMode(string_value);

	if (dictionary.Return("#ChangeFormationEnthalpy", string_vector))
		SetChangeFormationEnthalpy(string_vector);

	if (dictionary.Return("#FakeTemperatureIncrementThermalConductivity", double_value, string_value))
		SetChangeFakeTemperatureThermalConductivity(double_value, string_value);
}

void OpenSMOKE_Flame1D_DataManager::DefineFromFileOpposed(const std::string inputFile)
{
    OpenSMOKE_Dictionary_Flame1D dictionary;

	dictionary.PrepareForOpposedFlames();
    dictionary.ParseFile(inputFile);
	CheckDictionary(dictionary);
	CheckDictionaryForOpposedFlames(dictionary);
	LockForOpposedFlames();

	CopyInputFile(inputFile);
}

void OpenSMOKE_Flame1D_DataManager::DefineFromFilePremixed(const std::string inputFile)
{
    OpenSMOKE_Dictionary_Flame1D dictionary;

	dictionary.PrepareForPremixedFlames();
    dictionary.ParseFile(inputFile);
	CheckDictionary(dictionary);
	CheckDictionaryForPremixedFlames(dictionary);
	LockForPremixedFlames();

	CopyInputFile(inputFile);
}

void OpenSMOKE_Flame1D_DataManager::CopyInputFile(const std::string inputFile)
{
	const int SIZE = 400;
	char comment[SIZE];

	ifstream iFile;
	openInputFileAndControl(iFile, inputFile);
	
	input_file_lines.resize(0);
	input_file_lines.push_back("List of lines");

	while(!iFile.eof())
	{
		iFile.getline(comment, SIZE);
		input_file_lines.push_back(comment);
	}
	iFile.close();

	int number_of_lines = input_file_lines.size()-1;

	{
		bool indexLinePointsExist = false;
		for(int i=1;i<=number_of_lines;i++)
			if (StringFindSubString(input_file_lines[i], "#Points") == true)
			{
				indexLinePoints = i;
				indexLinePointsExist = true;
				break;
			}

		if (indexLinePointsExist == false)
			ErrorMessage("#Points line does not exist in input file");
	}

	{
		bool indexLineFuelVelocityExist = false;
		for(int i=1;i<=number_of_lines;i++)
			if (StringFindSubString(input_file_lines[i], "#FuelVelocity") == true)
			{
				indexLineFuelVelocity = i;
				indexLineFuelVelocityExist = true;
				break;
			}

		if (indexLineFuelVelocityExist == false)
			indexLineFuelVelocity = 0;
	}

	{
		bool indexLineOxidizerVelocityExist = false;
		for(int i=1;i<=number_of_lines;i++)
			if (StringFindSubString(input_file_lines[i], "#OxidizerVelocity") == true)
			{
				indexLineOxidizerVelocity = i;
				indexLineOxidizerVelocityExist = true;
				break;
			}

		if (indexLineOxidizerVelocityExist == false)
			indexLineOxidizerVelocity = 0;
	}

	{
		bool indexLineEquivalenceRatioExist = false;
		for(int i=1;i<=number_of_lines;i++)
			if (StringFindSubString(input_file_lines[i], "#EquivalenceRatio") == true)
			{
				indexLineEquivalenceRatio = i;
				indexLineEquivalenceRatioExist = true;
				break;
			}

		if (indexLineEquivalenceRatioExist == false)
			indexLineEquivalenceRatio = 0;
	}

	{
		bool indexLineFlameSpeedIndexExist = false;
		for(int i=1;i<=number_of_lines;i++)
			if (StringFindSubString(input_file_lines[i], "#FlameSpeedIndex") == true)
			{		
				indexLineFlameSpeedIndex = i;
				indexLineFlameSpeedIndexExist = true;
				break;
			}

		if (indexLineFlameSpeedIndexExist == false)
			indexLineFlameSpeedIndex = 0;
	}

}

void OpenSMOKE_Flame1D_DataManager::LockForOpposedFlames()
{
	if (iAssignedFlameThickness == false)	wmix = 0.50*L;
	if (iAssignedFlamePosition  == false)	xcen = 0.50*L;
	if (iAssignedGrid			== false)	gridKind = "EQUISPACED";
	if (iAssignedPeaks			== false)	
	{
		ChangeDimensions(0, &xPeaks); 
		namePeaks.resize(0);
	}

	if (iAssignedInitialTemperatureProfile == false && iAssignedFixedTemperatureProfile == false)	
			iTemperatureProfile = 0;

	ChangeDimensions(0, &iX);
	for(int i=1;i<=mix->NumberOfSpecies();i++)
		if (XC[i]!=0. || XO[i]!=0.)	iX.Append(i);
}

void OpenSMOKE_Flame1D_DataManager::LockForPremixedFlames()
{
	if (MassFlowRate == 0.)
	{
		BzzVector x(mix->NumberOfSpecies());
	
		for(int i=1;i<=xC.Size();i++)
			x[mix->recognize_species(nameC[i])] = xC[i];
		double MWmix = mix->GetMWFromMoleFractions(x);
		double rho = P_Pascal*MWmix/Constants::R_J_kmol/TC;
		MassFlowRate = VC*rho*CrossSectionalArea;
	}

	if (iAssignedFlameThickness == false)	wmix = 2./1000.;
	if (iAssignedFlamePosition  == false)	xcen = 3./1000.;
	if (iAssignedGrid			== false)	gridKind = "STRETCHED";
	if (iAssignedPeaks			== false)	
	{
		ChangeDimensions(0, &xPeaks); 
		namePeaks.resize(0);
	}
	if (iAssignedInitialTemperatureProfile == false &&
		iAssignedFixedTemperatureProfile == false)	
			iTemperatureProfile = 0;

	ChangeDimensions(0, &iX);
	for(int i=1;i<=mix->NumberOfSpecies();i++)
		if (XC[i]!=0. || XO[i]!=0.)	iX.Append(i);
}

void OpenSMOKE_Flame1D_DataManager::SetOutputFolder(const std::string _outputFolderName)
{
	iUserDefinedFolderName = true;
	userDefinedFolderName = _outputFolderName;
}
void OpenSMOKE_Flame1D_DataManager::AssignPressure(const std::string units, const double value)
{
    P_Pascal = OpenSMOKE_Conversions::conversion_pressure(value, units);
	P_atm	 = P_Pascal / OpenSMOKE_Conversions::Pa_from_atm;
	P_bar	 = P_Pascal / OpenSMOKE_Conversions::Pa_from_bar;
}

void OpenSMOKE_Flame1D_DataManager::AssignDistance(const std::string units, const double value)
{
    L = OpenSMOKE_Conversions::conversion_length(value, units);
}

void OpenSMOKE_Flame1D_DataManager::AssignFuelVelocity(const std::string units, const double value)
{
    VC = OpenSMOKE_Conversions::conversion_velocity(value, units);
}

void OpenSMOKE_Flame1D_DataManager::AssignOxidizerVelocity(const std::string units, const double value)
{
    VO = OpenSMOKE_Conversions::conversion_velocity(value, units);
}

void OpenSMOKE_Flame1D_DataManager::AssignFuelTemperature(const std::string units, const double value)
{
    TC = OpenSMOKE_Conversions::conversion_temperature(value, units);
}

void OpenSMOKE_Flame1D_DataManager::AssignOxidizerTemperature(const std::string units, const double value)
{
    TO = OpenSMOKE_Conversions::conversion_temperature(value, units);
}

void OpenSMOKE_Flame1D_DataManager::AssignTemperaturePeak(const std::string units, const double value)
{
    Tpeak = OpenSMOKE_Conversions::conversion_temperature(value, units);
}

void OpenSMOKE_Flame1D_DataManager::AssignGeometry(const std::string string_value)
{
	if (string_value != "AXIS" && string_value != "PLANAR")
		ErrorMessage("#Geometry options: AXIS || PLANAR");

	geometry = string_value;
}

void OpenSMOKE_Flame1D_DataManager::AssignFuel(const std::string string_value)
{
	nameFuel = string_value;
	jFUEL = mix->recognize_species(nameFuel);
}

void OpenSMOKE_Flame1D_DataManager::AssignOxidizer(const std::string string_value)
{
	nameOxidizer = string_value;
	jO2 = mix->recognize_species(nameOxidizer);
}

void OpenSMOKE_Flame1D_DataManager::AssignInert(const std::string string_value)
{
	nameInert = string_value;
	jINERT = mix->recognize_species(nameInert);
}

void OpenSMOKE_Flame1D_DataManager::AssignOutputSpecies(const vector<string> string_vector)
{
	if (string_vector.size()==1 && string_vector[0] == "ALL")
	{
		int iOutput = mix->NumberOfSpecies();
		ChangeDimensions(iOutput, &iOut);
		nameOutput.resize(iOutput+1); 
		for(int i=1;i<=iOutput;i++)
		{
			nameOutput[i]	= mix->names[i];
			iOut[i]			= i;
		}
	}
	else
	{
		int iOutput = string_vector.size();
		ChangeDimensions(iOutput, &iOut);
		nameOutput.resize(iOutput+1);
		for(int i=1;i<=iOutput;i++)
		{
			nameOutput[i] = string_vector[i-1];
			iOut[i] = mix->recognize_species(string_vector[i-1]);
		}
	}
}

void OpenSMOKE_Flame1D_DataManager::AssignGridPoints(const int int_value)
{
	Np = int_value;
}

void OpenSMOKE_Flame1D_DataManager::AssignFuelMassFractions(const vector<string> _names, const vector<double> _values)
{
	int i;
	double MWmix;
	BzzVector y(mix->NumberOfSpecies());
	ChangeDimensions(mix->NumberOfSpecies(), &XC);
	
	for(i=1;i<=int(_names.size());i++)
		y[mix->recognize_species(_names[i-1])] = _values[i-1];
	mix->GetMWAndMoleFractionsFromMassFractions(MWmix, XC, y);

	// Assign composition
	int iC = _names.size();
	ChangeDimensions(iC, &xC); 
	nameC.resize(iC+1);

	for(i=1;i<=int(_names.size());i++)
	{
		nameC[i] = _names[i-1];
		xC[i]	 = XC[mix->recognize_species(_names[i-1])];
	}	

	ChangeDimensions(0, &iXC);
	for(i=1;i<=mix->NumberOfSpecies();i++)
		if (XC[i]!=0.)	iXC.Append(i);
}

void OpenSMOKE_Flame1D_DataManager::AssignFuelMassFractions(BzzVector &omega_fuel)
{
	int i;
	double MWmix;
	ChangeDimensions(mix->NumberOfSpecies(), &XC);
	mix->GetMWAndMoleFractionsFromMassFractions(MWmix, XC, omega_fuel);

	// Assign composition
	int iC = 0;
	for(i=1;i<=mix->NumberOfSpecies();i++)
		if (omega_fuel[i] > 0.) iC++;
	
	ChangeDimensions(iC, &xC); 
	nameC.resize(iC+1);

	iC = 0;
	for(i=1;i<=mix->NumberOfSpecies();i++)
		if (omega_fuel[i] > 0.) 
		{
			iC++;
			nameC[iC] = mix->names[i];
			xC[iC]	 = XC[i];
		}	

	ChangeDimensions(0, &iXC);
	for(i=1;i<=mix->NumberOfSpecies();i++)
		if (XC[i]!=0.)	iXC.Append(i);
}

void OpenSMOKE_Flame1D_DataManager::AssignFuelMoleFractions(const vector<string> _names, const vector<double> _values)
{
	int i;
	int iC = _names.size();
	ChangeDimensions(iC, &xC);
	ChangeDimensions(mix->NumberOfSpecies(), &XC);
	nameC.resize(iC+1);

	for(i=1;i<=int(_names.size());i++)
	{
		nameC[i] = _names[i-1];
		xC[i]	 = _values[i-1];
		XC[mix->recognize_species(_names[i-1])] = _values[i-1];
	}	

	ChangeDimensions(0, &iXC);
	for(i=1;i<=mix->NumberOfSpecies();i++)
		if (XC[i]!=0.)	iXC.Append(i); 
}
	
void OpenSMOKE_Flame1D_DataManager::AssignOxidizerMassFractions(const vector<string> _names, const vector<double> _values)
{
	int i;
	double MWmix;
	BzzVector y(mix->NumberOfSpecies());
	ChangeDimensions(mix->NumberOfSpecies(), &XO);
	
	for(i=1;i<=int(_names.size());i++)
		y[mix->recognize_species(_names[i-1])] = _values[i-1];
	mix->GetMWAndMoleFractionsFromMassFractions(MWmix, XO, y);

	// Assign composition
	int iO = _names.size();
	ChangeDimensions(iO, &xO); 
	nameO.resize(iO+1);

	for(i=1;i<=int(_names.size());i++)
	{
		nameO[i] = _names[i-1];
		xO[i]	 = XO[mix->recognize_species(_names[i-1])];
	}

	ChangeDimensions(0, &iXO);
	for(i=1;i<=mix->NumberOfSpecies();i++)
		if (XO[i]!=0.)	iXO.Append(i); 	
}
    
void OpenSMOKE_Flame1D_DataManager::AssignOxidizerMoleFractions(const vector<string> _names, const vector<double> _values)
{
	int i;
	int iO = _names.size();
	ChangeDimensions(iO, &xO); 
	ChangeDimensions(mix->NumberOfSpecies(), &XO);
	nameO.resize(iO+1);

	for(i=1;i<=int(_names.size());i++)
	{
		nameO[i] = _names[i-1];
		xO[i]	 = _values[i-1];
		XO[mix->recognize_species(_names[i-1])] = _values[i-1];
	}

	ChangeDimensions(0, &iXO);
	for(i=1;i<=mix->NumberOfSpecies();i++)
		if (XO[i]!=0.)	iXO.Append(i); 	
}

void OpenSMOKE_Flame1D_DataManager::SetPeaks(const vector<string> _names, const vector<double> _values)
{
	int nPeaks = _names.size();	
	ChangeDimensions(nPeaks, &iPeaks); 
	ChangeDimensions(nPeaks, &xPeaks); 
	namePeaks.resize(nPeaks+1);
	
	for(int i=1;i<=nPeaks;i++)
	{
		namePeaks[i] = _names[i-1];
		xPeaks[i]	 = _values[i-1];
		iPeaks[i]	 = mix->recognize_species(_names[i-1]);
	}

	iAssignedPeaks = true;
}

void OpenSMOKE_Flame1D_DataManager::SetFixedTemperatureProfile(const std::string string_value)
{
	iTemperatureProfile = 1;
	ud_temperature_profile_file_name = string_value;
	iAssignedFixedTemperatureProfile = true;
}

void OpenSMOKE_Flame1D_DataManager::SetFixedTemperatureProfileProvisional(const std::string string_value)
{
	ud_temperature_profile_file_name = string_value;
	iAssignedFixedTemperatureProfileProvisional = true;
}

void OpenSMOKE_Flame1D_DataManager::SetInitialTemperatureProfile(const std::string string_value)
{
	iTemperatureProfile = 2;
	ud_temperature_profile_file_name = string_value;
	iAssignedInitialTemperatureProfile = true;
}

void OpenSMOKE_Flame1D_DataManager::SetFlameThickness(const std::string units, const double value)
{
    wmix = OpenSMOKE_Conversions::conversion_length(value, units);
	iAssignedFlameThickness = true;
}

void OpenSMOKE_Flame1D_DataManager::SetFlamePosition(const std::string units, const double value)
{
    xcen = OpenSMOKE_Conversions::conversion_length(value, units);
	iAssignedFlamePosition = true;
}

void OpenSMOKE_Flame1D_DataManager::SetGrid(const std::string string_value)
{
	if (string_value != "CENTERED" && string_value != "EQUISPACED" &&
		string_value != "STRETCHED" && string_value != "STRETCHED_POOL_FIRE" && string_value != "USER" &&
		string_value != "STRETCHED_STAGNATION")
		ErrorMessage("#Geometry options: EQUISPACED || CENTERED || STRETCHED || STRETCHED_POOL_FIRE || USER || STRETCHED_STAGNATION");

	gridKind = string_value;
	iAssignedGrid = true;
}

void OpenSMOKE_Flame1D_DataManager::SetVideoSteps(const int int_value)
{
	nStepsVideo = int_value;
	iAssignedVideoSteps = true;
}

void OpenSMOKE_Flame1D_DataManager::SetFileSteps(const int int_value)
{
	nStepsFile = int_value;
	iAssignedFileSteps = true;
}

void OpenSMOKE_Flame1D_DataManager::SetBackupSteps(const int int_value)
{
	nStepsBackUp = int_value;
	iAssignedBackupSteps = true;
}

void OpenSMOKE_Flame1D_DataManager::SetNlsRelativeTolerance(const double double_value)
{
	rel_nls_Tolerances = double_value;
	iAssignedNlsRelativeTolerance = true;
}

void OpenSMOKE_Flame1D_DataManager::SetNlsAbsoluteTolerance(const double double_value)
{
	abs_nls_Tolerances = double_value;
	iAssignedNlsAbsoluteTolerance = true;
}

void OpenSMOKE_Flame1D_DataManager::SetDaeRelativeTolerance(const double double_value)
{
	rel_dae_Tolerances = double_value;
	iAssignedDaeRelativeTolerance = true;
}

void OpenSMOKE_Flame1D_DataManager::SetDaeAbsoluteTolerance(const double double_value)
{
	abs_dae_Tolerances = double_value;
	iAssignedDaeAbsoluteTolerance = true;
}

void OpenSMOKE_Flame1D_DataManager::SetMaximumIntegrationOrder(const int int_value)
{
	max_integration_order = int_value;
}

void OpenSMOKE_Flame1D_DataManager::SetInitialTimeStep(const std::string units, const double value)
{
	initial_time_step = OpenSMOKE_Conversions::conversion_time(value, units);
}

void OpenSMOKE_Flame1D_DataManager::SetFuelRadialGradient(const std::string units, const double value)
{
	radialGradientC = OpenSMOKE_Conversions::conversion_frequency(value, units);
	iAssignedFuelRadialGradient = true;
}

void OpenSMOKE_Flame1D_DataManager::SetOxidizerRadialGradient(const std::string units, const double value)
{
	radialGradientO = OpenSMOKE_Conversions::conversion_frequency(value, units);
	iAssignedOxidizerRadialGradient = true;
}

void OpenSMOKE_Flame1D_DataManager::SetEnvironmentTemperature(const std::string units, const double value)
{
	environmentTemperature = OpenSMOKE_Conversions::conversion_temperature(value, units);
	iAssignedEnvironmentTemperature = true;
}

void OpenSMOKE_Flame1D_DataManager::SetGasRadiation()
{
	iGasRadiation = true;
	iAssignedGasRadiation = true;
}

void OpenSMOKE_Flame1D_DataManager::SetDepositionWall()
{
	iDepositionWall = true;
	//iThermophoreticEffect = true;
}

void OpenSMOKE_Flame1D_DataManager::SetUncorrectDiffusionFormulation()
{
	iCorrectDiffusionFormulation = false;
}

void OpenSMOKE_Flame1D_DataManager::SetSootRadiation(const std::string value)
{
	iAssignedSootRadiation = true;
	if (value == "none")			iRadiativeSootModel = RADIATIVE_SOOT_MODEL_NONE;
	else if (value == "fluent")		iRadiativeSootModel = RADIATIVE_SOOT_MODEL_FLUENT;
	else if (value == "bilger")		iRadiativeSootModel = RADIATIVE_SOOT_MODEL_BILGER;
	else if (value == "widmann")	iRadiativeSootModel = RADIATIVE_SOOT_MODEL_WIDMANN;
	else ErrorMessage("Wrong #SootRadiation option");
}

void OpenSMOKE_Flame1D_DataManager::SetGridRefineGradient(const int int_value)
{
	nDiff = int_value;
	iAssignedGridRefineGradient = true;
}

void OpenSMOKE_Flame1D_DataManager::SetGridRefineCurvature(const int int_value)
{
	nGrad = int_value;
	iAssignedGridRefineCurvature = true;
}

void OpenSMOKE_Flame1D_DataManager::SetGridRefineGradientStep(const double double_value)
{
	deltaDiff = double_value;
	iAssignedGridRefineGradientStep = true;
}

void OpenSMOKE_Flame1D_DataManager::SetGridRefineCurvatureStep(const double double_value)
{
	deltaGrad = double_value;
	iAssignedGridRefineCurvatureStep = true;
}

void OpenSMOKE_Flame1D_DataManager::SetDerivativeG(const char char_value)
{
	if (char_value != 'U' && char_value != 'C' && 
		char_value != 'B' && char_value != 'F')
			ErrorMessage("#DerivativeG options: U || C || B || F");

	iDerG = char_value;
	iAssignedDerivativeG = true;
}

void OpenSMOKE_Flame1D_DataManager::SetDerivativeT(const char char_value)
{
	if (char_value != 'U' && char_value != 'C' && 
		char_value != 'B' && char_value != 'F')
			ErrorMessage("#DerivativeT options: U || C || B || F");

	iDerT = char_value;
	iAssignedDerivativeT = true;
}

void OpenSMOKE_Flame1D_DataManager::SetDerivativeW(const char char_value)
{
	if (char_value != 'U' && char_value != 'C' && 
		char_value != 'B' && char_value != 'F')
			ErrorMessage("#DerivativeW options: U || C || B || F");

	iDerW = char_value;
	iAssignedDerivativeW = true;
}

void OpenSMOKE_Flame1D_DataManager::SetSensitivityOnFile(const vector<string> _names)
{
	GiveMeIndicesAndNames(*mix, _names, index_sensitivity, names_sensitivity);
	iAssignedSensitivity = true;
}

void OpenSMOKE_Flame1D_DataManager::SetFormationRatesOnFile(const vector<string> _names)
{
	GiveMeIndicesAndNames(*mix, _names, index_formation_rates, names_formation_rates);
	iAssignedFormationRates = true;
}

void OpenSMOKE_Flame1D_DataManager::SetVerboseROPAOnFile(const vector<string> _names)
{
	GiveMeIndicesAndNames(*mix, _names, index_ROPA, names_ROPA);
	iAssignedROPA			= true;
	iVerboseAssignedROPA	= true;
}

void OpenSMOKE_Flame1D_DataManager::SetROPAOnFile()
{
	iAssignedROPA			= true;
}

void OpenSMOKE_Flame1D_DataManager::SetElementFluxAnalysisOnFile(const vector<string> _names)
{
	element_flux_analysis = new OpenSMOKE_ElementFluxAnalysis();
	element_flux_analysis->Initialize(mix, _names);
}

void OpenSMOKE_Flame1D_DataManager::SetGlobalElementFluxAnalysisOnFile(const vector<string> _names)
{
	element_flux_analysis_manager = new OpenSMOKE_ElementFluxAnalysisManager();
	element_flux_analysis_manager->Initialize(_names);
	iAssignedGlobalElementFluxAnalysis = true;
}


void OpenSMOKE_Flame1D_DataManager::SetExperimentOnFile(const vector<string> _names)
{
	GiveMeIndicesAndNames(*mix, _names, index_Experiment, names_Experiment);
	iAssignedExperiment = true;
}

void OpenSMOKE_Flame1D_DataManager::SetReactionRatesOnFile(const vector<string> _names)
{
	if (_names[0] == "ALL")
	{
		ChangeDimensions(mix->NumberOfReactions(), &index_reaction_rates);
		for(int k=1;k<=mix->NumberOfReactions();k++)
		{
			index_reaction_rates[k] = k;
			names_reaction_rates.push_back(mix->strReaction[index_reaction_rates[k]]);
		}
	}
	else
	{
		ChangeDimensions(_names.size(), &index_reaction_rates);
		for(int k=0;k<int(_names.size());k++)
		{
			stringstream index_string;
			index_string << _names[k];
			index_string >> index_reaction_rates[k+1];
			if (index_reaction_rates[k+1] > mix->NumberOfReactions())
				ErrorMessage("The requested reaction for the reaction rates analysis is not available: " + _names[k]);
			names_reaction_rates.push_back(mix->strReaction[index_reaction_rates[k+1]]);
		}
	}
	
	iAssignedReactionRates = true;
}

void OpenSMOKE_Flame1D_DataManager::SetBINDensities(const vector<string> _names)
{
	if (_names.size() != 4)
		ErrorMessage("Wrong definition of BIN densities");

	bin_index_zero = atoi(_names[0].c_str());
	bin_index_final = atoi(_names[2].c_str());
	bin_density_A = atof(_names[1].c_str());
	bin_density_B = atof(_names[3].c_str());
}



void OpenSMOKE_Flame1D_DataManager::AssignMassFlowRate(const std::string units, const double value)
{
	MassFlowRate = OpenSMOKE_Conversions::conversion_massFlowRate(value, units);
}

void OpenSMOKE_Flame1D_DataManager::AssignInletVelocity(const std::string units, const double value)
{
	AssignFuelVelocity(units, value);
}
	
void OpenSMOKE_Flame1D_DataManager::AssignInletTemperature(const std::string units, const double value)
{
	AssignFuelTemperature(units, value);
}
	
void OpenSMOKE_Flame1D_DataManager::AssignCrossSection(const std::string units, const double value)
{
	CrossSectionalArea = OpenSMOKE_Conversions::conversion_area(value, units);
}

void OpenSMOKE_Flame1D_DataManager::AssignInletMassFractions(const vector<string> names, const vector<double> values)
{
	AssignFuelMassFractions(names, values);
}
 
void OpenSMOKE_Flame1D_DataManager::AssignInletMoleFractions(const vector<string> names, const vector<double> values)
{
	AssignFuelMoleFractions(names, values);
}

void OpenSMOKE_Flame1D_DataManager::SetOutletTemperature(const std::string units, const double value)
{
	AssignOxidizerTemperature(units, value);
	iAssignedOutletTemperature = true;
}

void OpenSMOKE_Flame1D_DataManager::SetOutletMassFractions(const vector<string> names, const vector<double> values)
{
	AssignOxidizerMassFractions(names, values);
	iAssignedOutletMassFractions = true;
}

void OpenSMOKE_Flame1D_DataManager::SetOutletMoleFractions(const vector<string> names, const vector<double> values)
{
	AssignOxidizerMoleFractions(names, values);
	iAssignedOutletMoleFractions = true;
}

void OpenSMOKE_Flame1D_DataManager::SetFlameSpeedIndex(const int int_value)
{
	iFixedTemperature = int_value;
	iAssignedFlameSpeedIndex = true;
}

void OpenSMOKE_Flame1D_DataManager::SetFlameSpeedTemperature(const std::string units, const double value)
{
	fixedTemperature = OpenSMOKE_Conversions::conversion_temperature(value, units);
	iAssignedFlameSpeedTemperature = true;
}

void OpenSMOKE_Flame1D_DataManager::SetStretchingFactor(const double double_value)
{
	alfa = double_value;
	iAssignedStretchingFactor = true;
}

void OpenSMOKE_Flame1D_DataManager::SetCrossSectionProfile(const std::string string_value)
{
	iAreaProfile = 1;
	ud_cross_section_profile_file_name = string_value;
	iAssignedCrossSectionProfile = true;
}

void OpenSMOKE_Flame1D_DataManager::SetCrossSectionMIT()
{
	iAreaProfile = 2;
	iAssignedCrossSectionMIT = true;
}

void OpenSMOKE_Flame1D_DataManager::SetSoretEffect()
{
	iSoretEffect = true;
}

void OpenSMOKE_Flame1D_DataManager::SetPhysicalSootDiffusionCoefficients(const int value)
{
	iPhysicalSootDiffusionCoefficients = value;
}

void OpenSMOKE_Flame1D_DataManager::SetThermophoreticEffect()
{
	iThermophoreticEffect = true;
}

void OpenSMOKE_Flame1D_DataManager::SetLewisNumbers()
{
	iLewisNumbers = true;
}

void OpenSMOKE_Flame1D_DataManager::SetUnityLewisNumbers()
{
	iUnityLewisNumbers = true;
}

void OpenSMOKE_Flame1D_DataManager::SetUserDefinedLewisNumbers(const vector<string> names, const vector<double> values)
{
	iUserDefinedLewisNumbers = true;
	
	ChangeDimensions(mix->NumberOfSpecies(), &user_defined_lewis_numbers);
	user_defined_lewis_numbers = 1.;

	for(int i=1;i<=names.size();i++)
	{
		if (names[i-1] == "default" || names[i-1] == "DEFAULT" || names[i-1] == "Default")
			user_defined_lewis_numbers = values[i-1];
	}

	for(int i=1;i<=names.size();i++)
	{
		if (names[i-1] == "default" || names[i-1] == "DEFAULT" || names[i-1] == "Default")
		{
			continue;
		}
		else
		{
			int index = mix->recognize_species(names[i-1]);
			user_defined_lewis_numbers[index]	 = values[i-1];
		}
	}
}

void OpenSMOKE_Flame1D_DataManager::SetVerboseMixtureProperties()
{
	iVerboseMixtureProperties = true;
}

void OpenSMOKE_Flame1D_DataManager::SetVerboseFluxes()
{
	iVerboseFluxes = true;
}

void OpenSMOKE_Flame1D_DataManager::SetCorrectionReactionRates(const double value)
{
	iCorrectionReactionRates	= true;
	correctionReactionRates		= value;
}

void OpenSMOKE_Flame1D_DataManager::SetTurbulentDiffusivity(const double value)
{
	iTurbulentDiffusivity = true;
	diffusivityEnhancingFactor = value;
}

void OpenSMOKE_Flame1D_DataManager::AssignEquivalenceRatioForPremixedFlames(const double value)
{
	iEquivalenceRatioForPremixedFlames = true;
	equivalence_ratio = value;
}

void OpenSMOKE_Flame1D_DataManager::AssignFuelMassFractionsForPremixedFlames(const vector<string> names, const vector<double> values)
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
	AssignFuelMoleFractionsForPremixedFlames(_names, _values);
}
 
void OpenSMOKE_Flame1D_DataManager::AssignFuelMoleFractionsForPremixedFlames(const vector<string> names, const vector<double> values)
{
	fuel_names.resize(names.size()+1);
	ChangeDimensions(names.size(), &moles_fuel);
	for(int i=1;i<=int(names.size());i++)
	{
		fuel_names[i] = names[i-1];
		moles_fuel[i] = values[i-1];
	}
}

void OpenSMOKE_Flame1D_DataManager::AssignOxidizerMassFractionsForPremixedFlames(const vector<string> names, const vector<double> values)
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
	AssignOxidizerMoleFractionsForPremixedFlames(_names, _values);
}
 
void OpenSMOKE_Flame1D_DataManager::AssignOxidizerMoleFractionsForPremixedFlames(const vector<string> names, const vector<double> values)
{
	oxidizer_names.resize(names.size()+1);
	ChangeDimensions(names.size(), &moles_oxidizer);
	for(int i=1;i<=int(names.size());i++)
	{
		oxidizer_names[i] = names[i-1];
		moles_oxidizer[i] = values[i-1];
	}
}

void OpenSMOKE_Flame1D_DataManager::AssignInletStreamCompositionForPremixedFlames()
{
	if (iEquivalenceRatioForPremixedFlames == true)
	{
		OpenSMOKE_GasStream gas_stream;
		gas_stream.AssignKineticScheme(*flame->mix);
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

		AssignFuelMoleFractions(names, values);
	}
}

void OpenSMOKE_Flame1D_DataManager::SetSingleContributions(const vector<string> _names)
{
	iSingleContributions = true;
	GiveMeIndicesAndNames(*mix, _names, index_SingleContributions, names_SingleContributions);
}

void OpenSMOKE_Flame1D_DataManager::SetChangeFrequencyFactor(const vector<string> _names)
{
	const int n = _names.size();
	for(int j=1;j<=n;j+=2)
	{
		cout << "Change frequency factor: " << atoi( _names[j-1].c_str() ) << " " << atof( _names[j].c_str() )*100. << "%" << endl;
		mix->kinetics.ChangeFrequencyFactor(atoi(_names[j-1].c_str()), atof(_names[j].c_str()));
	}
}

void OpenSMOKE_Flame1D_DataManager::SetChangeActivationEnergy(const vector<string> _names)
{
	const int n = _names.size();
	for(int j=1;j<=n;j+=3)
	{
		if (_names[j+1] != "cal/mol")
			ErrorMessage("Wrong units in #ChangeActivationEnergy option!");
		mix->kinetics.ChangeActivationEnergy(atoi(_names[j-1].c_str()), atof( _names[j].c_str() ));				// cal/mol
	}
}

void OpenSMOKE_Flame1D_DataManager::SetChangeDiffusivity(const vector<string> _names)
{
	iCorrectionDiffusivity = true;

	const int n = _names.size();
	for(int j=1;j<=n;j+=3)
	{
		cout << "Change diffusivity: " << _names[j-1].c_str() << " " << _names[j].c_str() << " " << atof( _names[j+1].c_str() )*100. << "%" << endl;
		index_correction_diffusivity.Append(mix->recognize_species(_names[j-1]));
		index_correction_diffusivity.Append(mix->recognize_species(_names[j]));
		correction_diffusivity.Append(1.+atof(_names[j+1].c_str()));
	}
}

void OpenSMOKE_Flame1D_DataManager::SetChangeFakeTemperatureThermalConductivity(const double value, const std::string units)
{
	iFakeTemperatureThermalConductivity = true;
	fakeTemperatureThermalConductivityIncrement = OpenSMOKE_Conversions::conversion_temperature(value, units);
}

void OpenSMOKE_Flame1D_DataManager::SetChangeFormationEnthalpy(const vector<string> _names)
{
	iCorrectionFormationEnthalpy = true;

	const int n = _names.size();
	for(int j=1;j<=n;j+=3)
	{
		cout	<< "Change formation enthalpy: " << _names[j-1].c_str() << " " 
			<< OpenSMOKE_Conversions::conversion_specificEnergyMolar(atof( _names[j].c_str() ), _names[j+1]) << " J/kmol" << endl;
		index_correction_formation_enthalpy.Append(mix->recognize_species(_names[j-1]));
		correction_formation_enthalpy.Append( OpenSMOKE_Conversions::conversion_specificEnergyMolar(atof( _names[j].c_str() ), _names[j+1]) );
	}
}

void OpenSMOKE_Flame1D_DataManager::SetLennardJonesMode(const std::string option)
{
	if (option == "ALL")
		sensitivityLennardJonesMode = OPENSMOKE_FITTING_ALL;
	else if (option == "ONLY_DIFFUSIVITIES")
		sensitivityLennardJonesMode = OPENSMOKE_FITTING_ONLY_DIFFUSIVITIES;
	else if (option == "ONLY_CONDUCTIVITY")
		sensitivityLennardJonesMode = OPENSMOKE_FITTING_ONLY_CONDUCTIVITY;
	else if (option == "ONLY_VISCOSITY")
		sensitivityLennardJonesMode = OPENSMOKE_FITTING_ONLY_VISCOSITY;
	else
		ErrorMessage("Wrong LennardJonesMode option...");
}


void OpenSMOKE_Flame1D_DataManager::SetAdaptiveGrid(const std::string name)
{
	if		(name=="GRADIENT")		kind_of_adaptive_grid = ADAPTIVE_GRADIENT;
	else if (name=="CHEMKIN")		kind_of_adaptive_grid = ADAPTIVE_CHEMKIN;
	else if (name=="CHEMKIN_SQRT")	kind_of_adaptive_grid = ADAPTIVE_CHEMKIN_SQRT;
	else if (name=="CHEMKIN_POW")	kind_of_adaptive_grid = ADAPTIVE_CHEMKIN_POW;
	else if (name=="EISEMAN")		kind_of_adaptive_grid = ADAPTIVE_EISEMAN;
	else ErrorMessage("Wrong kind of adaptive grid...");
}

void OpenSMOKE_Flame1D_DataManager::SetAdaptiveGridCoefficients(const vector<string> values)
{
	if (values.size()!=3)	ErrorMessage("Wrong number of adaptive grid coefficients...");
	
	iAssignedAdaptiveGridCoefficients = true;
	adaptive_grid_alfa  = atof(values[0].c_str()); 
	adaptive_grid_beta  = atof(values[1].c_str()); 
	adaptive_grid_gamma = atof(values[2].c_str()); 
}

void OpenSMOKE_Flame1D_DataManager::SetPoolFireGridOptions(const vector<string> values)
{
	if (values.size()!=8)	ErrorMessage("Wrong number of pool fire grid coefficients...");
	
	if (values[0] != "alfa_fuel")			ErrorMessage("Expected alfa_fuel - Found: " + values[0]);
	if (values[2] != "alfa_oxidizer")		ErrorMessage("Expected alfa_oxidizer - Found: " + values[2]);
	if (values[4] != "point_fraction")		ErrorMessage("Expected point_fraction - Found: " + values[4]);
	if (values[6] != "distance_fraction")	ErrorMessage("Expected distance_fraction - Found: " + values[6]);
	poolfire_grid_alfa_fuel         = atof(values[1].c_str()); 
	poolfire_grid_alfa_oxidizer     = atof(values[3].c_str()); 
	poolfire_grid_point_fraction    = atof(values[5].c_str()); 
	poolfire_grid_distance_fraction = atof(values[7].c_str()); 
}


void OpenSMOKE_Flame1D_DataManager::SetPoolFire(const vector<string> values)
{
	// Assign Pool Fire composition
	{
		vector<string> names;	names.push_back(nameFuel);
		vector<double> values;	values.push_back(1.);
		AssignFuelMassFractions(names, values);
	}

	int nsize = values.size();


	if (values[1] == "Equilibrium")		iPoolFire = POOL_FIRE_EQUILIBRIUM;
	if (values[1] == "Tboiling")		iPoolFire = POOL_FIRE_TASSIGNED;
	if (values[1] == "LiquidPool")		iPoolFire = POOL_FIRE_LIQUIDPOOL;

	OpenSMOKE_LiquidProperties_Database *liquid_properties = new OpenSMOKE_LiquidProperties_Database();
	liquid_properties->ReadFromFolder(values[0]);

	pool_fire_liquid_species = new OpenSMOKE_LiquidSpecies();
	pool_fire_liquid_species->SetName(nameFuel);
	pool_fire_liquid_species->SetProperties(*liquid_properties);
	pool_fire_liquid_species->Summary();
	pool_fire_temperature = pool_fire_liquid_species->TNormalBoiling(P_Pascal);

	if (iPoolFire == POOL_FIRE_TASSIGNED)
	{
		if (values[1] == "Tboiling")
			{
				// Do nothing
			}
			else
			{
				if (nsize==2)
					ErrorMessage("#PoolFire option requires the pool temperature when non-equilibrium model is adopted!");
				else
					pool_fire_temperature = OpenSMOKE_Conversions::conversion_temperature(atof(values[1].c_str()), values[2]);
			}		
	}

	if (iPoolFire == POOL_FIRE_EQUILIBRIUM)
		if (nsize > 2)	ErrorMessage("#PoolFire option requires only 2 arguments when Equilibrium model is adopted!");

	if (iPoolFire == POOL_FIRE_LIQUIDPOOL)
	{
		if (nsize != 6)	ErrorMessage("#PoolFire option requires 6 arguments when LiquidPool model is adopted!");
		pool_fire_feed_temperature = OpenSMOKE_Conversions::conversion_temperature(atof(values[2].c_str()), values[3]);
		pool_fire_depth = OpenSMOKE_Conversions::conversion_length(atof(values[4].c_str()), values[5]);
	}

	// First guess values
	{
		double dT_over_dx = 7.e5;	// [K/m]
		double rho_gas = P_Pascal/Constants::R_J_kmol/pool_fire_temperature*mix->M(jFUEL);
		mix->SpeciesConductivityFromFitting(pool_fire_temperature);
		double lambda_gas = mix->lambda[jFUEL];
		double DHvap = pool_fire_liquid_species->Hv(pool_fire_temperature);
		AssignFuelTemperature("K", pool_fire_temperature);
		AssignFuelVelocity("m/s", lambda_gas/rho_gas/DHvap*dT_over_dx);
	}
}		

void OpenSMOKE_Flame1D_DataManager::SetRobustTemperature()
{
	iRobustTemperature = true;
}

void OpenSMOKE_Flame1D_DataManager::SetPoolFireCorrectionFactorVaporizationHeat(const double value)
{
	correctionFactorVaporizationHeat = value;
}

void OpenSMOKE_Flame1D_DataManager::SetPoolFireCorrectionFactorVaporPressure(const double value)
{
	correctionFactorVaporPressure = value;
}

void OpenSMOKE_Flame1D_DataManager::SetPoolFireCorrectionFactorThermalConductivity(const double value)
{
	correctionFactorThermalConductivity = value;
}

void OpenSMOKE_Flame1D_DataManager::SetPoolFireCorrectionFactorSpecificHeat(const double value)
{
	correctionFactorSpecificHeat = value;
}
