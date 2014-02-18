/***************************************************************************
 *   Copyright (C) 2007 by Alberto Cuoci   	                               *
 *   alberto.cuoci@polimi.it   						                       *
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

#include "qmom/OpenSMOKE_QMOM_Module.h"


void OpenSMOKE_QMOM_Module::setupGasMixture(OpenSMOKE_ReactingGas &mix)
{
	int j;

	// Initialize indices
	jC2H2	= 0;
	jO2		= 0;
	jOH		= 0;

	// Recognize species
	for(j=1;j<=mix.NumberOfSpecies();j++)
	{
		if(mix.names[j] == "C2H2")
			jC2H2 = j;
		if(mix.names[j] == "O2")
			jO2 = j;
		if(mix.names[j] == "OH")
			jOH = j;
	}

	// Checking
	if (jC2H2==0)
		MessageError("It was impossible to recognize C2H2 in the kinetic scheme!");
	if (jO2==0)
		MessageError("It was impossible to recognize O2 in the kinetic scheme!");
	if (jOH==0)
		MessageError("It was impossible to recognize OH in the kinetic scheme!");
}


void OpenSMOKE_QMOM_Module::initialize(int _N, int _iFractalDimension, 
									  int _iNucleation, int _iGrowth, int _iOxidation, 
									  int _iAggregation, double _L0, 
									  int _iNucleationDistribution, double _seed_density,
									  int _iTfluctuations)
{
	// --------------------------------
	// DQMOM - Options and Data
	// --------------------------------	
	N	= _N;
	L0	= _L0;
	iFractalDimension	= _iFractalDimension;
	iNucleation			= _iNucleation;
	iGrowth				= _iGrowth;
	iOxidation			= _iOxidation;
	iAggregation		= _iAggregation;
	seed_density		= _seed_density;
	iNucleationDistribution = _iNucleationDistribution;
	epsilon				= 2.e-9*L0;
	iTfluctuations		= _iTfluctuations;

	finalize_initialization();
}

void OpenSMOKE_QMOM_Module::finalize_initialization()
{
	// --------------------------------
	// DQMOM - Memory allocation
	// --------------------------------	
	ChangeDimensions(N, &w);
	ChangeDimensions(N, &csi);
	ChangeDimensions(2*N, &moments);
	ChangeDimensions(N, &L);

	// --------------------------------
	// DQMOM - Physical Model Setup
	// --------------------------------	
	qmom.setup(N);
	models.setup(N);

	// ---------------------------------------------------------
	// DQMOM - Nucleation Distribution
	// ---------------------------------------------------------
	if (iNucleationDistribution == 0) // Uniform
	{
		moments[1] = 1.;
		for(int i=2;i<=N*2;i++)
			moments[i] = pow(2.*L0, i-1)/double(i);
		qmom.giveWeigthsAndAbscissasFromMoments( moments, w, L);
		L /= (2.*L0);							
	}
	else 
		MessageError("Only the uniform distribution is available in this version!");

	// -------------------------------------------------------------
	// SOOT Models - Physical Model Setup
	// -------------------------------------------------------------
	// Checking
	if (iFractalDimension<0 || iFractalDimension>1)
		MessageError("Wrong choice of FractalDimension Option!!!");
	if (iNucleation<1 || iNucleation>7)
		MessageError("Wrong choice of nucleation model!!!");
	if ( iGrowth<0 || iGrowth>7   )
		MessageError("Wrong choice of growth model!!!");
	if (iOxidation<0 || iOxidation>3)
		MessageError("Wrong choice of oxidation model!!!");
	if (iAggregation<0 || iAggregation>1)
		MessageError("Wrong choice of aggregation model!!!");
	if (seed_density<=1.e0)
		MessageError("Wrong choice of seed density!!!");
	if (iNucleationDistribution!=0)
		MessageError("Wrong choice of nucleation distribution!!!");
	if (iTfluctuations<0 || iTfluctuations>1)
		MessageError("Wrong choice of temperature fluctuations!!!");

	// Initialize
	soot.initialize(N,	iFractalDimension, 
						iNucleation, iGrowth, iOxidation, iAggregation, L0, iTfluctuations);
}

void OpenSMOKE_QMOM_Module::updateMoments(BzzVector &moments)
{
	// ---------------------------------------------------------------
	// Gordon Algorithm to obtain weights and abscissas
	// ---------------------------------------------------------------		
	qmom.giveWeigthsAndAbscissasFromMoments( moments, w, csi);
		
	// ---------------------------------------------------------------
	// QMOM Updating
	// ---------------------------------------------------------------
	w   *= seed_density;
	csi *= epsilon;
	models.update( w, csi );
}

void OpenSMOKE_QMOM_Module::updateData(double &T, double &P, double &rho, 
									  double &mu, double &tr,
									  double &omega_O2, double &omega_C2H2, 
									  double &omega_OH, double &SNTV)
{

	// ---------------------------------------------------------------
	// SootModels Updating
	// ---------------------------------------------------------------
	soot.update(w, csi);
	soot.assign_physical_properties (	T, P, rho, mu, tr, 
										omega_O2, omega_C2H2, omega_OH, SNTV);
}

BzzVector OpenSMOKE_QMOM_Module::calculateSources()
{	
	// 1. Nucleation phenomena
	// ---------------------------------------------------------------
	double epsilon;
	double J;
	soot.nucleation_rate(J, epsilon);
	models.nucleation(epsilon, J, L);

	// 2. Growth phenomena
	// ---------------------------------------------------------------
	if (iGrowth!=0)
	{
		double Ggrowth, rGrowth;
		soot.growth_rate(Ggrowth, rGrowth);
		models.powerLawGrowth(Ggrowth, rGrowth);
	}

	// 3. Oxidation phenomena
	// ---------------------------------------------------------------
	if (iOxidation!=0)
	{
		double Goxidation, rOxidation;
		soot.oxydation_rate(Goxidation, rOxidation);
		models.powerLawGrowth(Goxidation, rOxidation);
	}
		
	// 4. Aggregation phenomena
	// ---------------------------------------------------------------
	if (iAggregation!=0)
	{
		BzzMatrix Kernel(N,N);
		soot.aggregation_kernel(Kernel);
		models.aggregationKernel_Soot(Kernel);
		models.homogeneousAggregation();
	}

	// Conversion of sources for moment equations
	// ---------------------------------------------------------------
	// The source terms are now normalized
	double coeff = seed_density;
	for(int j=1;j<=2*N;j++)
	{
		models.Source[j] /= coeff;
		coeff *= epsilon;
	}

	return models.Source;
}

void OpenSMOKE_QMOM_Module::MessageError(char *message)
{
	cout << endl;
	cout << "QMOM_Module Error: " << message << endl; 
	cout << "Press enter to continue..." << endl;
	getchar();
	exit(-1);
}

void OpenSMOKE_QMOM_Module::readFromFile(const string fileName)
{
	const int SIZE = 200;
	char comment[SIZE];
	ifstream fInput;

	openInputFileAndControl(fInput, fileName);

	// RIGHE DI COMMENTO
	fInput.getline(comment, SIZE);
	fInput.getline(comment, SIZE);
	fInput.getline(comment, SIZE);
	
	// NUMBER OF DELTA OF DIRAC
	fInput >> N;							fInput.getline(comment, SIZE);
	fInput >> iFractalDimension;			fInput.getline(comment, SIZE);
	fInput >> iNucleation;					fInput.getline(comment, SIZE);
	fInput >> iGrowth;						fInput.getline(comment, SIZE);
	fInput >> iOxidation;					fInput.getline(comment, SIZE);
	fInput >> iAggregation;					fInput.getline(comment, SIZE);
	fInput >> L0;							fInput.getline(comment, SIZE);
	fInput >> seed_density;					fInput.getline(comment, SIZE);
	fInput >> iNucleationDistribution;		fInput.getline(comment, SIZE);
	fInput >> iTfluctuations;				fInput.getline(comment, SIZE);
	
	soot.read_correction_coefficients_from_file(fInput);

 	fInput >> comment;						

	if (strcmp(comment, "END"))
	{
		cout << "ERROR: reading QMOM Input File!!" << comment << endl;
		exit(1);
	}

	fInput.close();

	// Derived Properties
	// ----------------------------------------------------------------------
	epsilon				= 2.e-9*L0;
	finalize_initialization();
}
