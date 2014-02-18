/***************************************************************************
 *   Copyright (C) 2006 by Alberto Cuoci   	                               *
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

#include "basic/OpenSMOKE_Constants.h"
#include "engine/OpenSMOKE_IdealGas.h"

const int			OpenSMOKE_IdealGas::N_LOOKUP	= 4000;								// [-]
const double		OpenSMOKE_IdealGas::logCATM		= log(101325./Constants::R_J_kmol);	// [kmol.K/m3]


OpenSMOKE_IdealGas::OpenSMOKE_IdealGas()
{
	TMIN		= 250.;
	TMAX		= 3200.;

	name_object   = "[not assigned]";
	name_place    = "[not assigned]";
	name_author   = "[not assigned]";
	building_date = GiveMeTimeAndDate();

	iMixViscositySimplified = 0;
}

void OpenSMOKE_IdealGas::ErrorMessage(const string message)
{
    cout << endl;
    cout << "Class:  OpenSMOKE_IdealGas"	<< endl;
    cout << "Object: " << name_object		<< endl;
    cout << "Error:  " << message           << endl;
    cout << "Press a key to continue... "   << endl;
    getchar();
    exit(-1);
}

void OpenSMOKE_IdealGas::WarningMessage(const string message)
{
    cout << endl;
    cout << "Class:  OpenSMOKE_IdealGas"		<< endl;
    cout << "Object: "		<< name_object	<< endl;
    cout << "Warning:  "	<< message      << endl;
    cout << "Press a key to continue... "   << endl;
    getchar();
}


void OpenSMOKE_IdealGas::SetMinimumTemperature(const double tmin)
{
	TMIN		= tmin;
	if (TMIN >= TMAX || TMIN <= 0.)
		ErrorMessage("Please check your minimum temperature!");
}

void OpenSMOKE_IdealGas::SetMaximumTemperature(const double tmax)
{
	TMAX		= tmax;
	if (TMIN >= TMAX || TMAX >= 6000.)
		ErrorMessage("Please check your maximum temperature!");
}

void OpenSMOKE_IdealGas::SetName(const string name)
{
	name_object = name;
}

void OpenSMOKE_IdealGas::SetAuthorName(const string name)
{
	name_author = name;
}

void OpenSMOKE_IdealGas::SetPlaceName(const string name)
{
	name_place = name;
}

void OpenSMOKE_IdealGas::Setup(BzzLoad &binaryFile)
{
	int k;
	char dummy[Constants::NAME_SIZE];

	cout << "--------------------------------------------------------------------------" << endl;
	cout << "                 THERMODYNAMICS AND TRANSPORT PROPERTIES                  " << endl;
	cout << "--------------------------------------------------------------------------" << endl;

	binaryFile.fileLoad.read((char*) dummy, sizeof(dummy));
	string version = dummy;
	if (version != "V090524")
		ErrorMessage("This version of thermodynamics.bin file is not supported: " + version);
	cout << "Version: " << version << endl;

	// Reading information
	binaryFile.fileLoad.read((char*) dummy, sizeof(dummy));
	cout << "Name:    "  << dummy << endl;
	binaryFile.fileLoad.read((char*) dummy, sizeof(dummy));
	cout << "Author:  "  << dummy << endl;
	binaryFile.fileLoad.read((char*) dummy, sizeof(dummy));
	cout << "Place:   "  << dummy << endl;
	binaryFile.fileLoad.read((char*) dummy, sizeof(dummy));
	cout << "Date:    "  << dummy << endl;

	// Minimum and maximum temperature
	binaryFile >> TMIN;			cout << "Minimum temperature: " << TMIN << " K" << endl;
	binaryFile >> TMAX;			cout << "Maximum temperature: " << TMAX << " K" << endl;
	SetMinimumTemperature(TMIN);
	SetMaximumTemperature(TMAX);

	// Reading number of Species
	binaryFile >> NC;

	// Lettura dei nomi delle specie e del numero di componenti
	char name[Constants::NAME_SIZE];
	names = new string[NC + 1];
	for(k=1;k<=NC;k++)
	{
		binaryFile.fileLoad.read((char*) name, sizeof(name));
		names[k] = name;
	}

	// Allocazione della memoria per le variabili
	Allocate();

	// Lettura delle informazioni termodinamiche
	for(k=1;k<=NC;k++)
	{
		binaryFile >> CpHT[k];
		binaryFile >> CpLT[k];
		binaryFile >> aDH[k];
		binaryFile >> bDH[k];
		binaryFile >> aDS[k];
		binaryFile >> bDS[k];
	}
	binaryFile >> M;
	binaryFile >> T1;
	binaryFile >> T2;
	binaryFile >> T3;

	// Read elements
	binaryFile >> elements;
	binaryFile >> m_elements;

	// Read from file: element list
	for(k=1;k<=elements.Rows();k++)
	{
		char name[Constants::NAME_SIZE];
		binaryFile.fileLoad.read((char*) name, sizeof(name));
		list_of_elements.push_back(name);
	}

	// Post Processing information
	for(int i = 1;i <= NC;i++)
		uM[i] = 1. / M[i];			// reciproco del peso molecolare di ciascun componente

	if (iMixViscositySimplified == false)
	{
		int i;

		// Coefficiente costante nella formula per il calcolo della viscosita di miscela
		i=1;
		for(k = 1;k <= NC;k++)
			for(int j = k+1;j <= NC;j++)
			{
				phi_eta_sup[i]=1./sqrt( 8.*(1.+M[k]/M[j]) );
				phi_eta_inf[i++]=1./sqrt( 8.*(1.+M[j]/M[k]) );
			}

		// Rapporto tra i pesi molecolari elevato alla 1/4
		i=1;
		for(k = 1;k <= NC;k++)
			for(int j = k+1;j <= NC;j++)
				PMRatio1su4[i++] = sqrt(sqrt(M[j]/M[k]));
	}
	else
	{
		// Rapporto tra i pesi molecolari elevato alla 1/2
		int i=1;
		for(k = 1;k <= NC;k++)
			for(int j = k+1;j <= NC;j++)
			{
				sqrtPMRatio_sup[i] = sqrt(M[j]/M[k]);
				sqrtPMRatio_inf[i++] = sqrt(M[k]/M[j]);
			}
	}

	// Soot and PAH manager
    soot_manager.recognizeSpecies(NC, names, M);
    pah_manager.recognizeSpecies(NC,  names, M);
}

// ************************************************************************************	//
// ************************************************************************************ //
//								PROPRIETA' DI MISCELA									//
// ************************************************************************************	//
// ************************************************************************************ //

double OpenSMOKE_IdealGas::MixViscosity_FromMolarFractions(BzzVector &x)
{
	int i, j, k;

	double etaMix = 0.;
	sumK=x;

	// Formula di Wilke - Journal of Chemical Physics 18:517 (1950)
	// Modificata da Bird, Stewart, Lightfoot - Transport phenomena (1960)
	// Riportata in: Reid, Prausnitz, Poling - The properties of gases and liquids, pag 407

	if(!iMixViscositySimplified)
	{
		i=1;
		for(k = 1;k <= NC;k++)
		{
			sqrtEta[k] = sqrt(eta[k]);
			usqrtEta[k] = 1./sqrtEta[k];
		}

		for(k = 1;k <= NC;k++)
			for(j = k+1;j <= NC;j++)
			{
				delta_phi[i] = sqrtEta[k]*usqrtEta[j]*PMRatio1su4[i];
				sumK[k] += x[j]*phi_eta_sup[i] * BzzPow2(1.+delta_phi[i]);	// F.(49)
				sumK[j] += x[k]*phi_eta_inf[i] * BzzPow2(1.+1./delta_phi[i]);	// F.(49)
				i++;
			}
	}

	// Formula di Herning and Zipperer
	// Riportata in: Reid, Prausnitz, Poling - The properties of gases and liquids, pag 410
	// - ?una formula semplificata della precedente ma molto meno precisa; ha il vantaggio
	//   di essere decisamente meno onerosa

	else
	{
		i=1;
		for(k = 1;k <= NC;k++)
			for(j = k+1;j <= NC;j++)
			{
				sumK[k] += x[j]*sqrtPMRatio_sup[i];	// F.(49)
				sumK[j] += x[k]*sqrtPMRatio_inf[i];	// F.(49)
				i++;
			}
	}

	for(k = 1;k <= NC;k++)
			etaMix += x[k]*eta[k]/sumK[k];				// F.(48)


	return etaMix;
}

double OpenSMOKE_IdealGas::MixConductivity_FromMolarFractions(BzzVector &x)
{
	// Calcolo della conducibilita della miscela
	// Formula di Mathur, Todor, Saxena - Molecular Physics 52:569 (1967)

	int k;
	double sum1 = 0.,
		   sum2 = 0.;

	for(k = 1;k <= NC;k++)
	{
		sum1+=x[k]*lambda[k];
		sum2+=x[k]/lambda[k];
	}

	return (0.50 * ( sum1 + 1./sum2 ));	// F.(50)
}

void OpenSMOKE_IdealGas::MixDiffusionCoefficient_FromMolarFractionsAndPM(BzzVector &Dkm, BzzVector &x)
{
	// Calcolo dei coefficienti di diffusione materiali di miscela F.(44)
	// Devono essere note le frazioni molari X
	int i,j,k;
	double	sum;
	double soglia = 1.e-14;

//	for(i = 1;i <= NC;i++)
//		x[i] = (x[i]<0.)? 0. : x[i];

	sum=1./(1. + soglia*NC);
	for(i = 1;i <= NC;i++)
		x_aux[i]= (x[i]+soglia)*sum;

	double massaMolareMedia=0.;
	for(i = 1;i <= NC;i++)
		massaMolareMedia+=x_aux[i]*M[i];

	// Calcolo del denominatore e del coefficiente di diffusione della miscela
	for(k = 1;k <= NC;k++)
	{
		sum = 0.;

		for(j = 1;j < k;j++)
			sum += x_aux[j] * Djk[j][k];

		for(j = k+1;j <= NC;j++)
			sum += x_aux[j] * Djk[k][j];

		Dkm[k] = (massaMolareMedia - x_aux[k]*M[k]) / (massaMolareMedia * sum);
	}

	x = x_aux;
}


// Entalpia e Entropia delle singole specie [J,mol,K]
void OpenSMOKE_IdealGas::SpeciesEnthalpyAndEntropy(double T)
{
	double logT	=	log(T);
	double uT	=	1./T;

	for(int i = 1;i <= NC;i++)
	{
		if(T > T2[i])	// High Temperature
		{
			h[i] = aDH[i][1]+T*(aDH[i][2]+T*(aDH[i][3]+T*(aDH[i][4]+T*aDH[i][5])))+aDH[i][6]*uT;
			s[i] = aDS[i][1]*logT +
				   T*( aDS[i][2] + T*( aDS[i][3] + T *(aDS[i][4] + T*aDS[i][5])))+ aDS[i][6];
		}

		else			// Low Temperature
		{
			h[i] = bDH[i][1]+T*(bDH[i][2]+T*(bDH[i][3]+T*(bDH[i][4]+T*bDH[i][5])))+bDH[i][6]*uT;
			s[i] = bDS[i][1]*logT +
				   T*( bDS[i][2] + T*( bDS[i][3] + T *(bDS[i][4] + T*bDS[i][5]))) + bDS[i][6];
		}
	}
}

// Entalpia delle singole specie [J/mol]
void OpenSMOKE_IdealGas::SpeciesEnthalpy(double T)
{
	double uT = 1./T;

	for(int i = 1;i <= NC;i++)
	{
		if(T > T2[i])	// High Temperature
			h[i] = aDH[i][1]+T*(aDH[i][2]+T*(aDH[i][3]+T*(aDH[i][4]+T*aDH[i][5]))) +
				   aDH[i][6]*uT;
		else			// Low Temperature
			h[i] = bDH[i][1]+T*(bDH[i][2]+T*(bDH[i][3]+T*(bDH[i][4]+T*bDH[i][5]))) +
				   bDH[i][6]*uT;
	}
}

// Entropia delle singole specie [J/mol.K]
void OpenSMOKE_IdealGas::SpeciesEntropy(double T)
{
	double logT=log(T);

	for(int i = 1;i <= NC;i++)
	{
		if(T > T2[i])	// High Temperature
			s[i] = aDS[i][1]*logT +
				   T*( aDS[i][2] + T*( aDS[i][3] + T *(aDS[i][4] + T*aDS[i][5])))+aDS[i][6];

		else			// Low Temperature
			s[i] = bDS[i][1]*logT +
				   T*( bDS[i][2] + T*( bDS[i][3] + T *(bDS[i][4] + T*bDS[i][5])))+bDS[i][6];
	}
}

void OpenSMOKE_IdealGas::SpeciesCp(double T)
{
	int i;

	for(i = 1;i <= NC;i++)
	{
		if(T > T2[i])
			Cp[i] = CpHT[i][1] + T*(CpHT[i][2] + T*(CpHT[i][3] + T*(CpHT[i][4] + T*CpHT[i][5])));

		else
			Cp[i] = CpLT[i][1] + T*(CpLT[i][2] + T*(CpLT[i][3] + T*(CpLT[i][4] + T*CpLT[i][5]))) ;
	}
}

void OpenSMOKE_IdealGas::SpeciesCv(void)
{
	for(int k = 1;k <= NC;k++)
		Cv[k] = Cp[k] - Constants::R_J_kmol*uM[k];
}

// ************************************************************************************	//
// ************************************************************************************ //
//								ALLOCAZIONI DI MEMORIA									//
// ************************************************************************************	//
// ************************************************************************************ //

void OpenSMOKE_IdealGas::Allocate(void)
{
	cout << "  Memory allocation..." << endl << endl;

	// Variabili termodinamiche
	aDH		=	new BzzVector [NC + 1];
	bDH		=   new BzzVector [NC + 1];
	aDS		=   new BzzVector [NC + 1];
	bDS		=   new BzzVector [NC + 1];

	CpHT	=	new BzzVector [NC + 1];
	CpLT	=	new BzzVector [NC + 1];

	for(int k=1;k<=NC;k++)
	{
		ChangeDimensions(6,&aDH[k]);
		ChangeDimensions(6,&bDH[k]);
		ChangeDimensions(6,&aDS[k]);
		ChangeDimensions(6,&bDS[k]);

		ChangeDimensions(5,&CpHT[k]);
		ChangeDimensions(5,&CpLT[k]);
	}

	// --------------------------------------------------------------
	// Molecular weigths
	// --------------------------------------------------------------
	ChangeDimensions(NC,&M);
	ChangeDimensions(NC,&uM);
	
	// --------------------------------------------------------------
	// Thermodynamics
	// --------------------------------------------------------------
	ChangeDimensions(NC,&Cp);
	ChangeDimensions(NC,&Cv);
	ChangeDimensions(NC,&h);
	ChangeDimensions(NC,&s);
	ChangeDimensions(NC,&T1);
	ChangeDimensions(NC,&T2);
	ChangeDimensions(NC,&T3);

	// --------------------------------------------------------------
	// Transport properties
	// --------------------------------------------------------------
	ChangeDimensions(NC,&lambda);
	ChangeDimensions(NC,&eta);
	ChangeDimensions(NC,&sumK);
	if (!iMixViscositySimplified)
	{
		ChangeDimensions(NC*(NC-1)/2,&PMRatio1su4);
		ChangeDimensions(NC*(NC-1)/2,&delta_phi);
		ChangeDimensions(NC*(NC-1)/2,&phi_eta_sup);
		ChangeDimensions(NC*(NC-1)/2,&phi_eta_inf);
		ChangeDimensions(NC,&sqrtEta);
		ChangeDimensions(NC,&usqrtEta);
	}
	else
	{
		ChangeDimensions(NC*(NC-1)/2,&sqrtPMRatio_inf);
		ChangeDimensions(NC*(NC-1)/2,&sqrtPMRatio_sup);
	}

	ChangeDimensions(NC, &x_aux);
	ChangeDimensions(NC,NC,&Djk);
	ChangeDimensions((NC*(NC-1))/2, &iD);
	int s=1;
	for (int j=1;j<=NC;j++)
		for (int k=j+1;k<=NC;k++)
			iD[s++]=k+(j-1)*NC;
}

void OpenSMOKE_IdealGas::Fitting(const string fileName, BzzLoad &binaryFile)
{
	double etaMin, etaMax;
	double lambdaMin, lambdaMax;
	double diffMin, diffMax;


	cout << "--------------------------------------------------------------------------" << endl;
	cout << "                     PROPERTIES FITTING                              " << endl;
	cout << "--------------------------------------------------------------------------" << endl;

	// Read From File
	binaryFile >> etaMin >> etaMax;
	binaryFile >> fittingEta;
	binaryFile >> lambdaMin >> lambdaMax;
	binaryFile >> fittingLambda;
	binaryFile >> diffMin >> diffMax;
	binaryFile >> fittingDbinary;
	cout << "--------------------------------------------------------------------------" << endl;
}

// ************************************************************************************	//
// ************************************************************************************ //
//										UTILITIES										//
// ************************************************************************************	//
// ************************************************************************************ //

double OpenSMOKE_IdealGas::GetMWFromMassFractions(BzzVector &y)
{
	return 1. / Dot(y, uM);
}

double OpenSMOKE_IdealGas::GetMWFromMoleFractions(BzzVector &x)
{
	return Dot(x, M);
}

void OpenSMOKE_IdealGas::GetMWAndMassFractionsFromMoleFractions(double &MWmix, BzzVector &y, BzzVector &x)
{
	ElementByElementProduct(x, M, &y);
	MWmix = y.GetSumElements();
	y /= MWmix;
}

void OpenSMOKE_IdealGas::GetMWAndMoleFractionsFromMassFractions(double &MWmix, BzzVector &x, BzzVector &y)
{
	ElementByElementProduct(y, uM, &x);
	MWmix = 1./x.GetSumElements();
	x *= MWmix;
}

void OpenSMOKE_IdealGas::GetMassFractionsFromMoleFractionsAndMW(BzzVector &y, BzzVector &x, const double MWmix)
{
	ElementByElementProduct(x, M, &y);
	y /= MWmix;
}

void OpenSMOKE_IdealGas::GetMoleFractionsFromMassFractionsAndMW(BzzVector &x, BzzVector &y, const double MWmix)
{
	ElementByElementProduct(y, uM, &x);
	x *= MWmix;
}


// ************************************************************************************	//
// ************************************************************************************ //
//										FITTING											//
// ************************************************************************************	//
// ************************************************************************************ //

void OpenSMOKE_IdealGas::SpeciesViscosityFromFitting(double T)
{
	// 1. Fitting polinomiale classico
	double logT=log(T);
	for (int j=1;j<=NC;j++)
	//	eta[j]=expVisc(fittingEta[j][1]+logT*(fittingEta[j][2]+logT*(fittingEta[j][3]+logT*fittingEta[j][4])));
		eta[j]=exp(fittingEta[j][1]+logT*(fittingEta[j][2]+logT*(fittingEta[j][3]+logT*fittingEta[j][4])));
}

void OpenSMOKE_IdealGas::SpeciesConductivityFromFitting(double T)
{
	double logT=log(T);
	for (int j=1;j<=NC;j++)
	//	lambda[j]=expCond(fittingLambda[j][1]+logT*(fittingLambda[j][2]+logT*(fittingLambda[j][3]+logT*fittingLambda[j][4])));
		lambda[j]=exp(fittingLambda[j][1]+logT*(fittingLambda[j][2]+logT*(fittingLambda[j][3]+logT*fittingLambda[j][4])));
}

void OpenSMOKE_IdealGas::SpeciesDiffusivityFromFitting(double T, double P)
{
	// Viene calcolata soltanto la meta' superiore di questa matrice, diagonale esclusa
	// si tratta infatti di una matrice simmetrica e la regola di miscela e' in grado
	// di sfruttare questa sola meta' superiore

	int *s=iD.GetHandle();

	double logT = log(T);
	for (int j=1;j<=NC;j++)
	for (int k=j+1;k<=NC;k++)
	{
		// Djk[j][k]=P/expDiff((fittingDbinary[*s][1]+logT*(fittingDbinary[*s][2]+logT*(fittingDbinary[*s][3]+logT*fittingDbinary[*s][4]))));
		Djk[j][k]=P/exp((fittingDbinary[*s][1]+logT*(fittingDbinary[*s][2]+logT*(fittingDbinary[*s][3]+logT*fittingDbinary[*s][4]))));
		s++;
	}
}

int OpenSMOKE_IdealGas::recognize_species(char* name)
{
	for(int i=1;i<=NC;i++)
		if (names[i] == name)
			return i;
	
	string dummy = name;
	ErrorMessage("This species is not included in the kinetic scheme: " + dummy);
	return -1;
}

int OpenSMOKE_IdealGas::recognize_species_without_exit(char* name)
{
	for(int i=1;i<=NC;i++)
		if (names[i] == name)
			return i;

	return 0;
}

int OpenSMOKE_IdealGas::recognize_species(string name)
{
	for(int i=1;i<=NC;i++)
		if (!name.compare(names[i]))
			return i;
	
	ErrorMessage("This species is not included in the kinetic scheme: " + name);
	return -1;
}


double GiveMeFunctionEnthalpy(BzzVector &A, BzzVector &B, const double T, const double Tlimit)
{
	if (T>=Tlimit)
		return T*(A[1]+T*(A[2]+T*(A[3]+T*(A[4]+T*A[5]))))+A[6];	// [J/kg]
	else
		return T*(B[1]+T*(B[2]+T*(B[3]+T*(B[4]+T*B[5]))))+B[6];	// [J/kg]
}

double OpenSMOKE_IdealGas::GetTemperatureFromMixEnthalpy(const double Hmix, const double MWmix, BzzVector &x)
{
	int i, k;
	const int MAX_ITERATIONS = 64;
	const double RELATIVE_ERROR = 1e-7;
	double TA, TB, TC;// TC_Old;
//	double fA, fB, fC;
	BzzVector A(6);
	BzzVector B(6);

	// Limit temperature
	double Tlimit = T2[1];

	// Coefficient for calculating the enthalpy [J, kg]
	for(k=1; k<=6; k++)	
		for(i=1; i<=NC; i++)	
			A[k] += aDH[i][k]*x[i];
	
	for(k=1; k<=6; k++)	
		for(i=1; i<=NC; i++)	
			B[k] += bDH[i][k]*x[i];

	A *= Constants::R_J_kmol / MWmix;
	B *= Constants::R_J_kmol / MWmix;

	// Find the temperature from the enthalpy
	TA		=  200.;
	TB		=  5000.;

	 // Newton iteration
	double dt;
	TC = 300.;
	for (int n = 0; n < 50; n++) 
	{
		cout << TC << endl;
		SpeciesCp(TC);
		double CpMix = MixCp_FromMoleFractions(x);

		double h0 = GiveMeFunctionEnthalpy(A,B, TC, Tlimit);
        dt = (Hmix - h0)/CpMix;

		// limit step size to 100 K
		if (dt > 100.0)			dt =  100.0;
		else if (dt < -100.0)	dt = -100.0; 
		TC += dt;
		if (fabs(dt) < 0.01)
			return TC;
	}

	ErrorMessage("Temperature from enthalpy: maximum number of iteration reached!");
	return -1;
}

int OpenSMOKE_IdealGas::NumberOfSpecies()
{
	return NC;
}

int	OpenSMOKE_IdealGas::NumberOfElements()
{
	return elements.Rows();
}

void OpenSMOKE_IdealGas::GetElementalMoleFractions(BzzVector &x, BzzVector &x_elemental)
{
	for(int k=1;k<=elements.Rows();k++)
	{
		BzzVector aux = elements.GetRow(k);
		x_elemental[k] = Dot(x, aux);
	}
	double sum = x_elemental.GetSumElements();
	x_elemental /= sum;
}

void OpenSMOKE_IdealGas::GetElementalMassFractions(BzzVector &omega, BzzVector &omega_elemental)
{
	for(int k=1;k<=elements.Rows();k++)
	{
		omega_elemental[k] = 0;
			for(int j=1;j<=NC;j++)
				omega_elemental[k] += omega[j]*elements[k][j]*m_elements[k]/M[j];
	}
	double sum = omega_elemental.GetSumElements();
	omega_elemental /= sum;
}

void OpenSMOKE_IdealGas::GetElementalMoleFractions(BzzMatrix &x, BzzMatrix &x_elemental)
{
	BzzVector x_vector;

	for(int j=1;j<=x.Rows();j++)
	{
		int k;
		x_vector = x.GetRow(j);
		for(k=1;k<=elements.Rows();k++)
		{
			BzzVector aux = elements.GetRow(k);
			x_elemental[j][k] = Dot(x_vector, aux);
		}

		double sum = x_elemental.GetRow(j).GetSumElements();
		for(k=1;k<=elements.Rows();k++)
			x_elemental[j][k] /= sum;
	}
}

void OpenSMOKE_IdealGas::GetElementalMassFractions(BzzMatrix &omega, BzzMatrix &omega_elemental)
{	
	for(int i=1;i<=omega.Rows();i++)
	{
		int k;
		for(k=1;k<=elements.Rows();k++)
		{
			omega_elemental[i][k] = 0;
				for(int j=1;j<=NC;j++)
					omega_elemental[i][k] += omega[i][j]*elements[k][j]*m_elements[k]/M[j];
		}
		double sum = omega_elemental.GetRow(i).GetSumElements();
		for(k=1;k<=elements.Rows();k++)
			omega_elemental[i][k] /= sum;
	}
}

double OpenSMOKE_IdealGas::GetMixtureFraction(BzzVector &omega, BzzVector &omegaFuel, BzzVector &omegaAir)
{
	int iH = recognize_element("h");
	int iO = recognize_element("o");
	int iC = recognize_element_without_error("c");

	double Z;

	if (iC>0)
	{
		Z = 2.*	(omega[iC] - omegaAir[iC])/m_elements[iC] + 
				(omega[iH] - omegaAir[iH])/(2.*m_elements[iH]) - 
				(omega[iO] - omegaAir[iO])/m_elements[iO];

		Z /= (	2.*(omegaFuel[iC] - omegaAir[iC])/m_elements[iC] + 
				(omegaFuel[iH] - omegaAir[iH])/(2.*m_elements[iH]) - 
				(omegaFuel[iO] - omegaAir[iO])/m_elements[iO]		);
	}
	else
	{
		Z =		(omega[iH] - omegaAir[iH])/(2.*m_elements[iH]) - 
				(omega[iO] - omegaAir[iO])/m_elements[iO];

		Z /= (	(omegaFuel[iH] - omegaAir[iH])/(2.*m_elements[iH]) - 
				(omegaFuel[iO] - omegaAir[iO])/m_elements[iO]		);
	}


	return Z;
}

int OpenSMOKE_IdealGas::recognize_element(const string name)
{
	for(int i=1;i<=NumberOfElements();i++)
		if (!name.compare(list_of_elements[i-1]))
			return i;

	ErrorMessage("This element is not included in the kinetic scheme: " + name);
	return -1;
}

int OpenSMOKE_IdealGas::recognize_element_without_error(const string name)
{
	for(int i=1;i<=NumberOfElements();i++)
		if (!name.compare(list_of_elements[i-1]))
			return i;
	return -1;
}

// *************************************************************************************** //
//								MIXTURE PROPERTIES (mass)                                  //
// *************************************************************************************** //

double OpenSMOKE_IdealGas::GetMixFormationEnthalpy_Mass(BzzVector &omega)
{
	return GetMixEnthalpy_Mass(Constants::T_Reference, omega);
}

double OpenSMOKE_IdealGas::GetMixFormationInternalEnergy_Mass(BzzVector &omega)
{
	return GetMixInternalEnergy_Mass(Constants::T_Reference, omega);
}

double OpenSMOKE_IdealGas::GetMixFormationEntropy_Mass(const double P_Pa, BzzVector &omega)
{
	return GetMixEntropy_Mass(Constants::T_Reference, P_Pa, omega);
}

double OpenSMOKE_IdealGas::GetMixFormationGibbsFreeEnergy_Mass(const double P_Pa, BzzVector &omega)
{
	return GetMixGibbsFreeEnergy_Mass(P_Pa, Constants::T_Reference, omega);
}

double OpenSMOKE_IdealGas::GetMixFormationHelmotzFreeEnergy_Mass(const double P_Pa, BzzVector &omega)
{
	return GetMixHelmotzFreeEnergy_Mass(P_Pa, Constants::T_Reference, omega);
}

// *************************************************************************************** //
//								MIXTURE PROPERTIES (mole)                                  //
// *************************************************************************************** //

double OpenSMOKE_IdealGas::GetMixFormationEnthalpy_Mole(BzzVector &x)
{
	return GetMixEnthalpy_Mole(Constants::T_Reference, x);
}

double OpenSMOKE_IdealGas::GetMixFormationInternalEnergy_Mole(BzzVector &x)
{
	return GetMixInternalEnergy_Mole(Constants::T_Reference, x);
}

double OpenSMOKE_IdealGas::GetMixFormationEntropy_Mole(const double P_Pa, BzzVector &x)
{
	return GetMixEntropy_Mole(Constants::T_Reference, P_Pa, x);
}

double OpenSMOKE_IdealGas::GetMixFormationGibbsFreeEnergy_Mole(const double P_Pa, BzzVector &x)
{
	return GetMixGibbsFreeEnergy_Mole(P_Pa, Constants::T_Reference, x);
}

double OpenSMOKE_IdealGas::GetMixFormationHelmotzFreeEnergy_Mole(const double P_Pa, BzzVector &x)
{
	return GetMixHelmotzFreeEnergy_Mole(P_Pa, Constants::T_Reference, x);
}

// *************************************************************************************** //
//								MIXTURE PROPERTIES (mass)                                  //
// *************************************************************************************** //

double OpenSMOKE_IdealGas::GetMixEnthalpy_Mass(const double T, BzzVector &omega)
{
	BzzVector H(NC);
	GetMixAveragedEnthalpy_Mass(H, T);	// [J/kg]
	return Dot(H, omega);				// [J/kg]
}

double OpenSMOKE_IdealGas::GetMixInternalEnergy_Mass(const double T, BzzVector &omega)
{
	BzzVector U(NC);
	GetMixAveragedInternalEnergy_Mass(U, T);	// [J/kg]
	return Dot(U, omega);						// [J/kg]
}

double OpenSMOKE_IdealGas::GetMixEntropy_Mass(const double P_Pa, const double T, BzzVector &omega)
{
	BzzVector S(NC);
	GetMixAveragedEntropy_Mass(S, P_Pa, T, omega);	// [J/kg/K]
	return Dot(S, omega);							// [J/kg/K]
}

double OpenSMOKE_IdealGas::GetMixGibbsFreeEnergy_Mass(const double P_Pa, const double T, BzzVector &omega)
{
	double Hmix = GetMixEnthalpy_Mass(T, omega);			// [J/kg]
	double Smix = GetMixEntropy_Mass(P_Pa, T, omega);		// [J/kg/K]
	
	return Hmix - T*Smix;									// [J/kg]
}

double OpenSMOKE_IdealGas::GetMixHelmotzFreeEnergy_Mass(const double P_Pa, const double T, BzzVector &omega)
{
	double Umix = GetMixInternalEnergy_Mass(T, omega);	// [J/kg]
	double Smix = GetMixEntropy_Mass(P_Pa, T, omega);	// [J/kg/K]

	return (Umix - T*Smix);								// [J/kg]
}


// *************************************************************************************** //
//								MIXTURE PROPERTIES (mole)                                  //
// *************************************************************************************** //

double OpenSMOKE_IdealGas::GetMixEnthalpy_Mole(const double T, BzzVector &x)
{
	BzzVector H(NC);
	GetMixAveragedEnthalpy_Mole(H, T);		// [J/kmol]
	return Dot(H, x);						// [J/kmol]
}

double OpenSMOKE_IdealGas::GetMixInternalEnergy_Mole(const double T, BzzVector &x)
{
	BzzVector U(NC);
	GetMixAveragedInternalEnergy_Mole(U, T);		// [J/kmol]
	return Dot(U, x);								// [J/kmol]
}

double OpenSMOKE_IdealGas::GetMixEntropy_Mole(const double P_Pa, const double T, BzzVector &x)
{
	BzzVector S(NC);
	GetMixAveragedEntropy_Mole(S, P_Pa, T, x);		// [J/kmol/K]
	return Dot(S, x);								// [J/kmol/K]
}

double OpenSMOKE_IdealGas::GetMixGibbsFreeEnergy_Mole(const double P_Pa, const double T, BzzVector &x)
{
	double Hmix = GetMixEnthalpy_Mole(T, x);		// [J/kmol]
	double Smix = GetMixEntropy_Mole(P_Pa, T, x);	// [J/kmol/K]
	
	return Hmix - T*Smix;							// [J/kmol]
}

double OpenSMOKE_IdealGas::GetMixHelmotzFreeEnergy_Mole(const double P_Pa, const double T, BzzVector &x)
{
	double Umix = GetMixInternalEnergy_Mole(T, x);	// [J/kmol]
	double Smix = GetMixEntropy_Mole(P_Pa, T, x);	// [J/kmol/K]

	return (Umix - T*Smix);							// [J/kmol]
}


// *************************************************************************************** //
//								STANDARD PROPERTIES (mass)                                 //
// *************************************************************************************** //

void OpenSMOKE_IdealGas::GetStandardEnthalpy_Mass(BzzVector &H, const double T)
{
	GetStandardEnthalpy_Mole(H, T);			// [J/kmol]
	ElementByElementProduct(H, uM, &H);		// [J/kg]
}

void OpenSMOKE_IdealGas::GetStandardEntropy_Mass(BzzVector &S, const double T)
{
	GetStandardEntropy_Mole(S, T);			// [J/kmol/K]
	ElementByElementProduct(S, uM, &S);		// [J/kg/K]
}

void OpenSMOKE_IdealGas::GetStandardInternalEnergy_Mass(BzzVector &U, const double T)
{
	GetStandardInternalEnergy_Mole(U, T);	// [J/kmol]
	ElementByElementProduct(U, uM, &U);		// [J/kg]
}

void OpenSMOKE_IdealGas::GetStandardGibbsFreeEnergy_Mass(BzzVector &G, const double T)
{
	GetStandardGibbsFreeEnergy_Mole(G, T);	// [J/kmol]
	ElementByElementProduct(G, uM, &G);		// [J/kg]
}

void OpenSMOKE_IdealGas::GetStandardHelmotzFreeEnergy_Mass(BzzVector &A, const double T)
{
	GetStandardHelmotzFreeEnergy_Mole(A, T);// [J/kmol]
	ElementByElementProduct(A, uM, &A);		// [J/kg]
}


// *************************************************************************************** //
//								STANDARD PROPERTIES (mole)                                 //
// *************************************************************************************** //

void OpenSMOKE_IdealGas::GetStandardEnthalpy_Mole(BzzVector &H, const double T)
{
	ChangeDimensions(NC, &H);
	for(int i = 1;i <= NC;i++)
	{
		if(T > T2[i])	// High Temperature
			H[i] = aDH[i][1] + T*(aDH[i][2]+T*(aDH[i][3]+T*(aDH[i][4]+T*aDH[i][5]))) + aDH[i][6]/T;
		else
			H[i] = bDH[i][1] + T*(bDH[i][2]+T*(bDH[i][3]+T*(bDH[i][4]+T*bDH[i][5]))) + bDH[i][6]/T;
	}

	H *= Constants::R_J_kmol*T;			// [J/kmol]
}

void OpenSMOKE_IdealGas::GetStandardEntropy_Mole(BzzVector &S, const double T)
{
	double logT = log(T);

	ChangeDimensions(NC, &S);
	for(int i = 1;i <= NC;i++)
	{
		if(T > T2[i])	// High Temperature
			S[i] = aDS[i][1]*logT +
				   T*( aDS[i][2] + T*( aDS[i][3] + T*(aDS[i][4] + T*aDS[i][5])))+aDS[i][6];
		else
			S[i] = bDS[i][1]*logT +
				   T*( bDS[i][2] + T*( bDS[i][3] + T*(bDS[i][4] + T*bDS[i][5])))+bDS[i][6];
	}

	S *= Constants::R_J_kmol;			// [J/kmol/K]
}

void OpenSMOKE_IdealGas::GetStandardInternalEnergy_Mole(BzzVector &U, const double T)
{
	BzzVector H(NC);
	
	GetStandardEnthalpy_Mole(H,T);		// [J/kmol]
	U  = H;								// [J/kmol]
	
	double additional = Constants::R_J_kmol*T;		// [J/kmol]
	for(int i=1;i<=NC;i++)	U[i] -= additional;		// [J/kmol]
}

void OpenSMOKE_IdealGas::GetStandardGibbsFreeEnergy_Mole(BzzVector &G, const double T)
{
	BzzVector H(NC);
	BzzVector S(NC);
	
	GetStandardEnthalpy_Mole(H,T);	// [J/kmol]
	GetStandardEntropy_Mole(S,T);	// [J/kmol/K]

	G  = H;
	G -= T*S;						// [J/kmol/K]
}

void OpenSMOKE_IdealGas::GetStandardHelmotzFreeEnergy_Mole(BzzVector &A, const double T)
{
	BzzVector U(NC);
	BzzVector S(NC);

	GetStandardInternalEnergy_Mole(U,T);	// [J/kmol]
	GetStandardEntropy_Mole(S,T);			// [J/kmol/K]

	A  = U;									// [J/kmol/K]
	A -= T*S;								// [J/kmol/K]
}


// *************************************************************************************** //
//								MIX-AVERAGED PROPERTIES (mass)                             //
// *************************************************************************************** //

void OpenSMOKE_IdealGas::GetMixAveragedEnthalpy_Mass(BzzVector &H, const double T)
{
	GetStandardEnthalpy_Mass(H, T);	// [J/kg/K]
}

void OpenSMOKE_IdealGas::GetMixAveragedEntropy_Mass(BzzVector &S, const double P_Pa, const double T, BzzVector &omega)
{
	double MWmix;
	BzzVector x(NC);
	GetMWAndMoleFractionsFromMassFractions(MWmix, x, omega);
	GetMixAveragedEntropy_Mole(S, P_Pa, T, x);					// [J/kmol/K]
	ElementByElementProduct(S, uM, &S);							// [J/kg/K]
}

void OpenSMOKE_IdealGas::GetMixAveragedInternalEnergy_Mass(BzzVector &U, const double T)
{
	GetStandardInternalEnergy_Mass(U, T);	// [J/kg/K]
}

void OpenSMOKE_IdealGas::GetMixAveragedGibbsFreeEnergy_Mass(BzzVector &G, const double P_Pa, const double T, BzzVector &omega)
{
	double MWmix;
	BzzVector x(NC);
	GetMWAndMoleFractionsFromMassFractions(MWmix, x, omega);
	GetMixAveragedGibbsFreeEnergy_Mole(G, P_Pa, T, x);			// [J/kmol]
	ElementByElementProduct(G, uM, &G);							// [J/kg]
}

void OpenSMOKE_IdealGas::GetMixAveragedHelmotzFreeEnergy_Mass(BzzVector &A, const double P_Pa, const double T, BzzVector &omega)
{
	double MWmix;
	BzzVector x(NC);
	GetMWAndMoleFractionsFromMassFractions(MWmix, x, omega);
	GetMixAveragedHelmotzFreeEnergy_Mole(A, P_Pa, T, x);		// [J/kmol]
	ElementByElementProduct(A, uM, &A);							// [J/kg]
}

// *************************************************************************************** //
//								MIX-AVERAGED PROPERTIES (mole)                             //
// *************************************************************************************** //

void OpenSMOKE_IdealGas::GetMixAveragedEnthalpy_Mole(BzzVector &H, const double T)
{
	GetStandardEnthalpy_Mole(H, T);	// [J/kmol]
}

void OpenSMOKE_IdealGas::GetMixAveragedEntropy_Mole(BzzVector &S, const double P_Pa, const double T, BzzVector &x)
{
	BzzVector add(NC);

	double additional = log(P_Pa/Constants::P_Reference);

	for(int i = 1;i <= NC;i++)
		if (x[i] > Constants::SmallMoleFraction)	add[i] = log(x[i]) + additional;	// [-]
		else										add[i] = additional;				// [-]
	add *= Constants::R_J_kmol;															// [J/kmol/K]

	GetStandardEntropy_Mole(S, T);	// [J/kmol/K]
	S -= add;						// [J/kmol/K]
}

void OpenSMOKE_IdealGas::GetMixAveragedInternalEnergy_Mole(BzzVector &U, const double T)
{
	GetStandardInternalEnergy_Mole(U, T);		// [J/kmol]
}

void OpenSMOKE_IdealGas::GetMixAveragedGibbsFreeEnergy_Mole(BzzVector &G, const double P_Pa, const double T, BzzVector &x)
{
	BzzVector H(NC);
	BzzVector S(NC);
	
	GetMixAveragedEnthalpy_Mole(H, T);			// [J/kmol]
	GetMixAveragedEntropy_Mole(S, P_Pa, T, x);	// [J/kmol/K]

	G  = H;										// [J/kmol]
	G -= T*S;									// [J/kmol]
}

void OpenSMOKE_IdealGas::GetMixAveragedHelmotzFreeEnergy_Mole(BzzVector &A, const double P_Pa, const double T, BzzVector &x)
{
	BzzVector U(NC);
	BzzVector S(NC);

	GetMixAveragedInternalEnergy_Mole(U,T);		// [J/kmol]
	GetMixAveragedEntropy_Mole(S, P_Pa, T, x);	// [J/kmol/K]

	A  = U;										// [J/kmol]
	A -= T*S;									// [J/kmol]
}

double OpenSMOKE_IdealGas::MixCp_FromMassFractions(BzzVector &y)
{
	double CpMix = Dot(Cp, y);
	return CpMix;
}

double OpenSMOKE_IdealGas::MixCp_FromMoleFractions(BzzVector &x)
{
	double MWmix;
	BzzVector y(NC);
	GetMWAndMassFractionsFromMoleFractions(MWmix, y, x);
	return Dot(Cp, y);
}

double OpenSMOKE_IdealGas::MixCv_FromCpMix(const double CpMix, const double MWmix)
{
	return CpMix - Constants::R_J_kmol/MWmix;
}
