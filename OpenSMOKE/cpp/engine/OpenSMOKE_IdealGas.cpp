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
#include "basic/OpenSMOKE_Conversions.h"
#include "engine/OpenSMOKE_IdealGas.h"
#include "addons/OpenSMOKE_PolimiSoot.h"

#include <iomanip>

const int			OpenSMOKE_IdealGas::N_LOOKUP	= 4000;								// [-]
const double		OpenSMOKE_IdealGas::logCATM		= log(101325./Constants::R_J_kmol);	// [kmol.K/m3]

void OpenSMOKE_IdealGas::SetASCII()
{
	binary_version_ = false;
}

OpenSMOKE_IdealGas::OpenSMOKE_IdealGas()
{
	TMIN		= 300.;
	TMAX		= 3600.; 

	name_object   = "[not assigned]";
	name_author   = "[not assigned]";
	building_date = GiveMeTimeAndDate();

	iMixViscositySimplified = 0;
	iTransportMode			= true;
	binary_version_			= true;
}

void OpenSMOKE_IdealGas::ErrorMessage(const std::string message)
{
    cout << endl;
    cout << "Class:  OpenSMOKE_IdealGas"	<< endl;
    cout << "Object: " << name_object		<< endl;
    cout << "Error:  " << message           << endl;
    cout << "Press a key to continue... "   << endl;
    getchar();
    exit(-1);
}

void OpenSMOKE_IdealGas::WarningMessage(const std::string message)
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

void OpenSMOKE_IdealGas::SetName(const std::string name)
{
	name_object = name;
}

void OpenSMOKE_IdealGas::SetAuthorName(const std::string name)
{
	name_author = name;
}


void OpenSMOKE_IdealGas::ReadVersion(BzzLoad &binaryFile)
{
	std::string version;
	
	cout << "Reading version (1): " << binary_version_ << endl;
	if (binary_version_ == true)	
	{
		cout << "Reading binary version " << Constants::NAME_SIZE << endl;
		char dummy[Constants::NAME_SIZE];
		binaryFile.fileLoad.read((char*) dummy, sizeof(dummy));
		version = dummy;

	}
	else
	{
		cout << "Reading ascii version " << Constants::NAME_SIZE << endl;
		binaryFile >> version;
		cout << "Reading ascii version " << version << endl;
	}
	
	if		(version == "V090524")	iVersion = V090524;	
	else if (version == "V090905")	iVersion = V090905;
	else if (version == "V101116")	iVersion = V101116;
	else	ErrorMessage("This kinetic file version is not supported: " + version);
	binaryFile.End();
	cout << "OK" <<endl;
}

bool OpenSMOKE_IdealGas::TransportMode() 
{
	return iTransportMode;
}

void OpenSMOKE_IdealGas::Setup(BzzLoad &binaryFile)
{
	int k;
	
	cout << "--------------------------------------------------------------------------" << endl;
	cout << "                 THERMODYNAMICS AND TRANSPORT PROPERTIES                  " << endl;
	cout << "--------------------------------------------------------------------------" << endl;

	std::string version;
	cout << "Reading version (2): " << binary_version_ << endl;
	if (binary_version_ == true)	
	{
		char dummy[Constants::NAME_SIZE];
		binaryFile.fileLoad.read((char*) dummy, sizeof(dummy));
		version = dummy;

		if		(version == "V090524")	iVersion = V090524;	
		else if (version == "V090905")	iVersion = V090905;
		else if (version == "V101116")	iVersion = V101116;
		else	ErrorMessage("This kinetic file version is not supported: " + version);
		cout << "Version: " << version << endl;

		// Reading information
		binaryFile.fileLoad.read((char*) dummy, sizeof(dummy));
		cout << "Name:    "  << dummy << endl;
		binaryFile.fileLoad.read((char*) dummy, sizeof(dummy));
		cout << "Author:  "  << dummy << endl;
		binaryFile.fileLoad.read((char*) dummy, sizeof(dummy));
		//cout << "Place:   "  << dummy << endl;
		binaryFile.fileLoad.read((char*) dummy, sizeof(dummy));
		cout << "Date:    "  << dummy << endl;
	}
	else
	{
		binaryFile >> version;
		
		if		(version == "V090524")	iVersion = V090524;	
		else if (version == "V090905")	iVersion = V090905;
		else if (version == "V101116")	iVersion = V101116;
		else	ErrorMessage("This kinetic file version is not supported: " + version);
		cout << "Version: " << version << endl;

		// Reading information
		std::string dummy;
		binaryFile >> dummy;
		cout << "Name:    "  << dummy << endl;
		binaryFile >> dummy;
		cout << "Author:  "  << dummy << endl;
		binaryFile >> dummy;
		//cout << "Place:   "  << dummy << endl;
		binaryFile >> dummy;
		cout << "Date:    "  << dummy << endl;
	}

	// Minimum and maximum temperature
	binaryFile >> TMIN;			cout << "Minimum temperature: " << TMIN << " K" << endl;
	binaryFile >> TMAX;			cout << "Maximum temperature: " << TMAX << " K" << endl;
	SetMinimumTemperature(TMIN);
	SetMaximumTemperature(TMAX);

	// Reading number of Species
	binaryFile >> NC;

	// Lettura dei nomi delle specie e del numero di componenti
	if (binary_version_ == true)	
	{
		char name[Constants::NAME_SIZE];
		names = new std::string[NC + 1];
		for(k=1;k<=NC;k++)
		{
			binaryFile.fileLoad.read((char*) name, sizeof(name));
			names[k] = name;
		}
	}
	else
	{
		names = new std::string[NC + 1];
		for(k=1;k<=NC;k++)
		{
			binaryFile >> names[k];
		}
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
	if (binary_version_ == true)	
	{
		for(k=1;k<=elements.Rows();k++)
		{
			char name[Constants::NAME_SIZE];
			binaryFile.fileLoad.read((char*) name, sizeof(name));
			list_of_elements.push_back(name);
		}
	}
	else
	{
		for(k=1;k<=elements.Rows();k++)
		{
			std::string name;
			binaryFile >> name;
			list_of_elements.push_back(name);
		}
	}

	// Post Processing information
	for(k=1;k<=NC;k++)
		uM[k] = 1./M[k];

	ChangeDimensions(m_elements.Size(),&um_elements);
	for(k=1;k<=m_elements.Size();k++)
		um_elements[k] = 1./m_elements[k];

	PostProcessMeanTemperatureForThermodynamicProperties();

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
//    pah_manager.recognizeSpecies(NC,  names, M);
}

void OpenSMOKE_IdealGas::SetPolimiSoot(	const unsigned int bin_index_zero, const double bin_density_A, 
					const unsigned int bin_index_final, const double bin_density_B,
					const double Df, const std::string bin_minimum_soot, const std::string bin_minimum_aggregates)
{
	polimiSoot = new OpenSMOKE_PolimiSoot;
	polimiSoot->Setup(*this, bin_minimum_soot, bin_minimum_aggregates, bin_index_zero, bin_density_A, bin_index_final, bin_density_B, Df);
}

void OpenSMOKE_IdealGas::SetPolimiSoot(	)
{
	polimiSoot = new OpenSMOKE_PolimiSoot;
	polimiSoot->Setup(*this, "BIN5", "BIN12");
}

void OpenSMOKE_IdealGas::PostProcessMeanTemperatureForThermodynamicProperties()
{
	ChangeDimensions(0, &Tintervals);
	BzzVector aux=T2;
	Sort(&aux);
	Tintervals.Append(aux[1]);
	for(int i=2;i<=NC;i++)
		if (aux[i] != Tintervals[Tintervals.Size()])
				Tintervals.Append(aux[i]);
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

	double massaMolareMedia=Dot(x_aux,M);
//	for(i = 1;i <= NC;i++)
//		massaMolareMedia+=x_aux[i]*M[i];

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

//	x = x_aux;
}

void OpenSMOKE_IdealGas::MixThermalDiffusionRatios(BzzVector &Teta, BzzVector &x)
{
	Teta = 0.;
	for (int k=1;k<=NC;k++)
		for (int j=1;j<=NC;j++)
			Teta[k] += Tetakj[k][j]*x[k]*x[j];
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

// Entalpia delle singole specie [J/kmol]
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

// Entropia delle singole specie [J/kmol/K]
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

void OpenSMOKE_IdealGas::SpeciesCp_ForcingCorrelation(BzzVector &T, const char low_or_high, BzzVector &cpi, BzzVector &d1cpi, BzzVector &d2cpi, BzzVector &d3cpi, BzzVector &d4cpi)
{
	ChangeDimensions(NC, &cpi);
	ChangeDimensions(NC, &d1cpi);
	ChangeDimensions(NC, &d2cpi);
	ChangeDimensions(NC, &d3cpi);
	ChangeDimensions(NC, &d4cpi);

	if (low_or_high == 'H')
	{
		for(int i = 1;i <= NC;i++)
		{
			cpi[i] = CpHT[i][1] + T[i]*(CpHT[i][2] + T[i]*(CpHT[i][3] + T[i]*(CpHT[i][4] + T[i]*CpHT[i][5])));
			d1cpi[i] = CpHT[i][2] + T[i]*(2.*CpHT[i][3] + T[i]*(3.*CpHT[i][4] + T[i]*4.*CpHT[i][5]));
			d2cpi[i] = 2.*CpHT[i][3] + T[i]*(6.*CpHT[i][4] + 12.*T[i]*CpHT[i][5]);
			d3cpi[i] = 6.*CpHT[i][4] + 24.*T[i]*CpHT[i][5];
			d4cpi[i] = 24.*CpHT[i][5];
		}
	}
	else if (low_or_high == 'L')
	{
		for(int i = 1;i <= NC;i++)
		{
			cpi[i] = CpLT[i][1] + T[i]*(CpLT[i][2] + T[i]*(CpLT[i][3] + T[i]*(CpLT[i][4] + T[i]*CpLT[i][5])));
			d1cpi[i] = CpLT[i][2] + T[i]*(2.*CpLT[i][3] + T[i]*(3.*CpLT[i][4] + T[i]*4.*CpLT[i][5]));
			d2cpi[i] = 2.*CpLT[i][3] + T[i]*(6.*CpLT[i][4] + 12.*T[i]*CpLT[i][5]);
			d3cpi[i] = 6.*CpLT[i][4] + 24.*T[i]*CpLT[i][5];
			d4cpi[i] = 24.*CpLT[i][5];
		}
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
	ChangeDimensions(NC, NC, &Tetakj);
}

void OpenSMOKE_IdealGas::Fitting(BzzLoad &binaryFile)
{
	double etaMin, etaMax;
	double lambdaMin, lambdaMax;
	double diffMin, diffMax;
	double tetaMin, tetaMax;

	if (iVersion == V101116)
	{
		int dummy_int;
		binaryFile >> dummy_int;
		if (dummy_int == 1)	iTransportMode = true;
		if (dummy_int == 0)	iTransportMode = false;
	}

	if (iTransportMode == true)
	{
		cout << " Reading fitting coefficients for transport properties" << endl;

		// Read From File
		binaryFile >> etaMin >> etaMax;
		binaryFile >> fittingEta;
		binaryFile >> lambdaMin >> lambdaMax;
		binaryFile >> fittingLambda;
		binaryFile >> diffMin >> diffMax;
		binaryFile >> fittingDbinary;
		
		if (iVersion == V090905 || iVersion == V101116)
		{
			binaryFile >> tetaMin >> tetaMax;
			binaryFile >> iThermalDiffusionRatios;
			binaryFile >> fittingTetaBinary;
		}
	}
	
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

void OpenSMOKE_IdealGas::GetMassFractionsFromMoleFractions(BzzVector &y, BzzVector &x)
{
	ElementByElementProduct(x, M, &y);
	const double MWmix = y.GetSumElements();
	y /= MWmix;
}

void OpenSMOKE_IdealGas::GetMoleFractionsFromMassFractions(BzzVector &x, BzzVector &y)
{
	ElementByElementProduct(y, uM, &x);
	const double MWmix = 1. / x.GetSumElements();
	x *= MWmix;
}

void OpenSMOKE_IdealGas::GetMWAndMoleFractionsFromMassFractions(BzzVector &MWmix, BzzMatrix &x, BzzMatrix &y)
{
	for(int i=1;i<=x.Rows();i++)
	{
		MWmix[i] = 0.;
		for(int j=1;j<=x.Columns();j++)
		{
			x[i][j]   = y[i][j]*uM[j];
			MWmix[i] += x[i][j];
		}
		MWmix[i] = 1./MWmix[i];

		for(int j=1;j<=x.Columns();j++)
			x[i][j] *= MWmix[i];
	}
}

void OpenSMOKE_IdealGas::GetMWAndMassFractionsFromMoleFractions(BzzVector &MWmix, BzzMatrix &y, BzzMatrix &x)
{
	for(int i=1;i<=x.Rows();i++)
	{
		MWmix[i] = 0.;
		for(int j=1;j<=x.Columns();j++)
		{
			y[i][j]   = x[i][j]*M[j];
			MWmix[i] += y[i][j];
		}

		for(int j=1;j<=x.Columns();j++)
			y[i][j] /= MWmix[i];
	}
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

void OpenSMOKE_IdealGas::CorrectSpeciesDiffusivityForSoot()
{
	/*
	unsigned int jBIN1A = 0;
	double MWBIN1A = 0.;

	for (int j = 1; j <= NC; j++)
	{
		if (names[j] == "BIN1A")
		{
			jBIN1A = j;
			MWBIN1A = M[j];
			break;
		}
	}

	if (jBIN1A > 0)
	{
		for (int j = 1; j <= NC; j++)
		{
			if (names[j].compare(0, 3, "BIN") == 0)
			{
				const double MWratio = M[j] / MWBIN1A;
				const double correctionCoefficient = pow(MWratio, -0.681);

				for (int k = j + 1; k <= NC; k++)
					Djk[j][k] = Djk[jBIN1A][k] * correctionCoefficient;
			}
			else
			{
				for (int k = j + 1; k <= NC; k++)
				{
					if (names[k].compare(0, 3, "BIN") == 0)
					{
						const double MWratio = M[k] / MWBIN1A;
						const double correctionCoefficient = pow(MWratio, -0.681);

						Djk[j][k] = Djk[j][jBIN1A] * correctionCoefficient;
					}
				}
			}
		}

	}
	*/
}

void OpenSMOKE_IdealGas::SpeciesThermalDiffusionRatiosFromFitting(double T)
{
	Tetakj = 0.;
	for (int i=1;i<=iThermalDiffusionRatios.Size();i++)
		for (int k=1;k<=NC;k++)
//		{	
//			int j = k + (i-1)*NC;
//			Tetakj[k][j] += fittingTetaBinary[j][1]+T*(fittingTetaBinary[j][2]+T*(fittingTetaBinary[j][3]+T*fittingTetaBinary[j][4]));
//		}
		{	
			int j = k + (i-1)*NC;
			Tetakj[iThermalDiffusionRatios[i]][k] += fittingTetaBinary[j][1]+T*(fittingTetaBinary[j][2]+T*(fittingTetaBinary[j][3]+T*fittingTetaBinary[j][4]));
		}
}

void OpenSMOKE_IdealGas::CorrectBinaryDiffusionCoefficients(const int j, const int k, const double correction)
{
	cout << "Indices: " << k << " " << j << endl;
	cout << "Indices: " << k+(j-1)*NC << " " << j+(k-1)*NC << endl;
	cout << "Before:  " << fittingDbinary[k+(j-1)*NC][1] << " " << fittingDbinary[j+(k-1)*NC][1] << endl;
	fittingDbinary[k+(j-1)*NC][1] += log(correction);
	fittingDbinary[j+(k-1)*NC][1] += log(correction);
	cout << "After:   " << fittingDbinary[k+(j-1)*NC][1] << " " << fittingDbinary[j+(k-1)*NC][1] << endl;
}

void OpenSMOKE_IdealGas::CorrectFormationEnthalpy(const int j, const double correction)
{
	cout << "Correct formation enthalpy " << names[j] << endl;
	cout << "Before: " << aDH[j][6] << " " << bDH[j][6] << endl;

	double c=correction/Constants::R_J_kmol;
	aDH[j][6] += c;
	bDH[j][6] += c;

	cout << "After: " << aDH[j][6] << " " << bDH[j][6] << endl;
}

int OpenSMOKE_IdealGas::recognize_species(char* name)
{
	for(int i=1;i<=NC;i++)
		if (names[i] == name)
			return i;
	
	std::string dummy = name;
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

int OpenSMOKE_IdealGas::recognize_species(const std::string name)
{
	for(int i=1;i<=NC;i++)
		if (!name.compare(names[i]))
			return i;
	
	ErrorMessage("This species is not included in the kinetic scheme: " + name);
	return -1;
}

int OpenSMOKE_IdealGas::recognize_species_without_exit(const std::string name)
{
	for(int i=1;i<=NC;i++)
		if (!name.compare(names[i]))
			return i;
	
	return 0;
}


double GiveMeFunctionEnthalpy(BzzVector &A, BzzVector &B, const double T, const double Tlimit)
{
	if (T>=Tlimit)
		return T*(A[1]+T*(A[2]+T*(A[3]+T*(A[4]+T*A[5]))))+A[6];	// [J/kg]
	else
		return T*(B[1]+T*(B[2]+T*(B[3]+T*(B[4]+T*B[5]))))+B[6];	// [J/kg]
}

int OpenSMOKE_IdealGas::NumberOfSpecies()
{
	return NC;
}

int	OpenSMOKE_IdealGas::NumberOfElements()
{
	return elements.Rows();
}

void OpenSMOKE_IdealGas::GetElementalMassFractionsFromElementalMoleFractionsAndMW(BzzVector &omega_elemental, BzzVector &x_elemental, const double MWmix)
{
	for(int k=1;k<=elements.Rows();k++)
		omega_elemental[k] = x_elemental[k]*m_elements[k]/MWmix;
}

void OpenSMOKE_IdealGas::GetElementalMoleFractionsFromElementalMassFractionsAndMW(BzzVector &x_elemental, BzzVector &omega_elemental, const double MWmix)
{
	for(int k=1;k<=elements.Rows();k++)
		x_elemental[k] = omega_elemental[k]*um_elements[k]*MWmix;
}

void OpenSMOKE_IdealGas::GetMWAndElementalMassFractionsFromElementalMoleFractions(double &MWmix, BzzVector &omega_elemental, BzzVector &x_elemental)
{
	MWmix = Dot(x_elemental, m_elements);
	for(int k=1;k<=elements.Rows();k++)
		omega_elemental[k] = x_elemental[k]*m_elements[k]/MWmix;
}

void OpenSMOKE_IdealGas::GetMWAndElementalMoleFractionsFromElementalMassFractions(double &MWmix, BzzVector &x_elemental, BzzVector &omega_elemental)
{
	MWmix = 1./Dot(omega_elemental, um_elements);
	for(int k=1;k<=elements.Rows();k++)
		x_elemental[k] = omega_elemental[k]*um_elements[k]*MWmix;
}


void OpenSMOKE_IdealGas::GetElementalMoleFractionsFromSpeciesMoleFractions(BzzVector &x_elemental, BzzVector &x)
{
	for(int k=1;k<=elements.Rows();k++)
	{
		BzzVector aux = elements.GetRow(k);
		x_elemental[k] = Dot(x, aux);
	}
	double sum = x_elemental.GetSumElements();
	x_elemental /= sum;
}

void OpenSMOKE_IdealGas::GetElementalMassFractionsFromSpeciesMassFractions(BzzVector &omega_elemental, BzzVector &omega)
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

void OpenSMOKE_IdealGas::GetElementalMassFractionsFromSpeciesMoleFractions(BzzVector &omega_elemental, BzzVector &x)
{
	double mw;
	BzzVector omega(x.Size());
	GetMWAndMassFractionsFromMoleFractions(mw, omega, x);
	GetElementalMassFractionsFromSpeciesMassFractions(omega_elemental, omega);
}


void OpenSMOKE_IdealGas::GetElementalMoleFractionsFromSpeciesMoleFractions(BzzMatrix &x_elemental, BzzMatrix &x)
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

void OpenSMOKE_IdealGas::GetElementalMassFractionsFromSpeciesMassFractions(BzzMatrix &omega_elemental, BzzMatrix &omega)
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
		if (iH>0)
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
			Z = 2.*	(omega[iC] - omegaAir[iC])/m_elements[iC] - 
					(omega[iO] - omegaAir[iO])/m_elements[iO];

			Z /= (	2.*(omegaFuel[iC] - omegaAir[iC])/m_elements[iC] - 
					(omegaFuel[iO] - omegaAir[iO])/m_elements[iO]		);
		}
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

int OpenSMOKE_IdealGas::recognize_element(const std::string name)
{
	for(int i=1;i<=NumberOfElements();i++)
		if (caseInsCompare(name, list_of_elements[i-1]) == true)
			return i;

	ErrorMessage("This element is not included in the kinetic scheme: " + name);
	return -1;
}

int OpenSMOKE_IdealGas::recognize_element_without_error(const std::string name)
{
	for(int i=1;i<=NumberOfElements();i++)
		if (caseInsCompare(name, list_of_elements[i-1]))
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
	GetStandardEnthalpy_Mass(H, T);	// [J/kg]
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

void OpenSMOKE_IdealGas::VerboseDataSpecies(ofstream &fOutput)
{
	int nIntervals = int((TMAX-TMIN)/100.);

	BzzVector T_Vector(nIntervals);
	BzzMatrix Cp_Matrix(nIntervals, NC);
	BzzMatrix H_Matrix(nIntervals, NC);
	BzzMatrix G_Matrix(nIntervals, NC);
	BzzMatrix S_Matrix(nIntervals,NC);

	BzzVector H_Row(NC);
	BzzVector S_Row(NC);
	BzzVector G_Row(NC);

	BzzVector H_Tref(NC);
	BzzVector G_Tref(NC);

	// Reference values
	{
		GetStandardEnthalpy_Mole(H_Tref, Constants::T_Reference);			// [J/kmol]
		GetStandardGibbsFreeEnergy_Mole(G_Tref, Constants::T_Reference);	// [J/kmol]
	}

	// Complete table
	for(int i=1;i<=nIntervals;i++)
	{
		T_Vector[i] = TMIN+(i-1)*100.;

		SpeciesCp(T_Vector[i]);
		Cp_Matrix.SetRow(i, Cp);

		GetStandardEnthalpy_Mole(H_Row, T_Vector[i]);			// [J/kmol]
		GetStandardGibbsFreeEnergy_Mole(G_Row, T_Vector[i]);	// [J/kmol]
		GetStandardEntropy_Mole(S_Row, T_Vector[i]);			// [J/kmol/K]

		H_Matrix.SetRow(i, H_Row);
		G_Matrix.SetRow(i, G_Row);
		S_Matrix.SetRow(i, S_Row);
	}

	for(int k=1;k<=NC;k++)
		VerboseDataSpecies(fOutput, k, T_Vector, Cp_Matrix, H_Matrix, G_Matrix, S_Matrix, H_Tref, G_Tref);
}

void OpenSMOKE_IdealGas::VerboseDataSpecies(ofstream &fOutput, const int k, BzzVector &T_Vector, BzzMatrix &Cp_Matrix,	BzzMatrix &H_Matrix, BzzMatrix &G_Matrix, BzzMatrix &S_Matrix,
																														BzzVector &H_Tref,	 BzzVector &G_Tref)
{
	int i;
	double conversion_h_and_g = OpenSMOKE_Conversions::J_from_kcal*1000.;
	double conversion_s       = OpenSMOKE_Conversions::J_from_kcal;
	double conversion_cp      = M[k]/OpenSMOKE_Conversions::J_from_kcal;

	fOutput << " ====================================================================================================" << endl;
	fOutput << "   THERMO TABLE FOR MOLECULE " << names[k] << " IN PHASE GAS                                         " << endl;
	fOutput << " ====================================================================================================" << endl;
	fOutput << "   Elemental Composition:" << endl;

	for(i=1;i<=elements.Rows();i++)
		if (elements[i][k] != 0)	fOutput << "     " << list_of_elements[i-1] << " " << elements[i][k] << endl;

	fOutput << "   Formation Enthalpy at 298K =          " << H_Tref[k]/conversion_h_and_g	<< " kcal/mol"	<< endl;
	fOutput << "   Formation Free Gibbs Energy at 298K = " << G_Tref[k]/conversion_h_and_g	<< " kcal/mol"	<< endl;
	fOutput << "   Molecular Weight =                    " << M[k]							<< " kg/kmol"	<< endl;
	fOutput << " -----------------------------------------------------------------------------------------------------" << endl;
	fOutput << "   Temperature       Cp             H             G            S              DH         DG                           " << endl;          
    fOutput << "       [K]       [cal/mol/K]    [kcal/mol]    [kcal/mol]  [cal/(mol-K)]   [kcal/mol]  [kcal/mol]                          " << endl;        
	fOutput << " -----------------------------------------------------------------------------------------------------" << endl;

	for(i=1;i<=T_Vector.Size();i++)
	{
		fOutput << "     " << setw(14) << left << setprecision(4) << fixed << T_Vector[i];
		fOutput            << setw(14) << left << setprecision(4) << fixed << Cp_Matrix[i][k]*conversion_cp;
		fOutput            << setw(14) << left << setprecision(4) << fixed << H_Matrix[i][k]/conversion_h_and_g;
		fOutput            << setw(14) << left << setprecision(4) << fixed << G_Matrix[i][k]/conversion_h_and_g;
		fOutput            << setw(14) << left << setprecision(4) << fixed << S_Matrix[i][k]/conversion_s;
		fOutput            << setw(14) << left << setprecision(4) << fixed << (H_Matrix[i][k]-H_Tref[k])/conversion_h_and_g;
		fOutput            << setw(14) << left << setprecision(4) << fixed << (G_Matrix[i][k]-G_Tref[k])/conversion_h_and_g;
		fOutput			   << endl;
	}

	fOutput << " -----------------------------------------------------------------------------------------------------" << endl;
}

void OpenSMOKE_IdealGas::VerboseDataSpeciesCheckConsistency(ofstream &fOutput)
{
	BzzVector cpiL, d1cpiL, d2cpiL, d3cpiL, d4cpiL;
	BzzVector cpiH, d1cpiH, d2cpiH, d3cpiH, d4cpiH;

	SpeciesCp_ForcingCorrelation(T2, 'L', cpiL, d1cpiL, d2cpiL, d3cpiL, d4cpiL);
	SpeciesCp_ForcingCorrelation(T2, 'H', cpiH, d1cpiH, d2cpiH, d3cpiH, d4cpiH);

	fOutput << "-------------------------------------------------------------------" << endl;
	fOutput << "     Consistency of specific heats (before corrections)            " << endl;
	fOutput << "-------------------------------------------------------------------" << endl;

	fOutput << setw(16) << left << "Name";
	fOutput << setw(12) << left << "T[K]";

	fOutput << setw(16) << left << "CpL";
	fOutput << setw(16) << left << "CpH";
	fOutput << setw(16) << left << "Err(%)";

	fOutput << setw(16) << left << "dCpL/dT";
	fOutput << setw(16) << left << "dCpH/dT";
	fOutput << setw(16) << left << "Err(%)";

	fOutput << endl;

	for (int k=1;k<=NC;k++)
	{
		fOutput << setw(16) << left << names[k];
		fOutput << setw(12) << setprecision(2) << fixed << T2[k];
		
		fOutput << setw(16) << setprecision(6) << left << scientific << cpiL[k];
		fOutput << setw(16) << setprecision(6) << left << scientific << cpiH[k];
		fOutput << setw(16) << setprecision(6) << left << scientific << fabs(cpiL[k]-cpiH[k])/cpiH[k]*100.;

		fOutput << setw(16) << setprecision(6) << left << scientific << d1cpiL[k];
		fOutput << setw(16) << setprecision(6) << left << scientific << d1cpiH[k];
		fOutput << setw(16) << setprecision(6) << left << scientific << fabs(d1cpiL[k]-d1cpiH[k])/d1cpiH[k]*100.;
		
		fOutput << endl;
	}
	fOutput << endl;

	BzzVector phiCorrectionsNorm2(NC);
	BzzVector phiCorrectionsNormI(NC);

	for (int k=1;k<=NC;k++)
	{
		if (fabs(cpiL[k]-cpiH[k])/cpiH[k]*100. < 0.1)
		{	
			BzzVector phi(5);

//			BzzVector b(2, cpiL[k]-cpiH[k], d1cpiL[k]-d1cpiH[k]);

//			BzzFactorizedLQ A(2,5, CpHT[k][1], CpHT[k][2]*T2[k], CpHT[k][3]*T2[k]*T2[k], CpHT[k][4]*T2[k]*T2[k]*T2[k], CpHT[k][5]*T2[k]*T2[k]*T2[k]*T2[k],
//						           0.        , CpHT[k][2],       CpHT[k][3]*T2[k]*2.   , CpHT[k][4]*T2[k]*T2[k]*3.   , CpHT[k][5]*T2[k]*T2[k]*T2[k]*4.);

			BzzVector b(2, cpiL[k]-cpiH[k]);
			BzzFactorizedLQ A(2,5, CpHT[k][1], CpHT[k][2]*T2[k], CpHT[k][3]*T2[k]*T2[k], CpHT[k][4]*T2[k]*T2[k]*T2[k], CpHT[k][5]*T2[k]*T2[k]*T2[k]*T2[k]);

			Solve(A, b, &phi);

			phiCorrectionsNormI[k] = phi.NormI();
			phiCorrectionsNorm2[k] = phi.Norm2();

			for(int i=1;i<=5;i++)
				CpHT[k][i] *= 1.+phi[i];
		}
	}
	fOutput << endl;

    fOutput << "-------------------------------------------------------------------" << endl;
	fOutput << "         Consistency of specific heats (corrections)               " << endl;
	fOutput << "-------------------------------------------------------------------" << endl;
	for (int k=1;k<=NC;k++)
	{
		fOutput << setw(16) << left << names[k];
		fOutput << setw(16) << left << phiCorrectionsNormI[k];
		fOutput << setw(16) << left << phiCorrectionsNorm2[k];
		fOutput << endl;
	}
	fOutput << endl;

	SpeciesCp_ForcingCorrelation(T2, 'L', cpiL, d1cpiL, d2cpiL, d3cpiL, d4cpiL);
	SpeciesCp_ForcingCorrelation(T2, 'H', cpiH, d1cpiH, d2cpiH, d3cpiH, d4cpiH);

	fOutput << "-------------------------------------------------------------------" << endl;
	fOutput << "     Consistency of specific heats (after corrections)             " << endl;
	fOutput << "-------------------------------------------------------------------" << endl;

	fOutput << setw(16) << left << "Name";
	fOutput << setw(12) << left << "T[K]";

	fOutput << setw(16) << left << "CpL";
	fOutput << setw(16) << left << "CpH";
	fOutput << setw(16) << left << "Err(%)";

	fOutput << setw(16) << left << "dCpL/dT";
	fOutput << setw(16) << left << "dCpL/dT";
	fOutput << setw(16) << left << "Err(%)";

	fOutput << endl;

	for (int k=1;k<=NC;k++)
	{
		fOutput << setw(16) << left << names[k];
		fOutput << setw(12) << setprecision(2) << fixed << T2[k];
		
		fOutput << setw(16) << setprecision(6) << left << scientific << cpiL[k];
		fOutput << setw(16) << setprecision(6) << left << scientific << cpiH[k];
		fOutput << setw(16) << setprecision(6) << left << scientific << fabs(cpiL[k]-cpiH[k])/cpiH[k]*100.;

		fOutput << setw(16) << setprecision(6) << left << scientific << d1cpiL[k];
		fOutput << setw(16) << setprecision(6) << left << scientific << d1cpiH[k];
		fOutput << setw(16) << setprecision(6) << left << scientific << fabs(d1cpiL[k]-d1cpiH[k])/d1cpiH[k]*100.;
		
		fOutput << endl;
	}
	fOutput << endl;
}


double OpenSMOKE_IdealGas::GetTemperatureFromMassEnthalpyAndMassFractions(const double TFirstGuess, const double Hmass, BzzVector &omega)
{
	double mw;
	BzzVector xmole(omega.Size());
	GetMWAndMoleFractionsFromMassFractions(mw, xmole, omega);
	return GetTemperatureFromMassEnthalpyAndMoleFractions(TFirstGuess, Hmass, xmole);
}

double OpenSMOKE_IdealGas::GetTemperatureFromMassEnthalpyAndMassFractions(const double TFirstGuess, const double Hmass, BzzVector &omega, const double dtmax)
{
	double mw;
	BzzVector xmole(omega.Size());
	GetMWAndMoleFractionsFromMassFractions(mw, xmole, omega);
	return GetTemperatureFromMassEnthalpyAndMoleFractions(TFirstGuess, Hmass, xmole, dtmax);
}

double OpenSMOKE_IdealGas::GetTemperatureFromMassEnthalpyAndMassFractions(const double TFirstGuess, const double Hmass, BzzVector &omega, 
	const double max_residual, const unsigned int max_newton_iterations, 
	const unsigned int max_searching_iterations, const double max_temperature_increment)
{
	double mw;
	BzzVector xmole(omega.Size());
	GetMWAndMoleFractionsFromMassFractions(mw, xmole, omega);
	return GetTemperatureFromMassEnthalpyAndMoleFractions(TFirstGuess, Hmass, xmole, max_residual, max_newton_iterations, max_searching_iterations, max_temperature_increment);
}

double OpenSMOKE_IdealGas::GetTemperatureFromMassEnthalpyAndMassFractions_NewVersion(const double TFirstGuess, const double Hmass, BzzVector &omega, 
	const double max_residual, const unsigned int max_newton_iterations, 
	const unsigned int max_searching_iterations, const double max_temperature_increment)
{
	double mw;
	BzzVector xmole(omega.Size());
	GetMWAndMoleFractionsFromMassFractions(mw, xmole, omega);
	return GetTemperatureFromMassEnthalpyAndMoleFractions_NewVersion(TFirstGuess, Hmass, xmole, max_residual, max_newton_iterations, max_searching_iterations, max_temperature_increment);
}

double OpenSMOKE_IdealGas::GetTemperatureFromMassEnthalpyAndMoleFractions(const double TFirstGuess, const double Hmass, BzzVector &xmole, const double dtmax)
{		
	double MW		= GetMWFromMoleFractions(xmole);	// [kg/kmol]
	
	double ratio=2.25;
	double TA = TFirstGuess;
	SpeciesCp(TA);
	double CpA = MixCp_FromMoleFractions(xmole);		// [J/kg/K]
	double HA  = GetMixEnthalpy_Mole(TA, xmole) / MW;	// [J/kg]
	double fA = HA-Hmass;
	double TB, CpB, HB, fB;
	double pow1 = 1;
	double incrementDT = 1.;

//	cout << TA << " " << CpA << " " << HA << " " << fA << endl;

	bool iFound = false;
	for(int i=1;i<=20;i++)
	{

		TB = TA + pow1*incrementDT;
		SpeciesCp(TB);
		CpB = MixCp_FromMoleFractions(xmole);		// [J/kg/K]
		HB  = GetMixEnthalpy_Mole(TB, xmole) / MW;	// [J/kg]
		fB  = HB-Hmass;

//		cout << TB << " " << CpB << " " << HB << " " << fB << " " << pow1 << " " << incrementDT << endl;
//		getchar();
		if (fA*fB <= 0.)
		{
			iFound = true;
			break;
		}		
		if (pow1<0)
			incrementDT*=ratio;

		pow1*=-1.;
	}

	if (iFound == false)
		ErrorMessage("Temperature from Enthalpy error (1)");

	for(int i=1;i<=500;i++)
	{
		if (i>=500)
			ErrorMessage("Temperature from Enthalpy error (2)");

		double TC = 0.50*(TA+TB);
		double CpC = MixCp_FromMoleFractions(xmole);		// [J/kg/K]
		double HC  = GetMixEnthalpy_Mole(TC, xmole) / MW;	// [J/kg]
		double fC  = HC-Hmass;

		if (fabs(TA-TB) < dtmax)
		{
	//		cout << TC << endl;
	//		cout << fA << " " << fB << " " << fC << endl;
			return TC;
		}

		if (fA*fC <= 0.)
		{
			TB = TC;
			CpB = CpC;
			HB = HC;
			fB = fC;
		}
		else
		{
			TA = TC;
			CpA = CpC;
			HA = HC;
			fA = fC;
		}
	}

	return -1;
}

double OpenSMOKE_IdealGas::GetTemperatureFromMassEnthalpyAndMoleFractions(const double TFirstGuess, const double Hmass, BzzVector &xmole)
{		
	double dt;
	double dtold	= 0.;
	double TE		= TFirstGuess;						// First guess temperature
	double TEold	= TFirstGuess;						// First guess temperature
	double MW		= GetMWFromMoleFractions(xmole);	// [kg/kmol]
	
	
	bool iConvergence = false;
	for (int n=0;n<500;n++) 
	{
		SpeciesCp(TE);
		
		double CpE = MixCp_FromMoleFractions(xmole);		// [J/kg/K]
		double HE  = GetMixEnthalpy_Mole(TE, xmole) / MW;	// [J/kg]
	
        dt = (Hmass - HE)/CpE;

//		cout << "n=" << n << " " << "TE=" << TE << " TEOld=" << TEold << " Diff=" << fabs(TE-TEold) << endl;

		if (n>1 && dt*dtold<=0.)
		{
			double xA, xB, fA, fB;

			if (TE>=TEold)	{ xA = TEold; fA=dtold; xB=TE; fB=dt;}
			else			{ xB = TEold; fB=dtold; xA=TE; fA=dt;}

			for (int j=1;j<=300;j++)
			{
				double m = (fA-fB)/(xA-xB);
				double q = fB-m*xB;
				TE = -q/m;

				SpeciesCp(TE);
		
				double CpE = MixCp_FromMoleFractions(xmole);		// [J/kg/K]
				double HE  = GetMixEnthalpy_Mole(TE, xmole) / MW;	// [J/kg]
	
				dt = (Hmass - HE)/CpE;

				if (dt*fA<=0.)	{xB=TE; fB=dt;}
				else			{xA=TE; fA=dt;}

//				cout << "  j=" << j << " " << "TE=" << TE << " TEOld" << TEold  << " Diff=" << fabs(TE-TEold) << endl;

				if (j>1)
					if (fabs(dt) <= 0.0001)
					{
						double m = (fA-fB)/(xA-xB);
						double q = fB-m*xB;
						TE = -q/m;
						
						iConvergence = true;
						break;
					}
			}
		}
		else
		{
			if (dt > 300.0)			dt =  300.0;
			else if (dt < -300.0)	dt = -300.0; 
		
			TEold = TE;
			dtold = dt;

			TE += dt;
			if (fabs(dt) <= 0.0001)
			{
				iConvergence = true;
				break;
			}

//			cout << "  n-j=" << n << " " << "TE=" << TE << " TEOld" << TEold  << " Diff=" << fabs(TE-TEold) << endl;
		}

		if (iConvergence == true)
			break;
	}

	if (iConvergence == false)
		ErrorMessage("Maximum number of iteration for calculating temperature from enthalpy...");

	return TE;
}

double OpenSMOKE_IdealGas::GetTemperatureFromMassEnthalpyAndMoleFractions(const double TFirstGuess, const double Hmass, BzzVector &xmole, 
	const double max_residual, const unsigned int max_newton_iterations, 
	const unsigned int max_searching_iterations, const double max_temperature_increment)
{		
	double TE, CpE, HE, residual;
	
	// Start calculations
	double TEold	= TFirstGuess;						// First guess temperature
	double MW		= GetMWFromMoleFractions(xmole);	// [kg/kmol]
	SpeciesCp(TEold);
	double CpEold = MixCp_FromMoleFractions(xmole);			// [J/kg/K]
	double HEold  = GetMixEnthalpy_Mole(TEold, xmole) / MW;	// [J/kg]
	double residualOld = (Hmass - HEold)/CpEold;			// [J/kg]

	// If the point is already OK
	if (fabs(residualOld) <= max_residual)
		return TEold;

	// Try Newton's method
	{
		unsigned int newton_iteration = 0;
		for (;;)
		{
			newton_iteration++;

			TE = TEold + residualOld;
			CpE = MixCp_FromMoleFractions(xmole);			// [J/kg/K]
			HE  = GetMixEnthalpy_Mole(TE, xmole) / MW;	// [J/kg]
			residual = (Hmass - HE)/CpE;

			if ( fabs(residual) > fabs(residualOld) )
			{
				residualOld = residual;
				CpEold = CpE;
				HEold = HE;
				break;
			}
			else
			{
				if (fabs(residual) < max_residual)
					return TE;
				residualOld = residual;
				CpEold = CpE;
				HEold = HE;
				TEold = TE;
			}

			if (newton_iteration > max_newton_iterations)
				ErrorMessage("Maximum number of newton's iterations (1) for calculating temperature from enthalpy.");
		}
	}

	// Search for a suitable interval
	{
		int iteration = 1;
		double delta = 10.;

		TE = TEold;
		for (;;)
		{ 
			CpE = MixCp_FromMoleFractions(xmole);			// [J/kg/K]

			TEold -= delta;
			if (TEold < 250.) TEold = 250.;
			HEold  = GetMixEnthalpy_Mole(TEold, xmole) / MW;	// [J/kg]
			CpEold = MixCp_FromMoleFractions(xmole);			// [J/kg/K]
			residualOld = (Hmass - HEold)/CpEold;

			TE    += delta;
			if (TE > 8000.) TE = 8000.;
			HE  = GetMixEnthalpy_Mole(TE, xmole) / MW;	// [J/kg]
			CpE = MixCp_FromMoleFractions(xmole);		// [J/kg/K]
			residual = (Hmass - HE)/CpE;

			//std::cout << TEold << " " << residualOld << " " << TE << " " << residual << endl;

			if (residual*residualOld <= 0.)
				break;

			if (iteration==max_searching_iterations)
				ErrorMessage("Search for a suitable interval: Unable to find the appropriate interval for calculating temperature from enthalpy.");
		}
	}

/*	else
	{
		TE = TEold - residualOld;
		CpE = MixCp_FromMoleFractions(xmole);			// [J/kg/K]
		HE  = GetMixEnthalpy_Mole(TE, xmole) / MW;	// [J/kg]
		residual = (Hmass - HE)/CpE;

		// Find interval
		unsigned int iteration = 0;
		for(;;)
		{
			iteration++;

			//cout << iteration << " " << TE << " " << residual << " " << TEold << " " << residualOld << std::endl;

			if (residual*residualOld <= 0.)
				break;
			else
			{
				if (residual < residualOld)
				{
					TE += max_temperature_increment;
					CpE = MixCp_FromMoleFractions(xmole);			// [J/kg/K]
					HE  = GetMixEnthalpy_Mole(TE, xmole) / MW;	// [J/kg]
					residual = (Hmass - HE)/CpE;
				}
				else
				{
					TEold += max_temperature_increment;
					CpEold = MixCp_FromMoleFractions(xmole);			// [J/kg/K]
					HEold  = GetMixEnthalpy_Mole(TEold, xmole) / MW;	// [J/kg]
					residualOld = (Hmass - HEold)/CpEold;
				}
			}
			
			if (iteration > max_searching_iterations)
				ErrorMessage("Unable to find the appropriate interval for calculating temperature from enthalpy.");
		}
*/
	// Reduce the interval
	{
		double XA, XB, fA, fB;

		if (residualOld > residual)	{ XA = TEold; fA=residualOld; XB=TE; fB=residual;}
		else						{ XB = TEold; fB=residualOld; XA=TE; fA=residual;}

	//	std::cout << "A: " << XA << " " << fA << std::endl;
	//	std::cout << "B: " << XB << " " << fB << std::endl;

		unsigned int iteration = 0;
		for (;;)
		{
			iteration++;

			double XC = (XA+XB)*0.50;
			SpeciesCp(XC);
			double CpC = MixCp_FromMoleFractions(xmole);		// [J/kg/K]
			double HC  = GetMixEnthalpy_Mole(XC, xmole) / MW;	// [J/kg]
			double residualC = (Hmass-HC)/CpC;

			cout << iteration << " " << XC << " " << residualC << " " << endl;

			if (fabs(residualC) <= 1.e-6)
			{
				TEold = XC;
				HEold = HC;
				residualOld = residualC;
				cout << "Break ..." << iteration << " " << endl;
				return XC;
				break;
			}

			if (fB*residualC>=0.)	{XB=XC; fB=residualC;}
			else					{XA=XC; fA=residualC;}

			if (iteration > max_searching_iterations)
				ErrorMessage("Unable to find the halving the interval for calculating temperature from enthalpy.");
		}
	}

	// Try Newton again
	{
		cout << "Again newton ..." << " " << endl;
		unsigned int newton_iteration = 0;
		residualOld = 0.;
		for (;;)
		{
			newton_iteration++;

			TE = TEold + residualOld;
			CpE = MixCp_FromMoleFractions(xmole);			// [J/kg/K]
			HE  = GetMixEnthalpy_Mole(TE, xmole) / MW;	// [J/kg]
			residual = (Hmass - HE)/CpE;

	//		if ( fabs(residual) > fabs(residualOld) )
	//		{
	//			residualOld = residual;
	//			CpEold = CpE;
	//			HEold = HE;
	//			if (newton_iteration > max_newton_iterations)
		//		ErrorMessage("Maximum number of newton's iterations (3) for calculating temperature from enthalpy.");
				//break;
		//	}
		//	else
			{
				if (fabs(residual) < max_residual)
					return TE;
				residualOld = residual;
				CpEold = CpE;
				HEold = HE;
				TEold = TE;
			}

			if (newton_iteration > max_newton_iterations)
				ErrorMessage("Maximum number of newton's iterations (2) for calculating temperature from enthalpy.");
		}
	}
}

double OpenSMOKE_IdealGas::GetTemperatureFromMassEnthalpyAndMoleFractions_NewVersion(
	const double TFirstGuess,	const double Hmass, BzzVector &xmole, 
	const double max_residual,	const unsigned int max_newton_iterations, 
	const unsigned int max_searching_iterations, const double max_temperature_increment)
{		
	// Start calculations
	double TStarting	= TFirstGuess;							// First guess temperature
	double MW			= GetMWFromMoleFractions(xmole);		// [kg/kmol]
	SpeciesCp(TStarting);
	double CpEStarting = MixCp_FromMoleFractions(xmole);			// [J/kg/K]
	double HEStarting  = GetMixEnthalpy_Mole(TStarting, xmole) / MW;	// [J/kg]
	double residualStarting = (Hmass - HEStarting)/CpEStarting;			// [J/kg]

	// If the point is already OK
	if (fabs(residualStarting) <= max_residual)
		return TStarting;

	// Try Newton's method
	{
		double TEold = TStarting;
		double residualOld = residualStarting;
		double CpEold, HEold;
		unsigned int newton_iteration = 0;
		for (;;)
		{
			newton_iteration++;

			const double TE = TEold + residualOld;
			const double CpE = MixCp_FromMoleFractions(xmole);		// [J/kg/K]
			const double HE  = GetMixEnthalpy_Mole(TE, xmole) / MW;	// [J/kg]
			const double residual = (Hmass - HE)/CpE;

			// In case of failure of Newton's method
			if ( fabs(residual) > fabs(residualOld) )
			{
				break;
			}
			// If the Newton's method works
			else
			{
				if (fabs(residual) <= max_residual)
					return TE;

				residualOld = residual;
				CpEold = CpE;
				HEold = HE;
				TEold = TE;
			}

			if (newton_iteration > max_newton_iterations)
				break;
		//		ErrorMessage("Maximum number of newton's iterations (1) for calculating temperature from enthalpy.");
		}
	} // end Newton's method

	// Search for a suitable interval
	double T1, T2;
	double residual1, residual2;
	{
		for (int j=1;j<=max_searching_iterations;j++)
		{ 
			T1 = TStarting - max_temperature_increment*j;
			if (T1 < 250.) T1 = 250.;
			const double HE1  = GetMixEnthalpy_Mole(T1, xmole) / MW;	// [J/kg]
			const double CpE1 = MixCp_FromMoleFractions(xmole);			// [J/kg/K]
			residual1 = (Hmass - HE1)/CpE1;

			T2 = TStarting + max_temperature_increment*j;
			if (T2 > 8000.) T2 = 8000.;
			const double HE2  = GetMixEnthalpy_Mole(T2, xmole) / MW;	// [J/kg]
			const double CpE2 = MixCp_FromMoleFractions(xmole);		// [J/kg/K]
			residual2 = (Hmass - HE2)/CpE2;

			if (residual1*residual2 <= 0.)
				break;

			if (j==max_searching_iterations)
				ErrorMessage("Search for a suitable interval: Unable to find the appropriate interval for calculating temperature from enthalpy.");
		}
	}

	// Reduce the interval
	{
		double XA, XB, fA, fB;

		XA = T1; 
		fA=residual1; 
		XB=T2; 
		fB=residual2;

		unsigned int iteration = 0;
		for (;;)
		{
			iteration++;

			double XC = (XA+XB)*0.50;
			SpeciesCp(XC);
			double CpC = MixCp_FromMoleFractions(xmole);		// [J/kg/K]
			double HC  = GetMixEnthalpy_Mole(XC, xmole) / MW;	// [J/kg]
			double residualC = (Hmass-HC)/CpC;

			if (fabs(residualC) <= max_residual)
			{
				return XC;
				break;
			}

			if (iteration == 100)
			{
				XA -= 0.01;
				XB += 0.01;
				continue;
			}

			if (iteration == 120)
				return XC;

			if (fB*residualC>=0.)	{XB=XC; fB=residualC;}
			else					{XA=XC; fA=residualC;}


		}
	}
}

double OpenSMOKE_IdealGas::GetTemperatureFromMassInternalEnergyAndMoleFractions(const double TFirstGuess, const double Umass, BzzVector &xmole)
{		
	double dt;
	double dtold	= 0.;
	double TE		= TFirstGuess;						// First guess temperature
	double TEold	= TFirstGuess;						// First guess temperature
	double MW		= GetMWFromMoleFractions(xmole);	// [kg/kmol]
	
	bool iConvergence = false;
	for (int n=0;n<50;n++) 
	{
		SpeciesCp(TE);
		
		double CpE = MixCp_FromMoleFractions(xmole);			// [J/kg/K]
		double UE  = GetMixInternalEnergy_Mole(TE, xmole) / MW;	// [J/kg]
	
        dt = (Umass-UE)/CpE;

		if (n>1 && dt*dtold<=0.)
		{
			double xA, xB, fA, fB;

			if (TE>=TEold)	{ xA = TEold; fA=dtold; xB=TE; fB=dt;}
			else			{ xB = TEold; fB=dtold; xA=TE; fA=dt;}

			for (int j=1;j<=30;j++)
			{
				double m = (fA-fB)/(xA-xB);
				double q = fB-m*xB;
				TE = -q/m;
				
				SpeciesCp(TE);
		
				double CpE = MixCp_FromMoleFractions(xmole);			// [J/kg/K]
				double UE  = GetMixInternalEnergy_Mole(TE, xmole) / MW;	// [J/kg]
	
				dt = (Umass - UE)/CpE;

				if (dt*fA<=0.)	{xB=TE; fB=dt;}
				else			{xA=TE; fA=dt;}

				if (j>1)
					if (fabs(dt) <= 0.001)
					{
						TE	= (xA+xB)/2.;
						iConvergence = true;
						break;
					}
			}
		}
		else
		{
			if (dt > 300.0)			dt =  300.0;
			else if (dt < -300.0)	dt = -300.0; 
		
			TEold = TE;
			dtold = dt;

			TE += dt;
			if (fabs(dt) <= 0.001)
			{
				iConvergence = true;
				break;
			}
		}

		if (iConvergence == true)
			break;
	}

	if (iConvergence == false)
		ErrorMessage("Maximum number of iteration for calculating temperature from internal energy...");

	return TE;
}

double OpenSMOKE_IdealGas::GetTemperatureFromMassEntropyAndMoleFractions(const double TFirstGuess, const double P_Pascal, const double Smass, BzzVector &xmole)
{		
	double dt;
	double dtold	= 0.;
	double TE		= TFirstGuess;						// First guess temperature
	double TEold	= TFirstGuess;						// First guess temperature
	double MW		= GetMWFromMoleFractions(xmole);	// [kg/kmol]
	
	bool iConvergence = false;
	for (int n=0;n<50;n++) 
	{
		SpeciesCp(TE);
		
		double CpE = MixCp_FromMoleFractions(xmole);				// [J/kg/K]
		double SE  = GetMixEntropy_Mole(P_Pascal, TE, xmole) / MW;	// [J/kg/K]
	
        dt = (Smass-SE)/CpE*TE;

		if (n>1 && dt*dtold<=0.)
		{
			double xA, xB, fA, fB;

			if (TE>=TEold)	{ xA = TEold; fA=dtold; xB=TE; fB=dt;}
			else			{ xB = TEold; fB=dtold; xA=TE; fA=dt;}

			for (int j=1;j<=30;j++)
			{
				double m = (fA-fB)/(xA-xB);
				double q = fB-m*xB;
				TE = -q/m;

				SpeciesCp(TE);
		
				double CpE = MixCp_FromMoleFractions(xmole);				// [J/kg/K]
				double SE  = GetMixEntropy_Mole(P_Pascal, TE, xmole) / MW;	// [J/kg/K]
	
				dt = (Smass-SE)/CpE*TE;

				if (dt*fA<=0.)	{xB=TE; fB=dt;}
				else			{xA=TE; fA=dt;}

				if (j>1)
					if (fabs(dt) <= 0.001)
					{
						double m = (fA-fB)/(xA-xB);
						double q = fB-m*xB;
						TE = -q/m;
						
						iConvergence = true;
						break;
					}
			}
		}
		else
		{
			if (dt > 300.0)			dt =  300.0;
			else if (dt < -300.0)	dt = -300.0; 
		
			TEold = TE;
			dtold = dt;

			TE += dt;
			if (fabs(dt) <= 0.001)
			{
				iConvergence = true;
				break;
			}
		}

		if (iConvergence == true)
			break;
	}

	if (iConvergence == false)
		ErrorMessage("Maximum number of iteration for calculating temperature from entropy...");

	return TE;
}

double OpenSMOKE_IdealGas::GetPressureFromMassEntropyAndMoleFractions(const double T, const double Smass, BzzVector &xmole)
{		
	double MW	   = GetMWFromMoleFractions(xmole);							// [kg/kmol]
	double S_Act   = Smass*MW;												// [J/kmol/K]
	double S_Ref   = GetMixEntropy_Mole(Constants::P_Reference, T, xmole);	// [J/kmol/K]
	
	return Constants::P_Reference*exp((-S_Act+S_Ref)/Constants::R_J_kmol);	// [Pa]			
}

void OpenSMOKE_IdealGas::Equilibrium_TP(BzzVector &xmole_E, double &N_E, const double T, const double P_Pa, BzzVector &x_elements, const bool iVerbose)
{
	OpenSMOKE_EquilibriumStanjan equilibrium;
	if (iVerbose == true) equilibrium.SetVerbose();
	equilibrium.Setup(this);	
	equilibrium.SetElementalComposition(x_elements);
	equilibrium.Equilibrate(T, P_Pa, xmole_E, N_E);
}

void OpenSMOKE_IdealGas::Equilibrium_HP(double &T_E, BzzVector &xmole_E, double &N_E, const double Hmass, const double P_Pa, BzzVector &x_elements, const bool iVerbose, int &flag)
{
//	ofstream fEq;openOutputFileAndControl(fEq, "fEq.out");
//	fEq.setf(ios::scientific);

	OpenSMOKE_EquilibriumStanjan equilibrium;
	if (iVerbose == true) equilibrium.SetVerbose();
	equilibrium.Setup(this);	
	equilibrium.SetElementalComposition(x_elements);
		
	double dt;
	double dtold	= 0.;
	double T_Eold	= T_E;		// First guess temperature
	
	bool iConvergence = false;
	for (int n=0;n<50;n++) 
	{
	//	cout << "Out It " << n << "Start: " << T_E << endl;

		SpeciesCp(T_E);
	
		equilibrium.Equilibrate(T_E, P_Pa, xmole_E, N_E);

		double Cp_E = MixCp_FromMoleFractions(xmole_E);			// [J/kg/K]
		double MW_E = GetMWFromMoleFractions(xmole_E);			// [kg/kmol]
		double H_E  = GetMixEnthalpy_Mole(T_E, xmole_E) / MW_E;	// [J/kg]
        dt = (Hmass-H_E)/Cp_E;
		
		if (n>1 && dt*dtold<=0.)
		{
	//		cout << "Entering dim " << T_E << " " << T_Eold << " " << dt << " " << dtold << endl;
			double xA, xB, fA, fB;

			if (T_E>=T_Eold)	{ xA = T_Eold; fA=dtold; xB=T_E; fB=dt;}
			else				{ xB = T_Eold; fB=dtold; xA=T_E; fA=dt;}

			double dtmin = -1.e16; double T_E_min;
			double dtmax =  1.e16; double T_E_max;
			for (int j=1;j<=20;j++)
			{
				T_E = (xA+xB)/2.;

				SpeciesCp(T_E);	
				equilibrium.Equilibrate(T_E, P_Pa, xmole_E, N_E);
				double Cp_E = MixCp_FromMoleFractions(xmole_E);			// [J/kg/K]
				double MW_E = GetMWFromMoleFractions(xmole_E);			// [kg/kmol]
				double H_E  = GetMixEnthalpy_Mole(T_E, xmole_E) / MW_E;	// [J/kg]
				dt = (Hmass - H_E)/Cp_E;
				if (dt>=0.	&& dt<dtmax)	{dtmax = dt; T_E_max=T_E;}
				if (dt<0.	&& dt>dtmin)	{dtmin = dt; T_E_min=T_E;}

//				cout << "     Dim " << j << " " << T_E << " " << dt << " " << fA << " " << fB << " " << xA << " " << xB << endl;
//				cout << "         " << Hmass << " " << H_E << " " << Cp_E << " " << MW_E << endl;
//				fEq << j << " " << T_E << " " << dt << " " << fA << " " << fB << " " << xA << " " << xB 
//								 << " " << Hmass << " " << H_E << " " << Cp_E << " " << MW_E << endl;

				if (fabs(xA-xB)<=0.01 && fabs(dt)<=0.001)
				{
					iConvergence = true;
					flag = 1;
					break;
				}

				if ( (dt*fA)<=0.)	{xB=T_E; fB=dt;}
				else				{xA=T_E; fA=dt;}
			}

			if (iConvergence == true)
				break;
			else
			{	
				T_E=0.50*(T_E_min+T_E_max);
				
				SpeciesCp(T_E);	
				equilibrium.Equilibrate(T_E, P_Pa, xmole_E, N_E);
				double Cp_E = MixCp_FromMoleFractions(xmole_E);			// [J/kg/K]
				double MW_E = GetMWFromMoleFractions(xmole_E);			// [kg/kmol]
				double H_E  = GetMixEnthalpy_Mole(T_E, xmole_E) / MW_E;	// [J/kg]
				dt = (Hmass - H_E)/Cp_E;
				
				iConvergence=true;
				flag = 2;
			}
		}
		else
		{
			// limit step size to 200 K
			if (dt > 100.0)			dt =  100.0 + 0.2*double(n);
			else if (dt < -100.0)	dt = -100.0 - 0.2*double(n); 
		
			T_Eold = T_E;
			dtold = dt;

			T_E += dt;

			if (fabs(dt) <= 0.001)
			{
				iConvergence = true;
				flag = 1;
				break;
			}
		}

		if (iConvergence == true)
			break;
	}
	
	if (iConvergence == false)
		ErrorMessage("It was impossible to solve the (H,P) equilibrium non linear system...");
}


void OpenSMOKE_IdealGas::Equilibrium_SP(double &T_E, BzzVector &xmole_E, double &N_E, const double Smass, const double P_Pa, BzzVector &x_elements, const bool iVerbose)
{
	OpenSMOKE_EquilibriumStanjan equilibrium;
	if (iVerbose == true) equilibrium.SetVerbose();	
	equilibrium.Setup(this);
	equilibrium.SetElementalComposition(x_elements);
		
	double dt;
	double dtold	= 0.;
	T_E		= 300.;				// First guess temperature
	double T_Eold	= T_E;		// First guess temperature
	
	bool iConvergence = false;
	for (int n=0;n<50;n++) 
	{
		SpeciesCp(T_E);
		
		equilibrium.Equilibrate(T_E, P_Pa, xmole_E, N_E);

		double Cp_E = MixCp_FromMoleFractions(xmole_E);			// [J/kg/K]
		double MW_E = GetMWFromMoleFractions(xmole_E);			// [kg/kmol]
		double S_E  = GetMixEntropy_Mole(P_Pa, T_E, xmole_E) / MW_E;	// [J/kg]
	
        dt = (Smass - S_E)/Cp_E*T_E;

		if (n>1 && dt*dtold<=0.)
		{
			double xA, xB, fA, fB;

			if (T_E>=T_Eold)	{ xA = T_Eold; fA=dtold; xB=T_E; fB=dt;}
			else				{ xB = T_Eold; fB=dtold; xA=T_E; fA=dt;}

			for (int j=1;j<=30;j++)
			{
				double m = (fA-fB)/(xA-xB);
				double q = fB-m*xB;
				T_E = -q/m;

				SpeciesCp(T_E);
		
				equilibrium.Equilibrate(T_E, P_Pa, xmole_E, N_E);

				double Cp_E = MixCp_FromMoleFractions(xmole_E);					// [J/kg/K]
				double MW_E = GetMWFromMoleFractions(xmole_E);					// [kg/kmol]
				double S_E  = GetMixEntropy_Mole(P_Pa, T_E, xmole_E) / MW_E;	// [J/kg]
	
				dt = (Smass - S_E)/Cp_E*T_E;

				if (dt*fA<=0.)	{xB=T_E; fB=dt;}
				else			{xA=T_E; fA=dt;}

				if (j>1)
					if (fabs(dt) <= 0.001)
					{
						double m = (fA-fB)/(xA-xB);
						double q = fB-m*xB;
						T_E = -q/m;

						iConvergence = true;
						break;
					}
			}
		}
		else
		{
			// limit step size to 100 K
			if (dt > 300.0)			dt =  300.0;
			else if (dt < -300.0)	dt = -300.0; 
		
			T_Eold = T_E;
			dtold = dt;

			T_E += dt;
			if (fabs(dt) <= 0.001)
			{
				iConvergence = true;
				//flag = 1;
				break;
			}
		}

		if (iConvergence == true)
			break;
	}	

	if (iConvergence == false)
		ErrorMessage("It was impossible to solve the (S,P) equilibrium non linear system...");

	equilibrium.Equilibrate(T_E, P_Pa, xmole_E, N_E);
}

void OpenSMOKE_IdealGas::Equilibrium_TRHO(double &P_E, BzzVector &xmole_E, double &N_E, const double T, const double rho, BzzVector &x_elements, const bool iVerbose)
{
	OpenSMOKE_EquilibriumStanjan equilibrium;
	if (iVerbose == true) equilibrium.SetVerbose();	
	equilibrium.Setup(this);
	equilibrium.SetElementalComposition(x_elements);
	
	double dp;
	double dpold	= 0.;
	P_E				= 101325.;	// First guess pressure
	double P_Eold	= P_E;		// First guess pressure
	
	bool iConvergence = false;
	bool iNewton = false;
	for (int n=0;n<500;n++) 
	{
		equilibrium.Equilibrate(T, P_E, xmole_E, N_E);

		double MW_E		= GetMWFromMoleFractions(xmole_E);		// [kg/kmol]
		double rho_E	= P_E*MW_E/Constants::R_J_kmol/T;		// [kg/m3]	

        dp = (rho - rho_E)*Constants::R_J_kmol*T/MW_E;

		if (n>1 && dp*dpold<=0.)
		{
			iNewton = true;

			double xA, xB, fA, fB;

			if (P_E>=P_Eold)	{ xA = P_Eold; fA=dpold; xB=P_E; fB=dp;}
			else				{ xB = P_Eold; fB=dpold; xA=P_E; fA=dp;}

			for (int j=1;j<=30;j++)
			{
				double m = (fA-fB)/(xA-xB);
				double q = fB-m*xB;
				P_E = -q/m;
		
				equilibrium.Equilibrate(T, P_E, xmole_E, N_E);

				double MW_E		= GetMWFromMoleFractions(xmole_E);		// [kg/kmol]
				double rho_E	= P_E*MW_E/Constants::R_J_kmol/T;		// [kg/m3]	
				
				dp = (rho - rho_E)*Constants::R_J_kmol*T/MW_E;

				if (dp*fA<=0.)	{xB=P_E; fB=dp;}
				else			{xA=P_E; fA=dp;}

				if (j>1)
					if (fabs(dp) <= 0.1)
					{
						double m = (fA-fB)/(xA-xB);
						double q = fB-m*xB;
						P_E = -q/m;
						
						iConvergence = true;
						break;
					}
			}
		}
		else
		{
			// limit step size to 1. bar
			if (dp > 1.e5)			dp =  1.e5;
			else if (dp < -1.e5)	dp = -1.e5; 
			if (P_E+dp <= 0.)		dp = -(P_E-1.);

			P_Eold = P_E;
			dpold = dp;

			P_E += dp;

			if ( (fabs(dp) <= 0.1) && n>1)
			{
				break;
				iConvergence = true;
			}
		}

		if (iConvergence == true)
			break;
	}	

	if (iConvergence == false)
		if (iNewton == true)	ErrorMessage("It was impossible to solve the (T,RHO) or (T,V) equilibrium non linear system (newton number of iterations)...");
		else					ErrorMessage("It was impossible to solve the (T,RHO) or (T,V) equilibrium non linear system (no newton)...");
		
	equilibrium.Equilibrate(T, P_E, xmole_E, N_E);
}

void OpenSMOKE_IdealGas::Equilibrium_URHO(double &T_E, double &P_E, BzzVector &xmole_E, double &N_E, const double Umass, const double rho, BzzVector &x_elements, const bool iVerbose)
{
	OpenSMOKE_Equilibrium_URHO_MyNonLinearSystem nls;
	BzzVector xMin(2, 200.,  Constants::P_Reference/1.e5);
	BzzVector xMax(2, 6000., Constants::P_Reference*1.e5);
	BzzVector xFirstGuess(2, Constants::T_Reference, Constants::P_Reference);

	nls.assignIdealGas(this, Umass, rho, x_elements, iVerbose);
	BzzNonLinearSystemObject o(xFirstGuess, &nls);

	o.SetMinimumConstraints(xMin);
	o.SetMaximumConstraints(xMax);

	char control = o();

	if (control<=1 || control>=7)
		ErrorMessage("It was impossible to solve the (U,RHO) or (U,V) equilibrium non linear system...");

	T_E		= nls.T_E;
	P_E		= nls.P_E;
	N_E		= nls.N_E;
	xmole_E = nls.xmole_E;
}

void OpenSMOKE_IdealGas::Equilibrium_SRHO(double &T_E, double &P_E, BzzVector &xmole_E, double &N_E, const double Smass, const double rho, BzzVector &x_elements, const bool iVerbose)
{
	OpenSMOKE_Equilibrium_SRHO_MyNonLinearSystem nls;
	BzzVector xMin(2, 200.,  Constants::P_Reference/1.e5);
	BzzVector xMax(2, 6000., Constants::P_Reference*1.e5);
	BzzVector xFirstGuess(2, Constants::T_Reference, Constants::P_Reference);

	nls.assignIdealGas(this, Smass, rho, x_elements, iVerbose);
	BzzNonLinearSystemObject o(xFirstGuess, &nls);

	o.SetMinimumConstraints(xMin);
	o.SetMaximumConstraints(xMax);

	char control = o();

	if (control<=1 || control>=7)
		ErrorMessage("It was impossible to solve the (S,RHO) or (S,V) equilibrium non linear system...");

	T_E		= nls.T_E;
	P_E		= nls.P_E;
	N_E		= nls.N_E;
	xmole_E = nls.xmole_E;
}

void OpenSMOKE_Equilibrium_URHO_MyNonLinearSystem::GetResiduals(BzzVector &x,BzzVector &f)
{
	T_E = x[1];
	P_E = x[2];

	equilibrium.Equilibrate(T_E, P_E, xmole_E, N_E);

	double MW_E		= ptMix->GetMWFromMoleFractions(xmole_E);				// [kg/kmol]
	double rho_E	= P_E*MW_E/Constants::R_J_kmol/T_E;						// [kg/m3]
	double U_E		= ptMix->GetMixInternalEnergy_Mole(T_E, xmole_E)/MW_E;	// [J/kg]
	
	f[1] = U_E-Umass;
	f[2] = rho_E-rho;
}

void OpenSMOKE_Equilibrium_URHO_MyNonLinearSystem::assignIdealGas(OpenSMOKE_IdealGas *mix, const double _Umass, const double _rho, BzzVector &_x_elements, const bool iVerbose)
{
	ptMix		=  mix;
	Umass		= _Umass;
	rho			= _rho;
	x_elements	= _x_elements;
	ChangeDimensions(ptMix->NumberOfSpecies(), &xmole_E);

	if (iVerbose == true) equilibrium.SetVerbose();	
	equilibrium.Setup(ptMix);
	equilibrium.SetElementalComposition(x_elements);
}

void OpenSMOKE_Equilibrium_SRHO_MyNonLinearSystem::GetResiduals(BzzVector &x,BzzVector &f)
{
	T_E = x[1];
	P_E = x[2];

	equilibrium.Equilibrate(T_E, P_E, xmole_E, N_E);

	double MW_E		= ptMix->GetMWFromMoleFractions(xmole_E);				// [kg/kmol]
	double rho_E	= P_E*MW_E/Constants::R_J_kmol/T_E;						// [kg/m3]
	double S_E		= ptMix->GetMixEntropy_Mole(P_E, T_E, xmole_E)/MW_E;	// [J/kg/K]
	
	f[1] = S_E-Smass;
	f[2] = rho_E-rho;
}

void OpenSMOKE_Equilibrium_SRHO_MyNonLinearSystem::assignIdealGas(OpenSMOKE_IdealGas *mix, const double _Smass, const double _rho, BzzVector &_x_elements, const bool iVerbose)
{
	ptMix		=  mix;
	Smass		= _Smass;
	rho			= _rho;
	x_elements	= _x_elements;
	ChangeDimensions(ptMix->NumberOfSpecies(), &xmole_E);

	if (iVerbose == true) equilibrium.SetVerbose();	
	equilibrium.Setup(ptMix);
	equilibrium.SetElementalComposition(x_elements);
}

BzzVector OpenSMOKE_IdealGas::GetMoleFractionsFromEquivalenceRatio(const double equivalence_ratio, const std::string fuel_name)
{
	double nC	= elements[recognize_element("c")][recognize_species(fuel_name)];
	double nH	= elements[recognize_element("h")][recognize_species(fuel_name)];

	double nO2		= (2*nC + 0.50*nH)/2.;
	double nN2		= 0.79/0.21*nO2;
	double nFuel	= equivalence_ratio;
	double n		= nO2 + nN2 + nFuel;

	BzzVector mole(NumberOfSpecies());	
	mole[recognize_species(fuel_name)]	= nFuel/n;
	mole[recognize_species("O2")]		= nO2/n;
	mole[recognize_species("N2")]		= nN2/n;

	return mole;
}

BzzVector OpenSMOKE_IdealGas::GetMoleFractionsFromEquivalenceRatio(const double equivalence_ratio,	const vector<string> fuel_names,		BzzVector &moles_fuel,
																									const vector<string> oxidizer_names,	BzzVector &moles_oxidizer)
{
	int j;
	int number_of_fuels		= moles_fuel.Size();
	int number_of_oxidizers	= moles_oxidizer.Size();

	BzzVector nC(number_of_fuels);
	BzzVector nH(number_of_fuels);
	BzzVector nO(number_of_fuels);
	BzzVector nOxidizers(number_of_oxidizers);

	int jC = recognize_element_without_error("c");
	int jH = recognize_element_without_error("h");
	int jO = recognize_element_without_error("o");
	
	if (jC>0)
		for(j=1;j<=number_of_fuels;j++)
			nC[j]	= elements[jC][recognize_species(fuel_names[j])];

	if (jH>0)
		for(j=1;j<=number_of_fuels;j++)
			nH[j]	= elements[jH][recognize_species(fuel_names[j])];

	if (jO>0)
		for(j=1;j<=number_of_fuels;j++)
			nO[j]	= elements[jO][recognize_species(fuel_names[j])];

	double nFuel	= equivalence_ratio*moles_fuel.GetSumElements();

	int jO2 = 0;
	for(j=1;j<=number_of_oxidizers;j++)
		if (oxidizer_names[j] == "O2")	{ jO2 = j; break;}

	double nO2		= (	2.*Dot(nC,moles_fuel) + 
						0.50*Dot(nH,moles_fuel) -
						Dot(nO,moles_fuel))/2.;

	for(j=1;j<=number_of_oxidizers;j++)
		nOxidizers[j] = moles_oxidizer[j]/moles_oxidizer[jO2]*nO2;
	
	double n		= nFuel + nOxidizers.GetSumElements();

	BzzVector mole(NumberOfSpecies());	
	for(j=1;j<=number_of_fuels;j++)
		mole[recognize_species(fuel_names[j])] += equivalence_ratio*moles_fuel[j]/n;
	for(j=1;j<=number_of_oxidizers;j++)
		mole[recognize_species(oxidizer_names[j])] += nOxidizers[j]/n;

	return mole;
}
