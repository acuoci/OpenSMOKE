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
#include "preprocessing/OpenSMOKE_PreProcessorIdealGas.h"
#include "preprocessing_chemkin/OpenSMOKE_CHEMKINInterpreter.h"
#include "preprocessing_chemkin/OpenSMOKE_CHEMKINInterpreter_TransportData.h";
#include "preprocessing_chemkin/OpenSMOKE_CHEMKINInterpreter_TransportSpecies.h";
#include <iomanip>

const double OpenSMOKE_PreProcessorIdealGas::BOLTZMANN		= Constants::kBoltzmann*1.e7;		// [erg/K]
const double OpenSMOKE_PreProcessorIdealGas::BOLTZMANN3		= BzzPow3(BOLTZMANN);				// [erg3/K3]
const double OpenSMOKE_PreProcessorIdealGas::MU_MIN			= 1.e-38;							// [TODO]
const double OpenSMOKE_PreProcessorIdealGas::ZROTA			= 0.5 * pow(Constants::pi, 1.5);	
const double OpenSMOKE_PreProcessorIdealGas::ZROTB			= 2.0 + 0.25 * BzzPow2(Constants::pi);
const double OpenSMOKE_PreProcessorIdealGas::ZROTC			= pow(Constants::pi, 1.5);
const double OpenSMOKE_PreProcessorIdealGas::CONST_2SUPI	= 2./Constants::pi;
const double OpenSMOKE_PreProcessorIdealGas::CONST_5SU3R	= 5./(3.*Constants::R_J_mol);

BzzVector TStar(37,
	.1,.2,.3,.4,.5,.6,.7,.8,.9,1.,1.2,1.4,1.6,1.8,2.,2.5,
	3.,3.5,4.,5.,6.,7.,8.,9.,10.,12.,14.,16.,18.,20.,25.,
	30.,35.,40.,50.,75.,100.);

BzzVector deltaStar(8,0.,.25,.5,.75,1.,1.5,2.,2.5);

BzzMatrix O37_8(37,8,
4.008 , 4.002 , 4.655 , 5.52  , 6.454 , 8.214 , 9.824 , 11.31,
3.130 , 3.164 , 3.355 , 3.721 , 4.198 , 5.23  , 6.225 , 7.160,
2.649 , 2.657 , 2.77  , 3.002 , 3.319 , 4.054 , 4.785 , 5.483,
2.314 , 2.32  , 2.402 , 2.572 , 2.812 , 3.386 , 3.972 , 4.539,
2.066 , 2.073 , 2.14  , 2.278 , 2.472 , 2.946 , 3.437 , 3.918,
1.877 , 1.885 , 1.944 , 2.06  , 2.225 , 2.628 , 3.054 , 3.747,
1.729 , 1.738 , 1.79  , 1.893 , 2.036 , 2.388 , 2.763 , 3.137,
1.6122, 1.622 , 1.67  , 1.76  , 1.886 , 2.198 , 2.535 , 2.872,
1.517 , 1.527 , 1.572 , 1.653 , 1.765 , 2.044 , 2.35  , 2.657,
1.44  , 1.45  , 1.49  , 1.564 , 1.665 , 1.917 , 2.196 , 2.4780,
1.3204, 1.33  , 1.364 , 1.425 , 1.51  , 1.72  , 1.956 , 2.199 ,
1.234 , 1.24  , 1.272 , 1.324 , 1.394 , 1.573 , 1.777 , 1.99  ,
1.168 , 1.176 , 1.202 , 1.246 , 1.306 , 1.46  , 1.64  , 1.827 ,
1.1166, 1.124 , 1.146 , 1.185 , 1.237 , 1.372 , 1.53  , 1.7   ,
1.075 , 1.082 , 1.102 , 1.135 , 1.181 , 1.3   , 1.441 , 1.592 ,
1.0006, 1.005 , 1.02  , 1.046 , 1.08  , 1.17  , 1.278 , 1.397 ,
 .95  ,  .9538,  .9656,  .9852, 1.012 , 1.082 , 1.168 , 1.265 ,
 .9131,  .9162,  .9256,  .9413,  .9626, 1.019 , 1.09  , 1.17  ,
 .8845,  .8871,  .8948,  .9076,  .9252,  .972 , 1.03  , 1.098 ,
 .8428,  .8446,  .850 ,  .859 ,  .8716,  .9053,  .9483,  .9984,
 .813 ,  .8142,  .8183,  .825 ,  .8344,  .8598,  .8927,  .9316,
 .7898,  .791 ,  .794 ,  .7993,  .8066,  .8265,  .8526,  .8836,
 .7711,  .772 ,  .7745,  .7788,  .7846,  .8007,  .822 ,  .8474,
 .7555,  .7562,  .7584,  .7619,  .7667,  .78  ,  .7976,  .8189,
 .7422,  .743 ,  .7446,  .7475,  .7515,  .7627,  .7776,  .796 ,
 .72022, .7206,  .722 ,  .7241,  .7271,  .7354,  .7464,  .76  ,
 .7025,  .703 ,  .704 ,  .7055,  .7078,  .7142,  .7228,  .7334,
 .68776, .688,   .6888,  .6901,  .6919,  .697 ,  .704 ,  .7125,
 .6751,  .6753,  .676 ,  .677 ,  .6785,  .6827,  .6884,  .6955,
 .664 ,  .6642,  .6648,  .6657,  .6669,  .6704,  .6752,  .681,
 .6414,  .6415,  .6418,  .6425,  .6433,  .6457,  .649 ,  .653 ,
 .6235,  .6236,  .6239,  .6243,  .6249,  .6267,  .629 ,  .632 ,
 .60882, .6089,  .6091,  .6094,  .61  ,  .6112,  .613 ,  .6154,
 .5964,  .5964,  .5966,  .597 ,  .5972,  .5983,  .600 ,  .6017,
 .5763,  .5763,  .5764,  .5766,  .5768,  .5775,  .5785,  .58  ,
 .5415,  .5415,  .5416,  .5416,  .5418,  .542 ,  .5424,  .543 ,
 .518 ,  .518 ,  .5182,  .5184,  .5184,  .5185,  .5186,  .5187);


BzzMatrix P37_8(37,8,
4.1   , 4.266 , 4.833 , 5.742 , 6.729 , 8.624 ,10.34  ,11.890 ,
3.263 , 3.305 , 3.516 , 3.914 , 4.433 , 5.57  , 6.637 , 7.618 ,
2.84  , 2.836 , 2.936 , 3.168 , 3.511 , 4.329 , 5.126 , 5.874 ,
2.531 , 2.522 , 2.586 , 2.749 , 3.004 , 3.64  , 4.282 , 4.895 ,
2.284 , 2.277 , 2.329 , 2.46  , 2.665 , 3.187 , 3.727 , 4.249 ,
2.084 , 2.081 , 2.13  , 2.243 , 2.417 , 2.862 , 3.329 , 3.786 ,
1.922 , 1.924 , 1.97  , 2.072 , 2.225 , 2.641 , 3.028 , 3.435 ,
1.7902, 1.795 , 1.84  , 1.934 , 2.07  , 2.417 , 2.788 , 3.156 ,
1.682 , 1.689 , 1.733 , 1.82  , 1.944 , 2.258 , 2.596 , 2.933 ,
1.593 , 1.60  , 1.644 , 1.725 , 1.84  , 2.124 , 2.435 , 2.746 ,
1.455 , 1.465 , 1.504 , 1.574 , 1.67  , 1.913 , 2.181 , 2.45  ,
1.355 , 1.365 , 1.4   , 1.461 , 1.544 , 1.754 , 1.989 , 2.228 ,
1.28  , 1.289 , 1.321 , 1.374 , 1.447 , 1.63  , 1.838 , 2.053 ,
1.222 , 1.231 , 1.26  , 1.306 , 1.37  , 1.532 , 1.718 , 1.912 ,
1.176 , 1.184 , 1.209 , 1.25  , 1.307 , 1.45  , 1.618 , 1.795 ,
1.0933, 1.1   , 1.119 , 1.15  , 1.193 , 1.304 , 1.435 , 1.578 ,
1.039 , 1.044 , 1.06  , 1.083 , 1.117 , 1.204 , 1.31  , 1.428 ,
 .9996, 1.004 , 1.016 , 1.035 , 1.062 , 1.133 , 1.22  , 1.32  ,
 .9699,  .9732,  .983 ,  .9991, 1.021 , 1.08  , 1.153 , 1.236 ,
 .9268,  .9291,  .936 ,  .9473,  .9628, 1.005 , 1.058 , 1.12,
.8962,  .8979,  .903 ,  .9114,  .923 ,  .9545,  .9955, 1.044 ,
.8727,  .8741,  .878 ,  .8845,  .8935,  .918 ,  .9505,  .9893,
.8538,  .8549,  .858 ,  .8632,  .8703,  .890 ,  .9164,  .9482,
.8379,  .8388,  .8414,  .8456,  .8515,  .868 ,  .8895,  .916 ,
.8243,  .8251,  .8273,  .8308,  .8356,  .8493,  .8676,  .89  ,
.8018,  .8024,  .8039,  .8065,  .810 ,  .820 ,  .8337,  .8504,
.7836,  .784 ,  .7852,  .7872,  .7899,  .7976,  .808 ,  .8212,
.7683,  .7687,  .7696,  .771 ,  .7733,  .7794,  .788 ,  .7983,
.7552,  .7554,  .7562,  .7575,  .7592,  .764 ,  .771 ,  .7797,
.7436,  .7438,  .7445,  .7455,  .747 ,  .7512,  .757 ,  .7642,
.71982, .72  ,  .7204,  .7211,  .7221,  .725 ,  .7289,  .7339,
.701 ,  .7011,  .7014,  .702 ,  .7026,  .7047,  .7076,  .7112,
.68545, .6855,  .686 ,  .686 ,  .6867,  .6883,  .6905,  .693 ,
.6723,  .6724,  .6726,  .673 ,  .6733,  .6745,  .676 ,  .6784,
.651 ,  .651 ,  .6512,  .6513,  .6516,  .6524,  .6534,  .6546,
.614 ,  .614 ,  .6143,  .6145,  .6147,  .6148,  .6148,  .6147,
.5887,  .5889,  .5894,  .59  ,  .5903,  .5901,  .5895,  .5885);

BzzVector FITASTAR(7, .1106910525E+01, -.7065517161E-02, -.1671975393E-01,  .1188708609E-01, .7569367323E-03, -.1313998345E-02,  .1720853282E-03);
BzzVector FITBSTAR(7, .1199673577E+01, -.1140928763E+00, -.2147636665E-02,  .2512965407E-01, -.3030372973E-02, -.1445009039E-02, .2492954809E-03);
BzzVector FITCSTAR(7, .8386993788E+00,  .4748325276E-01, .3250097527E-01, -.1625859588E-01, -.2260153363E-02,  .1844922811E-02, -.2115417788E-03);
 
void OpenSMOKE_PreProcessorIdealGas::ErrorMessage(const string message)
{
    cout << endl;
    cout << "Class:  OpenSMOKE_PreProcessorIdealGas"	<< endl;
    cout << "Object: " << name_object		<< endl;
    cout << "Error:  " << message           << endl;
    cout << "Press enter to continue... "   << endl;
    getchar();
    exit(-1);
}

void OpenSMOKE_PreProcessorIdealGas::WarningMessage(const string message)
{
    cout << endl;
    cout << "Class:	  OpenSMOKE_PreProcessorIdealGas"	<< endl;
    cout << "Object:  " << name_object		<< endl;
    cout << "Warning: " << message          << endl;
	cout << endl;
}

OpenSMOKE_PreProcessorIdealGas::OpenSMOKE_PreProcessorIdealGas()
{
	TMIN			= 300.;
	TMAX			= 3600.;
	fittingPoints	= 50;

	name_object   = "[not assigned]";
	name_author   = "[not assigned]";
	building_date = GiveMeTimeAndDate();
}

void OpenSMOKE_PreProcessorIdealGas::SetMinimumTemperature(const double tmin)
{
	TMIN		= tmin;
	if (TMIN >= TMAX || TMIN <= 0.)
		ErrorMessage("Please check your minimum temperature!");
}

void OpenSMOKE_PreProcessorIdealGas::SetMaximumTemperature(const double tmax)
{
	TMAX		= tmax;
	if (TMIN >= TMAX || TMAX >= 6000.)
		ErrorMessage("Please check your maximum temperature!");
}

void OpenSMOKE_PreProcessorIdealGas::SetName(const string name)
{
	name_object = name;
}

void OpenSMOKE_PreProcessorIdealGas::SetAuthorName(const string name)
{
	name_author = name;
}

void OpenSMOKE_PreProcessorIdealGas::WriteHeaderFile(BzzSave &binaryFile, BzzSave &asciiFile)
{
	string dummy;
	char name[Constants::NAME_SIZE];

	dummy = "V101116";
	strcpy(name, dummy.c_str());
	binaryFile.fileSave.write((char*) name, sizeof(name));
	asciiFile << name;
	dummy = name_object;
	strcpy(name, dummy.c_str());
	binaryFile.fileSave.write((char*) name, sizeof(name));
	asciiFile << "ASCIINameObject";
	dummy = name_author;
	strcpy(name, dummy.c_str());
	binaryFile.fileSave.write((char*) name, sizeof(name));
	asciiFile << "ASCIINameAuthor";
	dummy = name_author;
	strcpy(name, dummy.c_str());
	binaryFile.fileSave.write((char*) name, sizeof(name));
	asciiFile << "ASCIINameAuthor";
	dummy = building_date;
	strcpy(name, dummy.c_str());
	binaryFile.fileSave.write((char*) name, sizeof(name));
	asciiFile << "ASCIINameDate";

	binaryFile << TMIN;
	asciiFile  << TMIN;
	binaryFile << TMAX;
	asciiFile  << TMAX;
}

void OpenSMOKE_PreProcessorIdealGas::Setup(const string fileNames, const string fileTransport, const string fileThermo, const string fileElements, BzzSave &binaryFile, BzzSave &asciiFile)
{
	WriteHeaderFile(binaryFile, asciiFile);
	
	// Lettura dei nomi delle specie e del numero di componenti
	readNames(fileNames, binaryFile, asciiFile);

	// Allocazione della memoria per le variabili
	allocate();

	// Lettura delle informazioni termodinamiche
	readThermodynamics(fileThermo, binaryFile, asciiFile);

	// Lettura delle proprieta di trasporto
	readTransportProperties(fileTransport);

	// Lettura degli elementi
	ReadElements(fileElements, binaryFile, asciiFile);

	// Calcolo delle grandezze indipendenti dalla temperatura, pressione e composizione
	Initialize();
	InitializeTransportProperties();
}

void OpenSMOKE_PreProcessorIdealGas::SetupTransportSensitivity(const string pathName, const string fileNames, const string fileTransport, const string fileThermo, const string fileElements, BzzSave &binaryFile, BzzSave &asciiFile,const double eps)
{
	WriteHeaderFile(binaryFile, asciiFile);
	
	// Lettura dei nomi delle specie e del numero di componenti
	readNames(fileNames, binaryFile, asciiFile);

	// Allocazione della memoria per le variabili
	allocate();

	// Lettura delle informazioni termodinamiche
	readThermodynamics(fileThermo, binaryFile, asciiFile);

	// Lettura delle proprieta di trasporto
	readTransportProperties(fileTransport);

	// Lettura degli elementi
	ReadElements(fileElements, binaryFile, asciiFile);

	// Calcolo delle grandezze indipendenti dalla temperatura, pressione e composizione
	BzzSave fOutput;
	fOutput(pathName + "/TransportSensitivity.out");
	{
		cout << " ** TransportSensitivity.out file..." << endl;

		// Number of species
		fOutput << "NC";
		fOutput << NC;

		// Number of species
		fOutput << "EPSILON";
		fOutput << eps;

		// Species
		for(int i=1;i<=NC;i++)
			fOutput << names[i];
		
		// Parameters
		fOutput << epsylon_over_kb;
		fOutput << sigma;
		fOutput << mu;
		fOutput << alfa;
		fOutput << zRot298;


		// Base case
		fOutput << "BASE";
		Initialize();
		InitializeTransportProperties();
		Fitting(fittingPoints);
		WriteFittingTransportSensitivity(fOutput);

		// eps/k 1+eps
		fOutput << "EPSKPLUS";
		epsylon_over_kb *= (1.+eps);
		Initialize();
		InitializeTransportProperties();
		Fitting(fittingPoints);
		WriteFittingTransportSensitivity(fOutput);
		epsylon_over_kb /= (1.+eps);

		// eps/k 1-eps
		fOutput << "EPSKMINUS";
		epsylon_over_kb *= (1.-eps);
		Initialize();
		InitializeTransportProperties();
		Fitting(fittingPoints);
		WriteFittingTransportSensitivity(fOutput);
		epsylon_over_kb /= (1.-eps);

		// Sigma 1+eps
		fOutput << "SIGMAPLUS";
		sigma *= (1.+eps);
		Initialize();
		InitializeTransportProperties();
		Fitting(fittingPoints);
		WriteFittingTransportSensitivity(fOutput);
		sigma /= (1.+eps);

		// Sigma 1-eps
		fOutput << "SIGMAMINUS";
		sigma *= (1.-eps);
		Initialize();
		InitializeTransportProperties();
		Fitting(fittingPoints);
		WriteFittingTransportSensitivity(fOutput);
		sigma /= (1.-eps);

		// Mu 1+eps
		fOutput << "MUPLUS";
		mu *= (1.+eps);
		Initialize();
		InitializeTransportProperties();
		Fitting(fittingPoints);
		WriteFittingTransportSensitivity(fOutput);
		mu /= (1.+eps);

		// Mu 1-eps
		fOutput << "MUMINUS";
		mu *= (1.-eps);
		Initialize();
		InitializeTransportProperties();
		Fitting(fittingPoints);
		WriteFittingTransportSensitivity(fOutput);
		mu /= (1.-eps);

		// Alfa 1+eps
		fOutput << "ALFAPLUS";
		alfa *= (1.+eps);
		Initialize();
		InitializeTransportProperties();
		Fitting(fittingPoints);
		WriteFittingTransportSensitivity(fOutput);
		alfa /= (1.+eps);

		// Alfa 1-eps
		fOutput << "ALFAMINUS";
		alfa *= (1.-eps);
		Initialize();
		InitializeTransportProperties();
		Fitting(fittingPoints);
		WriteFittingTransportSensitivity(fOutput);
		alfa /= (1.-eps);

		// zRot 1+eps
		fOutput << "ZROTPLUS";
		zRot298 *= (1.+eps);
		Initialize();
		InitializeTransportProperties();
		Fitting(fittingPoints);
		WriteFittingTransportSensitivity(fOutput);
		zRot298 /= (1.+eps);

		// Alfa 1-eps
		fOutput << "ZROTMINUS";
		zRot298 *= (1.-eps);
		Initialize();
		InitializeTransportProperties();
		Fitting(fittingPoints);
		WriteFittingTransportSensitivity(fOutput);
		zRot298 /= (1.-eps);

		fOutput.End();
	}
}

// ************************************************************************************	//
// ************************************************************************************ //
//										INIZIALIZZAZIONI								//
// ************************************************************************************	//
// ************************************************************************************ //
//----------------------------------------------------------------------------------//
//										Initialize									//
//																					//
//		- allocazione della memoria e inizializzazione delle variabili NON			//
//		  dipendenti dalla temperatura												//
//																					//
//----------------------------------------------------------------------------------//
void OpenSMOKE_PreProcessorIdealGas::Initialize(void)
{
	int i;

	// Calcolo del potenziale di Lennard-Jones [erg]
	epsylon = OpenSMOKE_PreProcessorIdealGas::BOLTZMANN * epsylon_over_kb;
	for(i = 1;i <= NC;i++)
		KbSuepsylon[i] = 1./epsylon_over_kb[i];

	// Calcolo dei reciproci dei pesi molecolari [mol/g]
	for(i = 1;i <= NC;i++)
		uM[i] = 1. / M[i];							// reciproco del peso molecolare di ciascun componente

	// Controllo dei momenti dipolari
	for(i = 1;i <= NC;i++)
		if(mu[i] < OpenSMOKE_PreProcessorIdealGas::MU_MIN)  mu[i] = 0.;  // momento dipolare [sqrt(A3.erg)]


	// Momento dipolare ridotto per la molecola polare F.(13)
	for(i = 1;i <= NC;i++)
		muStar[i] = mu[i] / sqrt(epsylon[i] * BzzPow3(sigma[i])); //[-]
}

void OpenSMOKE_PreProcessorIdealGas::InitializeTransportProperties(void)
{
	cout << "       processing diffusion coefficients..." << endl;
	allocateDiffusivity();
	initializeDiffusivity();

	cout << "       processing viscosity coefficients..." << endl;
	allocateViscosity();
	initializeViscosity();

	cout << "       processing conductivity coefficients..." << endl;
	allocateConductivity();
	initializeConductivity();

	cout << "       processing thermal diffusion ratios..." << endl;
	allocateThermalDiffusionRatios();
}

void OpenSMOKE_PreProcessorIdealGas::initializeDiffusivity(void)
{
	int j,k;

	BzzVector alfaStar(NC);
	BzzMatrix mujkStar(NC, NC),
					epsylonjk(NC, NC),
					sigmajk(NC, NC),
					Mjk(NC, NC),
					mujk(NC, NC);


	// Calcolo della massa molecolare ridotta per tutte le coppie di specie
	// attenzione non ?la grandezza che corrisponde a quella del CHEMKIN
	for(j = 1;j <= NC;j++)
		for(k = 1;k <= NC;k++)
			Mjk[j][k] = 2. * M[j] * M[k]  / (M[j] + M[k]); // [kg/kmol]

	// Calcolo  dei potenziali, diametri e momenti effettivi di collisione

	// CASO 1: entrambe le sostanze polari o non polari

		// Calcolo dei potenziali effettivi di Lennard-Jones gi?divisi
		// per la costante di Boltzmann F.(5)
		for(j = 1;j <= NC;j++)
			for(k = 1;k <= NC;k++)
			{
				epsylonjkSuKb[j][k] = sqrt(epsylon_over_kb[j] * epsylon_over_kb[k]); //[K]
				KbSuepsylonjk[j][k] = 1./epsylonjkSuKb[j][k];				 //[1/K]
			}

		// Calcolo dei potenziali effettivi di Lennard-Jones
		// non divisi per la costante di Boltzmann F.(5)
		for(j = 1;j <= NC;j++)
			for(k = 1;k <= NC;k++)
				epsylonjk[j][k] = epsylonjkSuKb[j][k] * BOLTZMANN; // [erg]

		// Diametri effettivi di collisione di Lennard-Jones F.(6)
		for(j = 1;j <= NC;j++)
			for(k = 1;k <= NC;k++)
				sigmajk[j][k] = .5 * (sigma[j] + sigma[k]);	// [A]

		// Momento dipolare effettivo di collisione F.(7)
		for(j = 1;j <= NC;j++)
			for(k = 1;k <= NC;k++)
				mujk[j][k] = sqrt(mu[j] * mu[k]);	// [sqrt(A3.erg)]

	// CASO 2: una molecola polare e una non polare

	double xsi; // coefficiente correttivo
	for(j = 1;j <= NC;j++)
	{
		for(k = 1;k <= NC;k++)
		{
			if(mu[j] == 0. && mu[k] != 0.)
			{
				// a. polarizzabilit?ridotta per la molecola non polare F.(12)
				//    alfa misurato in A3 - sigma misurato in A
				//    non ?necessario adottare nessuna conversione
				alfaStar[j] = alfa[j] / BzzPow3(sigma[j]);	//[-]

				// c. Coefficiente correttivo F.(11)
				xsi = 1. + .25 * alfaStar[j] * BzzPow2(muStar[k]) * sqrt(epsylon[k] / epsylon[j]); //[-]

				// d. Correzioni
				epsylonjkSuKb[j][k] *= (xsi * xsi);		// F.(8)
				sigmajk[j][k] *= pow(xsi, -1./6.);		// F.(9)
			}

			else if(mu[k] == 0. && mu[j] != 0.)
			{

				// a. polarizzabilit?ridotta per la molecola non polare F.(12)
				alfaStar[k] = alfa[k] / BzzPow3(sigma[k]);	//[-]

				// c. Coefficiente correttivo F.(11)
				xsi = 1. + .25 * alfaStar[k] * BzzPow2(muStar[j]) * sqrt(epsylon[j] / epsylon[k]); //[-]

				// d. Correzioni
				epsylonjkSuKb[j][k] *= (xsi * xsi);		// F.(8)
				sigmajk[j][k] *= pow(xsi,-1./6.);		// F.(9)
			}
		}
	}

	// Momento dipolare effettivo di collisione
	// Vedi formula F.(2)
	for(j = 1;j <= NC;j++)
		for(k = 1;k <= NC;k++)
			mujkStar[j][k] = mujk[j][k]  / sqrt(epsylonjk[j][k] * BzzPow3(sigmajk[j][k])); // [-]


	// Calcolo del momento dipolare ridotto F.(15)
	for(j = 1;j <= NC;j++)
		for(k = 1;k <= NC;k++)
			deltajkStar[j][k] = .5 * BzzPow2(mujkStar[j][k]);


	// Coefficiente costante nella formula per il calcolo dei coefficienti di
	// diffusione binari
	for(k = 1;k <= NC;k++)
		for(j = 1;j <= NC;j++)
			coeff_Djk[j][k] = 2.6693e-7  / ( sqrt(Mjk[j][k])*BzzPow2(sigmajk[j][k]) );
}

void OpenSMOKE_PreProcessorIdealGas::initializeConductivity(void)
{
	int k;
	double aux;

	for(k = 1;k <= NC;k++)
	{
		// Funzione di Parker-Brau-Jonkman a 298K F.(33)
		aux = sqrt(epsylon_over_kb[k] / 298.); //[-]
		f298[k] = 1. + aux * ( ZROTA + aux * (ZROTB + aux * ZROTC)); // [-]
	}

	// Calcolo dei contributi traslazionali e rotazionali per il calore
	// specifico molare a volume costante [J/molK]
	cVtrans = 1.5 * Constants::R_J_mol; // F.(22) F.(25) F.(28)
	for(k = 1;k <= NC;k++)
	{
		if(shape_factor[k] == 0)				// specie monoatomica
			cVrot[k] = 0.;				// F.(?)

		else if(shape_factor[k] == 1)				// specie a molecola lineare
			cVrot[k] = Constants::R_J_mol;	// F.(23)

		else							// specie a molecola non lineare
			cVrot[k] = cVtrans[k];		// F.(26)
	}

	// Coefficiente costante nella formula per il calcolo dei coefficienti di
	// diffusione binari
	for(k = 1;k <= NC;k++)
			coeff_Dkk[k] = 2.6693e-7  / ( sqrt(M[k])*BzzPow2(sigma[k]) );
}

void OpenSMOKE_PreProcessorIdealGas::initializeViscosity(void)
{
	// Coefficiente costante nella formula per il calcolo della viscosit?dei singoli
	// componenti
	for(int k = 1;k <= NC;k++)
		coeff_eta[k] = 26.693e-7 * sqrt(M[k]) / BzzPow2(sigma[k]);
}

// ************************************************************************************	//
// ************************************************************************************ //
//										FUNZIONI INTERNE								//
// ************************************************************************************	//
// ************************************************************************************ //

//----------------------------------------------------------------------------------//
//									TemperaturaRidotta								//
//																					//
//		- calcolo delle temperature ridotte effettive								//
//																					//
//----------------------------------------------------------------------------------//

void OpenSMOKE_PreProcessorIdealGas::TemperaturaRidotta(double T)
{
	// - si tratta di una matrice quadrata simmetrica che viene utilizzata soltanto
	//   come coordinata di ingresso per il calcolo dell' integrale collisionale
	//   di diffusivit?binaria omega11
	// - viene riempita soltanto una met? la rimanente, per ragioni di simmetria,
	//   viene posta uguale alla prima (in realt?viene utilizzata soltanto la prima met?
	//   per il calcolo dell'integrale collisionale)

	int j,k;
	for(j = 1;j <= NC;j++)
		for(k = j;k <= NC;k++)
			TjkStar[j][k] = T * KbSuepsylonjk[j][k];	// F.(14)
}

//----------------------------------------------------------------------------------//
//							TemperaturaRidottaComponentePuro						//
//																					//
//		- calcolo delle temperature ridotte											//
//																					//
//----------------------------------------------------------------------------------//

void OpenSMOKE_PreProcessorIdealGas::TemperaturaRidottaComponentePuro(double T)
{
	// - si tratta di un vettore che viene utilizzato soltanto come coordinata di
	//   ingresso per il calcolo dell'integrale collisionale di viscosit?omega22

	for(int k = 1;k <= NC;k++)
		TkStar[k] = T * KbSuepsylon[k];									// F.(2)
}

//----------------------------------------------------------------------------------//
//										Omega11jk									//
//																					//
//		- calcolo dell'integrale collisionale per le diffusivit?				//
//																					//
//----------------------------------------------------------------------------------//

void OpenSMOKE_PreProcessorIdealGas::Omega11jk(void)
{
	// - si tratta di una matrice simmetrica
	// - viene riempita soltanto la met?superiore + la diagonale principale
	// - la met?inferiore dovrebbe essere posta uguale a quella superiore, ma ?una
	//   operazione inutile, tanto non verr?mai utilizzata
	// - questa matrice viene utilizzata soltanto per il calcolo dei coefficienti di
	//   diffusione binari Djk

	// - parametri di ingresso:
	//   TjkStar: matrice simmetrica [-]
	//   deltajkStar: matrice simmetrica [-]

	// - il calcolo viene eseguito facendo una interpolazione bidimensionale sulla
	//   matrice O37_8


	for(int j = 1;j <= NC;j++)
		for(int k = j+1;k <= NC;k++)
			omega11jk[j][k] = collisionIntegral_11(TjkStar[j][k],deltajkStar[j][k]);
}

void OpenSMOKE_PreProcessorIdealGas::Omega11kk(void)
{
	for(int k = 1;k <= NC;k++)
		omega11kk[k] = collisionIntegral_11(TkStar[k],deltakStar[k]);
}

double OpenSMOKE_PreProcessorIdealGas::collisionIntegral_11(double tjk, double djk)
{
	int itStar,idStar;

	double	udx21,udx321,
			dx31,dx32,
			dxx1,dxx2,
			x1,x2,x3,
			y1,y2,y3,
			a1,a2,a3;

	if ( djk < -0.00001 )
		cout << "WARNING: Diffusivity collision integral undefined (1)" << endl;

	if ( djk > 2.5 )
		cout << "WARNING: Diffusivity collision integral undefined (2)" << endl;

	if ( tjk < 0.09 )
	{
		//cout << "Diffusivity collision integral undefined (3)" << endl;
		tjk = 0.09;
	}

	if ( tjk > 500.)
		cout << "WARNING: Diffusivity collision integral undefined (4)" << endl;

	if ( fabs(djk)>1.e-5 && tjk>75.)
		cout << "WARNING: Diffusivity collision integral undefined (5)" << endl;


	// 1. Soluzione nel caso in cui TjkStar cada al di l?dell'ultimo valore tabulato
	if(tjk > TStar[36])
		return ( .623 + tjk * (-.136e-2 + tjk * (.346e-5 - tjk * .343e-8)) );

	// 2. Soluzione nel caso in cui deltajkStar sia pi piccolo dei valori tabulati
	if(fabs(djk) <= 1.e-5)
	{
		// fa l'interpolazione usando solo la prima colonna
		if(tjk <  TStar[2])
			itStar = 1;
		else
			itStar = LocateInSortedVector(TStar,tjk);

		// poi fa l'interpolazione
		x1 = TStar[itStar];
		x2 = TStar[itStar + 1];
		x3 = TStar[itStar + 2];
		udx21 = 1. / (x2 - x1);
		dx31 = x3 - x1;
		dx32 = x3 - x2;
		dxx1 = tjk - x1;
		dxx2 = tjk - x2;
		udx321 = 1. / (dx31 * dx32);

		// prima interpolazione rispetto a TStar
		a1 = O37_8[itStar][1];
		a2 = (O37_8[itStar + 1][1] - a1) * udx21;
		a3 = (O37_8[itStar + 2][1] - a1 - a2 * dx31) * udx321;

		return ( a1 + dxx1 * ( a2 + a3 * dxx2) );
	}

	// 3. Caso generale
	// deve individuare la posizione di *tjk sul vettore TStar
		if(tjk <  TStar[2])
			itStar = 1;
		else
			itStar = LocateInSortedVector(TStar,tjk);

		// deve individuare la posizione di *djk sul vettore deltaStar
		if(djk < deltaStar[2])
			idStar = 1;
		else if(djk > deltaStar[7])
			idStar = 6;
		else
			idStar = LocateInSortedVector(deltaStar,djk);

		// poi fa l'interpolazione
		x1 = TStar[itStar];
		x2 = TStar[itStar + 1];
		x3 = TStar[itStar + 2];
		udx21 = 1. / (x2 - x1);
		dx31 = x3 - x1;
		dx32 = x3 - x2;
		dxx1 = tjk - x1;
		dxx2 = tjk - x2;
		udx321 = 1. / (dx31 * dx32);

		// prima interpolazione rispetto a TStar
		a1 = O37_8[itStar][idStar];
		a2 = (O37_8[itStar + 1][idStar] - a1) * udx21;
		a3 = (O37_8[itStar + 2][idStar] - a1 - a2 * dx31) * udx321;
		y1 = a1 + dxx1 * ( a2 + a3 * dxx2);

		// seconda interpolazione rispetto a TStar
		a1 = O37_8[itStar][idStar + 1];
		a2 = (O37_8[itStar + 1][idStar + 1] - a1) * udx21;
		a3 = (O37_8[itStar + 2][idStar + 1] - a1 - a2 * dx31) * udx321;
		y2 = a1 + dxx1 * ( a2 + a3 * dxx2);

		// terza interpolazione rispetto a TStar
		a1 = O37_8[itStar][idStar + 2];
		a2 = (O37_8[itStar + 1][idStar + 2] - a1) * udx21;
		a3 = (O37_8[itStar + 2][idStar + 2] - a1 - a2 * dx31) * udx321;
		y3 = a1 + dxx1 * ( a2 + a3 * dxx2);

		// interpolazione rispetto a deltaStar
		x1 = deltaStar[idStar];
		x2 = deltaStar[idStar + 1];
		x3 = deltaStar[idStar + 2];
		udx21 = 1. / (x2 - x1);
		dx31 = x3 - x1;
		dx32 = x3 - x2;
		dxx1 = djk - x1;
		dxx2 = djk - x2;
		udx321 = 1. / (dx31 * dx32);
		a1 = y1;
		a2 = (y2 - a1) * udx21;
		a3 = (y3 - a1 - a2 * dx31) * udx321;

		return ( a1 + dxx1 * ( a2 + a3 * dxx2));

}

//----------------------------------------------------------------------------------//
//										Omega22k									//
//																					//
//		- calcolo dell'integrale collisionale per le viscosit?					//
//																					//
//----------------------------------------------------------------------------------//

void OpenSMOKE_PreProcessorIdealGas::Omega22k(void)
{
	// - ?un vettore che viene utilizzato soltanto per il calcolo della viscosit?
	//   dei singoli componenti

	// - parametri in ingresso
	//   TkStar: vettore [-]
	//   deltakStar: vettore [-]

	int itStar,idStar;

	int k;
	double udx21,dx31,dx32,udx321,dxx1,dxx2,x1,x2,x3,y1,y2,y3,a1,a2,a3;
	double *tk = TkStar.GetHandle();
	double *dk = deltakStar.GetHandle();
	double *ok = omega22k.GetHandle();

	for(k = 1;k <= NC;k++)
	{
		if( *dk < -0.00001 )
			cout << "Viscosity-Conductivity collision integral undefined (1) - Species: " << names[k] << endl;

		if( *dk > 2.5 )
			cout << "Viscosity-Conductivity collision integral undefined (2) - Species: " << names[k] << endl;

		if( *tk < 0.09 )
		{
			cout << "Viscosity-Conductivity collision integral undefined (3) - Species: " << names[k] << endl;
			*tk = 1.00; // e' una forzatura
		}
		if( *tk > 500. )
			cout << "Viscosity-Conductivity collision integral undefined (4) - Species: " << names[k] << endl;

		if ( fabs(*dk) > 1.e-5 && *tk > 75.)
			cout << "Viscosity-Conductivity collision integral undefined (5) - Species: " << names[k] << endl;

		if(*tk > TStar[36])
		{
			*ok = .703 + *tk * (-.146e-2 + *tk * (.357e-5 - *tk * .343e-8));
			tk++;
			dk++;
			ok++;
			continue;
		}
		if(fabs(*dk) <= 1.e-5)
		{
			// fa l'interpolazione usando solo la prima colonna
			if(*tk <  TStar[2])
				itStar = 1;
			else
				itStar = LocateInSortedVector(TStar,*tk);
			// poi fa l'interpolazione
			x1 = TStar[itStar];
			x2 = TStar[itStar + 1];
			x3 = TStar[itStar + 2];
			udx21 = 1. / (x2 - x1);
			dx31 = x3 - x1;
			dx32 = x3 - x2;
			dxx1 = *tk - x1;
			dxx2 = *tk - x2;
			udx321 = 1. / (dx31 * dx32);
			// prima interpolazione rispetto a TStar
			a1 = P37_8[itStar][1];
			a2 = (P37_8[itStar + 1][1] - a1) * udx21;
			a3 = (P37_8[itStar + 2][1] - a1 - a2 * dx31) * udx321;
			*ok = a1 + dxx1 * ( a2 + a3 * dxx2);
			tk++;
			dk++;
			ok++;
			continue;
		}

		// deve individuare la posizione di *tjk sul vettore TStar
		if(*tk <  TStar[2])
			itStar = 1;
		else
			itStar = LocateInSortedVector(TStar,*tk);
		// deve individuare la posizione di *djk sul vettore deltaStar
		if(*dk < deltaStar[2])
			idStar = 1;
		else if(*dk > deltaStar[7])
			idStar = 6;
		else
			idStar = LocateInSortedVector(deltaStar,*dk);
		// poi fa l'interpolazione
		x1 = TStar[itStar];
		x2 = TStar[itStar + 1];
		x3 = TStar[itStar + 2];
		udx21 = 1. / (x2 - x1);
		dx31 = x3 - x1;
		dx32 = x3 - x2;
		dxx1 = *tk - x1;
		dxx2 = *tk - x2;
		udx321 = 1. / (dx31 * dx32);
		// prima interpolazione rispetto a TStar
		a1 = P37_8[itStar][idStar];
		a2 = (P37_8[itStar + 1][idStar] - a1) * udx21;
		a3 = (P37_8[itStar + 2][idStar] - a1 - a2 * dx31) * udx321;
		y1 = a1 + dxx1 * ( a2 + a3 * dxx2);
		// seconda interpolazione rispetto a TStar
		a1 = P37_8[itStar][idStar + 1];
		a2 = (P37_8[itStar + 1][idStar + 1] - a1) * udx21;
		a3 = (P37_8[itStar + 2][idStar + 1] - a1 - a2 * dx31) * udx321;
		y2 = a1 + dxx1 * ( a2 + a3 * dxx2);
		// terza interpolazione rispetto a TStar
		a1 = P37_8[itStar][idStar + 2];
		a2 = (P37_8[itStar + 1][idStar + 2] - a1) * udx21;
		a3 = (P37_8[itStar + 2][idStar + 2] - a1 - a2 * dx31) * udx321;
		y3 = a1 + dxx1 * ( a2 + a3 * dxx2);
		// interpolazione rispetto a deltaStar
		x1 = deltaStar[idStar];
		x2 = deltaStar[idStar + 1];
		x3 = deltaStar[idStar + 2];
		udx21 = 1. / (x2 - x1);
		dx31 = x3 - x1;
		dx32 = x3 - x2;
		dxx1 = *dk - x1;
		dxx2 = *dk - x2;
		udx321 = 1. / (dx31 * dx32);
		a1 = y1;
		a2 = (y2 - a1) * udx21;
		a3 = (y3 - a1 - a2 * dx31) * udx321;
		*ok = a1 + dxx1 * ( a2 + a3 * dxx2);
		tk++;
		dk++;
		ok++;
	}
}

//----------------------------------------------------------------------------------//
//								CoefficienteDiDiffusioneBinario						//
//																					//
//		- calcolo della matrice delle diffusivit?binarie							//
//																					//
//----------------------------------------------------------------------------------//

void OpenSMOKE_PreProcessorIdealGas::CoefficienteDiDiffusioneBinario(double T, double P)
{
	// - ?una matrice simmetrica; viene calcolata soltanto la met?superiore,
	//   quella inferiore viene invece ricopiata
	// - la diagonale principale non viene calcolata qui; corrisponde ai coefficienti
	//   di autodiffusione, calcolati con un'altra function
	// - viene utilizzata soltanto per il calcolo delle diffusivit?in miscela

	// T = Temperatura [K]
	// P = Pressione [bar]
	// omega11jk = integrale collisionale [-]
	// Djk = coefficienti di diffusivit?binari [m2/s]

	int j,k;

	double T_P = pow(T,1.5) / P;

	for(j = 1;j <= NC;j++)
		for(k = j+1; k <= NC; k++)
			Djk[k][j] = Djk[j][k] = coeff_Djk[j][k] * T_P / omega11jk[j][k];
}

//----------------------------------------------------------------------------------//
//								CoefficienteDiAutoDiffusione						//
//																					//
//		- calcolo delle autodiffusivit?											//
//																					//
//----------------------------------------------------------------------------------//

void OpenSMOKE_PreProcessorIdealGas::CoefficienteDiAutoDiffusione(double T)
{
	// - ?la diagonale principale della matrice Djk
	// - viene utilizzata soltanto per calcolare la conducibilit?termica delle singole
	//   specie; non ?necessaria per il calcolo delle diffusivit?in miscela

	// - a rigore la T dovrebbe essere divisa per la pressione in bar; tuttavia il
	//   coefficiente di autodiffusione verr?poi moltiplicato per la densit?delle
	//   singole specie che prevedono la presenza della pressione a numeratore; per
	//   questo motivo si evita di usare la pressione in entrambi i casi

	int k;
	double T_P = pow(T,1.5);
	for(k = 1;k <= NC;k++)
		Dkk[k] = coeff_Dkk[k] * T_P / omega11kk[k];
}

//----------------------------------------------------------------------------------//
//								ViscositaDelComponente								//
//																					//
//		- calcolo delle viscosit?dei singoli componenti							//
//																					//
//----------------------------------------------------------------------------------//

void OpenSMOKE_PreProcessorIdealGas::ViscositaDelComponente(double T)
{
	// Calcolo della viscosit?del singolo componente secondo l'espressione
	// data dalla teoria cinetica standard dei gas F.(1) --> Kuo
	// T = temperatura [K]
	// omega22 = integrale collisionale [-]
	// eta = viscosit?dinamica [Pa.s]

	double sqrT = sqrt(T);
	for(int k = 1;k <= NC;k++)
		eta[k] = coeff_eta[k] * sqrT / omega22k[k];
}

//----------------------------------------------------------------------------------//
//									computeFT										//
//																					//
//		- calcolo della funzione di Parker-Brau-Jonkman F.(33)						//
//																					//
//----------------------------------------------------------------------------------//
void OpenSMOKE_PreProcessorIdealGas::computeFT(double T)
{
	// - viene utilizzata soltanto per il calcolo della conducibilit?termica delle
	//   singole specie

	double uT = 1./T;
	double aux;
	for(int k = 1;k <= NC;k++)
	{
		aux = sqrt(epsylon_over_kb[k] * uT);
		fT[k] = 1. + aux * ( ZROTA + aux * ( ZROTB + aux * ZROTC));
	}
}

//----------------------------------------------------------------------------------//
//									computeZrot										//
//																					//
//		- Calcolo del numero di collisione di rilassamento rotazionale F.(32)		//
//																					//
//----------------------------------------------------------------------------------//
void OpenSMOKE_PreProcessorIdealGas::computeZrot(void)
{
	// - viene utilizzato soltanto per il calcolo della conducibilit?termica
	//   delle singole specie

	for(int k = 1;k <= NC;k++)
		zRot[k] = zRot298[k] * f298[k] / fT[k];
}
//----------------------------------------------------------------------------------//
//									computeCvVib									//
//																					//
//		- Calcolo del contributo vibrazionale al CpV [J/molK]						//
//																					//
//----------------------------------------------------------------------------------//
void OpenSMOKE_PreProcessorIdealGas::computeCvVib(void)
{
	// - viene utilizzato soltanto per il calcolo della conducibilit?termica delle
	//   singole specie

	for(int k = 1;k <= NC;k++)
		if(shape_factor[k] != 0)
			cVvib[k] = Cv[k]*M[k]*1.e-3 - cVtrans[k] - cVrot[k];
}

//----------------------------------------------------------------------------------//
//						SpeciesDensity(double T, double P)							//
//																					//
//		- Calcolo delle densit?per le singole specie [kg/m3]						//
//																					//
//----------------------------------------------------------------------------------//
void OpenSMOKE_PreProcessorIdealGas::SpeciesDensity(double T)
{
	// - viene utilizzata soltanto per il calcolo della conducibilit?termica delle
	//   singole specie

	// - a rigore bisognerebbe moltiplicare per la pressione in bar; tuttavia il
	//   coefficiente di autodiffusione verr?poi moltiplicato per la densit?delle
	//   singole specie che prevedono la presenza della pressione a numeratore; per
	//   questo motivo si evita di usare la pressione in entrambi i casi

	// T = temperatura [K]
	// P = pressione [bar]

	double coeff= 100. / (Constants::R_J_mol * T);

	for(int k = 1;k <= NC;k++)
		rho[k] =  coeff * M[k];
}

//----------------------------------------------------------------------------------//
//							ConducibilitaComponente									//
//																					//
//		- Calcolo della conducibilit?termica per le singole specie [W/mK]			//
//																					//
//----------------------------------------------------------------------------------//
void OpenSMOKE_PreProcessorIdealGas::ConducibilitaComponente(void)
{
	// - ?necessario conoscere eta
	// - ?necessario aver calcolato prima: fT, zRot, rho, CvVib e Dkk

	int k;
	double A, B;

	for(k = 1;k <= NC;k++)
	{
		fVib[k] = rho[k] * Dkk[k] / eta[k];	// F.(19)

		A = 2.5-fVib[k];					// F.(20)
		B = zRot[k] + CONST_2SUPI*(CONST_5SU3R*cVrot[k] + fVib[k]);	// F.(21)
		fTrans[k] = 2.50*(1.-CONST_2SUPI*cVrot[k]/cVtrans[k]*A/B);		// F.(17)

		if(shape_factor[k] == 0)
			lambda[k] = eta[k] * uM[k] * (fTrans[k]*cVtrans[k]) ;	// F.(29)

		else
		{
			fRot[k] = fVib[k]*(1.+CONST_2SUPI*A/B);		// F.(18)
		double rdf=(fTrans[k]*cVtrans[k] + fRot[k]*cVrot[k] + fVib[k]*cVvib[k]);	// F.(16)

			lambda[k] = eta[k] * uM[k] * (fTrans[k]*cVtrans[k] + fRot[k]*cVrot[k] + fVib[k]*cVvib[k]);	// F.(16)
		}
	}

	lambda*=1000.;		// [W/mK]
}

// ************************************************************************************	//
// ************************************************************************************ //
//								PROPRIETA' SINGOLE SPECIE								//
// ************************************************************************************	//
// ************************************************************************************ //

void OpenSMOKE_PreProcessorIdealGas::SpeciesDiffusionCoefficient(double T, double P)
{
	// Diffusivit?in miscela [m2/s]
	// T = temperatura [K]
	// P = pressione [bar]
	// - temperature effettive ridotte
	// - momento dipolare effettivo ridotto
	// - integrale collisionale effettivo di diffusivit?

	// Diffusivit?e viscosit?sono del tutto indipendenti

	TemperaturaRidotta(T);
	Omega11jk();
	CoefficienteDiDiffusioneBinario(T,P);
}

void OpenSMOKE_PreProcessorIdealGas::SpeciesViscosity(double T)
{
	// Viscosit?[Pa.s]
	// T = temperatura [K]

	// - temperature ridotte per i singoli componenti
	// - integrale collisionale di viscosit?per i singoli componenti

	TemperaturaRidottaComponentePuro(T);
	Omega22k();
	ViscositaDelComponente(T);
}

void OpenSMOKE_PreProcessorIdealGas::SpeciesConductivity(double T)
{
	// Conducibilit?termica [W/mK]
	// T = temperatura [K]

	// - Calori specifici a volume costante [J/kgK]
	// - funzione di Jonkman alla temperatura T [-]
	// - zRot alla temperatura T [-]
	// - densit?dei singoli componenti [kg/m3]
	// - viscosit?dei singoli componenti [Pa.s]

	computeFT(T);
	computeZrot();
	computeCvVib();
	SpeciesDensity(T);
	Omega11kk();
	CoefficienteDiAutoDiffusione(T);
	ConducibilitaComponente();
}

// ************************************************************************************	//
// ************************************************************************************ //
//								LETTURA INFORMAZIONI DA FILE							//
// ************************************************************************************	//
// ************************************************************************************ //

void OpenSMOKE_PreProcessorIdealGas::readNames(const string fileName, BzzSave &binaryFile, BzzSave &asciiFile)
{
	cout << "Reading names ... ";

	int k,j;
	char name[Constants::NAME_SIZE];

	NC=0;

	// Lettura del numero di specie
	// l'if necessario per evitare che vengano prese in considerazione delle righe nulle
	// nel calcolo del numero di specie
	ifstream inputFile;
	openInputFileAndControl(inputFile, fileName);

	for(;;)
	{
		inputFile >> name;
		if (!strcmp("END",	name))
			break;
		NC++;
	}

	inputFile.close();


	// Dimensionamento dell'array con i nomi
	names = new string[NC + 1];

	// Lettura dei nomi delle specie
	ifstream inputFileControl;
	openInputFileAndControl(inputFileControl, fileName);

	k=0;
	do
	{
		inputFileControl >> name;
		if(strcmp(name,"\0"))
		//if(!inputFileControl.eof()!=0)
		{
			k++;
			names[k] = name;
		}
	} while (k<NC);

	inputFileControl.close();

	// Controllo che non ci siano specie inserite due volte
	for (k=1;k<=NC;k++)
		for (j=1;j<=NC;j++)
			if( (names[k] == names[j]) && (k!=j))
				ErrorMessage("Una specie e' stata inserita 2 volte! " + names[k]);

	WriteSpeciesNames(binaryFile, asciiFile);
}

void OpenSMOKE_PreProcessorIdealGas::WriteSpeciesNames(BzzSave &binaryFile, BzzSave &asciiFile)
{
	binaryFile << NC;
	asciiFile << NC;
	for(int k=1;k<=NC;k++)
	{
		char name[Constants::NAME_SIZE];
		strcpy(name, names[k].c_str());
		binaryFile.fileSave.write((char*) name, sizeof(name));
		asciiFile << name;
	}
}


void OpenSMOKE_PreProcessorIdealGas::readThermodynamics(const string fileName, BzzSave &binaryFile, BzzSave &asciiFile)
{
	cout << "Reading thermodynamics ... " << endl;

	int k;
	int iFound;
	char name[Constants::NAME_SIZE];

	ifstream inputFile;
	openInputFileAndControl(inputFile, fileName);


	double pm;
	double RGAS_su_PM;
	double a1,a2,a3,a4,a5,a6,a7;
	double b1,b2,b3,b4,b5,b6,b7;
	double t1, t2, t3;

	for(k=1;k<=NC;k++)
	{
		do
		{
			iFound=0;
			inputFile >> name >> pm;
			inputFile >> a1 >> a2 >> a3 >> a4 >> a5 >> a6 >> a7;
			inputFile >> b1 >> b2 >> b3 >> b4 >> b5 >> b6 >> b7;
			inputFile >> t1 >> t2 >> t3;

			if(names[k] == name)
			{
				M[k]=pm;
				RGAS_su_PM=Constants::R_J_mol/(pm*1.e-3);

				CpHT[k][1] = a1*RGAS_su_PM;		// [J/kgK]
				CpHT[k][2] = a2*RGAS_su_PM;
				CpHT[k][3] = a3*RGAS_su_PM;
				CpHT[k][4] = a4*RGAS_su_PM;
				CpHT[k][5] = a5*RGAS_su_PM;

				CpLT[k][1] = b1*RGAS_su_PM;
				CpLT[k][2] = b2*RGAS_su_PM;
				CpLT[k][3] = b3*RGAS_su_PM;
				CpLT[k][4] = b4*RGAS_su_PM;
				CpLT[k][5] = b5*RGAS_su_PM;

				aDH[k][1] = a1;				// L'entalpia viene calcolata in forma
				aDH[k][2] = .5 * a2;		// adimensionale, cio?dividendo per Rgas e
				aDH[k][3] = 1./3. * a3;		// per la temperatura in modo tale da
				aDH[k][4] = .25 * a4;		// poter calcolare pi velocemente le costanti
				aDH[k][5] = .2 * a5;		// di equilibrio
				aDH[k][6] = a6;

				bDH[k][1] = b1;
				bDH[k][2] = .5 * b2;
				bDH[k][3] = 1./3. * b3;
				bDH[k][4] = .25 * b4;
				bDH[k][5] = .2 * b5;
				bDH[k][6] = b6;

				aDS[k][1] = a1;				// L' entropia viene calcolata in questo modo
				aDS[k][2] = a2;				// cio?adimensionale perch?ci?velocizza
				aDS[k][3] = .5 * a3;		// il calcolo delle costanti di equilibrio
				aDS[k][4] = 1./3. * a4;
				aDS[k][5] = .25 * a5;
				aDS[k][6] = a7;

				bDS[k][1] = b1;
				bDS[k][2] = b2;
				bDS[k][3] = .5 * b3;
				bDS[k][4] = 1./3. * b4;
				bDS[k][5] = .25 * b5;
				bDS[k][6] = b7;

				T1[k] = t1;
				T2[k] = t2;
				T3[k] = t3;

				iFound=1;
				break;
			}
		} while (inputFile.eof()==0);

		if (iFound==0)
			ErrorMessage("The species could not be found: " + names[k]);

		inputFile.seekg(0);
	}
	inputFile.close();

	WriteThermodynamicData(binaryFile, asciiFile);
}

void OpenSMOKE_PreProcessorIdealGas::WriteThermodynamicData(BzzSave &binaryFile, BzzSave &asciiFile)
{
	for(int k=1;k<=NC;k++)
	{
		binaryFile << CpHT[k];
		binaryFile << CpLT[k];
		binaryFile << aDH[k];
		binaryFile << bDH[k];
		binaryFile << aDS[k];
		binaryFile << bDS[k];
	}
	binaryFile << M;
	binaryFile << T1;
	binaryFile << T2;
	binaryFile << T3;

	for(int k=1;k<=NC;k++)
	{
		asciiFile << CpHT[k];
		asciiFile << CpLT[k];
		asciiFile << aDH[k];
		asciiFile << bDH[k];
		asciiFile << aDS[k];
		asciiFile << bDS[k];
	}
	asciiFile << M;
	asciiFile << T1;
	asciiFile << T2;
	asciiFile << T3;
}

void OpenSMOKE_PreProcessorIdealGas::readTransportProperties(const string fileName)
{
	cout << "Reading transport properties ... " << endl;

	int k;
	ofstream fLog;
	openOutputFileAndControl(fLog, "Log.log");

	OpenSMOKE_CHEMKINInterpreter_TransportData transport;
	transport.ReadTransportData(fileName, &fLog);
	vector<string> list_of_names;
	list_of_names.push_back("list of names");
	for(k=1;k<=NC;k++)	
		list_of_names.push_back(names[k]);
	
	BzzVectorInt indices = transport.GiveMeSpeciesIndices(list_of_names);
	for(k=1;k<=NC;k++)	
	{
		int j=indices[k];

		shape_factor[k]=transport.species[j].shape_factor;			// indice 0=specie monoatomica 1=c. lineare 2=c.non lineare
		epsylon_over_kb[k]=transport.species[j].epsylon_over_kb;	// [K]
		sigma[k]=transport.species[j].sigma;						// [A]
		mu[k]=transport.species[j].mu* 1.e-6;						// [Debye --> sqrt(A3.erg)]
		alfa[k]=transport.species[j].alfa;							// [A3]
		zRot298[k]=transport.species[j].zRot298;					// [-]
	}

	fLog.close();
}


// ************************************************************************************	//
// ************************************************************************************ //
//								ALLOCAZIONI DI MEMORIA									//
// ************************************************************************************	//
// ************************************************************************************ //

void OpenSMOKE_PreProcessorIdealGas::allocate(void)
{
	// Allocazione di memoria per le variabili
	ChangeDimensions(NC,&shape_factor);
	ChangeDimensions(NC,&epsylon_over_kb);
	ChangeDimensions(NC,&KbSuepsylon);
	ChangeDimensions(NC,&sigma);
	ChangeDimensions(NC,&mu);
	ChangeDimensions(NC,&alfa);
	ChangeDimensions(NC,&zRot298);
	ChangeDimensions(NC,&M);
	ChangeDimensions(NC,&uM);
	ChangeDimensions(NC,&muStar);

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

	ChangeDimensions(NC,&T1);
	ChangeDimensions(NC,&T2);
	ChangeDimensions(NC,&T3);
	ChangeDimensions(NC,&h);
	ChangeDimensions(NC,&h_f);
	ChangeDimensions(NC,&s);
	ChangeDimensions(NC,&Cp);
	ChangeDimensions(NC,&Cv);
}

void OpenSMOKE_PreProcessorIdealGas::allocateDiffusivity(void)
{
	ChangeDimensions(NC,NC,&epsylonjkSuKb);
	ChangeDimensions(NC,NC,&KbSuepsylonjk);
	ChangeDimensions(NC,NC,&TjkStar);
	ChangeDimensions(NC,NC,&deltajkStar);
	ChangeDimensions(NC,NC,&omega11jk);
	ChangeDimensions(NC,NC,&coeff_Djk);
	ChangeDimensions(NC,NC,&Djk);
}

void OpenSMOKE_PreProcessorIdealGas::allocateViscosity(void)
{
	ChangeDimensions(NC,&deltakStar);
	ChangeDimensions(NC,&omega22k);
	ChangeDimensions(NC,&TkStar);
	ChangeDimensions(NC,&eta);
	ChangeDimensions(NC,&coeff_eta);
}

void OpenSMOKE_PreProcessorIdealGas::allocateConductivity(void)
{
	ChangeDimensions(NC,&fVib);
	ChangeDimensions(NC,&fRot);
	ChangeDimensions(NC,&fTrans);
	ChangeDimensions(NC,&cVtrans);
	ChangeDimensions(NC,&cVrot);
	ChangeDimensions(NC,&cVvib);
	ChangeDimensions(NC,&f298);
	ChangeDimensions(NC,&fT);
	ChangeDimensions(NC,&zRot);
	ChangeDimensions(NC,&rho);
	ChangeDimensions(NC,&Dkk);
	ChangeDimensions(NC,&coeff_Dkk);
	ChangeDimensions(NC,&omega11kk);
	ChangeDimensions(NC,&lambda);
}

void OpenSMOKE_PreProcessorIdealGas::allocateThermalDiffusionRatios(void)
{
	for(int j=1;j<=NC;j++)
		if (M[j] <= Constants::MaxMWThermalDiffusionRatios)	iThermalDiffusionRatios.Append(j);
	ChangeDimensions(iThermalDiffusionRatios.Size(), NC, &TetaStar);
}

// ************************************************************************************	//
// ************************************************************************************ //
//										FITTING											//
// ************************************************************************************	//
// ************************************************************************************ //

void OpenSMOKE_PreProcessorIdealGas::Fitting(int nPoints)
{
	viscosityFitting(nPoints, etaMin, etaMax);
	conductivityFitting(nPoints, lambdaMin, lambdaMax);
	diffusivityFitting(nPoints, diffMin, diffMax);
	thermalDiffusionRatiosFitting(nPoints, tetaMin, tetaMax);
}

void OpenSMOKE_PreProcessorIdealGas::WriteFitting(BzzSave &binaryFile, BzzSave &asciiFile)
{
	binaryFile << 1;

	binaryFile << etaMin << etaMax;
	binaryFile << fittingEta;

	binaryFile << lambdaMin << lambdaMax;
	binaryFile << fittingLambda;

	binaryFile << diffMin << diffMax;
	binaryFile << fittingDbinary;

	binaryFile << tetaMin << tetaMax;
	binaryFile << iThermalDiffusionRatios;
	binaryFile << fittingTetaBinary;

	asciiFile << 1;

	asciiFile << etaMin << etaMax;
	asciiFile << fittingEta;

	asciiFile << lambdaMin << lambdaMax;
	asciiFile << fittingLambda;

	asciiFile << diffMin << diffMax;
	asciiFile << fittingDbinary;

	asciiFile << tetaMin << tetaMax;
	asciiFile << iThermalDiffusionRatios;
	asciiFile << fittingTetaBinary;
}

void OpenSMOKE_PreProcessorIdealGas::WriteFittingTransportSensitivity(BzzSave &fOutput)
{
	fOutput << fittingEta;
	fOutput << fittingLambda;
	fOutput << fittingDbinary;
	fOutput << fittingTetaBinary;
}

void OpenSMOKE_PreProcessorIdealGas::viscosityFitting(int nPoints, double &etaMin, double &etaMax)
{
	cout << "       viscosity fitting...  ";

	double start = BzzGetCpuTime();

	int i, j;
	int nCoefficients=4;
	double T, logT;
	double deltaT=(TMAX-TMIN)/(nPoints-1);

	BzzLinearRegression linReg;

	BzzMatrix F(nPoints,nCoefficients);
	BzzVector y(nPoints);
	BzzMatrix matrixEta(NC, nPoints);
	ChangeDimensions(NC,nCoefficients,&fittingEta);


	// FITTING LOGARITMICO PER LA VISCOSITA'
	// - il polinomio e' esattamente quello che viene utilizzato nel CHEMKIN

	// 1. Calcolo delle proprieta'
	T=TMIN;
	for (i=1;i<=nPoints;i++)
	{
		SpeciesViscosity(T);

		logT = log(T);
		F[i][1]=1.;
		for (j=2;j<=nCoefficients;j++)
			F[i][j]=F[i][j-1]*logT;;

		matrixEta.SetColumn(i,eta);

		T+=deltaT;
	}

	// 2. Fitting sulle viscosita'
	for (j=1;j<=NC;j++)
	{
		y = matrixEta.GetRow(j);
		for (i=1;i<=nPoints;i++)
			y[i] = log(y[i]);

		linReg(F,y);
		fittingEta.SetRow(j,linReg.GetParameters());
	}

	// 3. Min and MAx values
	double *s=matrixEta.GetHandle();
	etaMax=0.;
	etaMin=1.e16;
	for (i=1;i<=nPoints;i++)
	for (j=1;j<=NC;j++)
	{
		if(etaMin>(*s) && (*s)>0.) etaMin=(*s);
		if(etaMax<(*s)) etaMax=(*s);
		*s++;
	}
	etaMax=log(etaMax)+0.01;
	etaMin=log(etaMin)-0.01;

	cout << "(" << BzzGetCpuTime() - start << " s)" << endl;
	//cout << etaMax << " " << etaMin << endl;
}

void OpenSMOKE_PreProcessorIdealGas::conductivityFitting(int nPoints, double &lambdaMin, double &lambdaMax)
{
	// FITTING LOGARITMICO PER LA CONDUCIBILITA' TERMICA
	// - il polinomio e' esattamente quello che viene utilizzato nel CHEMKIN

	cout << "       thermal conductivity fitting...  ";

	double start = BzzGetCpuTime();

	int i, j;
	int nCoefficients=4;
	double T, logT;
	double deltaT=(TMAX-TMIN)/(nPoints-1);

	BzzLinearRegression linReg;

	BzzMatrix F(nPoints,nCoefficients);
	BzzVector y(nPoints);

	BzzMatrix matrixLambda(NC, nPoints);
	ChangeDimensions(NC,nCoefficients,&fittingLambda);

	// 1. Calcolo delle propriet?
	T=TMIN;
	for (i=1;i<=nPoints;i++)
	{
		SpeciesCp(T);
		SpeciesCv();
		SpeciesViscosity(T);
		SpeciesConductivity(T);

		logT = log(T);
		F[i][1]=1.;
		for (j=2;j<=nCoefficients;j++)
			F[i][j]=F[i][j-1]*logT;;

		matrixLambda.SetColumn(i,lambda);

		T+=deltaT;
	}

	// 2. Fitting sulle conducibilit?termiche
	for (j=1;j<=NC;j++)
	{
		y=matrixLambda.GetRow(j);
		for (i=1;i<=nPoints;i++)
			y[i] = log(y[i]);

		linReg(F,y);
		fittingLambda.SetRow(j,linReg.GetParameters());
	}

	// 3. Min and Max values
	double *s=matrixLambda.GetHandle();
	lambdaMax=0.;
	lambdaMin=1.e16;
	for (i=1;i<=nPoints;i++)
	for (j=1;j<=NC;j++)
	{
		if(lambdaMin>(*s) && (*s)>0.) lambdaMin=(*s);
		if(lambdaMax<(*s)) lambdaMax=(*s);
		*s++;
	}

	lambdaMax=log(lambdaMax)+0.01;
	lambdaMin=log(lambdaMin)-0.01;

	cout << "(" << BzzGetCpuTime() - start << " s)" << endl;
	//cout << lambdaMax << " " << lambdaMin << endl;
}


void OpenSMOKE_PreProcessorIdealGas::diffusivityFitting(int nPoints, double &diffMin, double &diffMax)
{
	cout << "       diffusivity fitting...  ";

	double start = BzzGetCpuTime();

	int i, j, k;
	int nCoefficients=4;
	double T, logT;
	double deltaT=(TMAX-TMIN)/(nPoints-1);

	BzzLinearRegression linReg;

	BzzMatrix F(nPoints,nCoefficients);
	BzzVector y(nPoints);
	BzzVector b(nCoefficients);
	BzzMatrix matrixDbinary(NC*NC, nPoints);
	ChangeDimensions(NC*NC,nCoefficients,&fittingDbinary);

	// 1. Calcolo delle proprieta'
	T=TMIN;
	for (i=1;i<=nPoints;i++)
	{
		SpeciesDiffusionCoefficient(T,1.);

		logT = log(T);
		F[i][1]=1.;
		for (j=2;j<=nCoefficients;j++)
			F[i][j]=F[i][j-1]*logT;

		for (j=1;j<=NC;j++)
			for (k=1;k<=NC;k++)
				matrixDbinary[k+NC*(j-1)][i]=Djk[j][k];

		T+=deltaT;
	}

	// 2. Fitting sulle diffusivita' binarie
	for (j=1;j<=NC;j++)
		for (k=1;k<=NC;k++)
		{
			if (k!=j)
			{
				y=matrixDbinary.GetRow(k+(j-1)*NC);
				for (i=1;i<=nPoints;i++)
					y[i] = log(y[i]);

				linReg(F,y);
				b=linReg.GetParameters();
			}
			else b=0.;

			fittingDbinary.SetRow(k+(j-1)*NC,b);
		}

	// 3. Min and Max values
	double *s=matrixDbinary.GetHandle();
	diffMax=0.;
	diffMin=1.e16;
	for (i=1;i<=nPoints;i++)
		for (j=1;j<=NC;j++)
			for (k=1;k<=NC;k++)
			{
				if(diffMin>(*s) && (*s)>0.) diffMin=(*s);
				if(diffMax<(*s)) diffMax=(*s);
				*s++;
			}
	diffMax=log(diffMax)+0.01;
	diffMin=log(diffMin)-0.01;

	cout << "(" << BzzGetCpuTime() - start << " s)" << endl;
	//cout << diffMax << " " << diffMin << endl;
}

void OpenSMOKE_PreProcessorIdealGas::thermalDiffusionRatiosFitting(int nPoints, double &tetaMin, double &tetaMax)
{
	cout << "       thermal diffusion coefficients fitting...  ";

	double start = BzzGetCpuTime();

	int i, j, k;
	int nCoefficients=4;
	double T;
	double deltaT=(TMAX-TMIN)/(nPoints-1);

	BzzLinearRegression linReg;

	BzzMatrix F(nPoints,nCoefficients);
	BzzVector y(nPoints);
	BzzVector b(nCoefficients);
	BzzMatrix matrixTetaBinary(NC*iThermalDiffusionRatios.Size(), nPoints);
	ChangeDimensions(NC*iThermalDiffusionRatios.Size(),nCoefficients, &fittingTetaBinary);

	// 1. Calcolo delle proprieta'
	T=TMIN;
	for (i=1;i<=nPoints;i++)
	{
		SpeciesThermalDiffusionRatios(T);

		F[i][1]=1.;
		for (j=2;j<=nCoefficients;j++)
			F[i][j]=F[i][j-1]*T;

		for (k=1;k<=iThermalDiffusionRatios.Size();k++)
			for (j=1;j<=NC;j++)
					matrixTetaBinary[j+NC*(k-1)][i]=TetaStar[k][j];

		T+=deltaT;
	}

	// 2. Fitting sulle diffusivita' binarie
	for (k=1;k<=iThermalDiffusionRatios.Size();k++)
		for (j=1;j<=NC;j++)
		{
			if (iThermalDiffusionRatios[k] != j)
			{
				y=matrixTetaBinary.GetRow(j+NC*(k-1));

				linReg(F,y);
				b=linReg.GetParameters();
			}
			else b=0.;

			fittingTetaBinary.SetRow(j+(k-1)*NC, b);
		}

	// 3. Min and Max values
	double *s=matrixTetaBinary.GetHandle();
	tetaMax=0.;
	tetaMin=1.e16;
	for (i=1;i<=nPoints;i++)
		for (k=1;k<=iThermalDiffusionRatios.Size();k++)
			for (j=1;j<=NC;j++)
			{
				if(tetaMin >(*s)) tetaMin = (*s);
				if(tetaMax <(*s)) tetaMax = (*s);
				*s++;
			}
	tetaMax+=1e-10;
	tetaMin-=1e-10;

	cout << "(" << BzzGetCpuTime() - start << " s)" << endl;
	//cout << tetaMax << " " << tetaMin << endl;	
}

void OpenSMOKE_PreProcessorIdealGas::WriteViscosityFittingCoefficients(ofstream &fOutput)
{
	fOutput << " ----------------------------------------------------------------------------------------------------------------------------------" << endl;
	fOutput << "                       VISCOSITY FITTING COEFFICIENTS    " << endl;
	fOutput << endl;
	fOutput << "             mu = exp( A + B*logT + C*(logT)^2 + D*(logT)^3 )  [kg/m/s]      " << endl;
	fOutput << endl;
	fOutput << "    Species                         A                 B              C                 D              mu(298K)        mu(1000K)    " << endl;
	fOutput << " ----------------------------------------------------------------------------------------------------------------------------------" << endl;
	fOutput << endl;
   
	for(int i=1;i<=NC;i++)
	{
		double logT, mu298, mu1000;

		logT	= log(298.);
		mu298	= exp(fittingEta[i][1]+logT*(fittingEta[i][2]+logT*(fittingEta[i][3]+logT*fittingEta[i][4])));
		
		logT	= log(1000.);
		mu1000	= exp(fittingEta[i][1]+logT*(fittingEta[i][2]+logT*(fittingEta[i][3]+logT*fittingEta[i][4])));

		fOutput << right << setw(5) << i;
		fOutput << ". ";
		fOutput << setw(20) << left  << names[i];

		for(int j=1;j<=4;j++)
			fOutput << setw(17) << scientific << right << setprecision(6) << fittingEta[i][j];
		
		fOutput << setw(17) << scientific << right << setprecision(6) << mu298;
		fOutput << setw(17) << scientific << right << setprecision(6) << mu1000;
		
		fOutput << endl;
	}
	fOutput << endl;
	fOutput << " ----------------------------------------------------------------------------------------------------------------------------------" << endl;
	fOutput << endl << endl;
}

void OpenSMOKE_PreProcessorIdealGas::WriteConductivityFittingCoefficients(ofstream &fOutput)
{
	fOutput << " ----------------------------------------------------------------------------------------------------------------------------------" << endl;
	fOutput << "                     THERMAL CONDUCTIVITY FITTING COEFFICIENTS                                                                     " << endl;
	fOutput << endl;
	fOutput << "             lambda = exp( A + B*logT + C*(logT)^2 + D*(logT)^3 )  [W/m/K]                                                         " << endl;
	fOutput << endl;
	fOutput << "    Species                         A                 B              C                 D           lambda(298K)     lambda(1000K)  " << endl;
	fOutput << " ----------------------------------------------------------------------------------------------------------------------------------" << endl;
	fOutput << endl;
   
	for(int i=1;i<=NC;i++)
	{
		double logT, lambda298, lambda1000;

		logT		= log(298.);
		lambda298	= exp(fittingLambda[i][1]+logT*(fittingLambda[i][2]+logT*(fittingLambda[i][3]+logT*fittingLambda[i][4])));
		
		logT		= log(1000.);
		lambda1000	= exp(fittingLambda[i][1]+logT*(fittingLambda[i][2]+logT*(fittingLambda[i][3]+logT*fittingLambda[i][4])));

		fOutput << right << setw(5) << i;
		fOutput << ". ";
		fOutput << setw(20) << left  << names[i];

		for(int j=1;j<=4;j++)
			fOutput << setw(17) << scientific << right << setprecision(6) << fittingLambda[i][j];
		
		fOutput << setw(17) << scientific << right << setprecision(6) << lambda298;
		fOutput << setw(17) << scientific << right << setprecision(6) << lambda1000;
		
		fOutput << endl;
	}
	fOutput << endl;
	fOutput << " ----------------------------------------------------------------------------------------------------------------------------------" << endl;
	fOutput << endl << endl;
}

void OpenSMOKE_PreProcessorIdealGas::WriteDiffusivityFittingCoefficients(ofstream &fOutput)
{
	fOutput << " -------------------------------------------------------------------------------------------------------------------------------------------------------------" << endl;
	fOutput << "                     BINARY DIFFUSIVITY FITTING COEFFICIENTS                                                                                                  " << endl;
	fOutput << endl;
	fOutput << "             Djk = ( exp( A + B*logT + C*(logT)^2 + D*(logT)^3 ) ) / P    [m2/s]  (P in bar)                                                                  " << endl;
	fOutput << endl;
	fOutput << "             Species             Species                      A                 B              C                 D             Djk(298K)        Djk(1000K)    " << endl;
	fOutput << " -------------------------------------------------------------------------------------------------------------------------------------------------------------" << endl;
	fOutput << endl;
 
	BzzVectorInt iD(NC*(NC-1)/2);
	int m=1;
	int i;
	for (i=1;i<=NC;i++)
		for (int k=i+1;k<=NC;k++)
			iD[m++]=k+(i-1)*NC;

	double P=1.;
	int *s=iD.GetHandle();
	for (i=1;i<=NC;i++)
	{
		for (int k=i+1;k<=NC;k++)
		{
			double logT, Djk298, Djk1000;

			logT	= log(298.);
			Djk298	= exp(fittingDbinary[*s][1]+logT*(fittingDbinary[*s][2]+logT*(fittingDbinary[*s][3]+logT*fittingDbinary[*s][4]))) / P;
	
			logT	= log(1000.);
			Djk1000	= exp((fittingDbinary[*s][1]+logT*(fittingDbinary[*s][2]+logT*(fittingDbinary[*s][3]+logT*fittingDbinary[*s][4])))) / P;

			fOutput << right << setw(5) << i;
			fOutput << right << setw(5) << k;
			fOutput << right << setw(3) << "";
			fOutput << setw(20) << left  << names[i];
			fOutput << setw(20) << left  << names[k];

			for(int j=1;j<=4;j++)
				fOutput << setw(17) << scientific << right << setprecision(6) << fittingDbinary[*s][j];
			
			fOutput << setw(17) << scientific << right << setprecision(6) << Djk298;
			fOutput << setw(17) << scientific << right << setprecision(6) << Djk1000;
			
			fOutput << endl;

			s++;
		}

		fOutput << endl;
	}

	fOutput << " -------------------------------------------------------------------------------------------------------------------------------------------------------------" << endl;
	fOutput << endl << endl;
}

void OpenSMOKE_PreProcessorIdealGas::WriteThermalDiffusionRatiosFittingCoefficients(ofstream &fOutput)
{
	fOutput << " -------------------------------------------------------------------------------------------------------------------------------------------------------------" << endl;
	fOutput << "                   THERMAL DIFFUSION RATIOS FITTING COEFFICIENTS                                                                                                  " << endl;
	fOutput << endl;
	fOutput << "                       Teta_jk = A + B*T + C*T^2 + D*T^3   [-]                                                                  " << endl;
	fOutput << endl;
	fOutput << "             Species             Species                      A                 B              C                 D          Tetajk(298K)      Tetajk(1000K)   " << endl;
	fOutput << " -------------------------------------------------------------------------------------------------------------------------------------------------------------" << endl;
	fOutput << endl;
 
	for (int i=1;i<=iThermalDiffusionRatios.Size();i++)
	{
		for (int k=1;k<=NC;k++)
		{
			double T, Tetajk298, Tetajk1000;
	
			int j = k + (i-1)*NC;

			T = 298.;
			Tetajk298	= fittingTetaBinary[j][1]+T*(fittingTetaBinary[j][2]+T*(fittingTetaBinary[j][3]+T*fittingTetaBinary[j][4]));

			T = 1000.;
			Tetajk1000	= fittingTetaBinary[j][1]+T*(fittingTetaBinary[j][2]+T*(fittingTetaBinary[j][3]+T*fittingTetaBinary[j][4]));

			fOutput << right << setw(5) << i;
			fOutput << right << setw(5) << k;
			fOutput << right << setw(3) << "";
			fOutput << setw(20) << left  << names[iThermalDiffusionRatios[i]];
			fOutput << setw(20) << left  << names[k];

			for(int m=1;m<=4;m++)
				fOutput << setw(17) << scientific << right << setprecision(6) << fittingTetaBinary[j][m];
			
			fOutput << setw(17) << scientific << right << setprecision(6) << Tetajk298;
			fOutput << setw(17) << scientific << right << setprecision(6) << Tetajk1000;
			
			fOutput << endl;
		}

		fOutput << endl;
	}

	fOutput << " -------------------------------------------------------------------------------------------------------------------------------------------------------------" << endl;
	fOutput << endl << endl;


//	for(i=1;i<=20;i++)
//	{
//		double T = 300.+(i-1)*200.;
//		SpeciesThermalDiffusionRatios(T);
//		fOutput << T << "\t" << TetaStar[2][2] << endl;
//	}
}


void OpenSMOKE_PreProcessorIdealGas::WriteCpCoefficients(ofstream &fOutput)
{
	fOutput << " ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------" << endl;
	fOutput << "                     SPECIFIC HEAT COEFFICIENTS                                                                     " << endl;
	fOutput << endl;
	fOutput << "             Cp = A + B*T + C*T^2 + D*T^3 + E*T^4  [J/kg/K]                                                         " << endl;
	fOutput << endl;
	fOutput << "    Species                      Cp(298K)         Cp(1000K)          LT-HT            A(LT)            B(LT)            C(LT)            D(LT)            E(LT)             A(HT)           B(HT)            C(HT)            D(HT)           E(HT)      " << endl;
	fOutput << " ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------" << endl;
	fOutput << endl;
 
	for(int i=1;i<=NC;i++)
	{
		double T, Cp298, Cp1000;

		T = 298.;
		if(T>T2[i])	Cp298 = CpHT[i][1] + T*(CpHT[i][2] + T*(CpHT[i][3] + T*(CpHT[i][4] + T*CpHT[i][5])));
		else		Cp298 = CpLT[i][1] + T*(CpLT[i][2] + T*(CpLT[i][3] + T*(CpLT[i][4] + T*CpLT[i][5]))) ;

		T = 1000.;
		if(T>T2[i])	Cp1000 = CpHT[i][1] + T*(CpHT[i][2] + T*(CpHT[i][3] + T*(CpHT[i][4] + T*CpHT[i][5])));
		else		Cp1000 = CpLT[i][1] + T*(CpLT[i][2] + T*(CpLT[i][3] + T*(CpLT[i][4] + T*CpLT[i][5]))) ;

		fOutput << right << setw(5) << i;
		fOutput << ". ";
		fOutput << setw(20) << left  << names[i];

		fOutput << setw(17) << scientific << right << setprecision(6) << Cp298;
		fOutput << setw(17) << scientific << right << setprecision(6) << Cp1000;
		fOutput << setw(17) << scientific << right << setprecision(6) << T2[i];

		int j;
		for(j=1;j<=5;j++)
			fOutput << setw(17) << scientific << right << setprecision(6) << CpLT[i][j];

		for(j=1;j<=5;j++)
			fOutput << setw(17) << scientific << right << setprecision(6) << CpLT[i][j];
		
		fOutput << endl;
	}
	fOutput << endl;
	fOutput << " ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------" << endl;
	fOutput << endl << endl;
}

void OpenSMOKE_PreProcessorIdealGas::WriteEnthalpyCoefficients(ofstream &fOutput)
{
	fOutput << " -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------" << endl;
	fOutput << "                           ENTHALPY COEFFICIENTS                                                                     " << endl;
	fOutput << endl;
	fOutput << "             H/(RT) = A + B*T + C*T^2 + D*T^3 + E*T^4 + F/T  [-]                                                         " << endl;
	fOutput << endl;
	fOutput << "    Species                   H(298K)[J/kmol]  H(1000K)[J/kmol]      LT-HT            A(LT)            B(LT)            C(LT)            D(LT)            E(LT)           F(LT)             A(HT)            B(HT)            C(HT)           D(HT)             E(HT)            F(HT)      " << endl;
	fOutput << " -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------" << endl;
	fOutput << endl;
 
	for(int i=1;i<=NC;i++)
	{
		double T, H298, H1000;

		T = 298.;
		if (T>T2[i]) H298 = aDH[i][1] + T*(aDH[i][2]+T*(aDH[i][3]+T*(aDH[i][4]+T*aDH[i][5]))) + aDH[i][6]/T;
		else         H298 = bDH[i][1] + T*(bDH[i][2]+T*(bDH[i][3]+T*(bDH[i][4]+T*bDH[i][5]))) + bDH[i][6]/T;
		H298 *= Constants::R_J_kmol*T;	// [J/kmol]

		T = 1000.;
		if (T>T2[i]) H1000 = aDH[i][1] + T*(aDH[i][2]+T*(aDH[i][3]+T*(aDH[i][4]+T*aDH[i][5]))) + aDH[i][6]/T;
		else         H1000 = bDH[i][1] + T*(bDH[i][2]+T*(bDH[i][3]+T*(bDH[i][4]+T*bDH[i][5]))) + bDH[i][6]/T;
		H1000 *= Constants::R_J_kmol*T;	// [J/kmol]

		fOutput << right << setw(5) << i;
		fOutput << ". ";
		fOutput << setw(20) << left  << names[i];

		fOutput << setw(17) << scientific << right << setprecision(6) << H298;
		fOutput << setw(17) << scientific << right << setprecision(6) << H1000;
		fOutput << setw(17) << scientific << right << setprecision(6) << T2[i];

		int j;
		for(j=1;j<=6;j++)
			fOutput << setw(17) << scientific << right << setprecision(6) << bDH[i][j];

		for(j=1;j<=6;j++)
			fOutput << setw(17) << scientific << right << setprecision(6) << aDH[i][j];
		
		fOutput << endl;
	}
	fOutput << endl;
	fOutput << " -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------" << endl;
	fOutput << endl << endl;
}

void OpenSMOKE_PreProcessorIdealGas::WriteEntropyCoefficients(ofstream &fOutput)
{
	fOutput << " -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------" << endl;
	fOutput << "                             ENTROPY COEFFICIENTS                                                                     " << endl;
	fOutput << endl;
	fOutput << "             S/R = A*lnT + B*T + C*T^2 + D*T^3 + E*T^4 + F  [-]                                                         " << endl;
	fOutput << endl;
	fOutput << "    Species                 S(298K)[J/kmol/K] S(1000K)[J/kmol/K]     LT-HT            A(LT)            B(LT)            C(LT)            D(LT)            E(LT)           F(LT)             A(HT)            B(HT)            C(HT)           D(HT)             E(HT)            F(HT)      " << endl;
	fOutput << " -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------" << endl;
	fOutput << endl;
 
	for(int i=1;i<=NC;i++)
	{
		double T, logT, S298, S1000;

		T = 298.;
		logT = log(T);
		if(T > T2[i])	S298 = aDS[i][1]*logT+T*( aDS[i][2] + T*( aDS[i][3] + T*(aDS[i][4] + T*aDS[i][5])))+aDS[i][6];
		else			S298 = bDS[i][1]*logT+T*( bDS[i][2] + T*( bDS[i][3] + T*(bDS[i][4] + T*bDS[i][5])))+bDS[i][6];
		S298 *= Constants::R_J_kmol;	// [J/kmol/K]

		T = 1000.;
		logT = log(T);
		if(T > T2[i])	S1000 = aDS[i][1]*logT+T*( aDS[i][2] + T*( aDS[i][3] + T*(aDS[i][4] + T*aDS[i][5])))+aDS[i][6];
		else			S1000 = bDS[i][1]*logT+T*( bDS[i][2] + T*( bDS[i][3] + T*(bDS[i][4] + T*bDS[i][5])))+bDS[i][6];
		S1000 *= Constants::R_J_kmol;	// [J/kmol/K]

		fOutput << right << setw(5) << i;
		fOutput << ". ";
		fOutput << setw(20) << left  << names[i];

		fOutput << setw(17) << scientific << right << setprecision(6) << S298;
		fOutput << setw(17) << scientific << right << setprecision(6) << S1000;
		fOutput << setw(17) << scientific << right << setprecision(6) << T2[i];

		int j;
		for(j=1;j<=6;j++)
			fOutput << setw(17) << scientific << right << setprecision(6) << bDS[i][j];

		for(j=1;j<=6;j++)
			fOutput << setw(17) << scientific << right << setprecision(6) << aDS[i][j];
		
		fOutput << endl;
	}
	fOutput << endl;
	fOutput << " -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------" << endl;
	fOutput << endl << endl;
}



int OpenSMOKE_PreProcessorIdealGas::recognize_species(char* name)
{
	for(int i=1;i<=NC;i++)
		if (names[i] == name)
			return i;
	
	string dummy = name;
	ErrorMessage("This species is not included in the kinetic scheme: " + dummy);
	return -1;
}

int OpenSMOKE_PreProcessorIdealGas::recognize_species_without_exit(char* name)
{
	for(int i=1;i<=NC;i++)
		if (names[i] == name)
			return i;

	return 0;
}

int OpenSMOKE_PreProcessorIdealGas::recognize_species(string name)
{
	for(int i=1;i<=NC;i++)
		if (!name.compare(names[i]))
			return i;
	
	ErrorMessage("This species is not included in the kinetic scheme: " + name);
	return -1;
}


int OpenSMOKE_PreProcessorIdealGas::NumberOfSpecies()
{
	return NC;
}

int	OpenSMOKE_PreProcessorIdealGas::NumberOfElements()
{
	return elements_matrix.Rows();
}

void OpenSMOKE_PreProcessorIdealGas::ReadElements(const string fileElements, BzzSave &binaryFile, BzzSave &asciiFile)
{
	char comment[Constants::COMMENT_SIZE];
	double dummy_double;
	string dummy;
	string name;
	vector<string> elements_list;
	vector<double> m_elements_list;
	ifstream fInput;
	
	// Open input file: list of species and elements
	openInputFileAndControl(fInput, fileElements);

	// Reading the element list
	cout << "Open Elements..." << endl;
	fInput >> dummy;
	if (dummy != "ELEMENTS")	ErrorMessage("The Elements.bzz file must begin with ELEMENTS");

	for(;;)
	{
		fInput >> dummy;
		if (dummy == "END")	break;
		else 
		{
			fInput >> dummy_double;
			elements_list.push_back(dummy);
			m_elements_list.push_back(dummy_double);
		}
	}
	fInput.seekg(ios::beg);

	ChangeDimensions(elements_list.size(), NC, &elements_matrix);

	// Reading the species
	for(int i=1;i<=NC;i++)
	{
		bool iFound = false;

		fInput >> dummy;
		for(unsigned int k=1; k<=elements_list.size(); k++)
			fInput >> dummy >> dummy_double;
		fInput >> dummy;

		while(iFound == false)
		{
			fInput >> name;
			
			if (name == "END")	
				ErrorMessage("This species was not included in the Elements.bzz file: " + names[i]);

			if (name == names[i])
			{
				iFound = true;

				bool iEnd = false;
				while(iEnd == false)
				{
					string name_element;
					double number_element;

					fInput >> name_element;
					if( name_element != "/")
					{
						fInput >> number_element;
						int index_element = seek_index(name_element, elements_list);
						elements_matrix[index_element + 1][i] = number_element;
					}
					else
					{
						iEnd = true;
						fInput.getline(comment, Constants::COMMENT_SIZE);
					}
				}

				// Checking the molecular weigth
				double MW =0;
				for (unsigned int j=1;j<=elements_list.size();j++)
					MW += elements_matrix[j][i] * m_elements_list[j-1];
				if (MW > M[i]*(1.+1.e-10) || MW < M[i]*(1.-1.e-10))
					ErrorMessage("The elemental composition is incorrect: " + names[i]);
			}
			else fInput.getline(comment, Constants::COMMENT_SIZE);
		}

		fInput.seekg(ios::beg);
	}
	fInput.close();
	
	// Cutting the element matrix
	BzzVectorInt delete_list;
	for(unsigned int k=1; k<=elements_list.size(); k++)
	{
		BzzVector row = elements_matrix.GetRow(k);
		if (row.GetSumElements() == 0)
			delete_list.Append(k);
		else
		{
			list_of_elements.push_back(elements_list[k-1]);
			elements_mw.Append(m_elements_list[k-1]);
		}
	}	
	if (delete_list.Size() != 0)
		elements_matrix.DeleteRows(delete_list);

	WriteElementsData(binaryFile, asciiFile);
}

void OpenSMOKE_PreProcessorIdealGas::WriteElementsData(BzzSave &binaryFile, BzzSave &asciiFile)
{
	// Write on file: element matrix
	binaryFile << elements_matrix;
	binaryFile << elements_mw;
	asciiFile << elements_matrix;
	asciiFile << elements_mw;
	// Write on file: element list
	for(int k=1;k<=elements_matrix.Rows();k++)
	{
		char name[Constants::NAME_SIZE];
		strcpy(name, list_of_elements[k-1].c_str());
		binaryFile.fileSave.write((char*) name, sizeof(name));
		asciiFile << name;
	}
}

void OpenSMOKE_PreProcessorIdealGas::SpeciesCp(const double T)
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

void OpenSMOKE_PreProcessorIdealGas::SpeciesCv(void)
{
	for(int k = 1;k <= NC;k++)
		Cv[k] = Cp[k] - Constants::R_J_kmol*uM[k];
}


void OpenSMOKE_PreProcessorIdealGas::SetupFromCHEMKINInterpreter(OpenSMOKE_CHEMKINInterpreter *chemkin, const string file_name)
{
	int k;

	NC			  = chemkin->kinetics.species_list.size()-1;
	TMIN          = chemkin->TMIN;
	TMAX		  = chemkin->TMAX;
	fittingPoints = chemkin->FITTINGPOINTS;
	name_object   = chemkin->name_object;
	name_author   = chemkin->name_author;
	building_date = chemkin->building_date;

	// Species names
	names = new string[NC+1];
	for(k=1;k<=NC;k++)
		names[k] = chemkin->kinetics.species_list[k];

	// Allocazione della memoria per le variabili
	allocate();

	// Lettura delle informazioni termodinamiche
	for(k=1;k<=NC;k++)
	{
		aDH[k] = chemkin->thermo.aDH[k];
		bDH[k] = chemkin->thermo.bDH[k];
		aDS[k] = chemkin->thermo.aDS[k];
		bDS[k] = chemkin->thermo.bDS[k];

		CpHT[k] = chemkin->thermo.CpHT[k];
		CpLT[k] = chemkin->thermo.CpLT[k];
	}
	T1 = chemkin->thermo.T1;
	T2 = chemkin->thermo.T2;
	T3 = chemkin->thermo.T3;
	M = chemkin->thermo.M;

	// Lettura delle proprieta di trasporto
	if (chemkin->transport.IsActivated() == true)
	{
		shape_factor	=	chemkin->transport.shape_factor;
		epsylon_over_kb	=	chemkin->transport.epsylon_over_kb;
		sigma			=	chemkin->transport.sigma;			
		mu				=	chemkin->transport.mu;		
		alfa			=	chemkin->transport.alfa;				
		zRot298			=	chemkin->transport.zRot298;			
	}

	// Lettura degli elementi
	elements_matrix	 = chemkin->thermo.elements_matrix;
	elements_mw		 = chemkin->thermo.elements_mw;
	list_of_elements = chemkin->thermo.list_of_elements;

	// Calcolo delle grandezze indipendenti dalla temperatura, pressione e composizione
	Initialize();

	// Fitting
	if (chemkin->transport.IsActivated() == true)
	{
		InitializeTransportProperties();
		Fitting(fittingPoints);
	}

	// Write On Binary File
	BzzSave binaryFile;
	BzzSave asciiFile;
	binaryFile('*', file_name);
	string file_name_ascii = file_name + ".ascii";
	asciiFile(file_name_ascii);
	WriteHeaderFile(binaryFile, asciiFile);
	WriteSpeciesNames(binaryFile, asciiFile);
	WriteThermodynamicData(binaryFile, asciiFile);
	WriteElementsData(binaryFile, asciiFile);
	if (chemkin->transport.IsActivated() == true)	WriteFitting(binaryFile, asciiFile);
	else											{ binaryFile << 0; asciiFile << 0; }
	binaryFile.End();
	asciiFile.End();

	// Provisional ()
/*	ofstream fOpenSMOKE_R("Provisional.out", ios::out);
	fOpenSMOKE_R << NC << endl;
	for(int j=1;j<=NC;j++)
	{
		fOpenSMOKE_R << names[j] << endl;
		fOpenSMOKE_R << "[notassigned]" << endl;
		fOpenSMOKE_R << 'G' << endl;
		fOpenSMOKE_R << M[j] << endl;
		fOpenSMOKE_R << T2[j] << endl;
		if(names[j]=="OH")
			cout << fittingLambda[j][1] << " " << fittingLambda[j][2] << " " << fittingLambda[j][3] << " " << fittingLambda[j][4] << endl;
		for(int k=1;k<=4;k++)
			fOpenSMOKE_R << fittingEta[j][k] << endl;
		for(int k=1;k<=4;k++)
			fOpenSMOKE_R << fittingLambda[j][k] << endl;
		for(int kk=1;kk<=NC;kk++)
			for(int k=1;k<=4;k++)
				fOpenSMOKE_R << fittingDbinary[kk+(j-1)*NC][k] << endl;
		for(int k=1;k<=5;k++)
			fOpenSMOKE_R << CpHT[j][k] << endl;
		for(int k=1;k<=5;k++)
			fOpenSMOKE_R << CpLT[j][k] << endl;
		for(int k=1;k<=6;k++)
			fOpenSMOKE_R << aDH[j][k] << endl;
		for(int k=1;k<=6;k++)
			fOpenSMOKE_R << bDH[j][k] << endl;
		for(int k=1;k<=6;k++)
			fOpenSMOKE_R << aDS[j][k] << endl;
		for(int k=1;k<=6;k++)
			fOpenSMOKE_R << bDS[j][k] << endl;
	}
	fOpenSMOKE_R.close();*/
}


void OpenSMOKE_PreProcessorIdealGas::SpeciesThermalDiffusionRatios(const double T)
{
	int i,j;

	// DETERMINE A*, B*, AND C* FOR EACH SPECIES
	for(i=1;i<=iThermalDiffusionRatios.Size();i++)
		for(j=1;j<=NC;j++)
        {
			int k = iThermalDiffusionRatios[i];

			double TSLOG = log(T/epsylonjkSuKb[k][j]);
			double T1 = TSLOG;
			double T2 = TSLOG*T1;
			double T3 = TSLOG*T2;
			double T4 = TSLOG*T3;
			double T5 = TSLOG*T4;
			double T6 = TSLOG*T5;

			double ASTAR = FITASTAR[1]+FITASTAR[2]*T1+FITASTAR[3]*T2+FITASTAR[4]*T3+FITASTAR[5]*T4+FITASTAR[6]*T5+FITASTAR[7]*T6;
			double BSTAR = FITBSTAR[1]+FITBSTAR[2]*T1+FITBSTAR[3]*T2+FITBSTAR[4]*T3+FITBSTAR[5]*T4+FITBSTAR[6]*T5+FITBSTAR[7]*T6;
			double CSTAR = FITCSTAR[1]+FITCSTAR[2]*T1+FITCSTAR[3]*T2+FITCSTAR[4]*T3+FITCSTAR[5]*T4+FITCSTAR[6]*T5+FITCSTAR[7]*T6;

			TetaStar[i][j] = 7.5*(2.*ASTAR+5.)*(6.*CSTAR-5.)/
				             (ASTAR*(16.*ASTAR-12.*BSTAR+55.)) *
                             (M[k]-M[j])/(M[k]+M[j]);
		}
}

void OpenSMOKE_PreProcessorIdealGas::WriteNewOutputFile(ofstream &fOutput)
{
	fOutput.setf(ios::scientific);
	fOutput.precision(14);

	fOutput << NC << endl;

	for(int i=1;i<=NC;i++)
	{
		fOutput << names[i] << endl;
		fOutput << names[i] << endl;

		int count = 0;
		for(int k=1;k<=NumberOfElements();k++)
			if (elements_matrix[k][i] > 0.) count++;

		fOutput << count << endl;
		for(int k=1;k<=NumberOfElements();k++)
			if (elements_matrix[k][i] > 0.) fOutput << list_of_elements[k-1] << endl;
		
		fOutput << count << endl;
		for(int k=1;k<=NumberOfElements();k++)
			if (elements_matrix[k][i] > 0.) fOutput << elements_matrix[k][i] << endl;
		
		// **************************************************************************** //
		//							Thermodynamics      								//
		// **************************************************************************** //

		fOutput << CpHT[i][1]*M[i]/Constants::R_J_kmol << endl;
		fOutput << CpHT[i][2]*M[i]/Constants::R_J_kmol << endl;
		fOutput << CpHT[i][3]*M[i]/Constants::R_J_kmol << endl;
		fOutput << CpHT[i][4]*M[i]/Constants::R_J_kmol << endl;
		fOutput << CpHT[i][5]*M[i]/Constants::R_J_kmol << endl;
		fOutput <<  aDH[i][6] << endl;
		fOutput <<  aDS[i][6] << endl;

		fOutput << CpLT[i][1]*M[i]/Constants::R_J_kmol << endl;
		fOutput << CpLT[i][2]*M[i]/Constants::R_J_kmol << endl;
		fOutput << CpLT[i][3]*M[i]/Constants::R_J_kmol << endl;
		fOutput << CpLT[i][4]*M[i]/Constants::R_J_kmol << endl;
		fOutput << CpLT[i][5]*M[i]/Constants::R_J_kmol << endl;
		fOutput <<  bDH[i][6] << endl;
		fOutput <<  bDS[i][6] << endl;

		fOutput <<  T1[i] << endl;
		fOutput <<  T3[i] << endl;
		fOutput <<  T2[i] << endl;

		fOutput << M[i] << endl;

		// **************************************************************************** //
		//							Transport Properties								//
		// **************************************************************************** //
		
		// Molecular weight
		fOutput << M[i] << endl;

		// Conductivity
		for(int j=1;j<=4;j++)
			fOutput << fittingLambda[i][j] << endl;

		// Viscosity
		for(int j=1;j<=4;j++)
			fOutput << fittingEta[i][j] << endl;

		// Diffusivity
		fOutput << NC << endl;
		int iRow = (i-1)*NC;
		for(int k=1;k<=NC;k++)
			for(int j=1;j<=4;j++)
				fOutput << fittingDbinary[iRow+k][j] << endl;

		// Thermal diffusion ratios
		if (M[i] <= Constants::MaxMWThermalDiffusionRatios)
		{
			int index = 0;
			for (int k=1;k<=iThermalDiffusionRatios.Size();k++)
				if (iThermalDiffusionRatios[k] == i)
				{	
					index = k;
					break;
				}
				if (index == 0)
					ErrorMessage("Wrong Thermal Diffusion Ratios exporting");

			int iRow = (index-1)*NC;
			for(int k=1;k<=NC;k++)
				for(int j=1;j<=4;j++)
					fOutput << fittingTetaBinary[iRow+k][j] << endl;
		}
		else
		{
			for(int k=1;k<=4*NC;k++)
				fOutput << 0. << endl;
		}
	}
}