/***************************************************************************
 *   Copyright (C) 2006 by Alberto Cuoci   	   *
 *   alberto.cuoci@polimi.it   						   *
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

#ifndef OPENSMOKE_PREPROCESSORIDEALGAS_H
#define OPENSMOKE_PREPROCESSORIDEALGAS_H

#include "Linux_Windows.h"
#include "basic/OpenSMOKE_Utilities.h"

class OpenSMOKE_CHEMKINInterpreter;

class OpenSMOKE_PreProcessorIdealGas 
{
protected:
	int NC;

private:

	// Constants
	static const double BOLTZMANN;
	static const double BOLTZMANN3;
	static const double MU_MIN;
	static const double ZROTA;
	static const double ZROTB;
	static const double ZROTC;
	static const double CONST_2SUPI;
	static const double CONST_5SU3R;

	//-----------------------------------------------------------------------------------//
	//										VARIABILI									 //
	//-----------------------------------------------------------------------------------//

	BzzVector *CpHT;	// calore specifico a pressione costante (High Temperature)
	BzzVector *CpLT;	// calore specifico a pressione costante (Low Temperature)

	// Variabili di ingresso per la lettura dei dati relativi alle proprieta di trasporto
	BzzVector			epsylon_over_kb;	// Potenziale di Lennard-Jones (diviso per la costante di Boltzmann)
	BzzVector			sigma;			// Diametro di Lennard-Jones
	BzzVector			zRot298;		// Numero collisionale di rilassamento rotazionale a 298K
	BzzVector			mu;				// Momento dipolare
	BzzVector			alfa;			// Polarizzabilita
	BzzVectorInt		shape_factor;			// 0 = specie monoatomica 1 = catena lineare 2 = catena non lineare

	// Variabili utilizzzate solo per il calcolo delle conducibilita delle singole specie
	BzzVector			f298;			// Funzione di Parker-Brau-Jonkman a 298K
	BzzVector			fT;				// Funzione di Parker-Brau-Jonkman
	BzzVector			zRot;			// Numero collisionale di rilassamento rotazionale
	BzzVector			rho;			// Densit?del singolo componente
	BzzVector			Dkk;			// Coefficiente di autodiffusione
	BzzVector			fVib;
	BzzVector			fTrans;
	BzzVector			fRot;
	BzzVector			cVtrans;		// Contributo traslazionale al CV
	BzzVector			cVrot;			// Contributo rotazionale al CV
	BzzVector			cVvib;			// Contributo vibrazionale al CV


	// Variabili utilizzate solo per il calcolo delle diffusivita
	BzzMatrix			epsylonjkSuKb;	// Potenziale effettivo di Lennard-Jones (diviso per la costante di Boltzmann)
	BzzMatrix			KbSuepsylonjk;  // 1. / epsylonjkSuKb
	BzzMatrix			TjkStar;		// Temperatura effettiva ridotta
	BzzMatrix			deltajkStar;	// Momento dipolare ridotto effettivo (II)
	BzzMatrix			omega11jk;		// Integrale collisionale effettivo
	BzzMatrix			coeff_Djk;		// Coefficiente costante per il calcolo di Djk

	BzzVector			KbSuepsylon;	// Potenziale di Lennard-Jones (diviso per la costante di Boltzmann)
	BzzVector			TkStar;			// Temperatura ridotta
	BzzVector			deltakStar;		// Momento dipolare ridotto (II)
	BzzVector			omega22k;		// Integrale collisionale per la viscosita
	BzzVector			coeff_eta;

	// Serve sia per la diffusivita che per la viscosita
	BzzVector			epsylon;	// Potenziale di Lennard-Jones
	BzzVector			muStar;		// Momento dipolare ridotto
	BzzVector			coeff_Dkk;
	BzzVector			omega11kk;

	// Thermal diffusion ratios
	BzzVectorInt iThermalDiffusionRatios;
	BzzMatrix TetaStar;

	// initialize constructors
	void Initialize(void);
	void InitializeTransportProperties(void);
	double collisionIntegral_11(double tjk, double djk);


	// Funzioni per il calcolo delle diffusivita
	void TemperaturaRidottaComponentePuro(double T);
	void Omega11jk(void);
	void CoefficienteDiDiffusioneBinario(double T, double P);

	// Funzioni per il calcolo della viscosita
	void TemperaturaRidotta(double T);
	void Omega22k(void);
	void ViscositaDelComponente(double T);

	// Funzioni per il calcolo della conducibilita
	void Omega11kk(void);
	void CoefficienteDiAutoDiffusione(double T);
	void computeFT(double T);
	void computeZrot(void);
	void computeCvVib(void);
	void ConducibilitaComponente(void);

	// Matrici delle costanti di fitting
	BzzMatrix fittingEta;
	BzzMatrix fittingLambda;
	BzzMatrix fittingDbinary;
	BzzMatrix fittingTetaBinary;

	int fittingPoints;
	double etaMin, etaMax;
	double lambdaMin, lambdaMax;
	double diffMin, diffMax;
	double tetaMin, tetaMax;

	void viscosityFitting(int nPoints, double &etaMin, double &etaMax);
	void conductivityFitting(int nPoints, double &lambdaMin, double &lambdaMax);
	void diffusivityFitting(int nPoints, double &diffMin, double &diffMax);
	void thermalDiffusionRatiosFitting(int nPoints, double &tetaMin, double &tetaMax);

	// Specific heat functions
	void SpeciesCp(double T);
	void SpeciesCv(void);

	// Proprieta singole specie (trasporto)
	void SpeciesDiffusionCoefficient(double T, double P);
	void SpeciesViscosity(double T);
	void SpeciesConductivity(double T);
	void SpeciesDensity(double T);
	void SpeciesThermalDiffusionRatios(const double T);

	void readNames(const string fileName, BzzSave &binaryFile, BzzSave &asciiFile);
	void readThermodynamics(const string fileName, BzzSave &binaryFile, BzzSave &asciiFile);
	void readTransportProperties(const string fileName);

	void allocate(void);
	void allocateDiffusivity(void);
	void allocateViscosity(void);
	void allocateConductivity(void);
	void allocateThermalDiffusionRatios(void);

	void initializeDiffusivity(void);
	void initializeViscosity(void);
	void initializeConductivity(void);

public:

	// Variabili di ingresso per la lettura dei dati termodinamici
	BzzVector *aDH; // alta T
	BzzVector *bDH; // bassa T
	BzzVector *aDS; // alta T
	BzzVector *bDS; // bassa T
	BzzVector T1;
	BzzVector T2;
	BzzVector T3;

	// Proprieta delle singole specie, accessibili dall'utente
	BzzVector			Cp;				// Calore specifico molare a pressione costante
	BzzVector			Cv;				// Calore specifico molare a volume costante
	BzzVector			h;
	BzzVector			h_f;
	BzzVector			s;
	BzzVector			lambda;
	BzzVector			eta;
	BzzMatrix			Djk;

	// Inizializzazione della libreria
	void Setup(const string fileNames, const string fileTransport, const string fileThermo, const string fileElements, BzzSave &binaryFile, BzzSave &asciiFile);
	void Fitting(int nPoints);
	
	void SetupFromCHEMKINInterpreter(OpenSMOKE_CHEMKINInterpreter *chemkin, const string file_name);
	void SetupTransportSensitivity(const string pathName, const string fileNames, const string fileTransport, const string fileThermo, const string fileElements, BzzSave &binaryFile, BzzSave &asciiFile,const double eps);
	
	void WriteHeaderFile(BzzSave &binaryFile, BzzSave &asciiFile);
	void WriteSpeciesNames(BzzSave &binaryFile, BzzSave &asciiFile);
	void WriteThermodynamicData(BzzSave &binaryFile, BzzSave &asciiFile);
	void WriteElementsData(BzzSave &binaryFile, BzzSave &asciiFile);
	void WriteFitting(BzzSave &binaryFile, BzzSave &asciiFile);
	void WriteFakeFitting(BzzSave &binaryFile, BzzSave &asciiFile);
	void WriteFittingTransportSensitivity(BzzSave &fOutput);

	void WriteViscosityFittingCoefficients(ofstream &fOutput);
	void WriteConductivityFittingCoefficients(ofstream &fOutput);
	void WriteDiffusivityFittingCoefficients(ofstream &fOutput);
	void WriteThermalDiffusionRatiosFittingCoefficients(ofstream &fOutput);


	void WriteCpCoefficients(ofstream &fOutput);
	void WriteEnthalpyCoefficients(ofstream &fOutput);
	void WriteEntropyCoefficients(ofstream &fOutput);

	void WriteNewOutputFile(ofstream &fOutput);


	BzzVector	M;
	BzzVector	uM;				// Peso molecolare [g/mol]
	string *names;

	int recognize_species(char* name);
	int recognize_species_without_exit(char* name);
    int recognize_species(string name);


	int NumberOfSpecies();
	int	NumberOfElements();

	OpenSMOKE_PreProcessorIdealGas();

	void SetMinimumTemperature(const double tmin);
	void SetMaximumTemperature(const double tmax);
	void SetName(const string name);
	void SetAuthorName(const string name);

	string name_object;
	string name_author;
	string building_date;

	double TMIN;
	double TMAX;

	void ErrorMessage(const string message);
	void WarningMessage(const string message);

	void ReadElements(const string fileElements, BzzSave &binaryFile, BzzSave &asciiFile);

	// Element matrix
	vector<string>	list_of_elements;
	BzzVector elements_mw;
	BzzMatrix elements_matrix;

	// Elements
	BzzVector m_elements;
	int  recognize_element(const string name);
};

#endif // OPENSMOKE_PREPROCESSORIDEALGAS_H
