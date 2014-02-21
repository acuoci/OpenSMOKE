/***************************************************************************
 *   Copyright (C) 2006-2009 by Alberto Cuoci						       *
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

#include "Utilities.h"
#include "flamelet_group.h"
#include "engine/OpenSMOKE_ReactingGas.h"
#include "idealreactors/OpenSMOKE_GasStream.h"


void checkingWords(const string foundWord, const string expectedWord, const string fileName)
{
	if (foundWord != expectedWord)
	{
		cout << "ERROR in reading file: expected: " << expectedWord << "  found: " << foundWord << endl;
		cout << "Please check your Flamelet File: " << fileName << " !!" << endl;
		cout << "Press enter to continue..." << endl;
		getchar();
		exit(-1);
	}
}

void extract_name(char *word, int &compositionLabel, std::string &name)
{
	if (!strcmp(word, "massfraction-"), 13)
		compositionLabel = 'W';
	else if (!strcmp(word, "molefraction-"), 13)
		compositionLabel = 'X';
	else
	{
		cout << "ERROR in reading file: expected: massfraction- or molefraction- " << endl;
		cout << "       found: " << word << endl;
		cout << "Please check your Flamelet File!!" << endl;
		cout << "Press enter to continue..." << endl;
		getchar();
		exit(-1);
	}

	std::string stringa = word;
	name = stringa.substr(13, stringa.length()-13);
	name = uppercase(name);
}
	
void import_flamelet_file_from_FLUENT(OpenSMOKE_ReactingGas &mix, string fileName, flamelet_group &flameletGroup)
{
	int i,j;
	string word;
	char word_name_species[30];
	ifstream iFile;
	openInputFileAndControl(iFile, fileName);
	flameletGroup.initialize();
	
	int k=0;
	while (!iFile.eof())
	{
		k++;

		int numberOfSpecies;
		int gridPoints;
		double chi;
		double pressure_Pa;

		BzzVectorInt	compositionLabel;
		BzzVector csi;
		BzzVector temperature;
		BzzMatrix composition;
		vector<string> names;


		iFile >> word ;
		if (!strcmp(word.c_str(), "TEMPERATURE"))
			break;
		checkingWords(word, "HEADER", fileName);

		iFile >> word ;	
		checkingWords(word.c_str(), "CHI", fileName);
		iFile >> chi;

		iFile >> word ;
		checkingWords(word.c_str(), "NUMOFSPECIES", fileName);
		iFile >> numberOfSpecies;

		iFile >> word ;
		checkingWords(word.c_str(), "GRIDPOINTS", fileName);
		iFile >> gridPoints;

		iFile >> word ;
		checkingWords(word.c_str(), "PRESSURE", fileName);
		iFile >> pressure_Pa;

		iFile >> word ;
		checkingWords(word.c_str(), "BODY", fileName);

		cout << "Flamelet n. " << k << " - Scalar dissipation rate (st.): " << chi << " Hz" << endl;
		cout << "Number Of Species: " << numberOfSpecies << endl;
		cout << "Number Of Points:  " << gridPoints << endl;

		names.resize(numberOfSpecies+1);
		ChangeDimensions(numberOfSpecies, &compositionLabel);

		ChangeDimensions(gridPoints, &csi);
		ChangeDimensions(gridPoints, &temperature);
		ChangeDimensions(gridPoints, numberOfSpecies, &composition);

		iFile >> word ;
		checkingWords(word, "Z", fileName);
		for(i=1;i<=gridPoints;i++)
			iFile >> csi[i];

		iFile >> word ;
		checkingWords(word, "TEMPERATURE", fileName);
		for(i=1;i<=gridPoints;i++)
			iFile >> temperature[i];

		for(j=1;j<=numberOfSpecies;j++)
		{
			iFile >> word_name_species;
			extract_name(word_name_species, compositionLabel[j], names[j]);
			for(i=1;i<=gridPoints;i++)
				iFile >> composition[i][j];
		}

		// Checking on composition
		for(j=1;j<=numberOfSpecies;j++)
			if (compositionLabel[j] != 'W')
			{
				cout << "The composition must be in mass fraction!!!" << endl;
				cout << "Please correct your Flamelet.fla Input File!" << endl;
				cout << "Found: " << compositionLabel[j] << endl;
				cout << "Press enter to continue..." << endl;
				getchar();
				exit(-1);
			}

		cout << "Maximum temperature: " << temperature.Max() << " K" << endl;
		cout << endl;

		flamelet flame;
		flame.assign_flamelet(mix, names, pressure_Pa, chi, csi, temperature, composition);
		flameletGroup.append(flame);
		flameletGroup.reorder();
		flameletGroup.lock();
	}
}


void conversion_from_FlameCRECK_to_FLUENT(OpenSMOKE_ReactingGas &mix, string fileNameSource, ofstream &fOutput, bool verboseFormationRates)
{
	int i, j;
	const int SIZE = 30000;
	char commento[SIZE];

	char name_lower_case[50];
	vector<string> nameReducedSpecies;
	
	int nReducedSpecies;
	int nPoints;
	int iCount;

	double pressure;
	double chiSt;
	
	BzzVector coordinate;
	BzzVector T;
	BzzMatrix X;
	BzzMatrix Y;
	double dummy;


	// -------------------------------------------------------------------	
	// Number of Species and number of points
	// -------------------------------------------------------------------
	ifstream fInput;
	openInputFileAndControl(fInput, fileNameSource);
	fInput.getline(commento, SIZE);

		
	// Lettura del numero di punti e delle coordinate
	// ---------------------------------------------------------------------------
	do 
	{
		fInput >> dummy;
		fInput >> dummy; coordinate.Append(dummy);
		fInput.getline(commento, SIZE);
	} while (coordinate[coordinate.Size()]<1.00);
	fInput.close();
	nPoints = coordinate.Size();

	// Dimensionamento dei vettori
	// ---------------------------------------------------------------------------
	ChangeDimensions(nPoints, &T);
	ChangeDimensions(nPoints, mix.NumberOfSpecies(), &X);
	ChangeDimensions(nPoints, mix.NumberOfSpecies(), &Y);

	// -------------------------------------------------------------------	
	// Number of Species and number of points
	// -------------------------------------------------------------------
	nReducedSpecies = mix.NumberOfSpecies();
	nameReducedSpecies.push_back("species");
	for(j=1;j<=nReducedSpecies;j++)
		nameReducedSpecies.push_back(mix.names[j]);

	// Reading all the information
	// ---------------------------------------------------------------------------
	openInputFileAndControl(fInput, fileNameSource);
	fInput.setf(ios::scientific);

	// Lettura delle informazioni
	fInput.getline(commento, SIZE);
	for(i=1;i<=nPoints;i++)
	{
		fInput >> dummy;							// point index
		fInput >> dummy; coordinate[i] = dummy;		// coordinate
		fInput >> dummy; T[i] = dummy;				// temperature
		fInput >> dummy;							// enthalpy mole	
		fInput >> dummy;							// enthalpy mass
		fInput >> pressure;							// pressure
		fInput >> chiSt;							// chiSt [Hz]
		
		for(j=1;j<=mix.NumberOfSpecies();j++)	fInput >> X[i][j];
		for(j=1;j<=mix.NumberOfSpecies();j++)	fInput >> Y[i][j];

		fInput.getline(commento, SIZE);
	}
	fInput.close();

	
	// Write on file
	fOutput << "HEADER" << endl;
	fOutput << "CHI " << chiSt << endl;
	fOutput << "NUMOFSPECIES " << mix.NumberOfSpecies() << endl;
	fOutput << "GRIDPOINTS " << nPoints << endl;
	fOutput << "PRESSURE " << pressure << endl;
	fOutput << "BODY" << endl;
	
	fOutput << "Z" << endl;
	iCount = 1;
	for (i=nPoints;i>=1;i--)
	{
		fOutput << 1.-coordinate[i] << "   ";
		if (iCount==5 && i!=1)
		{	iCount = 0;	fOutput << endl;	}
		iCount++;
	}
	fOutput << endl;

	fOutput << "TEMPERATURE" << endl;
	iCount = 1;
	for (i=nPoints;i>=1;i--)
	{
		fOutput << T[i] << "   ";
		if (iCount==5 && i!=1)
		{	iCount = 0;	fOutput << endl;	}
		iCount++;
	}
	fOutput << endl;

	for(j=1;j<=mix.NumberOfSpecies();j++)
		for(int jR=1;jR<=nReducedSpecies;jR++)
			if ( mix.names[j] == nameReducedSpecies[jR])
			{
				strcpy(name_lower_case, "massfraction-");
				strcat(name_lower_case, mix.names[j].c_str());
				for(int k=0; k<50; k++)
					name_lower_case[k] = tolower(name_lower_case[k]);
				fOutput << name_lower_case << endl;
				iCount = 1;
				for (i=nPoints;i>=1;i--)
				{
					fOutput << std::max(Y[i][j], 1.e-48);
					if (j!=mix.NumberOfSpecies()) 
						fOutput << "   ";
					if (j==mix.NumberOfSpecies() && i!=1) 
						fOutput << "   ";
					if (iCount==5 && i!=1)
					{	iCount = 0;	fOutput << endl;	}
					iCount++;
				}
				if (j!=mix.NumberOfSpecies()) fOutput << endl;
			}


	if (verboseFormationRates == true)
	{
		BzzMatrix FR(nPoints, mix.NumberOfSpecies());
		{
			BzzVector mass_fractions(mix.NumberOfSpecies());
			BzzVector mole_fractions(mix.NumberOfSpecies());
			BzzVector concentrations(mix.NumberOfSpecies());
			BzzVector R(mix.NumberOfSpecies());
			double MWtot;
			for (i=nPoints;i>=1;i--)
			{
				
				mix.ComputeKineticParameters( T[i], log(T[i]), 1./T[i], pressure);

				Y.GetRow(i, &mass_fractions);
				mix.GetMWAndMoleFractionsFromMassFractions(MWtot, mole_fractions, mass_fractions);

				double cTot    = pressure  / (Constants::R_J_kmol*T[i]);
				concentrations = cTot*mole_fractions;

				mix.ComputeFromConcentrations( T[i], concentrations, cTot, &R);		// [kmol/m3/s]
				ElementByElementProduct(R, mix.M, &R);								// [kg/m3/s]

				FR.SetRow(i, R);

				cout << "Formation rates: " << i << " " << T[i] << " " << pressure << " " << R.Max() << " " << R.Min() << endl;
			}

		}

		fOutput << endl;
		for(j=1;j<=mix.NumberOfSpecies();j++)
			for(int jR=1;jR<=nReducedSpecies;jR++)
				if ( mix.names[j] == nameReducedSpecies[jR])
				{
					strcpy(name_lower_case, "form-rate-");
					strcat(name_lower_case, mix.names[j].c_str());
					for(int k=0; name_lower_case[k]; k++)
						name_lower_case[k] = tolower(name_lower_case[k]);
					fOutput << name_lower_case << endl;
					iCount = 1;
					for (i=nPoints;i>=1;i--)
					{
						fOutput << FR[i][j];
						if (j!=mix.NumberOfSpecies()) 
							fOutput << "   ";
						if (j==mix.NumberOfSpecies() && i!=1) 
							fOutput << "   ";
						if (iCount==5 && i!=1)
						{	iCount = 0;	fOutput << endl;	}
						iCount++;
					}
					if (j!=mix.NumberOfSpecies()) fOutput << endl;
				}
	}
	fOutput << endl;
}


void search_index_for_linear_interpolation(BzzVector &x, double mean, int &iMean, double &interpolationFactor)
{
	int N = x.Size();
	for(int i=2;i<=N;i++)
		if (x[i]>=mean)
		{
			iMean =i;
			interpolationFactor = (mean-x[iMean-1]) / (x[iMean]-x[iMean-1]);
			break;
		}
}

double interpolate_function(BzzVector &v, int iMean, double interpolationFactor)
{
	return v[iMean-1] +  (v[iMean]-v[iMean-1]) * interpolationFactor;
}

void double_the_vector(BzzVector &original, BzzVector &doublevector)
{
	int i;
	ChangeDimensions(original.Size()*2-1, &doublevector);

	for(i=1;i<=original.Size();i++)
		doublevector[2*i-1] = original[i];
	for(i=1;i<=original.Size()-1;i++)
		doublevector[2*i] = (original[i] + original[i+1])*0.50;
}

void double_the_vector(int nRefinements, BzzVector &original, BzzVector &doublevector)
{	
	if (nRefinements == 0)
		doublevector = original;
	else
	{
		BzzVector auxiliary;
		double_the_vector(original, doublevector);
		for(int i=2;i<=nRefinements;i++)
		{
			auxiliary = doublevector;
			double_the_vector(auxiliary, doublevector);
		}
	}
}


std::string lowercase (std::string command) // converts all names input to lower case
{
	int length;
	length = command.length();				// calculates command length
	for (int i=0; i<length; i++)			// converts letters in sequence to lower case stopping when it gets to the end
		command[i] = tolower(command[i]);
	return command;
} 

std::string uppercase (std::string command) // converts all names input to upper case
{
	int length;
	length = command.length();				// calculates command length
	for (int i=0; i<length; i++)			// converts letters in sequence to upper case stopping when it gets to the end
		command[i] = toupper(command[i]);
	return command;
} 

void give_K_Planck_SANDIA(double T, double &kH2O, double &kCO2, double &kCO, double &kCH4)
{
	double uT = 1000./T;

	// Carbon dioxide [1/(m.atm)]	
	kH2O = -0.23093 +uT*(-1.1239+uT*(9.4153 +uT*(-2.9988 +uT*( 0.51382 + uT*(-1.8684e-5)))));
	
	// Water [1/(m.atm)]
	kCO2 =  18.741  +uT*(-121.31+uT*(273.5  +uT*(-194.05 +uT*( 56.31   + uT*(-5.8169   )))));

	// Carbon moxide [1/(m.atm)]
	kCO=0.;
	kCH4=0.;
/*	if (T <= 750.)
		kCO = 4.7869 + T*(-0.06953 + T*(2.95775e-4 + T*(-4.25732e-7 + T*2.02894e-10)));
	else
		kCO = 10.09 + T*(-0.01183 + T*(4.7753e-6 + T*(-5.87209e-10 + T*-2.5334e-14)));

	// Methane [1/(m.atm)]
	kCH4 = 6.6334 + T*( -0.0035686 + T*(1.6682e-08 + T*(2.5611e-10 -2.6558e-14*T)));*/
}



void log_normal_distribution::construct_sequence(double _chiMean, BzzVector &discreteChi)
{
	int j;

	double sigma = 2.00;
	double sqrt_2 = sqrt(2.);

	double coefficientA = 1./(sqrt_2*sigma);
	double coefficientB = sigma/(2.*sqrt_2);

	chiMean = _chiMean+1.e-32;
	N = discreteChi.Size() + 1;
	ChangeDimensions(N, &chi);
	ChangeDimensions(N, &teta);
	ChangeDimensions(N, &error_function);
	ChangeDimensions(N-1, &integral);
	
	// Scalar Dissipation Rate Sequence
	chi[1] = 0.;
	for(j=2; j<=N-1; j++)
		chi[j] = 0.50*(discreteChi[j-1]+discreteChi[j]);
	chi[N] = discreteChi[N-1] + (discreteChi[N-1]-chi[N-1]);

	// Auxiliary Variable Sequence
	teta[1] = 0.;
	for(j=2; j<=N; j++)
		teta[j] = coefficientA*log(chi[j]/chiMean) - coefficientB;

	// Error Function Sequence
	error_function[1] = -1.;
	for(j=2; j<=N; j++)
		error_function[j] = BzzErf(teta[j]);

	// Integral Sequence
	for(j=1; j<=N-1; j++)
		integral[j]=0.50*(-error_function[j]+error_function[j+1]);

	double sumchecking = integral.GetSumElements();
	integral /= sumchecking;	// normalization (very important)
}

double log_normal_distribution::apply(BzzVector &Beta)
{
	double sum=0.;
	for(int j=1; j<=N-1; j++)
		sum+=Beta[j]*integral[j];
	return sum;
}


void read_SootProperties_from_FlameCRECK(const string fileNameSource, BzzVector &m0, BzzVector &fv)
{
	int i;
	const int SIZE = 30000;
	char commento[SIZE];
	int nPoints;
	BzzVector coordinate;
	double dummy;

	// -------------------------------------------------------------------	
	// Number of Species and number of points
	// -------------------------------------------------------------------
	ifstream fInput;
	openInputFileAndControl(fInput, fileNameSource);

	// Lettura riga intestazione
	fInput.getline(commento, SIZE);
		
	// Lettura del numero di punti e delle coordinate
	// ---------------------------------------------------------------------------
	i=0;
	do 
	{
		i++;
		fInput >> dummy; coordinate.Append(dummy);
		fInput.getline(commento, SIZE);
	} 
	while (coordinate[i]<1.00);
	fInput.close();

	nPoints = coordinate.Size();

	// Dimensionamento dei vettori
	// ---------------------------------------------------------------------------
	ChangeDimensions(nPoints, &m0);
	ChangeDimensions(nPoints, &fv);

	// Reading all the information
	// ---------------------------------------------------------------------------
	openInputFileAndControl(fInput, fileNameSource);
	fInput.setf(ios::scientific);

	// Lettura riga intestazione
	fInput.getline(commento, SIZE);

	// Lettura delle informazioni
	for(i=1;i<=nPoints;i++)
	{
		fInput >> coordinate[i];				// coordinate
		fInput >> dummy;						// phiN
		fInput >> dummy;						// phiM
		
		fInput >> m0[i];						// soot particle density [1./m3]	
		fInput >> fv[i];						// soot volume fraction [-]
		if (fv[i] < 1.e-14) fv[i] = 1.e-14;		// Correction

		fInput.getline(commento, SIZE);
	}
	fInput.close();

	Reverse(&m0);
	Reverse(&fv);
}

void extract_data_2d(char* source_file_name, char* dest_file_name, int xcolumn, int ycolumn)
{
	int nColumns;
	int nRows;
	double mixture_fraction;
	double variance;
	double dummy;
	double request = 0.;

	ifstream fInput(source_file_name, ios::in);
	ofstream fOutput(dest_file_name, ios::out);
	fOutput.setf(ios::scientific);

	string line;
	getline(fInput, line);	
	
	// TODO: recognize the number of columns in the data
	nColumns = 59;

	// Recognize the number of rows
	nRows = 0;
	for(;;)
	{	
		fInput >> dummy;
		getline(fInput, line);	
		if (dummy==0.) 
			nRows++;
		else
			break;
	}

	cout << nRows << endl;
	fInput.seekg(0);
	getline(fInput, line);	

	// Mixture Fraction or Variance
	if (xcolumn==1)
	{
		while(!fInput.eof())
		{
			fInput >> mixture_fraction;
			fInput >> dummy;
			fInput >> variance;
			if (variance == request)
			{
				fOutput << mixture_fraction << "\t";
				for(int j=4;j<=ycolumn-1;j++)
					fInput >> dummy;
				fInput >> dummy;
				fOutput << dummy << endl;
			}
			getline(fInput, line);	
		}

		fInput.close();
		fOutput.close();
	}

	else if (xcolumn==3)
	{
		while(!fInput.eof())
		{
			fInput >> mixture_fraction;
			if (mixture_fraction==request)
			{
				fInput >> dummy;
				fInput >> variance;
				fOutput << variance << "\t";
				for(int j=4;j<=ycolumn-1;j++)
					fInput >> dummy;
				fInput >> dummy;
				fOutput << dummy << endl;
			}
			getline(fInput, line);	
		}

		fInput.close();
		fOutput.close();
	}
	else
	{
		cout << "Error!";
		cout << "Press enter to continue..." << endl;
		getchar();
		exit(-1);
	}
	
}

void extract_data_2d(char* source_file_name, char* dest_file_name, int xcolumn)
{
	int nColumns;
	int nRows;
	double mixture_fraction;
	double variance;
	double dummy;
	double request = 0.;

	ifstream fInput(source_file_name, ios::in);
	ofstream fOutput(dest_file_name, ios::out);
	fOutput.setf(ios::scientific);

	string line;
	
	nColumns = 4;

	// Recognize the number of rows
	nRows = 0;
	for(;;)
	{	
		fInput >> dummy;
		getline(fInput, line);	
		if (dummy==0.) 
			nRows++;
		else
			break;
	}

	cout << nRows << endl;
	fInput.seekg(0);

	// Mixture Fraction or Variance
	if (xcolumn==1)
	{
		while(!fInput.eof())
		{
			fInput >> mixture_fraction;
			fInput >> variance;
			if (variance == request)
			{
				fOutput << mixture_fraction << "\t";
				fInput >> dummy;
				fOutput << dummy << "\t";
				fInput >> dummy;
				fOutput << dummy << endl;
			}
			getline(fInput, line);	
		}

		fInput.close();
		fOutput.close();
	}

	else if (xcolumn==2)
	{
		while(!fInput.eof())
		{
			fInput >> mixture_fraction;
			if (mixture_fraction==request)
			{
				fInput >> variance;
				fOutput << variance << "\t";
				fInput >> dummy;
				fOutput << dummy << "\t";
				fInput >> dummy;
				fOutput << dummy << endl;
			}
			getline(fInput, line);	
		}

		fInput.close();
		fOutput.close();
	}
	else
	{
		cout << "Error!";
		cout << "Press enter to continue..." << endl;
		getchar();
		exit(-1);
	}
}
