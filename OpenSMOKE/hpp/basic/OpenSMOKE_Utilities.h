/***************************************************************************
 *   Copyright (C) 2006-2008 by Alberto Cuoci   	                       *
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

#ifndef OPENSMOKE_UTILITIES
#define OPENSMOKE_UTILITIES

#include <vector>
#include "BzzMath.hpp"
using namespace std;

class OpenSMOKE_IdealGas;
class OpenSMOKE_ReactingGas;
class OpenSMOKE_RateOfProductionAnalysis;
class OpenSMOKE_SensitivityAnalysis;
class OpenSMOKE_Nu;

class OpenSMOKE_NuManager
{
public:

	OpenSMOKE_NuManager();
	void Clean();
	void Set(const int index, const std::string name);
	void Set(const double coefficient, const int reaction);
	void PrintOnFile(ofstream &fOutput);

	void SaveToBinaryFile(BzzSave &fSave);
	void LoadFromBinaryFile(BzzLoad &fLoad);

	int				nReactions;
	BzzVectorInt	iReactions;
	BzzVector		nuReactions;

private:

	int		index_species;
	std::string	name_species;
};

class OpenSMOKE_Nu
{
public:

	OpenSMOKE_Nu(OpenSMOKE_ReactingGas *mix);		// Default constructor
	OpenSMOKE_Nu(BzzLoad &fLoad);					// Default constructor
	OpenSMOKE_NuManager	*nu;						// Stoichiometric coefficients

	BzzVectorInt number_of_reactions;	// Number of reactions for each species 
	BzzVectorInt sparse_rows;			// Sparsity Pattern
	BzzVectorInt sparse_columns;		// Sparsity Pattern
	BzzVectorInt sparse_index;			// Sparsity Pattern

	void PrintOnFile(const std::string fileName);		// 
	void SaveToBinaryFile(BzzSave &fSave);		// 
	void SparsityPattern(BzzMatrixSparse &M);

private:

	int NC;
	int NR;

	void BuildNuMatrix();
	void SparsityPattern();
	void Initialize();
	OpenSMOKE_ReactingGas *mix;
};

class OpenSMOKE_RateOfProductionCoefficient
{
public:

	OpenSMOKE_RateOfProductionCoefficient() {};
	
	BzzVector *p;
	BzzVector *d;
	BzzVectorInt n;

	void Initialize(BzzVectorInt &n_);
	void Clean();
	void LoadFromBinaryFile(BzzLoad &fLoad);
	void SaveToBinaryFile(BzzSave &fSave);
	void LoadFromBinaryFile(BzzLoad &fLoad, BzzVectorInt &indices);
	void SaveToBinaryFile(BzzSave &fSave, BzzVectorInt &indices);
};

class OpenSMOKE_ChebishevPolynomialsReaction
{
public:
	OpenSMOKE_ChebishevPolynomialsReaction(void);
	void ReadFromFile(BzzLoad &fInput);
	double GiveMeKappa(const double T, const double P);
	void WriteCoefficientMatrix(ofstream &fOutput);

private:

	double GiveMePhi(const int n, const double x);
	void ErrorMessage(const std::string message);
	void WarningMessage(const std::string message);

private:

	int N;
	int M;

	BzzMatrix a;

	double Tmin;
	double Tmax;

	double Pmin;
	double Pmax;

	double log10_Pmin;
	double log10_Pmax;

	double conversion;

	BzzVector phi_n;
	BzzVector phi_m;
};

class OpenSMOKE_LogarithmicPressureReaction
{
public:
	OpenSMOKE_LogarithmicPressureReaction(void);
	void ReadFromFile(BzzLoad &fInput);
	double GiveMeKappa(const double T, const double Pressure);
	void WriteCoefficientMatrix(ofstream &fOutput);

private:

	void ErrorMessage(const std::string message);
	void WarningMessage(const std::string message);

private:

	int N;

	BzzVector P;
	BzzVector lnP;
	BzzVector lnA;
	BzzVector Beta;
	BzzVector EoverR;
};



// ************************************************************* //
//					Linear Interpolation Classes                 //
// ************************************************************* //

class linearInterpolationFixedStep
{
private:
	BzzVector x,y,dyStar;
	double xMin, udx;

public:

	double operator ()(double t);
	void operator ()(BzzVector &x, BzzVector &y);
};

class linearInterpolationFixedStep_2D
{
private:
	int nc;
	BzzVector x;
	BzzMatrix y,dyStar;
	double xMin, udx;

public:

	void operator ()(double t, BzzVector &f);
	void operator ()(BzzVector &x, BzzMatrix &y);
};

class LinearInterpolation
{
public:
    
    LinearInterpolation();
		void 		SetName(std::string name);
		double 	operator ()(double t);
		void 		operator ()(BzzVector &x, BzzVector &y);

private:
	
		BzzVector x,y;
		double xMin, xMax;
		int n;
	
		std::string name_object;
        
  	void ErrorMessage(const std::string message);
};

std::string GiveMeTimeAndDate();
int seek_index(std::string name, vector<string> list_names);

void GiveMeIndicesAndNames(OpenSMOKE_IdealGas &mix, const vector<string> _names, BzzVectorInt &index, vector<string> &names);


// **************************************************************************************************** //
//					Molecular Weight Classes                                                            //
// **************************************************************************************************** //

void massFractionsAndPMtot(BzzVector &molar, BzzVector &mass, double MW, BzzVector &PM);
void massFractionsAndPM(BzzVector &molar, BzzVector &mass, double MW, BzzVector &PM);
void molarFractionsAndPM(BzzVector &molar, BzzVector &mass, double MW, BzzVector &PM);

// **************************************************************************************************** //
//					Utilities                                                                           //
// **************************************************************************************************** //
void my_itoa(int value, char* str, int base);
int myMaxInt(int a, int b);

class BzzInverseErrorFunction
{
private:
	BzzVector c;
	BzzVector C;

public:
	BzzInverseErrorFunction();
	double at(const double x);
};

// **************************************************************************************************** //
//					Gaussian Function                                                                   //
// **************************************************************************************************** //
class gaussianFunction
{
private:
	double mu, sigma;
	double xCenter;

	double w, reductionCoefficient, peak;

public:
	double giveValue(double x);
	void setup(double xCenter, double w, double reductionCoefficient, double peak);
};



// **************************************************************************************************** //
//					File Manager                                                                        //
// **************************************************************************************************** //

void openInputFileAndControl(ifstream &inputFile, const char *fileName);
void openOutputFileAndControl(ofstream &inputFile, const char *fileName);
void openInputFileAndControl_BinaryMode(ifstream &outputFile, const char *fileName);
void openOutputFileAndControl_BinaryMode(ofstream &outputFile, const char *fileName);

void openInputFileAndControl(ifstream &inputFile, std::string fileName);
void openOutputFileAndControl(ofstream &inputFile, std::string fileName);
void openInputFileAndControl_BinaryMode(ifstream &outputFile, std::string fileName);
void openOutputFileAndControl_BinaryMode(ofstream &outputFile, std::string fileName);

// **************************************************************************************************** //
//					Grid Zone Class                                                                     //
// **************************************************************************************************** //
class gridZones
{
public:
	int nPoints;
	BzzVector x;

	double L;
	double xCenter;
	double wMix;
	double xLeft, xRight;

	void setup(int _nPoints, double  _L, double _xCenter, double _wMix);
	void setup(int _nPoints, double  _L, double _wMix);
};

// **************************************************************************************************** //
//					Useful Functions                                                                    //
// **************************************************************************************************** //
void functionRampa(LinearInterpolation &interpolatedProfile, gridZones &grid, double Left, double Right);
void functionGaussian(LinearInterpolation &interpolatedProfile, gridZones &grid, double peak, double reductionFactor);
void functionStep(LinearInterpolation &interpolatedProfile, gridZones &grid, double Left, double peak, double Right);
void functionStepTwin(LinearInterpolation &interpolatedProfile, gridZones &grid, double Left, double peak);


// **************************************************************************************************** //
//					Utilities BzzLibraries --> Lapack                                                   //
// **************************************************************************************************** //
void FromCToBzzVector(int N, BzzVector &B, double *L);
void FromBzzVectorToC(int N, BzzVector &B, double *L);
void FromBzzMatrixToC(int NR, int NC, BzzMatrix &B, double *L);
void FromCToBzzMatrix(int NR, int NC, BzzMatrix &B, double *L);



// **************************************************************************************************** //
//					Name List Class                                                                     //
// **************************************************************************************************** //
class nameList
{
public:
	int N;
	char **name;

	void setup(int _N);
};


// **************************************************************************************************** //
//					Spark Class                                                                         //
// **************************************************************************************************** //
class spark
{
private:
	double timeSpark;
	double density;
	double xCenter, xRight, xLeft;
	double width;
	int shape;
	gaussianFunction gaussianProfile;

public:
	double giveEnergy(double t, double x);
	void setup(char *shape, double time, double density, double xCenter, double width);
};

// **************************************************************************************************** //
//					OpenSMOKE logo                                                  //
// **************************************************************************************************** //

void OpenSMOKE_logo(std::string application_name, std::string version, std::string date);
void OpenSMOKE_logo(ofstream &fOutput, std::string application_name);


// **************************************************************************************************** //
//					Parser Class                                                    //
// **************************************************************************************************** //

#include "interfaces/SimpleOpt.h"

class ParserClass
{
    public:

    void setup(int _argc, char* _argv[], CSimpleOpt::SOption *_g_rgOptions);
    int parse(std::string label, std::string &argument);
    int parse(std::string label);

    void virtual ShowUsage();

private:

    int     argc;
    char*   argv[100];
    CSimpleOpt::SOption *g_rgOptions;
    void first_screening(CSimpleOpt args);
};

// *********************************************************************** //
//							Matrix Class                                   //
// *********************************************************************** //

template <class T>
class OpenSMOKE_Matrix
{
private :

    int rows, cols;
    T * element;
    void ErrorMessage(std::string message);

public :

    OpenSMOKE_Matrix (int r=0, int c=0);
    OpenSMOKE_Matrix (int r, int c, T value);
    OpenSMOKE_Matrix (const OpenSMOKE_Matrix<T> & m);


    ~OpenSMOKE_Matrix () { delete [] element; }

    int Rows () const       { return rows; }
    int Columns () const    { return cols; }

    T & operator ( ) (int i, int j) const;
    OpenSMOKE_Matrix<T> & operator = (const OpenSMOKE_Matrix<T> & m);

    void Show () const;
    void Resize(int r, int c);

};

template <class T>
OpenSMOKE_Matrix<T>::OpenSMOKE_Matrix (int r, int c)
{
    if (r<0 || c<0)
        ErrorMessage("Wrong indices");

    rows = r;
    cols = c;
    element = new T [r*c];
}

template <class T>
OpenSMOKE_Matrix<T>::OpenSMOKE_Matrix(int r, int c, T value)
{
    if (r<0 || c<0)
        ErrorMessage("Wrong indices");

    rows = r;
    cols = c;
    element = new T [r*c];

    for (int i=0; i<rows*cols; i++)
        element[i] = value;
}

template <class T>
void OpenSMOKE_Matrix<T>::Resize(int r, int c)
{
    if (r<0 || c<0)
        ErrorMessage("Wrong indices");

    rows = r;
    cols = c;
    delete[] element;

    element = new T [r*c];

    for (int i=0; i<rows*cols; i++)
        element[i] = 0.;
}

template <class T>
OpenSMOKE_Matrix<T>::OpenSMOKE_Matrix(const OpenSMOKE_Matrix<T> & m)
{
    rows = m.rows;
    cols = m.cols;
    element = new T [rows*cols];

    for (int i=0; i<rows*cols; i++)
        element[i] = m.element[i];
}

template <class T>
T & OpenSMOKE_Matrix<T>::operator () (int i, int j) const
{
    if (i<1 || i>rows || j<1 || j>cols)
    {
        cout << "Wrong indices";
        exit (1);
    }

    return element [ (i-1)*cols + j - 1 ];
}


template <class T>
void OpenSMOKE_Matrix<T>::Show () const
{
    for (int i=0; i<rows; i++)
    {
        for (int j=0; j<cols; j++)
            cout << element[(i*cols)+j] << '\t';
        cout << '\n';
    }
}

template <class T>
void OpenSMOKE_Matrix<T>::ErrorMessage(std::string message)
{
    cout << "FATAL ERROR: OpenSMOKE Matrix Class" << endl;
    cout << message << endl;
    cout << "Press a key to continue... " << endl;
    getchar();
    exit(-1);
}

std::string GiveMeFileNameFromFullPath(const std::string pathname);
std::string GiveMeFolderPathFromFullPath(const std::string pathname);

bool caseInsCharCompareN(char a, char b); 
bool caseInsCompare(const std::string& s1, const std::string& s2); 
bool CheckForCommentLine(std::string &line);
bool CheckForCommentLineFromStart(std::string &line);
bool CheckForBlankLine(const std::string line);
bool CheckForEndLine(std::string &line);
void StringSubstitution(std::string &original, const std::string old_string, const std::string new_string);

bool StringSubstitutionAll(std::string &s, const std::string p, const std::string new_string);
bool StringFindSubString(std::string &s, const std::string p);
void SeparateNumberFromString(const std::string expression, std::string &name, double &number);
void SeparateInstructions(const std::string line, vector<string> &instructions);
int SeparateNumberFromStringForElements(const std::string expression, std::string &name, double &number);

int StringFind(const std::string s, const std::string subs);
void CleanFromBlanks(std::string &name);
void CleanFromBlanksForThermo(std::string &name);
void TokenizeString(const std::string line, vector<string> &list);
void TokenizeString(char *line, vector<string> &list);

void WriteUnits_cm_mol_s(ofstream &fOutput, const double lambda);
void WriteUnits_m_kmol_s(ofstream &fOutput, const double lambda);

double KineticTemperature(const double Cc, const double Tmean, const double Tatt, const double n);

std::string GetLabelIndex(const int count);
std::string GetNumber(const int count);

void PrintTagOnGnuplotLabel(const int length, ofstream &fOut, const std::string tag, int &counter);
void PrintHeaderOnBinaryFile(BzzSave &fSave);
void PrintEndOnBinaryFile(BzzSave &fSave);
void CheckInBinaryFile(BzzLoad &fLoad, const std::string flag);
std::string NextInBinaryFile(BzzLoad &fLoad);
void PrepareLiquidPropertiesDictionary(const std::string tag, const std::string file_name, vector<string> &lines, BzzVectorInt &indexLines);
std::string StringToUpper(std::string strToConvert);
std::string StringToLower(std::string strToConvert);

#endif	// OPENSMOKE_UTILITIES

