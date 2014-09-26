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
#include <sstream>
#include <iomanip>
#include "basic/OpenSMOKE_Utilities.h"
#include "engine/OpenSMOKE_IdealGas.h"
#include "engine/OpenSMOKE_ReactingGas.h"
#include "basic/OpenSMOKE_Conversions.h"
#include "basic/OpenSMOKE_Constants.h"

int seek_index(std::string name, vector<string> list_names)
{
	for(int i=0; i<int(list_names.size()); i++)
		if (name == list_names[i])
			return i;

	cout << "Element was not found: " << name << endl;
	cout << "Press enter to continue... " << endl;
	getchar();
	exit(-1);
}

std::string GiveMeTimeAndDate()
{
	time_t current_time;
	tm *local_time;
	
	time(&current_time);
	local_time = localtime(&current_time);
	std::string time_string = asctime(local_time);

	return time_string;
}

void GiveMeIndicesAndNames(OpenSMOKE_IdealGas &mix, const vector<string> _names, BzzVectorInt &index, vector<string> &names)
{
	if (_names[0] == "ALL")
	{
		ChangeDimensions(mix.NumberOfSpecies(), &index);
		for(int k=1;k<=mix.NumberOfSpecies();k++)
		{
			index[k] = k;
			names.push_back(mix.names[k]);
		}
	}
	else
	{
		ChangeDimensions(_names.size(), &index);
		for(int k=0;k<int(_names.size());k++)
		{
			index[k+1] = mix.recognize_species(_names[k]);
			names.push_back(_names[k]);
		}
	}
}


void linearInterpolationFixedStep::operator () (BzzVector &X, BzzVector &Y)
{
	int i;
	int n;
	double dx;

	n=X.Size();
	if(n!=Y.Size()) BzzError("Errore in LinearInterpolationFixedStep\n Dim(x) != Dim(y)");
	
	dx = X[2]-X[1];
	udx=1./dx;
	xMin=X[1];

	for(i=2;i<=n;i++)
		if ( fabs( ((X[i]-X[i-1])-dx)/X[i])  > 1.e-9 ) BzzError("deltaX non costante\n");
	
	x = X;
	y = Y;
	ChangeDimensions(n-1,&dyStar);
	for(i=1;i<=n-1;i++)
		dyStar[i] = (y[i+1]-y[i])*udx;
}

double linearInterpolationFixedStep::operator() (double t)
{
	int i;
	i = int((t-xMin)*udx)+1;
	return y[i] + (t-x[i])* dyStar[i];
}


void linearInterpolationFixedStep_2D::operator () (BzzVector &X, BzzMatrix &Y)
{
	int i,j;
	int n;
	double dx;

	nc = Y.Columns();

	n=X.Size();
	if(n!=Y.Rows()) BzzError("Errore in LinearInterpolationFixedStep\n Dim(x) != Dim(y)");
	
	dx = X[2]-X[1];
	udx=1./dx;
	xMin=X[1];

	for(i=2;i<=n;i++)
		if ( fabs( ((X[i]-X[i-1])-dx)/X[i])  > 1.e-9 ) BzzError("deltaX non costante\n");
	
	x = X;
	y = Y;
	ChangeDimensions(n-1,nc,&dyStar);
	for(i=1;i<=nc;i++)
		for(j=1;j<=n-1;j++)
			dyStar[j][i] = (y[j+1][i]-y[j][i])*udx;
}

void linearInterpolationFixedStep_2D::operator() (double t, BzzVector &f)
{
	int i;
	double dt;
	i = int((t-xMin)*udx)+1;
	dt = (t-x[i]);

	for(int k=1;k<=nc;k++)
		f[k] =	y[i][k] + dt * dyStar[i][k];
}
/*
double BzzErfInv(const double x)
{
	double pi=acos(-1.);
	return sqrt(pi)* ( 0.50*x + 
		               1./24.*pi*BzzPow3(x) + 
					   7./960.*pi*pi*BzzPow5(x) + 
		               127./80640.*BzzPow3(pi)*BzzPow7(x) + 
		               4369./11612160.*BzzPow4(pi)*BzzPow9(x) + 
		               34807./364953600.*BzzPow5(pi)*BzzPow10(x)*x
					 );
}*/

BzzInverseErrorFunction::BzzInverseErrorFunction()
{
	ChangeDimensions(30, &c);
	ChangeDimensions(30, &C);

	c[1]  = 1./1.;
	c[2]  = 1./1.;
	c[3]  = 7./6.;
	c[4]  = 127./90.;
	c[5]  = 4369./2520.;
	c[6]  = 34807./16200.;
	c[7]  = 20036983./7484400.;
	c[8]  = 2280356863./681080400.;
	c[9]  = 49020204823./11675664000.;
	c[10] = 65967241200001./12504636144000.;
	c[11] = 15773461423793767./2375880867360000.;
	c[12] = 655889589032992201./78404068622880000.;
	c[13] = 94020690191035873697./8910391798788480000.;
	c[14] = 655782249799531714375489./49229914688306352000000.;
	c[15] = 44737200694996264619809969./2658415393168543008000000.;

	for(int k=15;k<=29;k++)
		for(int m=0;m<=k-1;m++)
			c[k+1] += c[m+1]*c[k-1-m+1]/(m+1.)/(2.*m+1.);

	double PI = acos(-1.);	 
	for(int j=1;j<=30;j++)
	{
		int k = j-1;
		C[j] = c[j]/(2.*k+1.)*pow(sqrt(PI)/2.,2.*k+1.);
	}
}

double BzzInverseErrorFunction::at(const double x)
{
	double x2 = x*x;
	return x*(C[1]+x2*(C[2]+x2*(C[3]+x2*(C[4]+x2*(C[5]+x2*(C[6]+x2*(C[7]+x2*(C[8]+x2*(C[9]+x2*(C[10]+x2*(C[11]+x2*(C[12]+x2*(C[13]+x2*(C[14]+x2*(C[15]+x2*(C[16]+x2*(C[17]+x2*(C[18]+x2*(C[19]+x2*(C[20]+x2*(C[21]+x2*(C[22]+x2*(C[23]+x2*(C[24]+x2*(C[25]+x2*(C[26]+x2*(C[27]+x2*(C[29]+x2*(C[29]+x2*(C[30]))))))))))))))))))))))))))))));
}

void massFractionsAndPM(BzzVector &molar, BzzVector &mass, double MW, BzzVector &PM)
{
	int j;
	int NC = PM.Size();

	MW = 0.;
	for(j=1;j<=NC;j++)
		MW += molar[j]*PM[j];
		
	for(j=1;j<=NC;j++)
		mass[j] = molar[j] * PM[j] / MW;
}

void massFractionsAndPMtot(BzzVector &molar, BzzVector &mass, double MW, BzzVector &PM)
{
	int j;
	int NC = PM.Size();

	MW = 0.;
		
	for(j=1;j<=NC;j++)
		MW += molar[j]*PM[j];
		
	for(j=1;j<=NC;j++)
		mass[j] = molar[j] * PM[j] / MW;
}

/**
	
 * Ansi C "itoa" based on Kernighan & Ritchie's "Ansi C"
	
 * with slight modification to optimize for specific architecture:
	
 */
	
void strreverse(char* begin, char* end) {
	
	char aux;
	
	while(end>begin)
	
		aux=*end, *end--=*begin, *begin++=aux;
	
}

void my_itoa(int value, char* str, int base) {
	
	static char num[] = "0123456789abcdefghijklmnopqrstuvwxyz";
	
	char* wstr=str;
	
	int sign;
	
	div_t res;
	

	
	// Validate base
	
	if (base<2 || base>35){ *wstr='\0'; return; }
	

	
	// Take care of sign
	
	if ((sign=value) < 0) value = -value;
	

	
	// Conversion. Number is reversed.
	
	do {
	
		res = div(value,base);
	
		*wstr++ = num[res.rem];
	
	}while(value=res.quot);
	
	if(sign<0) *wstr++='-';
	
	*wstr='\0';
	

	
	// Reverse std::string
	
	strreverse(str,wstr-1);
	
}

void molarFractionsAndPM(BzzVector &molar, BzzVector &mass, double MW, BzzVector &PM)
{
	int j;
	int NC = PM.Size();

	MW = 0.;
	for(j=1;j<=NC;j++)
		MW += mass[j] / PM[j];

	MW = 1./MW;
	for(j=1;j<=NC;j++)
		molar[j] = mass[j] * MW / PM[j];
}


int myMaxInt(int a, int b)
{
	return (a>=b? a : b);
}

void gaussianFunction::setup(double _xCenter, double _w, double _reductionCoefficient, double _peak)
{
	xCenter = _xCenter;
	w = _w;
	reductionCoefficient = _reductionCoefficient;
	peak = _peak;

	mu = xCenter;
	sigma = 1.;

	double test;
	do
	{
		sigma *= 0.99;
		test = giveValue(xCenter+0.50*w);
	} while (test>(peak/reductionCoefficient));

}

double gaussianFunction::giveValue(double x)
{
	return peak*exp( -0.50*BzzPow2((x-mu)/sigma) );
}


void spark::setup(char *__shape, double __time, double __density, double __xCenter, double __width)
{
	if (!strcmp(__shape, "PLATE"))
		shape = 1;
	if (!strcmp(__shape, "GAUSSIAN"))
		shape = 2;
	
	timeSpark = __time * 1.e-3;
	density = __density * 1.e6;
	xCenter = __xCenter * 1.e-2;
	width = __width * 1.e-2;

	xLeft = xCenter - 0.50*width;
	xRight = xCenter + 0.50*width;

	if (shape == 2) // Gaussian profile
		gaussianProfile.setup(xCenter, width, 10., density);

}

double spark::giveEnergy(double t, double x)
{
	if (t < timeSpark)
	{
		if (shape==1) // Plate profile
			if( x>=xLeft && x<=xRight && t<timeSpark)
				return density*t/timeSpark;
			else	return 0.;

		else		 // Gaussian Profile
			return gaussianProfile.giveValue(x)*t/timeSpark;
	}
	else return 0.;
}

void openInputFileAndControl(ifstream &inputFile, const char *fileName)
{
	inputFile.open(fileName, ios::in);
	if (!inputFile.is_open()) 
	{
		cout << "The file could not be opened! " <<  fileName << endl; 
		cout << "Press enter to continue... "    << endl;
		getchar();
		exit(-1);
	}
}

void openInputFileAndControl_BinaryMode(ifstream &inputFile, const char *fileName)
{
	inputFile.open(fileName, ios::in | ios::binary);
	if (!inputFile.is_open()) 
	{
		cout << "The file could not be opened! " <<  fileName << endl; 
		cout << "Press enter to continue... "    << endl;
		getchar();
		exit(1);
	}
}

void openOutputFileAndControl(ofstream &outputFile, const char *fileName)
{
	outputFile.open(fileName, ios::out);
	if (!outputFile.is_open()) 
	{
		cout << "The file could not be opened! " <<  fileName << endl; 
		cout << "Press enter to continue... "    << endl;
		getchar();
		exit(1);
	}
}

void openOutputFileAndControl_BinaryMode(ofstream &outputFile, const char *fileName)
{
	outputFile.open(fileName, ios::out | ios::binary);
	if (!outputFile.is_open()) 
	{
		cout << "The file could not be opened! " <<  fileName << endl; 
		cout << "Press enter to continue... "    << endl;
		getchar();
		exit(1);
	}
}




void openInputFileAndControl(ifstream &inputFile, std::string fileName)
{
	inputFile.open(fileName.c_str(), ios::in);
	if (!inputFile.is_open()) 
	{
		cout << "The file could not be opened! " <<  fileName << endl; 
		getchar();
		exit(0);
	}
}

void openInputFileAndControl_BinaryMode(ifstream &inputFile, std::string fileName)
{
	inputFile.open(fileName.c_str(), ios::in | ios::binary);
	if (!inputFile.is_open()) 
	{
		cout << "The file could not be opened! " <<  fileName << endl; 
		getchar();
		exit(0);
	}
}

void openOutputFileAndControl(ofstream &outputFile, std::string fileName)
{
	outputFile.open(fileName.c_str(), ios::out);
	if (!outputFile.is_open()) 
	{
		cout << "The file could not be opened! " <<  fileName << endl; 
		getchar();
		exit(0);
	}
}

void openOutputFileAndControl_BinaryMode(ofstream &outputFile, std::string fileName)
{
	outputFile.open(fileName.c_str(), ios::out | ios::binary);
	if (!outputFile.is_open()) 
	{
		cout << "The file could not be opened! " <<  fileName << endl; 
		getchar();
		exit(0);
	}
}


void gridZones::setup(int _nPoints, double  _L, double _xCenter, double _wMix)
{
	if (_xCenter == 0.)	_xCenter 	= 0.50*_L;
	if (_wMix == 0.)	_wMix 		= 0.50*_L;	
	
	nPoints 	= _nPoints;
	L 		= _L;
	xCenter 	= _xCenter;
	wMix 		= _wMix;

	ChangeDimensions(nPoints, &x);

	xLeft	= xCenter-0.50*wMix;
	xRight	= xCenter+0.50*wMix;

	if (xLeft<0. || xRight>L)
	{
		cout << "Bad choice for xCen(" << xCenter << ") and wMix(" << wMix <<")!!" << endl;
		exit(1);
	}

	x[1]=0.; 
	double dxGrid = L/(nPoints-1);
	for (int i=2;i<=nPoints;i++)
		x[i]=x[i-1]+dxGrid;
}

void gridZones::setup(int _nPoints, double  _L, double _wMix)
{
	nPoints = _nPoints;
	L 		= _L;
	wMix	= _wMix;

	ChangeDimensions(nPoints, &x);

	xLeft	= L-0.75*wMix;
	xRight	= L-0.50*wMix;

	if (xLeft<0. || xRight>L)
	{
		cout << "Bad choice for xCen(" << xCenter << ") and wMix(" << wMix <<")!!" << endl;
		getchar();
		exit(-1);
	}

	x[1]=0.; 
	double dxGrid = L/(nPoints-1);
	for (int i=2;i<=nPoints;i++)
		x[i]=x[i-1]+dxGrid;
}

void functionRampa(LinearInterpolation &interpolatedProfile, gridZones &grid, double Left, double Right)
{
	int i;
	int nGrid = grid.nPoints;
	BzzVector profile(nGrid);


	for (i=1;i<=nGrid;i++)
		if (grid.x[i]<=grid.xLeft) profile[i] = Left;
		else if(grid.x[i]>=grid.xRight)  profile[i] = Right;
		else profile[i] = Left + (Right-Left)/grid.wMix*(grid.x[i]-grid.xLeft);

	interpolatedProfile(grid.x, profile);
}

void functionGaussian(LinearInterpolation &interpolatedProfile, gridZones &grid, double peak, double reductionFactor)
{
	int i;
	double a, b, c;
	int nGrid = grid.nPoints;
	BzzVector profile(nGrid);

	a = peak;
	b = grid.xCenter;
	c = sqrt(-BzzPow2(0.50*grid.wMix)/log(1./reductionFactor));

	for (i=1;i<=nGrid;i++)
		profile[i] = a*exp(-BzzPow2(grid.x[i]-b)/(c*c));

	interpolatedProfile(grid.x, profile);
}

void functionStep(LinearInterpolation &interpolatedProfile, gridZones &grid, double yLeft, double peak, double yRight)
{
	int i;
	int nGrid = grid.nPoints;
	BzzVector profile(nGrid);

	for (i=1;i<=nGrid;i++)
		if(grid.x[i]< grid.xLeft)		profile[i] = yLeft;
		else if(grid.x[i]> grid.xRight)	profile[i] = yRight;
		else							profile[i] = peak;

	interpolatedProfile(grid.x, profile);
}

void functionStepTwin(LinearInterpolation &interpolatedProfile, gridZones &grid, double yLeft, double peak)
{
	int i;
	int nGrid = grid.nPoints;
	BzzVector profile(nGrid);

	for (i=1;i<=nGrid;i++)
		if(grid.x[i]< grid.xLeft)		profile[i] = yLeft;
		else if(grid.x[i]> grid.xRight)	profile[i] = peak;
		else							profile[i] = peak;

	interpolatedProfile(grid.x, profile);
}

void FromCToBzzVector(int N, BzzVector &B, double *L)
{
	int i;
	for ( i= 0; i< N; i++ ) 
		B[i+1] = L[i];
}

void FromBzzVectorToC(int N, BzzVector &B, double *L)
{
	int i;
	for ( i= 0; i< N; i++ ) 
		L[i] = B[i+1];
}


void FromBzzMatrixToC(int NR, int NC, BzzMatrix &B, double *L)
{
	int i, j;
    for ( j= 0; j< NC; j++ )
		for ( i= 0; i< NR; i++ )
			L[i+j*NR] = B[i+1][j+1];
}

void FromCToBzzMatrix(int NR, int NC, BzzMatrix &B, double *L)
{
	int i, j;
    for ( j= 0; j< NC; j++ )
		for ( i= 0; i< NR; i++ )
			B[i+1][j+1] = L[i+j*NR]; 
}

void nameList::setup(int _N)
{
	N = _N;
	name = new char*[N+1];
	for(int i=1;i<=N;i++)
		name[i] = new char[20];
};

// **************************************************************************************************** //
//					OpenSMOKE logo                                                  //
// **************************************************************************************************** //


void OpenSMOKE_logo(std::string application_name, std::string version, std::string date)
{
	int i;
	std::string additional_info = "Version " + version + " - " + date;

    cout << endl;
	cout << "---------------------------------------------------------------------" << endl;
    cout << endl;
    cout << endl;
    cout << "       ____                    _____ __  __  ____  _  ________ " << endl;
    cout << "      / __ \\                  / ____|  \\/  |/ __ \\| |/ /  ____|" << endl;
    cout << "     | |  | |_ __   ___ _ __ | (___ | \\  / | |  | | ' /| |__   " << endl;
    cout << "     | |  | | '_ \\ / _ \\ '_ \\ \\___ \\| |\\/| | |  | |  < |  __|  " << endl;
    cout << "     | |__| | |_) |  __/ | | |____) | |  | | |__| | . \\| |____ " << endl;
    cout << "      \\____/| .__/ \\___|_| |_|_____/|_|  |_|\\____/|_|\\_\\______|" << endl;
    cout << "            | |                                                " << endl;
    cout << "            |_|                                                " << endl;
    cout << endl;
    cout << endl;
	cout << "---------------------------------------------------------------------" << endl;

	for (i=1;i<=35-int(application_name.size())/2;i++)	cout << " ";
	cout << application_name << endl;
	for (i=1;i<=35-int(additional_info.size())/2;i++)	cout << " ";
	cout << additional_info << endl;

	cout << "     Department of Chemistry, Materials and Chemical Engineering     " << endl;
	cout << "                        Politecnico di Milano                        " << endl;
	cout << "                       alberto.cuoci@polimi.it                       " << endl;
	cout << "                     tiziano.faravelli@polimi.it                     " << endl;
	cout << "                     alessio.frassoldati@polimi.it                   " << endl;
	cout << "                        eliseo.ranzi@polimi.it                       " << endl;
	cout << "---------------------------------------------------------------------" << endl;
	cout << "" << endl;
}

void OpenSMOKE_logo(ofstream &fOutput, std::string application_name)
{
	int i;

	fOutput << "---------------------------------------------------------------------" << endl;
    fOutput << endl;
    fOutput << endl;
    fOutput << "       ____                    _____ __  __  ____  _  ________ " << endl;
    fOutput << "      / __ \\                  / ____|  \\/  |/ __ \\| |/ /  ____|" << endl;
    fOutput << "     | |  | |_ __   ___ _ __ | (___ | \\  / | |  | | ' /| |__   " << endl;
    fOutput << "     | |  | | '_ \\ / _ \\ '_ \\ \\___ \\| |\\/| | |  | |  < |  __|  " << endl;
    fOutput << "     | |__| | |_) |  __/ | | |____) | |  | | |__| | . \\| |____ " << endl;
    fOutput << "      \\____/| .__/ \\___|_| |_|_____/|_|  |_|\\____/|_|\\_\\______|" << endl;
    fOutput << "            | |                                                " << endl;
    fOutput << "            |_|                                                " << endl;
    fOutput << endl;
    fOutput << endl;
	fOutput << "---------------------------------------------------------------------" << endl;

	for (i=1;i<=35-int(application_name.size())/2;i++)	fOutput << " ";
	fOutput << application_name << endl;
 
	fOutput << "     Department of Chemistry, Materials and Chemical Engineering     " << endl;
	fOutput << "                        Politecnico di Milano                        " << endl;
	fOutput << "                       alberto.cuoci@polimi.it                       " << endl;
	fOutput << "                     tiziano.faravelli@polimi.it                     " << endl;
	fOutput << "                     alessio.frassoldati@polimi.it                   " << endl;
	fOutput << "                        eliseo.ranzi@polimi.it                       " << endl;
	fOutput << "---------------------------------------------------------------------" << endl;
	fOutput << "" << endl;
}

// **************************************************************************************************** //
//					Spark Class                                                                         //
// **************************************************************************************************** //

enum { OPT_HELP, OPT_FLAG, OPT_ARG };

void ParserClass::first_screening(CSimpleOpt args)
{
    while (args.Next())
    {
        if (args.LastError() == SO_SUCCESS)
        {
            if (args.OptionId() == OPT_HELP)
            {
                ShowUsage();
                cout << "Press enter to continue..." << endl;
                getchar();
                exit(0);
            }
        }
        else
        {
            std::string message = args.OptionText();
            cout << "Invalid argument: " + message << endl;
            cout << "Press enter to continue..." << endl;
            getchar();
            exit(-1);
        }
    }
}

void ParserClass::setup(int _argc, char* _argv[], CSimpleOpt::SOption *_g_rgOptions)
{
    argc = _argc;
    for(int i=0;i<argc;i++)   argv[i] = _argv[i];
    g_rgOptions = _g_rgOptions;

    CSimpleOpt args(argc, argv, g_rgOptions);
    first_screening(args);
}

int ParserClass::parse(std::string label, std::string &argument)
{
    CSimpleOpt args(argc, argv, g_rgOptions);

    while (args.Next())
    {
        if (!strcmp(args.OptionText(), label.c_str() ))
        {
            argument = args.OptionArg();
            return 1;
        }
    }

    return 0;
}

int ParserClass::parse(std::string label)
{
    CSimpleOpt args(argc, argv, g_rgOptions);

    while (args.Next())
    {
        if (!strcmp(args.OptionText(), label.c_str() ))
            return 1;
    }

    return 0;
}

// **************************************************************************************************** //
//					Linear Interpolation                                                                //
// **************************************************************************************************** //

LinearInterpolation::LinearInterpolation()
{
	name_object = "Default name";
}

void LinearInterpolation::SetName(std::string name)
{
	name_object = name;
}

void LinearInterpolation::ErrorMessage(const std::string message)
{
    cout << endl;
    cout << "Class:  Linear Interpolation"		<< endl;
    cout << "Object: " << name_object			<< endl;
    cout << message                             << endl;
    cout << "Press a key to continue... "       << endl;
    getchar();
    exit(-1);
}

void LinearInterpolation::operator () (BzzVector &X, BzzVector &Y)
{
	n=X.Size();
	if(n!=Y.Size()) ErrorMessage("Dim(x) != Dim(y)");
	
	xMin=X[1];
	xMax=X[n];
	
	x = X;
	y = Y;
}

double LinearInterpolation::operator() (double t)
{
	int i;
	double epsilon = 1.e-4;
	
	if (t == xMin ) return y[1];
	else if (t == xMax ) return y[n];
	else if (t>(xMax+epsilon*xMax) || t<(xMin-epsilon*xMax)) 	
	{
		cout.setf(ios::scientific);
		cout << "Error in Linear Interpolation: Out of Boundaries!!!" << endl;
		cout << "t = " << t << endl;
		cout << "xMin = " << xMin << endl;
		cout << "xMax = " << xMax << endl;
		cout << "xMinAcc = " << xMin-epsilon*xMax << endl;
		cout << "xMaxAcc = " << xMax+epsilon*xMax << endl;
		getchar();
		exit(1);
	}

	else
	{
		for(i=2;i<=n;i++)
			if (t<=x[i]) break;
	
		return y[i-1]+(y[i]-y[i-1])/(x[i]-x[i-1])*(t-x[i-1]);
	}
}




OpenSMOKE_NuManager::OpenSMOKE_NuManager()
{
	nReactions		= 0;
}

void OpenSMOKE_NuManager::Set(const int index, const std::string name)
{
	index_species	= index;
	name_species	= name;
}

void OpenSMOKE_NuManager::Set(const double coefficient, const int reaction)
{
	nReactions++;
	nuReactions.Append(coefficient);
	iReactions.Append(reaction);
}

void OpenSMOKE_NuManager::PrintOnFile(ofstream &fOutput)
{
	fOutput << "Species:            " << name_species << " (" << index_species << ")" << endl;
	fOutput << "Reactions:          " << nReactions         << endl;
	fOutput << "---------------------------------------------------" << endl;
	for(int i=1;i<=nReactions;i++)
		fOutput << iReactions[i] << "\t" << nuReactions[i] << endl;
	fOutput << endl;
}

void OpenSMOKE_NuManager::Clean()
{
	BzzVectorInt index;

	for(int i=1;i<=nReactions;i++)
		if (iReactions[i]>0)
		{
			for(int j=i+1;j<=nReactions;j++)
				if (iReactions[j] == iReactions[i])
				{
					nuReactions[i] += nuReactions[j];
					iReactions[j]	= 0;
					nuReactions[j]  = 0;
					index.Append(j);
				}
		}

	iReactions.DeleteElements(index);
	nuReactions.DeleteElements(index);
	nReactions = iReactions.Size();
}

void OpenSMOKE_NuManager::SaveToBinaryFile(BzzSave &fSave)
{
	char name[Constants::NAME_SIZE];
	strcpy(name, name_species.c_str());
	fSave.fileSave.write((char*) name, sizeof(name));

	fSave  << nReactions;
	fSave  << nuReactions;
	fSave  << iReactions;
	fSave  << index_species;
}

void OpenSMOKE_NuManager::LoadFromBinaryFile(BzzLoad &fLoad)
{
	char dummy[Constants::NAME_SIZE];
	fLoad.fileLoad.read((char*) dummy, sizeof(dummy));
	name_species = dummy;

	fLoad >> nReactions;
	fLoad >> nuReactions;
	fLoad >> iReactions;
	fLoad >> index_species;
}

OpenSMOKE_Nu::OpenSMOKE_Nu(OpenSMOKE_ReactingGas *_mix)
{
	int i;

	mix	= _mix;
	NC  =  mix->NumberOfSpecies();
	NR  =  mix->NumberOfReactions();

	nu = new OpenSMOKE_NuManager[NC+1];

	for(i =1;i<=NC;i++)	
		nu[i].Set(i, mix->names[i]);

	BuildNuMatrix();

	for(i=1;i<=NC;i++)
		nu[i].Clean();

	Initialize();
}

OpenSMOKE_Nu::OpenSMOKE_Nu(BzzLoad &fLoad)
{
	fLoad >> NC;
	fLoad >> NR;

	nu = new OpenSMOKE_NuManager[NC+1];

	for(int i=1;i<=NC;i++)
		nu[i].LoadFromBinaryFile(fLoad);

	Initialize();
}

void OpenSMOKE_Nu::SaveToBinaryFile(BzzSave &fSave)
{
	fSave << NC;
	fSave << NR;

	for(int i=1;i<=NC;i++)
		nu[i].SaveToBinaryFile(fSave);
}

void OpenSMOKE_Nu::Initialize()
{
	ChangeDimensions(NC, &number_of_reactions);
	for(int i=1;i<=NC;i++)
		number_of_reactions[i] = nu[i].nReactions;
	
	SparsityPattern();
}

void OpenSMOKE_Nu::BuildNuMatrix()
{
	int i,k;

	int*	jD1 = mix->kinetics.jDir1.GetHandle();
	int*	jD2 = mix->kinetics.jDir2.GetHandle();
	int*	jD3 = mix->kinetics.jDir3.GetHandle();
	int*	jD4 = mix->kinetics.jDir4.GetHandle();
	int*	jD5 = mix->kinetics.jDir5.GetHandle();
	double* vD5 = mix->kinetics.valDir5.GetHandle();

	int*	jIT1 = mix->kinetics.jInvTot1.GetHandle();
	int*	jIT2 = mix->kinetics.jInvTot2.GetHandle();
	int*	jIT3 = mix->kinetics.jInvTot3.GetHandle();
	int*	jIT4 = mix->kinetics.jInvTot4.GetHandle();
	int*	jIT5 = mix->kinetics.jInvTot5.GetHandle();
	double* vIT5 = mix->kinetics.valInvTot5.GetHandle();

	for(i = 1;i <= NC;i++)
	{
		for(k = 1;k <= mix->kinetics.numDir1[i];k++)
			nu[i].Set(-1., *jD1++);
		for(k = 1;k <= mix->kinetics.numDir2[i];k++)
			nu[i].Set(-2., *jD2++);
		for(k = 1;k <= mix->kinetics.numDir3[i];k++)
			nu[i].Set(-3., *jD3++);
		for(k = 1;k <= mix->kinetics.numDir4[i];k++)
			nu[i].Set(-0.50, *jD4++);
		for(k = 1;k <= mix->kinetics.numDir5[i];k++)
			nu[i].Set(-(*vD5++), *jD5++);

		for(k = 1;k <= mix->kinetics.numInvTot1[i];k++)
			nu[i].Set(1., *jIT1++); 
		for(k = 1;k <= mix->kinetics.numInvTot2[i];k++)
			nu[i].Set(2., *jIT2++);
		for(k = 1;k <= mix->kinetics.numInvTot3[i];k++)
			nu[i].Set(3., *jIT3++);
		for(k = 1;k <= mix->kinetics.numInvTot4[i];k++)
			nu[i].Set(0.50, *jIT4++);
		for(k = 1;k <= mix->kinetics.numInvTot5[i];k++)
			nu[i].Set(*vIT5++, *jIT5++);
	}
}

void OpenSMOKE_Nu::SparsityPattern()
{
	const int MeanNumberOfSpeciesPerReaction = 8;
	const int max_size = MeanNumberOfSpeciesPerReaction*NR;

	ChangeDimensions(max_size, &sparse_rows);
	ChangeDimensions(max_size, &sparse_columns);
	ChangeDimensions(max_size, &sparse_index);

	int count=0;
	for(int index=1;index<=NR;index++)
		for(int j=1;j<=NC;j++)
		{
			for(int i=1;i<=nu[j].nReactions;i++)
				if (nu[j].iReactions[i] == index)
				{
					count++;
					sparse_rows[count]		= j;
					sparse_columns[count]	= i;
					sparse_index[count]		= index;
					break;
				}
		}

	sparse_rows.DeleteLastNElements(sparse_rows.Size()-count);
	sparse_columns.DeleteLastNElements(sparse_columns.Size()-count);
	sparse_index.DeleteLastNElements(sparse_index.Size()-count);
}

void OpenSMOKE_Nu::SparsityPattern(BzzMatrixSparse &M)
{
	ChangeDimensions(NR, NC, &M);
	for(int i=1;i<=sparse_index.Size();i++)
		M(sparse_index[i],  sparse_rows[i]) = 1.;
}

void OpenSMOKE_Nu::PrintOnFile(const std::string fileName)
{
	ofstream fOutput;
	openOutputFileAndControl(fOutput, fileName);
	fOutput.setf(ios::scientific);

	for(int i=1;i<=NC;i++)
		nu[i].PrintOnFile(fOutput);

	fOutput.close();
}


void OpenSMOKE_RateOfProductionCoefficient::Initialize(BzzVectorInt &n_)
{
	n = n_;
	p = new BzzVector[n.Size()+1];
	d = new BzzVector[n.Size()+1];
	for(int i=1;i<=n.Size();i++)
	{
		ChangeDimensions(n[i], &p[i]);
		ChangeDimensions(n[i], &d[i]);
	}
}

void OpenSMOKE_RateOfProductionCoefficient::Clean()
{
	for(int i=1;i<=n.Size();i++)
	{
		p[i] = 0.;
		d[i] = 0.;
	}
}

void OpenSMOKE_RateOfProductionCoefficient::LoadFromBinaryFile(BzzLoad &fLoad)
{
	fLoad >> n;
	p = new BzzVector[n.Size()+1];
	d = new BzzVector[n.Size()+1];
	for(int i=1;i<=n.Size();i++)
	{
		fLoad >> p[i];
		fLoad >> d[i];
	}
}

void OpenSMOKE_RateOfProductionCoefficient::SaveToBinaryFile(BzzSave &fSave)
{
	fSave << n;
	for(int i=1;i<=n.Size();i++)
	{
		fSave << p[i];
		fSave << d[i];
	}
}

void OpenSMOKE_RateOfProductionCoefficient::LoadFromBinaryFile(BzzLoad &fLoad, BzzVectorInt &indices)
{
	fLoad >> n;
	p = new BzzVector[n.Size()+1];
	d = new BzzVector[n.Size()+1];

	for(int i=1;i<=n.Size();i++)
	{
		fLoad >> p[i];
		fLoad >> d[i];
	}
}

void OpenSMOKE_RateOfProductionCoefficient::SaveToBinaryFile(BzzSave &fSave, BzzVectorInt &indices)
{
	BzzVectorInt aux(indices.Size());
	for (int j=1;j<=indices.Size();j++)
		aux = n[indices[j]];
	
	fSave << aux;
	for(int i=1;i<=indices.Size();i++)
	{
		fSave << p[indices[i]];
		fSave << d[indices[i]];
	}
}


std::string GiveMeFileNameFromFullPath(const std::string pathname)
{
	size_t found;
	found = pathname.find_last_of("/\\");
	return pathname.substr(found+1);
}


std::string GiveMeFolderPathFromFullPath(const std::string pathname)
{
	size_t found;
	found = pathname.find_last_of("/\\");
	return pathname.substr(0,found);
}

bool caseInsCharCompareN(char a, char b) 
{
	return(toupper(a) == toupper(b));
}

bool caseInsCompare(const std::string& s1, const std::string& s2) 
{
	return( (s1.size( ) == s2.size( )) && equal(s1.begin( ), s1.end( ), s2.begin( ), caseInsCharCompareN));
}

bool CheckForCommentLine(std::string &line)
{
	if (line.size() > 0)
		for(int j=0;j<int(line.size());j++)
		{
			if( line.at(j) == '!' )	
			{
				for(int k=j;k<int(line.size());k++)	line.at(k) = ' ';
				return true;
			}
		}
	return false;
}

bool CheckForCommentLineFromStart(std::string &line)
{
	if (line.size() > 0)
		for(int j=0;j<int(line.size());j++)
		{
			if( line.at(j) != ' ')	
				if( line.at(j) == '!')	return true;
				else					return false;
		}

	return false;
}

bool CheckForEndLine(std::string &line)
{
	if (line.size() > 0)
	{
		std::string dummy;
		stringstream parsed_string(line);
		parsed_string >> dummy;
		if (caseInsCompare(dummy,"END") == true) 
			return true;
	}
	return false;
}


bool CheckForBlankLine(const std::string line)
{
	if (line.size() == 0)		return true;
	for(int j=0;j<int(line.size());j++)
		if (line.at(j) != ' ')	return false;
	return true;
}

void CleanFromBlanks(std::string &name)
{
	StringSubstitutionAll(name, " ", "");
}

void CleanFromBlanksForThermo(std::string &name)
{
	std::string name_buffer = name;
	for(int j=0;j<int(name_buffer.size());j++)
		if (name_buffer.at(j) == ' ')	
		{
			name = name_buffer.substr(0,j);
			return;
		}
}

void StringSubstitution(std::string &original, const std::string old_string, const std::string new_string)
{
	std::string::size_type i = original.find(old_string);
	if (i != std::string::npos)
	{
		original.erase(i, old_string.length( ));
		original.insert(i, new_string);
	}
}

bool StringSubstitutionAll(std::string &s, const std::string p, const std::string new_string)
{
	bool found = false;
	if (p!=new_string)
	{
		std::string::size_type n = p.length( );
		for (std::string::size_type i = s.find(p); i != std::string::npos; i = s.find(p))
		{
			s.erase(i, n);
			s.insert(i, new_string);
			found = true;
		}
	}
	return found;
}

int StringFind(const std::string s, const std::string subs)
{
	int count = 0;
	std::string::size_type i = s.find(subs);
	if (i != std::string::npos)
		count++;
	return count;
}

bool StringFindSubString(std::string &s, const std::string p)
{
	std::string::size_type n = p.length( );
	for (std::string::size_type i = s.find(p); i != std::string::npos; i = s.find(p))
		return true;
	return false;
}

bool IsANumber(const std::string expression, const int index)
{
	if ((expression.at(index) == '0') ||
		(expression.at(index) == '1') ||
		(expression.at(index) == '2') ||
		(expression.at(index) == '3') ||
		(expression.at(index) == '4') ||
		(expression.at(index) == '5') ||
		(expression.at(index) == '6') ||
		(expression.at(index) == '7') ||
		(expression.at(index) == '8') ||
		(expression.at(index) == '9') ||
		(expression.at(index) == '-') ||
		(expression.at(index) == '.') 
	   )
	   return true;

	return false;
}

void SeparateNumberFromString(const std::string expression, std::string &name, double &number)
{
	int index = -1;
	bool stringFound = false;
	
	for(int i=0;i<int(expression.size());i++)
		if (IsANumber(expression, i) == false)
		{
			index = i;
			stringFound = true;
			break;
		}

	if (stringFound == false)
	{
		name = "";
		number = atof(expression.c_str());
	}
	else if (index == 0)
	{
		name = expression;
		number = 1.;
	}
	else
	{
		name = expression.substr(index, expression.size());
		number = atof(expression.substr(0, index).c_str());
	}
}

int SeparateNumberFromStringForElements(const std::string expression, std::string &name, double &number)
{
	int index = -1;
	bool stringFound = false;
	
	// Locating the first number in the std::string
	for(int i=0;i<expression.size();i++)
		if (IsANumber(expression, i) == true)
		{
			index = i;
			stringFound = true;
			break;
		}

	if (stringFound == false)
	{
		name = "";
		number = 0;
		return 3;
	}
	else
	{
	//	cout << "expression: $$$"<< expression<<"$$$" << endl;
	//	cout << "expression: $$$"<< expression.substr(index, expression.size()) <<"$$$" << endl;
	//	cout << index << endl;
		number = atof( expression.substr(index, expression.size()).c_str() );
		name   = expression.substr(0, index);
		CleanFromBlanks(name);

		return 3;
	}

/*	std::string name_buf = expression.substr(0,2);
	if (name_buf[1] == ' ')	
		name = name_buf[0];
	else						name = name_buf;
	number = atof(expression.substr(2, expression.size()).c_str());

	return 3;

	if (expression[0] != ' ')
	{
		if (expression[1] != ' ')
		std::string name_buf = expression.substr(0,2);
	if (name_buf[1] == ' ')	
		name = name_buf[0];
	else						name = name_buf;
	number = atof(expression.substr(2, expression.size()).c_str());

	}

		if(expression[1] == )*/
}

// ---------------------------------------------------------------
// Line parsing
// ---------------------------------------------------------------
void SeparateInstructions(const std::string line, vector<string> &instructions)
{
	instructions.resize(0);
	instructions.push_back("instructions");
		
	std::string dummy;
	stringstream parsed_string(line);
	
	for(;;)
	{
		parsed_string >> dummy;
		if (parsed_string.fail())	break;
		instructions.push_back(dummy);
	}
}

OpenSMOKE_ChebishevPolynomialsReaction::OpenSMOKE_ChebishevPolynomialsReaction(void) 
{
}

void OpenSMOKE_ChebishevPolynomialsReaction::ErrorMessage(const std::string message)
{
    cout << endl;
    cout << "Class:  OpenSMOKE_ChebishevPolynomialsReaction"	<< endl;
    cout << "Error:  " << message								<< endl;
    cout << "Press a key to continue... "						<< endl;
    getchar();
    exit(-1);
}

void OpenSMOKE_ChebishevPolynomialsReaction::WarningMessage(const std::string message)
{
    cout << endl;
    cout << "Class:  OpenSMOKE_ChebishevPolynomialsReaction"	<< endl;
    cout << "Warning:  "	<< message							<< endl;
    cout << "Press a key to continue... "						<< endl;
    getchar();
}

void OpenSMOKE_ChebishevPolynomialsReaction::ReadFromFile(BzzLoad &fInput)
{
	fInput >> N;
	fInput >> M;

	ChangeDimensions(N,   &phi_n);
	ChangeDimensions(M,   &phi_m);
	ChangeDimensions(N,M, &a);

	for(int n=1;n<=N;n++)
		for(int m=1;m<=M;m++)
			fInput >> a[n][m];

	fInput >> conversion;
	fInput >> Tmin >> Tmax;
	fInput >> Pmin >> Pmax;

	Pmin = OpenSMOKE_Conversions::conversion_pressure(Pmin, "atm");
	Pmax = OpenSMOKE_Conversions::conversion_pressure(Pmax, "atm");
	
	if (Tmin >= Tmax) ErrorMessage("Tmin >= Tmax");
	if (Pmin >= Pmax) ErrorMessage("Pmin >= Pmax");

	log10_Pmin = log10(Pmin);
	log10_Pmax = log10(Pmax);
}


double OpenSMOKE_ChebishevPolynomialsReaction::GiveMeKappa(const double T, const double P)
{
//	if (T>Tmax)	{ cout << T << " " << Tmax << endl; ErrorMessage("Temperature is larger than Tmax");}
//	if (T<Tmin)	{ cout << T << " " << Tmin << endl; ErrorMessage("Temperature is smaller than Tmin");}
//	if (P>Pmax)	{ cout << P << " " << Pmax << endl; ErrorMessage("Pressure is larger than Pmax");}
//	if (P<Pmin)	{ cout << P << " " << Pmin << endl; ErrorMessage("Pressure is smaller than Pmin");}

	double Ttilde = (2./T-1./Tmin-1./Tmax) / 
			           (1./Tmax-1./Tmin);
	double Ptilde = (2.*log10(P)-log10_Pmin-log10_Pmax) / 
			               (log10_Pmax-log10_Pmin);

	int n,m;
	for(n=1;n<=N;n++)	phi_n[n] = GiveMePhi(n, Ttilde);
	for(m=1;m<=M;m++)	phi_m[m] = GiveMePhi(m, Ptilde);

	double sum = 0.;
	for(n=1;n<=N;n++)
		for(m=1;m<=M;m++)	sum += a[n][m]*phi_n[n]*phi_m[m];

	return pow(10., sum)/conversion;
}

double OpenSMOKE_ChebishevPolynomialsReaction::GiveMePhi(const int n, const double x)
{
	return cos(double(n-1.)*acos(x));
}

void OpenSMOKE_ChebishevPolynomialsReaction::WriteCoefficientMatrix(ofstream &fOutput)
{
	fOutput << "   Chebishev Coefficient Matrix (" << N << "x" << M << ")" << endl;
	for(int n=1;n<=N;n++)
	{
		fOutput << "     ";
		for(int m=1;m<=M;m++)	
			fOutput << setw(15) << setprecision(5) << scientific << right << a[n][m];
		fOutput << endl;
	}
}


OpenSMOKE_LogarithmicPressureReaction::OpenSMOKE_LogarithmicPressureReaction(void) 
{
}

void OpenSMOKE_LogarithmicPressureReaction::ErrorMessage(const std::string message)
{
    cout << endl;
    cout << "Class:  OpenSMOKE_LogarithmicPressureReaction"	<< endl;
    cout << "Error:  " << message								<< endl;
    cout << "Press a key to continue... "						<< endl;
    getchar();
    exit(-1);
}

void OpenSMOKE_LogarithmicPressureReaction::WarningMessage(const std::string message)
{
    cout << endl;
    cout << "Class:  OpenSMOKE_LogarithmicPressureReaction"	<< endl;
    cout << "Warning:  "	<< message							<< endl;
    cout << "Press a key to continue... "						<< endl;
    getchar();
}

void OpenSMOKE_LogarithmicPressureReaction::ReadFromFile(BzzLoad &fInput)
{
	fInput >> N;

	ChangeDimensions(N,   &P);
	ChangeDimensions(N,   &lnP);
	ChangeDimensions(N,   &lnA);
	ChangeDimensions(N,   &Beta);
	ChangeDimensions(N,   &EoverR);

	int n;
	for(n=1;n<=N;n++)	fInput >> P[n];
	for(n=1;n<=N;n++)	fInput >> lnA[n];
	for(n=1;n<=N;n++)	fInput >> Beta[n];
	for(n=1;n<=N;n++)	fInput >> EoverR[n];
	
	for(n=1;n<=N;n++)
	{
		P[n] = OpenSMOKE_Conversions::conversion_pressure(P[n], "atm");
		lnP[n] = log(P[n]);
		lnA[n] = log(lnA[n]);
		EoverR[n] = EoverR[n]/Constants::R_cal_mol;
	}	
}


double OpenSMOKE_LogarithmicPressureReaction::GiveMeKappa(const double T, const double Pressure)
{
	if		(Pressure<=P[1])	return exp(lnA[1]+Beta[1]*log(T)-EoverR[1]/T);
	else if (Pressure>=P[N])	return exp(lnA[N]+Beta[N]*log(T)-EoverR[N]/T);
	else
	{
		int n;
		for(n=1;n<=N-1;n++)
			if (Pressure<P[n+1])
				break;
			
		double lnKappaA = lnA[n]  +Beta[n]  *log(T)-EoverR[n]/T;
		double lnKappaB = lnA[n+1]+Beta[n+1]*log(T)-EoverR[n+1]/T;
		return	exp( lnKappaA+(lnKappaB-lnKappaA)*(log(Pressure)-lnP[n])/(lnP[n+1]-lnP[n]) );
	}
}

void OpenSMOKE_LogarithmicPressureReaction::WriteCoefficientMatrix(ofstream &fOutput)
{
	double conversion = OpenSMOKE_Conversions::conversion_pressure(1., "atm");

	fOutput << "   Logarithmic Pressure Interpolation" << endl;
	fOutput << "     P[atm]      A[kmol,m3,s]   Beta      E[cal/mol]" << endl;
	for(int n=1;n<=N;n++)
		fOutput << "     " << setw(8)  << setprecision(3) << left  << fixed << P[n]/conversion 
						   << setw(15) << scientific	  << right			<< exp(lnA[n])
						   << setw(9)  << setprecision(4) << right << fixed << Beta[n] 
						   << setw(14) << setprecision(2) << right << fixed << EoverR[n]*Constants::R_cal_mol 
						   << endl;
}

void WriteUnits_cm_mol_s(ofstream &fOutput, const double lambda)
{
				if (lambda == 1.)	fOutput << "    " << "[1/s]";
		else	if (lambda == 2.)	fOutput << "    " << "[cm3/mol/s]";
		else	if (lambda == 3.)	fOutput << "    " << "[cm6/mol2/s]";
		else	if (lambda == 4.)	fOutput << "    " << "[cm9/mol3/s]";
		else	if (lambda == 5.)	fOutput << "    " << "[cm12/mol4/s]";
		else	if (lambda == 6.)	fOutput << "    " << "[cm15/mol5/s]";
		else	if (lambda == 7.)	fOutput << "    " << "[cm18/mol6/s]";
		else						fOutput << "    " << "[(cm3/mol)^" << setw(4) << setprecision(2) << fixed << lambda-1. <<"/s]";
		fOutput << endl;
}

void WriteUnits_m_kmol_s(ofstream &fOutput, const double lambda)
{
				if (lambda == 1.)	fOutput << "    " << "[1/s]";
		else	if (lambda == 2.)	fOutput << "    " << "[m3/kmol/s]";
		else	if (lambda == 3.)	fOutput << "    " << "[m6/kmol2/s]";
		else	if (lambda == 4.)	fOutput << "    " << "[m9/kmol3/s]";
		else	if (lambda == 5.)	fOutput << "    " << "[m12/kmol4/s]";
		else	if (lambda == 6.)	fOutput << "    " << "[m15/kmol5/s]";
		else	if (lambda == 7.)	fOutput << "    " << "[m18/kmol6/s]";
		else						fOutput << "    " << "[(m3/kmol)^" << setw(4) << setprecision(2) << fixed << lambda-1. <<"/s]";
		fOutput << endl;
}

double KineticTemperature(const double Cc, const double Tmean, const double Tatt, const double n)
{
	if (n == 0.)
	{
		double CcStar = Tatt/Tmean-log(Cc);
		return Tatt/CcStar;
	}
	else
	{
		double Told   = Tmean;
	
		for(int k=1;k<=20;k++)
		{
			double CcStar = n*log(Tmean)-Tatt/Tmean+log(Cc);
			double T = Tatt/(-CcStar);
			
			if (fabs(Told-T)/Tmean < 0.001)
				return T;
			
			Told = T;
		}

		return 300.;
	}
}

void TokenizeString(const std::string line, vector<string> &list)
{
	char *cstr1;
	strcpy(cstr1, line.c_str());
	char delim[]=" ,";
	char *token;

	//In the first call to strtok, the first argument is the line to be tokenized
	token=strtok(cstr1, delim);
	list.push_back(token);

	//In subsequent calls to strtok, the first argument is NULL
	while((token=strtok(NULL, delim))!=NULL)
		list.push_back(token);

	for(int j=1;j<=int(list.size());j++)
		CleanFromBlanks(list[j-1]);
}

void TokenizeString(char *line, vector<string> &list)
{
	char delim[]=" ,";
	char *token;

	//In the first call to strtok, the first argument is the line to be tokenized
	token=strtok(line, delim);
	list.push_back(token);
	

	//In subsequent calls to strtok, the first argument is NULL
	while((token=strtok(NULL, delim))!=NULL)
		list.push_back(token);

	std::string prov;
	for(int i=1;i<=int(list[0].size());i++)
		if (list[0].at(i-1) != ' ' && list[0].at(i-1) != '\t')
			prov.push_back(list[0].at(i-1));
	list[0] = prov;

	for(int j=1;j<=int(list.size());j++)
		CleanFromBlanks(list[j-1]);

}

std::string GetLabelIndex(const int count)
{
	stringstream number;
	number << count;
	std::string label = "(" + number.str() + ")";
	return label;
}

std::string GetNumber(const int count)
{
	stringstream number;
	number << count;
	return number.str();
}

void PrintTagOnGnuplotLabel(const int length, ofstream &fOut, const std::string tag, int &counter)
{
	stringstream number;
	number << counter++;
	std::string label = tag + "(" + number.str() + ")";
	fOut << setw(length) << left << label;
}

void PrintHeaderOnBinaryFile(BzzSave &fSave)
{
	std::string dummy;
	char name[Constants::NAME_SIZE];

	dummy = "V20100417";
	strcpy(name, dummy.c_str());
	fSave.fileSave.write((char*) name, sizeof(name));
	
	std::string building_date = GiveMeTimeAndDate();
	dummy = building_date;
	strcpy(name, dummy.c_str());
	fSave.fileSave.write((char*) name, sizeof(name));
}

void PrintEndOnBinaryFile(BzzSave &fSave)
{
	std::string dummy;
	char name[Constants::NAME_SIZE];

	dummy = "END";
	strcpy(name, dummy.c_str());
	fSave.fileSave.write((char*) name, sizeof(name));
}

void CheckInBinaryFile(BzzLoad &fLoad, const std::string flag)
{
	char dummy[Constants::NAME_SIZE];
	fLoad.fileLoad.read((char*) dummy, sizeof(dummy));
	if (strcmp(dummy, flag.c_str()))	
	{
		cout << endl;
		cout << "File:   " << fLoad.GetFileName()		<< endl;
		cout << "Error:  " << "Expected: " << flag << " - Found: " << std::string(dummy)	<< endl;
		cout << "Press a key to continue... "		<< endl;
		getchar();
		exit(-1);
	}
}

std::string NextInBinaryFile(BzzLoad &fLoad)
{
	char dummy[Constants::NAME_SIZE];
	fLoad.fileLoad.read((char*) dummy, sizeof(dummy));
	return dummy;
}

void PrepareLiquidPropertiesDictionary(const std::string tag, const std::string file_name, vector<string> &lines, BzzVectorInt &indexLines)
{
	const int SIZE = 400;
	char comment[SIZE];

	int number_of_lines;
	int total_number_of_species;

	BzzVectorInt indexBlankLines;
	BzzVectorInt indexCommentLines;


	ifstream fInput;
	openInputFileAndControl(fInput, file_name);

	// ---------------------------------------------------------------
	// Reading lines
	// ---------------------------------------------------------------
	lines.push_back("List of lines");

	while(!fInput.eof())
	{
		fInput.getline(comment, SIZE);
		lines.push_back(comment);
	}
	fInput.close();

	number_of_lines = lines.size()-1;

	// ---------------------------------------------------------------
	// Parsing lines
	// ---------------------------------------------------------------
	int i;
	for(i=1;i<=number_of_lines;i++)
	{
			 if (CheckForBlankLine(lines[i])			== true)	indexBlankLines.Append(i);
		else if (CheckForCommentLineFromStart(lines[i])	== true)	indexCommentLines.Append(i);
		else if (CheckForEndLine(lines[i])				== true)	indexCommentLines.Append(i);		
		else
		{
			CheckForCommentLine(lines[i]);
			indexLines.Append(i);
		}
	}

	total_number_of_species = indexLines.Size();

	cout << endl;
	cout << " ----------------------------------------------------------------" << endl;
	cout << "                    " << tag										<< endl;
	cout << " ----------------------------------------------------------------" << endl;
	cout << "    Total number of full lines:    " << indexLines.Size()			<< endl;
	cout << "    Total number of blank lines:   " << indexBlankLines.Size()	<< endl;
	cout << "    Total number of comment lines: " << indexCommentLines.Size()	<< endl;
	cout << "    Total number of lines:         " << number_of_lines			<< endl;
	cout << "    Total number of species:       " << total_number_of_species	<< endl;
	cout << " ----------------------------------------------------------------" << endl;
	cout << endl;
}

std::string StringToUpper(std::string strToConvert)
{
	//change each element of the std::string to upper case
	for(unsigned int i=0;i<strToConvert.length();i++)
	{
		  strToConvert[i] = toupper(strToConvert[i]);
	}
	return strToConvert;//return the converted std::string
}

std::string StringToLower(std::string strToConvert)
{
	//change each element of the std::string to lower case
	for(unsigned int i=0;i<strToConvert.length();i++)
	{
		strToConvert[i] = tolower(strToConvert[i]);
	}
	return strToConvert;//return the converted std::string
}