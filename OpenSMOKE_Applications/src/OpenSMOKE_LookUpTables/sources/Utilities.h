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

#if !defined(MY_UTILITIES)
#define MY_UTILITIES

#include "basic/OpenSMOKE_Utilities.h"

class flamelet_group;

void checkingWords(const string foundWord, const string expectedWord, const string fileName);
void extract_name(char *word, int &compositionLabel, string &name);

void import_flamelet_file_from_FLUENT(OpenSMOKE_ReactingGas &mix, const string fileName, flamelet_group &flameletGroup);
void conversion_from_FlameCRECK_to_FLUENT(OpenSMOKE_ReactingGas &mix, string fileNameSource, ofstream &fOutput, bool verboseFormationRates);

void search_index_for_linear_interpolation(BzzVector &x, double mean, int &iMean, double &interpolationFactor);
double interpolate_function(BzzVector &v, int iMean, double interpolationFactor);
void double_the_vector(BzzVector &original, BzzVector &doublevector);
void double_the_vector(int nRefinements, BzzVector &original, BzzVector &doublevector);

string lowercase(string command); // converts all names input to lower case
string uppercase(string command); // converts all names input to upper case

void read_flamelet_library(char *fileName, flamelet_group &flameletGroup);

void give_K_Planck_SANDIA(double T, double &kH2O, double &kCO2, double &kCO, double &kCH4);

void read_SootProperties_from_FlameCRECK(const string fileNameSource, BzzVector &m0, BzzVector &fv);


// Extraction of data from the library
void extract_data_2d(char* source_file_name, char* dest_file_name, int xcolumn, int ycolumn);
void extract_data_2d(char* source_file_name, char* dest_file_name, int xcolumn);





class species_list
{
public:

	vector<string> names;
	BzzVector pm;

	int iO2;
	int iOH;
	int iO;
	int iC2H2;
	int iH2;
	int iC6H6;
	int iC6H5;
};

class log_normal_distribution
{
	double chiMean;
	int N;
	BzzVector chi;
	BzzVector teta;
	BzzVector error_function;
	BzzVector integral;

public:
	void	construct_sequence(double chiMean, BzzVector &discreteChi);
	double	apply(BzzVector &Beta);
};

#endif	// MY_UTILITIES