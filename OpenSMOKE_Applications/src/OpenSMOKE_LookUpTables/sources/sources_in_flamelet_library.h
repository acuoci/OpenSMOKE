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

#if !defined(SOURCES_IN_FLAMELET_LIBRARY)
#define SOURCES_IN_FLAMELET_LIBRARY

#include <string>
#include <iostream>
#include <fstream>
#include "BzzMath.hpp"

class sources_in_flamelet_library
{
public:
	BzzVector source;
	BzzMatrix source_mean;
	BzzMatrix source_mean_nofluctuations;
	int Nslices_MixtureFraction;
	int Nslices_Variance;

	void initialize(int _Nslices_MixtureFraction, int _Nslices_Variance);
	void print_on_file(ofstream &fOut, char *fileName, BzzVector &csi, BzzVector &v, double chi, const double pressure_Pa);
	void assign_source(BzzVector &_source);
	void assign_source_mean_nofluctuations(BzzVector &_source_mean_nofluctuations);

private:

};

class sourceField_flamelet_library
{
public:
	int nChi;
	double minimum_chi;
	double maximum_chi;

	BzzVector pressure_Pa;
	BzzVector chi;
	BzzVector max_normal_variance;
	BzzVectorInt nCsi;
	BzzVectorInt nVariance;

	BzzVector *mixture_fraction;
	BzzVector *normal_variance;
	BzzMatrix *source_pdf;
	BzzMatrix *source_nopdf;

	void read_sourceField_fromFile(const char* fileName);

	double find_interpolated_value_local(int iChi, double csiMean, double normalVariance);
	double find_interpolated_value_global(double csiMean, double normalVariance, double scalarDissipationRate);
	double find_interpolated_value_global_lognormal(double csiMean, double normalVariance, double scalarDissipationRate);

	double find_interpolated_value_global_lognormal_backup(double csiMean, double normalVariance, double scalarDissipationRate);
};


#endif // !defined(SOURCES_IN_FLAMELET_LIBRARY)