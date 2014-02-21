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

#include "sources_in_flamelet_library.h"
#include "Utilities.h"


void sources_in_flamelet_library::initialize(int _Nslices_MixtureFraction, int _Nslices_Variance)
{
	Nslices_MixtureFraction = _Nslices_MixtureFraction;
	Nslices_Variance = _Nslices_Variance;
	ChangeDimensions(Nslices_MixtureFraction, &source);
	ChangeDimensions(Nslices_MixtureFraction, Nslices_Variance, &source_mean);
	ChangeDimensions(Nslices_MixtureFraction, Nslices_Variance, &source_mean_nofluctuations);
}

void sources_in_flamelet_library::print_on_file(ofstream &fOut, char *fileName, BzzVector &csi, BzzVector &v, double chi, const double pressure_Pa)
{
	// Writing on file - Formatted File for Visualization and Post Processing
	{
		ofstream fOutput(fileName, ios::out);
		fOutput.setf(ios::scientific);
		for(int k=1;k<=Nslices_MixtureFraction;k++)	
		{
			for(int i=1;i<=Nslices_Variance;i++)
			{
				fOutput << csi[k]							<< "\t";
				fOutput << v[i]								<< "\t";
				
				if (source_mean[k][i]<1.e-48)
					fOutput << 1.e-48								<< "\t";
				else fOutput << source_mean[k][i]					<< "\t";
				if (source_mean_nofluctuations[k][i]<1.e-48)
					fOutput << 1.e-48								<< "\t";
				else fOutput << source_mean_nofluctuations[k][i]	<< "\t";

				fOutput << endl;
			}
			fOutput << endl;
		}
		fOutput.close();
	}

	int k, i;

	// Writing on File - Unformatted File for FLUENT
	fOut << "pressure "		<< pressure_Pa << endl;
	fOut << "chi "			<< chi << endl;
	fOut << "nCsi "			<< Nslices_MixtureFraction << endl;
	fOut << "nVariance "	<< Nslices_Variance << endl;
		
	// Mixture Fraction
	fOut << "mf" << endl;
	for(k=1;k<=Nslices_MixtureFraction;k++)	
		fOut << csi[k] << endl;

	// Variance
	fOut << "mfv" << endl;
	for(i=1;i<=Nslices_Variance;i++)
		fOut << v[i] << endl;

	// Source Fields
	fOut << "sources" << endl;
	for(k=1;k<=Nslices_MixtureFraction;k++)	
	{
		for(i=1;i<=Nslices_Variance;i++)
		{
			fOut << source_mean[k][i] << "\t";
			fOut << source_mean_nofluctuations[k][i] << "\t";
			fOut << endl;
		}
		fOut << endl;
	}
	fOut << endl;	
}


void sources_in_flamelet_library::assign_source(BzzVector &_source)
{ 
	source = _source;
}

void sources_in_flamelet_library::assign_source_mean_nofluctuations(BzzVector &_source_mean_nofluctuations)
{
	source_mean_nofluctuations = _source_mean_nofluctuations;
}


void sourceField_flamelet_library::read_sourceField_fromFile(const char* fileName)
{
	int i,j, k;
	char word[10];
	ifstream fInput;
	
	fInput.open(fileName, ios::in);

	fInput >> word;
	checkingWords(word, "nFlamelets", fileName);
	fInput >> nChi;

	ChangeDimensions(nChi, &pressure_Pa);
	ChangeDimensions(nChi, &chi);
	ChangeDimensions(nChi, &max_normal_variance);
	ChangeDimensions(nChi, &nCsi);
	ChangeDimensions(nChi, &nVariance);

	mixture_fraction	= new BzzVector[nChi+1];
	normal_variance		= new BzzVector[nChi+1];
	source_pdf			= new BzzMatrix[nChi+1];
	source_nopdf		= new BzzMatrix[nChi+1];

	for(k=1;k<=nChi;k++)
	{
		fInput >> word;
		checkingWords(word, "pressure", fileName);
		fInput >> pressure_Pa[k];

		fInput >> word;
		checkingWords(word, "chi", fileName);
		fInput >> chi[k];

		fInput >> word;
		checkingWords(word, "nCsi", fileName);
		fInput >> nCsi[k];

		fInput >> word;
		checkingWords(word, "nVariance", fileName);
		fInput >> nVariance[k];

		ChangeDimensions(nCsi[k], &mixture_fraction[k]);
		ChangeDimensions(nVariance[k], &normal_variance[k]);
		ChangeDimensions(nCsi[k], nVariance[k], &source_pdf[k]);
		ChangeDimensions(nCsi[k], nVariance[k], &source_nopdf[k]);

		fInput >> word;
		checkingWords(word, "mf", fileName);
		for(i=1;i<=nCsi[k];i++)
			fInput >> mixture_fraction[k][i];

		fInput >> word;
		checkingWords(word, "mfv", fileName);
		for(j=1;j<=nVariance[k];j++)
			fInput >> normal_variance[k][j];

		fInput >> word;
		checkingWords(word, "sources", fileName);
		for(i=1;i<=nCsi[k];i++)
			for(j=1;j<=nVariance[k];j++)
			{
				fInput >> source_pdf[k][i][j];
				fInput >> source_nopdf[k][i][j];
			}
		
		max_normal_variance[k] = normal_variance[k][nVariance[k]];
	}

	for(k=2;k<=nChi;k++)
		if (chi[k] <= chi[k-1])
		{
			cout << "ERROR in reading file: " << fileName << endl;
			cout << "Please check the order of the scalar dissipation rates!" << endl;
		}

	minimum_chi = chi[1];
	maximum_chi = chi[nChi];
}

double sourceField_flamelet_library::find_interpolated_value_local(int iChi, double csiMean, double normalVariance)
{
	double interp_A, interp_B;
	int iCsi, iVariance;
	double factorCsi, factorVariance;

	// Checking on normal variance
	if (normalVariance >= max_normal_variance[iChi])
	{
		cout << "WARNING: normal variance over the maximum value. Correction apported!";
		cout << " Actual: " << normalVariance;
		cout << " Max:    " << max_normal_variance[iChi] << endl;
		normalVariance = max_normal_variance[iChi]-1.e-6;
	}

	// Checking on the maximum value of the mixture fraction
	if (csiMean == 1.0)
		csiMean = 1.-1.e-10;

	// Interpolated Value
	search_index_for_linear_interpolation(mixture_fraction[iChi], csiMean, iCsi, factorCsi);
	search_index_for_linear_interpolation(normal_variance[iChi], normalVariance, iVariance, factorVariance);	

	// TODO Diversificare per pdf e no pdf
//	interp_A = interpolate_function(source_nopdf[iChi].GetRow(iCsi-1), iVariance, factorVariance);
//	interp_B = interpolate_function(source_nopdf[iChi].GetRow(iCsi), iVariance, factorVariance);
	
	BzzVector interp_a_vector = source_pdf[iChi].GetRow(iCsi-1);
	BzzVector interp_b_vector = source_pdf[iChi].GetRow(iCsi);

	interp_A = interpolate_function(interp_a_vector, iVariance, factorVariance);
	interp_B = interpolate_function(interp_b_vector, iVariance, factorVariance);

	return ( interp_A + (interp_B-interp_A)*factorCsi );
}


double sourceField_flamelet_library::find_interpolated_value_global(double csiMean, double normalVariance, double scalarDissipationRate)
{
	int iChi;
	double factorChi;

	if (nChi == 1)
		return find_interpolated_value_local(1, csiMean, normalVariance);
	
	// Search for the scalar dissipation rate
	else
	{
		// Checking on the scalar dissipation rate
		if (scalarDissipationRate < minimum_chi) scalarDissipationRate = minimum_chi;
		if (scalarDissipationRate > maximum_chi) scalarDissipationRate = maximum_chi;

		search_index_for_linear_interpolation(chi, scalarDissipationRate, iChi, factorChi);
		double interp_A = find_interpolated_value_local(iChi-1, csiMean, normalVariance);
		double interp_B = find_interpolated_value_local(iChi, csiMean, normalVariance);
		return ( interp_A + (interp_B-interp_A)*factorChi );
	}
}

double sourceField_flamelet_library::find_interpolated_value_global_lognormal(double csiMean, double normalVariance, double scalarDissipationRate)
{
	BzzVector interpolated_value(nChi);
	log_normal_distribution logpdf;
		
	logpdf.construct_sequence(scalarDissipationRate, chi);

	if (nChi < 4)
	{
		cout << "The number of scalr dissipation rate slices is too low for applying the log-normal distribution!" << endl;
		return 0.;
	}

	// Search for the scalar dissipation rate
	else	
	{	
		for(int j=1;j<=nChi;j++)
			interpolated_value[j] = find_interpolated_value_local(j, csiMean, normalVariance);
		return logpdf.apply(interpolated_value);
	}

	// TEST FOR LOG NORMAL DISTRIBUTION
	/*
	BzzVector test(10000);
	test[1]=1.e-12;
	test[2]=1.e-11;
	test[3]=1.e-10;
	test[4]=1.e-9;
	test[5]=1.e-8;
	test[6]=1.e-7;
	test[7]=1.e-6;
	test[8]=1.e-5;
	test[9]=1.e-4;
	test[10]=1.e-3;
	test[11]=1.e-2;
	test[12]=1.e-1;
	test[13]=2.e-1;
	test[14]=3.e-1;
	for(int i=15;i<=5000;i++)
		test[i] = test[i-1]+1.e-1;
	for(i=5001;i<=9000;i++)
		test[i] = test[i-1]+1.e0;
	for(i=9001;i<=10000;i++)
		test[i] = test[i-1]+1.e1;
	logpdf.construct_sequence(scalarDissipationRate, test);
	return 0;
	*/

}

double sourceField_flamelet_library::find_interpolated_value_global_lognormal_backup(double csiMean, double normalVariance, double scalarDissipationRate)
{
	BzzVector interpolated_value(nChi);
	log_normal_distribution logpdf;
		
	logpdf.construct_sequence(scalarDissipationRate, chi);

	if (nChi < 4)
	{
		cout << "The number of scalar dissipation rate slices is too low for applying the log-normal distribution!" << endl;
		return 0.;
	}

	// Search for the scalar dissipation rate
	else	
	{	
		cout.setf(ios::scientific);
		cout << "#Flame\t" << "Chi[Hz]\t" << "Source" << endl;
		for(int j=1;j<=nChi;j++)
		{
			interpolated_value[j] = find_interpolated_value_local(j, csiMean, normalVariance);
			cout << j << "\t" << chi[j] << "\t" << interpolated_value[j] << endl;
		}
		return logpdf.apply(interpolated_value);
	}
}
