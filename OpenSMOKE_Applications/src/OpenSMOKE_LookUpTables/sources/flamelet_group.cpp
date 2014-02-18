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

#include <sstream>
#include "flamelet_group.h"
#include "basic/OpenSMOKE_Utilities.h"
#include "engine/OpenSMOKE_ReactingGas.h"

void flamelet_group::ErrorMessage(const string message)
{
    cout << endl;
    cout << "Class:  OpenSMOKE_LookUp_Table_FlameletGroup"	<< endl;
    cout << "Object: " << name_object		<< endl;
    cout << "Error:  " << message           << endl;
    cout << "Press a key to continue... "   << endl;
    getchar();
    exit(-1);
}

void flamelet_group::WarningMessage(const string message)
{
    cout << endl;
    cout << "Class:  OpenSMOKE_LookUp_Table_FlameletGroup"	<< endl;
    cout << "Object: "		<< name_object	<< endl;
    cout << "Warning:  "	<< message      << endl;
    cout << "Press a key to continue... "   << endl;
    getchar();
}

void flamelet_group::initialize()
{
	nFlamelets = 0;
	flames = new flamelet[1];
}

void flamelet_group::append(flamelet newflame)
{
	int j;
	flamelet *flames_buf;
	
	flames_buf = new flamelet[nFlamelets+1];
	for(j=1;j<=nFlamelets;j++)
		flames_buf[j] = flames[j];

	nFlamelets++;
	flames = new flamelet[nFlamelets+1];
	for(j=1;j<=nFlamelets-1;j++)
		flames[j] = flames_buf[j];
	flames[nFlamelets] = newflame; 
}

void flamelet_group::reorder()
{
	int j;
	BzzVector chi(nFlamelets);
	BzzVectorInt ichi(nFlamelets);
	for(j=1;j<=nFlamelets;j++)
		chi[j] = flames[j].chi;
	Sort(&chi, &ichi);

	flamelet *flames_buf;
	
	flames_buf = new flamelet[nFlamelets+1];
	for(j=1;j<=nFlamelets;j++)
		flames_buf[j] = flames[j];

	flames = new flamelet[nFlamelets+1];
	for(j=1;j<=nFlamelets;j++)
		flames[j] = flames_buf[ichi[j]];
}


void flamelet_group::lock()
{
	int i,j;

	for(j=2;j<=nFlamelets;j++)
		if (flames[j].NC != flames[j-1].NC)
			ErrorMessage("The number of species must be the same for all the flamelets...");

	for(j=2;j<=nFlamelets;j++)
		for(i=1;i<=flames[1].NC;i++)
		if (flames[1].species.names[i] != flames[j].species.names[i])
			ErrorMessage("The order of species must be the same for all the flamelets...");
}


void flamelet_group::print_on_file(char *fileName)
{
	int j;
	ofstream fOutput;
	fOutput.open(fileName, ios::out);
	fOutput.setf(ios::scientific);

	fOutput << "Chi [Hz]" << "\t";
	fOutput << "Csi [-] " << "\t";
	fOutput << "T [K]   " << "\t";
	int count = 4; 
	for(j=1;j<=flames[1].NC;j++)
		fOutput << flames[1].species.names[j].c_str() << "\t[" << count++ << "]\t";
	fOutput << endl << endl;

	for(j=1;j<=nFlamelets;j++)
		flames[j].print_on_file_appending(fOutput);

	fOutput.close();
}

void flamelet_group::build_library(OpenSMOKE_ReactingGas &mix, const string option, const string file_name, int NVariance, double alfa)
{
	int j;
	char name[40];
	char number_char[3];

	ofstream fLibrary;
	openOutputFileAndControl(fLibrary, file_name);
	fLibrary.setf(ios::scientific);

	// Intestazione file library
	fLibrary << "nChi " << nFlamelets << endl << endl;

	// Beta PDF
	if (option == "BETA")
	{
		for (j=1;j<=nFlamelets;j++)
		{
			cout << endl;
			cout << "Flamelet no. " << j << " - Scalar dissipation rate (st.) = " << flames[j].chi << " Hz" << endl;
			cout << "----------------------------------------------------------" << endl;

			flames[j].prepare_library(NVariance, alfa);
			flames[j].apply_betaPDF(mix);

			strcpy(name, "Output/BetaPDF_");
			my_itoa(j, number_char, 10);
			strcat(name, number_char);
			strcat(name, ".out");
			flames[j].print_library(name);
			flames[j].print_library(fLibrary);
		}
	}


	// Clipped Gaussian PDF
	else if (option == "GAUSSIAN")
	{
		for (j=1;j<=nFlamelets;j++)
		{
			cout << endl;
			cout << "Flamelet no. " << j << " - Scalar dissipation rate (st.) = " << flames[j].chi << " Hz" << endl;
			cout << "----------------------------------------------------------" << endl;

			flames[j].prepare_library(NVariance, alfa);
			flames[j].apply_clippedGaussian();

			strcpy(name, "Output/ClippedGaussianPDF_");
			my_itoa(j, number_char, 10);
			strcat(name, number_char);
			strcat(name, ".out");
			flames[j].print_library(name);
			flames[j].print_library(fLibrary);
		}
	}

	fLibrary.close();
}

void flamelet_group::write_progress_variables(OpenSMOKE_ReactingGas &mix, const int nProgressVariable, vector<string>& progress_variable_name, BzzVectorInt* progress_variable_index, BzzVector* progress_variable_value)
{
	for (int j=1;j<=nFlamelets;j++)
		flames[j].write_progress_variables(mix, nProgressVariable, progress_variable_name, progress_variable_index, progress_variable_value);
}

void flamelet_group::apply_source_betaPDF(OpenSMOKE_ReactingGas &mix, bool iPerfectlyCorrelatedApproach, int nSources, 	BzzVector *m0_N, 	BzzVector *fv_N,
										  nucleation_models nucleation_model, growth_models growth_model, aggregation_models aggregation_model, oxidation_models oxidation_model)
{
	int i,j;
	string *file_names;
	ofstream *file_list;

	file_list	= new ofstream[nSources+1];
	file_names	= new string[nSources+1];

	string nucleation_model_string	= "None";
	string growth_model_string		= "None";
	string aggregation_model_string	= "None";
	string oxidation_model_string	= "None";

	     if (nucleation_model == NUCLEATION_LIU_2001) 		nucleation_model_string = "Liu2001";
	else if (nucleation_model == NUCLEATION_LIU_2002) 		nucleation_model_string = "Liu2002";
	else if (nucleation_model == NUCLEATION_LIU_2003) 		nucleation_model_string = "Liu2003";
	else if (nucleation_model == NUCLEATION_MOSS_1999) 		nucleation_model_string = "Moss1999";
	else if (nucleation_model == NUCLEATION_WEN_2003) 		nucleation_model_string = "Wen2003";
	else if (nucleation_model == NUCLEATION_LINDSTEDT_1994) nucleation_model_string = "Lindstedt1994";
	else if (nucleation_model == NUCLEATION_LEUNG_1991) 	nucleation_model_string = "Leung1991";											
	else if (nucleation_model == NUCLEATION_HALL_1997) 		nucleation_model_string = "Hall1997";											

	     if (growth_model == GROWTH_LIU_2001) 				growth_model_string = "Liu2001";
	else if (growth_model == GROWTH_LIU_2002) 				growth_model_string = "Liu2002";
	else if (growth_model == GROWTH_LIU_2003) 				growth_model_string = "Liu2003";
	else if (growth_model == GROWTH_MOSS_1999) 				growth_model_string = "Moss1999";
	else if (growth_model == GROWTH_WEN_2003) 				growth_model_string = "Wen2003";
	else if (growth_model == GROWTH_LINDSTEDT_1994) 		growth_model_string = "Lindstedt1994";
	else if (growth_model == GROWTH_LEUNG_1991) 			growth_model_string = "Leung1991";

	     if (oxidation_model == OXIDATION_LEE)				oxidation_model_string = "Lee";
	else if (oxidation_model == OXIDATION_NSC)				oxidation_model_string = "NSC";
	else if (oxidation_model == OXIDATION_NEOH)				oxidation_model_string = "Neoh";

	     if (aggregation_model == AGGREGATION_SMOLUCHOWSKI)	aggregation_model_string = "Smoluchowski";
	else if (aggregation_model == AGGREGATION_MOSS)			aggregation_model_string = "Moss";

	file_names[1] = "Nucleation_" + nucleation_model_string;
	file_names[2] = "Growth_" + growth_model_string;
	file_names[3] = "Aggregation_" + aggregation_model_string;
	file_names[4] = "Oxidation_" + oxidation_model_string;

	for (i=1;i<=nSources;i++)
	{
		char file_name[60];
		system("mkdir Sources");
		strcpy(file_name, "Sources/");
		strcat(file_name, file_names[i].c_str());
		strcat(file_name, ".source");
		openOutputFileAndControl(file_list[i], file_name);
		file_list[i].setf(ios::scientific);

		// Intestazione del file
		file_list[i] << "nFlamelets " << nFlamelets << endl << endl;
	}
	
	// Costruzione della libreria per i termini sorgenti
	for (j=1;j<=nFlamelets;j++)
	{	
		cout << "Flamelet no. " << j << " - " << "Scalar dissipation rate (st.): " << flames[j].chi << " Hz" << endl;
		cout << "------------------------------------------------------------------------" << endl;
		flames[j].apply_source_betaPDF(mix, iPerfectlyCorrelatedApproach, j, nSources, 
			                           file_list, file_names, m0_N[j], fv_N[j],
									   nucleation_model, growth_model, aggregation_model, oxidation_model);
	}

	// Chiusura dei files
	for (i=1;i<=nSources;i++)
		file_list[i].close();
}

void flamelet_group::read_flamelet_library(OpenSMOKE_ReactingGas &mix, const string fileName)
{
	int i,j,k;
	char word[30];
	vector<string> names;
	
	ifstream fInput;
	openInputFileAndControl(fInput, fileName);

	int number_of_points_mixture_fraction;
	int	number_of_points_variance;
	int	number_of_species;
	double scalar_dissipation_rate;
	double pressure_pa;

	BzzVector csi;
	BzzVector variance;
	BzzMatrix density;
	BzzMatrix temperature;
	BzzMatrix composition;
	BzzMatrix *w;
	BzzMatrix z_r;
	BzzMatrix zv2_r;

	fInput >> word;
	checkingWords(word, "nChi", fileName);
	fInput >> nFlamelets;

	flames = new flamelet[nFlamelets+1];

	cout << "Reading Flamelets... " << endl;
	cout << "-----------------------------------------------" << endl;
	for (int index = 1; index<=nFlamelets; index++)
	{
		fInput >> word;
		checkingWords(word, "pressure", fileName);
		fInput >> pressure_pa;

		fInput >> word;
		checkingWords(word, "chi", fileName);
		fInput >> scalar_dissipation_rate;

		fInput >> word;
		checkingWords(word, "nCsi", fileName);
		fInput >> number_of_points_mixture_fraction;
			
		fInput >> word;
		checkingWords(word, "nVariances", fileName);
		fInput >> number_of_points_variance;
		
		fInput >> word;
		checkingWords(word, "nSpecies", fileName);
		fInput >> number_of_species;
	
		// Memory Allocation
		names.resize(number_of_species+1);
		ChangeDimensions(number_of_points_mixture_fraction, &csi);
		ChangeDimensions(number_of_points_variance, &variance);
		ChangeDimensions(number_of_points_mixture_fraction, number_of_points_variance, &density);
		ChangeDimensions(number_of_points_mixture_fraction, number_of_points_variance, &z_r);
		ChangeDimensions(number_of_points_mixture_fraction, number_of_points_variance, &zv2_r);
		ChangeDimensions(number_of_points_mixture_fraction, number_of_points_variance, &temperature);
		w = new BzzMatrix[number_of_points_mixture_fraction+1];
		for(i=1;i<=number_of_points_mixture_fraction;i++)
			ChangeDimensions(number_of_points_variance, number_of_species, &w[i]);
		ChangeDimensions(number_of_points_mixture_fraction, number_of_species, &composition);

		for(i=1;i<=number_of_species;i++)
			fInput >> names[i];
	
		fInput >> word;
		checkingWords(word, "csi", fileName);
		for(i=1;i<=number_of_points_mixture_fraction;i++)
			fInput >> csi[i];

		fInput >> word;
		checkingWords(word, "variance", fileName);
		for(i=1;i<=number_of_points_variance;i++)
			fInput >> variance[i];

		fInput >> word;
		checkingWords(word, "density", fileName);
		for(k=1;k<=number_of_points_mixture_fraction;k++)
			for(j=1;j<=number_of_points_variance;j++)
				fInput >> density[k][j];

		fInput >> word;
		checkingWords(word, "z_reynolds", fileName);
		for(k=1;k<=number_of_points_mixture_fraction;k++)
			for(j=1;j<=number_of_points_variance;j++)
				fInput >> z_r[k][j];

		fInput >> word;
		checkingWords(word, "zv2_reynolds", fileName);
		for(k=1;k<=number_of_points_mixture_fraction;k++)
			for(j=1;j<=number_of_points_variance;j++)
				fInput >> zv2_r[k][j];

		fInput >> word;
		checkingWords(word, "temperature", fileName);
		for(k=1;k<=number_of_points_mixture_fraction;k++)
			for(j=1;j<=number_of_points_variance;j++)
				fInput >> temperature[k][j];
		
		fInput >> word;
		checkingWords(word, "massfractions", fileName);
		for(k=1;k<=number_of_points_mixture_fraction;k++)
		{
			int checking_number;
			cout << "Reading slice no. " << k << "..." <<endl;
			fInput >> checking_number;
			if (checking_number!=k)
				ErrorMessage("Error in importing the flamelet table: slice no. ");
	
			for(j=1;j<=number_of_points_variance;j++)
			{
				double sumchecking = 0.;
				for(i=1;i<=number_of_species;i++)
				{	
					fInput >> w[k][j][i];
					sumchecking += w[k][j][i];
				}
				if (sumchecking <= 0.999999 || sumchecking >= 1.00001)
					ErrorMessage("Error in reading the file: Sum mass fractions...");
			}
		}
		fInput >> word;
		checkingWords(word, "end", fileName);
	
		for(k=1;k<=number_of_points_mixture_fraction;k++)
			for(i=1;i<=number_of_species;i++)
				composition[k][i] = w[k][1][i];
	
		cout << " Flamelet no. " << index << endl;
		cout << " Scalar Dissipation Rate:  " << scalar_dissipation_rate << " Hz" << endl;
		cout << " Number of Species:  " << number_of_species << endl;
		cout << " Number of Slices - Mixture Fraction Space:  " << number_of_points_mixture_fraction << endl;
		cout << " Number of Slices - Normal Variance Space:  " << number_of_points_variance << endl;
	
		cout << " Generating the flamelet... " << endl; 
		flames[index].assign_flamelet(mix, names, pressure_pa, scalar_dissipation_rate, csi, temperature.GetColumn(1), composition);

		cout << " Importing the PDF... " << endl; 		
		flames[index].import_pdf(variance, density, temperature, w, z_r, zv2_r);

		cout << " Succesfully DONE!! " << endl << endl; 
	}
	
	lock();
}

#include "addons/OpenSMOKE_SootManager.h" 
void flamelet_group::soot_post_processing(OpenSMOKE_ReactingGas &mix)
{/*
	for (int index = 1; index<=nFlamelets; index++)
	{
		OpenSMOKE_SootManager soot_manager;
		string *names_of_species;
		names_of_species = new string[flames[index].NC+1];
		for (int j = 1; j<=flames[index].NC; j++)
			names_of_species[j] = uppercase(flames[index].species.names[j]);
		soot_manager.recognizeSpecies(flames[index].NC, names_of_species, flames[index].species.pm);

		BzzVector wGas(flames[index].NC);
		BzzVector xGas(flames[index].NC);
		BzzVector cGas(flames[index].NC);
		double rhoGas;
		double tGas;

		char sootfilename[200];
		char number[2];
		strcpy(sootfilename, "Output/Soot_");
		my_itoa(index, number, 10);
		strcat(sootfilename, number);
		strcat(sootfilename, ".out");
		ofstream fOut(sootfilename, ios::out);
		fOut.setf(ios::scientific);

		for(int k=1;k<=flames[index].Nslices_MixtureFraction;k++)	
		{
			for(int i=1;i<=flames[index].Nslices_Variance;i++)
			{
				flames[index].requests_for_soot( mix, k, i, wGas, xGas, cGas, rhoGas, tGas);
				soot_manager.small_calculateAll(wGas, xGas, cGas, rhoGas);
				soot_manager.large_calculateAll(wGas, xGas, cGas, rhoGas);
				
				fOut <<  flames[index].mixture_fraction_space[k] << "\t";
				fOut <<  flames[index].normal_variance_space[i] << "\t";
				
				fOut <<  soot_manager.large_omega_tot << "\t";
				fOut <<  soot_manager.large_x_tot << "\t";
				fOut <<  soot_manager.large_c_tot << "\t";
				fOut <<  soot_manager.large_fv_tot << "\t";

				fOut <<  soot_manager.small_omega_tot << "\t";
				fOut <<  soot_manager.small_x_tot << "\t";
				fOut <<  soot_manager.small_c_tot << "\t";
				fOut <<  soot_manager.small_fv_tot << "\t";

				fOut << endl;
			}
			fOut << endl;
		}
		fOut.close();

	} // end for cycle for flamelets
	*/
}



void flamelet_group::apply_soot()
{
	int i, j, k;
	int nSources = 8;
	std::string *file_names;
	ofstream *file_list;

	file_list	= new ofstream[nSources+1];
	file_names	= new std::string[nSources+1];

	file_names[1] = "LargeBins_fv";
	file_names[2] = "SmallBins_fv";
	file_names[3] = "LargeBins_omega";
	file_names[4] = "SmallBins_omega";
	file_names[5] = "LargeBins_x";
	file_names[6] = "SmallBins_x";
	file_names[7] = "LargeBins_c";
	file_names[8] = "SmallBins_c";


	for (i=1;i<=nSources;i++)
	{
		char file_name[60];
		strcpy(file_name, "Output/");
		strcat(file_name, file_names[i].c_str());
		strcat(file_name, ".soot");
		file_list[i].open(file_name, ios::out);
		file_list[i].setf(ios::scientific);

		// Intestazione del file
		file_list[i] << "nFlamelets " << nFlamelets << endl << endl;
	}


	// Costruzione della libreria per i termini sorgenti
	for (int index=1;index<=nFlamelets;index++)
	{	
		BzzVector mixture_fraction_space(flames[index].Nslices_MixtureFraction);
		BzzVector normal_variance_space(flames[index].Nslices_Variance);
		BzzMatrix large_omega_tot(flames[index].Nslices_MixtureFraction, flames[index].Nslices_Variance);
		BzzMatrix large_x_tot(flames[index].Nslices_MixtureFraction, flames[index].Nslices_Variance);
		BzzMatrix large_c_tot(flames[index].Nslices_MixtureFraction, flames[index].Nslices_Variance);
		BzzMatrix large_fv_tot(flames[index].Nslices_MixtureFraction, flames[index].Nslices_Variance);
		BzzMatrix small_omega_tot(flames[index].Nslices_MixtureFraction, flames[index].Nslices_Variance);
		BzzMatrix small_x_tot(flames[index].Nslices_MixtureFraction, flames[index].Nslices_Variance);
		BzzMatrix small_c_tot(flames[index].Nslices_MixtureFraction, flames[index].Nslices_Variance);
		BzzMatrix small_fv_tot(flames[index].Nslices_MixtureFraction, flames[index].Nslices_Variance);

		char sootfilename[200];
		char number[2];
		strcpy(sootfilename, "Output/Soot_");
		my_itoa(index, number, 10);
		strcat(sootfilename, number);
		strcat(sootfilename, ".out");
		ifstream fInput(sootfilename, ios::in);

		// Reading from Soot file
		for(k=1;k<=flames[index].Nslices_MixtureFraction;k++)	
		{
			for(int i=1;i<=flames[index].Nslices_Variance;i++)
			{
				fInput >>  mixture_fraction_space[k];
				fInput >>  normal_variance_space[i];
				
				fInput >>  large_omega_tot[k][i];
				fInput >>  large_x_tot[k][i];
				fInput >>  large_c_tot[k][i];
				fInput >>  large_fv_tot[k][i];

				fInput >>  small_omega_tot[k][i];
				fInput >>  small_x_tot[k][i];
				fInput >>  small_c_tot[k][i];
				fInput >>  small_fv_tot[k][i];
			}
		}

		// Writing on Soot file LargeBins And Small Bins
		for (i=1;i<=nSources;i++)
		{
			file_list[i] << "chi "				<< flames[index].chi << endl;
			file_list[i] << "nCsi "				<< flames[index].Nslices_MixtureFraction << endl;
			file_list[i] << "nVariance "		<< flames[index].Nslices_Variance << endl;
		
			// Mixture Fraction
			for(k=1;k<=flames[index].Nslices_MixtureFraction;k++)	
				file_list[i]  << mixture_fraction_space[k] << endl;
			file_list[i]  << endl;

			// Variance
			for(int j=1;j<=flames[index].Nslices_Variance;j++)
				file_list[i] << normal_variance_space[j] << endl;
			file_list[i] << endl;
		}

		// Source Fields
		for(k=1;k<=flames[index].Nslices_MixtureFraction;k++)	
		{
			for(i=1;i<=flames[index].Nslices_Variance;i++)
			{
				file_list[1] << large_fv_tot[k][i]		<< "\t";
				file_list[1] << large_fv_tot[k][i]		<< "\t";
				file_list[2] << small_fv_tot[k][i]		<< "\t";
				file_list[2] << small_fv_tot[k][i]		<< "\t";
				file_list[3] << large_omega_tot[k][i]	<< "\t";
				file_list[3] << large_omega_tot[k][i]	<< "\t";
				file_list[4] << small_omega_tot[k][i]	<< "\t";
				file_list[4] << small_omega_tot[k][i]	<< "\t";
				file_list[5] << large_x_tot[k][i]		<< "\t";
				file_list[5] << large_x_tot[k][i]		<< "\t";
				file_list[6] << small_x_tot[k][i]		<< "\t";
				file_list[6] << small_x_tot[k][i]		<< "\t";
				file_list[7] << large_c_tot[k][i]		<< "\t";
				file_list[7] << large_c_tot[k][i]		<< "\t";
				file_list[8] << small_c_tot[k][i]		<< "\t";
				file_list[8] << small_c_tot[k][i]		<< "\t";

				for (j=1;j<=nSources;j++)
					file_list[j]<< endl;
			}
			for (j=1;j<=nSources;j++)
				file_list[j]<< endl;
		}
		for (j=1;j<=nSources;j++)
			file_list[j]<< endl;
	}

	// Chiusura dei files
	for (i=1;i<=nSources;i++)
		file_list[i].close();
}


void flamelet_group::apply_source_betaPDF_SootProperties(OpenSMOKE_ReactingGas &mix, int nSources, BzzVector *m0, BzzVector *fv)
{
	int i,j;
	std::string *file_names;
	ofstream *file_list;

	file_list	= new ofstream[nSources+1];
	file_names	= new std::string[nSources+1];

	file_names[1] = "m0";
	file_names[2] = "fv";

	for (i=1;i<=nSources;i++)
	{
		char file_name[60];
		strcpy(file_name, "Output/");
		strcat(file_name, file_names[i].c_str());
		strcat(file_name, ".source");
		file_list[i].open(file_name, ios::out);
		file_list[i].setf(ios::scientific);

		// Intestazione del file
		file_list[i] << "nFlamelets " << nFlamelets << endl << endl;
	}

	// Costruzione della libreria per i termini sorgenti
	for (j=1;j<=nFlamelets;j++)
	{			
		// Construction of source term library
		cout << "Flamelet no. " << j << " - " << "Scalar Dissipation Rate (st.): " << flames[j].chi << " Hz" << endl;
		cout << "------------------------------------------------------------------------" << endl;
		flames[j].apply_source_betaPDF_SootProperties(j, nSources, file_list, file_names, m0[j], fv[j]);
	}

	// Chiusura dei files
	for (i=1;i<=nSources;i++)
		file_list[i].close();
}


void flamelet_group::print_on_file_FLUENT(string fileName, OpenSMOKE_ReactingGas &mix)
{
	int i,j,k;
	const double ZERO = 0.;
	double mw;

	ofstream fFLUENT;
	openOutputFileAndControl(fFLUENT, fileName);
	fFLUENT.setf(ios::scientific);

	int		NC = flames[1].NC;
	double	TF = flames[1].temperature[1];
	double	TO = flames[1].temperature[flames[1].NP];
	BzzVector mole_fractions(NC);
	BzzVector mass_fractions(NC);

	fFLUENT << "(0 \"Fluent6.3\")" << endl;
	fFLUENT << endl;

	fFLUENT << "(0 \"PDF Variables:\")" << endl;
	fFLUENT << "(100 (101)(" << endl;
	fFLUENT << "(pdf/includeequil? #t )" << endl;
	fFLUENT << "(pdf/frozenboundary? #t )" << endl;
	fFLUENT << "(pdf/multi? #f )" << endl;
	fFLUENT << "(pdf/adiabatic? #t )" << endl;
	fFLUENT << "(pdf/chemistry 3)" << endl;
	fFLUENT << "(pdf/empirialfuel? #f )" << endl;
	fFLUENT << "(pdf/empirialscnd? #f )" << endl;
	fFLUENT << "(pdf/parpremix #f )" << endl;
	fFLUENT << "(pdf/numspecies " << NC << ")" << endl;
	fFLUENT << "(pdf/initpressure 101325.000000)" << endl;
	fFLUENT << "(pdf/mintemp 298.000000)" << endl;
	fFLUENT << "(pdf/maxtempratio 1.000000)" << endl;
	fFLUENT << "(pdf/tloss 0.667000)" << endl;
	fFLUENT << "(pdf/tadd 0.250000)" << endl;
	fFLUENT << "(pdf/max-species "<< NC+2 <<")" << endl;
	fFLUENT << "(pdf/num-include-species 0)" << endl;
	fFLUENT << "(pdf/num-exclude-species "<< NC+1 <<")" << endl;
	fFLUENT << "(pdf/mole-or-mass? #t )" << endl;
//	fFLUENT << "(pdf/exclude-species NO NO2 N2O H2O(L) C(S) N NH NH2 NH3 NNH HCN HNO CN H2CN HCNN HCNO HOCN HNCO )" << endl;
	fFLUENT << "(pdf/exclude-species )" << endl;
	fFLUENT << "(pdf/include-species )" << endl;
	fFLUENT << "(pdf/stream1temp " << TF <<")" << endl;
	fFLUENT << "(pdf/stream2temp " << TO <<")" << endl;
	fFLUENT << "(pdf/stream3temp 300.000000)" << endl;
	fFLUENT << "(pdf/numf1point " << flames[1].Nslices_MixtureFraction << ")" << endl;
	fFLUENT << "(pdf/numf2point 21)" << endl;
	fFLUENT << "(pdf/richlimit1 1.000000)" << endl;
	fFLUENT << "(pdf/richlimit2 1.000000)" << endl;
	fFLUENT << "(pdf/numfvpoint " << flames[1].Nslices_Variance << ")" << endl;
	fFLUENT << "(pdf/numstream 2)" << endl;
	fFLUENT << "(pdf/betapdf? 1)" << endl;
	fFLUENT << "(pdf/empshfuel 50000000.000000)" << endl;
	fFLUENT << "(pdf/empcpfuel 1000.000000)" << endl;
	fFLUENT << "(pdf/empshscnd 50000000.000000)" << endl;
	fFLUENT << "(pdf/empcpscnd 1000.000000)" << endl;
	fFLUENT << "(pdf/numenthalpy 41)" << endl;
	fFLUENT << "(pdf/nne/cold 30)" << endl;
	fFLUENT << "(pdf/nne/hot 11)" << endl;
	fFLUENT << "(pdf/numfl " << nFlamelets << ")" << endl;
	fFLUENT << "(pdf/expinv 1.000000)" << endl;
	fFLUENT << "(pdf/scadmx 31.000000)" << endl;
	fFLUENT << "(pdf/flamelet-int-params " << nFlamelets - 1 << " " << flames[1].NP << ")" << endl;
	fFLUENT << "(pdf/flamelet-real-params 1.000000e-002 5.000000e+000 1.000000e+000 2.000000e+000 1.000000e-010 1.000000e-014 1.000000e-005 1.000000e+003)" << endl;
	fFLUENT << "(pdf/flamelet-mechanism-file-name C:\\Research\\My Projects\\KineticSchemes\\schemiFluent\\reduced25.che )" << endl;
	fFLUENT << "(pdf/flamelet-thermo-file-name C:\\Fluent.Inc\\fluent6.3.26\\cpropep\\data\\thermo.db )" << endl;
	fFLUENT << "))" << endl;
	fFLUENT << endl;

	// ------------------------------------------------------------------------------------------------------
	// Name of species
	// ------------------------------------------------------------------------------------------------------
	fFLUENT << "(0 \"Species Input Data:\")" << endl;
	fFLUENT << "(100 (102 " << NC <<")(" << endl;
	for(i=1;i<=NC;i++)
		fFLUENT << mix.names[i] << endl;
	fFLUENT << "))" << endl;
	fFLUENT << endl;

	// ------------------------------------------------------------------------------------------------------
	// Fuel Composition - Mole fractions
	// ------------------------------------------------------------------------------------------------------
	mass_fractions = flames[1].massfractions.GetRow(1);
	mix.GetMWAndMoleFractionsFromMassFractions(mw, mole_fractions, mass_fractions);
	fFLUENT << "(100 (106 " << NC <<")(" << endl;
	for(j=1;j<=NC;j++)
		fFLUENT << mole_fractions[j] << endl;
	fFLUENT << "))" << endl;
	fFLUENT << endl;

	// ------------------------------------------------------------------------------------------------------
	// Air Composition - Mole Fractions
	// ------------------------------------------------------------------------------------------------------
	mass_fractions = flames[1].massfractions.GetRow(flames[1].NP);
	mix.GetMWAndMoleFractionsFromMassFractions(mw, mole_fractions, mass_fractions);
	fFLUENT << "(100 (107 " << NC <<")(" << endl;
	for(j=1;j<=NC;j++)
		fFLUENT << mole_fractions[j] << endl;
	fFLUENT << "))" << endl;
	fFLUENT << endl;

	// ------------------------------------------------------------------------------------------------------
	// Air Composition - Secondary Stream
	// ------------------------------------------------------------------------------------------------------
	mass_fractions = flames[1].massfractions.GetRow(flames[1].NP);
	mix.GetMWAndMoleFractionsFromMassFractions(mw, mole_fractions, mass_fractions);
	fFLUENT << "(100 (108 " << NC <<")(" << endl;
	for(j=1;j<=NC;j++)
		fFLUENT << mole_fractions[j] << endl;
	fFLUENT << "))" << endl;
	fFLUENT << endl;

	// ------------------------------------------------------------------------------------------------------
	// TODO
	// ------------------------------------------------------------------------------------------------------
	fFLUENT << "(100 (109 " << flames[1].Nslices_MixtureFraction+1 <<")(" << endl;
	for(j=1;j<=flames[1].Nslices_MixtureFraction+1;j++)
		fFLUENT << ZERO << endl;
	fFLUENT << "))" << endl;
	fFLUENT << endl;

	// ------------------------------------------------------------------------------------------------------
	// MIXTURE FRACTION
	// ------------------------------------------------------------------------------------------------------
	fFLUENT << "(200 (201 " << flames[1].Nslices_MixtureFraction <<")(" << endl;
	for(j=1;j<=flames[1].Nslices_MixtureFraction;j++)
		fFLUENT << flames[1].mixture_fraction_space[j] << endl;
	fFLUENT << "))" << endl;
	fFLUENT << endl;

	// ------------------------------------------------------------------------------------------------------
	// MIXTURE FRACTION NORMALIZED VARIANCE
	// ------------------------------------------------------------------------------------------------------
	fFLUENT << "(200 (202 " << flames[1].Nslices_Variance <<")(" << endl;
	for(j=1;j<=flames[1].Nslices_Variance;j++)
		fFLUENT << 0.25*flames[1].normal_variance_space[j]/flames[1].normal_variance_space[flames[1].Nslices_Variance] << endl;
	fFLUENT << "))" << endl;
	fFLUENT << endl;

	// ------------------------------------------------------------------------------------------------------
	// SCALAR DISSIPATION RATES NORMALIZED
	// ------------------------------------------------------------------------------------------------------
	fFLUENT << "(200 (262 " << nFlamelets <<")(" << endl;
	for(j=1;j<=nFlamelets;j++)
		fFLUENT << flames[j].chi / flames[nFlamelets].chi << endl;
	fFLUENT << "))" << endl;
	fFLUENT << endl;

	// ------------------------------------------------------------------------------------------------------
	// TEMPERATURE
	// ------------------------------------------------------------------------------------------------------
	fFLUENT << "(200 (214 " << flames[1].Nslices_MixtureFraction << " " 
							<< flames[1].Nslices_Variance        << " "
							<< nFlamelets                        <<")(" << endl;
	for(k=1;k<=nFlamelets;k++)
		flames[k].print_on_file_appending_FLUENT(fFLUENT, -1, mix);
	fFLUENT << "))" << endl;
	fFLUENT << endl;

	// ------------------------------------------------------------------------------------------------------
	// DENSITY
	// ------------------------------------------------------------------------------------------------------
	fFLUENT << "(200 (216 " << flames[1].Nslices_MixtureFraction << " " 
							<< flames[1].Nslices_Variance        << " "
							<< nFlamelets                        <<")(" << endl;
	for(k=1;k<=nFlamelets;k++)
		flames[k].print_on_file_appending_FLUENT(fFLUENT, 0, mix);
	fFLUENT << "))" << endl;
	fFLUENT << endl;

	// ------------------------------------------------------------------------------------------------------
	// SPECIES MOLE FRACTIONS
	// ------------------------------------------------------------------------------------------------------
	fFLUENT << "(200 (218 " << flames[1].Nslices_MixtureFraction << " " 
							<< flames[1].Nslices_Variance        << " "
							<< nFlamelets				         << " "
							<< NC                                <<")(" << endl;
	for(j=1;j<=NC;j++)
		for(k=1;k<=nFlamelets;k++)
			flames[k].print_on_file_appending_FLUENT(fFLUENT, j, mix);
	fFLUENT << "))" << endl;
	fFLUENT << endl;

	// ------------------------------------------------------------------------------------------------------
	// INTERMEDIATE DATA
	// ------------------------------------------------------------------------------------------------------


	fFLUENT.close();
}

void flamelet_group::print_on_file_look_up_table(string folderName, OpenSMOKE_ReactingGas &mix)
{
	string command = "mkdir " + folderName;
	system(command.c_str());

	ofstream fLookUpTable;
	openOutputFileAndControl(fLookUpTable, folderName + "/LookUpTable.out");
	fLookUpTable.setf(ios::scientific);

	fLookUpTable << "Adiabatic"	<< endl;
	fLookUpTable << "nc "		<< mix.NumberOfSpecies() << endl;
	fLookUpTable << "chi "		<< nFlamelets << endl;
	for(int j=1;j<=nFlamelets;j++)
		fLookUpTable << flames[j].chi << endl;
	fLookUpTable.close();

	int k;
	for(k=1;k<=nFlamelets;k++)
	{
		cout << "Flamelet no. " << k << endl;

		int j;

		stringstream kstring;
		kstring << k;

		{
			ofstream fSingleFlamelet;
			openOutputFileAndControl(fSingleFlamelet, folderName + "/SR_" + kstring.str() + ".out");
			fSingleFlamelet.setf(ios::scientific);

			// 00. Mole fractions
			flames[k].update_mole_fractions(mix);

			fSingleFlamelet << "version " << "0.1" << endl;
			fSingleFlamelet << "nc		" << mix.NumberOfSpecies() <<  endl ;

			// 01. Mixture fraction
			fSingleFlamelet << "mf  " << flames[k].Nslices_MixtureFraction	<< endl;
			for(j=1;j<=flames[k].Nslices_MixtureFraction;j++)
				fSingleFlamelet << flames[k].mixture_fraction_space[j] << endl;

			// 02. Mixture fraction variance
			fSingleFlamelet << "mfv " << flames[k].Nslices_Variance			<< endl;
			for(j=1;j<=flames[k].Nslices_Variance;j++)
				fSingleFlamelet << flames[k].normal_variance_space[j] << endl;

			// 03. Mixture fraction
			fSingleFlamelet << "mf_reynolds"	<< endl;
			flames[k].print_on_file_appending_FLUENT(fSingleFlamelet, PDF_WRITE_MF_REYNOLDS, mix);		

			// 04. Mixture fraction variance
			fSingleFlamelet << "mfv_reynolds"	<< endl;
			flames[k].print_on_file_appending_FLUENT(fSingleFlamelet, PDF_WRITE_MFV_REYNOLDS, mix);		

			// 05. Temperature
			fSingleFlamelet << "temperature" << endl;
			flames[k].print_on_file_appending_FLUENT(fSingleFlamelet, PDF_WRITE_TEMPERATURE, mix);		
			
			// 06. Density
			fSingleFlamelet << "density"	<< endl;
			flames[k].print_on_file_appending_FLUENT(fSingleFlamelet, PDF_WRITE_DENSITY, mix);		

			// 07. Enthalpy
			fSingleFlamelet << "enthalpy"	<< endl;
			flames[k].print_on_file_appending_FLUENT(fSingleFlamelet, PDF_WRITE_ENTHALPY, mix);		

			// 08. Molecular weight
			fSingleFlamelet << "mw"	<< endl;
			flames[k].print_on_file_appending_FLUENT(fSingleFlamelet, PDF_WRITE_MOLECULAR_WEIGHT, mix);		

			// 09. Specific heat
			fSingleFlamelet << "cp"	<< endl;
			flames[k].print_on_file_appending_FLUENT(fSingleFlamelet, PDF_WRITE_SPECIFIC_HEAT, mix);		

			// 10. Thermal conductivity
			fSingleFlamelet << "lambda"	<< endl;
			flames[k].print_on_file_appending_FLUENT(fSingleFlamelet, PDF_WRITE_THERMAL_CONDUCTIVITY, mix);		

			// 11. Dynamic viscosity
			fSingleFlamelet << "mu"	<< endl;
			flames[k].print_on_file_appending_FLUENT(fSingleFlamelet, PDF_WRITE_DYNAMIC_VISCOSITY, mix);		

			// 12. Absorption coefficient
			fSingleFlamelet << "as"	<< endl;
			flames[k].print_on_file_appending_FLUENT(fSingleFlamelet, PDF_WRITE_ABSORPTION_COEFFICIENT, mix);	

			// 13. Chemical species
			for(j=1;j<=flames[k].NC;j++)
			{
				fSingleFlamelet << mix.names[j] << " " << mix.M[j] << endl;
				flames[k].print_on_file_appending_FLUENT(fSingleFlamelet, j, mix);	// species (mole fractions)
			}

			fSingleFlamelet.close();

		}

		{
			string dummy;
			const int SIZE = 40;
			char name[SIZE];

			BzzSave fSingleFlamelet('*', folderName + "/SR_" + kstring.str() + ".bin");

			// 00. Mole fractions
			flames[k].update_mole_fractions(mix);

			dummy = "version";
			strcpy(name, dummy.c_str());
			fSingleFlamelet.fileSave.write((char*) name, SIZE);

			dummy = "0.1";
			strcpy(name, dummy.c_str());
			fSingleFlamelet.fileSave.write((char*) name, SIZE);

			dummy = "nc";
			strcpy(name, dummy.c_str());
			fSingleFlamelet.fileSave.write((char*) name, SIZE);
			fSingleFlamelet << mix.NumberOfSpecies();

			// 01. Mixture fraction
			dummy = "mf";
			strcpy(name, dummy.c_str());
			fSingleFlamelet.fileSave.write((char*) name, SIZE);
			fSingleFlamelet << flames[k].Nslices_MixtureFraction;
			for(j=1;j<=flames[k].Nslices_MixtureFraction;j++)
				fSingleFlamelet << flames[k].mixture_fraction_space[j];

			// 02. Mixture fraction variance
			dummy = "mfv";
			strcpy(name, dummy.c_str());
			fSingleFlamelet.fileSave.write((char*) name, SIZE);
			fSingleFlamelet << flames[k].Nslices_Variance;
			for(j=1;j<=flames[k].Nslices_Variance;j++)
				fSingleFlamelet << flames[k].normal_variance_space[j];

			// 03. Mixture fraction
			dummy = "mf_reynolds";
			strcpy(name, dummy.c_str());
			fSingleFlamelet.fileSave.write((char*) name, SIZE);
			flames[k].print_on_file_appending_FLUENT(fSingleFlamelet, PDF_WRITE_MF_REYNOLDS, mix);		

			// 04. Mixture fraction variance
			dummy = "mfv_reynolds";
			strcpy(name, dummy.c_str());
			fSingleFlamelet.fileSave.write((char*) name, SIZE);
			flames[k].print_on_file_appending_FLUENT(fSingleFlamelet, PDF_WRITE_MFV_REYNOLDS, mix);		

			// 05. Temperature
			dummy = "temperature";
			strcpy(name, dummy.c_str());
			fSingleFlamelet.fileSave.write((char*) name, SIZE);
			flames[k].print_on_file_appending_FLUENT(fSingleFlamelet, PDF_WRITE_TEMPERATURE, mix);		
			
			// 06. Density
			dummy = "density";
			strcpy(name, dummy.c_str());
			fSingleFlamelet.fileSave.write((char*) name, SIZE);
			flames[k].print_on_file_appending_FLUENT(fSingleFlamelet, PDF_WRITE_DENSITY, mix);		

			// 07. Enthalpy
			dummy = "enthalpy";
			strcpy(name, dummy.c_str());
			fSingleFlamelet.fileSave.write((char*) name, SIZE);
			flames[k].print_on_file_appending_FLUENT(fSingleFlamelet, PDF_WRITE_ENTHALPY, mix);		

			// 08. Molecular weight
			dummy = "mw";
			strcpy(name, dummy.c_str());
			fSingleFlamelet.fileSave.write((char*) name, SIZE);
			flames[k].print_on_file_appending_FLUENT(fSingleFlamelet, PDF_WRITE_MOLECULAR_WEIGHT, mix);		

			// 09. Specific heat
			dummy = "cp";
			strcpy(name, dummy.c_str());
			fSingleFlamelet.fileSave.write((char*) name, SIZE);
			flames[k].print_on_file_appending_FLUENT(fSingleFlamelet, PDF_WRITE_SPECIFIC_HEAT, mix);		

			// 10. Thermal conductivity
			dummy = "lambda";
			strcpy(name, dummy.c_str());
			fSingleFlamelet.fileSave.write((char*) name, SIZE);
			flames[k].print_on_file_appending_FLUENT(fSingleFlamelet, PDF_WRITE_THERMAL_CONDUCTIVITY, mix);		

			// 11. Dynamic viscosity
			dummy = "mu";
			strcpy(name, dummy.c_str());
			fSingleFlamelet.fileSave.write((char*) name, SIZE);
			flames[k].print_on_file_appending_FLUENT(fSingleFlamelet, PDF_WRITE_DYNAMIC_VISCOSITY, mix);		

			// 12. Absorption coefficient
			dummy = "as";
			strcpy(name, dummy.c_str());
			fSingleFlamelet.fileSave.write((char*) name, SIZE);
			flames[k].print_on_file_appending_FLUENT(fSingleFlamelet, PDF_WRITE_ABSORPTION_COEFFICIENT, mix);	

			// 13. Chemical species
			for(j=1;j<=flames[k].NC;j++)
			{
				strcpy(name, mix.names[j].c_str());
				fSingleFlamelet.fileSave.write((char*) name, SIZE);
				fSingleFlamelet << mix.M[j];

				flames[k].print_on_file_appending_FLUENT(fSingleFlamelet, j, mix);	// species (mole fractions)
			}

			fSingleFlamelet.End();
		}
	}

	for(k=1;k<=nFlamelets;k++)
	{
		stringstream kstring;
		kstring << k;

		ofstream fSingleFlamelet;
		openOutputFileAndControl(fSingleFlamelet, folderName + "/SR_" + kstring.str() + ".gnuplot");
		fSingleFlamelet.setf(ios::scientific);

		fSingleFlamelet << "Z_Fav(1)          ";
		fSingleFlamelet << "ZV2N_Fav(2)       ";
		fSingleFlamelet << "Z_Rey(3)          ";
		fSingleFlamelet << "ZV2N_Rey(4)       ";
		fSingleFlamelet << "rho_Rey[kg/m3](5) ";
		fSingleFlamelet << "T_Fav[K](6)       ";
		int count = 7;
		for(int w=1;w<=flames[k].NC;w++)
			fSingleFlamelet << "W_" << mix.names[w] << "_Fav(" << count++ <<")\t";
		fSingleFlamelet << endl;

		for(int i=1;i<=flames[k].Nslices_Variance;i++)
		{
			for(int j=1;j<=flames[k].Nslices_MixtureFraction;j++)
			{
				fSingleFlamelet << flames[k].mixture_fraction_space[j] << "\t";
				fSingleFlamelet << flames[k].normal_variance_space[i]  << "\t";
				fSingleFlamelet << flames[k].z_r_library[j][i] << "\t";
				fSingleFlamelet << flames[k].zv2_r_library[j][i]  << "\t";
				fSingleFlamelet << flames[k].density_r_library[j][i]   << "\t";			
				fSingleFlamelet << flames[k].T_f_library[j][i]		   << "\t";
				for(int w=1;w<=flames[k].NC;w++)
					fSingleFlamelet << flames[k].w_f_library[j][i][w]	   << "\t";
				fSingleFlamelet << endl;
			}
			fSingleFlamelet << endl;
		}
		fSingleFlamelet.close();
	}
}
