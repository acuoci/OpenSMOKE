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

#if !defined(FLAMELET)
#define FLAMELET

#include "engine/OpenSMOKE_ReactingGas.h"
#include "sources_in_flamelet_library.h"
#include "Utilities.h"

class OpenSMOKE_ReactingGas;

enum  pdf_write_type {	PDF_WRITE_TEMPERATURE, PDF_WRITE_DENSITY, PDF_WRITE_MF_REYNOLDS, PDF_WRITE_MFV_REYNOLDS,
						PDF_WRITE_SPECIFIC_HEAT, PDF_WRITE_THERMAL_CONDUCTIVITY, PDF_WRITE_DYNAMIC_VISCOSITY,
						PDF_WRITE_ABSORPTION_COEFFICIENT, PDF_WRITE_ENTHALPY, PDF_WRITE_MOLECULAR_WEIGHT };
class flamelet
{
	friend class flamelet_group;
protected:

	BzzVector *z_r_library;
	BzzVector *zv2_r_library;
	BzzVector *density_r_library;
	BzzVector *T_f_library;
	BzzVector *mw_f_library;
	BzzMatrix *w_f_library;
	BzzMatrix *x_f_library;

private:
	
	BzzVector density;
	BzzVector udensity;
	BzzVector z_scaled;
	BzzVector zv2_scaled;


	BzzMatrix variance_space;
	BzzMatrix scaled_variance_space;

	void clean_massfractions();

	// BETA Probability Distribution Function
	void apply_betaPDF(int k, double mean);	
	void calculate_betaPDF_and_mean_values( double mean, double variance, 
											double &density_mean, double &temperature_mean, 
											BzzVector &massfractions_mean,
											double &z_r_mean, double &zv2_r_mean);


	void apply_source_betaPDF(int k, double mean, BzzVector &bzz_source, BzzMatrix &bzz_source_mean);
	void calculate_source_betaPDF(	double mean, double variance,
									double &source_mean, 
									BzzVector &bzz_source_mean);

	// CLIPPED GAUSSIAN Probability Distribution Function
	void apply_clippedGaussian(int k, double mean);	
	void calculate_clippedGaussian_and_mean_values( double mean, double variance, 
													double &density_mean, double &temperature_mean, 
													BzzVector &massfractions_mean,
													double &z_r_mean, double &zv2_r_mean);

	double maxNormalVariance;
	BzzVector auxiliary_NP;

public:

	BzzVector temperature;
	BzzMatrix massfractions;
	BzzVector csi;

	void apply_source_betaPDF(	OpenSMOKE_ReactingGas &mix, int index_in_flamelet_group, int nSources, ofstream *file_list, string *file_names);
	void apply_source_betaPDF(	OpenSMOKE_ReactingGas &mix, bool iPerfectlyCorrelatedApproach, int index_in_flamelet_group, int nSources, std::ofstream *file_list, std::string *file_names, BzzVector &m0_N, BzzVector &fv_N,
								nucleation_models nucleation_model, growth_models growth_model, aggregation_models aggregation_model, oxidation_models oxidation_model);

	void write_progress_variables(OpenSMOKE_ReactingGas &mix, const int nProgressVariable, vector<string>& progress_variable_name, BzzVectorInt* progress_variable_index, BzzVector* progress_variable_value);

	int NC;
	int NP;
	double chi;
	double pressure_Pa;

	species_list species;

	void assign_flamelet(	OpenSMOKE_ReactingGas &mix, vector<string> _names, const double _pressure_Pa, const double _chi, BzzVector _csi, 
							BzzVector _temperature, BzzMatrix _massfractions);

	// Library Preparation
	void prepare_library( int _Nslices_Variance, double alfa );
	void import_pdf(BzzVector &_normal_variance,  BzzMatrix &_density, 
					BzzMatrix &_temperature, BzzMatrix *_w, BzzMatrix &_z_r, BzzMatrix &_zv2_r);

	// BETA Probability Distribution Function
	void apply_betaPDF(OpenSMOKE_ReactingGas &mix);

	// CLIPPED GAUSSIAN Probability Distribution Function
	void apply_clippedGaussian();
	
	// Print Functions
	void print_on_file_appending(std::ofstream &fOutput);
	void print_library(char* fileName);
	void print_library(ofstream &fOut);

	sources_in_flamelet_library *sources;
	BzzVector mixture_fraction_space;
	BzzVector normal_variance_space;

	int Nslices_MixtureFraction;
	int Nslices_Variance;
	void requests_for_soot( OpenSMOKE_ReactingGas &mix, int csiRequest, int varianceRequest,
							BzzVector &wRequest, BzzVector &xRequest, 
							BzzVector &cRequest, double &rhoRequest, double &temperatureRequest);

	void apply_source_betaPDF_SootProperties(int index_in_flamelet_group, int nSources, std::ofstream *file_list, std::string *file_names,
											BzzVector &m0, BzzVector &fv);

	void print_on_file_appending_FLUENT(ofstream &fOutput, const int index, OpenSMOKE_ReactingGas &mix);
	void print_on_file_appending_FLUENT(ofstream &fOutput, const pdf_write_type index, OpenSMOKE_ReactingGas &mix);
	void print_on_file_appending_FLUENT(BzzSave &fOutput, const int index, OpenSMOKE_ReactingGas &mix);
	void print_on_file_appending_FLUENT(BzzSave &fOutput, const pdf_write_type index, OpenSMOKE_ReactingGas &mix);
	void update_mole_fractions(OpenSMOKE_ReactingGas &mix);



private:

	void ErrorMessage(const string message);
	void WarningMessage(const string message);
	string name_object;
};

#endif // !defined(FLAMELET)
