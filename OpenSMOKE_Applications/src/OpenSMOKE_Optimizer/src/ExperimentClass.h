/***************************************************************************
 *   Copyright (C) 2006-2008 by Alberto Cuoci							   *
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

#ifndef EXPERIMENTCLASS_H
#define EXPERIMENTCLASS_H

#include "OpenSMOKE.hpp"

class ExperimentClass
{
    public:

		ExperimentClass();
		void SetName(const string name);

		void   Setup(const string _folderName, const int _index);
        BzzVector GiveMeObjectiveFunction(ofstream &fLog, string *names,   BzzVector &numerical_time, BzzVector &numerical_space,
                                                        BzzVector &numerical_T, BzzVector &numerical_V, BzzMatrix &numerical_mass,
                                                        BzzMatrix &numerical_mole);

		BzzVector GiveMeObjectiveFunction( ofstream &fLog, string *names, double numerical_T, double numerical_V, 
										   BzzVector &numerical_mass, BzzVector &numerical_mole);

		BzzVector GiveMeObjectiveFunction( ofstream &fLog, BzzVector &numerical_phi, BzzVector &numerical_v);


        void PostProcessing(string flag, string *names, BzzVector &numerical_time,   BzzVector &numerical_space,
                                BzzVector &numerical_T, BzzVector &numerical_V,
                                BzzMatrix &numerical_mass, BzzMatrix &numerical_mole);
		void PostProcessing(string flag, BzzVector &numerical_phi,   BzzVector &numerical_T, BzzVector &numerical_v);

        void Gnuplot_plots();
        void Gnuplot_plots_online();
        void Latex_figure(OpenSMOKE_LatexInterface &latex);
		BzzVector ExtractTemperatures();
		BzzVector ExtractSupportTemperatures();

        string   kindOfReactor;
        string   folderName;

		void GiveMeRegressionFunction( string *names, 
									   BzzVector &numerical_time, BzzVector &numerical_space,
									   BzzVector &numerical_T, BzzVector &numerical_V, BzzMatrix &numerical_mole,			
									   BzzVector &required_x, BzzVector &expected_y);
		
		BzzVector Tau_History;
		BzzVector Eta_History;
		BzzVector T_History;
		BzzMatrix mole_History;

		int      nData;
		BzzMatrix *experimental_data;
    
	private:

        int     indexExperiment;
        int     indipendent_variable;

        
        string  *labels;
        string  *labels_print;
        string  *kind;
        string  *f;
        BzzVector   weights;
        BzzVectorInt      nSum;
        OpenSMOKE_Matrix<string> sum_labels;

       

		string name_object;
        void ErrorMessage(const string message);
        void WarningMessage(const string message);
        
		int  recognize_species(string *names, string label, int Nspecies);

};

#endif // EXPERIMENTCLASS_H
