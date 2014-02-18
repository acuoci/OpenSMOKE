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

#include "ExperimentClass.h"
#include "version.h"
#include <sstream>

void ExperimentClass::ErrorMessage(const string message)
{
    cout << endl;
	cout << "Class:  CRECK_Optimizer::ExperimentClass"	<< endl;
    cout << "Object: " << name_object					<< endl;
    cout << "Error:  " << message						<< endl;
    cout << "Press a key to continue... "				<< endl;
    getchar();
    exit(-1);
}

void ExperimentClass::WarningMessage(const string message)
{
    cout << endl;
	cout << "Class:  CRECK_Optimizer::ExperimentClass"	<< endl;
    cout << "Object: "		<< name_object				<< endl;
    cout << "Warning:  "	<< message					<< endl;
    cout << "Press a key to continue... "				<< endl;
    getchar();
}

ExperimentClass::ExperimentClass()
{
	name_object			= "[Name not assigned]";
}

void ExperimentClass::SetName(const string name)
{
	name_object = name;
}

void ExperimentClass::Setup(const string _folderName, const int _index)
{
    int i, k;
    ifstream fInput;
    string fileName;
    string x, y;
    string dummy;
    string x_unit;

    folderName  = _folderName;
    fileName    = folderName + "/Exp.inp";

	// Set object name
	stringstream number;
	number << _index;
	SetName("Experiment #"  + number.str());

    openInputFileAndControl(fInput, fileName.c_str());

    fInput >> dummy;
    if (dummy != "REACTOR")
        ErrorMessage("Expected REACTOR key word!");
    fInput >> kindOfReactor;
    if (kindOfReactor != "PFR" && kindOfReactor != "CFDF" && kindOfReactor != "CFDF-TWIN" && kindOfReactor != "FLAMELET" && kindOfReactor != "CSTR" && kindOfReactor != "FLAMESPEEDCURVE")
        ErrorMessage("Only PFR || CFDF || CFDF-TWIN || FLAMELET || CSTR || FLAMESPEEDCURVE type are allowed!");

    fInput >> dummy;
    if      (dummy == "TIME")		indipendent_variable = 1;
    else if (dummy == "SPACE")		indipendent_variable = 2;
    else if (dummy == "PHI")		indipendent_variable = 3;
	else    ErrorMessage("Only TIME || SPACE || PHI options are allowed!");
    fInput >> x_unit;


    fInput >> dummy;
    if (dummy != "DATA")
        ErrorMessage("Expected DATA key word!");
    fInput >> nData;

    labels = new string[nData+1];
    labels_print = new string[nData+1];
    kind   = new string[nData+1];
    f      = new string[nData+1];
    experimental_data = new BzzMatrix[nData+1];
    ChangeDimensions(nData, &weights);
    ChangeDimensions(nData, &nSum);
    sum_labels.Resize(nData, 10);

    BzzVector aux(2);
    for (i=1;i<=nData;i++)
    {
        fInput >> labels[i];
        labels_print[i] = labels[i];

        if (labels[i] == "SUM")
        {
            labels_print[i] = "SUM";

            fInput >> nSum[i];

            for(int j=1;j<=nSum[i];j++)
            {
                fInput >> dummy;
                sum_labels(i,j) = dummy;
                labels_print[i] += "-" + dummy;
            }
        }


        fInput >> kind[i];
        fInput >> weights[i];
        fInput >> f[i];

        if (kind[i] != "TEMPERATURE" && kind[i] != "MASS_FRACTION" && kind[i] != "MOLE_FRACTION" && kind[i] != "VELOCITY")
            ErrorMessage("Only TEMPERATURE || MASS_FRACTION || MOLE_FRACTION || VELOCITY options are allowed!");

        if (weights[i] < 0.)
            ErrorMessage("Only positive weigths are allowed!");

        if (f[i] != "SQMN" && f[i] != "SQM" && f[i] != "FABS" && f[i] != "FABSN" && f[i] != "CROSS" && f[i] != "MINIMUM" && f[i] != "CROSS_NEW")
            ErrorMessage("Only SQMN, SQM, FABS, FABSS, CROSS, CROSS_NEW, MINIMUM types are allowed!");

        if (f[i] == "MINIMUM" && kindOfReactor != "CFDF-TWIN")
            ErrorMessage("MINIMUM objective function can be used only for CFDF-TWIN reactors...");


        // Reading data
        k=0;
        for (;;)
        {
            fInput >> x;

            if ( x != "//")
            {
                k++;
                fInput >> y;

                if (indipendent_variable == 1)
                    aux[1] = OpenSMOKE_Conversions::conversion_time(atof(x.c_str()), x_unit);

                if (indipendent_variable == 2)
                    aux[1] = OpenSMOKE_Conversions::conversion_length(atof(x.c_str()), x_unit);

                if (indipendent_variable == 3)
                    aux[1] = atof(x.c_str());

                aux[2] = atof(y.c_str());
                experimental_data[i].AppendRow(aux);
            }
            else
                break;
        }
    }
}

BzzVector ExperimentClass::GiveMeObjectiveFunction( ofstream &fLog, string *names, BzzVector &numerical_time, BzzVector &numerical_space,
                                                    BzzVector &numerical_T, BzzVector &numerical_V, BzzMatrix &numerical_mass,
                                                    BzzMatrix &numerical_mole)
{
    int i, j;
    int Nspecies = numerical_mass.Columns();
    
//	BzzCubicSpline spline;
    
	LinearInterpolation spline;

    BzzVector      numerical_x;
    BzzVector      y_expected;

    if (indipendent_variable == 1)  numerical_x = numerical_time;
    if (indipendent_variable == 2)  numerical_x = numerical_space;

	BzzVectorInt iRemove;
	for(i=1;i<=numerical_x.Size()-1;i++)
		if (fabs(numerical_x[i] - numerical_x[i+1]) <= 1.e-16 )
			iRemove.Append(i+1);

	numerical_x.DeleteElements(iRemove);

    BzzVector fObjective(nData);
	fLog << name_object << endl;
    for (i=1;i<=nData;i++)
    {
        ChangeDimensions(experimental_data[i].Rows(), &y_expected);

        // Extracting data corresponding to the experimental data
        if (kind[i] == "TEMPERATURE")
        {
			numerical_T.DeleteElements(iRemove);

			fLog << "Temperature..." << endl;
			fLog << numerical_x.Size() << " " << numerical_T.Size() << endl;
			if (numerical_x.Size() != numerical_T.Size())
				fLog << "Size error" << endl;
			for (j=1;j<=numerical_T.Size();j++)
				fLog << "  " << numerical_x[j] << " " << numerical_T[j] << endl;


			spline(numerical_x, numerical_T);
		}
		else if (kind[i] == "VELOCITY")
		{
			numerical_V.DeleteElements(iRemove);
			spline(numerical_x, numerical_V);
		}
        else if (kind[i] == "MASS_FRACTION")
        {
            if (labels[i] != "SUM")
            {
                int index = recognize_species(names, labels[i], Nspecies);
                BzzVector y = numerical_mass.GetColumn(index);
				y.DeleteElements(iRemove);
                spline(numerical_x, y);
            }
            else
            {
                BzzVector y(numerical_x.Size());
                for(j=1;j<=nSum[i];j++)
                {
                    int index = recognize_species(names, sum_labels(i,j), Nspecies);
                    y += numerical_mass.GetColumn(index);
                }
				y.DeleteElements(iRemove);
                spline(numerical_x, y);
            }

        }
        else if (kind[i] == "MOLE_FRACTION")
        {
            if (labels[i] != "SUM")
            {
                int index = recognize_species(names, labels[i], Nspecies);
                BzzVector y = numerical_mole.GetColumn(index);
				y.DeleteElements(iRemove);
                spline(numerical_x, y);
            }
            else
            {
                BzzVector y(numerical_x.Size());
                for(j=1;j<=nSum[i];j++)
                {
                    int index = recognize_species(names, sum_labels(i,j), Nspecies);
                    y += numerical_mole.GetColumn(index);
                }
				y.DeleteElements(iRemove);
                spline(numerical_x, y);
            }
        }
        else
            ErrorMessage("Only TEMPERATURE || MASS_FRACTION || MOLE_FRACTION || VELOCITY options are possible!");

        // Extracting values corresponding to the experimental support point
        for (j=1;j<=experimental_data[i].Rows();j++)
        {
            y_expected[j] = spline(experimental_data[i][j][1]);

            #ifdef VERBOSE_FLAG
                cout    << experimental_data[i][j][1] << " "
                        << experimental_data[i][j][2] << " "
                        << y_expected[j] << endl;
            #endif
        }

        double sum = 0.;
        if (f[i] == "SQMN")
        {
            for (j=1;j<=experimental_data[i].Rows();j++)
                sum += BzzPow2( (y_expected[j]-experimental_data[i][j][2])/experimental_data[i][j][2] );
        }
        else if (f[i] == "SQM")
        {
            for (j=1;j<=experimental_data[i].Rows();j++)
                sum += BzzPow2( (y_expected[j]-experimental_data[i][j][2]) );
        }
        else if (f[i] == "FABSN")
        {
            for (j=1;j<=experimental_data[i].Rows();j++)
                sum += fabs( (y_expected[j]-experimental_data[i][j][2])/experimental_data[i][j][2] );
        }
        else if (f[i] == "FABS")
        {
            for (j=1;j<=experimental_data[i].Rows();j++)
                sum += fabs( y_expected[j]-experimental_data[i][j][2] );
        }
		else if (f[i] == "CROSS")
        {
            for (j=1;j<=experimental_data[i].Rows();j++)
                sum += BzzPow2( y_expected[j]-experimental_data[i][j][2] );
			sum /= experimental_data[i].Rows();
			sum /= BzzPow2(experimental_data[i].GetColumn(2).Max());
        }
		else if (f[i] == "CROSS_NEW")
        {
			fLog << "Prof: " << i << endl;
            for (j=1;j<=experimental_data[i].Rows();j++)
			{
				fLog << " Data: " << j << " " << y_expected[j] << " " << experimental_data[i][j][2] << endl;
				sum += BzzPow2( y_expected[j]-experimental_data[i][j][2] );
			}
			double denominator = experimental_data[i].Rows() * BzzPow2(experimental_data[i].GetColumn(2).Max());
			
			fLog << "Sum: " << sum << endl;
			fLog << "Den: " << denominator << endl;
			fLog << "Max: " << experimental_data[i].GetColumn(2).Max() << endl;
			fLog << "Row: " << experimental_data[i].Rows() << endl;
			fLog << endl;

			
			if (denominator <= 1.e-20)
			{
				cout << "Sum: " << sum << endl;
				cout << "Den: " << denominator << endl;
				cout << "Max: " << experimental_data[i].GetColumn(2).Max() << endl;
				cout << "Row: " << experimental_data[i].Rows() << endl;
				ErrorMessage("Denominator too small...");
			}

			sum /= denominator;
	    }
		else
            ErrorMessage("Only SQMN || SQM || FABSN || FABS || CROSS || CROSS_NEW options are possible!");

        fObjective[i] = sum*weights[i];

        #ifdef VERBOSE_FLAG
        cout << sum << " " << fObjective << endl;
        #endif
    }

	fLog << endl << "*********************************************************************" << endl;

    return fObjective;
}

BzzVector ExperimentClass::GiveMeObjectiveFunction( ofstream &fLog, BzzVector &numerical_phi, BzzVector &numerical_v)
{
    int i, j;
    
	LinearInterpolation spline;
    BzzVector y_expected;

    BzzVector fObjective(nData);
	
    for (i=1;i<=nData;i++)
    {
		cout << kind[i] << endl;
		if (kind[i] != "VELOCITY")	continue;

        ChangeDimensions(experimental_data[i].Rows(), &y_expected);

        // Extracting data corresponding to the experimental data
		spline(numerical_phi, numerical_v);

        // Extracting values corresponding to the experimental support point
        for (j=1;j<=experimental_data[i].Rows();j++)
            y_expected[j] = spline(experimental_data[i][j][1]);

        double sum = 0.;
        if (f[i] == "SQMN")
        {
            for (j=1;j<=experimental_data[i].Rows();j++)
                sum += BzzPow2( (y_expected[j]-experimental_data[i][j][2])/experimental_data[i][j][2] );

			for (j=1;j<=experimental_data[i].Rows();j++)
				fLog << experimental_data[i][j][1] << " " << y_expected[j] << " " << experimental_data[i][j][2] << " " << sum << endl;
        }
		else
            ErrorMessage("Only SQMN options are possible!");

        fObjective[i] = sum*weights[i];
    }

	fLog << endl << "*********************************************************************" << endl;

    return fObjective;
}

void ExperimentClass::GiveMeRegressionFunction(	string *names, 
												BzzVector &numerical_time, 
												BzzVector &numerical_space,
												BzzVector &numerical_T,
												BzzVector &numerical_V,
												BzzMatrix &numerical_mole,
												
												BzzVector &required_x,
												BzzVector &expected_y
											  )
{
    int i, j;
    BzzCubicSpline spline;
    BzzVector      numerical_x;
    BzzVector      y_expected;

    if (indipendent_variable == 1)  numerical_x = numerical_time;
    if (indipendent_variable == 2)  numerical_x = numerical_space;


    for (i=1;i<=nData;i++)
    {
        ChangeDimensions(required_x.Size(), &expected_y);

        // Extracting data corresponding to the experimental data
        if (kind[i] == "TEMPERATURE")		spline(numerical_x, numerical_T);
		
		else if (kind[i] == "VELOCITY")		spline(numerical_x, numerical_V);
        
		else if (kind[i] == "MOLE_FRACTION")
        {
			int Nspecies = numerical_mole.Columns();
            if (labels[i] != "SUM")
            {
                int index = recognize_species(names, labels[i], Nspecies);
                BzzVector y = numerical_mole.GetColumn(index);
                spline(numerical_x, y);
            }
            else
            {
                BzzVector y(numerical_x.Size());
                for(j=1;j<=nSum[i];j++)
                {
                    int index = recognize_species(names, sum_labels(i,j), Nspecies);
                    y += numerical_mole.GetColumn(index);
                }
                spline(numerical_x, y);
            }
        }
        else
            ErrorMessage("Only TEMPERATURE || MOLE_FRACTION || VELOCITY options are possible!");

        // Extracting values corresponding to the experimental support point
        for (j=1;j<=required_x.Size();j++)
            expected_y[j] = spline(required_x[j]);
	}
}

void ExperimentClass::PostProcessing(string flag, string *names,
                                        BzzVector &numerical_time,   BzzVector &numerical_space,
                                        BzzVector &numerical_T, BzzVector &numerical_V,
                                        BzzMatrix &numerical_mass, BzzMatrix &numerical_mole)
{
    int i, j;
    int Nspecies = numerical_mass.Columns();
    BzzMatrix yMatrix;
    BzzVector aux;
    BzzVector max_vector;
    BzzVector numerical_x;

    ofstream fExperimental;
    ofstream fNumerical;

    string folder = "Output_Optimization/";
    string experimental_name    = folder + folderName + "_experimental.out";
    string numerical_name       = folder + flag + "_" + folderName + "_numerical.out";

    openOutputFileAndControl(fExperimental, experimental_name.c_str());
    openOutputFileAndControl(fNumerical, numerical_name.c_str());
    fExperimental.setf(ios::scientific);
    fNumerical.setf(ios::scientific);

    if (indipendent_variable == 1) numerical_x = numerical_time;
    if (indipendent_variable == 2) numerical_x = numerical_space;

    for (i=1;i<=nData;i++)
    {
        // Extracting data corresponding to the experimental data
        if (kind[i] == "TEMPERATURE")
            yMatrix.AppendColumn(numerical_T);
        else if (kind[i] == "VELOCITY")
            yMatrix.AppendColumn(numerical_V);
        else if (kind[i] == "MASS_FRACTION")
        {
            if (labels[i] != "SUM")
            {
                int index = recognize_species(names, labels[i], Nspecies);
                aux = numerical_mass.GetColumn(index);
                yMatrix.AppendColumn(aux);
            }
            else
            {
                ChangeDimensions(numerical_mass.Rows(), &aux);
                for(j=1;j<=nSum[i];j++)
                {
                    int index = recognize_species(names, sum_labels(i,j), Nspecies);
                    aux += numerical_mass.GetColumn(index);
                }
                yMatrix.AppendColumn(aux);
            }
        }

        else if (kind[i] == "MOLE_FRACTION")
        {
            if (labels[i] != "SUM")
            {
                int index = recognize_species(names, labels[i], Nspecies);
                aux = numerical_mole.GetColumn(index);
                yMatrix.AppendColumn(aux);
            }
            else
            {
                ChangeDimensions(numerical_mole.Rows(), &aux);
                for(j=1;j<=nSum[i];j++)
                {
                    int index = recognize_species(names, sum_labels(i,j), Nspecies);
                    aux += numerical_mole.GetColumn(index);
                }
                yMatrix.AppendColumn(aux);
            }
        }
    }


    for (i=1;i<=nData;i++)
        max_vector.Append(experimental_data[i].Rows());

    for (i=1;i<=nData;i++)
    {
        fExperimental << "x"            << "\t\t";
        fExperimental << labels_print[i]      << "\t\t";

        fNumerical << "x"            << "\t\t";
        fNumerical << labels_print[i]      << "\t\t";
    }
    fExperimental   << endl;
    fNumerical      << endl;


    for (j=1;j<=max_vector.Max();j++)
    {
        for (i=1;i<=nData;i++)
        {
            if (j<=experimental_data[i].Rows())
                fExperimental   << experimental_data[i][j][1] << "\t"
                                << experimental_data[i][j][2] << "\t";
            else
                fExperimental   << "************" << "\t"
                                << "************" << "\t";
        }
        fExperimental << endl;
    }

    for (j=1;j<=numerical_x.Size();j++)
    {
        for (i=1;i<=nData;i++)
                fNumerical   << numerical_x[j]  << "\t"
                             << yMatrix[j][i]   << "\t";
        fNumerical << endl;
    }

    fExperimental.close();
    fNumerical.close();
}

void ExperimentClass::PostProcessing(string flag, BzzVector &numerical_phi,   BzzVector &numerical_T,
                                     BzzVector &numerical_v)
{
    ofstream fExperimental;
    ofstream fNumerical;

    string folder = "Output_Optimization/";
    string experimental_name    = folder + folderName + "_experimental.out";
    string numerical_name       = folder + flag + "_" + folderName + "_numerical.out";

    openOutputFileAndControl(fExperimental, experimental_name.c_str());
    openOutputFileAndControl(fNumerical, numerical_name.c_str());
    fExperimental.setf(ios::scientific);
    fNumerical.setf(ios::scientific);

    fExperimental	<< "phi"	<< "\t\t";
    fExperimental	<< "T[k]"	<< "\t\t";
    fExperimental	<< "phi"	<< "\t\t";
    fExperimental	<< "v[m/s]"	<< "\t\t";
    fExperimental   << endl;

    fNumerical	<< "phi"	<< "\t\t";
    fNumerical	<< "T[k]"	<< "\t\t";
    fNumerical	<< "phi"	<< "\t\t";
    fNumerical	<< "v[m/s]"	<< "\t\t";
    fNumerical  << endl;

	// Experimental
	BzzVector max_vector;
	for (int i=1;i<=nData;i++)
		max_vector.Append(experimental_data[i].Rows());
    for (int j=1;j<=max_vector.Max();j++)
    {
        for (int i=1;i<=nData;i++)
        {
            if (j<=experimental_data[i].Rows())
                fExperimental   << experimental_data[i][j][1] << "\t"
                                << experimental_data[i][j][2] << "\t";
            else
                fExperimental   << "************" << "\t"
                                << "************" << "\t";
        }
        fExperimental << endl;
    }

	// Numerical
    for (int j=1;j<=numerical_phi.Size();j++)
		fNumerical	<< numerical_phi[j]  << "\t"	<< numerical_T[j]   << "\t"
					<< numerical_phi[j]  << "\t"	<< numerical_v[j]   << "\t" 
					<< endl;

    fExperimental.close();
    fNumerical.close();
}

int ExperimentClass::recognize_species(string *names, string label, int Nspecies)
{
    for (int k=1;k<=Nspecies; k++)
        if (names[k] == label)
            return k;

    string message = "This species is not avaliable: " + label;
    ErrorMessage(message);

    return -1;
}

void ExperimentClass::Gnuplot_plots()
{
    for (int i=1;i<=nData;i++)
    {
        string fileFigure = folderName + "_" + labels_print[i];
        string fileData_1 = "Output_Optimization/" + folderName + "_experimental.out";
        string fileData_2 = "Output_Optimization/Start_" + folderName + "_numerical.out";
        string fileData_3 = "Output_Optimization/Final_" + folderName + "_numerical.out";

        OpenSMOKE_GnuPlotInterface gplot("Output_Optimization");
        if (indipendent_variable == 1)
            gplot.setPlot(folderName + " - " + labels_print[i], "time [s]", kind[i]);
        else
            gplot.setPlot(folderName + " - " + labels_print[i], "space [m]", kind[i]);
        gplot.setKey("Experimental", "Original", "Optimized");
        gplot.setKind('p', 'l', 'l');

		if (kindOfReactor == "PFR")		gplot.setLogScaleX();
		
		cout << "Gnuplot: " << fileFigure << endl;
        gplot.plot(fileFigure, fileData_1, fileData_2, fileData_3, (i-1)*2+1, (i-1)*2+2);
    }
}

void ExperimentClass::Gnuplot_plots_online()
{
    for (int i=1;i<=nData;i++)
    {
        string fileFigure = folderName + "_" + labels_print[i];
        string fileData_1 = "Output_Optimization/" + folderName + "_experimental.out";
        string fileData_2 = "Output_Optimization/Start_" + folderName + "_numerical.out";
        string fileData_3 = "Output_Optimization/online-post-processing_" + folderName + "_numerical.out";

        OpenSMOKE_GnuPlotInterface gplot("Output_Optimization");
        if (indipendent_variable == 1)
            gplot.setPlot(folderName + " - " + labels_print[i], "time [s]", kind[i]);
        else
            gplot.setPlot(folderName + " - " + labels_print[i], "space [m]", kind[i]);
        gplot.setKey("Experimental", "Original", "Optimized");
        gplot.setKind('p', 'l', 'l');

		if (kindOfReactor == "PFR")		gplot.setLogScaleX();
		
		cout << "Gnuplot: " << fileFigure << endl;
        gplot.plot(fileFigure, fileData_1, fileData_2, fileData_3, (i-1)*2+1, (i-1)*2+2);
    }
}


void ExperimentClass::Latex_figure(OpenSMOKE_LatexInterface &latex)
{
    for (int i=1;i<=nData;i++)
    {
        string fileName = "Output_Optimization/" + folderName + "_" + labels_print[i] + ".eps";
        string caption  = labels_print[i] + "(" + folderName + ")";
        latex.include_figure(fileName, caption);
    }
}

BzzVector ExperimentClass::ExtractTemperatures()
{
	for (int i=1;i<=nData;i++)
        if (kind[i] == "TEMPERATURE")
			return experimental_data[i].GetColumn(2);

	ErrorMessage("TEMPERATURE was not found!");
	return -1;
}

BzzVector ExperimentClass::ExtractSupportTemperatures()
{
	for (int i=1;i<=nData;i++)
        if (kind[i] == "TEMPERATURE")
			return experimental_data[i].GetColumn(1);

	ErrorMessage("TEMPERATURE (support) was not found!");
	return -1;
}

BzzVector ExperimentClass::GiveMeObjectiveFunction( ofstream &fLog, string *names, double numerical_T, double numerical_V, 
													BzzVector &numerical_mass, BzzVector &numerical_mole)
{
	int Nspecies = numerical_mass.Size();
    BzzVector fObjective(nData);
    
	for (int i=1;i<=nData;i++)
    {
		double y_expected;

		// Extracting data corresponding to the experimental data
        if (kind[i] == "TEMPERATURE")
			y_expected = numerical_T;
		else if (kind[i] == "VELOCITY")
			y_expected = numerical_V;
        else if (kind[i] == "MASS_FRACTION")
        {
            if (labels[i] != "SUM")
            {
                int index = recognize_species(names, labels[i], Nspecies);
                y_expected = numerical_mass[index];
            }
            else
            {
                y_expected = 0.;
                for(int j=1;j<=nSum[i];j++)
                {
                    int index = recognize_species(names, sum_labels(i,j), Nspecies);
                    y_expected += numerical_mass[index];
                }
            }
        }
        else if (kind[i] == "MOLE_FRACTION")
        {
            if (labels[i] != "SUM")
            {
                int index = recognize_species(names, labels[i], Nspecies);
				y_expected = numerical_mole[index];
            }
            else
            {
                y_expected = 0.;
                for(int j=1;j<=nSum[i];j++)
                {
                    int index = recognize_species(names, sum_labels(i,j), Nspecies);
                    y_expected += numerical_mole[index];
                }
            }
        }
        else
            ErrorMessage("Only TEMPERATURE || MASS_FRACTION || MOLE_FRACTION || VELOCITY options are possible!");

        double sum = 0.;
        if (f[i] == "SQMN")
                sum  = BzzPow2( (y_expected-experimental_data[i][1][2])/experimental_data[i][1][2] );
        else if (f[i] == "SQM")
                sum  = BzzPow2( (y_expected-experimental_data[i][1][2]) );
        else if (f[i] == "FABSN")
                sum  = fabs( (y_expected-experimental_data[i][1][2])/experimental_data[i][1][2] );
        else if (f[i] == "FABS")
                sum  = fabs( y_expected-experimental_data[i][1][2] );
		else if (f[i] == "CROSS")
        {
			sum  = BzzPow2( y_expected-experimental_data[i][1][2] );
			sum /= experimental_data[i].Rows();
			sum /= BzzPow2(experimental_data[i].GetColumn(2).Max());
        }
		else
            ErrorMessage("Only SQMN || SQM || FABSN || FABS || CROSS options are possible!");

        fObjective[i] = sum*weights[i];
    }

    return fObjective;
}