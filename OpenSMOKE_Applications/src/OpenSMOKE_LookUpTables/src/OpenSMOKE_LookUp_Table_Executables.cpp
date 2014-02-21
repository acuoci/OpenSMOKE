/***************************************************************************
 *   Copyright (C) 2006 by Alberto Cuoci								   *
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

#include "OpenSMOKE_LookUp_Table_Executables.h"

////////////////////////////////////////////////////////////////////////////////////////////////////////////
//									OpenSMOKE_LookUp_Table_Executables						              //
////////////////////////////////////////////////////////////////////////////////////////////////////////////

void OpenSMOKE_LookUp_Table_Executables::ErrorMessage(const string message)
{
    cout << endl;
    cout << "Class:  OpenSMOKE_LookUp_Table_Executables"	<< endl;
    cout << "Object: " << name_object						<< endl;
    cout << "Error:  " << message							<< endl;
    cout << "Press a key to continue... "					<< endl;
    getchar();
    exit(-1);
}

void OpenSMOKE_LookUp_Table_Executables::WarningMessage(const string message)
{
    cout << endl;
    cout << "Class:  OpenSMOKE_LookUp_Table_Executables"	<< endl;
    cout << "Object: "		<< name_object				<< endl;
    cout << "Warning:  "	<< message					<< endl;
    cout << "Press a key to continue... "				<< endl;
    getchar();
}

OpenSMOKE_LookUp_Table_Executables::OpenSMOKE_LookUp_Table_Executables()
{
	name_object						= "[not assigned]";

	iSootSourceTerms_Adiabatic		= false;
	iSootSourceTerms_NonAdiabatic	= false;

	iChiPDF				= 0;					// 0=delta dirac  1=log-normal distribution
	interactions_mode	= INTERACTIONS_NONE;	
}

void OpenSMOKE_LookUp_Table_Executables::SetSootClosureModel(const string _iSootMode)
{
		 if (_iSootMode == "CORRELATED")	interactions_mode	= INTERACTIONS_CORRELATED;
	else if (_iSootMode == "UNCORRELATED")	interactions_mode	= INTERACTIONS_UNCORRELATED;
	else if (_iSootMode == "BRUTAL")		interactions_mode	= INTERACTIONS_BRUTAL;
	else ErrorMessage("Wrong soot closure model...");
}

void OpenSMOKE_LookUp_Table_Executables::SetSootSourcesModel(const string _iSourcesModel)
{
		 if (_iSourcesModel == "ADIABATIC")		iSootSourceTerms_Adiabatic = true;
	else if (_iSourcesModel == "NON_ADIABATIC")	iSootSourceTerms_NonAdiabatic = true;
	else ErrorMessage("Wrong soot sources model...");
}

void OpenSMOKE_LookUp_Table_Executables::SetLogPDF()
{
	iChiPDF = 1;
}


void OpenSMOKE_LookUp_Table_Executables::DefineFromFile(const string inputFile)
{
 //   double			double_value;
    string			string_value;
 //   int				int_value;
	vector<string>  string_vector;

    OpenSMOKE_Dictionary_LookUp_Table_Executables dictionary;
    dictionary.ParseFile(inputFile);

	if (dictionary.Return("#SootClosure", string_value))
		SetSootClosureModel(string_value);

	if (dictionary.Return("#SootSources", string_value))
		SetSootSourcesModel(string_value);

	if (dictionary.Return("#LogPDF"))
		SetLogPDF();
}

void OpenSMOKE_LookUp_Table_Executables::Run()
{
	if (iSootSourceTerms_Adiabatic		== true)		Exe_SootSourceTerms_Adiabatic();
	if (iSootSourceTerms_NonAdiabatic	== true)		Exe_SootSourceTerms_NonAdiabatic();
}

void OpenSMOKE_LookUp_Table_Executables::Exe_SootSourceTerms_Adiabatic()
{
	// ----------------------------------------------------------------------------------
	// EXECUTABLE FILE FOR FLUENT - SOOT SOURCE TERMS ADIABATIC
	// ----------------------------------------------------------------------------------
	// 1. UNCORRELATED CLOSURE: the soot source terms for the corresponding transport
	//							equations are calculated using the uncorrelated closure
	//							model (see Brookes and Moss for further details)
	// 2.   CORRELATED CLOSURE: the soot source terms for the corresponding transport
	//							equations are calculated using the correlated closure
	//							model (see Brookes and Moss for further details); in this
	//                          it is necessary also to call the function for the
	//							calculation of normalized profiles

		SootSourceTerms_Adiabatic();
		exit(1);
}

void OpenSMOKE_LookUp_Table_Executables::Exe_SootSourceTerms_NonAdiabatic()
{
	// ----------------------------------------------------------------------------------
	// EXECUTABLE FILE FOR FLUENT - SOOT SOURCE TERMS WITH ENTHALPY DEFECT
	// ----------------------------------------------------------------------------------
	// 1. UNCORRELATED CLOSURE: the soot source terms for the corresponding transport
	//							equations are calculated using the uncorrelated closure
	//							model (see Brookes and Moss for further details)
	// 2.   CORRELATED CLOSURE: the soot source terms for the corresponding transport
	//							equations are calculated using the correlated closure
	//							model (see Brookes and Moss for further details); in this
	//                          it is necessary also to call the function for the
	//							calculation of normalized profiles

	SootSourceTerms_NonAdiabatic();
	exit(1);
}

void OpenSMOKE_LookUp_Table_Executables::Exe_SootProfiles_NonAdiabatic()
{
	// ----------------------------------------------------------------------------------
	// EXECUTABLE FILE FOR FLUENT - SOOT PREDICTIONS
	// ----------------------------------------------------------------------------------
	// 1. NORMAL PROFILES MODE: the function is used to obtain the values of the 
	//                          normalized profiles of soot properties (this step is 
	//							necessary when the Perfectly Correlated Approach (3)
	//							is adopted)
	// 2. BRUTAL FLAMELET MODE: the function is used to directly calculate the local
	//							values of soot properties; this case correspond to the 
	//							Brutal Flamelet Approach
	
	SootProfiles_NonAdiabatic();
	exit(1);
}


////////////////////////////////////////////////////////////////////////////////////////////////////////////
//											DICTIONARY										              //
////////////////////////////////////////////////////////////////////////////////////////////////////////////

OpenSMOKE_Dictionary_LookUp_Table_Executables::OpenSMOKE_Dictionary_LookUp_Table_Executables()
{
    SetupBase();

	Add("#SootSources",					'O', 'S', "Soot Sources: ADIABATIC || NON_ADIABATIC");
	Add("#SootClosure",					'O', 'S', "Soot closure model: CORRELATED || UNCORRELATED");
	Add("#LogPDF",						'O', 'N', "Logarithmic PDF for scalar dissipation rate");

    Lock();
}


void OpenSMOKE_LookUp_Table_Executables::SootSourceTerms_Adiabatic()
{
	string fileName;
	string rootName;

	int numberOfCells;
	
	char charNucleation[20];
	char charGrowth[20];
	char charAggregation[20];
	char charOxidation[20];
	string stringNucleation;
	string stringGrowth;
	string stringAggregation;
	string stringOxidation;
	
	int nInfo=4;

	double startTime;
	double endTime;
	ifstream fInput;
	ofstream fOutput;
	BzzVector info(nInfo);
	BzzMatrix infoMatrix;
	BzzVector Source;

	BzzLoad load;
	BzzSave save;
	
	startTime =	BzzGetCpuTime();

	// File Management
	fInput.open("FieldForFlamelet.bin", ios::in | ios::binary);
	if (!fInput.is_open()) 
		BzzError("File FieldForFlamelet.bin cannot be opened!!");

	// File Management - OutputFile
	if      (interactions_mode == INTERACTIONS_UNCORRELATED)	fileName = "Uncorrelated_Closure_Source_Terms.bin";
	else if (interactions_mode == INTERACTIONS_CORRELATED)		fileName = "Correlated_Closure_Source_Terms.bin";
	else ErrorMessage("Wrong Option in the extraction of Soot Source Terms with Enthalpy Defect!");

	fOutput.open(fileName.c_str(), ios::out | ios::binary);
	if (!fOutput.is_open()) 
		ErrorMessage("File FlameletSourceTerms.bin cannot be opened!!");

	// ----------------------------------------------------------------------------
	// Memory Allocation
	// ----------------------------------------------------------------------------
	if( !fInput.read( reinterpret_cast < char * > (&numberOfCells), sizeof(int)) )
		ErrorMessage("Error in file FieldForFlamelet.bin: Wrong number of cells!!");

	if( !fInput.read( reinterpret_cast < char * > (&charNucleation), sizeof(charNucleation)) )
		ErrorMessage("Error in file FieldForFlamelet.bin: Wrong iNucleation!!");
		stringNucleation = charNucleation;

	if( !fInput.read( reinterpret_cast < char * > (&charGrowth), sizeof(charGrowth)) )
		ErrorMessage("Error in file FieldForFlamelet.bin: Wrong iGrowth!!");
		stringGrowth = charGrowth;

	if( !fInput.read( reinterpret_cast < char * > (&charAggregation), sizeof(charAggregation)) )
		ErrorMessage("Error in file FieldForFlamelet.bin: Wrong iAggregation!!");
		stringAggregation = charAggregation;

	if( !fInput.read( reinterpret_cast < char * > (&charOxidation), sizeof(charOxidation)) )
		ErrorMessage("Error in file FieldForFlamelet.bin: Wrong iOxidation!!");
		stringOxidation = charOxidation;

	ChangeDimensions(numberOfCells, nInfo, &infoMatrix);
	ChangeDimensions(numberOfCells, &Source);

	sourceField_flamelet_library field_Nucleation;
	sourceField_flamelet_library field_Growth;
	sourceField_flamelet_library field_Aggregation;
	sourceField_flamelet_library field_Oxidation;
	
	// Print Info On Screen
	cout << "------------------------------------------------" << endl;
	cout << " Soot Source Term Calculations From FLUENT data " << endl;
	cout << "------------------------------------------------" << endl;
	cout << "Number Of CFD Cells:         " << numberOfCells << endl;
	cout << "Scalar Dissipation rate PDF: " << iChiPDF << endl;
	cout << "Calculation Mode:            " << interactions_mode << endl;


	// Opening Files
	if ( interactions_mode == INTERACTIONS_UNCORRELATED ) rootName = "Kinetics_Uncorrelated/H0/";
	if ( interactions_mode == INTERACTIONS_CORRELATED   ) rootName = "Kinetics_Correlated/H0/";

	if (stringNucleation != "None")
	{
		string list_names = rootName + "Nucleation_" + stringNucleation + ".source";
		field_Nucleation.read_sourceField_fromFile(list_names.c_str());
	}

	if (stringGrowth != "None")
	{
		string list_names = rootName + "Growth_" + stringGrowth + ".source";
		field_Growth.read_sourceField_fromFile(list_names.c_str());
	}

	if (stringAggregation != "None")
	{
		string list_names = rootName + "Aggregation_" + stringAggregation + ".source";
		field_Aggregation.read_sourceField_fromFile(list_names.c_str());
	}

	if (stringOxidation != "None")
	{
		string list_names = rootName + "Oxidation_" + stringOxidation + ".source";
		field_Oxidation.read_sourceField_fromFile(list_names.c_str());
	}

	// Recovering information from file - Binary Version
	for(int cell=1;cell<=numberOfCells;cell++)
	{	
		load.readVectorFromBinaryFile(fInput, nInfo, info);	// moments
		infoMatrix.SetRow(cell, info);
	}
	fInput.close();

	cout << "Max mixture fraction:          " << infoMatrix.GetColumn(1).Max() << endl;
	cout << "Min mixture fraction:          " << infoMatrix.GetColumn(1).Min() << endl;
	cout << "Max mixture fraction variance: " << infoMatrix.GetColumn(2).Max() << endl;
	cout << "Min mixture fraction variance: " << infoMatrix.GetColumn(2).Min() << endl;
	cout << "Max scalar dissipation rate:   " << infoMatrix.GetColumn(3).Max() << endl;
	cout << "Min scalar dissipation rate:   " << infoMatrix.GetColumn(3).Min() << endl;


	// Writing On File - Nucleation Source Term
	if (stringNucleation != "None")
	{
		cout << "Nucleation...";
		CalculateSourceTerms(Source, field_Nucleation, infoMatrix, numberOfCells);
		save.writeVectorOnBinaryFile(fOutput, Source);
		cout << " DONE!" << endl;
	}

	// Writing On File - Growth Source Term
	if (stringGrowth != "None")		
	{
		cout << "Growth...";
		CalculateSourceTerms(Source, field_Growth, infoMatrix, numberOfCells);
		save.writeVectorOnBinaryFile(fOutput, Source);
		cout << " DONE!" << endl;
	}

	// Writing On File - Aggregation Source Term
	if (stringAggregation != "None")				
	{
		cout << "Aggregation...";
		CalculateSourceTerms(Source, field_Aggregation, infoMatrix, numberOfCells);
		save.writeVectorOnBinaryFile(fOutput, Source);
		cout << " DONE!" << endl;
	}

	// Writing On File - Oxidation Source Term
	if (stringOxidation != "None")				
	{
		cout << "Oxidation...";
		CalculateSourceTerms(Source, field_Oxidation, infoMatrix, numberOfCells);
		save.writeVectorOnBinaryFile(fOutput, Source);
		cout << " DONE!" << endl;
	}

	// Closing the Output File
	fOutput.close();

	// Info On Screen
	endTime = BzzGetCpuTime() - startTime;
	cout << "Total time: " << endTime << " s" << endl;
}

void OpenSMOKE_LookUp_Table_Executables::SootSourceTerms_NonAdiabatic()
{
	// 1. UNCORRELATED MODE: in this case the function is used to calcultae the source
	//                       terms in the transport equations for soot using the uncorrelated approach;
	//                       therefore it is not necessary to call the function for the
	//                       normalization of soot properties profiles
	// 2. CORRELATED MODE:   in this case the function is used to calcultae the source
	//                       terms in the transport equations for soot using the correlated
	//                       approach; now it is not necessary to call the function fo the
	//                       normalization of soot properties profiles

	int j;
	string fileName;
	string rootName;
	vector<string> list_names;
	int numberOfCells;
	
	int iNucleation;
	int iGrowth;
	int iAggregation;
	int iOxidation;

	int nInfo				= 4;
	int nEnthalpyDefects	= 5;
	double defectStep		= 50.;							// [kJ/kg, heat loss]
	double defectMax = defectStep*(nEnthalpyDefects-1);		// [kJ/kg, heat loss]

	double startTime;
	double endTime;
	ifstream fInput;
	ofstream fOutput;
	BzzVector info(nInfo);
	BzzMatrix infoMatrix;
	BzzVector Source;

	BzzLoad load;
	BzzSave save;
	
	startTime =	BzzGetCpuTime();

	// File Management - File From CFD Simulation
	fInput.open("FieldForFlamelet.bin", ios::in | ios::binary);
	if (!fInput.is_open()) 
		BzzError("File FieldForFlamelet.bin cannot be opened!!");

	// File Management - OutputFile
		 if (interactions_mode == INTERACTIONS_UNCORRELATED)	fileName = "Uncorrelated_Closure_Source_Terms.bin";
	else if (interactions_mode == INTERACTIONS_CORRELATED)		fileName = "Correlated_Closure_Source_Terms.bin";
	fOutput.open(fileName.c_str(), ios::out | ios::binary);
	if (!fOutput.is_open())	ErrorMessage(fileName + "file cannot be opened!!");

	// Memory Allocation
	if( !fInput.read( reinterpret_cast < char * > (&numberOfCells), sizeof(int)) )
		ErrorMessage("Error in file FieldForFlamelet.bin: Wrong number of cells!!");
	if( !fInput.read( reinterpret_cast < char * > (&iNucleation), sizeof(int)) )
		ErrorMessage("Error in file FieldForFlamelet.bin: Wrong iNucleation!!");
	if( !fInput.read( reinterpret_cast < char * > (&iGrowth), sizeof(int)) )
		ErrorMessage("Error in file FieldForFlamelet.bin: Wrong iGrowth!!");
	if( !fInput.read( reinterpret_cast < char * > (&iAggregation), sizeof(int)) )
		ErrorMessage("Error in file FieldForFlamelet.bin: Wrong iAggregation!!");
	if( !fInput.read( reinterpret_cast < char * > (&iOxidation), sizeof(int)) )
		ErrorMessage("Error in file FieldForFlamelet.bin: Wrong iOxidation!!");

	// Info :	1=mixture fraction 2=variance of mixture fraction
	//			3=scalar dissipation rate 4=enthalpy defect [ positive if heat loss, kJ/kg]
	ChangeDimensions(numberOfCells, nInfo, &infoMatrix);
	
	ChangeDimensions(numberOfCells, &Source);

	// Opening Library for nucleation rate
	sourceField_flamelet_library *field_Nucleation;
	sourceField_flamelet_library *field_Growth;
	sourceField_flamelet_library *field_Oxidation;
	sourceField_flamelet_library *field_Aggregation;

	field_Nucleation	= new sourceField_flamelet_library[nEnthalpyDefects];
	field_Growth		= new sourceField_flamelet_library[nEnthalpyDefects];
	field_Oxidation		= new sourceField_flamelet_library[nEnthalpyDefects];
	field_Aggregation	= new sourceField_flamelet_library[nEnthalpyDefects];

	if (interactions_mode == INTERACTIONS_UNCORRELATED) rootName = "Kinetics_Uncorrelated/H";
	if (interactions_mode == INTERACTIONS_CORRELATED)	rootName = "Kinetics_Correlated/H";

	if (iNucleation != 0)
	{
		list_names = construct_name_list(rootName, "Nucleation", iNucleation);
		for (j=0;j<nEnthalpyDefects;j++)
			field_Nucleation[j].read_sourceField_fromFile(list_names[j].c_str());
	}

	if (iGrowth != 0)
	{
		list_names = construct_name_list(rootName, "Growth", iGrowth);
		for (j=0;j<nEnthalpyDefects;j++)
			field_Growth[j].read_sourceField_fromFile(list_names[j].c_str());
	}

	if (iAggregation != 0)
	{
		list_names = construct_name_list(rootName, "Aggregation", iAggregation);
		for (j=0;j<nEnthalpyDefects;j++)
			field_Aggregation[j].read_sourceField_fromFile(list_names[j].c_str());
	}

	if (iOxidation != 0)
	{
		list_names = construct_name_list(rootName, "Oxidation", iOxidation);
		for (j=0;j<nEnthalpyDefects;j++)
			field_Oxidation[j].read_sourceField_fromFile(list_names[j].c_str());
	}

	// Recovering information from file - Binary Version
	for(int cell=1;cell<=numberOfCells;cell++)
	{	
		load.readVectorFromBinaryFile(fInput, nInfo, info);	
		infoMatrix.SetRow(cell, info);
	}
	fInput.close();

	// Print Info On Screen
	cout << "------------------------------------------------" << endl;
	cout << " Soot Source Term Calculations From FLUENT data " << endl;
	cout << "------------------------------------------------" << endl;
	cout << "Number Of CFD Cells:         " << numberOfCells << endl;
	cout << "Scalar Dissipation rate PDF: " << iChiPDF << endl;
	cout << "Calculation Mode:            " << interactions_mode << endl;


	// Writing On File - Nucleation Source Term
	if (iNucleation != 0)
	{
		cout << "Nucleation...";
		CalculateSourceTerms(Source, field_Nucleation, infoMatrix, numberOfCells, nEnthalpyDefects, defectStep);
		save.writeVectorOnBinaryFile(fOutput, Source);
		cout << " DONE!" << endl;
	}

	// Writing On File - Growth Source Term
	if (iGrowth != 0)		
	{
		cout << "Growth...";
		CalculateSourceTerms(Source, field_Growth, infoMatrix, numberOfCells, nEnthalpyDefects, defectStep);
		save.writeVectorOnBinaryFile(fOutput, Source);
		cout << " DONE!" << endl;
	}

	// Writing On File - Aggregation Source Term
	if (iAggregation != 0)				
	{
		cout << "Aggregation... ";
		CalculateSourceTerms(Source, field_Aggregation, infoMatrix, numberOfCells, nEnthalpyDefects, defectStep);
		save.writeVectorOnBinaryFile(fOutput, Source);
		cout << " DONE!" << endl;
	}

	// Writing On File - Oxidation Source Term
	if (iOxidation != 0)				
	{
		cout << "Oxidation... ";
		CalculateSourceTerms(Source, field_Oxidation, infoMatrix, numberOfCells, nEnthalpyDefects, defectStep);
		save.writeVectorOnBinaryFile(fOutput, Source);
		cout << " DONE!" << endl;
	}

	// Closure Of Output File
	fOutput.close();

	// Information On Screen
	endTime = BzzGetCpuTime() - startTime;
	cout << "Total time: " << endTime << " s" << endl;
}



void OpenSMOKE_LookUp_Table_Executables::SootProfiles_NonAdiabatic()
{
	// 1. NORMAL PROFILES MODE
	//	  To obtain the normalized values of soot properties from the data extracted from the CFD code; 
	//    in this case the file FieldForFlamelet.bin must be read from the external CFD code and the
	//    local values of m0_N and fv_N are written on the file
	// 2. BRUTAL FLAMELET (Approach 4, according to Brookes and Moss, 1999)
	//	  To directly obtain the soot properties mean values, without normalizing the profiles
	//    In this case the the file FieldForFlamelet.bin must be read from the 
	//    external CFD code and the local values of m0 and fv are written on the file

	int numberOfCells;
	int cell;

	int nInfo=4;
	int nEnthalpyDefects=5;
	double defectStep = 50.;								// [kJ/kg, heat loss]
	double defectMax = defectStep*(nEnthalpyDefects-1);		// [kJ/kg, heat loss]

	double startTime;
	double endTime;
	double dummy;
	ifstream fInput;
	ofstream fOutput;
	BzzVector info(nInfo);
	BzzMatrix infoMatrix;
	BzzVector Source;

	BzzLoad load;
	BzzSave save;
	
	startTime =	BzzGetCpuTime();

	// ---------------------------------------------------
	// File Management - Open CFD File
	// ---------------------------------------------------
	fInput.open("FieldForFlamelet.bin", ios::in | ios::binary);
	if (!fInput.is_open()) 
		BzzError("File FieldForFlamelet.bin cannot be opened!!");

	// ---------------------------------------------------
	// File Management - Open Output File
	// ---------------------------------------------------
	// Perfect Correlation (with normalized profiles, see Brookes and Moss 1999)
	if (interactions_mode == INTERACTIONS_CORRELATED)
		fOutput.open("NormalProfiles.bin", ios::out | ios::binary);

	// Perfect Uncorrelation (see Brookes and Moss, 1999)
	else if (interactions_mode == INTERACTIONS_BRUTAL)
		fOutput.open("Brutal_Flamelet.bin", ios::out | ios::binary);
	
	if (!fOutput.is_open()) 
		ErrorMessage("File FlameletSootResults.bin cannot be opened!!");

	// ----------------------------------------------------------------------------
	// Memory Allocation
	// ----------------------------------------------------------------------------
	if( !fInput.read( reinterpret_cast < char * > (&numberOfCells), sizeof(int)) )
		ErrorMessage("Error in file FieldForFlamelet.bin: Wrong number of cells!!");
	if( !fInput.read( reinterpret_cast < char * > (&dummy), sizeof(int)) )
		ErrorMessage("Error in file FieldForFlamelet.bin: Wrong iNucleation!!");
	if( !fInput.read( reinterpret_cast < char * > (&dummy), sizeof(int)) )
		ErrorMessage("Error in file FieldForFlamelet.bin: Wrong iGrowth!!");
	if( !fInput.read( reinterpret_cast < char * > (&dummy), sizeof(int)) )
		ErrorMessage("Error in file FieldForFlamelet.bin: Wrong iAggregation!!");
	if( !fInput.read( reinterpret_cast < char * > (&dummy), sizeof(int)) )
		ErrorMessage("Error in file FieldForFlamelet.bin: Wrong iOxidation!!");

	ChangeDimensions(numberOfCells, nInfo, &infoMatrix);
	ChangeDimensions(numberOfCells, &Source);

	// ---------------------------------------------------------------
	// Opening Soot Flamelet Library
	// ---------------------------------------------------------------
	sourceField_flamelet_library *field_A;
	sourceField_flamelet_library *field_B;
	field_A	= new sourceField_flamelet_library[nEnthalpyDefects];
	field_B	= new sourceField_flamelet_library[nEnthalpyDefects];

	if (interactions_mode == INTERACTIONS_CORRELATED)
	{
		field_A[0].read_sourceField_fromFile(	"Kinetics_Correlated/H0/NormalProfiles/m0_N.source");
		field_A[1].read_sourceField_fromFile(	"Kinetics_Correlated/H50/NormalProfiles/m0_N.source");
		field_A[2].read_sourceField_fromFile(	"Kinetics_Correlated/H100/NormalProfiles/m0_N.source");
		field_A[3].read_sourceField_fromFile(	"Kinetics_Correlated/H150/NormalProfiles/m0_N.source");
		field_A[4].read_sourceField_fromFile(	"Kinetics_Correlated/H200/NormalProfiles/m0_N.source");

		field_B[0].read_sourceField_fromFile(	"Kinetics_Correlated/H0/NormalProfiles/fv_N.source");
		field_B[1].read_sourceField_fromFile(	"Kinetics_Correlated/H50/NormalProfiles/fv_N.source");
		field_B[2].read_sourceField_fromFile(	"Kinetics_Correlated/H100/NormalProfiles/fv_N.source");
		field_B[3].read_sourceField_fromFile(	"Kinetics_Correlated/H150/NormalProfiles/fv_N.source");
		field_B[4].read_sourceField_fromFile(	"Kinetics_Correlated/H200/NormalProfiles/fv_N.source");
	}

	if (interactions_mode == INTERACTIONS_BRUTAL)
	{
		field_A[0].read_sourceField_fromFile(	"Kinetics_Brutal/LargeBins_H0.soot");
		field_A[1].read_sourceField_fromFile(	"Kinetics_Brutal/LargeBins_H50.soot");
		field_A[2].read_sourceField_fromFile(	"Kinetics_Brutal/LargeBins_H100.soot");
		field_A[3].read_sourceField_fromFile(	"Kinetics_Brutal/LargeBins_H150.soot");
		field_A[4].read_sourceField_fromFile(	"Kinetics_Brutal/LargeBins_H200.soot");

		field_B[0].read_sourceField_fromFile(	"Kinetics_Brutal/SmallBins_H0.soot");
		field_B[1].read_sourceField_fromFile(	"Kinetics_Brutal/SmallBins_H50.soot");
		field_B[2].read_sourceField_fromFile(	"Kinetics_Brutal/SmallBins_H100.soot");
		field_B[3].read_sourceField_fromFile(	"Kinetics_Brutal/SmallBins_H150.soot");
		field_B[4].read_sourceField_fromFile(	"Kinetics_Brutal/SmallBins_H200.soot");
	}

	// ---------------------------------------------------------------
	// Recovering information from file - Binary Version
	// ---------------------------------------------------------------
	for(cell=1;cell<=numberOfCells;cell++)
	{	
		load.readVectorFromBinaryFile(fInput, nInfo, info);
		infoMatrix.SetRow(cell, info);
	}
	fInput.close();

	// Print Info On Screen
	cout << "-----------------------------------------------" << endl;
	cout << "  Soot Calculations From FLUENT data           " << endl;
	cout << "-----------------------------------------------" << endl;
	cout << "Number Of CFD Cells:         " << numberOfCells << endl;
	cout << "Scalar Dissipation rate PDF: " << iChiPDF << endl;
	cout << "Calculation Mode:            " << interactions_mode << endl;

	// Field A: (Normalized Profiles = m0  -  Brutal Flamelet = LargeBins)
	CalculateSourceTerms(Source, field_A, infoMatrix, numberOfCells, nEnthalpyDefects, defectStep);
	save.writeVectorOnBinaryFile(fOutput, Source);
		
	if (interactions_mode == INTERACTIONS_CORRELATED)	cout << "Max Normalized Soot Particle Number Density = " << Source.Max() << endl;
	if (interactions_mode == INTERACTIONS_BRUTAL)		cout << "Max Volume Fraction for Large Bins = "  << Source.Max() << endl;

	
	// -------------------------------------------------------------------------------
	// Field B: (Normalized fdfles = fv  -  Brutal Flamelet = SmallBins)
	// -------------------------------------------------------------------------------	
	CalculateSourceTerms(Source, field_B, infoMatrix, numberOfCells, nEnthalpyDefects, defectStep);
	save.writeVectorOnBinaryFile(fOutput, Source);

	if (interactions_mode == INTERACTIONS_CORRELATED)	cout << "Max Normalized Soot Volume Fraction = " << Source.Max() << endl;	
	if (interactions_mode == INTERACTIONS_BRUTAL)		cout << "Max Volume Fraction for Small Bins = "  << Source.Max() << endl;

	// ---------------------------------------------------------------
	// Close the Output File
	// ---------------------------------------------------------------
	fOutput.close();

	// ---------------------------------------------------------------
	// Information on screen
	// ---------------------------------------------------------------
	endTime = BzzGetCpuTime() - startTime;
	cout << "Total time: " << endTime << " s" << endl;
}


void OpenSMOKE_LookUp_Table_Executables::CalculateSourceTerms(	BzzVector &Source, sourceField_flamelet_library *field, 
																BzzMatrix &infoMatrix, int numberOfCells, int nEnthalpyDefects, double defectStep)
{
	int index;
	double sourceA, sourceB;
	
	if (iChiPDF == 0)
	{
		for(int cell=1;cell<=numberOfCells;cell++)	
		{
			index = int(infoMatrix[cell][4]/defectStep);
			
			if (infoMatrix[cell][4]<=0.)		// Heat Gain
				Source[cell] =  field[0].find_interpolated_value_global
								(infoMatrix[cell][1], infoMatrix[cell][2], infoMatrix[cell][3]); 
		
			else if (index<(nEnthalpyDefects-1))	// linear interpolation
			{
				sourceA = field[index].find_interpolated_value_global
						 (infoMatrix[cell][1], infoMatrix[cell][2], infoMatrix[cell][3]); 
				sourceB = field[index+1].find_interpolated_value_global
						 (infoMatrix[cell][1], infoMatrix[cell][2], infoMatrix[cell][3]);

				double xA = index*defectStep;
				double xB = (index+1)*defectStep;
				Source[cell] =  linear_interpolation_for_enthalpy_defect(infoMatrix[cell][4], sourceA, sourceB, xA, xB);
			}
						
			else if (index>=(nEnthalpyDefects-1))	// Heat Loss
			{
				sourceA = field[nEnthalpyDefects-2].find_interpolated_value_global
						  (infoMatrix[cell][1], infoMatrix[cell][2], infoMatrix[cell][3]); 
				sourceB = field[nEnthalpyDefects-1].find_interpolated_value_global
						  (infoMatrix[cell][1], infoMatrix[cell][2], infoMatrix[cell][3]);

				double xA = (nEnthalpyDefects-2)*defectStep;
				double xB = (nEnthalpyDefects-1)*defectStep;
				Source[cell] =  linear_extrapolation_for_enthalpy_defect(infoMatrix[cell][4], sourceA, sourceB, xA, xB);
			}
		}
	}

	else if (iChiPDF == 1)
	{
		for(int cell=1;cell<=numberOfCells;cell++)	
		{
			index = int(infoMatrix[cell][4]/defectStep);
			
			if (infoMatrix[cell][4]<=0.)		// Heat Gain
				Source[cell] =  field[0].find_interpolated_value_global_lognormal
								(infoMatrix[cell][1], infoMatrix[cell][2], infoMatrix[cell][3]); 
		
			else if (index<(nEnthalpyDefects-1))	// linear interpolation
			{
				sourceA = field[index].find_interpolated_value_global_lognormal
						 (infoMatrix[cell][1], infoMatrix[cell][2], infoMatrix[cell][3]); 
				sourceB = field[index+1].find_interpolated_value_global_lognormal
						 (infoMatrix[cell][1], infoMatrix[cell][2], infoMatrix[cell][3]);

				double xA = index*defectStep;
				double xB = (index+1)*defectStep;
				Source[cell] =  linear_interpolation_for_enthalpy_defect(infoMatrix[cell][4], sourceA, sourceB, xA, xB);
			}
						
			else if (index>=(nEnthalpyDefects-1))	// Heat Loss
			{
				sourceA = field[nEnthalpyDefects-2].find_interpolated_value_global_lognormal
						  (infoMatrix[cell][1], infoMatrix[cell][2], infoMatrix[cell][3]); 
				sourceB = field[nEnthalpyDefects-1].find_interpolated_value_global_lognormal
						  (infoMatrix[cell][1], infoMatrix[cell][2], infoMatrix[cell][3]);

				double xA = (nEnthalpyDefects-2)*defectStep;
				double xB = (nEnthalpyDefects-1)*defectStep;
				Source[cell] =  linear_extrapolation_for_enthalpy_defect(infoMatrix[cell][4], sourceA, sourceB, xA, xB);
			}
		}
	}
}

	
void OpenSMOKE_LookUp_Table_Executables::CalculateSourceTerms(	BzzVector &Source, sourceField_flamelet_library &field, 
																BzzMatrix &infoMatrix, int numberOfCells)
{
	if (iChiPDF == 0)		
		for(int cell=1;cell<=numberOfCells;cell++)	
			Source[cell] = field.find_interpolated_value_global
							(infoMatrix[cell][1], infoMatrix[cell][2], infoMatrix[cell][3]);
			
	else if (iChiPDF == 1)		
		for(int cell=1;cell<=numberOfCells;cell++)	
			Source[cell] = field.find_interpolated_value_global_lognormal
							(infoMatrix[cell][1], infoMatrix[cell][2], infoMatrix[cell][3]);
}


vector<string> OpenSMOKE_LookUp_Table_Executables::construct_name_list(const string root, const string phenomenon, const int iModel)
{
	char number[4];
	BzzVectorInt enthalpy_list(5, 0, 50, 100, 150, 200);
	int number_of_enthalpy_defect = enthalpy_list.Size();
	
	vector<string> list_of_names;
	list_of_names.resize(number_of_enthalpy_defect);


	for (int i=0;i<number_of_enthalpy_defect;i++)
	{
		list_of_names[i] =  root;
			my_itoa(enthalpy_list[i+1], number, 10);
		
		list_of_names[i] += number;
		list_of_names[i] += "/";
		list_of_names[i] += phenomenon;
		list_of_names[i] += "/";
		list_of_names[i] += phenomenon;
		list_of_names[i] += "_";
			my_itoa(iModel, number, 10);
		list_of_names[i] += number;
		list_of_names[i] += ".source";
	}

	return list_of_names;
}

double OpenSMOKE_LookUp_Table_Executables::linear_interpolation_for_enthalpy_defect(double &x, double &sourceA, double &sourceB, double &xA, double &xB)
{
	return sourceA + (x-xA)*(sourceB-sourceA)/(xB-xA);
}

double OpenSMOKE_LookUp_Table_Executables::linear_extrapolation_for_enthalpy_defect(double &x, double &sourceA, double &sourceB, double &xA, double &xB)
{
	double source = sourceB + (sourceB-sourceA)/(xB-xA)*(x-xB);
	if (source>=0.) return source;
	else ErrorMessage("The extrapolation could be wrong!");
}
