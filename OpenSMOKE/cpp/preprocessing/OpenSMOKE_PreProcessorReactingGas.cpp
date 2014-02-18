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

#include "preprocessing/OpenSMOKE_PreProcessorReactingGas.h"
 
int	OpenSMOKE_PreProcessorReactingGas::NumberOfReactions()
{
	return NR;
}

void OpenSMOKE_PreProcessorReactingGas::PreProcessing(const string pathName, const string fileNameTransport, const string fileNameElements, const string option)
{
	string fileNameStoichiometry	= pathName + "/stoichiometry.bzz";
	string fileNameThermodynamics	= pathName + "/termo.bzz";
	string fileNameReactions		= pathName + "/reactions.bzz";
	string fileNameNames			= pathName + "/names.bzz";
	
	kinetics.reactionRates = this;

	// 2. Lettura del numero di componenti e del numero di reazioni	
	BzzLoad inputFile;
	inputFile(fileNameStoichiometry);
	inputFile >> NC >> NR;
	inputFile.End();

	// 3. Lettura dello schema cinetico e dei parametri cinetici
	kinetics.readFromFile(fileNameStoichiometry, fileNameReactions);

	// 4. Lettura delle informazioni termodinamiche
	BzzSave binaryFile;
	BzzSave asciiFile;
	binaryFile('*', "idealgas.bin");
	asciiFile("idealgas.ascii");
	Setup(fileNameNames, fileNameTransport, fileNameThermodynamics, fileNameElements, binaryFile, asciiFile);
	kinetics.CheckingStoichiometry();
	Fitting(50);
	WriteFitting(binaryFile,asciiFile);
	binaryFile.End();
	asciiFile.End();

	// 4.1 Sparsity
	kinetics.SparsityStructures(fileNameStoichiometry);
}

void OpenSMOKE_PreProcessorReactingGas::PreProcessingTransportSensitivity(const string pathName, const string fileNameTransport, const string fileNameElements, const double eps)
{
	string fileNameStoichiometry	= pathName + "/ascii/stoichiometry.bzz";
	string fileNameThermodynamics	= pathName + "/ascii/termo.bzz";
	string fileNameReactions		= pathName + "/ascii/reactions.bzz";
	string fileNameNames			= pathName + "/ascii/names.bzz";
	
	kinetics.reactionRates = this;

	// 2. Lettura del numero di componenti e del numero di reazioni	
	BzzLoad inputFile;
	inputFile(fileNameStoichiometry);
	inputFile >> NC >> NR;
	inputFile.End();

	// 3. Lettura dello schema cinetico e dei parametri cinetici
	//kinetics.readFromFile(fileNameStoichiometry, fileNameReactions);

	// 4. Lettura delle informazioni termodinamiche
	BzzSave binaryFile;
	BzzSave asciiFile;
	binaryFile('*', pathName +"/idealgas.dummy");
	asciiFile(pathName +"/idealgas.dummy.ascii");
	SetupTransportSensitivity(pathName, fileNameNames, fileNameTransport, fileNameThermodynamics, fileNameElements, binaryFile, asciiFile, eps);
}

void OpenSMOKE_LennardJonesSensitivityCoefficientsManager::ErrorMessage(const string message)
{
    cout << endl;
    cout << "Class:  OpenSMOKE_LennardJonesSensitivityCoefficientsManager"		<< endl;
    cout << "Object: " << name_object		<< endl;
    cout << "Error:  " << message           << endl;
    cout << "Press a key to continue... "   << endl;
    getchar();
    exit(-1);
}

void OpenSMOKE_LennardJonesSensitivityCoefficientsManager::WarningMessage(const string message)
{
    cout << endl;
    cout << "Class:  OpenSMOKE_LennardJonesSensitivityCoefficientsManager"		<< endl;
    cout << "Object: "		<< name_object	<< endl;
    cout << "Warning:  "	<< message      << endl;
    cout << "Press a key to continue... "   << endl;
    getchar();
} 

OpenSMOKE_LennardJonesSensitivityCoefficientsManager::OpenSMOKE_LennardJonesSensitivityCoefficientsManager()
{
	name_object = "[not assigned]";
	fittingMode = OPENSMOKE_FITTING_ALL; 
}

void OpenSMOKE_LennardJonesSensitivityCoefficientsManager::SetName(const string name)
{
	name_object = name;
}


void OpenSMOKE_LennardJonesSensitivityCoefficientsManager::Setup(const string file_name, fittingCoefficientSensitivityMode _fittingMode)
{
	string dummy;

	fittingMode = _fittingMode;

	BzzLoad fInput(file_name);

	fInput >> dummy;
	if (dummy!="NC")
		ErrorMessage("Expected NC - Found: " + dummy);
	fInput >> NC;

	fInput >> dummy;
	if (dummy!="EPSILON")
		ErrorMessage("Expected EPSILON - Found: " + dummy);
	fInput >> epsilon;

	names.resize(NC+1);
	for(int i=1;i<=NC;i++)
		fInput >> names[i];

	fInput >> epsylon_over_kb;
	fInput >> sigma;
	fInput >> mu;
	fInput >> alfa;
	fInput >> zRot298;

	fInput >> dummy;
	if (dummy!="BASE")
		ErrorMessage("Expected BASE - Found: " + dummy);
	fInput >> fittingEta_Base;
	fInput >> fittingLambda_Base;
	fInput >> fittingDbinary_Base;
	fInput >> fittingTetaBinary_Base;

	cout << "Reading Epsilon/K..." << endl;

	fInput >> dummy;
	if (dummy!="EPSKPLUS")
		ErrorMessage("Expected EPSKPLUS - Found: " + dummy);
	fInput >> fittingEta_KE_Plus;
	fInput >> fittingLambda_KE_Plus;
	fInput >> fittingDbinary_KE_Plus;
	fInput >> fittingTetaBinary_KE_Plus;

	fInput >> dummy;
	if (dummy!="EPSKMINUS")
		ErrorMessage("Expected EPSKMINUS - Found: " + dummy);
	fInput >> fittingEta_KE_Minus;
	fInput >> fittingLambda_KE_Minus;
	fInput >> fittingDbinary_KE_Minus;
	fInput >> fittingTetaBinary_KE_Minus;

	cout << "Reading Sigma..." << endl;

	fInput >> dummy;
	if (dummy!="SIGMAPLUS")
		ErrorMessage("Expected SIGMAPLUS - Found: " + dummy);
	fInput >> fittingEta_Sigma_Plus;
	fInput >> fittingLambda_Sigma_Plus;
	fInput >> fittingDbinary_Sigma_Plus;
	fInput >> fittingTetaBinary_Sigma_Plus;

	fInput >> dummy;
	if (dummy!="SIGMAMINUS")
		ErrorMessage("Expected SIGMAMINUS - Found: " + dummy);
	fInput >> fittingEta_Sigma_Minus;
	fInput >> fittingLambda_Sigma_Minus;
	fInput >> fittingDbinary_Sigma_Minus;
	fInput >> fittingTetaBinary_Sigma_Minus;

	cout << "Reading Mu..." << endl;

	fInput >> dummy;
	if (dummy!="MUPLUS")
		ErrorMessage("Expected MUPLUS - Found: " + dummy);
	fInput >> fittingEta_Mu_Plus;
	fInput >> fittingLambda_Mu_Plus;
	fInput >> fittingDbinary_Mu_Plus;
	fInput >> fittingTetaBinary_Mu_Plus;

	fInput >> dummy;
	if (dummy!="MUMINUS")
		ErrorMessage("Expected MUMINUS - Found: " + dummy);
	fInput >> fittingEta_Mu_Minus;
	fInput >> fittingLambda_Mu_Minus;
	fInput >> fittingDbinary_Mu_Minus;
	fInput >> fittingTetaBinary_Mu_Minus;


	cout << "Reading Alfa..." << endl;

	fInput >> dummy;
	if (dummy!="ALFAPLUS")
		ErrorMessage("Expected ALFAPLUS - Found: " + dummy);
	fInput >> fittingEta_Alfa_Plus;
	fInput >> fittingLambda_Alfa_Plus;
	fInput >> fittingDbinary_Alfa_Plus;
	fInput >> fittingTetaBinary_Alfa_Plus;

	fInput >> dummy;
	if (dummy!="ALFAMINUS")
		ErrorMessage("Expected ALFAMINUS - Found: " + dummy);
	fInput >> fittingEta_Alfa_Minus;
	fInput >> fittingLambda_Alfa_Minus;
	fInput >> fittingDbinary_Alfa_Minus;
	fInput >> fittingTetaBinary_Alfa_Minus;

	cout << "Reading zRot298..." << endl;

	fInput >> dummy;
	if (dummy!="ZROTPLUS")
		ErrorMessage("Expected ZROTPLUS - Found: " + dummy);
	fInput >> fittingEta_zRot298_Plus;
	fInput >> fittingLambda_zRot298_Plus;
	fInput >> fittingDbinary_zRot298_Plus;
	fInput >> fittingTetaBinary_zRot298_Plus;

	fInput >> dummy;
	if (dummy!="ZROTMINUS")
		ErrorMessage("Expected ZROTMINUS - Found: " + dummy);
	fInput >> fittingEta_zRot298_Minus;
	fInput >> fittingLambda_zRot298_Minus;
	fInput >> fittingDbinary_zRot298_Minus;
	fInput >> fittingTetaBinary_zRot298_Minus;

	fInput.End();

	ChangeDimensions(NC, &delta_epsylon_over_kb);
	ChangeDimensions(NC, &delta_sigma);
	ChangeDimensions(NC, &delta_mu);
	ChangeDimensions(NC, &delta_alfa);
	ChangeDimensions(NC, &delta_zRot298);
	const double _2epsilon = 2.*epsilon;
	for(int i=1;i<=NC;i++)
	{
		if (epsylon_over_kb[i]>0.)	delta_epsylon_over_kb[i] = _2epsilon*epsylon_over_kb[i];
		else						delta_epsylon_over_kb[i] = 1.e16;
	}
	for(int i=1;i<=NC;i++)
	{
		if (sigma[i]>0.)	delta_sigma[i] = _2epsilon*sigma[i];
		else				delta_sigma[i] = 1.e16;
	}
	for(int i=1;i<=NC;i++)
	{
		if (mu[i]>0.)	delta_mu[i] = _2epsilon*mu[i];
		else			delta_mu[i] = 1.e16;
	}
	for(int i=1;i<=NC;i++)
	{
		if (alfa[i]>0.)	delta_alfa[i] = _2epsilon*alfa[i];
		else			delta_alfa[i] = 1.e16;
	}
	for(int i=1;i<=NC;i++)
	{
		if (zRot298[i]>0.)	delta_zRot298[i] = _2epsilon*zRot298[i];
		else				delta_zRot298[i] = 1.e16;
	}
}

void OpenSMOKE_LennardJonesSensitivityCoefficientsManager::GetFittingCoefficients(const fittingCoefficientExtraction fitting, const int j, OpenSMOKE_ReactingGas &mix)
{
	BzzVector aux;
	if (fitting == OPENSMOKE_FITTING_BASE)
	{
		#if LINUX_SO==1
			aux = fittingEta_Base.GetRow(j);
			mix.fittingEta.SetRow(j, aux);
			aux = fittingLambda_Base.GetRow(j);
			mix.fittingLambda.SetRow(j, aux);

			unsigned int index = (j-1)*NC;
			for (unsigned int k=index+1;k<=index+NC;k++)
			{
				aux = fittingDbinary_Base.GetRow(k);
				mix.fittingDbinary.SetRow(k, aux);
			}
		#else
			mix.fittingEta.SetRow(j, fittingEta_Base.GetRow(j));
			mix.fittingLambda.SetRow(j, fittingLambda_Base.GetRow(j));

			unsigned int index = (j-1)*NC;
			for (unsigned int k=index+1;k<=index+NC;k++)
				mix.fittingDbinary.SetRow(k, fittingDbinary_Base.GetRow(k));
		#endif
	}

	else if (fitting == OPENSMOKE_FITTING_KE_PLUS)
	{
		#if LINUX_SO==1
			aux = fittingEta_KE_Plus.GetRow(j);
			mix.fittingEta.SetRow(j, aux);
			aux = fittingLambda_KE_Plus.GetRow(j);
			mix.fittingLambda.SetRow(j, aux);

			unsigned int index = (j-1)*NC;
			for (unsigned int k=index+1;k<=index+NC;k++)
			{
				aux = fittingDbinary_KE_Plus.GetRow(k);
				mix.fittingDbinary.SetRow(k, aux);
			}
		#else
			mix.fittingEta.SetRow(j, fittingEta_KE_Plus.GetRow(j));
			mix.fittingLambda.SetRow(j, fittingLambda_KE_Plus.GetRow(j));

			unsigned int index = (j-1)*NC;
			for (unsigned int k=index+1;k<=index+NC;k++)
				mix.fittingDbinary.SetRow(k, fittingDbinary_KE_Plus.GetRow(k));
		#endif
	}

	else if (fitting == OPENSMOKE_FITTING_KE_MINUS)
	{
		#if LINUX_SO==1
			aux = fittingEta_KE_Minus.GetRow(j);
			mix.fittingEta.SetRow(j, aux);
			aux = fittingLambda_KE_Minus.GetRow(j);
			mix.fittingLambda.SetRow(j, aux);

			unsigned int index = (j-1)*NC;
			for (unsigned int k=index+1;k<=index+NC;k++)
			{
				aux = fittingDbinary_KE_Minus.GetRow(k);
				mix.fittingDbinary.SetRow(k, aux);
			}
		#else
			mix.fittingEta.SetRow(j, fittingEta_KE_Minus.GetRow(j));
			mix.fittingLambda.SetRow(j, fittingLambda_KE_Minus.GetRow(j));

			unsigned int index = (j-1)*NC;
			for (unsigned int k=index+1;k<=index+NC;k++)
				mix.fittingDbinary.SetRow(k, fittingDbinary_KE_Minus.GetRow(k));
		#endif
	}

	else if (fitting == OPENSMOKE_FITTING_SIGMA_PLUS)
	{
		#if LINUX_SO==1
			aux = fittingEta_Sigma_Plus.GetRow(j);
			mix.fittingEta.SetRow(j, aux);
			aux = fittingLambda_Sigma_Plus.GetRow(j);
			mix.fittingLambda.SetRow(j, aux);

			unsigned int index = (j-1)*NC;
			for (unsigned int k=index+1;k<=index+NC;k++)
			{
				aux = fittingDbinary_Sigma_Plus.GetRow(k);
				mix.fittingDbinary.SetRow(k, aux);
			}
		#else
			mix.fittingEta.SetRow(j, fittingEta_Sigma_Plus.GetRow(j));
			mix.fittingLambda.SetRow(j, fittingLambda_Sigma_Plus.GetRow(j));

			unsigned int index = (j-1)*NC;
			for (unsigned int k=index+1;k<=index+NC;k++)
				mix.fittingDbinary.SetRow(k, fittingDbinary_Sigma_Plus.GetRow(k));
		#endif
	}

	else if (fitting == OPENSMOKE_FITTING_SIGMA_MINUS)
	{
		#if LINUX_SO==1
			aux = fittingEta_Sigma_Minus.GetRow(j);
			mix.fittingEta.SetRow(j, aux);
			aux = fittingLambda_Sigma_Minus.GetRow(j);
			mix.fittingLambda.SetRow(j, aux);

			unsigned int index = (j-1)*NC;
			for (unsigned int k=index+1;k<=index+NC;k++)
			{
				aux = fittingDbinary_Sigma_Minus.GetRow(k);
				mix.fittingDbinary.SetRow(k, aux);
			}
		#else
			mix.fittingEta.SetRow(j, fittingEta_Sigma_Minus.GetRow(j));
			mix.fittingLambda.SetRow(j, fittingLambda_Sigma_Minus.GetRow(j));

			unsigned int index = (j-1)*NC;
			for (unsigned int k=index+1;k<=index+NC;k++)
				mix.fittingDbinary.SetRow(k, fittingDbinary_Sigma_Minus.GetRow(k));
		#endif
	}

	else if (fitting == OPENSMOKE_FITTING_MU_PLUS)
	{
		#if LINUX_SO==1
			aux = fittingEta_Mu_Plus.GetRow(j);
			mix.fittingEta.SetRow(j, aux);
			aux = fittingLambda_Mu_Plus.GetRow(j);
			mix.fittingLambda.SetRow(j, aux);
			
			unsigned int index = (j-1)*NC;
			for (unsigned int k=index+1;k<=index+NC;k++)
			{
				aux = fittingDbinary_Mu_Plus.GetRow(k);
				mix.fittingDbinary.SetRow(k, aux);
			}
		#else
			mix.fittingEta.SetRow(j, fittingEta_Mu_Plus.GetRow(j));
			mix.fittingLambda.SetRow(j, fittingLambda_Mu_Plus.GetRow(j));

			unsigned int index = (j-1)*NC;
			for (unsigned int k=index+1;k<=index+NC;k++)
				mix.fittingDbinary.SetRow(k, fittingDbinary_Mu_Plus.GetRow(k));
		#endif
	}

	else if (fitting == OPENSMOKE_FITTING_MU_MINUS)
	{
		#if LINUX_SO==1
			aux = fittingEta_Mu_Minus.GetRow(j);
			mix.fittingEta.SetRow(j, aux);
			aux = fittingLambda_Mu_Minus.GetRow(j);
			mix.fittingLambda.SetRow(j, aux);
			
			unsigned int index = (j-1)*NC;
			for (unsigned int k=index+1;k<=index+NC;k++)
			{			
				aux = fittingDbinary_Mu_Minus.GetRow(k);
				mix.fittingDbinary.SetRow(k, aux);
			}	
		#else
			mix.fittingEta.SetRow(j, fittingEta_Mu_Minus.GetRow(j));
			mix.fittingLambda.SetRow(j, fittingLambda_Mu_Minus.GetRow(j));

			unsigned int index = (j-1)*NC;
			for (unsigned int k=index+1;k<=index+NC;k++)
				mix.fittingDbinary.SetRow(k, fittingDbinary_Mu_Minus.GetRow(k));
		#endif
	}

	else if (fitting == OPENSMOKE_FITTING_ALFA_PLUS)
	{
		#if LINUX_SO==1
			aux = fittingEta_Alfa_Plus.GetRow(j);
			mix.fittingEta.SetRow(j, aux);
			aux = fittingLambda_Alfa_Plus.GetRow(j);
			mix.fittingLambda.SetRow(j, aux);
			
			unsigned int index = (j-1)*NC;
			for (unsigned int k=index+1;k<=index+NC;k++)
			{			
				aux = fittingDbinary_Alfa_Plus.GetRow(k);
				mix.fittingDbinary.SetRow(k, aux);
			}	
		#else
			mix.fittingEta.SetRow(j, fittingEta_Alfa_Plus.GetRow(j));
			mix.fittingLambda.SetRow(j, fittingLambda_Alfa_Plus.GetRow(j));

			unsigned int index = (j-1)*NC;
			for (unsigned int k=index+1;k<=index+NC;k++)
				mix.fittingDbinary.SetRow(k, fittingDbinary_Alfa_Plus.GetRow(k));
		#endif
	}

	else if (fitting == OPENSMOKE_FITTING_ALFA_MINUS)
	{
		#if LINUX_SO==1
			aux = fittingEta_Alfa_Minus.GetRow(j);
			mix.fittingEta.SetRow(j, aux);
			aux = fittingLambda_Alfa_Minus.GetRow(j);
			mix.fittingLambda.SetRow(j, aux);
			
			unsigned int index = (j-1)*NC;
			for (unsigned int k=index+1;k<=index+NC;k++)
			{			
				aux = fittingDbinary_Alfa_Minus.GetRow(k);
				mix.fittingDbinary.SetRow(k, aux);
			}	
		#else
			mix.fittingEta.SetRow(j, fittingEta_Alfa_Minus.GetRow(j));
			mix.fittingLambda.SetRow(j, fittingLambda_Alfa_Minus.GetRow(j));

			unsigned int index = (j-1)*NC;
			for (unsigned int k=index+1;k<=index+NC;k++)
				mix.fittingDbinary.SetRow(k, fittingDbinary_Alfa_Minus.GetRow(k));
		#endif
	}

	else if (fitting == OPENSMOKE_FITTING_ZROT_PLUS)
	{
		#if LINUX_SO==1
			aux = fittingEta_zRot298_Plus.GetRow(j);
			mix.fittingEta.SetRow(j, aux);
			aux = fittingLambda_zRot298_Plus.GetRow(j);
			mix.fittingLambda.SetRow(j, aux);
			
			unsigned int index = (j-1)*NC;
			for (unsigned int k=index+1;k<=index+NC;k++)
			{			
				aux = fittingDbinary_zRot298_Plus.GetRow(k);
				mix.fittingDbinary.SetRow(k, aux);
			}	
		#else
			mix.fittingEta.SetRow(j, fittingEta_zRot298_Plus.GetRow(j));
			mix.fittingLambda.SetRow(j, fittingLambda_zRot298_Plus.GetRow(j));

			unsigned int index = (j-1)*NC;
			for (unsigned int k=index+1;k<=index+NC;k++)
				mix.fittingDbinary.SetRow(k, fittingDbinary_zRot298_Plus.GetRow(k));
		#endif
	}

	else if (fitting == OPENSMOKE_FITTING_ZROT_MINUS)
	{
		#if LINUX_SO==1
			aux = fittingEta_zRot298_Minus.GetRow(j);
			mix.fittingEta.SetRow(j, aux);
			aux = fittingLambda_zRot298_Minus.GetRow(j);
			mix.fittingLambda.SetRow(j, aux);
			
			unsigned int index = (j-1)*NC;
			for (unsigned int k=index+1;k<=index+NC;k++)
			{			
				aux = fittingDbinary_zRot298_Minus.GetRow(k);
				mix.fittingDbinary.SetRow(k, aux);
			}
		#else
			mix.fittingEta.SetRow(j, fittingEta_zRot298_Minus.GetRow(j));
			mix.fittingLambda.SetRow(j, fittingLambda_zRot298_Minus.GetRow(j));

			unsigned int index = (j-1)*NC;
			for (unsigned int k=index+1;k<=index+NC;k++)
				mix.fittingDbinary.SetRow(k, fittingDbinary_zRot298_Minus.GetRow(k));
		#endif
	}

	else	ErrorMessage("No fitting option available...");



	if ( (fittingMode != OPENSMOKE_FITTING_ALL) && (fitting != OPENSMOKE_FITTING_BASE) )
	{
		#if LINUX_SO==1
			if (fittingMode == OPENSMOKE_FITTING_ONLY_VISCOSITY)
			{
				aux = fittingLambda_Base.GetRow(j);
				mix.fittingLambda.SetRow(j, aux);
				unsigned int index = (j-1)*NC;
				for (unsigned int k=index+1;k<=index+NC;k++)
				{			
					aux = fittingDbinary_Base.GetRow(k);
					mix.fittingDbinary.SetRow(k, aux);
				}
			}
			else if (fittingMode == OPENSMOKE_FITTING_ONLY_CONDUCTIVITY)
			{
				aux = fittingEta_Base.GetRow(j);
				mix.fittingEta.SetRow(j, aux);
				unsigned int index = (j-1)*NC;
				for (unsigned int k=index+1;k<=index+NC;k++)
				{			
					aux = fittingDbinary_Base.GetRow(k);
					mix.fittingDbinary.SetRow(k, aux);
				}
			}
			else if (fittingMode == OPENSMOKE_FITTING_ONLY_DIFFUSIVITIES)
			{
				aux = fittingEta_Base.GetRow(j);
				mix.fittingEta.SetRow(j, aux);
				aux = fittingLambda_Base.GetRow(j);
				mix.fittingLambda.SetRow(j, aux);
			}	
		#else
			if (fittingMode == OPENSMOKE_FITTING_ONLY_VISCOSITY)
			{
				mix.fittingLambda.SetRow(j, fittingLambda_Base.GetRow(j));
				unsigned int index = (j-1)*NC;
				for (unsigned int k=index+1;k<=index+NC;k++)
					mix.fittingDbinary.SetRow(k, fittingDbinary_Base.GetRow(k));
			}
			else if (fittingMode == OPENSMOKE_FITTING_ONLY_CONDUCTIVITY)
			{
				mix.fittingEta.SetRow(j, fittingEta_Base.GetRow(j));
				unsigned int index = (j-1)*NC;
				for (unsigned int k=index+1;k<=index+NC;k++)
					mix.fittingDbinary.SetRow(k, fittingDbinary_Base.GetRow(k));
			}
			else if (fittingMode == OPENSMOKE_FITTING_ONLY_DIFFUSIVITIES)
			{
				mix.fittingEta.SetRow(j, fittingEta_Base.GetRow(j));
				mix.fittingLambda.SetRow(j, fittingLambda_Base.GetRow(j));
			}			
		#endif
	}
}
