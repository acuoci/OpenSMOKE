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

#include <vector>
#include "basic/OpenSMOKE_Utilities.h"
#include "engine/OpenSMOKE_ReactingGas.h"
#include "addons/OpenSMOKE_SensitivityAnalysis.h"


OpenSMOKE_SensitivityAnalysis::OpenSMOKE_SensitivityAnalysis()
{
}

void OpenSMOKE_SensitivityAnalysis::Initialize(OpenSMOKE_ReactingGas *_mix, BzzVectorInt &_indices, bool _iTemperature)
{
	int i,k;

	mix				= _mix;
	iTemperature	= _iTemperature;

	NC = mix->NumberOfSpecies();			// Number of Species
	NR = mix->NumberOfReactions();			// Number of reactions

	NP = NR;								// Number of parameters
	
	if (iTemperature == true)
	{
		NV = NC+1;									// Number of variables
		indexSpecies     = 1;						// index of species
		indexTemperature = NC+1;					// index of temperature
	}
	
	if (iTemperature == false)
	{
		NV = NC;					// Number of variables
		indexSpecies     = 1;		// index of species
		indexTemperature = 0;		// index of temperature
	}

	ChangeDimensions(NV, NP, &JAlfa);
	ChangeDimensions(NV, NP, &S);
	ChangeDimensions(NV, NP, &Somega);
	ChangeDimensions(NV, NP, &Sx);
	ChangeDimensions(NV, NV, &Jacobian);
	ChangeDimensions(NV, &I);				// Identity matrix (Diagonal)
	ChangeDimensions(NV, NP, &Sold);		// Old sensitivity coefficient matrix
	ChangeDimensions(NV,	 &scaling);		// Scaling factors


	M  = mix->M;							// Molecular weigths [kg/kmol]
	uM = mix->uM;							// uMolecular weigths [kmol/kg]
	A  = mix->kinetics.exp_k01;				// Frequency factors [kmol, m3, s]

	nu = new OpenSMOKE_NuManager[NC+1];

	for(i =1;i<=NC;i++)	
		nu[i].Set(i, mix->names[i]);

	BuildNuMatrix(&mix->kinetics);
	
	for(k=1;k<=NC;k++)
		nu[k].Clean();
	
	indices_print_species = _indices;

	threshold_normalization = 1e-10;
	iImplicit				= true;
	I						= 1.0;
}

void OpenSMOKE_SensitivityAnalysis::BuildNuMatrix(OpenSMOKE_Kinetics *kinetics)
{
	int i,k;

	int*	jD1 = kinetics->jDir1.GetHandle();
	int*	jD2 = kinetics->jDir2.GetHandle();
	int*	jD3 = kinetics->jDir3.GetHandle();
	int*	jD4 = kinetics->jDir4.GetHandle();
	int*	jD5 = kinetics->jDir5.GetHandle();
	double* vD5 = kinetics->valDir5.GetHandle();

	int*	jIT1 = kinetics->jInvTot1.GetHandle();
	int*	jIT2 = kinetics->jInvTot2.GetHandle();
	int*	jIT3 = kinetics->jInvTot3.GetHandle();
	int*	jIT4 = kinetics->jInvTot4.GetHandle();
	int*	jIT5 = kinetics->jInvTot5.GetHandle();
	double* vIT5 = kinetics->valInvTot5.GetHandle();
	
	for(i = 1;i <= NC;i++)
	{
		for(k = 1;k <= kinetics->numDir1[i];k++)
			nu[i].Set(-1., *jD1++);
		for(k = 1;k <= kinetics->numDir2[i];k++)
			nu[i].Set(-2., *jD2++);
		for(k = 1;k <= kinetics->numDir3[i];k++)
			nu[i].Set(-4., *jD3++);
		for(k = 1;k <= kinetics->numDir4[i];k++)
			nu[i].Set(-0.50, *jD4++);
		for(k = 1;k <= kinetics->numDir5[i];k++)
			nu[i].Set(-(*vD5++), *jD5++);

		for(k = 1;k <= kinetics->numInvTot1[i];k++)
			nu[i].Set(1., *jIT1++);
		for(k = 1;k <= kinetics->numInvTot2[i];k++)
			nu[i].Set(2., *jIT2++);
		for(k = 1;k <= kinetics->numInvTot3[i];k++)
			nu[i].Set(3., *jIT3++);
		for(k = 1;k <= kinetics->numInvTot4[i];k++)
			nu[i].Set(0.50, *jIT4++);
		for(k = 1;k <= kinetics->numInvTot5[i];k++)
			nu[i].Set(*vIT5++, *jIT5++);
	}
}

void OpenSMOKE_SensitivityAnalysis::BuildJAlfaMatrix(BzzVector &r, const double T)
{
/*	int i,j;

	JAlfa = 0.;

	mix->GiveMe_Jalfa_Species(JAlfa, nu, indexSpecies);
	if (iTemperature == true)
		mix->GiveMe_Jalfa_Temperature(JAlfa, T, indexTemperature);

	for(i=1;i<=NV;i++)
		for(j=1;j<=NP;j++)
			JAlfa[i][j] /= scaling[i]; */
}


void OpenSMOKE_SensitivityAnalysis::GiveMe_SensitivityCoefficients()
{
	BzzFactorizedGauss AFactorized(Jacobian);
	Solve(&AFactorized, -JAlfa, &S);
}

void OpenSMOKE_SensitivityAnalysis::GiveMe_SensitivityCoefficients(const double deltat)
{
	// Implicit Method
	if (iImplicit == true)
	{
		JAlfa *= deltat;
		Sum(JAlfa, Sold, &Sold);
		Jacobian *= -deltat;
		Sum(I,Jacobian, &Jacobian);

		BzzFactorizedGauss AFactorized(Jacobian);
		Solve(&AFactorized, Sold, &S);
	}
/*	else
	// Explicit Method
	{
		Product(JacobianOld, Sold, &S);
		Sum(JacobianOld, JalfaOld, &JacobianOld);

		JacobianOld	*= deltat;
				
		Sum(Sold, JacobianOld, &S);
	}*/

	// Update old value
	Sold = S;
}

void OpenSMOKE_SensitivityAnalysis::Normalize_SensitivityCoefficients(BzzVector &omega, BzzVector &x, const double MWmix, const double T)
{
	int k;

	Somega	= 0.;
	Sx		= 0.;

	for(k=1;k<=NC;k++)
	{
		if (omega[k] > threshold_normalization)
			for(int i=1;i<=NP;i++)
				Somega[k][i] = S[k][i] * A[i]/omega[k];
	}

	
	for(k=1;k<=NC;k++)
	{
		if (omega[k] > threshold_normalization)
			for(int i=1;i<=NP;i++)
			{
				double sum = 0.;
				for(int j=1;j<=NC;j++)
					sum += uM[j]*S[j][i];
				sum *= A[i]*MWmix;
		
				Sx[k][i] = Somega[k][i] - sum;
			}
	}

	if (iTemperature == true)
		for(int i=1;i<=NP;i++)
			Sx[indexTemperature][i] = A[i] / T * S[indexTemperature][i];
}

void OpenSMOKE_SensitivityAnalysis::PrintOnFile_SensitivityCoefficients(ofstream &fOutput, int count)
{
	for(int i=1;i<=indices_print_species.Size();i++)
	{
		int k = indices_print_species[i];
		for (int jj=1;jj<=NR;jj++)
			fOutput << mix->names[k] << "_R" << jj << "(" << count++ << ")\t";
	}

	if (iTemperature == true)
		for (int jj=1;jj<=NR;jj++)
			fOutput << "T_R" << jj << "(" << count++ << ")\t";

	fOutput << endl;
	fOutput << endl;
}

void OpenSMOKE_SensitivityAnalysis::PrintOnFile_SensitivityCoefficients(ofstream &fOutput)
{
	for(int i=1;i<=indices_print_species.Size();i++)
	{
		int k = indices_print_species[i];
		for (int jj=1;jj<=NR;jj++)
			fOutput << Sx[k][jj] << "\t";
	}

	if (iTemperature == true)
		for (int jj=1;jj<=NR;jj++)
			fOutput << Sx[indexTemperature][jj] << "\t";

	fOutput << endl;
}