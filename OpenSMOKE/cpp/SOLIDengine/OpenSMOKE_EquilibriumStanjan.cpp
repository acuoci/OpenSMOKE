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

#include "Linux_Windows.h"
#include <iomanip>

#if LINUX_SO==0
	#include "nr.h"
#endif

#include "basic/OpenSMOKE_Constants.h"
#include "engine/OpenSMOKE_IdealGas.h"
#include "engine/OpenSMOKE_EquilibriumStanjan.h"

OpenSMOKE_EquilibriumStanjan *ptEquilibrium;

const double OpenSMOKE_EquilibriumStanjan::threshold_element	=	1.e-16;

// ---------------------------------------------------------------------------------------------
// Non linear system solver
// ---------------------------------------------------------------------------------------------
void MyNonLinearSystem_EquilibriumStanjan::GetResiduals(BzzVector &lambda, BzzVector &f)
{
	//
}

void MyNonLinearSystem_EquilibriumStanjan::AssignEquilibrium(OpenSMOKE_EquilibriumStanjan *equilibrium)
{
	ptEquilibrium = equilibrium;
}

void MyNonLinearSystem_EquilibriumStanjan::ObjectBzzPrint(void)
{
}

void OpenSMOKE_EquilibriumStanjan::ErrorMessage(const string message)
{
    cout << endl;
    cout << "Class:  OpenSMOKE_EquilibriumStanjan"	<< endl;
    cout << "Object: " << name_object				<< endl;
    cout << "Error:  " << message					<< endl;
    cout << "Press enter to continue... "			<< endl;
    getchar();
    exit(-1);
}

void OpenSMOKE_EquilibriumStanjan::WarningMessageStrong(const string message)
{
    cout << endl;
    cout << "Class:	  OpenSMOKE_EquilibriumStanjan"	<< endl;
    cout << "Object:  " << name_object				<< endl;
    cout << "Warning: " << message					<< endl;
    cout << "Press enter to continue..."			<< endl;
	getchar();
	cout << endl;
}

void OpenSMOKE_EquilibriumStanjan::WarningMessageSoft(const string message)
{
    cout << endl;
    cout << "Class:	  OpenSMOKE_EquilibriumStanjan"	<< endl;
    cout << "Object:  " << name_object				<< endl;
    cout << "Warning: " << message					<< endl;
	cout << endl;
}

OpenSMOKE_EquilibriumStanjan::OpenSMOKE_EquilibriumStanjan()
{
	name_object	= "[Name not assigned]";
	assignedElementalComposition = false;
	iVerboseFile				 = true;

	countMode1		= 0;
	countMode2		= 0;
	countMode2bis	= 0;
	countMode3		= 0;

	NmaxIterationsMode1 = 100;
	NmaxIterationsMode2 = 10;
	NmaxIterations      = 10;

	toleranceMode1		= 1.e-3;
	toleranceMode2		= 1.e-4;
	toleranceSteeping	= 1.e-6;

	iSolidPhase  = false;
	iLiquidPhase = false;
	iMultiPhase  = false;
}

void OpenSMOKE_EquilibriumStanjan::SetName(string name)
{
	name_object = name;
}

void OpenSMOKE_EquilibriumStanjan::Setup(OpenSMOKE_IdealGas *ideal_gas)
{
	int k;

	gas = ideal_gas;

	NC	= gas->NumberOfSpecies();
	NE	= gas->NumberOfElements();
	E	= gas->elements;
	
	ChangeDimensions(NC, &gtilde);
	ChangeDimensions(NC, &indexSpecies_Original);
	ChangeDimensions(NE, &indexElement_Original);
	for(k=1;k<=NC;k++)	
		indexSpecies_Original[k] = k;
	for(k=1;k<=NE;k++)	
		indexElement_Original[k] = k;

	jCarbon		= gas->recognize_species_without_exit("CSOLID");
	jWater		= gas->recognize_species_without_exit("H2OLIQUID");
	
	if (iVerboseFile == true)
	{
		openOutputFileAndControl(fLog, "Equilibrium.log");
		fLog.setf(ios::scientific);
	}
}

void OpenSMOKE_EquilibriumStanjan::SetElementalCompositionFromSpeciesMoleFractions(BzzVector &_xInitial)
{
	// Elemental composition
	Product(E, _xInitial, &p);
	
	// Reduced variables: if some elements are zero
	NER	= NE;	// number of elements
	ER	= E;	// matrix of elements
	pR	= p;	// elemental compositions

	// Additional reduced variables
	indexIncludeElement = indexElement_Original;;
	for(int k=NE;k>=1;k--)
		if (p[k] <= threshold_element)
		{
			NER--;
			ER.DeleteRow(k);
			pR.DeleteElement(k);
			indexIncludeElement.DeleteElement(k);
		}

	// Additional reduced variables
	indexIncludeSpecies = indexSpecies_Original;
	ChangeDimensions(0, &indexNoSpecies);
	if (NE != NER)
	{	for(int i=NC;i>=1;i--)
		{
			for(int k=NE;k>=1;k--)
				if (p[k] <= threshold_element)
					if(E[k][i] != 0.)
					{
						indexNoSpecies.Append(i);
						indexIncludeSpecies.DeleteElement(i);
						break;
					}
		}

		// indexNoSpecies is in the opposite order
		for(i=1;i<=indexNoSpecies.Size();i++)
			ER.DeleteColumn(indexNoSpecies[i]);
	}

	// First guess solution preparation
	NCR = indexIncludeSpecies.Size();
	ChangeDimensions(NER, &lambdaFinalSolutionR);
	ChangeDimensions(NCR, &gtildeR);
	ChangeDimensions(NCR, &jPhase);

	ChangeDimensions(NER, &H);
	ChangeDimensions(NER, &f);
	ChangeDimensions(NER, &deltalambda);
	ChangeDimensions(NER, NER, &Q);

	ChangeDimensions(NCR, &BR);
	ChangeDimensions(NER, &bR);

	// For initializing the problem
	ERT = ER;
	Transpose(&ERT);

	Product(ERT, ER, &ERT_x_ER);
	ERT_x_ER_Factorized = ERT_x_ER;
	Product(ER, ERT, &ER_x_ERT);
	ER_x_ERT_Factorized = ER_x_ERT;

	// 
	ERT_Row = new BzzVector[NCR+1];
	for(int j=1;j<=NCR;j++)
	{
		ChangeDimensions(NER, &ERT_Row[j]);
		ERT_Row[j] = ERT.GetRow(j);
	}

	// Number of phases
	nPhases = 1;

	// Checking for solid phase
	if (jCarbon>0)	
	{
		iSolidPhase = true;
		for(j=1; j<=indexNoSpecies.Size();j++)
			if (indexNoSpecies[j] == jCarbon) 
			{
				iSolidPhase = false;
				break;
			}
		for(j=1; j<=NCR;j++)
			if (gas->names[indexIncludeSpecies[j]] == "CSOLID")
			{
				jCarbonR = j; 
				break;
			}
	}
	iSolidPhase = false;
	if (iSolidPhase == true) nPhases++;

	// Checking for liquid phase
	if (jWater>0)	
	{
		iLiquidPhase = true;
		for(j=1; j<=indexNoSpecies.Size();j++)
			if (indexNoSpecies[j] == jWater) 
			{
				iLiquidPhase = false;
				break;
			}
		for(j=1; j<=NCR;j++)
			if (gas->names[indexIncludeSpecies[j]] == "H2OLIQUID")
			{
				jWaterR = j; 
				break;
			}
	}
	if (iLiquidPhase == true) nPhases++;

	// Phase of each species
	jPhase = 1;
	if (iSolidPhase  == true) jPhase[jCarbon] = 2;
	if (iLiquidPhase == true) jPhase[jWater]  = 3;
	if (iSolidPhase == true || iLiquidPhase == true) iMultiPhase = true;

	// Matrices
	ChangeDimensions(NER, nPhases, &D);
	ChangeDimensions(NER, nPhases, &EE);
	ChangeDimensions(NER+nPhases, NER+nPhases, &QExtended_Factorized);
	ChangeDimensions(NER+nPhases, &bExtended);
	ChangeDimensions(NER+nPhases, &deltaExtended);


	// Assigned elemental composition
	assignedElementalComposition = true;
}

void OpenSMOKE_EquilibriumStanjan::Equilibrate(const double _T, const double _P_Pa, BzzVector &xFinal)
{
	// Get Free Gibbs Energy for the species[-]
	GetGtilde(_T);

	// Get Lambda First Guess
	GetFirstGuess(_T, _P_Pa);
}

void OpenSMOKE_EquilibriumStanjan::GetLambdaFirstGuessFromLeastSquares(BzzVector &xFirstGuessR, BzzVector &lambdaFirstGuessR)
{
	int i;

	// Build System
	for(i=1;i<=NCR;i++) 
		BR[i] = gtildeR[i] + log(xFirstGuessR[i]);

	// Solving Least squares problem
	Product(ER, BR, &bR);
	Solve(ER_x_ERT_Factorized, bR, &lambdaFirstGuessR);

	// Write on File
	if (iVerboseFile == true)
	{
		int j;

		fLog << "--------------------------------------------------------" << endl;
		fLog << "                First Guess (Least Squares)             " << endl;
		fLog << "--------------------------------------------------------" << endl;
		for(j=1;j<=NER;j++)
			fLog << setw(13) <<  " lambda_"	<< gas->list_of_elements[indexIncludeElement[j]-1] << "\t" << lambdaFirstGuessR[j] << endl;
		fLog << endl << endl;

		fLog << "--------------------------------------------------------" << endl;
		fLog << "                First Guess (Least Squares)             " << endl;
		fLog << "--------------------------------------------------------" << endl;		
		for(j=1;j<=NCR;j++)
			fLog << setw(16) << gas->names[indexIncludeSpecies[j]] << "\t" << xFirstGuessR[j] << endl;
		fLog << endl << endl;
	}
}

void OpenSMOKE_EquilibriumStanjan::GetLambdaFirstGuessFromLinearProgramming(double &NtotFirstGuess, BzzVector &xFirstGuessR, BzzVector &lambdaFirstGuessR)
{
	int i;
	BzzVector NFirstGuessMin(NCR);
	BzzVectorInt indexBase(NER);
	BzzMatrix G(NER,NER);
	BzzFactorizedGauss G_Factorized;

	// Composition from linear programming
	LinearProgramming(NFirstGuessMin);
	
	// Composition from linear programming
	NtotFirstGuess = NFirstGuessMin.GetSumElements();
	for(i=1;i<=NCR;i++)
		xFirstGuessR[i] = NFirstGuessMin[i]/NtotFirstGuess;
	
	// Finding base species
	int countBase = 0;
	for(i=1;i<=NCR;i++)
		if (xFirstGuessR[i] > threshold_element)
		{
			countBase++;
			indexBase[countBase] = i;
		}

	// Checking number of base species
	if (countBase == NER)	
	{
		// Assembling linear system
		for(i=1;i<=NER;i++)	G.SetRow(i, ER.GetColumn(indexBase[i]));
		G_Factorized = G;
		for(i=1;i<=NER;i++) lambdaFirstGuessR[i] = log(xFirstGuessR[indexBase[i]])+gtildeR[indexBase[i]];
	
		// Solving linear system
		Solve(G_Factorized, &lambdaFirstGuessR);
	}
	else
		WarningMessageStrong("The number of base species is lower than the number of atomic species!");

	if (iVerboseFile == true)
	{
		int j;

		fLog << "----------------------------------------------------------" << endl;
		fLog << " First Guess: Linear Programming                          " << endl;
		fLog << "----------------------------------------------------------" << endl;
		fLog << " lambda"	<< "\t";	for(j=1;j<=NER;j++)	fLog << lambdaFirstGuessR[j]	<< "\t"; fLog << endl;
		for(j=1;j<=NCR;j++)
			fLog << j << " " << gas->names[indexIncludeSpecies[j]] << "\t" << xFirstGuessR[j] << endl;
		fLog << endl;
	}
}

void OpenSMOKE_EquilibriumStanjan::GetFirstGuess(const double T, const double P_Pa)
{
	int i;

	// First Guess Lambda
	BzzVector NtotFirstGuess(nPhases);
	BzzVector lambdaFirstGuessR(NER);
	BzzVector xFirstGuessR(NCR); 

	// Pope method
	{
		BzzVector NFirstGuessMin(NCR); 
		BzzVector NFirstGuessMinMax(NCR); 
		LinearProgramming(NFirstGuessMin);
		LinearProgrammingMinMax(NFirstGuessMinMax);
		Blending(NtotFirstGuess, xFirstGuessR, NFirstGuessMin, NFirstGuessMinMax);
		GetLambdaFirstGuessFromLeastSquares(xFirstGuessR, lambdaFirstGuessR);
	}

	// Stanjan method
	{
//		GetLambdaFirstGuessFromLinearProgramming(NtotFirstGuess, xFirstGuessR, lambdaFirstGuessR);
	}

	BzzVector	Ntot	= NtotFirstGuess;
	BzzVector	lambdaR	= lambdaFirstGuessR;
	BzzVector	xR		= xFirstGuessR;

	exit_status flag;
	
	Solution(Ntot, lambdaR, xR);
	
	// Solution
	FinalSummary(Ntot, xR, T, P_Pa);
}

void OpenSMOKE_EquilibriumStanjan::GetGtilde(const double T)
{
	gas->GetStandardGibbsFreeEnergy_Mass(gtilde, T);
	gtilde /= (Constants::R_J_kmol*T);
	ElementByElementProduct(gtilde, gas->M, &gtilde);

	for (int i=1;i<=NCR;i++)
		gtildeR[i] = gtilde[indexIncludeSpecies[i]];
}



void OpenSMOKE_EquilibriumStanjan::Mode1(BzzVector &Ntot, BzzVector &lambda, BzzVector &x, exit_status &flag)
{
	countMode1++;
	
	int		k,j;
	double	Beta;
	double	deltas;
	double	Wold;
	double	W;
	BzzVector	Zold(nPhases);
	BzzVector	Z(nPhases);

	// Z and W
	if (nPhases == 1)
	{
		Zold[1] = x.GetSumElements();
		Wold    = Ntot[1]*(Zold[1]-1.) - Dot(lambda, pR);
	}

	if (nPhases == 2)
	{
		Zold[1] = x.GetSumElements();
		Zold[1] -= x[jCarbonR];
		Zold[2]  = x[jCarbonR];
		Wold     = Ntot[1]*(Zold[1]-1.) + Ntot[2]*(Zold[2]-1.) - Dot(lambda, pR);
	}

	// Vector H
	GiveMeH(Ntot,x);

	// Beta
	Beta = -sqrt(Dot(H,H));

	if(Beta < -threshold_element*100.)	// otherwise convergence is verified
	{
		// Matrix Q
		GiveMeQ(Ntot,x);

		// Vector f
		for(k=1;k<=NER;k++)
			f[k] = H[k]/Beta;

		// DeltaS
		double sum=0.;
		for(k=1;k<=NER;k++)
			for(j=1;j<=NER;j++)
				sum += Q[j][k]*f[j]*f[k];
		deltas = -Beta/sum;

		// lambda
		for(j=1;j<=NER;j++)
			lambda[j] += f[j]*deltas;
	}
		
	// x
	for(j=1;j<=NCR;j++)
		x[j] = exp(-gtildeR[j] + Dot( ERT_Row[j],lambda) );

	// Z and W
	if (nPhases == 1)
	{
		Z[1] = x.GetSumElements();
		W    = Ntot[1]*(Z[1]-1.) - Dot(lambda, pR);
	}

	if (nPhases == 2)
	{
		Z[1]  = x.GetSumElements();
		Z[1] -= x[jCarbonR];
		Z[2]  = x[jCarbonR];
		W     = Ntot[1]*(Z[1]-1.) + Ntot[2]*(Z[2]-1.) - Dot(lambda, pR);
	}

	if (iVerboseFile == true && nPhases == 1)
	{
		fLog << "----------------------------------------------------------------------"	<< endl;
		fLog << " Mode 1" << "\t" << "#" << countMode1 << "                            " 	<< endl;
		fLog << "----------------------------------------------------------------------"	<< endl;
		fLog << " W\t" << Wold << "\t" << W    << "\t" << W-Wold << "\t" << (W-Wold)/W*100. << "%" << endl;
		fLog << " Z\t" << Zold[1] << "\t" << Z[1]    << "\t" << Z[1]-Zold[1] << "\t" << (Z[1]-Zold[1])/Z[1]*100. << "%" << endl;
		fLog << " N\t" << Ntot[1] << "\t" << Ntot[1] << "\t" << 0.     << "\t" << 0.              << "%" << endl;
		for(j=1;j<=NER;j++)
			fLog << setw(8) <<  " lambda_"	<< gas->list_of_elements[indexIncludeElement[j]-1] << "\t\t" << lambda[j] << endl;
		fLog << endl << endl;
	}

	if (iVerboseFile == true && nPhases == 2)
	{
		fLog << "----------------------------------------------------------------------"	<< endl;
		fLog << " Mode 1" << "\t" << "#" << countMode1 << "                            " 	<< endl;
		fLog << "----------------------------------------------------------------------"	<< endl;
		fLog << " W\t" << Wold << "\t" << W    << "\t" << W-Wold << "\t" << (W-Wold)/W*100. << "%" << endl;
		fLog << " ZGas\t" << Zold[1] << "\t" << Z[1]    << "\t" << Z[1]-Zold[1] << "\t" << (Z[1]-Zold[1])/Z[1]*100. << "%" << endl;
		fLog << " ZSol\t" << Zold[2] << "\t" << Z[2]    << "\t" << Z[2]-Zold[2] << "\t" << (Z[2]-Zold[2])/Z[2]*100. << "%" << endl;
		fLog << " Ngas\t" << Ntot[1] << "\t" << Ntot[1] << "\t" << 0.     << "\t" << 0.              << "%" << endl;
		fLog << " Nsol\t" << Ntot[2] << "\t" << Ntot[2] << "\t" << 0.     << "\t" << 0.              << "%" << endl;
		fLog << " Ntot\t" << Ntot[1]+Ntot[2] << "\t" << Ntot[1]+Ntot[2] << "\t" << 0.     << "\t" << 0.              << "%" << endl;
		for(j=1;j<=NER;j++)
			fLog << setw(8) <<  " lambda_"	<< gas->list_of_elements[indexIncludeElement[j]-1] << "\t\t" << lambda[j] << endl;
		fLog << endl << endl;
	}

	if ( fabs((W-Wold)/W) < toleranceMode1 )	{ flag =  CONVERGED;	return; }	// Convergence was reached
	if ( W < Wold )								{ flag =  GO_ON;		return; }	
	if ( W > Wold )								{ flag =  FAILED;		return; }	
}


void OpenSMOKE_EquilibriumStanjan::Mode2(BzzVector &Ntot, BzzVector &lambda, BzzVector &x, exit_status &flag)
{
	countMode2++;

	int j;
	double W;
	double Wold;
	BzzVector Z(nPhases);
	BzzVector Zold(nPhases);

	// Z and W
	if (nPhases == 1)
	{
		Zold[1] = x.GetSumElements();
		Wold    = Ntot[1]*(Zold[1]-1.) - Dot(lambda, pR);
	}

	if (nPhases == 2)
	{
		Zold[1] = x.GetSumElements();
		Zold[1] -= x[jCarbonR];
		Zold[2]  = x[jCarbonR];
		Wold     = Ntot[1]*(Zold[1]-1.) + Ntot[2]*(Zold[2]-1.) - Dot(lambda, pR);
	}

	// Vector H
	GiveMeH(Ntot,x);

	// Matrix Q
	GiveMeQ(Ntot,x);

	// Linear System Solution
	Q_Factorized = Q;
	Solve(Q_Factorized, -H, &deltalambda);

	// lambdaR
	for(j=1;j<=NER;j++)
		lambda[j] += deltalambda[j];

	// x
	for(j=1;j<=NCR;j++)
		x[j] = exp( -gtildeR[j] + Dot(ERT_Row[j],lambda) );

	// Z and W
	if (nPhases == 1)
	{
		Z[1] = x.GetSumElements();
		W    = Ntot[1]*(Z[1]-1.) - Dot(lambda, pR);
	}

	if (nPhases == 2)
	{
		Z[1]  = x.GetSumElements();
		Z[1] -= x[jCarbonR];
		Z[2]  = x[jCarbonR];
		W     = Ntot[1]*(Z[1]-1.) + Ntot[2]*(Z[2]-1.) - Dot(lambda, pR);
	}

	if (iVerboseFile == true && nPhases == 1)
	{
		fLog << "----------------------------------------------------------------------"	<< endl;
		fLog << " Mode 2" << "\t" << "#" << countMode2 << "                            " 	<< endl;
		fLog << "----------------------------------------------------------------------"	<< endl;
		fLog << " W\t" << Wold << "\t" << W    << "\t" << W-Wold << "\t" << (W-Wold)/W*100. << "%" << endl;
		fLog << " Z\t" << Zold[1] << "\t" << Z[1]    << "\t" << Z[1]-Zold[1] << "\t" << (Z[1]-Zold[1])/Z[1]*100. << "%" << endl;
		fLog << " N\t" << Ntot[1] << "\t" << Ntot[1] << "\t" << 0.     << "\t" << 0.              << "%" << endl;
		for(j=1;j<=NER;j++)
			fLog << setw(8) <<  " lambda_"	<< gas->list_of_elements[indexIncludeElement[j]-1] << "\t\t" << lambda[j] << endl;
		fLog << endl << endl;
	}

	if (iVerboseFile == true && nPhases == 2)
	{
		fLog << "----------------------------------------------------------------------"	<< endl;
		fLog << " Mode 2" << "\t" << "#" << countMode2 << "                            " 	<< endl;
		fLog << "----------------------------------------------------------------------"	<< endl;
		fLog << " W\t" << Wold << "\t" << W    << "\t" << W-Wold << "\t" << (W-Wold)/W*100. << "%" << endl;
		fLog << " ZGas\t" << Zold[1] << "\t" << Z[1]    << "\t" << Z[1]-Zold[1] << "\t" << (Z[1]-Zold[1])/Z[1]*100. << "%" << endl;
		fLog << " ZSol\t" << Zold[2] << "\t" << Z[2]    << "\t" << Z[2]-Zold[2] << "\t" << (Z[2]-Zold[2])/Z[2]*100. << "%" << endl;
		fLog << " Ngas\t" << Ntot[1] << "\t" << Ntot[1] << "\t" << 0.     << "\t" << 0.              << "%" << endl;
		fLog << " Nsol\t" << Ntot[2] << "\t" << Ntot[2] << "\t" << 0.     << "\t" << 0.              << "%" << endl;
		fLog << " Ntot\t" << Ntot[1]+Ntot[2] << "\t" << Ntot[1]+Ntot[2] << "\t" << 0.     << "\t" << 0.              << "%" << endl;
		for(j=1;j<=NER;j++)
			fLog << setw(8) <<  " lambda_"	<< gas->list_of_elements[indexIncludeElement[j]-1] << "\t\t" << lambda[j] << endl;
		fLog << endl << endl;
	}

	if ( fabs((W-Wold)/W) < toleranceMode2 )	{ flag = CONVERGED; return; }	// Convergence was reached
	if ( W > Wold*1.01 )						{ flag = FAILED;	return; }	
	if ( W < Wold*1.01 )						{ flag = GO_ON;		return; }	
}

void OpenSMOKE_EquilibriumStanjan::Mode2bis(BzzVector &Ntot, BzzVector &lambda, BzzVector &x, exit_status &flag)
{
	countMode2bis++;

	int j;
	double alfa;
	double deltas;
	double W;
	double Wold;
	BzzVector V(nPhases);
	BzzVector r(nPhases);
	BzzMatrix A(nPhases, nPhases);
	BzzVector Z(nPhases);
	BzzVector Zold(nPhases);

	// Z and W
	if (nPhases == 1)
	{
		Zold[1] = x.GetSumElements();
		V[1]    = Zold[1]-1.;
		Wold    = Ntot[1]*(Zold[1]-1.) - Dot(lambda, pR);
	}

	if (nPhases == 2)
	{
		Zold[1]  = x.GetSumElements();
		Zold[1] -= x[jCarbonR];
		V[1]     = Zold[1]-1.;
		Zold[2]  = x[jCarbonR];
		V[2]     = Zold[2]-1.;
		Wold     = Ntot[1]*V[1] + Ntot[2]*V[2] - Dot(lambda, pR);
	}

	// alfa
	alfa = sqrt(Dot(V,V));

	// r
	for(j=1;j<=nPhases;j++)
		r[j] = V[j]/alfa;

	// Vector D
	if (nPhases == 1)
	{
		BzzMatrix mX(NCR,1);
		mX.SetColumn(1, x);
		Product(ER, mX, &D);
	}
	else
	{
		BzzMatrix mX(NCR,2);
		mX.SetColumn(1, x);
		mX[jCarbonR][1] = 0.;
		mX[jCarbonR][2] = x[jCarbonR];
		Product(ER, mX, &D);
	}

	// Q
	GiveMeQ(Ntot,x);

	// E
	Q_Factorized = Q;
	Solve(Q_Factorized, -D, &EE);

	// A
	for(int m=1;m<=nPhases;m++)
		for(int n=1;n<=nPhases;n++)
			for(int j=1;j<=NER;j++)
				A[m][n] += D[j][m]*EE[j][n];

	// DeltaS
	double sum = 0.;
	for(m=1;m<=nPhases;m++)
		for(int n=1;n<=nPhases;n++)
			sum += A[n][m]*r[n]*r[m];
	deltas = -alfa/sum;

	// Ntot
	r *= deltas;
	BzzVector deltar(nPhases);
	BzzVector NtotNew(nPhases);
	for (int k=1;k<=100;k++)
	{
		deltar = r;
		r /= pow(2.,double(k-1));
		NtotNew = Ntot;
		NtotNew +=deltar;
		if (NtotNew[2]>=0.)
			break;
	}
	Ntot = NtotNew;

	// x
	for(j=1;j<=NCR;j++)
		x[j] = exp(-gtildeR[j] + Dot( ERT_Row[j],lambda) );

	// Z and W
	if (nPhases == 1)
	{
		Z[1] = x.GetSumElements();
		W    = Ntot[1]*(Z[1]-1.) - Dot(lambda, pR);
	}

	if (nPhases == 2)
	{
		Z[1]  = x.GetSumElements();
		Z[1] -= x[jCarbonR];
		Z[2]  = x[jCarbonR];
		W     = Ntot[1]*(Z[1]-1.) + Ntot[2]*(Z[2]-1.) - Dot(lambda, pR);
	}

	if (iVerboseFile == true && nPhases == 1)
	{
		fLog << "----------------------------------------------------------------------"	<< endl;
		fLog << " Mode 2bis" << "\t" << "#" << countMode2bis << "                            " 	<< endl;
		fLog << "----------------------------------------------------------------------"	<< endl;
		fLog << " W\t"<< Wold			<< "\t" << W		<< "\t" << W-Wold		<< "\t" << (W-Wold)/W*100.						<< "%" << endl;
		fLog << " Z\t" << Zold[1]		<< "\t" << Z[1]		<< "\t" << Z[1]-Zold[1] << "\t" << (Z[1]-Zold[1])/Z[1]*100.				<< "%" << endl;
		fLog << " N\t" << Ntot[1]-r[1]	<< "\t" << Ntot[1]	<< "\t" << r[1]			<< "\t" << r[1]/Ntot[1]*100.					<< "%" << endl;
		for(j=1;j<=NER;j++)
			fLog << setw(8) <<  " lambda_"	<< gas->list_of_elements[indexIncludeElement[j]-1] << "\t\t" << lambda[j] << endl;
		fLog << endl << endl;
	}

	if (iVerboseFile == true && nPhases == 2)
	{
		fLog << "----------------------------------------------------------------------"	<< endl;
		fLog << " Mode 2bis" << "\t" << "#" << countMode2bis << "                            " 	<< endl;
		fLog << "----------------------------------------------------------------------"	<< endl;
		fLog << " W\t"	  << Wold						<< "\t" << W				<< "\t" << W-Wold		<< "\t" << (W-Wold)/W*100.						<< "%" << endl;
		fLog << " ZGas\t" << Zold[1]					<< "\t" << Z[1]				<< "\t" << Z[1]-Zold[1] << "\t" << (Z[1]-Zold[1])/Z[1]*100.				<< "%" << endl;
		fLog << " ZSol\t" << Zold[2]					<< "\t" << Z[2]				<< "\t" << Z[2]-Zold[2] << "\t" << (Z[2]-Zold[2])/Z[2]*100.				<< "%" << endl;
		fLog << " Ngas\t" << Ntot[1]-r[1]				<< "\t" << Ntot[1]			<< "\t" << r[1]			<< "\t" << r[1]/Ntot[1]*100.					<< "%" << endl;
		fLog << " Nsol\t" << Ntot[2]-r[2]				<< "\t" << Ntot[2]			<< "\t" << r[2]			<< "\t" << r[2]/Ntot[2]*100.					<< "%" << endl;
		fLog << " Ntot\t" << Ntot[1]+Ntot[2]-r[1]-r[2]	<< "\t" << Ntot[1]+Ntot[2]	<< "\t" << r[1]+r[2]	<< "\t" << (r[1]+r[2])/(Ntot[1]+Ntot[2])*100.	<< "%" << endl;
		for(j=1;j<=NER;j++)
			fLog << setw(8) <<  " lambda_"	<< gas->list_of_elements[indexIncludeElement[j]-1] << "\t\t" << lambda[j] << endl;
		fLog << endl << endl;
	}

	if ( fabs(W-Wold) < 1.e-7 )			{ flag = CONVERGED; return; }	// Convergence was reached
	if ( W < Wold*1.01 )				{ flag = GO_ON;		return; }	
	if ( W > Wold*1.01 )				{ flag = FAILED;	return; }	
}


void OpenSMOKE_EquilibriumStanjan::Mode3(BzzVector &Ntot, BzzVector &lambda, BzzVector &x, exit_status &flag)
{
	countMode3++;

	int k,j;
	
	double W;
	double Wold;
	BzzVector Z(nPhases);
	BzzVector Zold(nPhases);

	// Z and W
	if (nPhases == 1)
	{
		Zold[1] = x.GetSumElements();
		Wold    = Ntot[1]*(Zold[1]-1.) - Dot(lambda, pR);
	}

	if (nPhases == 2)
	{
		Zold[1] = x.GetSumElements();
		Zold[1] -= x[jCarbonR];
		Zold[2]  = x[jCarbonR];
		Wold     = Ntot[1]*(Zold[1]-1.) + Ntot[2]*(Zold[2]-1.) - Dot(lambda, pR);
	}

	// Vector H
	GiveMeH(Ntot,x);

	// Vector D
	if (nPhases == 1)
	{
		BzzMatrix mX(NCR,1);
		mX.SetColumn(1, x);
		Product(ER, mX, &D);
	}
	else
	{
		BzzMatrix mX(NCR,2);
		mX.SetColumn(1, x);
		mX[jCarbonR][1] = 0.;
		mX[jCarbonR][2] = x[jCarbonR];
		Product(ER, mX, &D);
	}

	// Q
	GiveMeQ(Ntot, x);

	// Linear System
	for(j=1;j<=NER;j++)
		for(k=1;k<=NER;k++)
			QExtended_Factorized[j][k] = Q[j][k];
	for(int m=1;m<=nPhases;m++)
		for(j=1;j<=NER;j++)
		{
			QExtended_Factorized[NER+m][j] = D[j][m];
			QExtended_Factorized[j][NER+m] = D[j][m];
		}
	for(j=1;j<=NER;j++)
		bExtended[j] = -H[j];
	for(m=1;m<=nPhases;m++)
		bExtended[NER+m] = 1.-Zold[m];

	Solve(QExtended_Factorized, bExtended, &deltaExtended);

	// lambda
	for(j=1;j<=NER;j++)
		lambda[j] += deltaExtended[j];
	for(m=1;m<=nPhases;m++)
		Ntot[m] += deltaExtended[NER+m];

	// x
	for(j=1;j<=NCR;j++)
		x[j] = exp( -gtildeR[j] + Dot(ERT_Row[j],lambda) );

	// Z and W
	if (nPhases == 1)
	{
		Z[1] = x.GetSumElements();
		W    = Ntot[1]*(Z[1]-1.) - Dot(lambda, pR);
	}

	if (nPhases == 2)
	{
		Z[1]  = x.GetSumElements();
		Z[1] -= x[jCarbonR];
		Z[2]  = x[jCarbonR];
		W     = Ntot[1]*(Z[1]-1.) + Ntot[2]*(Z[2]-1.) - Dot(lambda, pR);
	}


	if (iVerboseFile == true && nPhases == 1)
	{
		BzzVector NtotOld = Ntot;
		BzzVector dNtot(nPhases);
		for(m=1;m<=nPhases;m++)	NtotOld[m] -= deltaExtended[NER+m];
		for(m=1;m<=nPhases;m++)	dNtot[m]    = deltaExtended[NER+m];

		fLog << "----------------------------------------------------------------------"	<< endl;
		fLog << " Mode 3" << "\t" << "#" << countMode3 << "                            " 	<< endl;
		fLog << "----------------------------------------------------------------------"	<< endl;
		fLog << " W\t"<< Wold			<< "\t" << W		<< "\t" << W-Wold		<< "\t" << (W-Wold)/W*100.				<< "%" << endl;
		fLog << " Z\t" << Zold[1]		<< "\t" << Z[1]		<< "\t" << Z[1]-Zold[1] << "\t" << (Z[1]-Zold[1])/Z[1]*100.		<< "%" << endl;
		fLog << " N\t" << NtotOld[1]	<< "\t" << Ntot[1]	<< "\t" << dNtot[1]		<< "\t" << dNtot[1]/Ntot[1]*100.		<< "%" << endl;
		for(j=1;j<=NER;j++)
			fLog << setw(8) <<  " lambda_"	<< gas->list_of_elements[indexIncludeElement[j]-1] << "\t\t" << lambda[j] << endl;
		fLog << endl << endl;
	}

	if (iVerboseFile == true && nPhases == 2)
	{
		BzzVector NtotOld = Ntot;
		BzzVector dNtot(nPhases);
		for(m=1;m<=nPhases;m++)	NtotOld[m] -= deltaExtended[NER+m];
		for(m=1;m<=nPhases;m++)	dNtot[m]    = deltaExtended[NER+m];
		
		fLog << "----------------------------------------------------------------------"	<< endl;
		fLog << " Mode 3" << "\t" << "#" << countMode3 << "                            " 	<< endl;
		fLog << "----------------------------------------------------------------------"	<< endl;
		fLog << " W\t"	  << Wold					<< "\t" << W				<< "\t" << W-Wold				<< "\t" << (W-Wold)/W*100.								<< "%" << endl;
		fLog << " ZGas\t" << Zold[1]				<< "\t" << Z[1]				<< "\t" << Z[1]-Zold[1]			<< "\t" << (Z[1]-Zold[1])/Z[1]*100.						<< "%" << endl;
		fLog << " ZSol\t" << Zold[2]				<< "\t" << Z[2]				<< "\t" << Z[2]-Zold[2]			<< "\t" << (Z[2]-Zold[2])/Z[2]*100.						<< "%" << endl;
		fLog << " Ngas\t" << NtotOld[1]				<< "\t" << Ntot[1]			<< "\t" << dNtot[1]				<< "\t" << dNtot[1]/Ntot[1]*100.						<< "%" << endl;
		fLog << " Nsol\t" << NtotOld[2]				<< "\t" << Ntot[2]			<< "\t" << dNtot[2]				<< "\t" << dNtot[1]/Ntot[2]*100.						<< "%" << endl;
		fLog << " Ntot\t" << NtotOld[1]+NtotOld[2]	<< "\t" << Ntot[1]+Ntot[2]	<< "\t" << dNtot[1]+dNtot[2]	<< "\t" << (dNtot[1]+dNtot[1])/(Ntot[1]+Ntot[2])*100.	<< "%" << endl;
		for(j=1;j<=NER;j++)
			fLog << setw(8) <<  " lambda_"	<< gas->list_of_elements[indexIncludeElement[j]-1] << "\t\t" << lambda[j] << endl;
		fLog << endl << endl;
	}

	if ( fabs(W-Wold) < 1.e-7 )			{ flag = CONVERGED; return; }	// Convergence was reached
	if ( W < Wold )						{ flag = GO_ON;		return; }	
	if ( W > Wold )						{ flag = FAILED;	return; }	
}


void OpenSMOKE_EquilibriumStanjan::GiveMeQ(BzzVector Ntot, BzzVector &x)
{
	Q = 0.;
	for(int i=1;i<=NER;i++)
		for(int k=1;k<=NER;k++)
			for(int j=1;j<=NCR;j++)
				Q[i][k] += Ntot[jPhase[j]]*ER[i][j]*ER[k][j]*x[j];
}

void OpenSMOKE_EquilibriumStanjan::GiveMeH(BzzVector Ntot, BzzVector &x)
{
	Product(ER, x, &H);
	H *= Ntot[1];
	H -= pR;
	
	if(nPhases == 2)
	{
		BzzVector aux = ER.GetColumn(jCarbonR);
		for(int j=1;j<=NER;j++)
		{
		//	cout << H[j] << "\t";
			H[j] += (Ntot[2]-Ntot[1])*aux[j];
		//	cout << H[j] << endl;
		}
	//	getchar();
	}
}

void OpenSMOKE_EquilibriumStanjan::LinearProgramming(BzzVector &Nsolution)
{
	int i,j;
	int icase;

	int m1 = 0;
	int m2 = 0;
	int m3 = NER;
	int N  = NCR;
	int M  = m1+m2+m3;

	Mat_IO_DP a(M+2, N+1);
	Vec_O_INT izrov(N);
	Vec_O_INT iposv(M);

	// Initializing matrix
	for (j=1;j<=M+1;j++)			
		for (i=1;i<=N+1;i++)				
			a[j-1][i-1] = 0.;

	// Assembling problem
	for (i=2;i<=N+1;i++)	a[0][i-1] = -gtildeR[i-1];	// Objective function
	for (i=2;i<=M+1;i++)	a[i-1][0] = pR[i-1];		// Constraints
	
	for (j=2;j<=M+1;j++)								// Constraints
		for (i=2;i<=N+1;i++)				
			a[j-1][i-1] = -ER[j-1][i-1];

	// Solving Linear Programming Problem
	NR::simplx(a, m1, m2, m3, icase, izrov, iposv);

	// Check Solution
	if (icase != 0)	ErrorMessage("Failure in Linear Programming (Min)!");

	// Solution
	for(j=0;j<=M-1;j++)	if (iposv[j] < N)	Nsolution[iposv[j]+1] = a[j+1][0];
	for(j=0;j<=N-1;j++)	if (izrov[j] < N)	Nsolution[izrov[j]+1] = 0.;

	// Write on File
	if (iVerboseFile == true)
	{
		fLog << "--------------------------------------------------------" << endl;
		fLog << "                Linear Programming (Min)                " << endl;
		fLog << "--------------------------------------------------------" << endl;

		double Ntot = Nsolution.GetSumElements();
		for(j=1;j<=N;j++)
			if (Nsolution[j] > threshold_element)
				fLog << setw(16) << gas->names[indexIncludeSpecies[j]] << "\t" << Nsolution[j] <<  "\t" << Nsolution[j]/Ntot << endl;
		
		fLog << "--------------------------------------------------------" << endl;
		fLog << setw(16) << " " << "\t" << Ntot << "\t" << 1. << endl;
		fLog << endl << endl;
	}
}

void OpenSMOKE_EquilibriumStanjan::LinearProgrammingMinMax(BzzVector &Nsolution)
{
	int i,j;
	int icase;

	BzzVector xSolution(NCR+1);
	BzzMatrix Amatrix(NCR, NCR+1);
	for (i=1;i<=NCR;i++)	Amatrix[i][i]		=	-1.;
	for (i=1;i<=NCR;i++)	Amatrix[i][NCR+1]	=	 1.;

	int m1 = NCR;
	int m2 = 0;
	int m3 = NER;

	int N = NCR+1;
	int M = m1+m2+m3;

	Mat_IO_DP a(M+2, N+1);
	Vec_O_INT izrov(N);
	Vec_O_INT iposv(M);

	for (i=1;i<=M+2;i++)
		for (j=1;j<=N+1;j++)
				a[i-1][j-1] = 0.;

	a[0][(N+1)-1] = 1.;
	for (i=2;i<=m1+1;i++)	a[i-1][0] = 0;
	for (i=m1+2;i<=M+1;i++)	a[i-1][0] = pR[i-(m1+2)+1];

	for (j=2;j<=m1+1;j++)
		for (i=2;i<=N+1;i++)
			a[j-1][i-1] = -Amatrix[j-1][i-1];

	for (j=m1+2;j<=M+1;j++)
		for (i=2;i<=N;i++)
			a[j-1][i-1] = -ER[j-(m1+2)+1][i-1];

	NR::simplx(a, m1, m2, m3, icase, izrov, iposv);

	// Check Solution
	if (icase != 0)	ErrorMessage("Failure in Linear Programming (Min-Max)!");

	// Solution
	for(j=0;j<=M-1;j++)	if (iposv[j] < N)	xSolution[iposv[j]+1] = a[j+1][0];
	for(j=0;j<=N-1;j++)	if (izrov[j] < N)	xSolution[izrov[j]+1] = 0.;
	for(j=1;j<=N-1;j++)	Nsolution[j] = xSolution[j];

	// Write on File
	if (iVerboseFile == true)
	{
		fLog << "--------------------------------------------------------" << endl;
		fLog << "              Linear Programming (Min-Max)              " << endl;
		fLog << "--------------------------------------------------------" << endl;

		double Ntot = Nsolution.GetSumElements();
		for(j=1;j<=N;j++)
			if (Nsolution[j] > threshold_element)
				fLog << setw(16) << gas->names[indexIncludeSpecies[j]] << "\t" << Nsolution[j] <<  "\t" << Nsolution[j]/Ntot << endl;
				
		fLog << "--------------------------------------------------------" << endl;
		fLog << setw(16) << " " << "\t" << Ntot << "\t" << 1. << endl;
		fLog << endl << endl;
	}
}

void OpenSMOKE_EquilibriumStanjan::Blending(BzzVector &NtotFirstGuess, BzzVector &xFirstGuess, BzzVector &Nmin, BzzVector &Nmm)
{
	int j;
	double f=0.10;

	for (j=1;j<=NCR;j++)
		xFirstGuess[j] = (1.-f)*Nmin[j] + f*Nmm[j];

	NtotFirstGuess[1] = xFirstGuess.GetSumElements();
	
	if (iSolidPhase == true) 
	{
		NtotFirstGuess[2]  = xFirstGuess[jCarbonR];	// Solid Phase
		NtotFirstGuess[1] -= NtotFirstGuess[2];		// Gas Phase
		xFirstGuess		  /= NtotFirstGuess[1];		// Gas Phase
		xFirstGuess[jCarbonR] = 1.0;				// Solid Phase
	}
	else
		xFirstGuess		  /= NtotFirstGuess[1];
}

exit_status OpenSMOKE_EquilibriumStanjan::ModuleMode1(BzzVector &Ntot, BzzVector &lambda, BzzVector &x)
{
	exit_status flag;
	
	BzzVector lambdaNew	= lambda;	// Copy of initial values
	BzzVector xNew		= x;		// Copy of initial values
	
	for (int k=1;k<=NmaxIterationsMode1;k++)
	{		
		// Iterations on Mode1
		Mode1(Ntot, lambdaNew, xNew, flag);

		// 2. No more changes
		if (flag == CONVERGED)
		{
			lambda	= lambdaNew;
			x		= xNew;
			return CONVERGED;
		}

		// The solution is worse than the previous one
		if (flag == FAILED)									
		{
		}
	}
	
	lambda	= lambdaNew;
	x		= xNew;

	return GO_ON;
}

exit_status OpenSMOKE_EquilibriumStanjan::ModuleMode2(BzzVector &Ntot, BzzVector &lambda, BzzVector &x)
{
	exit_status flag;
	
	BzzVector lambdaNew	= lambda;	// Copy of initial values
	BzzVector xNew		= x;		// Copy of initial values
	
	for (int k=1;k<=NmaxIterationsMode2;k++)
	{		
		// Iterations on Mode1
		Mode2(Ntot, lambdaNew, xNew, flag);

		// 2. No more changes
		if (flag == CONVERGED)
		{
			lambda	= lambdaNew;
			x		= xNew;
			return CONVERGED;
		}

		// The solution is worse than the previous one
		if (flag == FAILED)									
		{
			fLog << "FAIL" << endl; 
			WarningMessageStrong("Mode2: FAILED!");
			return FAILED;
		}
	}
	
	lambda	= lambdaNew;
	x		= xNew;

	return GO_ON;
}

exit_status OpenSMOKE_EquilibriumStanjan::Solution(BzzVector &Ntot, BzzVector &lambda, BzzVector &x)
{
	exit_status flag;
	
	ModuleMode1(Ntot, lambda, x);			// Steeping descent
	ModuleMode2(Ntot, lambda, x);			// Newton descent
	Mode2bis(Ntot, lambda, x, flag);		// Adjusting number of moles

	BzzVector NtotOld = Ntot;
	for (int k=1;k<=NmaxIterations;k++)
	{
		ModuleMode1(Ntot, lambda, x);			// Steeping descent
		ModuleMode2(Ntot, lambda, x);			// Newton descent
		Mode2bis(Ntot, lambda, x, flag);		// Adjusting number of moles

		if (fabs(Ntot[1]-NtotOld[1])/Ntot[1] <= toleranceSteeping)
			break;
		NtotOld = Ntot;
	}
	Mode3(Ntot, lambda, x, flag);
	Mode3(Ntot, lambda, x, flag);
	Mode3(Ntot, lambda, x, flag);
	
	return GO_ON;
}

void OpenSMOKE_EquilibriumStanjan::FinalSummary(BzzVector &Ntot, BzzVector &xR, const double P_Pa, const double T)
{
	int j;
	double MWmix;
	double cTot;

	ChangeDimensions(NC, &NGas);
	ChangeDimensions(NC, &xGas);
	ChangeDimensions(NC, &omegaGas);
	ChangeDimensions(NC, &cGas);

	// Total concentration [kmol/m3]
	cTot = P_Pa/(Constants::R_J_kmol*T);

	// Gas Phase
	{
		xGas = 0.;
		for(j=1;j<=NCR;j++)
			xGas[indexIncludeSpecies[j]] = xR[j];
		if (iSolidPhase == true) xGas[jCarbon] = 0.;
		gas->GetMWAndMassFractionsFromMoleFractions(MWmix, omegaGas, xGas);
		cGas = cTot*xGas;
		NGas = Ntot[1]*xGas;
	}
	
	if (iVerboseFile == true)
	{
		//double Nsolid = 
		if (iSolidPhase == true)
		{
			fLog << "---------------------------------------------------------------------------------------" << endl; 
			fLog << setw(16) << "Sol. Phase  Name" << left << setw(15) << "  N[kmol]" << setw(15) << " x" << setw(15) << "omega" << setw(15) << "c[kmol/m3]" << setw(15) << "g[-]" << endl;
			fLog << "---------------------------------------------------------------------------------------" << endl; 
			fLog << right << setw(16)			<< gas->names[jCarbon]	<< "\t" 
						  << Ntot[2]			<< "\t"
						  << 1.					<< "\t"
						  << 1.					<< "\t"
						  << "############"		<< "\t"
						  << gtilde[jCarbon]	<< "\t"
						  << endl << endl;;
		}

		fLog << "---------------------------------------------------------------------------------------" << endl; 
			fLog << setw(16) << "Gas Phase   Name" << left << setw(15) << "  N[kmol]" << setw(15) << " x" << setw(15) << "omega" << setw(15) << "c[kmol/m3]" << setw(15) << "g[-]" << endl;
			fLog << "---------------------------------------------------------------------------------------" << endl; 
		for(j=1;j<=NC;j++)
			fLog << right << setw(16)	<< gas->names[j]	<< "\t" 
										<< NGas[j]			<< "\t"
										<< xGas[j]			<< "\t"
										<< omegaGas[j]		<< "\t"
										<< cGas[j]			<< "\t"
										<< gtilde[j]		<< "\t"
										<< endl;


	}

}
