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
#include "nr.h"

#include "basic/OpenSMOKE_Constants.h"
#include "engine/OpenSMOKE_IdealGas.h"
#include "engine/OpenSMOKE_EquilibriumStanjan.h"

OpenSMOKE_EquilibriumStanjan *ptEquilibrium;

const double OpenSMOKE_EquilibriumStanjan::threshold_element	=	1.e-16;

void OpenSMOKE_EquilibriumStanjan::ErrorMessage(const std::string message)
{
    cout << endl;
    cout << "Class:  OpenSMOKE_EquilibriumStanjan"	<< endl;
    cout << "Object: " << name_object				<< endl;
    cout << "Error:  " << message					<< endl;
    cout << "Press enter to continue... "			<< endl;
    getchar();
    exit(-1);
}

void OpenSMOKE_EquilibriumStanjan::WarningMessageStrong(const std::string message)
{
    cout << endl;
    cout << "Class:	  OpenSMOKE_EquilibriumStanjan"	<< endl;
    cout << "Object:  " << name_object				<< endl;
    cout << "Warning: " << message					<< endl;
    cout << "Press enter to continue..."			<< endl;
	getchar();
	cout << endl;
}

void OpenSMOKE_EquilibriumStanjan::WarningMessageSoft(const std::string message)
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
	iVerboseFile				 = false;

	NmaxIterationsMode1 = 500;
	NmaxIterationsMode2 = 30;
	NmaxIterations      = 30;

	toleranceMode1		= 1.e-7;
	toleranceMode2		= 1.e-8;
	toleranceSteeping	= 1.e-13;

	iAlreadySolved		= false;
}

void OpenSMOKE_EquilibriumStanjan::Reset()
{
	iAlreadySolved = false;
}

void OpenSMOKE_EquilibriumStanjan::SetName(std::string name)
{
	name_object = name;
}

void OpenSMOKE_EquilibriumStanjan::SetVerbose()
{
	iVerboseFile = true;
}

void OpenSMOKE_EquilibriumStanjan::UnsetVerbose()
{
	iVerboseFile = false;
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

	if (iVerboseFile == true)
	{
		openOutputFileAndControl(fLog, "Equilibrium.log");
		fLog.setf(ios::scientific);
	}
}

void OpenSMOKE_EquilibriumStanjan::SetElementalCompositionFromSpeciesMoleFractions(BzzVector &_xInitial)
{
	BzzVector _xInitialElemental(NE);
	Product(E, _xInitial, &_xInitialElemental);
	SetElementalComposition(_xInitialElemental);
	Reset();
}

void OpenSMOKE_EquilibriumStanjan::SetElementalComposition(BzzVector &_xInitialElemental)
{
	// Elemental composition
	p = _xInitialElemental;
	
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
	{	
		int i;
		for(i=NC;i>=1;i--)
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
	NCR	= indexIncludeSpecies.Size();

	ReorderingSpecies();

	// Memory Allocation
	ChangeDimensions(NER, &lambdaFinalSolutionR);
	ChangeDimensions(NCR, &gtildeR);

	ChangeDimensions(NER, &H);
	ChangeDimensions(NER, &f);
	ChangeDimensions(NER, &deltalambda);
	ChangeDimensions(NER, NER, &Q);

	ChangeDimensions(NER,&D);
	ChangeDimensions(NER,&EE);
	ChangeDimensions(NER+1, NER+1, &QExtended_Factorized);
	ChangeDimensions(NER+1, &bExtended);
	ChangeDimensions(NER+1, &deltaExtended);

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

	// Assigned elemental composition
	assignedElementalComposition = true;

	// Reset
	Reset();
}

void OpenSMOKE_EquilibriumStanjan::ReorderingSpecies()
{
	NCR_Gas		= NCR;
}

void OpenSMOKE_EquilibriumStanjan::Equilibrate(const double _T, const double _P_Pa, BzzVector &xFinal, double &NFinal)
{
	countMode1		= 0;
	countMode2		= 0;
	countMode2bis	= 0;
	countMode3		= 0;

	// Get Free Gibbs Energy for the species[-]
	GetGtilde(_T, _P_Pa);

	// Get Lambda First Guess
	GetFirstGuess(_T, _P_Pa);

	// Final
	xFinal = xGas;
	NFinal = Database_NtotFirstGuess;

	// Already Solved
	iAlreadySolved = false;
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
	if (iAlreadySolved == false)
	{
		// First Guess Lambda
		double NtotFirstGuess; 
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

		double	Ntot;
		Ntot				= NtotFirstGuess;
		BzzVector	lambdaR	= lambdaFirstGuessR;
		BzzVector	xR		= xFirstGuessR;

		// Solution
		Solution(Ntot, lambdaR, xR);

		// Final Summary
		FinalSummary(Ntot, xR, T, P_Pa);

		Database_NtotFirstGuess			= Ntot;
		Database_lambdaFirstGuessR		= lambdaR;
		Database_xFirstGuessR			= xR;
	}
	else
	{
		double	Ntot;
		Ntot				= Database_NtotFirstGuess;
		BzzVector	lambdaR	= Database_lambdaFirstGuessR;
		BzzVector	xR		= Database_xFirstGuessR;

		// Solution
		Solution(Ntot, lambdaR, xR);

		// Final Summary
		FinalSummary(Ntot, xR, T, P_Pa);

		Database_NtotFirstGuess			= Ntot;
		Database_lambdaFirstGuessR		= lambdaR;
		Database_xFirstGuessR			= xR;
	}
}

void OpenSMOKE_EquilibriumStanjan::GetGtilde(const double T, const double P_Pa)
{
	int i;

	gas->GetStandardGibbsFreeEnergy_Mass(gtilde, T);
	gtilde /= (Constants::R_J_kmol*T);
	ElementByElementProduct(gtilde, gas->M, &gtilde);

	for (i=1;i<=NCR;i++)
		gtildeR[i] = gtilde[indexIncludeSpecies[i]];

	double additional_term = log(P_Pa/Constants::P_Reference);
	for (i=1;i<=NCR;i++)
		gtildeR[i] += additional_term;
}



void OpenSMOKE_EquilibriumStanjan::Mode1(double &Ntot, BzzVector &lambda, BzzVector &x, exit_status &flag)
{
	countMode1++;

	int    maxNumberOfSubIterations = 50;
	double maxIncreasingFactor		= 1.01;
	double reductionFactor			= 0.25;
	double minBeta					= 1.e-20;

	int		iter = 0;
	int		k,j;
	double	Beta;
	double	deltas;
	double	Wold;
	double	W;
	double  Zold;
	double  Z;
	BzzVector lambdaNew(NER);


	// Z and W
	GiveMeZandW(x, lambda, Ntot, Zold, Wold);	
	
	// Vector H
	GiveMeH(Ntot,x);

	// Beta
	Beta = -sqrt(Dot(H,H));

	// Perturbation (in case of first iteration)
	if (Beta >= -minBeta && countMode1 == 1)
	{
		Ntot *= 1.001;	
		GiveMeZandW(x, lambda, Ntot, Zold, Wold);	
		GiveMeH(Ntot,x);
		Beta = -sqrt(Dot(H,H));
	}

	if(Beta < -minBeta)	// otherwise convergence is verified
	{
		// Matrix Q
		GiveMeQ(Ntot,x);

		// Vector f
		for(k=1;k<=NER;k++)
			f[k] = H[k]/Beta;

		// DeltaS
		double csi=0.;
		for(k=1;k<=NER;k++)
			for(j=1;j<=NER;j++)
				csi += Q[j][k]*f[j]*f[k];
		deltas = -Beta/csi;

		for (iter=1;iter<=maxNumberOfSubIterations;iter++)
		{
			// lambda
			for(j=1;j<=NER;j++)
				lambdaNew[j] = lambda[j] + f[j]*deltas;

			// x
			for(j=1;j<=NCR;j++)
				x[j] = exp(-gtildeR[j] + Dot( ERT_Row[j],lambdaNew) );

			// Z and W
			GiveMeZandW(x, lambdaNew, Ntot, Z, W);

			if (W > maxIncreasingFactor*Wold)	
				deltas *= reductionFactor;		// Reject
			else
			{
				lambda = lambdaNew;				// Accept
				break;
			}
		}
	}
	else
	{
		deltas	= 0.;
		Z		= Zold;
		W		= Wold;
	}
		
	if ( fabs((W-Wold)/W) < toleranceMode1 )	flag =  CONVERGED;		// Convergence was reached
	else if ( W < Wold )						flag =  GO_ON;	
	else if ( W >= Wold )						flag =  FAILED;	


	if (iVerboseFile == true)
	{
		std::string status;
		if (flag == CONVERGED)	status = "CONVERGED";
		if (flag == GO_ON)		status = "GO_ON";
		if (flag == FAILED)		status = "FAILED";

		fLog << "----------------------------------------------------------------------"	<< endl;
		fLog << " Mode 1" << "\t" << "#" << countMode1 << " - " << status					<< endl;
		fLog << "----------------------------------------------------------------------"	<< endl;
		fLog << "  Number of iterations: " << iter << "/" << maxNumberOfSubIterations << endl;
		fLog << "  Beta:     " << setw(16) << right << Beta   << endl;
		fLog << "  DeltaS:   " << setw(16) << right << deltas << endl;
		for(j=1;j<=NER;j++)
		{
			fLog <<  "  H_" << setw(6) << left << gas->list_of_elements[indexIncludeElement[j]-1] << H[j] << "   ";
			fLog <<  "  lambda_" << setw(10) << left << gas->list_of_elements[indexIncludeElement[j]-1] << lambda[j] << endl;
		}
		fLog << endl;

		fLog << "  W   " << setw(16) << right << Wold << setw(16) << right << W    << setw(16) << right << W-Wold << setw(16) << right << (W-Wold)/W*100. << "%" << endl;
		fLog << "  Z   " << setw(16) << right << Zold << setw(16) << right << Z    << setw(16) << right << Z-Zold << setw(16) << right << (Z-Zold)/Z*100. << "%" << endl;	
		fLog << "  N   " << setw(16) << right << Ntot << setw(16) << right << Ntot << setw(16) << right << 0.	  << setw(16) << right << 0.			  << "%" << endl;
	
		fLog << endl << endl;
	}

	return;
}


void OpenSMOKE_EquilibriumStanjan::Mode2(double &Ntot, BzzVector &lambda, BzzVector &x, exit_status &flag)
{
	countMode2++;

	int j;
	double W;
	double Wold;
	double Zold;
	double Z;

	// Z and W
	GiveMeZandW(x, lambda, Ntot, Zold, Wold);

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
	fLog << x.GetSumElements() << endl;
	
	// Z and W
	GiveMeZandW(x, lambda, Ntot, Z, W);

	if (iVerboseFile == true)
	{
		fLog << "----------------------------------------------------------------------"	<< endl;
		fLog << " Mode 2" << "\t" << "#" << countMode2 << "                            " 	<< endl;
		fLog << "----------------------------------------------------------------------"	<< endl;
//		fLog << "  Number of iterations: " << iter << "/" << maxNumberOfSubIterations << endl;
//		fLog << "  Beta:     " << setw(16) << right << Beta   << endl;
//		fLog << "  DeltaS:   " << setw(16) << right << deltas << endl;
		for(j=1;j<=NER;j++)
		{
			fLog <<  "  H_" << setw(6) << left << gas->list_of_elements[indexIncludeElement[j]-1] << H[j] << "   ";
			fLog <<  "  dlambda_" << setw(6) << left << gas->list_of_elements[indexIncludeElement[j]-1] << deltalambda[j] << "   ";
			fLog <<  "  lambda_"  << setw(10) << left << gas->list_of_elements[indexIncludeElement[j]-1] << lambda[j] << endl;
		}
		fLog << endl;

		fLog << " W\t" << Wold << "\t" << W    << "\t" << W-Wold << "\t" << (W-Wold)/W*100. << "%" << endl;
		
		fLog << " Z" << "\t" << Zold << "\t" << Z    << "\t" << Z-Zold << "\t" << (Z-Zold)/Z*100. << "%" << endl;
		
		fLog << " N" << "\t" << Ntot << "\t" << Ntot << "\t" << 0.			<< "\t" << 0.						<< "%" << endl;
		
		for(j=1;j<=NER;j++)
			fLog << setw(8) <<  " lambda_"	<< gas->list_of_elements[indexIncludeElement[j]-1] << "\t\t" << lambda[j] << endl;
		
		fLog << endl << endl;
	}

	if ( fabs((W-Wold)/W) < toleranceMode2 )	{ flag = CONVERGED; return; }	// Convergence was reached
	if ( W > Wold*1.001 )						{ flag = FAILED;	return; }	
	if ( W < Wold*1.001 )						{ flag = GO_ON;		return; }	
}

void OpenSMOKE_EquilibriumStanjan::Mode3(double &Ntot, BzzVector &lambda, BzzVector &x, exit_status &flag)
{
	countMode3++;

	int k,j;
	
	double W;
	double Wold;
	double Zold;
	double Z;

	// Z and W
	GiveMeZandW(x, lambda, Ntot, Zold, Wold);

	// Vector H
	GiveMeH(Ntot,x);

	// Vector D
	GiveMeD(x);

	// Q
	GiveMeQ(Ntot, x);

	// Linear System
	for(j=1;j<=NER;j++)
		for(k=1;k<=NER;k++)
			QExtended_Factorized[j][k] = Q[j][k];
		for(j=1;j<=NER;j++)
		{
			QExtended_Factorized[NER+1][j] = D[j];
			QExtended_Factorized[j][NER+1] = D[j];
		}
	for(j=1;j<=NER;j++)
		bExtended[j] = -H[j];
	bExtended[NER+1] = 1.-Zold;

	Solve(QExtended_Factorized, bExtended, &deltaExtended);

	// lambda
	for(j=1;j<=NER;j++)
		lambda[j] += deltaExtended[j];
	Ntot += deltaExtended[NER+1];

	// x
	for(j=1;j<=NCR;j++)
		x[j] = exp( -gtildeR[j] + Dot(ERT_Row[j],lambda) );

	// Z and W
	GiveMeZandW(x, lambda, Ntot, Z, W);

	if (iVerboseFile == true)
	{
		double NtotOld = Ntot;
		double dNtot;
		NtotOld -= deltaExtended[NER+1];
		dNtot    = deltaExtended[NER+1];

		fLog << "----------------------------------------------------------------------"	<< endl;
		fLog << " Mode 3" << "\t" << "#" << countMode3 << "                            " 	<< endl;
		fLog << "----------------------------------------------------------------------"	<< endl;
		fLog << " W\t"<< Wold			<< "\t" << W		<< "\t" << W-Wold		<< "\t" << (W-Wold)/W*100.				<< "%" << endl;
		
		fLog << " Z" << "\t" << Zold	<< "\t" << Z		<< "\t" << Z-Zold << "\t" << (Z-Zold)/Z*100.		<< "%" << endl;

		fLog << " N" << "\t" << NtotOld	<< "\t" << Ntot	<< "\t" << dNtot	<< "\t" << dNtot/Ntot*100.		<< "%" << endl;
		
		for(j=1;j<=NER;j++)
			fLog << setw(8) <<  " lambda_"	<< gas->list_of_elements[indexIncludeElement[j]-1] << "\t\t" << lambda[j] << endl;
		
		fLog << endl << endl;
	}

	if ( fabs(W-Wold) < 1.e-7 )			{ flag = CONVERGED; return; }	// Convergence was reached
	if ( W < Wold )						{ flag = GO_ON;		return; }	
	if ( W > Wold )						{ flag = FAILED;	return; }	
}


void OpenSMOKE_EquilibriumStanjan::GiveMeQ(double &Ntot, BzzVector &x)
{
	Q = 0.;
	for(int i=1;i<=NER;i++)
		for(int k=1;k<=NER;k++)
			for(int j=1;j<=NCR_Gas;j++)
				Q[i][k] += ER[i][j]*ER[k][j]*x[j];
	Q *= Ntot;
}

void OpenSMOKE_EquilibriumStanjan::GiveMeH(double &Ntot, BzzVector &x)
{
	H=0.;
	for(int i=1;i<=NER;i++)
		for(int j=1;j<=NCR_Gas;j++)
			H[i] += Ntot*ER[i][j]*x[j];
	H -= pR;
}

void OpenSMOKE_EquilibriumStanjan::GiveMeD(BzzVector &x)
{
	D =0.;
	for(int i=1;i<=NER;i++)
		for(int j=1;j<=NCR_Gas;j++)
			D[i] += ER[i][j]*x[j];
}

void OpenSMOKE_EquilibriumStanjan::GiveMeA()
{
	A = Dot(D,EE);
}

void OpenSMOKE_EquilibriumStanjan::GiveMeZandW(BzzVector &x, BzzVector &lambda, 
											   double &Ntot, double &Z, double &W)	
{
	Z =x.GetSumElements();
	W = Ntot*(Z-1.)-Dot(lambda, pR);
}

void OpenSMOKE_EquilibriumStanjan::GiveMeZandWandV(BzzVector &x, BzzVector &lambda, double &Ntot,
												   double &Z, double &W, double &V)	
{
	GiveMeZandW(x, lambda, Ntot, Z, W);
	V = Z-1.;
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

		for(j=1;j<=NCR;j++)
			if (Nsolution[j] > threshold_element)
				fLog << setw(16) << gas->names[indexIncludeSpecies[j]] << "\t" << Nsolution[j] <<  "\t" << Nsolution[j]/Ntot << endl;
				
		fLog << "--------------------------------------------------------" << endl;
		fLog << setw(16) << " " << "\t" << Ntot << "\t" << 1. << endl;
		fLog << endl << endl;
	}
}

void OpenSMOKE_EquilibriumStanjan::Blending(double &NtotFirstGuess, BzzVector &xFirstGuess, BzzVector &Nmin, BzzVector &Nmm)
{
	int j;
	double f=0.10;

	for (j=1;j<=NCR;j++)
		xFirstGuess[j] = (1.-f)*Nmin[j] + f*Nmm[j];

	NtotFirstGuess = 0.;
	for (j=1;j<=NCR_Gas;j++)	
		NtotFirstGuess += xFirstGuess[j];

	for (j=1;j<=NCR_Gas;j++)	
		xFirstGuess[j] /= NtotFirstGuess;
}

exit_status OpenSMOKE_EquilibriumStanjan::ModuleMode1(double &Ntot, BzzVector &lambda, BzzVector &x)
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

exit_status OpenSMOKE_EquilibriumStanjan::ModuleMode2(double &Ntot, BzzVector &lambda, BzzVector &x)
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
		//	fLog << "FAIL" << endl; 
		//	WarningMessageStrong("Mode2: FAILED!");
			return FAILED;
		}
	}
	
	lambda	= lambdaNew;
	x		= xNew;

	return GO_ON;
}

exit_status OpenSMOKE_EquilibriumStanjan::Solution(double &Ntot, BzzVector &lambda, BzzVector &x)
{
	exit_status flag;

//	ModuleMode1(Ntot, lambda, x);			// Steeping descent
//	ModuleMode2(Ntot, lambda, x);			// Newton descent
//	ModuleMode2bis(Ntot, lambda, x);		// Steeping ascent

	double NtotOld = Ntot;
	for (int k=1;k<=NmaxIterations;k++)
	{
		ModuleMode1(Ntot, lambda, x);			// Steeping descent
		ModuleMode2(Ntot, lambda, x);			// Newton descent
		ModuleMode2bis(Ntot, lambda, x);		// Steeping ascent

		if (fabs(Ntot-NtotOld)/Ntot <= toleranceSteeping)
			break;
		NtotOld = Ntot;
	}

//	for (int k=1;k<=10;k++)
//		Mode3(Ntot, lambda, x, flag);
	
	return GO_ON;
}

void OpenSMOKE_EquilibriumStanjan::FinalSummary(double &Ntot, BzzVector &xR, const double P_Pa, const double T)
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
		for(j=1;j<=NCR_Gas;j++)
			xGas[indexIncludeSpecies[j]] = xR[j];
		gas->GetMWAndMassFractionsFromMoleFractions(MWmix, omegaGas, xGas);
		cGas = cTot*xGas;
		NGas = Ntot*xGas;
	}
	
	if (iVerboseFile == true)
	{ 
		fLog << "---------------------------------------------------------------------------------------" << endl; 
			fLog << setw(16) << "Gas Phase   Name" << left << setw(15) << "  N[kmol]" << setw(15) << " x" << setw(15) << "omega" << setw(15) << "c[kmol/m3]" << setw(15) << "g[-]" << endl;
			fLog << "---------------------------------------------------------------------------------------" << endl; 
		for(j=1;j<=NCR_Gas;j++)
			fLog << right << setw(16)	<< gas->names[indexIncludeSpecies[j]]	<< "\t" 
										<< NGas[indexIncludeSpecies[j]]			<< "\t"
										<< xGas[indexIncludeSpecies[j]]			<< "\t"
										<< omegaGas[indexIncludeSpecies[j]]		<< "\t"
										<< cGas[indexIncludeSpecies[j]]			<< "\t"
										<< gtilde[indexIncludeSpecies[j]]		<< "\t"
										<< endl;


	}

}

void OpenSMOKE_EquilibriumStanjan::Mode2bis(double &Ntot, BzzVector &lambda, BzzVector &x, exit_status &flag)
{
	countMode2bis++;

	int j;
	double alfa;
	double deltaN;
	double W;
	double Wold;
	double Zold;
	double Z;
	double V;
	double r;

	const double alfa_min = 1.e-16;

	// Z and W
	Zold = x.GetSumElements();
	V    = Zold-1.;
	Wold = Ntot*V - Dot(lambda, pR);

	// alfa
	alfa = sqrt(V*V);

	if (alfa > alfa_min)	
	{	
		// r
		r = V/alfa;

		// Q
		GiveMeQ(Ntot,x);

		// Matrix EE
		GiveMeD(x);
		Q_Factorized = Q;
		Solve(Q_Factorized, -D, &EE);

		// DeltaN
		GiveMeA();
		deltaN = -alfa/(A*r*r);

		// Ntot
		Ntot += r*deltaN;

		// x
		for(j=1;j<=NCR_Gas;j++)
			x[j] = exp(-gtildeR[j] + Dot( ERT_Row[j],lambda) );

		// Z and W
		Z = x.GetSumElements();
		W = Ntot*(Z-1.) - Dot(lambda, pR);
	}
	else
	{
		W = Wold;
		Z = Zold;
	}

	if (iVerboseFile == true)
	{
		fLog << "----------------------------------------------------------------------"	<< endl;
		fLog << " Mode 2bis - SP" << "\t" << "#" << countMode2bis << "                            " 	<< endl;
		fLog << "----------------------------------------------------------------------"	<< endl;
		fLog << "r " << r << " alfa " << alfa << " V " << V << endl;
		fLog << " W\t"<< Wold			<< "\t" << W		<< "\t" << W-Wold		<< "\t" << (W-Wold)/W*100.		<< "%" << endl;
		
		fLog << " Z" << 1 << "\t" << Zold << "\t" << Z    << "\t" << Z-Zold << "\t" << (Z-Zold)/Z*100. << "%" << endl;	
		fLog << " N" << 1 << "\t" << Ntot-r*deltaN<< "\t" << Ntot << "\t" << r*deltaN << "\t" << r*deltaN/Ntot*100.	<< "%" << endl;

		for(j=1;j<=NER;j++)
			fLog << setw(8) <<  " lambda_"	<< gas->list_of_elements[indexIncludeElement[j]-1] << "\t\t" << lambda[j] << endl;

		fLog << endl << endl;
	}

	if ( fabs(W-Wold) < 1.e-7 )			{ flag = CONVERGED; return; }	// Convergence was reached
	if ( W < Wold )						{ flag = GO_ON;		return; }	
	if ( W > Wold )						{ flag = FAILED;	return; }	
}

exit_status OpenSMOKE_EquilibriumStanjan::ModuleMode2bis(double &Ntot, BzzVector &lambda, BzzVector &x)
{
	exit_status flag;

	for(int k=1;k<=1;k++)
		Mode2bis(Ntot, lambda, x, flag);		// No solid phase: Adjusting number of moles

	return flag;
}