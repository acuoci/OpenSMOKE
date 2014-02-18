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
#if LINUX_SO==0
	#include "nr.h"
#endif

#include "basic/OpenSMOKE_Constants.h"
#include "engine/OpenSMOKE_IdealGas.h"
#include "engine/OpenSMOKE_Equilibrium.h"

OpenSMOKE_Equilibrium *ptEquilibrium;

const double OpenSMOKE_Equilibrium::threshold_element	=	1.e-16;

// ---------------------------------------------------------------------------------------------
// Non linear system solver
// ---------------------------------------------------------------------------------------------
void MyNonLinearSystem_Equilibrium::GetResiduals(BzzVector &lambda, BzzVector &f)
{
	ptEquilibrium->nls(lambda, f);
}

void MyNonLinearSystem_Equilibrium::AssignEquilibrium(OpenSMOKE_Equilibrium *equilibrium)
{
	ptEquilibrium = equilibrium;
}

void MyNonLinearSystem_Equilibrium::ObjectBzzPrint(void)
{
}

void OpenSMOKE_Equilibrium::ErrorMessage(const string message)
{
    cout << endl;
    cout << "Class:  OpenSMOKE_Equilibrium"	<< endl;
    cout << "Object: " << name_object		<< endl;
    cout << "Error:  " << message           << endl;
    cout << "Press enter to continue... "   << endl;
    getchar();
    exit(-1);
}

void OpenSMOKE_Equilibrium::WarningMessage(const string message)
{
    cout << endl;
    cout << "Class:	  OpenSMOKE_Equilibrium"	<< endl;
    cout << "Object:  " << name_object			<< endl;
    cout << "Warning: " << message				<< endl;
	cout << endl;
}

OpenSMOKE_Equilibrium::OpenSMOKE_Equilibrium()
{
	name_object				= "[Name not assigned]";
}

void OpenSMOKE_Equilibrium::SetName(string name)
{
	name_object = name;
}

void OpenSMOKE_Equilibrium::Setup(OpenSMOKE_IdealGas *ideal_gas)
{
	gas = ideal_gas;

	n			= gas->NumberOfSpecies();
	m_Original	= gas->NumberOfElements();

	ChangeDimensions(n, &Z);
	ChangeDimensions(n, &ZFirstGuess);
	ChangeDimensions(n, &q);
	ChangeDimensions(n, &gtilde);
	ChangeDimensions(n, &B);
	ChangeDimensions(n, &N_x_lambda_Reduced);
	ChangeDimensions(m_Original, &alfa_Original);

	NT_Original = gas->elements;
	N_Original  = NT_Original;	
	Transpose(&N_Original);

	BzzMatrix NT_x_N_Original;
	Product(NT_Original, N_Original, &NT_x_N_Original);
	NT_x_N_Factored_Original = NT_x_N_Original;
	
	ChangeDimensions(n, &indexSpecies_Original);
	for(int k=1;k<=n;k++)
		indexSpecies_Original[k] = k;
}

void OpenSMOKE_Equilibrium::Equilibrate(const double _T, const double _P_Pa, BzzVector &_xInitial, BzzVector &xFinal)
{
	// Initial composition
	xInitial = _xInitial;
	TInitial = _T;
	PInitial = _P_Pa;

	// Elemental composition
	Product(NT_Original, xInitial, &alfa_Original);
	
	// Reduced variables: if some elements are zero
	m_Reduced		= m_Original;		// number of elements
	NT_Reduced		= NT_Original;		// matrix of elements
	alfa_Reduced	= alfa_Original;	// matrix of elements

	// Additional reduced variables
	for(int k=m_Original;k>=1;k--)
		if (alfa_Original[k] <= threshold_element)
		{
			m_Reduced--;
			NT_Reduced.DeleteRow(k);
			alfa_Reduced.DeleteElement(k);
		}

	// Additional reduced variables
	indexIncludeSpecies = indexSpecies_Original;
	ChangeDimensions(0, &indexNoSpecies);
	if (m_Original != m_Reduced)
		for(int i=n;i>=1;i--)
		{
			for(int k=m_Original;k>=1;k--)
				if (alfa_Original[k] <= threshold_element)
					if(NT_Original[k][i] != 0.)
					{
						indexNoSpecies.Append(i);
						indexIncludeSpecies.DeleteElement(i);
						break;
					}
		}
	
	// First guess solution preparation
	ChangeDimensions(m_Reduced, &b);
	ChangeDimensions(m_Reduced, &lambdaFirstGuess);
	ChangeDimensions(m_Reduced, &lambdaSolution);

	N_Reduced  = NT_Reduced;	
	Transpose(&N_Reduced);

	BzzMatrix NT_x_N_Reduced;
	Product(NT_Reduced, N_Reduced, &NT_x_N_Reduced);
	NT_x_N_Factored_Reduced = NT_x_N_Reduced;

	// First Guess Solution
	GetGtilde(_T);
	GetFirstGuess();

	// Prov
//	BzzVector xFirstGuess(n);
//	Product(N_Reduced, lambdaFirstGuess, &xFirstGuess);
	//

	// Solution
//	lambdaFirstGuess = 0.;
	Run();
	xFinal = Z;

	int i;

	BzzVector x_elemental(m_Original);
	gas->GetElementalMoleFractions(xInitial, x_elemental);

	Product(N_Reduced, lambdaFirstGuess, &N_x_lambda_Reduced);
	for (i=1; i<=indexIncludeSpecies.Size(); i++)
	{
		int k				= indexIncludeSpecies[i];
		ZFirstGuess[k]	= exp(-gtilde[k] + N_x_lambda_Reduced[k]);
	}
	
	for(i=1;i<=m_Original;i++)
		cout << gas->list_of_elements[i-1] << "\t" << x_elemental[i] << endl;
	
	for(i=1;i<=m_Reduced;i++)
		cout << i << "\t" << lambdaFirstGuess[i] << "\t" << lambdaSolution[i] << endl;

	for(i=1;i<=n;i++)
		cout << gas->names[i] << "\t" << q[i] << "\t" << ZFirstGuess[i] << "\t" << xFinal[i] << endl;

	NT_Original.BzzPrint("NT_Original");
	N_Original.BzzPrint("N_Original");
	NT_x_N_Factored_Original.BzzPrint("NT_x_N Original");
	
	alfa_Original.BzzPrint("ALFA_Original");
	alfa_Reduced.BzzPrint("ALFA_Reduced");
	
	indexIncludeSpecies.BzzPrint("includeSpecies");
	indexNoSpecies.BzzPrint("excludeSpecies");
	
	NT_Reduced.BzzPrint("NT_Reduced");
	N_Reduced.BzzPrint("N_Reduced");
	NT_x_N_Factored_Reduced.BzzPrint("NT_x_N Reduced");
	
	gtilde.BzzPrint("gtilde");
	B.BzzPrint("B");
	b.BzzPrint("bVector");
	lambdaFirstGuess.BzzPrint("lambda first guess");
//	xFirstGuess.BzzPrint("xFirstGuess");
	lambdaSolution.BzzPrint("lambda solution");
	N_x_lambda_Reduced.BzzPrint("Nxlambda");

	BzzVector prov = -gtilde + N_x_lambda_Reduced;
	prov.BzzPrint("exp");
	xFinal.BzzPrint("xFinal");
}

void OpenSMOKE_Equilibrium::GetFirstGuess(BzzVector &x)
{
}

void OpenSMOKE_Equilibrium::GetFirstGuess()
{
	int i;

	// E1
	double g = gtilde.GetSumElements() / double(n);
	for(i=1;i<=n;i++)	q[i] = exp(-gtilde[i]/g);
	double qsum = q.GetSumElements();
	q /= qsum;

	// E2
	q = 1./double(n);

	// E3
//	q = xInitial;
//	Sum(&q,1e-3);
//	double sum = q.GetSumElements();
//	q /= sum;

	for(i=1;i<=n;i++) 
		B[i] = gtilde[i] + log(q[i]);
//		B[i] = 1./double(n);
//		B[i] = exp(gtilde[i] + log_x);

	Product(NT_Reduced, B, &b);
	b.BzzPrint("NT*B");
	Solve(NT_x_N_Factored_Reduced, b, &lambdaFirstGuess);

	BzzFactorizedGauss N_x_NT_Factored_Original;
	BzzMatrix auxx(n,n);
	Product(N_Original, NT_Original, &auxx);
	N_x_NT_Factored_Original = auxx;
	cout << "A" << endl;
	BzzVector NFirstGuess(n);
	BzzVector aux(n);
	Product(N_Reduced,alfa_Reduced, &aux);
	cout << "B" << endl;
	Solve(N_x_NT_Factored_Original, aux, &NFirstGuess);
	cout << "C" << endl;
	double Ntot = NFirstGuess.GetSumElements();

	cout << "Ntot" << endl;
	for(i=1;i<=n;i++)
		cout << i << "  x " << q[i] << "  N " << NFirstGuess[i] << endl;
	getchar();

	for(i=1;i<=100;i++)
	{
		Mode1(Ntot, lambdaFirstGuess, q);
		Mode1(Ntot, lambdaFirstGuess, q);
		Mode1(Ntot, lambdaFirstGuess, q);
		Mode2(Ntot, lambdaFirstGuess, q);
		Mode3(Ntot, lambdaFirstGuess, q);
		Mode4(Ntot, lambdaFirstGuess, q);

		getchar();
	}

	{
		BzzVector h0(gas->NumberOfSpecies());
		BzzVector s0(gas->NumberOfSpecies());
		BzzVector h(gas->NumberOfSpecies());
		BzzVector s(gas->NumberOfSpecies());

		gas->GetStandardEnthalpy_Mole(h0, Constants::T_Reference);		// [J/kmol]
		gas->GetStandardEntropy_Mole(s0, Constants::T_Reference);		// [J/kmol/K]
		gas->GetStandardEnthalpy_Mole(h, TInitial);						// [J/kmol]
		gas->GetStandardEntropy_Mole(s, TInitial);						// [J/kmol/K]
		gas->GetStandardGibbsFreeEnergy_Mole(gtilde, TInitial);			// [J/kmol]
		gtilde /= (Constants::R_J_kmol*TInitial);						// [-]

		for(int j=1;j<=gas->NumberOfSpecies();j++)
			cout	<< gas->names[j]								<< "\t"
					<< gas->M[j]									<< "\t" 
					<< h[j]/1.e6/OpenSMOKE_Conversions::J_from_cal	<< "\t" 
					<< s[j]/1.e3/OpenSMOKE_Conversions::J_from_cal	<< "\t"
					<< gtilde[j]									<< "\t"
					<< endl;

	}

	{/*
		BzzVector x0;
		BzzVector s;
		BzzMatrixDoubleSparse E(NT_Original.Rows(), NT_Original.Columns());
		BzzVector e;

		x0 = xInitial;
		x0[1] = 0.50;;
		x0[2] = 0.50;;
		x0[3] = 0.25;;
		x0[4] = 0.00;;
		s  = gtilde; 
		e  = alfa_Original;
		
		for (int i=1;i<=NT_Original.Rows();i++)
			for (int k=1;k<=NT_Original.Columns();k++)
				if (NT_Original[i][k]!=0) E(i,k) = NT_Original[i][k];
		
		BzzVectorInt iL(xInitial.Size());
		BzzVectorInt iU(xInitial.Size());
		BzzVector xL(xInitial.Size());
		BzzVector xU(xInitial.Size());
		for (int j=1;j<=xInitial.Size();j++)
		{
			iL[j] = j; xL[j] = 0.;
			iU[j] = j; xU[j] = 1e16;
		}
		E.BzzPrint("E");
		e.BzzPrint("e");
		x0.BzzPrint("x0");
		s.BzzPrint("s");
		s.BzzPrint("s");
		iL.BzzPrint("iL");
		iU.BzzPrint("iU");
		xL.BzzPrint("xL");
		xU.BzzPrint("xU");
	
		cout << "Sono 1" << endl;
		BzzLinearProgrammingDouble lp(&x0,&s,&E,&e,0,0,0,0,&iL,&xL,&iU,&xU);
		cout << "Sono 2" << endl;
		lp();
		BzzPause();
		exit(-1);*/
/*
	BzzVector x0(5);
	BzzVector s(5,-10.,-9.,-2.,3.,1.);
	BzzMatrixDoubleSparse E(4,5);
	E(1,1) = 4.; E(1,2) = 2.; // ie = -2 Constraint incompatible
	// the projection of x1 = x2 = 1
	// is y1 = 1.2 y2 = 1.1
	// and 4 + 2 < 7
	
	E(2,1) = 1.; E(2,2) = 2.; E(2,5) = -3.; // ie = -1 Constraint incompatible
	// the projection of x1 .1 x2 = .3 x5 = 1.
	// is y1 = -.307 y2 = -.514 y5 = 1.22
	// and .1 + 2*.3 - 5 > -5
	
	E(3,2) = 2.; E(3,3) = 1.;  E(3,4) = 8.; // ie = 0
	// the projection of x2 .3 x3 = 0. x4 = 1.
	// is y2 = .46 y3 = .08 y4 = 1.

	E(4,1) = 1.; // ie = 0
	BzzVector e(4,7.,-5.,9.,.5);
	BzzVector xL(5);
	BzzVector xU(5);
	xL[1] = .1;
	xL[2] = .3;
	xL[4] = .9;
	xU = 1.;
//	xL[3] = -BZZ_BIG;
//	xU[3] = BZZ_BIG;
	BzzVectorInt iL(5,1,2,3,4,5);
	BzzVectorInt iU(5,1,2,3,4,5);
	
	BzzLinearProgrammingDouble lp(&x0,&s,&E,&e,0,0,0,0,&iL,&xL,&iU,&xU);
	lp();
	lp.BzzPrint("Results");
	BzzPause();*/
	}

	{
		int i,j;
		int icase;
#if LINUX_SO == 0		
		int N = NT_Original.Columns();
		int M = NT_Original.Rows();
		Mat_IO_DP a(M+2, N+1);
		Vec_O_INT izrov(N);
		Vec_O_INT iposv(M);
		int m1 = 0;
		int m2 = 0;
		int m3 = M;


		for (i=2;i<=N+1;i++)	a[0][i-1] = -gtilde[i-1];
		for (i=2;i<=M+1;i++)	a[i-1][0] = alfa_Original[i-1];
		
		for (j=2;j<=M+1;j++)
			for (i=2;i<=N+1;i++)
				a[j-1][i-1] = -NT_Original[j-1][i-1];


		//cout << "AStarting " << endl;
		//for (j=1;j<=M+2;j++)
		//{
		//	for (i=1;i<=N+1;i++)
		//		cout << a[j-1][i-1] << "\t";
		//	cout << endl;
		//}
/*
		int N = 4;
		int M = 4;
		Mat_IO_DP a(M+2, N+1);	
		int m1 = 2;
		int m2 = 1;
		int m3 = 1;
		int icase;
		Vec_O_INT izrov(N+1);
		Vec_O_INT iposv(M+1);
		
		cout << "Matrix..." << endl;
		i = 0; j = 0;
		
		a[i][j++] = 0;		a[i][j++] = 1;	a[i][j++] = 1;	a[i][j++] = 3;	a[i][j++] = -0.50;	cout << "Matrix line 1 done..." << endl;	i++; j=0;		
		a[i][j++] = 740;	a[i][j++] = -1;	a[i][j++] = 0;	a[i][j++] = -2;	a[i][j++] = 0;		cout << "Matrix line 2 done..." << endl;	i++; j=0;	
		a[i][j++] = 0;		a[i][j++] = 0;	a[i][j++] = -2;	a[i][j++] = 0;	a[i][j++] = 7;		cout << "Matrix line 3 done..." << endl;	i++; j=0;	
		a[i][j++] = 0.5;	a[i][j++] = 0;	a[i][j++] = -1;	a[i][j++] = 1;	a[i][j++] = -2;		cout << "Matrix line 4 done..." << endl;	i++; j=0;	
		a[i][j++] = 9;		a[i][j++] = -1;	a[i][j++] = -1;	a[i][j++] = -1;	a[i][j++] = -1;		cout << "Matrix line 5 done..." << endl;	i++; j=0;	
		a[i][j++] = 0;		a[i][j++] = 0;	a[i][j++] = 0;	a[i][j++] = 0;	a[i][j++] = 0;		cout << "Matrix line 6 done..." << endl;	i++; j=0;	
		cout << "End Matrix..." << endl;
*/

		NR::simplx(a, m1, m2, m3, icase, izrov, iposv);

		cout << "Status " << icase << endl;
//		cout << "A " << endl;
//		for (j=1;j<=M+2;j++)
//		{
//			for (i=1;i<=N+1;i++)
//				cout << a[j-1][i-1] << "\t";
//			cout << endl;
//		}
//		cout << "izrov " << endl;
//		for (j=1;j<=N;j++)
//			cout << izrov[j-1] << "\t";
//		cout << endl << endl;
//		cout << "iposv " << endl;
//		for (j=1;j<=M;j++)
//			cout << iposv[j-1] << "\t";
//
//		cout << "End simplex..." << endl;
//
		// Solution
		BzzVector xStarting(N);
		for(j=0;j<=M-1;j++)	if (iposv[j] < N)	xStarting[iposv[j]+1] = a[j+1][0];
		for(j=0;j<=N-1;j++)	if (izrov[j] < N)	xStarting[izrov[j]+1] = 0.;
		double sum = xStarting.GetSumElements();
		xStarting /= sum;

		cout << "Initial: " << Dot(xInitial, gtilde)	<< endl;
		cout << "Final:   " << -a[0][0]					<< endl;
		for(j=1;j<=N;j++)
			cout << "x" << j << "\t= " << xInitial[j] << "\t" << xStarting[j] << "\t" << gas->names[j] <<  endl;
#endif

	}
}

void OpenSMOKE_Equilibrium::GetGtilde(const double T)
{
	gas->GetStandardGibbsFreeEnergy_Mass(gtilde, T);
	gtilde /= (Constants::R_J_kmol*T);
	ElementByElementProduct(gtilde, gas->M, &gtilde);
}

void OpenSMOKE_Equilibrium::Run()
{
	double startTime = BzzGetCpuTime();

	nls_object.AssignEquilibrium(this);
	BzzNonLinearSystemObject o(lambdaFirstGuess, &nls_object);
	o();
	o.GetSolution(&lambdaSolution, &residuals);

	double endTime = BzzGetCpuTime();

	double maxResidual	= residuals.MaxAbs();
	double meanResidual = 0.;
	for (int i=1;i<=m_Reduced;i++)
		meanResidual += fabs(residuals[i]);
	meanResidual /= m_Reduced;

	cout << "Solution time: " << endTime - startTime	<< " s"	<< endl;
	cout << "Mean residual: " << meanResidual			<< endl;
	cout << "Max. residual: " << maxResidual			<< endl;
}

void OpenSMOKE_Equilibrium::nls(BzzVector &lambda, BzzVector &f)
{
	int i, k;

	Product(N_Reduced, lambda, &N_x_lambda_Reduced);
	for (i=1; i<=indexIncludeSpecies.Size(); i++)
	{
		k = indexIncludeSpecies[i];
		Z[k] = exp(-gtilde[k] + N_x_lambda_Reduced[k]);
	}

	for (i=1;i<=indexNoSpecies.Size();i++)
		Z[indexNoSpecies[i]] = 0.;

	BzzVector aux = NT_Reduced.GetRow(m_Reduced);
	double coefficient = Dot( aux, Z);
	for(k=1;k<=m_Reduced-1;k++)
	{
		BzzVector aux = NT_Reduced.GetRow(k);
		f[k] = alfa_Reduced[m_Reduced]*Dot(aux, Z) - alfa_Reduced[k]*coefficient;
	}

	f[m_Reduced] = Z.GetSumElements() - 1.0;
}

double OpenSMOKE_Equilibrium::AdiabaticTemperature(const double Tinitial, BzzVector &xInitial)
{
/*	double Hfixed;
	double Hadiabatic;
	double Tadiabatic;
	double MW;

	BzzVector omegaInitial(n);
	BzzVector omegaAdiabatic(n);
	BzzVector xAdiabatic(n);

	gas->GetMWAndMassFractionsFromMoleFractions(MW, omegaInitial, xInitial);		// [kg/kmol]
	Hfixed = gas->GetMixEnthalpy_Mass(Tinitial, omegaInitial);		// [J/kg]
	

	// loop to estimate T
	double T0;
	{
		double Tmin = 200.;
		double Tmax = 4000.;
		double slope, dT;
		double Hhigh, Hlow, Hval;
		

		Hhigh	= gas->GetMixEnthalpy_Mass(Tmax, omegaInitial);		// [J/kg]
		Hlow	= gas->GetMixEnthalpy_Mass(Tmin, omegaInitial);		// [J/kg]

		// start with T at the midpoint of the range
		T0 = 0.505*(Tmin + Tmax);
		
		// loop up to 5 times
		for (int it = 0; it < 5; it++) 
		{
			// set the composition and get p1
			Hval	= gas->GetMixEnthalpy_Mass(T0, omegaInitial);		// [J/kg]

	        // If this value of p1 is greater than the specified
	        // property value, then the current temperature is too
	        // high. Use it as the new upper bound. Otherwise, it
	        // is too low, so use it as the new lower bound.
	        if (Hval > Hfixed) 
	        { 
	            Tmax = T0; 
		        Hhigh = Hval; 
			}
			else 
			{
				Tmin = T0; 
				Hlow = Hval; 
			}

			// Determine the new T estimate by linearly intepolation
			// between the upper and lower bounds
			slope = (Hhigh - Hlow)/(Tmax - Tmin);
			dT = (Hfixed - Hlow)/slope;

			// If within 100 K, terminate the search
			if (fabs(dT) < 100.0) break;

			// update the T estimate
			T0 = Tmin + dT;
		}
	}

	Tadiabatic = T0;
	cout << "TFirstGuess: " << Tadiabatic << endl;
	
	xAdiabatic = xInitial;
	for(int i=1;i<=64;i++)
	{
		cout << "it: " << i << endl;
		cout << "TA " << Tadiabatic  << endl;
		Equilibrate(Tadiabatic, PInitial, xInitial, xAdiabatic);
		
	//	for(int j=1;j<=n;j++)
	//		cout << gas->names[j] << "\t" << xInitial[j] << "\t" << xAdiabatic[j] << endl;

		gas->GetMWAndMassFractionsFromMoleFractions(MW, omegaAdiabatic, xAdiabatic);	// [kg/kmol]
		
		
		getchar();

		exit(-1);

		Hadiabatic = gas->GetMixEnthalpy_Mass(Tadiabatic, omegaAdiabatic);			// [J/kg]

		cout << "Hf " << Hfixed  << endl;
		cout << "HA " << Hadiabatic  << endl;

		if (fabs(Hadiabatic-Hfixed)/max(fabs(Hfixed), fabs(Hadiabatic)) <= 1.e-5)
			return Tadiabatic;


		double deltaT = gas->GetTemperatureFromMixEnthalpy(Hfixed, MW, xAdiabatic) - Tadiabatic;
		if (deltaT >  100.)	deltaT = 100.;
		if (deltaT < -100.)	deltaT = -100.;
		Tadiabatic += deltaT;
		cout << "TA " << Tadiabatic  << endl;
	}
*/
	return -1;
}


void OpenSMOKE_Equilibrium::Mode1(const double Ntot, BzzVector &lambda, BzzVector &x)
{
	int k,i,j;

	int NE = NT_Reduced.Rows();
	int NS = NT_Reduced.Columns();
	
	double Beta;
	double deltas;
	double W;

	BzzVector	H(NE);
	BzzVector	f(NE);
	BzzVector	deltalambda(NE);
	BzzMatrix	Q(NE,NE);

	cout << "Mode 1" << endl;
	cout << "--------------------------------------" << endl;	
	for(j=1;j<=NE;j++)
		cout << "lambda" << j << "_Before " << lambda[j] << endl;

	for(j=1;j<=NS;j++)
		cout << "x" << j << "_Before " << x[j] << endl;

	// W
	W = Ntot*(x.GetSumElements()-1) - Dot(lambda, alfa_Reduced);
	cout << "W: " << W << endl;

	// Vector H
	Product(NT_Reduced, x, &H);
	H *= Ntot;
	H -= alfa_Reduced;

	for(j=1;j<=NE;j++)
	{
		cout << "H " << j << H[j]			<< endl;
		cout << "p " << j << alfa_Reduced[j] << endl;
	}

	// Q
	for(i=1;i<=NE;i++)
		for(k=1;k<=NE;k++)
		{
			Q[i][k] = 0.;
			for(j=1;j<=NS;j++)
				Q[i][k] += NT_Reduced[i][j]*NT_Reduced[k][j]*x[j]; 
		}

	// Beta
	Beta = -sqrt(Dot(H,H));
	cout << "Beta: " << Beta << endl;

	// f
	for(k=1;k<=NE;k++)
		f[k] = H[k]/Beta;

	// DeltaS
	double sum=0.;
	for(k=1;k<=NE;k++)
		for(i=1;i<=NE;i++)
			sum += Q[i][k]*f[i]*f[k];
	deltas = -Beta/sum;

	cout << "deltas: " << deltas << endl;

	// lambda
	for(j=1;j<=NE;j++)
		lambda[j] += f[j]*deltas;

	// x
	
	for(j=1;j<=NS;j++)
	{
		BzzVector aux = N_Reduced.GetRow(j);
		x[j] = exp(-gtilde[j] + Dot( aux,lambda) );
	}

	for(j=1;j<=NE;j++)
		cout << "lambda" << j << "_After " << lambda[j] << endl;

	for(j=1;j<=NS;j++)
		cout << "x" << j << "_After " << x[j] << endl;

	W = Ntot*(x.GetSumElements()-1) - Dot(lambda, alfa_Reduced);
	cout << "W: " << W << endl;
	cout << "N: " << Ntot << endl;
	cout << "Z: " << x.GetSumElements() << endl;
}

void OpenSMOKE_Equilibrium::Mode2(const double Ntot, BzzVector &lambda, BzzVector &x)
{
	int k,i,j;

	int NE = NT_Reduced.Rows();
	int NS = NT_Reduced.Columns();
	
	double W;

	BzzVector	H(NE);
	BzzVector	deltalambda(NE);
	BzzMatrix	Q(NE,NE);

	cout << "Mode 2" << endl;
	cout << "--------------------------------------" << endl;

	for(j=1;j<=NE;j++)
		cout << "lambda" << j << "_Before " << lambda[j] << endl;

	for(j=1;j<=NS;j++)
		cout << "x" << j << "_Before " << x[j] << endl;

	// W
	W = Ntot*(x.GetSumElements()-1) - Dot(lambda, alfa_Reduced);
	cout << "W: " << W << endl;

	// Vector H
	Product(NT_Reduced, x, &H);
	H *= Ntot;
	H -= alfa_Reduced;

	for(j=1;j<=NE;j++)
	{
		cout << "H" << j << H[j]			<< endl;
		cout << "p" << j << alfa_Reduced[j] << endl;
	}

	// Q
	for(i=1;i<=NE;i++)
		for(k=1;k<=NE;k++)
		{
			Q[i][k] = 0.;
			for(j=1;j<=NS;j++)
				Q[i][k] += NT_Reduced[i][j]*NT_Reduced[k][j]*x[j]; 
		}

	BzzFactorizedGauss Q_Factored;
	Q_Factored = Q;
	Solve(Q_Factored, -H, &deltalambda);

	// lambda
	for(j=1;j<=NE;j++)
		lambda[j] += deltalambda[j];

	// x
	for(j=1;j<=NS;j++)
	{
		BzzVector aux = N_Reduced.GetRow(j);
		x[j] = exp(-gtilde[j] + Dot( aux,lambda) );
	}

	for(j=1;j<=NE;j++)
		cout << "lambda" << j << "_After " << lambda[j] << endl;

	for(j=1;j<=NS;j++)
		cout << "x" << j << "_After " << x[j] << endl;

	W = Ntot*(x.GetSumElements()-1) - Dot(lambda, alfa_Reduced);
	cout << "W: " << W << endl;
	cout << "N: " << Ntot << endl;
	cout << "Z: " << x.GetSumElements() << endl;
}

void OpenSMOKE_Equilibrium::Mode3(double &Ntot, BzzVector &lambda, BzzVector &x)
{
	int k,i,j;

	int NE = NT_Reduced.Rows();
	int NS = NT_Reduced.Columns();
	int NF = 1;

	double alfa;
	double deltas;
	double W;
	double V;
	double r;


	BzzVector	D(NE);
	BzzMatrix	Q(NE,NE);
	BzzMatrix	E(NE,NF);
	BzzMatrix	A(NF,NF);
	BzzVector	deltalambda(NE);

	cout << "Mode 3" << endl;
	cout << "--------------------------------------" << endl;	
	for(j=1;j<=NE;j++)
		cout << "lambda" << j << "_Before " << lambda[j] << endl;

	for(j=1;j<=NS;j++)
		cout << "x" << j << "_Before " << x[j] << endl;

	// W
	W = Ntot*(x.GetSumElements()-1) - Dot(lambda, alfa_Reduced);
	cout << "W: "		<< W	<< endl;
	cout << "Ntot: "	<< Ntot << endl;

	// Vector H
	Product(NT_Reduced, x, &D);
	for(j=1;j<=NE;j++)
		cout << "D " << j << D[j]			<< endl;

	// V
	V = x.GetSumElements()-1.;

	// Q
	for(i=1;i<=NE;i++)
		for(k=1;k<=NE;k++)
		{
			Q[i][k] = 0.;
			for(j=1;j<=NS;j++)
				Q[i][k] += NT_Reduced[i][j]*NT_Reduced[k][j]*x[j]; 
		}

	// E
	BzzFactorizedGauss Q_Factored;
	Q_Factored = Q;
	Solve(Q_Factored, -D, &E);

	// A
	A[1][1] = 0.;
	for(k=1;k<=NE;k++)
		A[1][1] += D[k]*E[k][1]; 

	// alfa
	alfa = -sqrt(V*V);
	cout << "alfa: " << alfa << endl;

	// r
	r = V/alfa;

	// DeltaS
	deltas = -alfa/(A[1][1]*r*r);
	cout << "deltas: " << deltas << endl;

	// lambda
	Ntot += deltas;

	// x
	for(j=1;j<=NS;j++)
	{
		BzzVector aux = N_Reduced.GetRow(j);
		x[j] = exp(-gtilde[j] + Dot( aux,lambda) );
	}

	for(j=1;j<=NE;j++)
		cout << "lambda" << j << "_After " << lambda[j] << endl;

	for(j=1;j<=NS;j++)
		cout << "x" << j << "_After " << x[j] << endl;

	W = Ntot*(x.GetSumElements()-1) - Dot(lambda, alfa_Reduced);
	cout << "W: "		<< W	<< endl;
	cout << "Ntot: "	<< Ntot << endl;
	cout << "Z: "		<< x.GetSumElements() << endl;
}

void OpenSMOKE_Equilibrium::Mode4(double &Ntot, BzzVector &lambda, BzzVector &x)
{
	int k,i,j;

	int NE = NT_Reduced.Rows();
	int NS = NT_Reduced.Columns();
	
	double W;

	BzzVector	H(NE);
	BzzVector	D(NE);
	BzzVector	delta(NE+1);
	BzzMatrix	Q(NE,NE);

	cout << "Mode 4" << endl;
	cout << "--------------------------------------" << endl;

	// W
	W = Ntot*(x.GetSumElements()-1) - Dot(lambda, alfa_Reduced);

	// Vector H
	Product(NT_Reduced, x, &H);
	H *= Ntot;
	H -= alfa_Reduced;

	// Vector D
	Product(NT_Reduced, x, &D);

	// Q
	for(i=1;i<=NE;i++)
		for(k=1;k<=NE;k++)
		{
			Q[i][k] = 0.;
			for(j=1;j<=NS;j++)
				Q[i][k] += NT_Reduced[i][j]*NT_Reduced[k][j]*x[j]; 
		}

	BzzFactorizedGauss Q_Factored(NE+1, NE+1);
	for(i=1;i<=NE;i++)
		for(k=1;k<=NE;k++)
			Q_Factored[i][k] = Q[i][k];
	for(i=1;i<=NE;i++)
	{
		Q_Factored[NE+1][i] = D[i];
		Q_Factored[i][NE+1] = D[i];
	}
	BzzVector b(NE+1);
	for(i=1;i<=NE;i++)
		b[i] = -H[i];
	b[NE+1] = 1.-x.GetSumElements();

	Solve(Q_Factored, b, &delta);

	// lambda
	for(j=1;j<=NE;j++)
		lambda[j] += delta[j];
	Ntot += delta[NE+1];

	// x
	for(j=1;j<=NS;j++)
	{
		BzzVector aux = N_Reduced.GetRow(j);
		x[j] = exp(-gtilde[j] + Dot( aux,lambda) );
	}
	W = Ntot*(x.GetSumElements()-1.) - Dot(lambda, alfa_Reduced);

	for(j=1;j<=NE;j++)
		cout << "lambda"	<< j << "\t" << lambda[j]	<< endl;

	for(j=1;j<=NS;j++)
		cout << "x"			<< j << "\t" << x[j]		<< endl;

	for(j=1;j<=NE;j++)
		cout << "H"			<< j << "\t" << H[j]		<< endl;

	
	cout << "W:\t" << W << endl;
	cout << "N:\t" << Ntot << endl;
	cout << "Z:\t" << x.GetSumElements() << endl;
}
