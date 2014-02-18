/***************************************************************************
 *   Copyright (C) 2006 by Alberto Cuoci   	   *
 *   alberto.cuoci@polimi.it   						   *
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

#include "engine/OpenSMOKE_ReactingGas.h"
#include "engine/OpenSMOKE_Kinetics.h"
#include "basic/OpenSMOKE_Constants.h"
#include "basic/OpenSMOKE_Utilities.h"
#include <vector>
#include <iomanip>


const int		OpenSMOKE_Kinetics::NUM_MAX					= 40000;
const int		OpenSMOKE_Kinetics::NUM_MAX_CHEBISHEV		= 10000;
const int		OpenSMOKE_Kinetics::NUM_MAX_LOGPRESSURE		= 10000;
const double	OpenSMOKE_Kinetics::S1						= 1.001;	
const double	OpenSMOKE_Kinetics::S2						= S1 * S1;	
const double	OpenSMOKE_Kinetics::S3						= S2 * S1;	
const double	OpenSMOKE_Kinetics::SSQ						= sqrt(S1);	

void OpenSMOKE_Kinetics::ErrorMessage(const string message)
{
    cout << endl;
    cout << "Class:  OpenSMOKE_Kinetics"		<< endl;
    cout << "Object: " << name_object			<< endl;
    cout << "Error:  " << message				<< endl;
    cout << "Press enter to continue... "		<< endl;
    getchar();
    exit(-1);
}

void OpenSMOKE_Kinetics::WarningMessage(const string message)
{
    cout << endl;
    cout << "Class:	  OpenSMOKE_Kinetics"		<< endl;
    cout << "Object:  " << name_object			<< endl;
    cout << "Warning: " << message				<< endl;
	cout << endl;
}

OpenSMOKE_Kinetics::OpenSMOKE_Kinetics()
{
	name_object	= "[Name not assigned]";
}

void OpenSMOKE_Kinetics::SetName(const string name)
{
	name_object = name;
}


// Costruzione della matrice GAMMA (NC x NR)
BzzMatrix OpenSMOKE_Kinetics::constructGamma(const string fileSt)
{
	int i, j;
	double val;
	int numComponents, numReactions;
	
	cout << "Gamma Matrix construction... ";
 
	ifstream inputFile;
	openInputFileAndControl(inputFile, fileSt);

	inputFile >> numComponents >> numReactions;

	BzzMatrix gamma;
	ChangeDimensions(numComponents, numReactions, &gamma);

	while(1)
	{
		inputFile >> i >> j >> val;
	
		if (inputFile.eof()!=0) break;

		gamma[i][j] += val;
	}

	inputFile.close();

	cout << " DONE" << endl;

	return gamma;
}

//----------------------------------------------------------------------------------//
//																					//
//						   computeDirectAndInverse									//
//																					//
//----------------------------------------------------------------------------------//	

// Calcolo delle velocita di reazione dirette e inverse (in realt?si 
// tratta soltanto della produttoria delle concentrazioni, non ancora 
// moltiplicata per alcuna costante cinetica)

void OpenSMOKE_Kinetics::ComputeDirectAndInverse(BzzVector& c, BzzVector &rDirC, BzzVector &rInvC)
{
	rDirC =  1.;
	rInvC = -1.;

//	if (numGlobal != 0)	
//		GlobalReactionRates(c, rDirC, rInvC);

	if (numGlobal != NR)
	{
		int i,k;
		double c1,c2,c3,csq;

		int *jD1, *jD2, *jD3, *jD4, *jD5;
		double *vD5;

		int *jIE1, *jIE2, *jIE3, *jIE4, *jIE5;
		double *vIE5;

		jD1 = lambda_jDir1.GetHandle();
		jD2 = lambda_jDir2.GetHandle();
		jD3 = lambda_jDir3.GetHandle();
		jD4 = lambda_jDir4.GetHandle();
		jD5 = lambda_jDir5.GetHandle();
		vD5 = lambda_valDir5.GetHandle();

		jIE1 = lambda_jInvEq1.GetHandle();
		jIE2 = lambda_jInvEq2.GetHandle();
		jIE3 = lambda_jInvEq3.GetHandle();
		jIE4 = lambda_jInvEq4.GetHandle();
		jIE5 = lambda_jInvEq5.GetHandle();
		vIE5 = lambda_valInvEq5.GetHandle();

		
		for(i = 1;i <= NC;i++)
		{
			c1 = c[i];
			c2 = c1 * c1;
			c3 = c2 * c1;
			if(lambda_numDir4[i] != 0)
				csq = sqrt(c1);
			for(k = 1;k <= lambda_numDir1[i];k++)
			{
				rDirC[*jD1] *= c1;
				jD1++;
			}
			for(k = 1;k <= lambda_numDir2[i];k++)
			{
				rDirC[*jD2] *= c2;
				jD2++;
			}
			for(k = 1;k <= lambda_numDir3[i];k++)
			{
				rDirC[*jD3] *= c3;
				jD3++;
			}
			for(k = 1;k <= lambda_numDir4[i];k++)
			{
				rDirC[*jD4] *= csq;
				jD4++;
			}
			for(k = 1;k <= lambda_numDir5[i];k++)
			{
				rDirC[*jD5] *= pow(c1,*vD5);
				jD5++;
				vD5++;
			}

			for(k = 1;k <= lambda_numInvEq1[i];k++)
			{
				rInvC[*jIE1] *= c1;
				jIE1++;
			}
			for(k = 1;k <= lambda_numInvEq2[i];k++)
			{
				rInvC[*jIE2] *= c2;
				jIE2++;
			}
			for(k = 1;k <= lambda_numInvEq3[i];k++)
			{
				rInvC[*jIE3] *= c3;
				jIE3++;
			}
			for(k = 1;k <= lambda_numInvEq4[i];k++)
			{
				rInvC[*jIE4] *= csq;
				jIE4++;
			}
			for(k = 1;k <= lambda_numInvEq5[i];k++)
			{
				rInvC[*jIE5] *= pow(c1,*vIE5);
				jIE5++;
				vIE5++;
			}
		}
	}

	// Negative frequency factors
	for(int j=1;j<=negativeSigns.Size();j++)
	{
		rDirC[negativeSigns[j]] *= -1.;
		rInvC[negativeSigns[j]] *= -1.;
	}

	if (numGlobal != 0)	
		GlobalReactionRates(c, rDirC, rInvC);
}

void OpenSMOKE_Kinetics::GlobalReactionRates(BzzVector& c, BzzVector &rDirC, BzzVector &rInvC)
{
	int i;
	
	double Cstar = 1.e-8;
	double ALFA  = 1.e-5;
	double H = 1.50*log(ALFA/(1.-ALFA));
	double K = 2.00*log((1.-ALFA)/ALFA) / Cstar;
	double delta = 1.e9;

	// Direct reactions
	for(i=1;i<=numGlobal;i++)
	{
		rDirC[iGlobal[i]] = 1.;

		for(int j=1;j<=iGlobalDirect[i].Size();j++)
		{
			double C		= c[iGlobalDirect[i][j]];
			double lambda	= lambdaGlobalDirect[i][j];

			if (lambda>=1.)
			{
				rDirC[iGlobal[i]] *= pow(C, lambda);
			}
			else
			{
				double m = (tanh(K*C + H)+1.)/2.;	//con questi coefficienti inizia a decrescere per c<1.5e-5 e finisce per 9.e-6
				double gamma =	m*pow(C+m/delta,lambda) + 
								(1-m)*pow(Cstar,lambda-1.)*C;
				rDirC[iGlobal[i]] *= gamma;
			}
		}
	}

	// Inverse reactions
	for(i=1;i<=numGlobal;i++)
	{
		rInvC[iGlobal[i]] = -1.;

		for(int j=1;j<=iGlobalInverse[i].Size();j++)
		{
			double C		= c[iGlobalInverse[i][j]];
			double lambda	= lambdaGlobalInverse[i][j];
/*
			cout << "  j:      " << j << endl;
			cout << "  r:      " << iGlobal[i] << endl;
			cout << "  ind:    " << iGlobalInverse[i][j] << endl;
			cout << "  C:      " << c[iGlobalInverse[i][j]] << endl;
			cout << "  lambda: " << lambdaGlobalInverse[i][j] << endl;
			cout << "  rInvI:    " << rInvC[iGlobal[i]] << endl;
*/
			if (lambda>=1.)
			{
				rInvC[iGlobal[i]] *= pow(C, lambda);
			}
			else
			{
				double m = (tanh(K*C + H)+1.)/2.;	//con questi coefficienti inizia a decrescere per c<1.5e-5 e finisce per 9.e-6
				double gamma =	m*pow(C+m/delta,lambda) + 
								(1-m)*pow(Cstar,lambda-1.)*C;
				rInvC[iGlobal[i]] *= gamma;
			}
		}
	}
}

//----------------------------------------------------------------------------------//
//						          computeFallOff									//
//----------------------------------------------------------------------------------//	
void OpenSMOKE_Kinetics::FallOff(double T, BzzVector &coeffM, BzzVector &k1, BzzVector &k2, BzzVector &logFcent, BzzVector &coeffFallOff)
{
	int j,k;
	double Pr, nTroe,cTroe,sTroe,logPr,xSRI, wF;

	coeffFallOff = k1;
	
	for(k = 1;k <= numFallOff;k++)
	{
		j=iFallOff[k];

		Pr = k1[j] * coeffM[j] / k2[j];
		
		switch(jThirdBody[j])
		{
			case 2: // Lindeman
				wF = 1.;
				break;
			case 3: // Troe
				logPr = log10(Pr);
				nTroe = 0.75 - 1.27 * logFcent[j];
				cTroe = -0.4 - 0.67 * logFcent[j];
				sTroe = logPr + cTroe;
				wF = pow(10.,logFcent[j] /(1. + BzzPow2(sTroe / (nTroe - 0.14 * sTroe))));
				break;
			case 4:	// SRI
				logPr = log10(Pr);
				xSRI = 1. / (1. + BzzPow2(logPr));
				wF = pow(logFcent[j],xSRI) * dPressure[j];
				if(ePressure[j] != 0.)
					wF *= pow(T,ePressure[j]);
				break;
			default: // Lindeman
				wF = 1.;
				break;
			}
	
		coeffFallOff[j] = k2[j] * (Pr / (1. + Pr)) * wF;
	}
}

void OpenSMOKE_Kinetics::ChemicallyActivatedBimolecularReactions(double T, BzzVector &coeffM, BzzVector &k1, BzzVector &k2, BzzVector &logFcent, BzzVector &coeffFallOff)
{
	int j,k;
	double Pr, nTroe,cTroe,sTroe,logPr,xSRI, wF;

	// coeffFallOff = k1; Viene gia' fatto in FallOff(...)

	for(k = 1;k <= numCABR;k++)
	{
		j=iCABR[k];

		Pr = k1[j] * coeffM[j] / k2[j];
		
		switch(jThirdBody[j])
		{
			case 5: // Lindeman
				wF = 1.;
				break;
			case 6: // Troe
				logPr = log10(Pr);
				nTroe = 0.75 - 1.27 * logFcent[j];
				cTroe = -0.4 - 0.67 * logFcent[j];
				sTroe = logPr + cTroe;
				wF = pow(10.,logFcent[j] /(1. + BzzPow2(sTroe / (nTroe - 0.14 * sTroe))));
				break;
			case 7:
				logPr = log10(Pr);
				xSRI = 1. / (1. + BzzPow2(logPr));
				wF = pow(logFcent[j],xSRI) * dPressure[j];
				if(ePressure[j] != 0.)
					wF *= pow(T,ePressure[j]);
				break;
			default: // Lindeman
				wF = 1.;
				break;
			}
	
		coeffFallOff[j] = k1[j] /(1. + Pr) * wF;
	}
}

//----------------------------------------------------------------------------------//
//						          computeThirdBody									//
//----------------------------------------------------------------------------------//	

// Calcola il coefficiente correttivo Gamma per tutte quelle reazioni che prevedono
// la presenza di un terzo corpo (anche quelle di FallOff)
void OpenSMOKE_Kinetics::ThirdBody(BzzVector &c, double cTot, BzzVector &coeffM)
{	
	int j, k;

	for(int s = 1;s <= numThirdBody;s++)
		{
			j=iThirdBody[s];

			if(jNumEfficiency[j] != -1) // Reazioni con terzo corpo "Convenzionali" 
			{
				coeffM[j] = cTot; // questa espressione corrisponde a quella 
							      // che si avrebbe se tutte le efficienze fossero unitarie

				// quel termine viene corretto se esistono delle efficienze con valori
				// diversi da quelli unitari (il -1 che dovrebbe comparire ?gi?stato sottratto)
				if(jNumEfficiency[j] != 0) 
					for(k = 1;k <= jNumEfficiency[j];k++)
					coeffM[j] += c[iEfficiency[j][k]] * efficiency[j][k];
			}

			else // Reazioni con terzo corpo "Non Convenzionali"
			coeffM[j] = c[iEfficiency[j][1]] * efficiency[j][1];
		}
}

//----------------------------------------------------------------------------------//
//																					//
//						   reactionsWithEquilibrium									//
//																					//
//----------------------------------------------------------------------------------//	

void OpenSMOKE_Kinetics::reactionsWithEquilibrium(BzzVector &rDirC, BzzVector &rInvC,BzzVector &uKeq,BzzVector &r)
{
	int j,k;

	r=rDirC;

	for(k = 1;k <= numEquilibrium;k++)
	{
		j = reactionWithEquil[k];
		r[j] += rInvC[j] * uKeq[j];
	}
}

void OpenSMOKE_Kinetics::reactionsWithEquilibriumForwardAndBackward(BzzVector &rBackward,BzzVector &uKeq)
{
	int j,k;
	for(k = 1;k <= numEquilibrium;k++)
	{
		j = reactionWithEquil[k];
		rBackward[j] *= uKeq[j];
	}
	for(k=1;k<=NR;k++)
	{
		if (IsAReversibleReaction(k) == 0)
			rBackward[k] = 0.;
	}
}


//----------------------------------------------------------------------------------//
//																					//
//						   reactionsWithThirdBody									//
//																					//
//----------------------------------------------------------------------------------//	

void OpenSMOKE_Kinetics::reactionsWithThirdBody(BzzVector &coeffM,BzzVector &r)
{
	int j,k;

	for(k = 1;k <= numThirdBodyOnly;k++)
	{
		j=iThirdBodyOnly[k];
		r[j] *= coeffM[j];
	}	
}

void OpenSMOKE_Kinetics::reactionsWithThirdBodyForwardAndBackward(BzzVector &coeffM, BzzVector &rForward, BzzVector &rBackward)
{
	int j,k;

	for(k = 1;k <= numThirdBodyOnly;k++)
	{
		j=iThirdBodyOnly[k];
		rForward[j] *= coeffM[j];
		rBackward[j] *= coeffM[j];
	}	
}

//----------------------------------------------------------------------------------//
//							 compositionReactionRates								//
//----------------------------------------------------------------------------------//	

void OpenSMOKE_Kinetics::compositionReactionRates(BzzVector &r, BzzVector *R)
{
	int i,k;

	int* jD1 = jDir1.GetHandle();
	int* jD2 = jDir2.GetHandle();
	int* jD3 = jDir3.GetHandle();
	int* jD4 = jDir4.GetHandle();
	int* jD5 = jDir5.GetHandle();
	double* vD5 = valDir5.GetHandle();

	int* jIT1 = jInvTot1.GetHandle();
	int* jIT2 = jInvTot2.GetHandle();
	int* jIT3 = jInvTot3.GetHandle();
	int* jIT4 = jInvTot4.GetHandle();
	int* jIT5 = jInvTot5.GetHandle();
	double* vIT5 = valInvTot5.GetHandle();
	
	*R = 0.;

	for(i = 1;i <= NC;i++)
		{
		for(k = 1;k <= numDir1[i];k++)
			(*R)[i] -= r[*jD1++];
		for(k = 1;k <= numDir2[i];k++)
			(*R)[i] -= (r[*jD2] + r[*jD2++]);
		for(k = 1;k <= numDir3[i];k++)
			(*R)[i] -= (r[*jD3] + r[*jD3] + r[*jD3++]);
		for(k = 1;k <= numDir4[i];k++)
			(*R)[i] -= .5 * r[*jD4++];
		for(k = 1;k <= numDir5[i];k++)
			(*R)[i] -= (*vD5++) * r[*jD5++];

		for(k = 1;k <= numInvTot1[i];k++)
			(*R)[i] += r[*jIT1++];
		for(k = 1;k <= numInvTot2[i];k++)
			(*R)[i] += (r[*jIT2] + r[*jIT2++]);
		for(k = 1;k <= numInvTot3[i];k++)
			(*R)[i] += (r[*jIT3] + r[*jIT3] + r[*jIT3++]);
		for(k = 1;k <= numInvTot4[i];k++)
			(*R)[i] += .5 * r[*jIT4++];
		for(k = 1;k <= numInvTot5[i];k++)
			(*R)[i] += (*vIT5++) * r[*jIT5++];
		}
}

//----------------------------------------------------------------------------------//
//																					//
//						   computeKineticParameters									//
//																					//
//----------------------------------------------------------------------------------//	

void OpenSMOKE_Kinetics::ComputeKineticParameters(BzzVector &reactionDH, BzzVector &reactionDS, double T, double P_Pa, double logT, double uT, double loguRT, BzzVector &uKeq, BzzVector &k1, BzzVector &k2, BzzVector &logFcent)
{
	int k,j;

	for(k = 1;k <= numEquilibrium;k++)
	{
		j = reactionWithEquil[k];
		uKeq[j] = exp(-reactionDS[j] + reactionDH[j] - loguRT * sumNuij[j]);
	}

	// Conventional reactions
	for(j = 1;j <= NR;j++)
		k1[j] = exp(k01[j] + beta1[j] * logT + E1[j] * uT);

	// Landau-Teller reactions
	for(k=1; k<=numLandauTeller;k++)
	{
		j=iLandauTeller[k];
		k1[j] *= exp(LandauTeller[k][1]/pow(T,1./3.) + LandauTeller[k][2]/pow(T,2./3.));
	}

	// Janev-Langer reactions
	for(k=1; k<=numJanevLanger;k++)
	{
		j=iJanevLanger[k];
		double coeffJanevLanger = JanevLanger[k][1]+logT*(JanevLanger[k][2]+logT*(JanevLanger[k][3]+logT*(JanevLanger[k][4]+logT*(JanevLanger[k][5]+logT*(JanevLanger[k][6]+logT*(JanevLanger[k][7]+logT*(JanevLanger[k][8]+logT*JanevLanger[k][9])))))));
		k1[j] = exp( k01[j] + beta1[j]*logT +E1[j]*uT  + coeffJanevLanger );
	}

	// Power-Series reactions
	for(k=1; k<=numPowerSeries;k++)
	{
		j=iPowerSeries[k];
		k1[j]  = exp(k01[j] + beta1[j]*logT);
		k1[j] *= exp( PowerSeries[k][1]/T +	PowerSeries[k][2]/T/T + PowerSeries[k][3]/T/T/T + PowerSeries[k][4]/T/T/T/T);
	}

	// Chebishev-Polynomials reactions
	for(k=1; k<=numChebishev;k++)
	{
		j=iChebishev[k];
		k1[j]  = ChebishevPolynomials[k].GiveMeKappa(T, P_Pa);
	}

	// Logarithmic-Pressure reactions
	for(k=1; k<=numLogarithmicPressure;k++)
	{
		j=iLogarithmicPressure[k];
		k1[j]  = LogarithmicPressure[k].GiveMeKappa(T, P_Pa);
	}

	// Landau-Teller reactions
	for(k=1; k<=numCollisionEfficiency;k++)
	{
		j=iCollisionEfficiency[k];
		double gamma = Min(1., exp(k01[j] + beta1[j] * logT + E1[j] * uT));
		k1[j] = gamma*CollisionEfficiency[k]*sqrt(T);
	}

	// Fall-Off Reactions
	double Fcent,e1Troe,e2Troe,e3Troe;	
	for(k = 1;k <= numFallOff;k++)
	{
		j=iFallOff[k];

		k2[j] = exp(k02[j] + beta2[j] * logT + E2[j] * uT);

		switch(jThirdBody[j])
		{
			case 2: // Lindeman
				break;
			case 3: // Troe
				e1Troe = e2Troe = 0.;
				e3Troe = 0.;
				if(bPressure[j] != 0.)
					e1Troe = exp(-T / bPressure[j]);
				if(cPressure[j] != 0.)
					e2Troe = exp(-T / cPressure[j]);
				if(dPressure[j] != 0.)
					e3Troe = exp(-dPressure[j] / T);
				Fcent = (1. - aPressure[j]) * e1Troe + aPressure[j] * e2Troe + e3Troe;
				if(Fcent < 1.e-300)
						logFcent[j] = -300.;
				else
					logFcent[j] = log10(Fcent);
				break;
			case 4:	// SRI
				e1Troe = 1.;
				e2Troe = 0.;
				if(bPressure[j] != 0.)
					e1Troe = exp(-bPressure[j] / T);
				if(cPressure[j] != 0.)
					e2Troe = exp(-T / cPressure[j]);
				logFcent[j] = (aPressure[j] * e1Troe + e2Troe);
				break;
			default: // Lindeman
				break;
		}
	}

	// CAB Reactions
	for(k = 1;k <= numCABR;k++)
	{
		j=iCABR[k];

		k2[j] = exp(k02[j] + beta2[j] * logT + E2[j] * uT);

		switch(jThirdBody[j])
		{
			case 5: // Lindeman
				break;
			case 6: // Troe
				e1Troe = e2Troe = 0.;
				e3Troe = 0.;
				if(bPressure[j] != 0.)
					e1Troe = exp(-T / bPressure[j]);
				if(cPressure[j] != 0.)
					e2Troe = exp(-T / cPressure[j]);
				if(dPressure[j] != 0.)
					e3Troe = exp(-dPressure[j] / T);
				Fcent = (1. - aPressure[j]) * e1Troe + aPressure[j] * e2Troe + e3Troe;
				if(Fcent < 1.e-300)
						logFcent[j] = -300.;
				else
					logFcent[j] = log10(Fcent);
				break;
			case 7:
				e1Troe = 1.;
				e2Troe = 0.;
				if(bPressure[j] != 0.)
					e1Troe = exp(-bPressure[j] / T);
				if(cPressure[j] != 0.)
					e2Troe = exp(-T / cPressure[j]);
				logFcent[j] = aPressure[j] * e1Troe + e2Troe;
				break;
			default: // Lindeman
				break;
		}
	}

	// TAR reactions
	for(k=1; k<=numTAR;k++)
	{
		j=iTAR[k];
		k1[j]  = exp( TAR_series[k][1] + TAR_series[k][2]*TAR_omegaC + TAR_series[k][3]*TAR_omegaC*TAR_omegaC + 
			          TAR_series[k][4]*log(T)-1./T*(TAR_series[k][5]-TAR_series[k][6]*TAR_omegaC) );

//		cout << "OmegaC:     " << TAR_omegaC << endl;
//		cout << "Reazione:   " << j << endl;
//		cout << "Temp.:      " << T << endl;
//		for (int hh=1;hh<=6;hh++)
//			cout << "Parameters: " << TAR_series[k][hh] << endl;
//		cout << "Kappa:      " << k1[j] << endl;
//		getchar();
	}
}

void OpenSMOKE_Kinetics::UpdateOmegaC(const double omegaC)
{
	TAR_omegaC = omegaC;
}

//----------------------------------------------------------------------------------//
//																					//
//						      computeDHandDSreaction								//
//																					//
//----------------------------------------------------------------------------------//	

void OpenSMOKE_Kinetics::ComputeDHandDSreaction(BzzVector &componentDH, BzzVector &componentDS,BzzVector &reactionDH,BzzVector &reactionDS)
{
	int i, k;

	int *jD1 = jDir1.GetHandle();
	int *jD2 = jDir2.GetHandle();
	int *jD3 = jDir3.GetHandle();
	int *jD4 = jDir4.GetHandle();
	int *jD5 = jDir5.GetHandle();
	double *vD5 = valDir5.GetHandle();

	int *jIT1 = jInvTot1.GetHandle();
	int *jIT2 = jInvTot2.GetHandle();
	int *jIT3 = jInvTot3.GetHandle();
	int *jIT4 = jInvTot4.GetHandle();
	int *jIT5 = jInvTot5.GetHandle();
	double *vIT5 = valInvTot5.GetHandle();
	
	reactionDH = 0.;
	reactionDS = 0.;
			
	for(i = 1;i <= NC;i++)
	{
		for(k = 1;k <= numDir1[i];k++)
		{
			reactionDH[*jD1] -= componentDH[i];
			reactionDS[*jD1] -= componentDS[i];
			jD1++;
		}
		for(k = 1;k <= numDir2[i];k++)
		{
			reactionDH[*jD2] -= (componentDH[i] + componentDH[i]);
			reactionDS[*jD2] -= (componentDS[i] + componentDS[i]);
			jD2++;
		}
		for(k = 1;k <= numDir3[i];k++)
		{
			reactionDH[*jD3] -= (componentDH[i] + componentDH[i] + componentDH[i]);
			reactionDS[*jD3] -= (componentDS[i] + componentDS[i] + componentDS[i]);
			jD3++;
		}
		for(k = 1;k <= numDir4[i];k++)
		{
			reactionDH[*jD4] -= (0.5 * componentDH[i]);
			reactionDS[*jD4] -= (0.5 * componentDS[i]);
			jD4++;
		}
		for(k = 1;k <= numDir5[i];k++)
		{
			reactionDH[*jD5] -= (*vD5 * componentDH[i]);
			reactionDS[*jD5] -= (*vD5 * componentDS[i]);
			jD5++;
			vD5++;
		}

		for(k = 1;k <= numInvTot1[i];k++)
		{
			reactionDH[*jIT1] += componentDH[i];
			reactionDS[*jIT1] += componentDS[i];
			jIT1++;
		}
		for(k = 1;k <= numInvTot2[i];k++)
		{
			reactionDH[*jIT2] += (componentDH[i] + componentDH[i]);
			reactionDS[*jIT2] += (componentDS[i] + componentDS[i]);
			jIT2++;
		}
		for(k = 1;k <= numInvTot3[i];k++)
		{
			reactionDH[*jIT3] += (componentDH[i] + componentDH[i] + componentDH[i]);
			reactionDS[*jIT3] += (componentDS[i] + componentDS[i] + componentDS[i]);
			jIT3++;
		}
		for(k = 1;k <= numInvTot4[i];k++)
		{
			reactionDH[*jIT4] += (0.5 * componentDH[i]);
			reactionDS[*jIT4] += (0.5 * componentDS[i]);
			jIT4++;
		}
		for(k = 1;k <= numInvTot5[i];k++)
		{
			reactionDH[*jIT5] += (*vIT5 * componentDH[i]);
			reactionDS[*jIT5] += (*vIT5 * componentDS[i]);
			jIT5++;
			vIT5++;
		}
	}
}



void OpenSMOKE_Kinetics::GetDerivativesC(double T, double cTot, BzzMatrix *dRC, BzzVector &cRes, BzzVector &R)
{
	ChangeDimensions(NC,NC,dRC);
	BzzVector wR(NC);

	int *jD1 = jDir1.GetHandle();
	int *jD2 = jDir2.GetHandle();
	int *jD3 = jDir3.GetHandle();
	int *jD4 = jDir4.GetHandle();
	int *jD5 = jDir5.GetHandle();
	double *vD5 = valDir5.GetHandle();

	int *jIE1 = jInvEq1.GetHandle();
	int *jIE2 = jInvEq2.GetHandle();
	int *jIE3 = jInvEq3.GetHandle();
	int *jIE4 = jInvEq4.GetHandle();
	int *jIE5 = jInvEq5.GetHandle();
	double *vIE5 = valInvEq5.GetHandle();


	double dc,udc;
	int j,k,kd;

	mc = cRes;
	mR = R;

	mr = reactionRates->r;
	mrDirC = reactionRates->rDirC;
	mrInvC = reactionRates->rInvC;

	// derivata rispetto a ckd
	for(kd = 1;kd <= NC;kd++)
	{
		// if(kd==1) cout << "Jacobian evaluation" << endl;
		// Solo se la concentrazione e' praticamente nulla
		if(mc[kd] <= 1.e-100)
		{
			cRes = mc;
			dc = 1.e-10 + 1.e-12 * cRes[kd];
			udc = 3.e-8 * Max(cRes[kd],1./dc);
			udc = Max(udc,1./dc);
			udc = Max(udc,1.e-19);
			dc = Min(udc,.001 + .001 * cRes[kd]);
			cRes[kd] += dc;
			udc = 1. / dc;

			reactionRates->ComputeFromConcentrations(T, cRes, cTot, &wR);
		
			for(j = 1;j <= mc.Size();j++)
				(*dRC)[j][kd] = (wR[j] - mR[j]) * udc;
			for(k = 1;k <= numDir1[kd];k++)
				jD1++;
			for(k = 1;k <= numDir2[kd];k++)
				jD2++;
			for(k = 1;k <= numDir3[kd];k++)
				jD3++;
			for(k = 1;k <= numDir4[kd];k++)
				jD4++;
			for(k = 1;k <= numDir5[kd];k++)
				{jD5++;vD5++;}

			for(k = 1;k <= numInvEq1[kd];k++)
				jIE1++;
			for(k = 1;k <= numInvEq2[kd];k++)
				jIE2++;
			for(k = 1;k <= numInvEq3[kd];k++)
				jIE3++;
			for(k = 1;k <= numInvEq4[kd];k++)
				jIE4++;
			for(k = 1;k <= numInvEq5[kd];k++)
				{jIE5++;vIE5++;}
			continue;
		}

		cRes = mc;
		dc = (S1 - 1.) * mc[kd];
		if(dc == 0.) BzzError("TODO dc = 0. R derivatives");
		cRes[kd] += dc;
		udc = 1. / dc;

		// 1. CALCOLO DELLE REAZIONI DIRETTE E INVERSE
		reactionRates->rDirC = mrDirC;
		reactionRates->rInvC = mrInvC;		
		for(k = 1;k <= numDir1[kd];k++)
			reactionRates->rDirC[*jD1++] *= S1;
		for(k = 1;k <= numDir2[kd];k++)
			reactionRates->rDirC[*jD2++] *= S2;
		for(k = 1;k <= numDir3[kd];k++)
			reactionRates->rDirC[*jD3++] *= S3;
		for(k = 1;k <= numDir4[kd];k++)
			reactionRates->rDirC[*jD4++] *= SSQ;
		for(k = 1;k <= numDir5[kd];k++)
			reactionRates->rDirC[*jD5++] *= pow(S1,*vD5++);

		for(k = 1;k <= numInvEq1[kd];k++)
			reactionRates->rInvC[*jIE1++] *= S1;
		for(k = 1;k <= numInvEq2[kd];k++)
			reactionRates->rInvC[*jIE2++] *= S2;
		for(k = 1;k <= numInvEq3[kd];k++)
			reactionRates->rInvC[*jIE3++] *= S3;
		for(k = 1;k <= numInvEq4[kd];k++)
			reactionRates->rInvC[*jIE4++] *= SSQ;
		for(k = 1;k <= numInvEq5[kd];k++)
			reactionRates->rInvC[*jIE5++] *= pow(S1,*vIE5++);

		// 2. CONTRIBUTI DALLE REAZIONI CON TERZO CORPO
		ThirdBody(cRes, cTot, reactionRates->coeffM);

		// 3. CONTRIBUTI DERIVANTI DALLE REAZIONI DI FALLOFF
		FallOff(T, reactionRates->coeffM, reactionRates->k1, reactionRates->k2,
		       reactionRates->logFcent, reactionRates->coeffFallOff);

		// 3a. CONTRIBUTI DERIVANTI DALLE REAZIONI DI FALLOFF
		ChemicallyActivatedBimolecularReactions(T, reactionRates->coeffM, reactionRates->k1, reactionRates->k2,
		       reactionRates->logFcent, reactionRates->coeffFallOff);
	
		// 4. ASEMBLAGGIO DEI DIVERSI CONTRIBUTI
		reactionsWithEquilibrium(reactionRates->rDirC, reactionRates->rInvC, reactionRates->uKeq, reactionRates->r);
		reactionsWithThirdBody(reactionRates->coeffM, reactionRates->r);
		
		for(j = 1;j <= NR;j++)
			reactionRates->r[j] *= reactionRates->coeffFallOff[j];

		// 5. COSTRUZIONE DEI CONTRIBUTI PER LE SINGOLE SPECIE
		compositionReactionRates(reactionRates->r, &wR);


		//---------------------------------------------------------------------------
		// Calcolo dello Jacobiano
		//---------------------------------------------------------------------------
		Difference(&wR,mR);
		Product(udc,&wR);
		dRC->SetColumn(kd,wR);
	}

}

void OpenSMOKE_Kinetics::readFromFileBinary(const string fileKin)
{
	int j;

	// Recupero del numero di reazioni e del numero di specie
	// -------------------------------------------------------------------------------
	NR = reactionRates->NumberOfReactions();
	NC = reactionRates->NumberOfSpecies();


	// Inizializzazione dei contatori per l'analisi dei tipi di reazioni che
	// compaiono nello schema cinetico
	// -------------------------------------------------------------------------------
	numThirdBodyOnly		= 0;
	numFallOff				= 0;
	numThirdBody			= 0;
	numConventional			= 0;
	numCABR					= 0;
	numLandauTeller			= 0;
	numJanevLanger			= 0;
	numPowerSeries			= 0;
	numChebishev			= 0;
	numLogarithmicPressure	= 0;
	numIrreversible			= 0;
	numGlobal				= 0;
	numCollisionEfficiency	= 0;
	numTAR					= 0;

	// Dimensionamento delle matrici e dei vettori
	// -------------------------------------------------------------------------------
	iEfficiency	 =	 new BzzVectorInt [NR + 1];
	efficiency	 =	 new BzzVector [NR + 1];

	ChangeDimensions(NR,&jEquil);
	ChangeDimensions(NR,&jThirdBody);
	ChangeDimensions(NR,&jNumEfficiency);
	ChangeDimensions(NR,&beta1);
	ChangeDimensions(NR,&beta2);
	ChangeDimensions(NR,&E1);
	ChangeDimensions(NR,&E2);
	ChangeDimensions(NR,&k01);
	ChangeDimensions(NR,&exp_k01);
	ChangeDimensions(NR,&k02);
	ChangeDimensions(NR,&exp_k02);
	ChangeDimensions(NR,&aPressure);
	ChangeDimensions(NR,&bPressure);
	ChangeDimensions(NR,&cPressure);
	ChangeDimensions(NR,&dPressure);
	ChangeDimensions(NR,&ePressure);
	ChangeDimensions(NUM_MAX,&iFallOff);
	ChangeDimensions(NUM_MAX,&iCABR);
	ChangeDimensions(NUM_MAX,&iThirdBody);
	ChangeDimensions(NUM_MAX,&iThirdBodyOnly);
	ChangeDimensions(NUM_MAX,&iLandauTeller);
	ChangeDimensions(NUM_MAX,&iPowerSeries);
	ChangeDimensions(NUM_MAX,&iJanevLanger);
	ChangeDimensions(NUM_MAX,&iChebishev);
	ChangeDimensions(NUM_MAX,&iLogarithmicPressure);
	ChangeDimensions(NUM_MAX,&iCollisionEfficiency);
	ChangeDimensions(NUM_MAX,&iTAR);

	OpenSMOKE_ChebishevPolynomialsReaction* ChebishevPolynomialsProvisional;
	ChebishevPolynomialsProvisional = new OpenSMOKE_ChebishevPolynomialsReaction[NUM_MAX_CHEBISHEV+1];
	OpenSMOKE_LogarithmicPressureReaction* LogarithmicPressureProvisional;
	LogarithmicPressureProvisional = new OpenSMOKE_LogarithmicPressureReaction[NUM_MAX_LOGPRESSURE+1];


	// 1. Reading reaction data
	// -------------------------------------------------------------------------------
	BzzLoad inputFile;
	if (reactionRates->binary_version_ == true)	inputFile('*', fileKin);
	else										inputFile(fileKin);

	inputFile >> NC;
	inputFile >> NR;

	for(j = 1;j <= NR;j++)
	{
		inputFile  >> jEquil[j] >> jThirdBody[j] >> jNumEfficiency[j];
		
		if(jNumEfficiency[j] == -1)
		{
			ChangeDimensions(1,&iEfficiency[j]);
			ChangeDimensions(1,&efficiency[j]);
			inputFile >> iEfficiency[j][1] >> efficiency[j][1];
		}
		
		// Reazioni con il terzo corpo: lettura delle efficienze
		else if(jNumEfficiency[j] > 0)
		{
			ChangeDimensions(jNumEfficiency[j],&iEfficiency[j]);
			ChangeDimensions(jNumEfficiency[j],&efficiency[j]);
			for(int ke = 1;ke <= jNumEfficiency[j];ke++)
			{
				inputFile  >> iEfficiency[j][ke] >> efficiency[j][ke];
				efficiency[j][ke]-=1.;
			}
		}

		// Aggiornamento delle variabili per il conteggio e l'analisi dei tipi di
		// reazione coinvolte nello schema cinetico
		if(jThirdBody[j]==0 && jNumEfficiency[j]==0) 
			numConventional++;
		else
		{
			// - tutte le reazioni di FallOff-CABR sono automaticamente anche di terzo corpo
			//   e dunque si rende comunque necessario il calcolo del coefficiente correttivo coeffM

			// la reazione e' (solo) di terzo corpo (ma la dip. dalla temperatura puo' essere LT, JAN, ecc.)
			if( jThirdBody[j]==1  || jThirdBody[j]==10 || 
				jThirdBody[j]==11 || jThirdBody[j]==12 ||
				jThirdBody[j]==13 || jThirdBody[j]==14  )
			{
				numThirdBodyOnly++;
				iThirdBodyOnly[numThirdBodyOnly]=j;
			}

			// reazioni con terzo corpo (sia quelle pure che di falloff)
			if(jThirdBody[j]<100) 
			{  
				numThirdBody++;
				iThirdBody[numThirdBody]=j;
			}

			// in questo caso la reazione e' sicuramente di FallOff e quindi in realta anche di terzo corpo
			if(jThirdBody[j]>=2 && jThirdBody[j]<=4) 
			{
				numFallOff++;
				iFallOff[numFallOff]=j;
			}

			if(jThirdBody[j]>=5 && jThirdBody[j]<=7) 
			{
				numCABR++;
				iCABR[numCABR]=j;
			}

			if(jThirdBody[j]==10 || jThirdBody[j]==100) 
			{
				numLandauTeller++;
				iLandauTeller[numLandauTeller]=j;
			}

			if(jThirdBody[j]==11 || jThirdBody[j]==110) 
			{
				numJanevLanger++;
				iJanevLanger[numJanevLanger]=j;
			}

			if(jThirdBody[j]==12 || jThirdBody[j]==120) 
			{
				numPowerSeries++;
				iPowerSeries[numPowerSeries]=j;
			}

			if(jThirdBody[j]==160) 
			{
				numTAR++;
				iTAR[numTAR]=j;
			}

			if(jThirdBody[j]==13 || jThirdBody[j]==130) 
			{
				numChebishev++;
				iChebishev[numChebishev]=j;
			}

			if(jThirdBody[j]==14 || jThirdBody[j]==140) 
			{
				numLogarithmicPressure++;
				iLogarithmicPressure[numLogarithmicPressure]=j;
			}

			if(jThirdBody[j]==15 || jThirdBody[j]==150) 
			{
				numCollisionEfficiency++;
				iCollisionEfficiency[numCollisionEfficiency]=j;
			}
		}

		// Reading Kinetic data
		inputFile  >> k01[j] >> beta1[j] >> E1[j];
		
		if (k01[j] < 0.)
		{
			k01[j] = -k01[j];
			negativeSigns.Append(j);
		}

		exp_k01[j] = k01[j];
		k01[j] = log(k01[j]);
		E1[j] = -E1[j] /  Constants::R_cal_mol;

		// In case of Fall-Off reactions
		if(jThirdBody[j]>=2 && jThirdBody[j]<=7)
		{
			inputFile  >> k02[j] >> beta2[j] >> E2[j];
			exp_k02[j] = k02[j];
			k02[j] = log(k02[j]);
			E2[j] = -E2[j] /  Constants::R_cal_mol;
					
			inputFile  >> aPressure[j] >> bPressure[j] >> cPressure[j] >> dPressure[j] >> ePressure[j];
		}

		// Landau-Teller
		if(jThirdBody[j]==10 || jThirdBody[j]==100)
		{
			BzzVector aux(2);
			inputFile  >> aux[1] >> aux[2];
			LandauTeller.AppendRow(aux);
		}

		// Janev-Langer
		if(jThirdBody[j]==11 || jThirdBody[j]==110)
		{
			BzzVector aux(9);
			inputFile  >> aux[1] >> aux[2] >> aux[3] >> aux[4] >> aux[5] >> aux[6] >> aux[7] >> aux[8] >> aux[9];
			JanevLanger.AppendRow(aux);
		}

		// Power-Series
		if(jThirdBody[j]==12 || jThirdBody[j]==120)
		{
			BzzVector aux(4);
			inputFile  >> aux[1] >> aux[2] >> aux[3] >> aux[4];
			PowerSeries.AppendRow(aux);
		}

		// Power-Series
		if(jThirdBody[j]==160)
		{
			BzzVector aux(6);
			inputFile  >> aux[1] >> aux[2] >> aux[3] >> aux[4] >> aux[5] >> aux[6];
			TAR_series.AppendRow(aux);
		}

		// Chebishev
		if(jThirdBody[j]==13 || jThirdBody[j]==130)
			ChebishevPolynomialsProvisional[numChebishev].ReadFromFile(inputFile);

		// Chebishev
		if(jThirdBody[j]==14 || jThirdBody[j]==140)
			LogarithmicPressureProvisional[numLogarithmicPressure].ReadFromFile(inputFile);

		// Landau-Teller
		if(jThirdBody[j]==15 || jThirdBody[j]==150)
		{
			double aux;
			inputFile  >> aux;
			CollisionEfficiency.Append(aux);
		}
	}

	// 2. Cleaning matrices
	// -------------------------------------------------------------------------------
	iFallOff.DeleteLastNElements(NUM_MAX-numFallOff);
	iCABR.DeleteLastNElements(NUM_MAX-numCABR);
	iThirdBody.DeleteLastNElements(NUM_MAX-numThirdBody);
	iThirdBodyOnly.DeleteLastNElements(NUM_MAX-numThirdBodyOnly);
	iLandauTeller.DeleteLastNElements(NUM_MAX-numLandauTeller);
	iJanevLanger.DeleteLastNElements(NUM_MAX-numJanevLanger);
	iPowerSeries.DeleteLastNElements(NUM_MAX-numPowerSeries);
	iTAR.DeleteLastNElements(NUM_MAX-numTAR);
	iChebishev.DeleteLastNElements(NUM_MAX-numChebishev);
	iCollisionEfficiency.DeleteLastNElements(NUM_MAX-numCollisionEfficiency);

	ChebishevPolynomials = new OpenSMOKE_ChebishevPolynomialsReaction[numChebishev+1];
	for(j=1;j<=numChebishev;j++)
		ChebishevPolynomials[j] = ChebishevPolynomialsProvisional[j];
	LogarithmicPressure = new OpenSMOKE_LogarithmicPressureReaction[numLogarithmicPressure+1];
	for(j=1;j<=numLogarithmicPressure;j++)
		LogarithmicPressure[j] = LogarithmicPressureProvisional[j];

	// 3. Equilibrium reactions
	// -------------------------------------------------------------------------------
	numEquilibrium = 0;
	for(j = 1;j <= NR;j++)
	{
		if(jEquil[j] == 1)	numEquilibrium++;
		if(jEquil[j] == 0)	numIrreversible++;
	}

	ChangeDimensions(numEquilibrium,&reactionWithEquil);	
	int neq = 0;
	for(j = 1;j <= NR;j++)
	{
		if(jEquil[j] == 1)
			reactionWithEquil[++neq] = j;
	}


	// 4. Importing stoichiometry data
	// -------------------------------------------------------------------------------
	readStoichiometricFileBinary(inputFile);

	// 5. Reaction names
	// -------------------------------------------------------------------------------
	// Allocating memory
	reactionRates->strReaction = new char* [NR + 1];
	for(j =1;j<=NR;j++)
		reactionRates->strReaction[j] = new char [Constants::REACTION_NAME_SIZE];

	// Lettura dei nomi delle specie e del numero di componenti
	for(j =1;j<=NR;j++)
	{
		if (reactionRates->binary_version_ == true)	
		{
			char name[Constants::REACTION_NAME_SIZE];
			inputFile.fileLoad.read((char*) name, sizeof(name));
			strcpy(reactionRates->strReaction[j], name);
		}
		else
		{
			inputFile.fileLoad >> reactionRates->strReaction[j];
		}
	}

	// Writing additional info on stoichiometry
	if (reactionRates->iVersion == V090905 || reactionRates->iVersion == V101116)
	{
		inputFile >> forwardOrders;
		inputFile >> backwardOrders;
		inputFile >> numGlobal;
		inputFile >> reactionRates->indexSoot;

	//	if (numGlobal == 0)
	//	{	
			lambda_jDir1	 = jDir1;
			lambda_jDir2	 = jDir2;
			lambda_jDir3	 = jDir3;
			lambda_jDir4	 = jDir4;
			lambda_jDir5	 = jDir5;
			lambda_valDir5	 = valDir5;

			lambda_jInvEq1   = jInvEq1;
			lambda_jInvEq2   = jInvEq2;
			lambda_jInvEq3   = jInvEq3;
			lambda_jInvEq4   = jInvEq4;
			lambda_jInvEq5   = jInvEq5;
			lambda_valInvEq5 = valInvEq5;

			lambda_numDir1	 = numDir1;
			lambda_numDir2	 = numDir2;
			lambda_numDir3	 = numDir3;
			lambda_numDir4	 = numDir4;
			lambda_numDir5	 = numDir5;

			lambda_numInvEq1 = numInvEq1;
			lambda_numInvEq2 = numInvEq2;
			lambda_numInvEq3 = numInvEq3;
			lambda_numInvEq4 = numInvEq4;
			lambda_numInvEq5 = numInvEq5;
	//	}
		//else
		if (numGlobal > 0)
		{
			iGlobalDirect       = new BzzVectorInt[numGlobal+1];
			iGlobalInverse      = new BzzVectorInt[numGlobal+1];
			lambdaGlobalDirect  = new BzzVector[numGlobal+1];
			lambdaGlobalInverse = new BzzVector[numGlobal+1];

		/*	inputFile >> lambda_numDir1;
			inputFile >> lambda_numDir2;
			inputFile >> lambda_numDir3;
			inputFile >> lambda_numDir4;
			inputFile >> lambda_numDir5;

			inputFile >> lambda_numInvEq1;
			inputFile >> lambda_numInvEq2;
			inputFile >> lambda_numInvEq3;
			inputFile >> lambda_numInvEq4;
			inputFile >> lambda_numInvEq5;

			inputFile >> lambda_jDir1;
			inputFile >> lambda_jDir2;
			inputFile >> lambda_jDir3;
			inputFile >> lambda_jDir4;
			inputFile >> lambda_jDir5;
			inputFile >> lambda_valDir5;

			inputFile >> lambda_jInvEq1;
			inputFile >> lambda_jInvEq2;
			inputFile >> lambda_jInvEq3;
			inputFile >> lambda_jInvEq4;
			inputFile >> lambda_jInvEq5;
			inputFile >> lambda_valInvEq5;
		*/
			inputFile >> iGlobal;

			for (int n=1;n<=numGlobal;n++)	
			{	
				inputFile >> iGlobalDirect[n];
				inputFile >> lambdaGlobalDirect[n];
				
				int tag;
				inputFile >> tag;
				if (tag == 1)
				{
					inputFile >> iGlobalInverse[n];
					inputFile >> lambdaGlobalInverse[n];
				}
			}
		}

		if (reactionRates->indexSoot != 0)	reactionRates->iSootMode = true;
	}
	
	// 6. End
	// -------------------------------------------------------------------------------
	inputFile.End();

	
	// 5. Riepilogo delle informazioni a video
	// -------------------------------------------------------------------------------
	cout << "--------------------------------------------------------------------------" << endl;
	cout << "                       KINETIC SCHEME SUMMARY                             " << endl;
	cout << "--------------------------------------------------------------------------" << endl;
	cout << "  Total number of species                   = " << NC						 << endl;
	cout << "  Total number of reactions                 = " << NR						 << endl;
	cout << "    Number of irreversible reactions        = " << numIrreversible			 << endl;
	cout << "    Number of equilibrium reactions         = " << numEquilibrium			 << endl;
	cout << "    Number of conventional reactions        = " << numConventional			 << endl;
	cout << "    Number of third body reactions (total)  = " << numThirdBody		 	 << endl;
	cout << "    Number of third body reactions          = " << numThirdBodyOnly		 << endl;
	cout << "    Number of Fall-Off reactions            = " << numFallOff				 << endl;
	cout << "    Number of C.A.B. reactions              = " << numCABR					 << endl;
	cout << "    Number of Landau-Teller reactions       = " << numLandauTeller			 << endl;
	cout << "    Number of Janev-Langer reactions        = " << numJanevLanger			 << endl;
	cout << "    Number of Power-Series reactions        = " << numPowerSeries			 << endl;
	cout << "    Number of Chebishev-Pol. reactions      = " << numChebishev			 << endl;
	cout << "    Number of Log-Pressure reactions        = " << numLogarithmicPressure	 << endl;
	cout << "    Number of Global reactions              = " << numGlobal				 << endl;
	cout << "    Number of Bimolecular Coll. reactions   = " << numCollisionEfficiency	 << endl;
	cout << "    Number of TAR reactions                 = " << numTAR	                 << endl;
	cout << "--------------------------------------------------------------------------" << endl;
	cout << endl;

	// TODO
	TAR_omegaC = 0.;
};


void OpenSMOKE_Kinetics::readStoichiometricFileBinary(BzzLoad &inputFile)
{
	inputFile >> numDir1;
	inputFile >> numDir2;
	inputFile >> numDir3;
	inputFile >> numDir4;
	inputFile >> numDir5;

	inputFile >> numInvTot1;
	inputFile >> numInvTot2;
	inputFile >> numInvTot3;
	inputFile >> numInvTot4;
	inputFile >> numInvTot5;

	inputFile >> numInvEq1;
	inputFile >> numInvEq2;
	inputFile >> numInvEq3;
	inputFile >> numInvEq4;
	inputFile >> numInvEq5;

	inputFile >> jDir1;
	inputFile >> jDir2;
	inputFile >> jDir3;
	inputFile >> jDir4;
	inputFile >> jDir5;
	inputFile >> valDir5;

	inputFile >> jInvTot1;
	inputFile >> jInvTot2;
	inputFile >> jInvTot3;
	inputFile >> jInvTot4;
	inputFile >> jInvTot5;
	inputFile >> valInvTot5;

	inputFile >> jInvEq1;
	inputFile >> jInvEq2;
	inputFile >> jInvEq3;
	inputFile >> jInvEq4;
	inputFile >> jInvEq5;
	inputFile >> valInvEq5;

	inputFile >> sumNuij;
}

void OpenSMOKE_Kinetics::GiveMeIndexOfSpeciesInEachReaction(const string fileName, BzzVectorIntArray &indices)
{
	int i,j,k;

	BzzMatrix gamma;

	gamma = constructGamma(fileName);

	// Forward reactions
	// -----------------------------------------
	for(i=1;i<=NR;i++)
	{
		BzzVectorInt index;
		for(j=1;j<=NC;j++)
			if (gamma[j][i] < 0) 
				index.Append(j);

		for(k = 1;k <= numEquilibrium;k++)
		{	
			if (reactionWithEquil[k]==i)
			{
				for(j=1;j<=NC;j++)
					if (gamma[j][i] > 0) 
						index.Append(j);
				break;
			}
		}

		indices.SetVector(i,index);
	}
}

int OpenSMOKE_Kinetics::IsAReversibleReaction(const int index)
{
	for (int j=1;j<=numEquilibrium;j++)
		if (reactionWithEquil[j] == index)
			return j;
	return 0;
}
int OpenSMOKE_Kinetics::IsAFallOffReaction(const int index)
{
	for (int j=1;j<=numFallOff;j++)
		if (iFallOff[j] == index)
			return j;
	return 0;
}

int OpenSMOKE_Kinetics::IsACABReaction(const int index)
{
	for (int j=1;j<=numCABR;j++)
		if (iCABR[j] == index)
			return j;
	return 0;
}

int OpenSMOKE_Kinetics::IsALandauTellerReaction(const int index)
{
	for (int j=1;j<=numLandauTeller;j++)
		if (iLandauTeller[j] == index)
			return j;
	return 0;
}

int OpenSMOKE_Kinetics::IsAJanevLangerReaction(const int index)
{
	for (int j=1;j<=numJanevLanger;j++)
		if (iJanevLanger[j] == index)
			return j;
	return 0;
}

int OpenSMOKE_Kinetics::IsAPowerSeriesReaction(const int index)
{
	for (int j=1;j<=numPowerSeries;j++)
		if (iPowerSeries[j] == index)
			return j;
	return 0;
}

int OpenSMOKE_Kinetics::IsATARReaction(const int index)
{
	for (int j=1;j<=numTAR;j++)
		if (iTAR[j] == index)
			return j;
	return 0;
}

int OpenSMOKE_Kinetics::IsAChebishevPolynomialsReaction(const int index)
{
	for (int j=1;j<=numChebishev;j++)
		if (iChebishev[j] == index)
			return j;
	return 0;
}

int OpenSMOKE_Kinetics::IsALogarithmicPressureReaction(const int index)
{
	for (int j=1;j<=numLogarithmicPressure;j++)
		if (iLogarithmicPressure[j] == index)
			return j;
	return 0;
}

int OpenSMOKE_Kinetics::IsACollisionEfficiencyReaction(const int index)
{
	for (int j=1;j<=numCollisionEfficiency;j++)
		if (iCollisionEfficiency[j] == index)
			return j;
	return 0;
}

void OpenSMOKE_Kinetics::KineticExpressionString(ofstream &fOutput, const int j)
{
	int jReversibleReaction		= IsAReversibleReaction(j);
	int jFallOffReaction		= IsAFallOffReaction(j);
	int jCABReaction			= IsACABReaction(j);
	int jLandauTeller			= IsALandauTellerReaction(j);
	int jCollisionEfficiency	= IsACollisionEfficiencyReaction(j);
	int jPowerSeries			= IsAPowerSeriesReaction(j);
	int jChebishev				= IsAChebishevPolynomialsReaction(j);
	int jLogarithmicPressure	= IsALogarithmicPressureReaction(j);
	int jJanevLanger			= IsAJanevLangerReaction(j);
	int jTAR					= IsATARReaction(j);

	fOutput << "   Change in moles in the reaction = " << sumNuij[j] << endl;
	fOutput << "   Reaction order (forward)        = " << forwardOrders[j] << endl;
	if (jReversibleReaction != 0)
		fOutput << "   Reaction order (backward)       = " << backwardOrders[j] << endl;

	// --------------------------------------------------------------------------------------------------------------------------------------------
	// Kinetic expression - High Pressure for Fall-Off
	// --------------------------------------------------------------------------------------------------------------------------------------------
	if (jFallOffReaction != 0)
	{
		fOutput << "   Fall-off Reaction" << endl;
		fOutput << "   kLow = " << setw(6) << setprecision(4) << scientific << exp(k01[j]);
		if (beta1[j] != 0) fOutput << " T^" << setw(4) << setprecision(3) << fixed << beta1[j];
		fOutput << " exp(" << setw(6) << setprecision(2) << fixed << E1[j]*Constants::R_cal_mol << "/RT)";	
		WriteUnits_m_kmol_s(fOutput, forwardOrders[j]);

		fOutput << "   kLow = " << setw(6) << setprecision(4) << scientific << exp(k01[j])*pow(1.e3,forwardOrders[j]-1);
		if (beta1[j] != 0) fOutput << " T^" << setw(4) << setprecision(3) << fixed << beta1[j];
		fOutput << " exp(" << setw(6) << setprecision(2) << fixed << E1[j]*Constants::R_cal_mol << "/RT)";	
		WriteUnits_cm_mol_s(fOutput, forwardOrders[j]);

		fOutput << "   kHigh = " << setw(6) << setprecision(4) << scientific << exp(k02[j])*pow(1.e3,forwardOrders[j]-1.);
		if (beta1[j] != 0) fOutput << " T^" << setw(4) << setprecision(3) << fixed << beta2[j];
		fOutput << " exp(" << setw(6) << setprecision(2) << fixed << E2[j]*Constants::R_cal_mol << "/RT)";	
		WriteUnits_cm_mol_s(fOutput, forwardOrders[j]-1.);

		fOutput << "   kHigh = " << setw(6) << setprecision(4) << scientific << exp(k02[j])*pow(1.e3,forwardOrders[j]-2.);
		if (beta2[j] != 0) fOutput << " T^" << setw(4) << setprecision(3) << fixed << beta2[j];
		fOutput << " exp(" << setw(6) << setprecision(2) << fixed << E2[j]*Constants::R_cal_mol << "/RT)";	
		WriteUnits_cm_mol_s(fOutput, forwardOrders[j]-1.);
	}

	// --------------------------------------------------------------------------------------------------------------------------------------------
	// Kinetic expression - Low Pressure for Chemically Activated Based Reactions
	// --------------------------------------------------------------------------------------------------------------------------------------------
	else if (jCABReaction != 0)
	{
		fOutput << "   Chemically Activated Bimolecular Reaction" << endl;
		fOutput << "   kHigh = " << setw(6) << setprecision(4) << scientific << exp(k01[j])*pow(1.e3,forwardOrders[j]-2.);
		if (beta1[j] != 0) fOutput << " T^" << setw(4) << setprecision(3) << fixed << beta1[j];
		fOutput << " exp(" << setw(6) << setprecision(2) << fixed << E1[j]*Constants::R_cal_mol << "/RT)";	
		WriteUnits_cm_mol_s(fOutput, forwardOrders[j]-2.);

		fOutput << "   kHigh = " << setw(6) << setprecision(4) << scientific << exp(k01[j])*pow(1.e3,forwardOrders[j]-3.);
		if (beta1[j] != 0) fOutput << " T^" << setw(4) << setprecision(3) << fixed << beta1[j];
		fOutput << " exp(" << setw(6) << setprecision(2) << fixed << E1[j]*Constants::R_cal_mol << "/RT)";	
		WriteUnits_cm_mol_s(fOutput, forwardOrders[j]-2.);

		fOutput << "   kLow  = " << setw(6) << setprecision(4) << scientific << exp(k02[j])*pow(1.e3,forwardOrders[j]-1.);
		if (beta1[j] != 0) fOutput << " T^" << setw(4) << setprecision(3) << fixed << beta2[j];
		fOutput << " exp(" << setw(6) << setprecision(2) << fixed << E2[j]*Constants::R_cal_mol << "/RT)";	
		WriteUnits_m_kmol_s(fOutput, forwardOrders[j]-1.);

		fOutput << "   kLow  = " << setw(6) << setprecision(4) << scientific << exp(k02[j])*pow(1.e3,forwardOrders[j]-2.);
		if (beta2[j] != 0) fOutput << " T^" << setw(4) << setprecision(3) << fixed << beta2[j];
		fOutput << " exp(" << setw(6) << setprecision(2) << fixed << E2[j]*Constants::R_cal_mol << "/RT)";	
		WriteUnits_cm_mol_s(fOutput, forwardOrders[j]-1.);
	}

	else if (jLandauTeller != 0)
	{
		fOutput << "   Landau-Teller reaction" << endl;
		fOutput << "   k = " << setw(6) << setprecision(4) << scientific << exp(k01[j]);
		if (beta1[j] != 0) fOutput << " T^" << setw(4) << setprecision(3) << fixed << beta1[j];
		fOutput << " exp("  << setw(6) << setprecision(2) << fixed << E1[j]*Constants::R_cal_mol << "/RT)";
		fOutput << " exp( " << setprecision(3) << LandauTeller[jLandauTeller][1] << "/T^(1/3) + " 
						    << setprecision(3) << LandauTeller[jLandauTeller][2] << "/T^(2/3) )";
		WriteUnits_m_kmol_s(fOutput, forwardOrders[j]);

		fOutput << "   k = " << setw(6) << setprecision(4) << scientific << exp(k01[j])*pow(1.e3,forwardOrders[j]-1);
		if (beta1[j] != 0) fOutput << " T^" << setw(4) << setprecision(3) << fixed << beta1[j];
		fOutput << " exp("  << setw(6) << setprecision(2) << fixed << E1[j]*Constants::R_cal_mol << "/RT)";	
		fOutput << " exp( " << setprecision(3) << LandauTeller[jLandauTeller][1] << "/T^(1/3) + " 
						    << setprecision(3) << LandauTeller[jLandauTeller][2] << "/T^(2/3) )";
		WriteUnits_cm_mol_s(fOutput, forwardOrders[j]);
	}

	else if (jPowerSeries != 0)
	{
		fOutput << "   Power-series expansion reaction" << endl;
		fOutput << "   k = " << setw(6) << setprecision(4) << scientific << exp(k01[j]);
		if (beta1[j] != 0) fOutput << " T^" << setw(4) << setprecision(3) << fixed << beta1[j];
		fOutput << " exp( " << setprecision(3) << PowerSeries[jPowerSeries][1] << "/T + " 
						    << setprecision(3) << PowerSeries[jPowerSeries][2] << "/T^2 +"
						    << setprecision(3) << PowerSeries[jPowerSeries][3] << "/T^3 +"
						    << setprecision(3) << PowerSeries[jPowerSeries][4] << "/T^4 )";
		WriteUnits_m_kmol_s(fOutput, forwardOrders[j]);

		fOutput << "   k = " << setw(6) << setprecision(4) << scientific << exp(k01[j])*pow(1.e3,forwardOrders[j]-1);
		if (beta1[j] != 0) fOutput << " T^" << setw(4) << setprecision(3) << fixed << beta1[j];
		fOutput << " exp( " << setprecision(3) << PowerSeries[jPowerSeries][1] << "/T + " 
						    << setprecision(3) << PowerSeries[jPowerSeries][2] << "/T^2 +"
						    << setprecision(3) << PowerSeries[jPowerSeries][3] << "/T^3 +"
						    << setprecision(3) << PowerSeries[jPowerSeries][4] << "/T^4 )";
		WriteUnits_cm_mol_s(fOutput, forwardOrders[j]);
	}

	else if (jTAR != 0)
	{
		fOutput << "   TAR reaction" << endl;
		fOutput << "   k = " 
				<< " exp( " << setprecision(3) << TAR_series[jTAR][1] << " + "
							<< setprecision(3) << TAR_series[jTAR][2] << "*omC + " 
							<< setprecision(3) << TAR_series[jTAR][3] << "*omC^2 + " 
							<< setprecision(3) << TAR_series[jTAR][4] << "*ln(T) - 1/T*( "
							<< setprecision(3) << TAR_series[jTAR][5] << " + "
							<< setprecision(3) << TAR_series[jTAR][6] << " *omegaC ) )";

		WriteUnits_m_kmol_s(fOutput, forwardOrders[j]);

		// TODO
	}

	else if (jJanevLanger != 0)
	{
		fOutput << "   Janev-Langer Reaction" << endl;
		fOutput << "   k = " << setw(6) << setprecision(4) << scientific << exp(k01[j]);
		if (beta1[j] != 0) fOutput << " T^" << setw(4) << setprecision(3) << fixed << beta1[j];
		fOutput << " exp(" << setw(6) << setprecision(2) << fixed << E1[j]*Constants::R_cal_mol << "/RT)";	
		fOutput << " exp( " << setprecision(3) << JanevLanger[jJanevLanger][1] << " lnT + " 
						    << setprecision(3) << JanevLanger[jJanevLanger][2] << " (lnT)^2 +"
						    << setprecision(3) << JanevLanger[jJanevLanger][3] << " (lnT)^3 +"
						    << setprecision(3) << JanevLanger[jJanevLanger][4] << " (lnT)^4 )";
		WriteUnits_m_kmol_s(fOutput, forwardOrders[j]);

		fOutput << "   k = " << setw(6) << setprecision(4) << scientific << exp(k01[j])*pow(1.e3,forwardOrders[j]-1);
		if (beta1[j] != 0) fOutput << " T^" << setw(4) << setprecision(3) << fixed << beta1[j];
		fOutput << " exp(" << setw(6) << setprecision(2) << fixed << E1[j]*Constants::R_cal_mol << "/RT)";	
		fOutput << " exp( " << setprecision(3) << JanevLanger[jJanevLanger][1] << " lnT + " 
						    << setprecision(3) << JanevLanger[jJanevLanger][2] << " (lnT)^2 +"
						    << setprecision(3) << JanevLanger[jJanevLanger][3] << " (lnT)^3 +"
						    << setprecision(3) << JanevLanger[jJanevLanger][4] << " (lnT)^4 )";
		WriteUnits_cm_mol_s(fOutput, forwardOrders[j]);
	}

	else if (jChebishev != 0)
	{
		fOutput << "   Chebishev Polynomial Expansion Reaction" << endl;
		fOutput << "   k = Chebishev Polynomial Expansion";
		WriteUnits_m_kmol_s(fOutput, forwardOrders[j]);

		fOutput << "   k = Chebishev Polynomial Expansion";
		WriteUnits_cm_mol_s(fOutput, forwardOrders[j]);

		ChebishevPolynomials[jChebishev].WriteCoefficientMatrix(fOutput);
	}

	else if (jLogarithmicPressure != 0)
	{
		fOutput << "   Logarithmic Pressure Interpolation Reaction" << endl;

		fOutput << "   k = A T^Beta exp(-E/RT)";
		WriteUnits_m_kmol_s(fOutput, forwardOrders[j]);

		fOutput << "   k = A T^Beta exp(-E/RT)";
		WriteUnits_cm_mol_s(fOutput, forwardOrders[j]);

		LogarithmicPressure[jLogarithmicPressure].WriteCoefficientMatrix(fOutput);
	}

	else if (jCollisionEfficiency != 0)
	{
		fOutput << "   Bimolecular Reaction with Collisional Efficiency" << endl;

		fOutput << "   k = " << setw(6) << setprecision(4) << scientific << exp(k01[j]);
		if (beta1[j] != 0) fOutput << " T^" << setw(4) << setprecision(3) << fixed << beta1[j];
		fOutput << " exp(" << setw(6) << setprecision(2) << fixed << E1[j]*Constants::R_cal_mol << "/RT)";	
		
		WriteUnits_m_kmol_s(fOutput, forwardOrders[j]);

		fOutput << "   k = " << setw(6) << setprecision(4) << scientific << exp(k01[j])*pow(1.e3,forwardOrders[j]-1);
		if (beta1[j] != 0) fOutput << " T^" << setw(4) << setprecision(3) << fixed << beta1[j];
		fOutput << " exp(" << setw(6) << setprecision(2) << fixed << E1[j]*Constants::R_cal_mol << "/RT)";	

		WriteUnits_cm_mol_s(fOutput, forwardOrders[j]);
	}

	else
	{
		fOutput << "   k = " << setw(6) << setprecision(4) << scientific << exp(k01[j]);
		if (beta1[j] != 0) fOutput << " T^" << setw(4) << setprecision(3) << fixed << beta1[j];
		fOutput << " exp(" << setw(6) << setprecision(2) << fixed << E1[j]*Constants::R_cal_mol << "/RT)";	
		
		WriteUnits_m_kmol_s(fOutput, forwardOrders[j]);

		fOutput << "   k = " << setw(6) << setprecision(4) << scientific << exp(k01[j])*pow(1.e3,forwardOrders[j]-1);
		if (beta1[j] != 0) fOutput << " T^" << setw(4) << setprecision(3) << fixed << beta1[j];
		fOutput << " exp(" << setw(6) << setprecision(2) << fixed << E1[j]*Constants::R_cal_mol << "/RT)";	

		WriteUnits_cm_mol_s(fOutput, forwardOrders[j]);
	}

	// --------------------------------------------------------------------------------------------------------------------------------------------
	// Reverse reaction
	// --------------------------------------------------------------------------------------------------------------------------------------------
	if (jReversibleReaction != 0)
	{
		fOutput << "   Reverse reaction units:";
		WriteUnits_m_kmol_s(fOutput, backwardOrders[j]);
	}

	if (jFallOffReaction != 0 || jCABReaction !=0)
	{
		if (jThirdBody[j] == 3 || jThirdBody[j] == 6)	
		{
			fOutput << "   Troe Parameters:  " << endl;
			fOutput << "    a    = " << scientific << setprecision(4) << aPressure[j] << endl;
			fOutput << "    T*** = " << scientific << setprecision(4) << bPressure[j] << " K" << endl;
			fOutput << "    T*   = " << scientific << setprecision(4) << cPressure[j] << " K" << endl;
			fOutput << "    T**  = " << scientific << setprecision(4) << dPressure[j] << " K" << endl;
		}
		
		if (jThirdBody[j] == 4 || jThirdBody[j] == 7)	
		{
			fOutput << "   SRI Parameters:  " << endl;
			fOutput << "    a = " << scientific << setprecision(4) << aPressure[j] << endl;
			fOutput << "    b = " << scientific << setprecision(4) << bPressure[j] << " K" << endl;
			fOutput << "    c = " << scientific << setprecision(4) << cPressure[j] << " K" << endl;
			fOutput << "    d = " << scientific << setprecision(4) << dPressure[j] << endl;
			fOutput << "    e = " << scientific << setprecision(4) << ePressure[j] << endl;
		}
	}

	if (jNumEfficiency[j] > 0)
	{
			fOutput << "   Three-body efficiencies:  " << endl;
			for(int k=1;k<=jNumEfficiency[j];k++)
				fOutput << "    " << setw(20) << left << reactionRates->names[iEfficiency[j][k]]
								  << setw(6) << setprecision(3) << left << efficiency[j][k]+1. << endl;
	}

	for(int i=1;i<=iGlobal.Size();i++)
	if (iGlobal[i] == j)	
	{
		int k;

		fOutput << "   Reaction orders (forward)" << endl;
		for(k=1;k<=iGlobalDirect[i].Size();k++)
			fOutput << "     "	<< setw(20) << left << reactionRates->names[iGlobalDirect[i][k]] 
								<< setw(6) << setprecision(3) << left << lambdaGlobalDirect[i][k] << endl;

		if (jReversibleReaction != 0)
		{
			fOutput << "   Reaction orders (backward)" << endl;
			for(k=1;k<=iGlobalInverse[i].Size();k++)
				fOutput << "     "	<< setw(20) << left << reactionRates->names[iGlobalInverse[i][k]] 
									<< setw(6) << setprecision(3) << left << lambdaGlobalInverse[i][k] << endl;
		}
	}
}

double OpenSMOKE_Kinetics::GiveMe_A(const int j)
{
	return exp_k01[j];
}

double OpenSMOKE_Kinetics::GiveMe_E(const int j)
{
	return -E1[j]*Constants::R_cal_mol * Constants::Conversion_cal_J;
}

double OpenSMOKE_Kinetics::GiveMe_Beta(const int j)
{
	return beta1[j];
}

double OpenSMOKE_Kinetics::GiveMe_ForwardOrder(const int j)
{
	return forwardOrders[j];
}

void OpenSMOKE_Kinetics::UpdateKineticParameters(const int iModel, solid_regression_parameters &parameters, BzzVector &b, BzzMatrixInt &user)
{
	int j;

	if (parameters == SOLID_REGRESSION_ALL)
	{
		if (iModel == 0)
		{
			for(j=1;j<=NR;j++)	{ exp_k01[j] = b[j]; k01[j] = log(exp_k01[j]); }
			for(j=1;j<=NR;j++)	beta1[j] = b[NR+j];
			for(j=1;j<=NR;j++)	E1[j]	= -b[2*NR+j]/Constants::R_cal_mol/Constants::Conversion_cal_J;
		}
		else if (iModel == 1)
		{
			for(j=1;j<=NR;j++)	{ exp_k01[j] = exp(b[j]); k01[j] = log(exp_k01[j]); }
			for(j=1;j<=NR;j++)	beta1[j] = b[NR+j];
			for(j=1;j<=NR;j++)	E1[j]	= -(b[2*NR+j]*Constants::R_J_kmol)/Constants::R_cal_mol/Constants::Conversion_cal_J;
		}
		else if (iModel == 2)
		{
			ErrorMessage("Model 2 not yet implemented...");
		}
	}
	
	else if (parameters == SOLID_REGRESSION_A)
	{
		if (iModel == 0)		for(j=1;j<=NR;j++)	{ exp_k01[j] = b[j]; k01[j] = log(exp_k01[j]); }
		else if (iModel == 1)	for(j=1;j<=NR;j++)	{ exp_k01[j] = exp(b[j]); k01[j] = log(exp_k01[j]); }
		else if (iModel == 2)	ErrorMessage("Model 2 not yet implemented...");
	}
	
	else if (parameters == SOLID_REGRESSION_BETA)
	{
		if (iModel == 0)		for(j=1;j<=NR;j++)	beta1[j] = b[j];
		else if (iModel == 1)	for(j=1;j<=NR;j++)	beta1[j] = b[j];
		else if (iModel == 2)	ErrorMessage("Model 2 not yet implemented...");
	}
	
	else if (parameters == SOLID_REGRESSION_E)
	{
		if (iModel == 0)		for(j=1;j<=NR;j++)	E1[j]	= -b[j]/Constants::R_cal_mol/Constants::Conversion_cal_J;
		else if (iModel == 1)	for(j=1;j<=NR;j++)	E1[j]	= -(b[j]*Constants::R_J_kmol)/Constants::R_cal_mol/Constants::Conversion_cal_J;
		else if (iModel == 2)	ErrorMessage("Model 2 not yet implemented...");
	}

	else if (parameters == SOLID_REGRESSION_A_E)
	{
		if (iModel == 0)
		{
			for(j=1;j<=NR;j++)	{ exp_k01[j] = b[j]; k01[j] = log(exp_k01[j]); }
			for(j=1;j<=NR;j++)	E1[j] = -b[NR+j]/Constants::R_cal_mol/Constants::Conversion_cal_J;
		}
		else if (iModel == 1)
		{
			for(j=1;j<=NR;j++)	{ exp_k01[j] = exp(b[j]); k01[j] = log(exp_k01[j]); }
			for(j=1;j<=NR;j++)	E1[j] = -(b[NR+j]*Constants::R_J_kmol)/Constants::R_cal_mol/Constants::Conversion_cal_J;
		}
		else if (iModel == 2)
		{
			ErrorMessage("Model 2 not yet implemented...");
		}
	}

	else if (parameters == SOLID_REGRESSION_A_BETA)
	{
		if (iModel == 0)
		{
			for(j=1;j<=NR;j++)	{ exp_k01[j] = b[j]; k01[j] = log(exp_k01[j]); }
			for(j=1;j<=NR;j++)	beta1[j] = b[NR+j];
		}
		else if (iModel == 1)
		{
			for(j=1;j<=NR;j++)	{ exp_k01[j] = exp(b[j]); k01[j] = log(exp_k01[j]); }
			for(j=1;j<=NR;j++)	beta1[j] = b[NR+j];
		}
		else if (iModel == 2)
		{
			ErrorMessage("Model 2 not yet implemented...");
		}
	}
	
	else if (parameters == SOLID_REGRESSION_BETA_E)
	{
		if (iModel == 0)
		{
			for(j=1;j<=NR;j++)	beta1[j] = b[j];
			for(j=1;j<=NR;j++)	E1[j]	= -b[NR+j]/Constants::R_cal_mol/Constants::Conversion_cal_J;
		}
		else if (iModel == 1)
		{
			for(j=1;j<=NR;j++)	beta1[j] = b[j];
			for(j=1;j<=NR;j++)	E1[j]	= -(b[NR+j]*Constants::R_J_kmol)/Constants::R_cal_mol/Constants::Conversion_cal_J;
		}
		else if (iModel == 2)
		{
			ErrorMessage("Model 2 not yet implemented...");
		}
	}

	else if (parameters == SOLID_REGRESSION_USERDEFINED)
	{
		for(j=1;j<=b.Size();j++)
		{
			int jReaction = user[j][1];

			if (iModel == 0)
			{
				     if (user[j][2] == 1)	{ exp_k01[jReaction] = b[j]; k01[jReaction] = log(exp_k01[jReaction]); }
				else if (user[j][2] == 2)	beta1[jReaction] = b[j];
				else if (user[j][2] == 3)	E1[jReaction]	= -b[j]/Constants::R_cal_mol/Constants::Conversion_cal_J;
				else ErrorMessage("Only A || E || Beta can be optimized...");
			}
			else if (iModel == 1)
			{
				     if (user[j][2] == 1)	{ exp_k01[jReaction] = exp(b[j]); k01[jReaction] = log(exp_k01[jReaction]); }
				else if (user[j][2] == 2)	beta1[jReaction] = b[j];
				else if (user[j][2] == 3)	E1[jReaction]	= -(b[j]*Constants::R_J_kmol)/Constants::R_cal_mol/Constants::Conversion_cal_J;
				else ErrorMessage("Only A || E || Beta can be optimized...");
			}
			else if (iModel == 2)
			{
				ErrorMessage("Model 2 not yet implemented...");
			}
		}
	}
}

void OpenSMOKE_Kinetics::UpdateOptimizerParameters(const int iModel, solid_regression_parameters &parameters, BzzVector &b, BzzMatrixInt &user)
{
	int j;

	// First guess values
	if (parameters == SOLID_REGRESSION_ALL)
	{
		if (iModel == 0)
		{
			for(j=1;j<=NR;j++)	b[j]						= exp_k01[j];
			for(j=1;j<=NR;j++)	b[NR+j]	= beta1[j];		
			for(j=1;j<=NR;j++)	b[2*NR+j]	= -E1[j]*Constants::Conversion_cal_J*Constants::R_cal_mol;		
		}
		else if (iModel == 1)
		{
			for(j=1;j<=NR;j++)	b[j]						= log(exp_k01[j]);
			for(j=1;j<=NR;j++)	b[NR+j]	= beta1[j];
			for(j=1;j<=NR;j++)	b[2*NR+j]	= -(E1[j]*Constants::Conversion_cal_J/Constants::R_J_kmol)*Constants::R_cal_mol;
		}
		else if (iModel == 2)
		{
			ErrorMessage("Model 2 not yet implemented...");
		}
	}
	
	else if (parameters == SOLID_REGRESSION_A)
	{
		if (iModel == 0)		for(j=1;j<=NR;j++)	b[j]	= exp_k01[j];
		else if (iModel == 1)	for(j=1;j<=NR;j++)	b[j]	= log(exp_k01[j]);
		else if (iModel == 2)	ErrorMessage("Model 2 not yet implemented...");
	}

	else if (parameters == SOLID_REGRESSION_BETA)
	{
		if (iModel == 0)		for(j=1;j<=NR;j++)	b[j]	= beta1[j];
		else if (iModel == 1)	for(j=1;j<=NR;j++)	b[j]	= beta1[j];
		else if (iModel == 2)	ErrorMessage("Model 2 not yet implemented...");
	}

	else if (parameters == SOLID_REGRESSION_E)
	{
		if (iModel == 0)		for(j=1;j<=NR;j++)	b[j]	= -E1[j]*Constants::Conversion_cal_J*Constants::R_cal_mol;			
		else if (iModel == 1)	for(j=1;j<=NR;j++)	b[j]	= -E1[j]*Constants::Conversion_cal_J/Constants::R_J_kmol*Constants::R_cal_mol;
		else if (iModel == 2)	ErrorMessage("Model 2 not yet implemented...");
	}

	else if (parameters == SOLID_REGRESSION_A_E)
	{
		if (iModel == 0)
		{
			for(j=1;j<=NR;j++)	b[j]						= exp_k01[j];
			for(j=1;j<=NR;j++)	b[NR+j]	= -E1[j]*Constants::Conversion_cal_J*Constants::R_cal_mol;	
		}
		else if (iModel == 1)
		{
			for(j=1;j<=NR;j++)	b[j]						= log(exp_k01[j]);
			for(j=1;j<=NR;j++)	b[NR+j]	= -E1[j]*Constants::Conversion_cal_J/Constants::R_J_kmol*Constants::R_cal_mol;
		}
		else if (iModel == 2)
		{
			ErrorMessage("Model 2 not yet implemented...");
		}
	}
	else if (parameters == SOLID_REGRESSION_A_BETA)
	{
		if (iModel == 0)
		{
			for(j=1;j<=NR;j++)	b[j]						= exp_k01[j];
			for(j=1;j<=NR;j++)	b[NR+j]	= beta1[j];		
		}
		else if (iModel == 1)
		{
			for(j=1;j<=NR;j++)	b[j]						= log(exp_k01[j]);
			for(j=1;j<=NR;j++)	b[NR+j]	= beta1[j];
		}
		else if (iModel == 2)
		{
			ErrorMessage("Model 2 not yet implemented...");
		}
	}

	else if (parameters == SOLID_REGRESSION_BETA_E)
	{
		if (iModel == 0)
		{
			for(j=1;j<=NR;j++)	b[j]						= beta1[j];		
			for(j=1;j<=NR;j++)	b[NR+j]	= -E1[j]*Constants::Conversion_cal_J*Constants::R_cal_mol;			
		}
		else if (iModel == 1)
		{
			for(j=1;j<=NR;j++)	b[j]						= beta1[j];
			for(j=1;j<=NR;j++)	b[NR+j]	= -E1[j]*Constants::Conversion_cal_J/Constants::R_J_kmol*Constants::R_cal_mol;
		}
		else if (iModel == 2)
		{
			ErrorMessage("Model 2 not yet implemented...");
		}
	}

	else if (parameters == SOLID_REGRESSION_USERDEFINED)
	{
		for(j=1;j<=b.Size();j++)
		{
			int jReaction = user[j][1];

			if (iModel == 0)
			{
				     if (user[j][2] == 1)	b[j] = exp_k01[jReaction];
				else if (user[j][2] == 2)	b[j] = beta1[jReaction];	
				else if (user[j][2] == 3)	b[j] = -E1[jReaction]*Constants::Conversion_cal_J*Constants::R_cal_mol;	
				else ErrorMessage("Only A || E || Beta can be optimized...");
			}
			else if (iModel == 1)
			{
				     if (user[j][2] == 1)	b[j] = log(exp_k01[jReaction]);
				else if (user[j][2] == 2)	b[j] = beta1[jReaction];	
				else if (user[j][2] == 3)	b[j] = -E1[jReaction]*Constants::Conversion_cal_J/Constants::R_J_kmol*Constants::R_cal_mol;
				else ErrorMessage("Only A || E || Beta can be optimized...");
			}
			else if (iModel == 2)
			{
				ErrorMessage("Model 2 not yet implemented...");
			}
		}
	}
}


BzzVectorInt OpenSMOKE_Kinetics::ReactionIndices(const vector<string> list_of_species)
{
	int* jD1 = jDir1.GetHandle();
	int* jD2 = jDir2.GetHandle();
	int* jD3 = jDir3.GetHandle();
	int* jD4 = jDir4.GetHandle();
	int* jD5 = jDir5.GetHandle();
	double* vD5 = valDir5.GetHandle();

	int* jIT1 = jInvTot1.GetHandle();
	int* jIT2 = jInvTot2.GetHandle();
	int* jIT3 = jInvTot3.GetHandle();
	int* jIT4 = jInvTot4.GetHandle();
	int* jIT5 = jInvTot5.GetHandle();
	double* vIT5 = valInvTot5.GetHandle();

	BzzMatrix gamma(NC, NR);

	int i,j,k;
	for(i = 1;i <= NC;i++)
	{
		for(k = 1;k <= numDir1[i];k++)
			gamma[i][*jD1++] += -1.;
		for(k = 1;k <= numDir2[i];k++)
			gamma[i][*jD2++] += -2.;
		for(k = 1;k <= numDir3[i];k++)
			gamma[i][*jD3++] += -3.;
		for(k = 1;k <= numDir4[i];k++)
			gamma[i][*jD4++] += -0.5;
		for(k = 1;k <= numDir5[i];k++)
			gamma[i][*jD5++] += -(*vD5++);

		for(k = 1;k <= numInvTot1[i];k++)
			gamma[i][*jIT1++] += 1.;
		for(k = 1;k <= numInvTot2[i];k++)
			gamma[i][*jIT2++] += 2.;
		for(k = 1;k <= numInvTot3[i];k++)
			gamma[i][*jIT3++] += 3.;
		for(k = 1;k <= numInvTot4[i];k++)
			gamma[i][*jIT4++] += 0.5;
		for(k = 1;k <= numInvTot5[i];k++)
			gamma[i][*jIT5++] += (*vIT5++);
	}

	BzzVectorInt indices;
	for(i=1;i<=int(list_of_species.size());i++)
	{
		int index = reactionRates->recognize_species(list_of_species[i-1]);
		for(j=1;j<=NR;j++)
			if(gamma[index][j] != 0.)
			{
				bool iPresent = false;
				for(int k=1;k<=indices.Size();k++)
					if (indices[k] == j)
					{	
						iPresent = true;
						break;
					}
				if (iPresent == false)
					indices.Append(j);
			}
	}

	return indices;
}

void OpenSMOKE_Kinetics::ChangeFrequencyFactor(const int j, const double variation)
{
	k01[j]  = exp(k01[j]);
	k01[j] *= (1.+variation);
	exp_k01[j] = k01[j];
	k01[j] = log(k01[j]);
}

void OpenSMOKE_Kinetics::ChangeActivationEnergy(const int j, const double new_value_cal_mol)
{
	cout << "Change activation energy for reaction #" << j << endl;
	cout << " " << reactionRates->strReaction[j] << endl;
	cout << " Original value: " << -E1[j]*Constants::R_cal_mol << " cal/mol" <<  endl;
	cout << " Modified value: " << new_value_cal_mol << " cal/mol" << endl;

	E1[j] = -new_value_cal_mol /  Constants::R_cal_mol;
}

void OpenSMOKE_Kinetics::StoichiometricMatrix(BzzMatrix &gamma)
{
	ChangeDimensions(NC, NR, &gamma);

	int i,k;

	int* jD1 = jDir1.GetHandle();
	int* jD2 = jDir2.GetHandle();
	int* jD3 = jDir3.GetHandle();
	int* jD4 = jDir4.GetHandle();
	int* jD5 = jDir5.GetHandle();
	double* vD5 = valDir5.GetHandle();

	int* jIT1 = jInvTot1.GetHandle();
	int* jIT2 = jInvTot2.GetHandle();
	int* jIT3 = jInvTot3.GetHandle();
	int* jIT4 = jInvTot4.GetHandle();
	int* jIT5 = jInvTot5.GetHandle();
	double* vIT5 = valInvTot5.GetHandle();
	
	for(i = 1;i <= NC;i++)
	{
		for(k = 1;k <= numDir1[i];k++)
			gamma[i][*jD1++] -= 1.;
		for(k = 1;k <= numDir2[i];k++)
			gamma[i][*jD2++] -= 2.;
		for(k = 1;k <= numDir3[i];k++)
			gamma[i][*jD3++] -= 3.;
		for(k = 1;k <= numDir4[i];k++)
			gamma[i][*jD4++] -= .5;
		for(k = 1;k <= numDir5[i];k++)
			gamma[i][*jD5++] -= (*vD5++);

		for(k = 1;k <= numInvTot1[i];k++)
			gamma[i][*jIT1++] += 1.;
		for(k = 1;k <= numInvTot2[i];k++)
			gamma[i][*jIT2++] += 2.;
		for(k = 1;k <= numInvTot3[i];k++)
			gamma[i][*jIT3++] += 3.;
		for(k = 1;k <= numInvTot4[i];k++)
			gamma[i][*jIT4++] += .5;
		for(k = 1;k <= numInvTot5[i];k++)
			gamma[i][*jIT5++] += (*vIT5++);
	}
}

BzzVector OpenSMOKE_Kinetics::GiveMeSumNuDirect()
{
	BzzVector sumNuDirect(NR);

	int i, k;

	int *jD1 = jDir1.GetHandle();
	int *jD2 = jDir2.GetHandle();
	int *jD3 = jDir3.GetHandle();
	int *jD4 = jDir4.GetHandle();
	int *jD5 = jDir5.GetHandle();
	double *vD5 = valDir5.GetHandle();
			
	for(i = 1;i <= NC;i++)
	{
		for(k = 1;k <= numDir1[i];k++)
		{
			sumNuDirect[*jD1] += 1.;
			jD1++;
		}
		for(k = 1;k <= numDir2[i];k++)
		{
			sumNuDirect[*jD2] += 2.;
			jD2++;
		}
		for(k = 1;k <= numDir3[i];k++)
		{
			sumNuDirect[*jD3] += 3.;
			jD3++;
		}
		for(k = 1;k <= numDir4[i];k++)
		{
			sumNuDirect[*jD4] += 0.5;
			jD4++;
		}
		for(k = 1;k <= numDir5[i];k++)
		{
			sumNuDirect[*jD5] += *vD5;
			jD5++;
			vD5++;
		}
	}

	return sumNuDirect;
}