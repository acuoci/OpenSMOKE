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


const int		OpenSMOKE_Kinetics::NUM_MAX					= 20000;
const double	OpenSMOKE_Kinetics::S1						= 1.001;	
const double	OpenSMOKE_Kinetics::S2						= S1 * S1;	
const double	OpenSMOKE_Kinetics::S3						= S2 * S1;	
const double	OpenSMOKE_Kinetics::SSQ					= sqrt(S1);	

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
	int i,k;
	double c1,c2,c3,csq;

	int *jD1, *jD2, *jD3, *jD4, *jD5;
	double *vD5;

	int *jIE1, *jIE2, *jIE3, *jIE4, *jIE5;
	double *vIE5;

	jD1 = jDir1.GetHandle();
	jD2 = jDir2.GetHandle();
	jD3 = jDir3.GetHandle();
	jD4 = jDir4.GetHandle();
	jD5 = jDir5.GetHandle();
	vD5 = valDir5.GetHandle();

	jIE1 = jInvEq1.GetHandle();
	jIE2 = jInvEq2.GetHandle();
	jIE3 = jInvEq3.GetHandle();
	jIE4 = jInvEq4.GetHandle();
	jIE5 = jInvEq5.GetHandle();
	vIE5 = valInvEq5.GetHandle();

	rDirC = 1.;
	rInvC = -1.;
	
	for(i = 1;i <= NC;i++)
	{
		c1 = c[i];
		c2 = c1 * c1;
		c3 = c2 * c1;
		if(numDir4[i] != 0)
			csq = sqrt(c1);
		for(k = 1;k <= numDir1[i];k++)
			{
			rDirC[*jD1] *= c1;
			jD1++;
			}
		for(k = 1;k <= numDir2[i];k++)
			{
			rDirC[*jD2] *= c2;
			jD2++;
			}
		for(k = 1;k <= numDir3[i];k++)
			{
			rDirC[*jD3] *= c3;
			jD3++;
			}
		for(k = 1;k <= numDir4[i];k++)
			{
			rDirC[*jD4] *= csq;
			jD4++;
			}
		for(k = 1;k <= numDir5[i];k++)
			{
			rDirC[*jD5] *= pow(c1,*vD5);
			jD5++;
			vD5++;
			}

		for(k = 1;k <= numInvEq1[i];k++)
			{
			rInvC[*jIE1] *= c1;
			jIE1++;
			}
		for(k = 1;k <= numInvEq2[i];k++)
			{
			rInvC[*jIE2] *= c2;
			jIE2++;
			}
		for(k = 1;k <= numInvEq3[i];k++)
			{
			rInvC[*jIE3] *= c3;
			jIE3++;
			}
		for(k = 1;k <= numInvEq4[i];k++)
			{
			rInvC[*jIE4] *= csq;
			jIE4++;
			}
		for(k = 1;k <= numInvEq5[i];k++)
			{
			rInvC[*jIE5] *= pow(c1,*vIE5);
			jIE5++;
			vIE5++;
			}
		}



/*
	rDirC = 1.;
	rInvC = -1.;
	
	for(i = 1;i <= NC;i++)
	{
		c1 = c[i];
		c2 = c1 * c1;
		c3 = c2 * c1;
		if(numDir4[i] != 0 || numInvEq4[i] != 0)
			csq = sqrt(c1);
		
		for(k = 1;k <= numDir1[i];k++)
		{
			rDirC[*jD1] *= c1;
			jD1++;
		}
		for(k = 1;k <= numDir2[i];k++)
		{
			rDirC[*jD2] *= c2;
			jD2++;
		}
		for(k = 1;k <= numDir3[i];k++)
		{
			rDirC[*jD3] *= c3;
			jD3++;
		}
		for(k = 1;k <= numDir4[i];k++)
		{
			rDirC[*jD4] *= csq;
			jD4++;
		}
		for(k = 1;k <= numDir5[i];k++)
		{
			rDirC[*jD5] *= pow(c1,*vD5);
			jD5++;
			vD5++;
		}

		for(k = 1;k <= numInvEq1[i];k++)
		{
			rInvC[*jIE1] *= c1;
			jIE1++;
		}
		for(k = 1;k <= numInvEq2[i];k++)
		{
			rInvC[*jIE2] *= c2;
			jIE2++;
		}
		for(k = 1;k <= numInvEq3[i];k++)
		{
			rInvC[*jIE3] *= c3;
			jIE3++;
		}
		for(k = 1;k <= numInvEq4[i];k++)
		{
			rInvC[*jIE4] *= csq;
			jIE4++;
		}
		for(k = 1;k <= numInvEq5[i];k++)
		{
				rInvC[*jIE5] *= pow(c1,*vIE5);
			jIE5++;
			vIE5++;
		}
	}
*/
}

//----------------------------------------------------------------------------------//
//																					//
//						          computeFallOff									//
//																					//
//----------------------------------------------------------------------------------//	

void OpenSMOKE_Kinetics::fallOff(double T, BzzVector &coeffM, BzzVector &k1, BzzVector &k2, BzzVector &logFcent, BzzVector &coeffFallOff)
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
			case 4:
				logPr = log10(Pr);
				xSRI = 1. / (1. + BzzPow2(logPr));
				wF = pow(logFcent[j],xSRI) * dFall[j];
				if(eFall[j] != 0.)
					wF *= pow(T,eFall[j]);
				break;
			default: // Lindeman
				wF = 1.;
				break;
			}
	
		coeffFallOff[j] = k2[j] * (Pr / (1. + Pr)) * wF;
	}
}



//----------------------------------------------------------------------------------//
//																					//
//						          computeThirdBody									//
//																					//
//----------------------------------------------------------------------------------//	

// Calcola il coefficiente correttivo Gamma per tutte quelle reazioni che prevedono
// la presenza di un terzo corpo (anche quelle di FallOff)
void OpenSMOKE_Kinetics::thirdBody(BzzVector &c, double cTot, BzzVector &coeffM)
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

//----------------------------------------------------------------------------------//
//																					//
//							 reactionsWithFallOff									//
//																					//
//----------------------------------------------------------------------------------//	
/*
void OpenSMOKE_Kinetics::reactionsWithFallOff(BzzVector &coeffFallOff,BzzVector &r)
{
	int k, j;

	for(k = 1;k <= numFallOff;k++)
	{
		j=iFallOff[k];
		r[j]*=coeffFallOff[j];
	}
}
*/
//----------------------------------------------------------------------------------//
//																					//
//							 compositionReactionRates								//
//																					//
//----------------------------------------------------------------------------------//	

void OpenSMOKE_Kinetics::compositionReactionRates(BzzVector &r, BzzVector *R)
{
	int i,k;
/*
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
			{(*R)[i] -= (r[*jD2] + r[*jD2]);*jD2++;}
		for(k = 1;k <= numDir3[i];k++)
			{(*R)[i] -= (r[*jD3] + r[*jD3] + r[*jD3]); *jD3++;}
		for(k = 1;k <= numDir4[i];k++)
			(*R)[i] -= 0.5 * r[*jD4++];
		for(k = 1;k <= numDir5[i];k++)
			(*R)[i] -= (*vD5++) * r[*jD5++];

		for(k = 1;k <= numInvTot1[i];k++)
			(*R)[i] += r[*jIT1++];
		for(k = 1;k <= numInvTot2[i];k++)
			{(*R)[i] += (r[*jIT2] + r[*jIT2]);*jIT2++;}
		for(k = 1;k <= numInvTot3[i];k++)
			{(*R)[i] += (r[*jIT3] + r[*jIT3] + r[*jIT3]); *jIT3++;}
		for(k = 1;k <= numInvTot4[i];k++)
			(*R)[i] += 0.5 * r[*jIT4++];
		for(k = 1;k <= numInvTot5[i];k++)
			(*R)[i] += (*vIT5++) * r[*jIT5++];
	}
*/
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

void OpenSMOKE_Kinetics::ComputeKineticParameters(BzzVector &reactionDH, BzzVector &reactionDS, double T, double logT, double uT, double loguRT, BzzVector &uKeq, BzzVector &k1, BzzVector &k2, BzzVector &logFcent)
{
	int k,j;

	for(k = 1;k <= numEquilibrium;k++)
	{
		j = reactionWithEquil[k];
		uKeq[j] = exp(-reactionDS[j] + reactionDH[j] - loguRT * sumNuij[j]);
	}

	for(j = 1;j <= NR;j++)
		k1[j] = exp(k01[j] + beta1[j] * logT + E1[j] * uT);

	
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
				if(bFall[j] != 0.)
					e1Troe = exp(-T / bFall[j]);
				if(cFall[j] != 0.)
					e2Troe = exp(-T / cFall[j]);
				if(dFall[j] != 0.)
					e3Troe = exp(-dFall[j] / T);
				Fcent = (1. - aFall[j]) * e1Troe + aFall[j] * e2Troe + e3Troe;
				if(Fcent < 1.e-300)
						logFcent[j] = -300.;
				else
					logFcent[j] = log10(Fcent);
				break;
			case 4:
				e1Troe = 1.;
				e2Troe = 0.;
				if(bFall[j] != 0.)
					e1Troe = exp(-bFall[j] / T);
				if(cFall[j] != 0.)
					e2Troe = exp(-T / cFall[j]);
				logFcent[j] = aFall[j] * e1Troe + e2Troe;
				break;
			default: // Lindeman
				break;
		}
	}
	
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



void OpenSMOKE_Kinetics::GetDerivativesC(double T, double cTot, BzzMatrix *dRC, 
									   BzzVector &cRes, BzzVector &R)
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
		// Solo se la concentrazione ?praticamente nulla
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
		thirdBody(cRes, cTot, reactionRates->coeffM);

		// 3. CONTRIBUTI DERIVANTI DALLE REAZIONI DI FALLOFF
		fallOff(T, reactionRates->coeffM, reactionRates->k1, reactionRates->k2,
		       reactionRates->logFcent, reactionRates->coeffFallOff);
  
		// 4. ASEMBLAGGIO DEI DIVERSI CONTRIBUTI
		reactionsWithEquilibrium(reactionRates->rDirC, reactionRates->rInvC, reactionRates->uKeq, reactionRates->r);
		reactionsWithThirdBody(reactionRates->coeffM, reactionRates->r);
		//reactionsWithFallOff(reactionRates->coeffFallOff, reactionRates->r);
		
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
	numThirdBodyOnly=0;
	numFallOff=0;
	numThirdBody=0;
	numConventional=0;

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
	ChangeDimensions(NR,&aFall);
	ChangeDimensions(NR,&bFall);
	ChangeDimensions(NR,&cFall);
	ChangeDimensions(NR,&dFall);
	ChangeDimensions(NR,&eFall);
	ChangeDimensions(NUM_MAX,&iFallOff);
	ChangeDimensions(NUM_MAX,&iThirdBody);
	ChangeDimensions(NUM_MAX,&iThirdBodyOnly);


	// 1. Reading reaction data
	// -------------------------------------------------------------------------------
	BzzLoad inputFile;
	inputFile('*', fileKin);

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
		else if(jNumEfficiency[j]  > 0)
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
		if(jThirdBody[j]==0 && jNumEfficiency[j]==0) numConventional++;
		else
		{
			// - tutte le reazioni di FallOff sono automaticamente anche di terzo corpo
			//   e dunque si rende comunque necessario il calcolo del coefficiente correttivo coeffM

			// la reazione e' solo di terzo corpo
			if(jThirdBody[j]==1)
			{
				numThirdBodyOnly++;
				iThirdBodyOnly[numThirdBodyOnly]=j;
			}

			// reazioni con terzo corpo (sia quelle pure che di falloff)
			if(jThirdBody[j]!=0) 
			{  
				numThirdBody++;
				iThirdBody[numThirdBody]=j;
			}

			// in questo caso la reazione e' sicuramente di FallOff e quindi in realta anche di terzo corpo
			if(jThirdBody[j]>=2) 
			{
				numFallOff++;
				iFallOff[numFallOff]=j;
			}

		}

		// Reading Kinetic data
		inputFile  >> k01[j] >> beta1[j] >> E1[j];
		exp_k01[j] = k01[j];
		k01[j] = log(k01[j]);
		E1[j] = -E1[j] /  Constants::R_cal_mol;

		// In case of Fall-Off reactions
		if(jThirdBody[j] >= 2)
		{
			inputFile  >> k02[j] >> beta2[j] >> E2[j];
			k02[j] = log(k02[j]);
			E2[j] = -E2[j] /  Constants::R_cal_mol;
					
			inputFile  >> aFall[j] >> bFall[j] >> cFall[j] >> dFall[j] >> eFall[j];
		}
	}

	// 2. Cleaning matrices
	// -------------------------------------------------------------------------------
	iFallOff.DeleteLastNElements(NUM_MAX-numFallOff);
	iThirdBody.DeleteLastNElements(NUM_MAX-numThirdBody);
	iThirdBodyOnly.DeleteLastNElements(NUM_MAX-numThirdBodyOnly);


	// 3. Equilibrium reactions
	// -------------------------------------------------------------------------------
	numEquilibrium = 0;
	for(j = 1;j <= NR;j++)
	{
		if(jEquil[j] == 1) numEquilibrium++;
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
		char name[Constants::REACTION_NAME_SIZE];
		inputFile.fileLoad.read((char*) name, sizeof(name));
		strcpy(reactionRates->strReaction[j], name);
		// cout << j << " " << reactionRates->strReaction[j] << endl;
	}
	
	// 6. End
	// -------------------------------------------------------------------------------
	inputFile.End();

	
	// 5. Riepilogo delle informazioni a video
	// -------------------------------------------------------------------------------
	cout << "--------------------------------------------------------------------------" << endl;
	cout << "                       KINETIC SCHEME SUMMARY                             " << endl;
	cout << "--------------------------------------------------------------------------" << endl;
	cout << "  Total number of species                   = " << NC			 << endl;
	cout << "  Total number of reactions                 = " << NR			 << endl;
	cout << "    Number of irreversible reactions        = " << numIrreversible				 << endl;
	cout << "    Number of equilibrium reactions         = " << numEquilibrium				 << endl;
	cout << "    Number of conventional reactions        = " << numConventional		 << endl;
	cout << "    Number of third body reactions (total)  = " << numThirdBody		 	 << endl;
	cout << "    Number of third body reactions          = " << numThirdBodyOnly		 << endl;
	cout << "    Number of Fall-Off reactions            = " << numFallOff				 << endl;
	cout << "--------------------------------------------------------------------------" << endl;
	cout << endl;
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