/***************************************************************************
 *   Copyright (C) 2003-2008 by                                            *
 *   Alberto Cuoci		                                                   *
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

#include "symbolickinetics/polimi_h2conox_1101/OpenSMOKE_SymbolicKinetics_Polimi_H2CONOX_1101.h"

void OpenSMOKE_SymbolicKinetics_Polimi_H2CONOX_1101::giveJacobian(BzzVector &c, BzzMatrix &J) 
{

	// ============================================================ 
	// ===== DERIVATIVES FOR FALLOFF REACTIONS ==================== 
	// ============================================================ 
	dCFOdM3 =	coeffFallOff3*coeffFallOff3/BzzPow2(coeffFallOff3+k3*coeffM3); 
	dCFOdM14 =	coeffFallOff14*coeffFallOff14/BzzPow2(coeffFallOff14+k14*coeffM14); 
	dCFOdM18 =	coeffFallOff18*coeffFallOff18/BzzPow2(coeffFallOff18+k18*coeffM18); 
	dCFOdM45 =	coeffFallOff45*coeffFallOff45/BzzPow2(coeffFallOff45+k45*coeffM45); 
	dCFOdM107 =	coeffFallOff107*coeffFallOff107/BzzPow2(coeffFallOff107+k107*coeffM107); 
	dCFOdM161 =	coeffFallOff161*coeffFallOff161/BzzPow2(coeffFallOff161+k161*coeffM161); 


	if (iAccurateJacobian == 1)
	{
		double num3 = -0.40  + 0.43429448*lnPr3-0.670000*logFcent3; 
		double den3 =  0.806 - 0.0608012 *lnPr3-1.176199*logFcent3; 
		dwFdM3 = -2.30258509299 * wF3 * logFcent3 * num3/(coeffM3 * den3*den3) * (0.86858896+0.12160240*num3/den3) / BzzPow2(1+BzzPow2(num3/den3)); 

		double num14 = -0.40  + 0.43429448*lnPr14-0.670000*logFcent14; 
		double den14 =  0.806 - 0.0608012 *lnPr14-1.176199*logFcent14; 
		dwFdM14 = -2.30258509299 * wF14 * logFcent14 * num14/(coeffM14 * den14*den14) * (0.86858896+0.12160240*num14/den14) / BzzPow2(1+BzzPow2(num14/den14)); 

		dwFdM18 = 0.; 
		dwFdM45 = 0.; 
		dwFdM107 = 0.; 
		double num161 = -0.40  + 0.43429448*lnPr161-0.670000*logFcent161; 
		double den161 =  0.806 - 0.0608012 *lnPr161-1.176199*logFcent161; 
		dwFdM161 = -2.30258509299 * wF161 * logFcent161 * num161/(coeffM161 * den161*den161) * (0.86858896+0.12160240*num161/den161) / BzzPow2(1+BzzPow2(num161/den161)); 

	}
	else
	{
		dwFdM3 = 0.; 
		dwFdM14 = 0.; 
		dwFdM18 = 0.; 
		dwFdM45 = 0.; 
		dwFdM107 = 0.; 
		dwFdM161 = 0.; 
	}


	// ============================================================ 
	// ===== DERIVATIVES FOR EVERY REACTION ======================= 
	// ============================================================ 
	d1d6 = 	   c[24]*k1;
	d1d22 = 	   -c[27]*k1*uK1;
	d1d24 = 	   c[6]*k1;
	d1d27 = 	   -c[22]*k1*uK1;

	d2d13 = 	   c[22]*k2;
	d2d22 = 	   c[13]*k2;
	d2d24 = 	   -c[27]*k2*uK2;
	d2d27 = 	   -c[24]*k2*uK2;

	sigma3 =	   CFO3*dwFdM3+dCFOdM3*wF3;
	d3d1 = 	   rFlat3*sigma3*8.0E-1;
	d3d2 = 	   rFlat3*sigma3*8.0E-1;
	d3d3 = 	   rFlat3*sigma3*1.26;
	d3d4 = 	   rFlat3*sigma3*1.0;
	d3d5 = 	   rFlat3*sigma3*1.0;
	d3d6 = 	   c[24]*coeffFallOff3*k3;
	d3d7 = 	   rFlat3*sigma3*1.0;
	d3d8 = 	   rFlat3*sigma3*1.0;
	d3d9 = 	   rFlat3*sigma3*1.0;
	d3d10 = 	   rFlat3*sigma3*1.0;
	d3d11 = 	   rFlat3*sigma3*1.0;
	d3d12 = 	   rFlat3*sigma3*1.0;
	d3d13 = 	   rFlat3*sigma3*2.5;
	d3d14 = 	   rFlat3*sigma3*1.0;
	d3d15 = 	   rFlat3*sigma3*1.8E1;
	d3d16 = 	   rFlat3*sigma3*1.0;
	d3d17 = 	   rFlat3*sigma3*1.0;
	d3d18 = 	   rFlat3*sigma3*1.0;
	d3d19 = 	   rFlat3*sigma3*1.2;
	d3d20 = 	   rFlat3*sigma3*2.4;
	d3d21 = 	   rFlat3*sigma3*1.0;
	d3d22 = 	   rFlat3*sigma3*1.0;
	d3d23 = 	   rFlat3*sigma3*1.0;
	d3d24 = 	   rFlat3*sigma3*1.0+c[6]*coeffFallOff3*k3;
	d3d25 = 	   rFlat3*sigma3*1.0;
	d3d26 = 	   rFlat3*sigma3*1.0;
	d3d27 = 	   rFlat3*sigma3*1.0;
	d3d28 = 	   rFlat3*sigma3*1.0-coeffFallOff3*k3*uK3;
	d3d29 = 	   rFlat3*sigma3*1.0;
	d3d30 = 	   rFlat3*sigma3*1.0;
	d3d31 = 	   rFlat3*sigma3*1.0;
	d3d32 = 	   rFlat3*sigma3*1.0;

	d4d6 = 	   k4*(c[24]*c[6]*2.0-c[28]*uK4);
	d4d24 = 	   (c[6]*c[6])*k4;
	d4d28 = 	   -c[6]*k4*uK4;

	d5d6 = 	   -c[15]*k5*uK5;
	d5d15 = 	   -c[6]*k5*uK5;
	d5d27 = 	   c[28]*k5;
	d5d28 = 	   c[27]*k5;

	d6d24 = 	   c[28]*k6;
	d6d27 = 	   c[27]*k6*uK6*-2.0;
	d6d28 = 	   c[24]*k6;

	d7d6 = 	   -c[27]*k7*uK7;
	d7d22 = 	   c[28]*k7;
	d7d27 = 	   -c[6]*k7*uK7;
	d7d28 = 	   c[22]*k7;

	d8d15 = 	   -c[22]*k8*uK8;
	d8d22 = 	   -c[15]*k8*uK8;
	d8d27 = 	   c[27]*k8*2.0;

	d9d1 = 	   rFlat9*5.0E-1;
	d9d2 = 	   rFlat9*5.0E-1;
	d9d3 = 	   rFlat9*1.0;
	d9d4 = 	   rFlat9*1.0;
	d9d5 = 	   rFlat9*1.0;
	d9d6 = 	   rFlat9*1.0;
	d9d7 = 	   rFlat9*1.0;
	d9d8 = 	   rFlat9*1.0;
	d9d9 = 	   rFlat9*1.0;
	d9d10 = 	   rFlat9*1.0;
	d9d11 = 	   rFlat9*1.0;
	d9d12 = 	   rFlat9*1.0;
	d9d13 = 	   rFlat9*2.5+coeffM9*k9;
	d9d14 = 	   rFlat9*1.0;
	d9d15 = 	   rFlat9*1.2E1;
	d9d16 = 	   rFlat9*1.0;
	d9d17 = 	   rFlat9*1.0;
	d9d18 = 	   rFlat9*1.0;
	d9d19 = 	   rFlat9*1.9;
	d9d20 = 	   rFlat9*3.8;
	d9d21 = 	   rFlat9*1.0;
	d9d22 = 	   rFlat9*1.0;
	d9d23 = 	   rFlat9*1.0;
	d9d24 = 	   rFlat9*1.0-c[24]*coeffM9*k9*uK9*2.0;
	d9d25 = 	   rFlat9*1.0;
	d9d26 = 	   rFlat9*1.0;
	d9d27 = 	   rFlat9*1.0;
	d9d28 = 	   rFlat9*1.0;
	d9d29 = 	   rFlat9*1.0;
	d9d30 = 	   rFlat9*1.0;
	d9d31 = 	   rFlat9*1.0;
	d9d32 = 	   rFlat9*1.0;

	d10d1 = 	   rFlat10*2.0E-1;
	d10d2 = 	   rFlat10*2.0E-1;
	d10d3 = 	   rFlat10*1.0;
	d10d4 = 	   rFlat10*1.0;
	d10d5 = 	   rFlat10*1.0;
	d10d6 = 	   rFlat10*1.0+coeffM10*k10;
	d10d7 = 	   rFlat10*1.0;
	d10d8 = 	   rFlat10*1.0;
	d10d9 = 	   rFlat10*1.0;
	d10d10 = 	   rFlat10*1.0;
	d10d11 = 	   rFlat10*1.0;
	d10d12 = 	   rFlat10*1.0;
	d10d13 = 	   rFlat10*2.5;
	d10d14 = 	   rFlat10*1.0;
	d10d15 = 	   rFlat10*1.2E1;
	d10d16 = 	   rFlat10*1.0;
	d10d17 = 	   rFlat10*1.0;
	d10d18 = 	   rFlat10*1.0;
	d10d19 = 	   rFlat10*1.9;
	d10d20 = 	   rFlat10*3.8;
	d10d21 = 	   rFlat10*1.0;
	d10d22 = 	   rFlat10*1.0-c[22]*coeffM10*k10*uK10*2.0;
	d10d23 = 	   rFlat10*1.0;
	d10d24 = 	   rFlat10*1.0;
	d10d25 = 	   rFlat10*1.0;
	d10d26 = 	   rFlat10*1.0;
	d10d27 = 	   rFlat10*1.0;
	d10d28 = 	   rFlat10*1.0;
	d10d29 = 	   rFlat10*1.0;
	d10d30 = 	   rFlat10*1.0;
	d10d31 = 	   rFlat10*1.0;
	d10d32 = 	   rFlat10*1.0;

	d11d1 = 	   rFlat11*1.0;
	d11d2 = 	   rFlat11*1.0;
	d11d3 = 	   rFlat11*1.0;
	d11d4 = 	   rFlat11*1.0;
	d11d5 = 	   rFlat11*1.0;
	d11d6 = 	   rFlat11*1.0;
	d11d7 = 	   rFlat11*1.0;
	d11d8 = 	   rFlat11*1.0;
	d11d9 = 	   rFlat11*1.0;
	d11d10 = 	   rFlat11*1.0;
	d11d11 = 	   rFlat11*1.0;
	d11d12 = 	   rFlat11*1.0;
	d11d13 = 	   rFlat11*2.0;
	d11d14 = 	   rFlat11*1.0;
	d11d15 = 	   rFlat11*1.6E1-coeffM11*k11*uK11;
	d11d16 = 	   rFlat11*1.0;
	d11d17 = 	   rFlat11*1.0;
	d11d18 = 	   rFlat11*1.0;
	d11d19 = 	   rFlat11*1.0;
	d11d20 = 	   rFlat11*1.9;
	d11d21 = 	   rFlat11*1.0;
	d11d22 = 	   rFlat11*1.0;
	d11d23 = 	   rFlat11*1.0;
	d11d24 = 	   rFlat11*1.0+c[27]*coeffM11*k11;
	d11d25 = 	   rFlat11*1.0;
	d11d26 = 	   rFlat11*1.0;
	d11d27 = 	   rFlat11*1.0+c[24]*coeffM11*k11;
	d11d28 = 	   rFlat11*1.0;
	d11d29 = 	   rFlat11*1.0;
	d11d30 = 	   rFlat11*1.0;
	d11d31 = 	   rFlat11*1.0;
	d11d32 = 	   rFlat11*1.0;

	d12d6 = 	   -c[13]*k12*uK12;
	d12d13 = 	   -c[6]*k12*uK12;
	d12d24 = 	   c[28]*k12;
	d12d28 = 	   c[24]*k12;

	d13d6 = 	   -c[16]*k13*uK13;
	d13d16 = 	   -c[6]*k13*uK13;
	d13d28 = 	   c[28]*k13*2.0;

	sigma14 =	   CFO14*dwFdM14+dCFOdM14*wF14;
	d14d1 = 	   rFlat14*sigma14*7.0E-1;
	d14d2 = 	   rFlat14*sigma14*7.0E-1;
	d14d3 = 	   rFlat14*sigma14*1.0;
	d14d4 = 	   rFlat14*sigma14*1.0;
	d14d5 = 	   rFlat14*sigma14*1.0;
	d14d6 = 	   rFlat14*sigma14*1.0;
	d14d7 = 	   rFlat14*sigma14*1.0;
	d14d8 = 	   rFlat14*sigma14*1.0;
	d14d9 = 	   rFlat14*sigma14*1.0;
	d14d10 = 	   rFlat14*sigma14*1.0;
	d14d11 = 	   rFlat14*sigma14*1.0;
	d14d12 = 	   rFlat14*sigma14*1.0;
	d14d13 = 	   rFlat14*sigma14*2.0;
	d14d14 = 	   rFlat14*sigma14*1.0;
	d14d15 = 	   rFlat14*sigma14*6.0;
	d14d16 = 	   rFlat14*sigma14*1.0-coeffFallOff14*k14*uK14;
	d14d17 = 	   rFlat14*sigma14*1.0;
	d14d18 = 	   rFlat14*sigma14*1.0;
	d14d19 = 	   rFlat14*sigma14*1.5;
	d14d20 = 	   rFlat14*sigma14*2.0;
	d14d21 = 	   rFlat14*sigma14*1.0;
	d14d22 = 	   rFlat14*sigma14*1.0;
	d14d23 = 	   rFlat14*sigma14*1.0;
	d14d24 = 	   rFlat14*sigma14*1.0;
	d14d25 = 	   rFlat14*sigma14*1.0;
	d14d26 = 	   rFlat14*sigma14*1.0;
	d14d27 = 	   rFlat14*sigma14*1.0+c[27]*coeffFallOff14*k14*2.0;
	d14d28 = 	   rFlat14*sigma14*1.0;
	d14d29 = 	   rFlat14*sigma14*1.0;
	d14d30 = 	   rFlat14*sigma14*1.0;
	d14d31 = 	   rFlat14*sigma14*1.0;
	d14d32 = 	   rFlat14*sigma14*1.0;

	d15d1 = 	   rFlat15;
	d15d2 = 	   rFlat15;
	d15d3 = 	   rFlat15;
	d15d4 = 	   rFlat15;
	d15d5 = 	   rFlat15;
	d15d6 = 	   rFlat15;
	d15d7 = 	   rFlat15;
	d15d8 = 	   rFlat15;
	d15d9 = 	   rFlat15;
	d15d10 = 	   rFlat15;
	d15d11 = 	   rFlat15;
	d15d12 = 	   rFlat15;
	d15d13 = 	   rFlat15;
	d15d14 = 	   rFlat15;
	d15d15 = 	   rFlat15;
	d15d16 = 	   rFlat15;
	d15d17 = 	   rFlat15;
	d15d18 = 	   rFlat15;
	d15d19 = 	   rFlat15;
	d15d20 = 	   rFlat15;
	d15d21 = 	   rFlat15;
	d15d22 = 	   rFlat15+c[27]*coeffM15*k15;
	d15d23 = 	   rFlat15;
	d15d24 = 	   rFlat15;
	d15d25 = 	   rFlat15;
	d15d26 = 	   rFlat15;
	d15d27 = 	   rFlat15+c[22]*coeffM15*k15;
	d15d28 = 	   rFlat15-coeffM15*k15*uK15;
	d15d29 = 	   rFlat15;
	d15d30 = 	   rFlat15;
	d15d31 = 	   rFlat15;
	d15d32 = 	   rFlat15;

	d16d6 = 	   c[19]*k16;
	d16d19 = 	   c[6]*k16;
	d16d20 = 	   -c[22]*k16*uK16;
	d16d22 = 	   -c[20]*k16*uK16;

	d17d6 = 	   c[32]*k17;
	d17d19 = 	   -c[28]*k17*uK17;
	d17d28 = 	   -c[19]*k17*uK17;
	d17d32 = 	   c[6]*k17;

	sigma18 =	   CFO18*dwFdM18+dCFOdM18*wF18;
	d18d1 = 	   rFlat18*sigma18*1.0;
	d18d2 = 	   rFlat18*sigma18*5.0E-1;
	d18d3 = 	   rFlat18*sigma18*1.0;
	d18d4 = 	   rFlat18*sigma18*1.0;
	d18d5 = 	   rFlat18*sigma18*1.0;
	d18d6 = 	   rFlat18*sigma18*1.0;
	d18d7 = 	   rFlat18*sigma18*1.0;
	d18d8 = 	   rFlat18*sigma18*1.0;
	d18d9 = 	   rFlat18*sigma18*1.0;
	d18d10 = 	   rFlat18*sigma18*1.0;
	d18d11 = 	   rFlat18*sigma18*1.0;
	d18d12 = 	   rFlat18*sigma18*1.0;
	d18d13 = 	   rFlat18*sigma18*2.0;
	d18d14 = 	   rFlat18*sigma18*1.0;
	d18d15 = 	   rFlat18*sigma18*1.2E1;
	d18d16 = 	   rFlat18*sigma18*1.0;
	d18d17 = 	   rFlat18*sigma18*1.0;
	d18d18 = 	   rFlat18*sigma18*1.0;
	d18d19 = 	   rFlat18*sigma18*1.5+c[22]*coeffFallOff18*k18;
	d18d20 = 	   rFlat18*sigma18*2.0-coeffFallOff18*k18*uK18;
	d18d21 = 	   rFlat18*sigma18*1.0;
	d18d22 = 	   rFlat18*sigma18*1.0+c[19]*coeffFallOff18*k18;
	d18d23 = 	   rFlat18*sigma18*1.0;
	d18d24 = 	   rFlat18*sigma18*1.0;
	d18d25 = 	   rFlat18*sigma18*1.0;
	d18d26 = 	   rFlat18*sigma18*1.0;
	d18d27 = 	   rFlat18*sigma18*1.0;
	d18d28 = 	   rFlat18*sigma18*1.0;
	d18d29 = 	   rFlat18*sigma18*1.0;
	d18d30 = 	   rFlat18*sigma18*1.0;
	d18d31 = 	   rFlat18*sigma18*1.0;
	d18d32 = 	   rFlat18*sigma18*1.0;

	d19d19 = 	   c[27]*k19;
	d19d20 = 	   -c[24]*k19*uK19;
	d19d24 = 	   -c[20]*k19*uK19;
	d19d27 = 	   c[19]*k19;

	d20d19 = 	   c[27]*k20;
	d20d20 = 	   -c[24]*k20*uK20;
	d20d24 = 	   -c[20]*k20*uK20;
	d20d27 = 	   c[19]*k20;

	d21d19 = 	   c[28]*k21;
	d21d20 = 	   -c[27]*k21*uK21;
	d21d27 = 	   -c[20]*k21*uK21;
	d21d28 = 	   c[19]*k21;

	d22d13 = 	   -c[20]*k22*uK22;
	d22d15 = 	   c[19]*k22;
	d22d19 = 	   c[15]*k22;
	d22d20 = 	   -c[13]*k22*uK22;

	d23d1 = 	   rFlat23*1.0;
	d23d2 = 	   rFlat23*1.0;
	d23d3 = 	   rFlat23*1.0;
	d23d4 = 	   rFlat23*1.0;
	d23d5 = 	   rFlat23*1.0;
	d23d6 = 	   rFlat23*1.0;
	d23d7 = 	   rFlat23*1.0;
	d23d8 = 	   rFlat23*1.0;
	d23d9 = 	   rFlat23*1.0;
	d23d10 = 	   rFlat23*1.0;
	d23d11 = 	   rFlat23*1.0;
	d23d12 = 	   rFlat23*1.0;
	d23d13 = 	   rFlat23*1.9;
	d23d14 = 	   rFlat23*1.0;
	d23d15 = 	   rFlat23*5.0;
	d23d16 = 	   rFlat23*1.0;
	d23d17 = 	   rFlat23*1.0;
	d23d18 = 	   rFlat23*1.0;
	d23d19 = 	   rFlat23*1.9-c[24]*coeffM23*k23*uK23;
	d23d20 = 	   rFlat23*3.0;
	d23d21 = 	   rFlat23*1.0;
	d23d22 = 	   rFlat23*1.0;
	d23d23 = 	   rFlat23*1.0;
	d23d24 = 	   rFlat23*1.0-c[19]*coeffM23*k23*uK23;
	d23d25 = 	   rFlat23*1.0;
	d23d26 = 	   rFlat23*1.0;
	d23d27 = 	   rFlat23*1.0;
	d23d28 = 	   rFlat23*1.0;
	d23d29 = 	   rFlat23*1.0;
	d23d30 = 	   rFlat23*1.0;
	d23d31 = 	   rFlat23*1.0;
	d23d32 = 	   rFlat23*1.0+coeffM23*k23;

	d24d20 = 	   -c[24]*k24*uK24;
	d24d22 = 	   c[32]*k24;
	d24d24 = 	   -c[20]*k24*uK24;
	d24d32 = 	   c[22]*k24;

	d25d13 = 	   -c[19]*k25*uK25;
	d25d19 = 	   -c[13]*k25*uK25;
	d25d24 = 	   c[32]*k25;
	d25d32 = 	   c[24]*k25;

	d26d15 = 	   -c[19]*k26*uK26;
	d26d19 = 	   -c[15]*k26*uK26;
	d26d27 = 	   c[32]*k26;
	d26d32 = 	   c[27]*k26;

	d27d16 = 	   -c[19]*k27*uK27;
	d27d19 = 	   -c[16]*k27*uK27;
	d27d28 = 	   c[32]*k27;
	d27d32 = 	   c[28]*k27;

	d28d28 = 	   c[32]*k28;
	d28d32 = 	   c[28]*k28;

	d29d13 = 	   -c[27]*k29*uK29;
	d29d15 = 	   c[24]*k29;
	d29d24 = 	   c[15]*k29;
	d29d27 = 	   -c[13]*k29*uK29;

	d30d15 = 	   -c[27]*k30*uK30;
	d30d16 = 	   c[24]*k30;
	d30d24 = 	   c[16]*k30;
	d30d27 = 	   -c[15]*k30*uK30;

	d31d13 = 	   -c[28]*k31*uK31;
	d31d16 = 	   c[24]*k31;
	d31d24 = 	   c[16]*k31;
	d31d28 = 	   -c[13]*k31*uK31;

	d32d1 = 	   rFlat32;
	d32d2 = 	   rFlat32;
	d32d3 = 	   rFlat32;
	d32d4 = 	   rFlat32;
	d32d5 = 	   rFlat32;
	d32d6 = 	   rFlat32;
	d32d7 = 	   rFlat32;
	d32d8 = 	   rFlat32;
	d32d9 = 	   rFlat32;
	d32d10 = 	   rFlat32;
	d32d11 = 	   rFlat32;
	d32d12 = 	   rFlat32;
	d32d13 = 	   rFlat32;
	d32d14 = 	   rFlat32;
	d32d15 = 	   rFlat32;
	d32d16 = 	   rFlat32;
	d32d17 = 	   rFlat32+coeffM32*k32;
	d32d18 = 	   rFlat32;
	d32d19 = 	   rFlat32;
	d32d20 = 	   rFlat32;
	d32d21 = 	   rFlat32;
	d32d22 = 	   rFlat32;
	d32d23 = 	   rFlat32;
	d32d24 = 	   rFlat32-c[29]*coeffM32*k32*uK32;
	d32d25 = 	   rFlat32;
	d32d26 = 	   rFlat32;
	d32d27 = 	   rFlat32;
	d32d28 = 	   rFlat32;
	d32d29 = 	   rFlat32-c[24]*coeffM32*k32*uK32;
	d32d30 = 	   rFlat32;
	d32d31 = 	   rFlat32;
	d32d32 = 	   rFlat32;

	d33d13 = 	   -c[25]*k33*uK33;
	d33d24 = 	   c[29]*k33;
	d33d25 = 	   -c[13]*k33*uK33;
	d33d29 = 	   c[24]*k33;

	d34d8 = 	   -c[24]*k34*uK34;
	d34d22 = 	   c[29]*k34;
	d34d24 = 	   -c[8]*k34*uK34;
	d34d29 = 	   c[22]*k34;

	d35d22 = 	   c[29]*k35;
	d35d25 = 	   -c[27]*k35*uK35;
	d35d27 = 	   -c[25]*k35*uK35;
	d35d29 = 	   c[22]*k35;

	d36d4 = 	   -c[13]*k36*uK36;
	d36d13 = 	   -c[4]*k36*uK36;
	d36d22 = 	   c[29]*k36;
	d36d29 = 	   c[22]*k36;

	d37d6 = 	   c[29]*k37;
	d37d8 = 	   -c[27]*k37*uK37;
	d37d27 = 	   -c[8]*k37*uK37;
	d37d29 = 	   c[6]*k37;

	d38d15 = 	   -c[25]*k38*uK38;
	d38d25 = 	   -c[15]*k38*uK38;
	d38d27 = 	   c[29]*k38;
	d38d29 = 	   c[27]*k38;

	d39d27 = 	   -c[30]*k39*uK39;
	d39d28 = 	   c[29]*k39;
	d39d29 = 	   c[28]*k39;
	d39d30 = 	   -c[27]*k39*uK39;

	d40d6 = 	   -c[29]*k40*uK40;
	d40d22 = 	   c[30]*k40;
	d40d29 = 	   -c[6]*k40*uK40;
	d40d30 = 	   c[22]*k40;

	d41d3 = 	   -(c[24]*c[24])*k41*uK41;
	d41d21 = 	   c[29]*k41;
	d41d24 = 	   c[24]*c[3]*k41*uK41*-2.0;
	d41d29 = 	   c[21]*k41;

	d42d14 = 	   -c[24]*k42*uK42;
	d42d24 = 	   -c[14]*k42*uK42;
	d42d25 = 	   c[29]*k42;
	d42d29 = 	   c[25]*k42;

	d43d13 = 	   -c[14]*k43*uK43;
	d43d14 = 	   -c[13]*k43*uK43;
	d43d29 = 	   c[29]*k43*2.0;

	d44d17 = 	   -c[25]*k44*uK44;
	d44d25 = 	   -c[17]*k44*uK44;
	d44d29 = 	   c[29]*k44*2.0;

	sigma45 =	   CFO45*dwFdM45+dCFOdM45*wF45;
	d45d1 = 	   rFlat45*sigma45*1.0;
	d45d2 = 	   rFlat45*sigma45*1.0;
	d45d3 = 	   rFlat45*sigma45*2.5;
	d45d4 = 	   rFlat45*sigma45*1.0;
	d45d5 = 	   rFlat45*sigma45*1.0;
	d45d6 = 	   rFlat45*sigma45*1.0;
	d45d7 = 	   rFlat45*sigma45*1.0;
	d45d8 = 	   rFlat45*sigma45*1.0;
	d45d9 = 	   rFlat45*sigma45*1.0;
	d45d10 = 	   rFlat45*sigma45*1.0;
	d45d11 = 	   rFlat45*sigma45*1.0;
	d45d12 = 	   rFlat45*sigma45*1.0;
	d45d13 = 	   rFlat45*sigma45*1.0;
	d45d14 = 	   rFlat45*sigma45*1.0;
	d45d15 = 	   rFlat45*sigma45*5.0;
	d45d16 = 	   rFlat45*sigma45*1.0;
	d45d17 = 	   rFlat45*sigma45*1.0E1;
	d45d18 = 	   rFlat45*sigma45*1.0-coeffFallOff45*k45*uK45;
	d45d19 = 	   rFlat45*sigma45*1.0;
	d45d20 = 	   rFlat45*sigma45*1.0;
	d45d21 = 	   rFlat45*sigma45*1.0;
	d45d22 = 	   rFlat45*sigma45*1.0;
	d45d23 = 	   rFlat45*sigma45*1.0;
	d45d24 = 	   rFlat45*sigma45*1.0;
	d45d25 = 	   rFlat45*sigma45*1.0;
	d45d26 = 	   rFlat45*sigma45*1.0;
	d45d27 = 	   rFlat45*sigma45*1.0;
	d45d28 = 	   rFlat45*sigma45*1.0;
	d45d29 = 	   rFlat45*sigma45*1.0+c[29]*coeffFallOff45*k45*2.0;
	d45d30 = 	   rFlat45*sigma45*1.0;
	d45d31 = 	   rFlat45*sigma45*1.0;
	d45d32 = 	   rFlat45*sigma45*1.0;

	d46d4 = 	   c[29]*k46;
	d46d26 = 	   -c[27]*k46*uK46;
	d46d27 = 	   -c[26]*k46*uK46;
	d46d29 = 	   c[4]*k46;

	d47d3 = 	   -c[15]*k47*uK47;
	d47d4 = 	   c[29]*k47;
	d47d15 = 	   -c[3]*k47*uK47;
	d47d29 = 	   c[4]*k47;

	d48d5 = 	   -c[15]*k48*uK48;
	d48d7 = 	   c[29]*k48;
	d48d15 = 	   -c[5]*k48*uK48;
	d48d29 = 	   c[7]*k48;

	d49d4 = 	   -c[30]*k49*uK49;
	d49d7 = 	   c[29]*k49;
	d49d29 = 	   c[7]*k49;
	d49d30 = 	   -c[4]*k49*uK49;

	d50d13 = 	   -c[21]*k50*uK50;
	d50d21 = 	   -c[13]*k50*uK50;
	d50d24 = 	   c[25]*k50;
	d50d25 = 	   c[24]*k50;

	d51d4 = 	   -c[24]*k51*uK51;
	d51d22 = 	   c[25]*k51;
	d51d24 = 	   -c[4]*k51*uK51;
	d51d25 = 	   c[22]*k51;

	d52d8 = 	   -c[24]*k52*uK52;
	d52d24 = 	   -c[8]*k52*uK52;
	d52d25 = 	   c[27]*k52;
	d52d27 = 	   c[25]*k52;

	d53d15 = 	   -c[21]*k53*uK53;
	d53d21 = 	   -c[15]*k53*uK53;
	d53d25 = 	   c[27]*k53;
	d53d27 = 	   c[25]*k53;

	d54d4 = 	   -c[13]*k54*uK54;
	d54d13 = 	   -c[4]*k54*uK54;
	d54d25 = 	   c[27]*k54;
	d54d27 = 	   c[25]*k54;

	d55d6 = 	   c[25]*k55;
	d55d8 = 	   -c[22]*k55*uK55;
	d55d22 = 	   -c[8]*k55*uK55;
	d55d25 = 	   c[6]*k55;

	d56d4 = 	   -c[27]*k56*uK56;
	d56d6 = 	   c[25]*k56;
	d56d25 = 	   c[6]*k56;
	d56d27 = 	   -c[4]*k56*uK56;

	d57d3 = 	   -(c[24]*c[24])*k57*uK57;
	d57d24 = 	   c[24]*c[3]*k57*uK57*-2.0;
	d57d25 = 	   c[25]*k57*2.0;

	d58d3 = 	   -c[24]*k58*uK58;
	d58d21 = 	   c[25]*k58;
	d58d24 = 	   -c[3]*k58*uK58;
	d58d25 = 	   c[21]*k58;

	d59d3 = 	   -c[27]*k59*uK59;
	d59d4 = 	   c[25]*k59;
	d59d25 = 	   c[4]*k59;
	d59d27 = 	   -c[3]*k59*uK59;

	d60d5 = 	   -c[27]*k60*uK60;
	d60d7 = 	   c[25]*k60;
	d60d25 = 	   c[7]*k60;
	d60d27 = 	   -c[5]*k60*uK60;

	d61d4 = 	   -c[24]*k61*uK61;
	d61d21 = 	   c[27]*k61;
	d61d24 = 	   -c[4]*k61*uK61;
	d61d27 = 	   c[21]*k61;

	d62d4 = 	   -c[22]*k62*uK62;
	d62d6 = 	   c[21]*k62;
	d62d21 = 	   c[6]*k62;
	d62d22 = 	   -c[4]*k62*uK62;

	d63d3 = 	   -c[22]*k63*uK63;
	d63d4 = 	   c[21]*k63;
	d63d21 = 	   c[4]*k63;
	d63d22 = 	   -c[3]*k63*uK63;

	d64d13 = 	   -c[31]*k64*uK64;
	d64d18 = 	   c[24]*k64;
	d64d24 = 	   c[18]*k64;
	d64d31 = 	   -c[13]*k64*uK64;

	d65d14 = 	   -c[15]*k65*uK65;
	d65d15 = 	   -c[14]*k65*uK65;
	d65d18 = 	   c[22]*k65;
	d65d22 = 	   c[18]*k65;

	d66d15 = 	   -c[31]*k66*uK66;
	d66d18 = 	   c[27]*k66;
	d66d27 = 	   c[18]*k66;
	d66d31 = 	   -c[15]*k66*uK66;

	d67d17 = 	   -c[31]*k67*uK67;
	d67d18 = 	   c[29]*k67;
	d67d29 = 	   c[18]*k67;
	d67d31 = 	   -c[17]*k67*uK67;

	d68d1 = 	   rFlat68;
	d68d2 = 	   rFlat68;
	d68d3 = 	   rFlat68;
	d68d4 = 	   rFlat68;
	d68d5 = 	   rFlat68;
	d68d6 = 	   rFlat68;
	d68d7 = 	   rFlat68;
	d68d8 = 	   rFlat68;
	d68d9 = 	   rFlat68;
	d68d10 = 	   rFlat68;
	d68d11 = 	   rFlat68;
	d68d12 = 	   rFlat68;
	d68d13 = 	   rFlat68;
	d68d14 = 	   rFlat68-c[24]*coeffM68*k68*uK68;
	d68d15 = 	   rFlat68;
	d68d16 = 	   rFlat68;
	d68d17 = 	   rFlat68;
	d68d18 = 	   rFlat68;
	d68d19 = 	   rFlat68;
	d68d20 = 	   rFlat68;
	d68d21 = 	   rFlat68;
	d68d22 = 	   rFlat68;
	d68d23 = 	   rFlat68;
	d68d24 = 	   rFlat68-c[14]*coeffM68*k68*uK68;
	d68d25 = 	   rFlat68;
	d68d26 = 	   rFlat68;
	d68d27 = 	   rFlat68;
	d68d28 = 	   rFlat68;
	d68d29 = 	   rFlat68;
	d68d30 = 	   rFlat68;
	d68d31 = 	   rFlat68+coeffM68*k68;
	d68d32 = 	   rFlat68;

	d69d24 = 	   c[31]*k69;
	d69d29 = 	   c[29]*k69*uK69*-2.0;
	d69d31 = 	   c[24]*k69;

	d70d14 = 	   -c[27]*k70*uK70;
	d70d22 = 	   c[31]*k70;
	d70d27 = 	   -c[14]*k70*uK70;
	d70d31 = 	   c[22]*k70;

	d71d8 = 	   -c[29]*k71*uK71;
	d71d22 = 	   c[31]*k71;
	d71d29 = 	   -c[8]*k71*uK71;
	d71d31 = 	   c[22]*k71;

	d72d14 = 	   -c[15]*k72*uK72;
	d72d15 = 	   -c[14]*k72*uK72;
	d72d27 = 	   c[31]*k72;
	d72d31 = 	   c[27]*k72;

	d73d8 = 	   -c[17]*k73*uK73;
	d73d17 = 	   -c[8]*k73*uK73;
	d73d27 = 	   c[31]*k73;
	d73d31 = 	   c[27]*k73;

	d74d14 = 	   -c[29]*k74*uK74;
	d74d25 = 	   c[31]*k74;
	d74d29 = 	   -c[14]*k74*uK74;
	d74d31 = 	   c[25]*k74;

	d75d1 = 	   rFlat75*1.0;
	d75d2 = 	   rFlat75*1.0;
	d75d3 = 	   rFlat75*2.0;
	d75d4 = 	   rFlat75*1.0;
	d75d5 = 	   rFlat75*1.0;
	d75d6 = 	   rFlat75*2.0;
	d75d7 = 	   rFlat75*1.0;
	d75d8 = 	   rFlat75*1.0;
	d75d9 = 	   rFlat75*1.0;
	d75d10 = 	   rFlat75*1.0;
	d75d11 = 	   rFlat75*1.0;
	d75d12 = 	   rFlat75*1.0;
	d75d13 = 	   rFlat75*2.0;
	d75d14 = 	   rFlat75*1.0+coeffM75*k75;
	d75d15 = 	   rFlat75*1.5E1;
	d75d16 = 	   rFlat75*1.0;
	d75d17 = 	   rFlat75*1.0;
	d75d18 = 	   rFlat75*1.0;
	d75d19 = 	   rFlat75*1.0;
	d75d20 = 	   rFlat75*1.0;
	d75d21 = 	   rFlat75*1.0;
	d75d22 = 	   rFlat75*1.0;
	d75d23 = 	   rFlat75*1.0;
	d75d24 = 	   rFlat75*1.0-c[26]*coeffM75*k75*uK75;
	d75d25 = 	   rFlat75*1.0;
	d75d26 = 	   rFlat75*1.0-c[24]*coeffM75*k75*uK75;
	d75d27 = 	   rFlat75*1.0;
	d75d28 = 	   rFlat75*1.0;
	d75d29 = 	   rFlat75*1.0;
	d75d30 = 	   rFlat75*1.0;
	d75d31 = 	   rFlat75*1.0;
	d75d32 = 	   rFlat75*1.0;

	d76d13 = 	   -c[26]*k76*uK76;
	d76d14 = 	   c[24]*k76;
	d76d24 = 	   c[14]*k76;
	d76d26 = 	   -c[13]*k76*uK76;

	d77d4 = 	   -c[29]*k77*uK77;
	d77d14 = 	   c[22]*k77;
	d77d22 = 	   c[14]*k77;
	d77d29 = 	   -c[4]*k77*uK77;

	d78d14 = 	   c[22]*k78;
	d78d22 = 	   c[14]*k78;
	d78d26 = 	   -c[27]*k78*uK78;
	d78d27 = 	   -c[26]*k78*uK78;

	d79d14 = 	   c[27]*k79;
	d79d15 = 	   -c[26]*k79*uK79;
	d79d26 = 	   -c[15]*k79*uK79;
	d79d27 = 	   c[14]*k79;

	d80d4 = 	   c[14]*k80;
	d80d5 = 	   -c[29]*k80*uK80;
	d80d14 = 	   c[4]*k80;
	d80d29 = 	   -c[5]*k80*uK80;

	d81d14 = 	   c[25]*k81;
	d81d25 = 	   c[14]*k81;
	d81d26 = 	   -c[29]*k81*uK81;
	d81d29 = 	   -c[26]*k81*uK81;

	d82d14 = 	   c[29]*k82;
	d82d17 = 	   -c[26]*k82*uK82;
	d82d26 = 	   -c[17]*k82*uK82;
	d82d29 = 	   c[14]*k82;

	d83d3 = 	   -c[24]*k83*uK83;
	d83d24 = 	   -c[3]*k83*uK83;
	d83d26 = 	   k83;

	d84d3 = 	   -c[13]*k84*uK84;
	d84d13 = 	   -c[3]*k84*uK84;
	d84d24 = 	   c[26]*k84;
	d84d26 = 	   c[24]*k84;

	d85d5 = 	   -c[24]*k85*uK85;
	d85d22 = 	   c[26]*k85;
	d85d24 = 	   -c[5]*k85*uK85;
	d85d26 = 	   c[22]*k85;

	d86d4 = 	   -c[25]*k86*uK86;
	d86d22 = 	   c[26]*k86;
	d86d25 = 	   -c[4]*k86*uK86;
	d86d26 = 	   c[22]*k86;

	d87d3 = 	   -c[15]*k87*uK87;
	d87d15 = 	   -c[3]*k87*uK87;
	d87d26 = 	   c[27]*k87;
	d87d27 = 	   c[26]*k87;

	d88d3 = 	   -c[28]*k88*uK88;
	d88d6 = 	   c[26]*k88;
	d88d26 = 	   c[6]*k88;
	d88d28 = 	   -c[3]*k88*uK88;

	d89d3 = 	   -c[24]*c[6]*k89*uK89;
	d89d6 = 	   k89*(c[26]-c[24]*c[3]*uK89);
	d89d24 = 	   -c[3]*c[6]*k89*uK89;
	d89d26 = 	   c[6]*k89;

	d90d3 = 	   -c[29]*k90*uK90;
	d90d25 = 	   c[26]*k90;
	d90d26 = 	   c[25]*k90;
	d90d29 = 	   -c[3]*k90*uK90;

	d91d3 = 	   -c[17]*k91*uK91;
	d91d17 = 	   -c[3]*k91*uK91;
	d91d26 = 	   c[29]*k91;
	d91d29 = 	   c[26]*k91;

	d92d3 = 	   -c[8]*k92*uK92;
	d92d4 = 	   c[26]*k92;
	d92d8 = 	   -c[3]*k92*uK92;
	d92d26 = 	   c[4]*k92;

	d93d1 = 	   rFlat93;
	d93d2 = 	   rFlat93;
	d93d3 = 	   rFlat93;
	d93d4 = 	   rFlat93;
	d93d5 = 	   rFlat93-c[24]*coeffM93*k93*uK93;
	d93d6 = 	   rFlat93;
	d93d7 = 	   rFlat93;
	d93d8 = 	   rFlat93;
	d93d9 = 	   rFlat93+coeffM93*k93;
	d93d10 = 	   rFlat93;
	d93d11 = 	   rFlat93;
	d93d12 = 	   rFlat93;
	d93d13 = 	   rFlat93;
	d93d14 = 	   rFlat93;
	d93d15 = 	   rFlat93;
	d93d16 = 	   rFlat93;
	d93d17 = 	   rFlat93;
	d93d18 = 	   rFlat93;
	d93d19 = 	   rFlat93;
	d93d20 = 	   rFlat93;
	d93d21 = 	   rFlat93;
	d93d22 = 	   rFlat93;
	d93d23 = 	   rFlat93;
	d93d24 = 	   rFlat93-c[5]*coeffM93*k93*uK93;
	d93d25 = 	   rFlat93;
	d93d26 = 	   rFlat93;
	d93d27 = 	   rFlat93;
	d93d28 = 	   rFlat93;
	d93d29 = 	   rFlat93;
	d93d30 = 	   rFlat93;
	d93d31 = 	   rFlat93;
	d93d32 = 	   rFlat93;

	d94d1 = 	   rFlat94;
	d94d2 = 	   rFlat94;
	d94d3 = 	   rFlat94-c[27]*coeffM94*k94*uK94;
	d94d4 = 	   rFlat94;
	d94d5 = 	   rFlat94;
	d94d6 = 	   rFlat94;
	d94d7 = 	   rFlat94;
	d94d8 = 	   rFlat94;
	d94d9 = 	   rFlat94+coeffM94*k94;
	d94d10 = 	   rFlat94;
	d94d11 = 	   rFlat94;
	d94d12 = 	   rFlat94;
	d94d13 = 	   rFlat94;
	d94d14 = 	   rFlat94;
	d94d15 = 	   rFlat94;
	d94d16 = 	   rFlat94;
	d94d17 = 	   rFlat94;
	d94d18 = 	   rFlat94;
	d94d19 = 	   rFlat94;
	d94d20 = 	   rFlat94;
	d94d21 = 	   rFlat94;
	d94d22 = 	   rFlat94;
	d94d23 = 	   rFlat94;
	d94d24 = 	   rFlat94;
	d94d25 = 	   rFlat94;
	d94d26 = 	   rFlat94;
	d94d27 = 	   rFlat94-c[3]*coeffM94*k94*uK94;
	d94d28 = 	   rFlat94;
	d94d29 = 	   rFlat94;
	d94d30 = 	   rFlat94;
	d94d31 = 	   rFlat94;
	d94d32 = 	   rFlat94;

	d95d5 = 	   -c[13]*k95*uK95;
	d95d9 = 	   c[24]*k95;
	d95d13 = 	   -c[5]*k95*uK95;
	d95d24 = 	   c[9]*k95;

	d96d9 = 	   c[24]*k96;
	d96d24 = 	   c[9]*k96;
	d96d26 = 	   -c[27]*k96*uK96;
	d96d27 = 	   -c[26]*k96*uK96;

	d97d5 = 	   -c[27]*k97*uK97;
	d97d9 = 	   c[22]*k97;
	d97d22 = 	   c[9]*k97;
	d97d27 = 	   -c[5]*k97*uK97;

	d98d6 = 	   -c[26]*k98*uK98;
	d98d9 = 	   c[22]*k98;
	d98d22 = 	   c[9]*k98;
	d98d26 = 	   -c[6]*k98*uK98;

	d99d5 = 	   -c[15]*k99*uK99;
	d99d9 = 	   c[27]*k99;
	d99d15 = 	   -c[5]*k99*uK99;
	d99d27 = 	   c[9]*k99;

	d100d9 = 	   c[27]*k100;
	d100d26 = 	   -c[28]*k100*uK100;
	d100d27 = 	   c[9]*k100;
	d100d28 = 	   -c[26]*k100*uK100;

	d101d4 = 	   c[9]*k101;
	d101d5 = 	   -c[8]*k101*uK101;
	d101d8 = 	   -c[5]*k101*uK101;
	d101d9 = 	   c[4]*k101;

	d102d4 = 	   c[9]*k102;
	d102d7 = 	   -c[26]*k102*uK102;
	d102d9 = 	   c[4]*k102;
	d102d26 = 	   -c[7]*k102*uK102;

	d103d7 = 	   c[9]*k103;
	d103d9 = 	   c[7]*k103;
	d103d23 = 	   -c[26]*k103*uK103;
	d103d26 = 	   -c[23]*k103*uK103;

	d104d5 = 	   -c[10]*k104*uK104;
	d104d7 = 	   c[9]*k104;
	d104d9 = 	   c[7]*k104;
	d104d10 = 	   -c[5]*k104*uK104;

	d105d4 = 	   c[28]*k105;
	d105d7 = 	   -c[27]*k105*uK105;
	d105d27 = 	   -c[7]*k105*uK105;
	d105d28 = 	   c[4]*k105;

	d106d1 = 	   rFlat106*1.0;
	d106d2 = 	   rFlat106*1.0;
	d106d3 = 	   rFlat106*1.7;
	d106d4 = 	   rFlat106*1.0+c[22]*coeffM106*k106;
	d106d5 = 	   rFlat106*1.0;
	d106d6 = 	   rFlat106*1.5;
	d106d7 = 	   rFlat106*1.0-coeffM106*k106*uK106;
	d106d8 = 	   rFlat106*1.0;
	d106d9 = 	   rFlat106*1.0;
	d106d10 = 	   rFlat106*1.0;
	d106d11 = 	   rFlat106*1.0;
	d106d12 = 	   rFlat106*1.0;
	d106d13 = 	   rFlat106*1.0;
	d106d14 = 	   rFlat106*1.0;
	d106d15 = 	   rFlat106*1.0E1;
	d106d16 = 	   rFlat106*1.0;
	d106d17 = 	   rFlat106*1.0;
	d106d18 = 	   rFlat106*1.0;
	d106d19 = 	   rFlat106*1.0;
	d106d20 = 	   rFlat106*1.0;
	d106d21 = 	   rFlat106*1.0;
	d106d22 = 	   rFlat106*1.0+c[4]*coeffM106*k106;
	d106d23 = 	   rFlat106*1.0;
	d106d24 = 	   rFlat106*1.0;
	d106d25 = 	   rFlat106*1.0;
	d106d26 = 	   rFlat106*1.0;
	d106d27 = 	   rFlat106*1.0;
	d106d28 = 	   rFlat106*1.0;
	d106d29 = 	   rFlat106*1.0;
	d106d30 = 	   rFlat106*1.0;
	d106d31 = 	   rFlat106*1.0;
	d106d32 = 	   rFlat106*1.0;

	sigma107 =	   CFO107*dwFdM107+dCFOdM107*wF107;
	d107d1 = 	   rFlat107*sigma107;
	d107d2 = 	   rFlat107*sigma107;
	d107d3 = 	   rFlat107*sigma107;
	d107d4 = 	   rFlat107*sigma107+c[27]*coeffFallOff107*k107;
	d107d5 = 	   rFlat107*sigma107;
	d107d6 = 	   rFlat107*sigma107;
	d107d7 = 	   rFlat107*sigma107;
	d107d8 = 	   rFlat107*sigma107;
	d107d9 = 	   rFlat107*sigma107;
	d107d10 = 	   rFlat107*sigma107-coeffFallOff107*k107*uK107;
	d107d11 = 	   rFlat107*sigma107;
	d107d12 = 	   rFlat107*sigma107;
	d107d13 = 	   rFlat107*sigma107;
	d107d14 = 	   rFlat107*sigma107;
	d107d15 = 	   rFlat107*sigma107;
	d107d16 = 	   rFlat107*sigma107;
	d107d17 = 	   rFlat107*sigma107;
	d107d18 = 	   rFlat107*sigma107;
	d107d19 = 	   rFlat107*sigma107;
	d107d20 = 	   rFlat107*sigma107;
	d107d21 = 	   rFlat107*sigma107;
	d107d22 = 	   rFlat107*sigma107;
	d107d23 = 	   rFlat107*sigma107;
	d107d24 = 	   rFlat107*sigma107;
	d107d25 = 	   rFlat107*sigma107;
	d107d26 = 	   rFlat107*sigma107;
	d107d27 = 	   rFlat107*sigma107+c[4]*coeffFallOff107*k107;
	d107d28 = 	   rFlat107*sigma107;
	d107d29 = 	   rFlat107*sigma107;
	d107d30 = 	   rFlat107*sigma107;
	d107d31 = 	   rFlat107*sigma107;
	d107d32 = 	   rFlat107*sigma107;

	d108d4 = 	   c[32]*k108;
	d108d8 = 	   -c[19]*k108*uK108;
	d108d19 = 	   -c[8]*k108*uK108;
	d108d32 = 	   c[4]*k108;

	d109d1 = 	   rFlat109*1.0;
	d109d2 = 	   rFlat109*1.0;
	d109d3 = 	   rFlat109*1.0;
	d109d4 = 	   rFlat109*1.0+c[24]*coeffM109*k109;
	d109d5 = 	   rFlat109*1.0;
	d109d6 = 	   rFlat109*1.0;
	d109d7 = 	   rFlat109*1.0;
	d109d8 = 	   rFlat109*1.0-coeffM109*k109*uK109;
	d109d9 = 	   rFlat109*1.0;
	d109d10 = 	   rFlat109*1.0;
	d109d11 = 	   rFlat109*1.0;
	d109d12 = 	   rFlat109*1.0;
	d109d13 = 	   rFlat109*1.25;
	d109d14 = 	   rFlat109*1.0;
	d109d15 = 	   rFlat109*4.1;
	d109d16 = 	   rFlat109*1.0;
	d109d17 = 	   rFlat109*1.0;
	d109d18 = 	   rFlat109*1.0;
	d109d19 = 	   rFlat109*1.0;
	d109d20 = 	   rFlat109*1.0;
	d109d21 = 	   rFlat109*1.0;
	d109d22 = 	   rFlat109*1.0;
	d109d23 = 	   rFlat109*1.0;
	d109d24 = 	   rFlat109*1.0+c[4]*coeffM109*k109;
	d109d25 = 	   rFlat109*1.0;
	d109d26 = 	   rFlat109*1.0;
	d109d27 = 	   rFlat109*1.0;
	d109d28 = 	   rFlat109*1.0;
	d109d29 = 	   rFlat109*1.0;
	d109d30 = 	   rFlat109*1.0;
	d109d31 = 	   rFlat109*1.0;
	d109d32 = 	   rFlat109*1.0;

	d110d4 = 	   -c[13]*k110*uK110;
	d110d8 = 	   c[24]*k110;
	d110d13 = 	   -c[4]*k110*uK110;
	d110d24 = 	   c[8]*k110;

	d111d4 = 	   -c[27]*k111*uK111;
	d111d8 = 	   c[22]*k111;
	d111d22 = 	   c[8]*k111;
	d111d27 = 	   -c[4]*k111*uK111;

	d112d4 = 	   -c[15]*k112*uK112;
	d112d8 = 	   c[27]*k112;
	d112d15 = 	   -c[4]*k112*uK112;
	d112d27 = 	   c[8]*k112;

	d113d4 = 	   -c[28]*k113*uK113;
	d113d6 = 	   c[8]*k113;
	d113d8 = 	   c[6]*k113;
	d113d28 = 	   -c[4]*k113*uK113;

	d114d4 = 	   -c[17]*k114*uK114;
	d114d8 = 	   c[29]*k114;
	d114d17 = 	   -c[4]*k114*uK114;
	d114d29 = 	   c[8]*k114;

	d115d4 = 	   c[8]*k115;
	d115d5 = 	   -c[27]*k115*uK115;
	d115d8 = 	   c[4]*k115;
	d115d27 = 	   -c[5]*k115*uK115;

	d116d4 = 	   -c[10]*k116*uK116;
	d116d7 = 	   c[8]*k116;
	d116d8 = 	   c[7]*k116;
	d116d10 = 	   -c[4]*k116*uK116;

	d117d5 = 	   -c[15]*k117*uK117;
	d117d8 = 	   c[8]*k117*2.0;
	d117d15 = 	   -c[5]*k117*uK117;

	d118d7 = 	   c[13]*k118;
	d118d10 = 	   -c[24]*k118*uK118;
	d118d13 = 	   c[7]*k118;
	d118d24 = 	   -c[10]*k118*uK118;

	d119d7 = 	   -c[27]*k119*uK119;
	d119d10 = 	   c[22]*k119;
	d119d22 = 	   c[10]*k119;
	d119d27 = 	   -c[7]*k119*uK119;

	d120d7 = 	   -c[15]*k120*uK120;
	d120d10 = 	   c[27]*k120;
	d120d15 = 	   -c[7]*k120*uK120;
	d120d27 = 	   c[10]*k120;

	d121d7 = 	   -c[29]*k121*uK121;
	d121d10 = 	   c[25]*k121;
	d121d25 = 	   c[10]*k121;
	d121d29 = 	   -c[7]*k121*uK121;

	d122d7 = 	   -c[17]*k122*uK122;
	d122d10 = 	   c[29]*k122;
	d122d17 = 	   -c[7]*k122*uK122;
	d122d29 = 	   c[10]*k122;

	d123d4 = 	   -c[15]*c[7]*k123*uK123;
	d123d7 = 	   -c[15]*c[4]*k123*uK123;
	d123d10 = 	   c[10]*k123*2.0;
	d123d15 = 	   -c[4]*c[7]*k123*uK123;

	d124d1 = 	   rFlat124;
	d124d2 = 	   rFlat124;
	d124d3 = 	   rFlat124;
	d124d4 = 	   rFlat124;
	d124d5 = 	   rFlat124;
	d124d6 = 	   rFlat124;
	d124d7 = 	   rFlat124;
	d124d8 = 	   rFlat124-c[24]*coeffM124*k124*uK124;
	d124d9 = 	   rFlat124;
	d124d10 = 	   rFlat124;
	d124d11 = 	   rFlat124;
	d124d12 = 	   rFlat124;
	d124d13 = 	   rFlat124;
	d124d14 = 	   rFlat124;
	d124d15 = 	   rFlat124;
	d124d16 = 	   rFlat124;
	d124d17 = 	   rFlat124;
	d124d18 = 	   rFlat124;
	d124d19 = 	   rFlat124;
	d124d20 = 	   rFlat124;
	d124d21 = 	   rFlat124;
	d124d22 = 	   rFlat124;
	d124d23 = 	   rFlat124;
	d124d24 = 	   rFlat124-c[8]*coeffM124*k124*uK124;
	d124d25 = 	   rFlat124;
	d124d26 = 	   rFlat124;
	d124d27 = 	   rFlat124;
	d124d28 = 	   rFlat124;
	d124d29 = 	   rFlat124;
	d124d30 = 	   rFlat124+coeffM124*k124;
	d124d31 = 	   rFlat124;
	d124d32 = 	   rFlat124;

	d125d8 = 	   -c[13]*k125*uK125;
	d125d13 = 	   -c[8]*k125*uK125;
	d125d24 = 	   c[30]*k125;
	d125d30 = 	   c[24]*k125;

	d126d24 = 	   c[30]*k126;
	d126d27 = 	   -c[29]*k126*uK126;
	d126d29 = 	   -c[27]*k126*uK126;
	d126d30 = 	   c[24]*k126;

	d127d8 = 	   -c[27]*k127*uK127;
	d127d22 = 	   c[30]*k127;
	d127d27 = 	   -c[8]*k127*uK127;
	d127d30 = 	   c[22]*k127;

	d128d8 = 	   -c[15]*k128*uK128;
	d128d15 = 	   -c[8]*k128*uK128;
	d128d27 = 	   c[30]*k128;
	d128d30 = 	   c[27]*k128;

	d129d4 = 	   c[30]*k129;
	d129d8 = 	   c[8]*k129*uK129*-2.0;
	d129d30 = 	   c[4]*k129;

	d130d8 = 	   -c[17]*k130*uK130;
	d130d17 = 	   -c[8]*k130*uK130;
	d130d29 = 	   c[30]*k130;
	d130d30 = 	   c[29]*k130;

	d131d7 = 	   c[30]*k131;
	d131d8 = 	   -c[10]*k131*uK131;
	d131d10 = 	   -c[8]*k131*uK131;
	d131d30 = 	   c[7]*k131;

	d132d7 = 	   -c[27]*k132*uK132;
	d132d23 = 	   c[24]*k132;
	d132d24 = 	   c[23]*k132;
	d132d27 = 	   -c[7]*k132*uK132;

	d133d6 = 	   -c[7]*k133*uK133;
	d133d7 = 	   -c[6]*k133*uK133;
	d133d22 = 	   c[23]*k133;
	d133d23 = 	   c[22]*k133;

	d134d7 = 	   -c[28]*k134*uK134;
	d134d23 = 	   c[27]*k134;
	d134d27 = 	   c[23]*k134;
	d134d28 = 	   -c[7]*k134*uK134;

	d135d6 = 	   -c[27]*c[7]*k135*uK135;
	d135d7 = 	   -c[27]*c[6]*k135*uK135;
	d135d23 = 	   c[28]*k135;
	d135d27 = 	   -c[6]*c[7]*k135*uK135;
	d135d28 = 	   c[23]*k135;

	d136d4 = 	   -c[6]*c[7]*k136*uK136;
	d136d6 = 	   -c[4]*c[7]*k136*uK136;
	d136d7 = 	   k136*(c[23]-c[4]*c[6]*uK136);
	d136d23 = 	   c[7]*k136;

	d137d4 = 	   -c[6]*k137*uK137;
	d137d6 = 	   -c[4]*k137*uK137;
	d137d23 = 	   k137;

	d138d6 = 	   -(c[7]*c[7])*k138*uK138;
	d138d7 = 	   c[6]*c[7]*k138*uK138*-2.0;
	d138d23 = 	   c[23]*k138*2.0;

	d139d4 = 	   -c[20]*c[24]*k139*uK139;
	d139d7 = 	   c[32]*k139;
	d139d20 = 	   -c[24]*c[4]*k139*uK139;
	d139d24 = 	   -c[20]*c[4]*k139*uK139;
	d139d32 = 	   c[7]*k139;

	d140d4 = 	   -c[27]*k140*uK140;
	d140d7 = 	   c[24]*k140;
	d140d24 = 	   c[7]*k140;
	d140d27 = 	   -c[4]*k140*uK140;

	d141d4 = 	   -c[6]*k141*uK141;
	d141d6 = 	   -c[4]*k141*uK141;
	d141d7 = 	   c[22]*k141;
	d141d22 = 	   c[7]*k141;

	d142d1 = 	   rFlat142;
	d142d2 = 	   rFlat142;
	d142d3 = 	   rFlat142;
	d142d4 = 	   rFlat142;
	d142d5 = 	   rFlat142;
	d142d6 = 	   rFlat142;
	d142d7 = 	   rFlat142+c[22]*coeffM142*k142;
	d142d8 = 	   rFlat142;
	d142d9 = 	   rFlat142;
	d142d10 = 	   rFlat142;
	d142d11 = 	   rFlat142;
	d142d12 = 	   rFlat142;
	d142d13 = 	   rFlat142;
	d142d14 = 	   rFlat142;
	d142d15 = 	   rFlat142;
	d142d16 = 	   rFlat142;
	d142d17 = 	   rFlat142;
	d142d18 = 	   rFlat142;
	d142d19 = 	   rFlat142;
	d142d20 = 	   rFlat142;
	d142d21 = 	   rFlat142;
	d142d22 = 	   rFlat142+c[7]*coeffM142*k142;
	d142d23 = 	   rFlat142-coeffM142*k142*uK142;
	d142d24 = 	   rFlat142;
	d142d25 = 	   rFlat142;
	d142d26 = 	   rFlat142;
	d142d27 = 	   rFlat142;
	d142d28 = 	   rFlat142;
	d142d29 = 	   rFlat142;
	d142d30 = 	   rFlat142;
	d142d31 = 	   rFlat142;
	d142d32 = 	   rFlat142;

	d143d6 = 	   -c[10]*k143*uK143;
	d143d7 = 	   c[28]*k143;
	d143d10 = 	   -c[6]*k143*uK143;
	d143d28 = 	   c[7]*k143;

	d144d4 = 	   c[4]*c[6]*k144*uK144*-2.0;
	d144d6 = 	   -(c[4]*c[4])*k144*uK144;
	d144d7 = 	   c[7]*k144*2.0;

	d145d4 = 	   -c[23]*k145*uK145;
	d145d7 = 	   c[7]*k145*2.0;
	d145d23 = 	   -c[4]*k145*uK145;

	d146d4 = 	   -c[23]*k146*uK146;
	d146d7 = 	   c[7]*k146*2.0;
	d146d23 = 	   -c[4]*k146*uK146;

	d147d4 = 	   -c[20]*k147*uK147;
	d147d7 = 	   c[19]*k147;
	d147d19 = 	   c[7]*k147;
	d147d20 = 	   -c[4]*k147*uK147;

	d148d7 = 	   c[32]*k148;
	d148d10 = 	   -c[19]*k148*uK148;
	d148d19 = 	   -c[10]*k148*uK148;
	d148d32 = 	   c[7]*k148;

	d149d4 = 	   -c[19]*k149*uK149;
	d149d19 = 	   -c[4]*k149*uK149;
	d149d20 = 	   c[21]*k149;
	d149d21 = 	   c[20]*k149;

	d150d1 = 	   rFlat150*1.0;
	d150d2 = 	   rFlat150*1.0;
	d150d3 = 	   rFlat150*1.7-c[22]*coeffM150*k150*uK150;
	d150d4 = 	   rFlat150*1.0;
	d150d5 = 	   rFlat150*1.0+coeffM150*k150;
	d150d6 = 	   rFlat150*1.4;
	d150d7 = 	   rFlat150*1.0;
	d150d8 = 	   rFlat150*1.0;
	d150d9 = 	   rFlat150*1.0;
	d150d10 = 	   rFlat150*1.0;
	d150d11 = 	   rFlat150*1.0;
	d150d12 = 	   rFlat150*1.0;
	d150d13 = 	   rFlat150*1.0;
	d150d14 = 	   rFlat150*1.0;
	d150d15 = 	   rFlat150*1.2E1;
	d150d16 = 	   rFlat150*1.0;
	d150d17 = 	   rFlat150*1.0;
	d150d18 = 	   rFlat150*1.0;
	d150d19 = 	   rFlat150*1.5;
	d150d20 = 	   rFlat150*3.0;
	d150d21 = 	   rFlat150*1.0;
	d150d22 = 	   rFlat150*1.0-c[3]*coeffM150*k150*uK150;
	d150d23 = 	   rFlat150*1.0;
	d150d24 = 	   rFlat150*1.0;
	d150d25 = 	   rFlat150*1.0;
	d150d26 = 	   rFlat150*1.0;
	d150d27 = 	   rFlat150*1.0;
	d150d28 = 	   rFlat150*1.0;
	d150d29 = 	   rFlat150*1.0;
	d150d30 = 	   rFlat150*1.0;
	d150d31 = 	   rFlat150*1.0;
	d150d32 = 	   rFlat150*1.0;

	d151d3 = 	   -c[27]*k151*uK151;
	d151d5 = 	   c[24]*k151;
	d151d24 = 	   c[5]*k151;
	d151d27 = 	   -c[3]*k151*uK151;

	d152d3 = 	   -c[27]*k152*uK152;
	d152d5 = 	   c[24]*k152;
	d152d24 = 	   c[5]*k152;
	d152d27 = 	   -c[3]*k152*uK152;

	d153d4 = 	   c[25]*k153;
	d153d5 = 	   -c[24]*k153*uK153;
	d153d24 = 	   -c[5]*k153*uK153;
	d153d25 = 	   c[4]*k153;

	d154d4 = 	   c[4]*k154*uK154*-2.0;
	d154d5 = 	   c[22]*k154;
	d154d22 = 	   c[5]*k154;

	d155d3 = 	   -c[6]*k155*uK155;
	d155d5 = 	   c[22]*k155;
	d155d6 = 	   -c[3]*k155*uK155;
	d155d22 = 	   c[5]*k155;

	d156d3 = 	   -c[28]*k156*uK156;
	d156d5 = 	   c[27]*k156;
	d156d27 = 	   c[5]*k156;
	d156d28 = 	   -c[3]*k156*uK156;

	d157d7 = 	   -c[13]*k157*uK157;
	d157d11 = 	   c[24]*k157;
	d157d13 = 	   -c[7]*k157*uK157;
	d157d24 = 	   c[11]*k157;

	d158d7 = 	   -c[27]*k158*uK158;
	d158d11 = 	   c[22]*k158;
	d158d22 = 	   c[11]*k158;
	d158d27 = 	   -c[7]*k158*uK158;

	d159d7 = 	   -c[15]*k159*uK159;
	d159d11 = 	   c[27]*k159;
	d159d15 = 	   -c[7]*k159*uK159;
	d159d27 = 	   c[11]*k159;

	d160d10 = 	   -k160*uK160;
	d160d11 = 	   k160;

	sigma161 =	   CFO161*dwFdM161+dCFOdM161*wF161;
	d161d1 = 	   rFlat161*sigma161*1.0;
	d161d2 = 	   rFlat161*sigma161*1.0;
	d161d3 = 	   rFlat161*sigma161*1.0;
	d161d4 = 	   rFlat161*sigma161*1.0;
	d161d5 = 	   rFlat161*sigma161*1.0;
	d161d6 = 	   rFlat161*sigma161*1.0;
	d161d7 = 	   rFlat161*sigma161*1.0+c[27]*coeffFallOff161*k161;
	d161d8 = 	   rFlat161*sigma161*1.0;
	d161d9 = 	   rFlat161*sigma161*1.0;
	d161d10 = 	   rFlat161*sigma161*1.0;
	d161d11 = 	   rFlat161*sigma161*1.0;
	d161d12 = 	   rFlat161*sigma161*1.0-coeffFallOff161*k161*uK161;
	d161d13 = 	   rFlat161*sigma161*1.0;
	d161d14 = 	   rFlat161*sigma161*1.0;
	d161d15 = 	   rFlat161*sigma161*5.0;
	d161d16 = 	   rFlat161*sigma161*1.0;
	d161d17 = 	   rFlat161*sigma161*1.0;
	d161d18 = 	   rFlat161*sigma161*1.0;
	d161d19 = 	   rFlat161*sigma161*1.0;
	d161d20 = 	   rFlat161*sigma161*1.0;
	d161d21 = 	   rFlat161*sigma161*1.0;
	d161d22 = 	   rFlat161*sigma161*1.0;
	d161d23 = 	   rFlat161*sigma161*1.0;
	d161d24 = 	   rFlat161*sigma161*1.0;
	d161d25 = 	   rFlat161*sigma161*1.0;
	d161d26 = 	   rFlat161*sigma161*1.0;
	d161d27 = 	   rFlat161*sigma161*1.0+c[7]*coeffFallOff161*k161;
	d161d28 = 	   rFlat161*sigma161*1.0;
	d161d29 = 	   rFlat161*sigma161*1.0;
	d161d30 = 	   rFlat161*sigma161*1.0;
	d161d31 = 	   rFlat161*sigma161*1.0;
	d161d32 = 	   rFlat161*sigma161*1.0;

	d162d12 = 	   c[27]*k162;
	d162d15 = 	   -c[23]*k162*uK162;
	d162d23 = 	   -c[15]*k162*uK162;
	d162d27 = 	   c[12]*k162;

	d163d13 = 	   c[24]*k163;
	d163d24 = 	   c[13]*k163;

	d164d13 = 	   c[29]*k164;
	d164d29 = 	   c[13]*k164;

	d165d15 = 	   c[27]*k165;
	d165d27 = 	   c[15]*k165;

	d166d15 = 	   c[28]*k166;
	d166d28 = 	   c[15]*k166;

	d167d15 = 	   c[29]*k167;
	d167d29 = 	   c[15]*k167;

	d168d7 = 	   c[16]*k168;
	d168d16 = 	   c[7]*k168;

	d169d16 = 	   c[27]*k169;
	d169d27 = 	   c[16]*k169;

	d170d16 = 	   c[22]*k170;
	d170d22 = 	   c[16]*k170;

	d171d16 = 	   c[28]*k171;
	d171d28 = 	   c[16]*k171;

	d172d16 = 	   c[29]*k172;
	d172d29 = 	   c[16]*k172;

	d173d6 = 	   c[17]*k173;
	d173d17 = 	   c[6]*k173;

	d174d17 = 	   c[24]*k174;
	d174d24 = 	   c[17]*k174;

	d175d17 = 	   c[27]*k175;
	d175d27 = 	   c[17]*k175;

	d176d17 = 	   c[22]*k176;
	d176d22 = 	   c[17]*k176;

	d177d17 = 	   c[28]*k177;
	d177d28 = 	   c[17]*k177;

	d178d17 = 	   c[29]*k178;
	d178d29 = 	   c[17]*k178;

	// ============================================================ 
	// ===== JACOBIAN MAXTRIX ===================================== 
	// ============================================================ 
	J[1][1] =	   0.0;
	J[1][2] =	   0.0;
	J[1][3] =	   0.0;
	J[1][4] =	   0.0;
	J[1][5] =	   0.0;
	J[1][6] =	   0.0;
	J[1][7] =	   0.0;
	J[1][8] =	   0.0;
	J[1][9] =	   0.0;
	J[1][10] =	   0.0;
	J[1][11] =	   0.0;
	J[1][12] =	   0.0;
	J[1][13] =	   0.0;
	J[1][14] =	   0.0;
	J[1][15] =	   0.0;
	J[1][16] =	   0.0;
	J[1][17] =	   0.0;
	J[1][18] =	   0.0;
	J[1][19] =	   0.0;
	J[1][20] =	   0.0;
	J[1][21] =	   0.0;
	J[1][22] =	   0.0;
	J[1][23] =	   0.0;
	J[1][24] =	   0.0;
	J[1][25] =	   0.0;
	J[1][26] =	   0.0;
	J[1][27] =	   0.0;
	J[1][28] =	   0.0;
	J[1][29] =	   0.0;
	J[1][30] =	   0.0;
	J[1][31] =	   0.0;
	J[1][32] =	   0.0;

	J[2][1] =	   0.0;
	J[2][2] =	   0.0;
	J[2][3] =	   0.0;
	J[2][4] =	   0.0;
	J[2][5] =	   0.0;
	J[2][6] =	   0.0;
	J[2][7] =	   0.0;
	J[2][8] =	   0.0;
	J[2][9] =	   0.0;
	J[2][10] =	   0.0;
	J[2][11] =	   0.0;
	J[2][12] =	   0.0;
	J[2][13] =	   0.0;
	J[2][14] =	   0.0;
	J[2][15] =	   0.0;
	J[2][16] =	   0.0;
	J[2][17] =	   0.0;
	J[2][18] =	   0.0;
	J[2][19] =	   0.0;
	J[2][20] =	   0.0;
	J[2][21] =	   0.0;
	J[2][22] =	   0.0;
	J[2][23] =	   0.0;
	J[2][24] =	   0.0;
	J[2][25] =	   0.0;
	J[2][26] =	   0.0;
	J[2][27] =	   0.0;
	J[2][28] =	   0.0;
	J[2][29] =	   0.0;
	J[2][30] =	   0.0;
	J[2][31] =	   0.0;
	J[2][32] =	   0.0;

	J[3][1] =	   d150d1+d94d1;
	J[3][2] =	   d150d2+d94d2;
	J[3][3] =	   d150d3+d151d3+d152d3+d155d3+d156d3+d41d3+d47d3+d57d3+d58d3+d59d3+d63d3+d83d3+d84d3+d87d3+d88d3+d89d3+d90d3+d91d3+d92d3+d94d3;
	J[3][4] =	   d150d4+d47d4+d59d4+d63d4+d92d4+d94d4;
	J[3][5] =	   d150d5+d151d5+d152d5+d155d5+d156d5+d94d5;
	J[3][6] =	   d150d6+d155d6+d88d6+d89d6+d94d6;
	J[3][7] =	   d150d7+d94d7;
	J[3][8] =	   d150d8+d92d8+d94d8;
	J[3][9] =	   d150d9+d94d9;
	J[3][10] =	   d150d10+d94d10;
	J[3][11] =	   d150d11+d94d11;
	J[3][12] =	   d150d12+d94d12;
	J[3][13] =	   d150d13+d84d13+d94d13;
	J[3][14] =	   d150d14+d94d14;
	J[3][15] =	   d150d15+d47d15+d87d15+d94d15;
	J[3][16] =	   d150d16+d94d16;
	J[3][17] =	   d150d17+d91d17+d94d17;
	J[3][18] =	   d150d18+d94d18;
	J[3][19] =	   d150d19+d94d19;
	J[3][20] =	   d150d20+d94d20;
	J[3][21] =	   d150d21+d41d21+d58d21+d63d21+d94d21;
	J[3][22] =	   d150d22+d155d22+d63d22+d94d22;
	J[3][23] =	   d150d23+d94d23;
	J[3][24] =	   d150d24+d151d24+d152d24+d41d24+d57d24+d58d24+d83d24+d84d24+d89d24+d94d24;
	J[3][25] =	   d150d25+d57d25+d58d25+d59d25+d90d25+d94d25;
	J[3][26] =	   d150d26+d83d26+d84d26+d87d26+d88d26+d89d26+d90d26+d91d26+d92d26+d94d26;
	J[3][27] =	   d150d27+d151d27+d152d27+d156d27+d59d27+d87d27+d94d27;
	J[3][28] =	   d150d28+d156d28+d88d28+d94d28;
	J[3][29] =	   d150d29+d41d29+d47d29+d90d29+d91d29+d94d29;
	J[3][30] =	   d150d30+d94d30;
	J[3][31] =	   d150d31+d94d31;
	J[3][32] =	   d150d32+d94d32;

	J[4][1] =	   -d106d1-d107d1-d109d1;
	J[4][2] =	   -d106d2-d107d2-d109d2;
	J[4][3] =	   -d106d3-d107d3-d109d3-d47d3-d59d3-d63d3-d92d3;
	J[4][4] =	   -d101d4-d102d4-d105d4-d106d4-d107d4-d108d4-d109d4+d110d4+d111d4+d112d4+d113d4+d114d4-d115d4+d116d4+d123d4-d129d4+d136d4+d137d4+d139d4+d140d4+d141d4+d144d4*2.0+d145d4+d146d4+d147d4+d149d4-d153d4+d154d4*2.0+d36d4-d46d4-d47d4+d49d4+d51d4+d54d4+d56d4-d59d4+d61d4+d62d4-d63d4+d77d4-d80d4+d86d4-d92d4;
	J[4][5] =	   -d101d5-d106d5-d107d5-d109d5-d115d5-d153d5+d154d5*2.0-d80d5;
	J[4][6] =	   -d106d6-d107d6-d109d6+d113d6+d136d6+d137d6+d141d6+d144d6*2.0+d56d6+d62d6;
	J[4][7] =	   -d102d7-d105d7-d106d7-d107d7-d109d7+d116d7+d123d7+d136d7+d139d7+d140d7+d141d7+d144d7*2.0+d145d7+d146d7+d147d7+d49d7;
	J[4][8] =	   -d101d8-d106d8-d107d8-d108d8-d109d8+d110d8+d111d8+d112d8+d113d8+d114d8-d115d8+d116d8-d129d8-d92d8;
	J[4][9] =	   -d101d9-d102d9-d106d9-d107d9-d109d9;
	J[4][10] =	   -d106d10-d107d10-d109d10+d116d10+d123d10;
	J[4][11] =	   -d106d11-d107d11-d109d11;
	J[4][12] =	   -d106d12-d107d12-d109d12;
	J[4][13] =	   -d106d13-d107d13-d109d13+d110d13+d36d13+d54d13;
	J[4][14] =	   -d106d14-d107d14-d109d14+d77d14-d80d14;
	J[4][15] =	   -d106d15-d107d15-d109d15+d112d15+d123d15-d47d15;
	J[4][16] =	   -d106d16-d107d16-d109d16;
	J[4][17] =	   -d106d17-d107d17-d109d17+d114d17;
	J[4][18] =	   -d106d18-d107d18-d109d18;
	J[4][19] =	   -d106d19-d107d19-d108d19-d109d19+d147d19+d149d19;
	J[4][20] =	   -d106d20-d107d20-d109d20+d139d20+d147d20+d149d20;
	J[4][21] =	   -d106d21-d107d21-d109d21+d149d21+d61d21+d62d21-d63d21;
	J[4][22] =	   -d106d22-d107d22-d109d22+d111d22+d141d22+d154d22*2.0+d36d22+d51d22+d62d22-d63d22+d77d22+d86d22;
	J[4][23] =	   -d106d23-d107d23-d109d23+d136d23+d137d23+d145d23+d146d23;
	J[4][24] =	   -d106d24-d107d24-d109d24+d110d24+d139d24+d140d24-d153d24+d51d24+d61d24;
	J[4][25] =	   -d106d25-d107d25-d109d25-d153d25+d51d25+d54d25+d56d25-d59d25+d86d25;
	J[4][26] =	   -d102d26-d106d26-d107d26-d109d26-d46d26+d86d26-d92d26;
	J[4][27] =	   -d105d27-d106d27-d107d27-d109d27+d111d27+d112d27-d115d27+d140d27-d46d27+d54d27+d56d27-d59d27+d61d27;
	J[4][28] =	   -d105d28-d106d28-d107d28-d109d28+d113d28;
	J[4][29] =	   -d106d29-d107d29-d109d29+d114d29+d36d29-d46d29-d47d29+d49d29+d77d29-d80d29;
	J[4][30] =	   -d106d30-d107d30-d109d30-d129d30+d49d30;
	J[4][31] =	   -d106d31-d107d31-d109d31;
	J[4][32] =	   -d106d32-d107d32-d108d32-d109d32+d139d32;

	J[5][1] =	   -d150d1+d93d1;
	J[5][2] =	   -d150d2+d93d2;
	J[5][3] =	   -d150d3-d151d3-d152d3-d155d3-d156d3+d93d3;
	J[5][4] =	   d101d4+d115d4-d150d4+d153d4-d154d4+d80d4+d93d4;
	J[5][5] =	   d101d5+d104d5+d115d5+d117d5-d150d5-d151d5-d152d5+d153d5-d154d5-d155d5-d156d5+d48d5+d60d5+d80d5+d85d5+d93d5+d95d5+d97d5+d99d5;
	J[5][6] =	   -d150d6-d155d6+d93d6;
	J[5][7] =	   d104d7-d150d7+d48d7+d60d7+d93d7;
	J[5][8] =	   d101d8+d115d8+d117d8-d150d8+d93d8;
	J[5][9] =	   d101d9+d104d9-d150d9+d93d9+d95d9+d97d9+d99d9;
	J[5][10] =	   d104d10-d150d10+d93d10;
	J[5][11] =	   -d150d11+d93d11;
	J[5][12] =	   -d150d12+d93d12;
	J[5][13] =	   -d150d13+d93d13+d95d13;
	J[5][14] =	   -d150d14+d80d14+d93d14;
	J[5][15] =	   d117d15-d150d15+d48d15+d93d15+d99d15;
	J[5][16] =	   -d150d16+d93d16;
	J[5][17] =	   -d150d17+d93d17;
	J[5][18] =	   -d150d18+d93d18;
	J[5][19] =	   -d150d19+d93d19;
	J[5][20] =	   -d150d20+d93d20;
	J[5][21] =	   -d150d21+d93d21;
	J[5][22] =	   -d150d22-d154d22-d155d22+d85d22+d93d22+d97d22;
	J[5][23] =	   -d150d23+d93d23;
	J[5][24] =	   -d150d24-d151d24-d152d24+d153d24+d85d24+d93d24+d95d24;
	J[5][25] =	   -d150d25+d153d25+d60d25+d93d25;
	J[5][26] =	   -d150d26+d85d26+d93d26;
	J[5][27] =	   d115d27-d150d27-d151d27-d152d27-d156d27+d60d27+d93d27+d97d27+d99d27;
	J[5][28] =	   -d150d28-d156d28+d93d28;
	J[5][29] =	   -d150d29+d48d29+d80d29+d93d29;
	J[5][30] =	   -d150d30+d93d30;
	J[5][31] =	   -d150d31+d93d31;
	J[5][32] =	   -d150d32+d93d32;

	J[6][1] =	   -d10d1-d3d1;
	J[6][2] =	   -d10d2-d3d2;
	J[6][3] =	   -d10d3+d155d3-d3d3-d88d3;
	J[6][4] =	   -d10d4-d113d4+d136d4+d137d4+d141d4+d144d4-d3d4-d56d4-d62d4;
	J[6][5] =	   -d10d5+d155d5-d3d5;
	J[6][6] =	   -d10d6-d113d6+d12d6+d13d6+d133d6+d135d6+d136d6+d137d6+d138d6+d141d6+d143d6+d144d6+d155d6-d16d6-d17d6-d173d6-d1d6-d37d6-d3d6+d40d6-d4d6-d55d6-d56d6+d5d6-d62d6+d7d6-d88d6+d98d6;
	J[6][7] =	   -d10d7+d133d7+d135d7+d136d7+d138d7+d141d7+d143d7+d144d7-d3d7;
	J[6][8] =	   -d10d8-d113d8-d37d8-d3d8-d55d8;
	J[6][9] =	   -d10d9-d3d9+d98d9;
	J[6][10] =	   -d10d10+d143d10-d3d10;
	J[6][11] =	   -d10d11-d3d11;
	J[6][12] =	   -d10d12-d3d12;
	J[6][13] =	   -d10d13+d12d13-d3d13;
	J[6][14] =	   -d10d14-d3d14;
	J[6][15] =	   -d10d15-d3d15+d5d15;
	J[6][16] =	   -d10d16+d13d16-d3d16;
	J[6][17] =	   -d10d17-d173d17-d3d17;
	J[6][18] =	   -d10d18-d3d18;
	J[6][19] =	   -d10d19-d16d19-d17d19-d3d19;
	J[6][20] =	   -d10d20-d16d20-d3d20;
	J[6][21] =	   -d10d21-d3d21-d62d21;
	J[6][22] =	   -d10d22+d133d22+d141d22+d155d22-d16d22-d1d22-d3d22+d40d22-d55d22-d62d22+d7d22+d98d22;
	J[6][23] =	   -d10d23+d133d23+d135d23+d136d23+d137d23+d138d23-d3d23;
	J[6][24] =	   -d10d24+d12d24-d1d24-d3d24-d4d24;
	J[6][25] =	   -d10d25-d3d25-d55d25-d56d25;
	J[6][26] =	   -d10d26-d3d26-d88d26+d98d26;
	J[6][27] =	   -d10d27+d135d27-d1d27-d37d27-d3d27-d56d27+d5d27+d7d27;
	J[6][28] =	   -d10d28-d113d28+d12d28+d13d28+d135d28+d143d28-d17d28-d3d28-d4d28+d5d28+d7d28-d88d28;
	J[6][29] =	   -d10d29-d37d29-d3d29+d40d29;
	J[6][30] =	   -d10d30-d3d30+d40d30;
	J[6][31] =	   -d10d31-d3d31;
	J[6][32] =	   -d10d32-d17d32-d3d32;

	J[7][1] =	   d106d1-d142d1-d161d1;
	J[7][2] =	   d106d2-d142d2-d161d2;
	J[7][3] =	   d106d3-d142d3-d161d3;
	J[7][4] =	   d102d4+d105d4+d106d4-d116d4+d123d4-d139d4-d140d4-d141d4-d142d4-d144d4*2.0-d145d4*2.0-d146d4*2.0-d147d4-d161d4-d49d4;
	J[7][5] =	   -d104d5+d106d5-d142d5-d161d5-d48d5-d60d5;
	J[7][6] =	   d106d6+d133d6+d135d6+d138d6*2.0-d141d6-d142d6-d143d6-d144d6*2.0-d161d6;
	J[7][7] =	   d102d7-d103d7-d104d7+d105d7+d106d7-d116d7-d118d7+d119d7+d120d7+d121d7+d122d7+d123d7-d131d7+d132d7+d133d7+d134d7+d135d7+d138d7*2.0-d139d7-d140d7-d141d7-d142d7-d143d7-d144d7*2.0-d145d7*2.0-d146d7*2.0-d147d7-d148d7+d157d7+d158d7+d159d7-d161d7-d168d7-d48d7-d49d7-d60d7;
	J[7][8] =	   d106d8-d116d8-d131d8-d142d8-d161d8;
	J[7][9] =	   d102d9-d103d9-d104d9+d106d9-d142d9-d161d9;
	J[7][10] =	   -d104d10+d106d10-d116d10-d118d10+d119d10+d120d10+d121d10+d122d10+d123d10-d131d10-d142d10-d143d10-d148d10-d161d10;
	J[7][11] =	   d106d11-d142d11+d157d11+d158d11+d159d11-d161d11;
	J[7][12] =	   d106d12-d142d12-d161d12;
	J[7][13] =	   d106d13-d118d13-d142d13+d157d13-d161d13;
	J[7][14] =	   d106d14-d142d14-d161d14;
	J[7][15] =	   d106d15+d120d15+d123d15-d142d15+d159d15-d161d15-d48d15;
	J[7][16] =	   d106d16-d142d16-d161d16-d168d16;
	J[7][17] =	   d106d17+d122d17-d142d17-d161d17;
	J[7][18] =	   d106d18-d142d18-d161d18;
	J[7][19] =	   d106d19-d142d19-d147d19-d148d19-d161d19;
	J[7][20] =	   d106d20-d139d20-d142d20-d147d20-d161d20;
	J[7][21] =	   d106d21-d142d21-d161d21;
	J[7][22] =	   d106d22+d119d22+d133d22-d141d22-d142d22+d158d22-d161d22;
	J[7][23] =	   -d103d23+d106d23+d132d23+d133d23+d134d23+d135d23+d138d23*2.0-d142d23-d145d23*2.0-d146d23*2.0-d161d23;
	J[7][24] =	   d106d24-d118d24+d132d24-d139d24-d140d24-d142d24+d157d24-d161d24;
	J[7][25] =	   d106d25+d121d25-d142d25-d161d25-d60d25;
	J[7][26] =	   d102d26-d103d26+d106d26-d142d26-d161d26;
	J[7][27] =	   d105d27+d106d27+d119d27+d120d27+d132d27+d134d27+d135d27-d140d27-d142d27+d158d27+d159d27-d161d27-d60d27;
	J[7][28] =	   d105d28+d106d28+d134d28+d135d28-d142d28-d143d28-d161d28;
	J[7][29] =	   d106d29+d121d29+d122d29-d142d29-d161d29-d48d29-d49d29;
	J[7][30] =	   d106d30-d131d30-d142d30-d161d30-d49d30;
	J[7][31] =	   d106d31-d142d31-d161d31;
	J[7][32] =	   d106d32-d139d32-d142d32-d148d32-d161d32;

	J[8][1] =	   d109d1+d124d1;
	J[8][2] =	   d109d2+d124d2;
	J[8][3] =	   d109d3+d124d3+d92d3;
	J[8][4] =	   d101d4+d108d4+d109d4-d110d4-d111d4-d112d4-d113d4-d114d4-d115d4-d116d4+d124d4+d129d4*2.0+d92d4;
	J[8][5] =	   d101d5+d109d5-d115d5-d117d5*2.0+d124d5;
	J[8][6] =	   d109d6-d113d6+d124d6+d37d6+d55d6;
	J[8][7] =	   d109d7-d116d7+d124d7+d131d7;
	J[8][8] =	   d101d8+d108d8+d109d8-d110d8-d111d8-d112d8-d113d8-d114d8-d115d8-d116d8-d117d8*2.0+d124d8+d125d8+d127d8+d128d8+d129d8*2.0+d130d8+d131d8+d34d8+d37d8+d52d8+d55d8+d71d8+d73d8+d92d8;
	J[8][9] =	   d101d9+d109d9+d124d9;
	J[8][10] =	   d109d10-d116d10+d124d10+d131d10;
	J[8][11] =	   d109d11+d124d11;
	J[8][12] =	   d109d12+d124d12;
	J[8][13] =	   d109d13-d110d13+d124d13+d125d13;
	J[8][14] =	   d109d14+d124d14;
	J[8][15] =	   d109d15-d112d15-d117d15*2.0+d124d15+d128d15;
	J[8][16] =	   d109d16+d124d16;
	J[8][17] =	   d109d17-d114d17+d124d17+d130d17+d73d17;
	J[8][18] =	   d109d18+d124d18;
	J[8][19] =	   d108d19+d109d19+d124d19;
	J[8][20] =	   d109d20+d124d20;
	J[8][21] =	   d109d21+d124d21;
	J[8][22] =	   d109d22-d111d22+d124d22+d127d22+d34d22+d55d22+d71d22;
	J[8][23] =	   d109d23+d124d23;
	J[8][24] =	   d109d24-d110d24+d124d24+d125d24+d34d24+d52d24;
	J[8][25] =	   d109d25+d124d25+d52d25+d55d25;
	J[8][26] =	   d109d26+d124d26+d92d26;
	J[8][27] =	   d109d27-d111d27-d112d27-d115d27+d124d27+d127d27+d128d27+d37d27+d52d27+d73d27;
	J[8][28] =	   d109d28-d113d28+d124d28;
	J[8][29] =	   d109d29-d114d29+d124d29+d130d29+d34d29+d37d29+d71d29;
	J[8][30] =	   d109d30+d124d30+d125d30+d127d30+d128d30+d129d30*2.0+d130d30+d131d30;
	J[8][31] =	   d109d31+d124d31+d71d31+d73d31;
	J[8][32] =	   d108d32+d109d32+d124d32;

	J[9][1] =	   -d93d1-d94d1;
	J[9][2] =	   -d93d2-d94d2;
	J[9][3] =	   -d93d3-d94d3;
	J[9][4] =	   -d101d4-d102d4-d93d4-d94d4;
	J[9][5] =	   -d101d5-d104d5-d93d5-d94d5-d95d5-d97d5-d99d5;
	J[9][6] =	   -d93d6-d94d6-d98d6;
	J[9][7] =	   -d102d7-d103d7-d104d7-d93d7-d94d7;
	J[9][8] =	   -d101d8-d93d8-d94d8;
	J[9][9] =	   -d100d9-d101d9-d102d9-d103d9-d104d9-d93d9-d94d9-d95d9-d96d9-d97d9-d98d9-d99d9;
	J[9][10] =	   -d104d10-d93d10-d94d10;
	J[9][11] =	   -d93d11-d94d11;
	J[9][12] =	   -d93d12-d94d12;
	J[9][13] =	   -d93d13-d94d13-d95d13;
	J[9][14] =	   -d93d14-d94d14;
	J[9][15] =	   -d93d15-d94d15-d99d15;
	J[9][16] =	   -d93d16-d94d16;
	J[9][17] =	   -d93d17-d94d17;
	J[9][18] =	   -d93d18-d94d18;
	J[9][19] =	   -d93d19-d94d19;
	J[9][20] =	   -d93d20-d94d20;
	J[9][21] =	   -d93d21-d94d21;
	J[9][22] =	   -d93d22-d94d22-d97d22-d98d22;
	J[9][23] =	   -d103d23-d93d23-d94d23;
	J[9][24] =	   -d93d24-d94d24-d95d24-d96d24;
	J[9][25] =	   -d93d25-d94d25;
	J[9][26] =	   -d100d26-d102d26-d103d26-d93d26-d94d26-d96d26-d98d26;
	J[9][27] =	   -d100d27-d93d27-d94d27-d96d27-d97d27-d99d27;
	J[9][28] =	   -d100d28-d93d28-d94d28;
	J[9][29] =	   -d93d29-d94d29;
	J[9][30] =	   -d93d30-d94d30;
	J[9][31] =	   -d93d31-d94d31;
	J[9][32] =	   -d93d32-d94d32;

	J[10][1] =	   d107d1;
	J[10][2] =	   d107d2;
	J[10][3] =	   d107d3;
	J[10][4] =	   d107d4+d116d4-d123d4*2.0;
	J[10][5] =	   d104d5+d107d5;
	J[10][6] =	   d107d6+d143d6;
	J[10][7] =	   d104d7+d107d7+d116d7+d118d7-d119d7-d120d7-d121d7-d122d7-d123d7*2.0+d131d7+d143d7+d148d7+d168d7;
	J[10][8] =	   d107d8+d116d8+d131d8;
	J[10][9] =	   d104d9+d107d9;
	J[10][10] =	   d104d10+d107d10+d116d10+d118d10-d119d10-d120d10-d121d10-d122d10-d123d10*2.0+d131d10+d143d10+d148d10+d160d10;
	J[10][11] =	   d107d11+d160d11;
	J[10][12] =	   d107d12;
	J[10][13] =	   d107d13+d118d13;
	J[10][14] =	   d107d14;
	J[10][15] =	   d107d15-d120d15-d123d15*2.0;
	J[10][16] =	   d107d16+d168d16;
	J[10][17] =	   d107d17-d122d17;
	J[10][18] =	   d107d18;
	J[10][19] =	   d107d19+d148d19;
	J[10][20] =	   d107d20;
	J[10][21] =	   d107d21;
	J[10][22] =	   d107d22-d119d22;
	J[10][23] =	   d107d23;
	J[10][24] =	   d107d24+d118d24;
	J[10][25] =	   d107d25-d121d25;
	J[10][26] =	   d107d26;
	J[10][27] =	   d107d27-d119d27-d120d27;
	J[10][28] =	   d107d28+d143d28;
	J[10][29] =	   d107d29-d121d29-d122d29;
	J[10][30] =	   d107d30+d131d30;
	J[10][31] =	   d107d31;
	J[10][32] =	   d107d32+d148d32;

	J[11][1] =	   0.0;
	J[11][2] =	   0.0;
	J[11][3] =	   0.0;
	J[11][4] =	   0.0;
	J[11][5] =	   0.0;
	J[11][6] =	   0.0;
	J[11][7] =	   -d157d7-d158d7-d159d7;
	J[11][8] =	   0.0;
	J[11][9] =	   0.0;
	J[11][10] =	   -d160d10;
	J[11][11] =	   -d157d11-d158d11-d159d11-d160d11;
	J[11][12] =	   0.0;
	J[11][13] =	   -d157d13;
	J[11][14] =	   0.0;
	J[11][15] =	   -d159d15;
	J[11][16] =	   0.0;
	J[11][17] =	   0.0;
	J[11][18] =	   0.0;
	J[11][19] =	   0.0;
	J[11][20] =	   0.0;
	J[11][21] =	   0.0;
	J[11][22] =	   -d158d22;
	J[11][23] =	   0.0;
	J[11][24] =	   -d157d24;
	J[11][25] =	   0.0;
	J[11][26] =	   0.0;
	J[11][27] =	   -d158d27-d159d27;
	J[11][28] =	   0.0;
	J[11][29] =	   0.0;
	J[11][30] =	   0.0;
	J[11][31] =	   0.0;
	J[11][32] =	   0.0;

	J[12][1] =	   d161d1;
	J[12][2] =	   d161d2;
	J[12][3] =	   d161d3;
	J[12][4] =	   d161d4;
	J[12][5] =	   d161d5;
	J[12][6] =	   d161d6;
	J[12][7] =	   d161d7;
	J[12][8] =	   d161d8;
	J[12][9] =	   d161d9;
	J[12][10] =	   d161d10;
	J[12][11] =	   d161d11;
	J[12][12] =	   d161d12-d162d12;
	J[12][13] =	   d161d13;
	J[12][14] =	   d161d14;
	J[12][15] =	   d161d15-d162d15;
	J[12][16] =	   d161d16;
	J[12][17] =	   d161d17;
	J[12][18] =	   d161d18;
	J[12][19] =	   d161d19;
	J[12][20] =	   d161d20;
	J[12][21] =	   d161d21;
	J[12][22] =	   d161d22;
	J[12][23] =	   d161d23-d162d23;
	J[12][24] =	   d161d24;
	J[12][25] =	   d161d25;
	J[12][26] =	   d161d26;
	J[12][27] =	   d161d27-d162d27;
	J[12][28] =	   d161d28;
	J[12][29] =	   d161d29;
	J[12][30] =	   d161d30;
	J[12][31] =	   d161d31;
	J[12][32] =	   d161d32;

	J[13][1] =	   -d9d1;
	J[13][2] =	   -d9d2;
	J[13][3] =	   d84d3-d9d3;
	J[13][4] =	   d110d4+d36d4+d54d4-d9d4;
	J[13][5] =	   d95d5-d9d5;
	J[13][6] =	   d12d6-d9d6;
	J[13][7] =	   -d118d7+d157d7-d9d7;
	J[13][8] =	   d110d8+d125d8-d9d8;
	J[13][9] =	   d95d9-d9d9;
	J[13][10] =	   -d118d10-d9d10;
	J[13][11] =	   d157d11-d9d11;
	J[13][12] =	   -d9d12;
	J[13][13] =	   d110d13-d118d13+d12d13+d125d13+d157d13-d164d13+d22d13+d25d13+d29d13-d2d13+d31d13+d33d13+d36d13+d43d13+d50d13+d54d13+d64d13+d76d13+d84d13+d95d13-d9d13;
	J[13][14] =	   d43d14+d76d14-d9d14;
	J[13][15] =	   d22d15+d29d15-d9d15;
	J[13][16] =	   d31d16-d9d16;
	J[13][17] =	   d174d17-d9d17;
	J[13][18] =	   d64d18-d9d18;
	J[13][19] =	   d22d19+d25d19-d9d19;
	J[13][20] =	   d22d20-d9d20;
	J[13][21] =	   d50d21-d9d21;
	J[13][22] =	   -d2d22+d36d22-d9d22;
	J[13][23] =	   -d9d23;
	J[13][24] =	   d110d24-d118d24+d12d24+d125d24+d157d24+d174d24+d25d24+d29d24-d2d24+d31d24+d33d24+d50d24+d64d24+d76d24+d84d24+d95d24-d9d24;
	J[13][25] =	   d33d25+d50d25+d54d25-d9d25;
	J[13][26] =	   d76d26+d84d26-d9d26;
	J[13][27] =	   d29d27-d2d27+d54d27-d9d27;
	J[13][28] =	   d12d28+d31d28-d9d28;
	J[13][29] =	   -d164d29+d33d29+d36d29+d43d29-d9d29;
	J[13][30] =	   d125d30-d9d30;
	J[13][31] =	   d64d31-d9d31;
	J[13][32] =	   d25d32-d9d32;

	J[14][1] =	   d68d1-d75d1;
	J[14][2] =	   d68d2-d75d2;
	J[14][3] =	   d68d3-d75d3;
	J[14][4] =	   d68d4-d75d4-d77d4-d80d4;
	J[14][5] =	   d68d5-d75d5-d80d5;
	J[14][6] =	   d68d6-d75d6;
	J[14][7] =	   d68d7-d75d7;
	J[14][8] =	   d68d8-d75d8;
	J[14][9] =	   d68d9-d75d9;
	J[14][10] =	   d68d10-d75d10;
	J[14][11] =	   d68d11-d75d11;
	J[14][12] =	   d68d12-d75d12;
	J[14][13] =	   d43d13+d68d13-d75d13-d76d13;
	J[14][14] =	   d42d14+d43d14+d65d14+d68d14+d70d14+d72d14+d74d14-d75d14-d76d14-d77d14-d78d14-d79d14-d80d14-d81d14-d82d14;
	J[14][15] =	   d65d15+d68d15+d72d15-d75d15-d79d15;
	J[14][16] =	   d68d16-d75d16;
	J[14][17] =	   d68d17-d75d17-d82d17;
	J[14][18] =	   d65d18+d68d18-d75d18;
	J[14][19] =	   d68d19-d75d19;
	J[14][20] =	   d68d20-d75d20;
	J[14][21] =	   d68d21-d75d21;
	J[14][22] =	   d65d22+d68d22+d70d22-d75d22-d77d22-d78d22;
	J[14][23] =	   d68d23-d75d23;
	J[14][24] =	   d42d24+d68d24-d75d24-d76d24;
	J[14][25] =	   d42d25+d68d25+d74d25-d75d25-d81d25;
	J[14][26] =	   d68d26-d75d26-d76d26-d78d26-d79d26-d81d26-d82d26;
	J[14][27] =	   d68d27+d70d27+d72d27-d75d27-d78d27-d79d27;
	J[14][28] =	   d68d28-d75d28;
	J[14][29] =	   d42d29+d43d29+d68d29+d74d29-d75d29-d77d29-d80d29-d81d29-d82d29;
	J[14][30] =	   d68d30-d75d30;
	J[14][31] =	   d68d31+d70d31+d72d31+d74d31-d75d31;
	J[14][32] =	   d68d32-d75d32;

	J[15][1] =	   d11d1;
	J[15][2] =	   d11d2;
	J[15][3] =	   d11d3+d47d3+d87d3;
	J[15][4] =	   d11d4+d112d4+d123d4+d47d4;
	J[15][5] =	   d11d5+d117d5+d48d5+d99d5;
	J[15][6] =	   d11d6+d5d6;
	J[15][7] =	   d11d7+d120d7+d123d7+d159d7+d48d7;
	J[15][8] =	   d11d8+d112d8+d117d8+d128d8;
	J[15][9] =	   d11d9+d99d9;
	J[15][10] =	   d11d10+d120d10+d123d10;
	J[15][11] =	   d11d11+d159d11;
	J[15][12] =	   d11d12+d162d12;
	J[15][13] =	   d11d13-d22d13-d29d13;
	J[15][14] =	   d11d14+d65d14+d72d14+d79d14;
	J[15][15] =	   d11d15+d112d15+d117d15+d120d15+d123d15+d128d15+d159d15+d162d15-d166d15-d167d15-d22d15+d26d15-d29d15+d30d15+d38d15+d47d15+d48d15+d53d15+d5d15+d65d15+d66d15+d72d15+d79d15+d87d15+d8d15+d99d15;
	J[15][16] =	   d11d16+d169d16+d30d16;
	J[15][17] =	   d11d17+d175d17;
	J[15][18] =	   d11d18+d65d18+d66d18;
	J[15][19] =	   d11d19-d22d19+d26d19;
	J[15][20] =	   d11d20-d22d20;
	J[15][21] =	   d11d21+d53d21;
	J[15][22] =	   d11d22+d65d22+d8d22;
	J[15][23] =	   d11d23+d162d23;
	J[15][24] =	   d11d24-d29d24+d30d24;
	J[15][25] =	   d11d25+d38d25+d53d25;
	J[15][26] =	   d11d26+d79d26+d87d26;
	J[15][27] =	   d11d27+d112d27+d120d27+d128d27+d159d27+d162d27+d169d27+d175d27+d26d27-d29d27+d30d27+d38d27+d53d27+d5d27+d66d27+d72d27+d79d27+d87d27+d8d27+d99d27;
	J[15][28] =	   d11d28-d166d28+d5d28;
	J[15][29] =	   d11d29-d167d29+d38d29+d47d29+d48d29;
	J[15][30] =	   d11d30+d128d30;
	J[15][31] =	   d11d31+d66d31+d72d31;
	J[15][32] =	   d11d32+d26d32;

	J[16][1] =	   d14d1;
	J[16][2] =	   d14d2;
	J[16][3] =	   d14d3;
	J[16][4] =	   d14d4;
	J[16][5] =	   d14d5;
	J[16][6] =	   d13d6+d14d6;
	J[16][7] =	   d14d7-d168d7;
	J[16][8] =	   d14d8;
	J[16][9] =	   d14d9;
	J[16][10] =	   d14d10;
	J[16][11] =	   d14d11;
	J[16][12] =	   d14d12;
	J[16][13] =	   d14d13-d31d13;
	J[16][14] =	   d14d14;
	J[16][15] =	   d14d15+d166d15-d30d15;
	J[16][16] =	   d13d16+d14d16-d168d16-d169d16-d170d16-d172d16+d27d16-d30d16-d31d16;
	J[16][17] =	   d14d17+d177d17;
	J[16][18] =	   d14d18;
	J[16][19] =	   d14d19+d27d19;
	J[16][20] =	   d14d20;
	J[16][21] =	   d14d21;
	J[16][22] =	   d14d22-d170d22;
	J[16][23] =	   d14d23;
	J[16][24] =	   d14d24-d30d24-d31d24;
	J[16][25] =	   d14d25;
	J[16][26] =	   d14d26;
	J[16][27] =	   d14d27-d169d27-d30d27;
	J[16][28] =	   d13d28+d14d28+d166d28+d177d28+d27d28-d31d28;
	J[16][29] =	   d14d29-d172d29;
	J[16][30] =	   d14d30;
	J[16][31] =	   d14d31;
	J[16][32] =	   d14d32+d27d32;

	J[17][1] =	   -d32d1;
	J[17][2] =	   -d32d2;
	J[17][3] =	   -d32d3+d91d3;
	J[17][4] =	   d114d4-d32d4;
	J[17][5] =	   -d32d5;
	J[17][6] =	   -d173d6-d32d6;
	J[17][7] =	   d122d7-d32d7;
	J[17][8] =	   d114d8+d130d8-d32d8+d73d8;
	J[17][9] =	   -d32d9;
	J[17][10] =	   d122d10-d32d10;
	J[17][11] =	   -d32d11;
	J[17][12] =	   -d32d12;
	J[17][13] =	   d164d13-d32d13;
	J[17][14] =	   -d32d14+d82d14;
	J[17][15] =	   d167d15-d32d15;
	J[17][16] =	   d172d16-d32d16;
	J[17][17] =	   d114d17+d122d17+d130d17-d173d17-d174d17-d175d17-d176d17-d177d17-d32d17+d44d17+d67d17+d73d17+d82d17+d91d17;
	J[17][18] =	   -d32d18+d67d18;
	J[17][19] =	   -d32d19;
	J[17][20] =	   -d32d20;
	J[17][21] =	   -d32d21;
	J[17][22] =	   -d176d22-d32d22;
	J[17][23] =	   -d32d23;
	J[17][24] =	   -d174d24-d32d24;
	J[17][25] =	   -d32d25+d44d25;
	J[17][26] =	   -d32d26+d82d26+d91d26;
	J[17][27] =	   -d175d27-d32d27+d73d27;
	J[17][28] =	   -d177d28-d32d28;
	J[17][29] =	   d114d29+d122d29+d130d29+d164d29+d167d29+d172d29-d32d29+d44d29+d67d29+d82d29+d91d29;
	J[17][30] =	   d130d30-d32d30;
	J[17][31] =	   -d32d31+d67d31+d73d31;
	J[17][32] =	   -d32d32;

	J[18][1] =	   d45d1;
	J[18][2] =	   d45d2;
	J[18][3] =	   d45d3;
	J[18][4] =	   d45d4;
	J[18][5] =	   d45d5;
	J[18][6] =	   d45d6;
	J[18][7] =	   d45d7;
	J[18][8] =	   d45d8;
	J[18][9] =	   d45d9;
	J[18][10] =	   d45d10;
	J[18][11] =	   d45d11;
	J[18][12] =	   d45d12;
	J[18][13] =	   d45d13-d64d13;
	J[18][14] =	   d45d14-d65d14;
	J[18][15] =	   d45d15-d65d15-d66d15;
	J[18][16] =	   d45d16;
	J[18][17] =	   d45d17-d67d17;
	J[18][18] =	   d45d18-d64d18-d65d18-d66d18-d67d18;
	J[18][19] =	   d45d19;
	J[18][20] =	   d45d20;
	J[18][21] =	   d45d21;
	J[18][22] =	   d45d22-d65d22;
	J[18][23] =	   d45d23;
	J[18][24] =	   d45d24-d64d24;
	J[18][25] =	   d45d25;
	J[18][26] =	   d45d26;
	J[18][27] =	   d45d27-d66d27;
	J[18][28] =	   d45d28;
	J[18][29] =	   d45d29-d67d29;
	J[18][30] =	   d45d30;
	J[18][31] =	   d45d31-d64d31-d66d31-d67d31;
	J[18][32] =	   d45d32;

	J[19][1] =	   -d18d1+d23d1;
	J[19][2] =	   -d18d2+d23d2;
	J[19][3] =	   -d18d3+d23d3;
	J[19][4] =	   d108d4-d147d4+d149d4-d18d4+d23d4;
	J[19][5] =	   -d18d5+d23d5;
	J[19][6] =	   -d16d6+d17d6-d18d6+d23d6;
	J[19][7] =	   -d147d7+d148d7-d18d7+d23d7;
	J[19][8] =	   d108d8-d18d8+d23d8;
	J[19][9] =	   -d18d9+d23d9;
	J[19][10] =	   d148d10-d18d10+d23d10;
	J[19][11] =	   -d18d11+d23d11;
	J[19][12] =	   -d18d12+d23d12;
	J[19][13] =	   -d18d13-d22d13+d23d13+d25d13;
	J[19][14] =	   -d18d14+d23d14;
	J[19][15] =	   -d18d15-d22d15+d23d15+d26d15;
	J[19][16] =	   -d18d16+d23d16+d27d16;
	J[19][17] =	   -d18d17+d23d17;
	J[19][18] =	   -d18d18+d23d18;
	J[19][19] =	   d108d19-d147d19+d148d19+d149d19-d16d19+d17d19-d18d19-d19d19-d20d19-d21d19-d22d19+d23d19+d25d19+d26d19+d27d19;
	J[19][20] =	   -d147d20+d149d20-d16d20-d18d20-d19d20-d20d20-d21d20-d22d20+d23d20;
	J[19][21] =	   d149d21-d18d21+d23d21;
	J[19][22] =	   -d16d22-d18d22+d23d22;
	J[19][23] =	   -d18d23+d23d23;
	J[19][24] =	   -d18d24-d19d24-d20d24+d23d24+d25d24;
	J[19][25] =	   -d18d25+d23d25;
	J[19][26] =	   -d18d26+d23d26;
	J[19][27] =	   -d18d27-d19d27-d20d27-d21d27+d23d27+d26d27;
	J[19][28] =	   d17d28-d18d28-d21d28+d23d28+d27d28;
	J[19][29] =	   -d18d29+d23d29;
	J[19][30] =	   -d18d30+d23d30;
	J[19][31] =	   -d18d31+d23d31;
	J[19][32] =	   d108d32+d148d32+d17d32-d18d32+d23d32+d25d32+d26d32+d27d32;

	J[20][1] =	   d18d1;
	J[20][2] =	   d18d2;
	J[20][3] =	   d18d3;
	J[20][4] =	   d139d4+d147d4-d149d4+d18d4;
	J[20][5] =	   d18d5;
	J[20][6] =	   d16d6+d18d6;
	J[20][7] =	   d139d7+d147d7+d18d7;
	J[20][8] =	   d18d8;
	J[20][9] =	   d18d9;
	J[20][10] =	   d18d10;
	J[20][11] =	   d18d11;
	J[20][12] =	   d18d12;
	J[20][13] =	   d18d13+d22d13;
	J[20][14] =	   d18d14;
	J[20][15] =	   d18d15+d22d15;
	J[20][16] =	   d18d16;
	J[20][17] =	   d18d17;
	J[20][18] =	   d18d18;
	J[20][19] =	   d147d19-d149d19+d16d19+d18d19+d19d19+d20d19+d21d19+d22d19;
	J[20][20] =	   d139d20+d147d20-d149d20+d16d20+d18d20+d19d20+d20d20+d21d20+d22d20+d24d20;
	J[20][21] =	   -d149d21+d18d21;
	J[20][22] =	   d16d22+d18d22+d24d22;
	J[20][23] =	   d18d23;
	J[20][24] =	   d139d24+d18d24+d19d24+d20d24+d24d24;
	J[20][25] =	   d18d25;
	J[20][26] =	   d18d26;
	J[20][27] =	   d18d27+d19d27+d20d27+d21d27;
	J[20][28] =	   d18d28+d21d28+d28d28;
	J[20][29] =	   d18d29;
	J[20][30] =	   d18d30;
	J[20][31] =	   d18d31;
	J[20][32] =	   d139d32+d18d32+d24d32+d28d32;

	J[21][1] =	   0.0;
	J[21][2] =	   0.0;
	J[21][3] =	   -d41d3-d58d3-d63d3;
	J[21][4] =	   -d149d4-d61d4-d62d4-d63d4;
	J[21][5] =	   0.0;
	J[21][6] =	   -d62d6;
	J[21][7] =	   0.0;
	J[21][8] =	   0.0;
	J[21][9] =	   0.0;
	J[21][10] =	   0.0;
	J[21][11] =	   0.0;
	J[21][12] =	   0.0;
	J[21][13] =	   d50d13;
	J[21][14] =	   0.0;
	J[21][15] =	   d53d15;
	J[21][16] =	   0.0;
	J[21][17] =	   0.0;
	J[21][18] =	   0.0;
	J[21][19] =	   -d149d19;
	J[21][20] =	   -d149d20;
	J[21][21] =	   -d149d21-d41d21+d50d21+d53d21-d58d21-d61d21-d62d21-d63d21;
	J[21][22] =	   -d62d22-d63d22;
	J[21][23] =	   0.0;
	J[21][24] =	   -d41d24+d50d24-d58d24-d61d24;
	J[21][25] =	   d50d25+d53d25-d58d25;
	J[21][26] =	   0.0;
	J[21][27] =	   d53d27-d61d27;
	J[21][28] =	   0.0;
	J[21][29] =	   -d41d29;
	J[21][30] =	   0.0;
	J[21][31] =	   0.0;
	J[21][32] =	   0.0;

	J[22][1] =	   d10d1*2.0-d106d1-d142d1-d15d1+d150d1-d18d1;
	J[22][2] =	   d10d2*2.0-d106d2-d142d2-d15d2+d150d2-d18d2;
	J[22][3] =	   d10d3*2.0-d106d3-d142d3-d15d3+d150d3-d155d3-d18d3+d63d3;
	J[22][4] =	   d10d4*2.0-d106d4-d111d4-d141d4-d142d4-d15d4+d150d4-d154d4-d18d4-d36d4-d51d4+d62d4+d63d4-d77d4-d86d4;
	J[22][5] =	   d10d5*2.0-d106d5-d142d5-d15d5+d150d5-d154d5-d155d5-d18d5-d85d5-d97d5;
	J[22][6] =	   d10d6*2.0-d106d6-d133d6-d141d6-d142d6-d15d6+d150d6-d155d6+d16d6-d18d6+d1d6-d40d6+d55d6+d62d6-d7d6-d98d6;
	J[22][7] =	   d10d7*2.0-d106d7-d119d7-d133d7-d141d7-d142d7-d15d7+d150d7-d158d7-d18d7;
	J[22][8] =	   d10d8*2.0-d106d8-d111d8-d127d8-d142d8-d15d8+d150d8-d18d8-d34d8+d55d8-d71d8;
	J[22][9] =	   d10d9*2.0-d106d9-d142d9-d15d9+d150d9-d18d9-d97d9-d98d9;
	J[22][10] =	   d10d10*2.0-d106d10-d119d10-d142d10-d15d10+d150d10-d18d10;
	J[22][11] =	   d10d11*2.0-d106d11-d142d11-d15d11+d150d11-d158d11-d18d11;
	J[22][12] =	   d10d12*2.0-d106d12-d142d12-d15d12+d150d12-d18d12;
	J[22][13] =	   d10d13*2.0-d106d13-d142d13-d15d13+d150d13-d18d13-d2d13-d36d13;
	J[22][14] =	   d10d14*2.0-d106d14-d142d14-d15d14+d150d14-d18d14-d65d14-d70d14-d77d14-d78d14;
	J[22][15] =	   d10d15*2.0-d106d15-d142d15-d15d15+d150d15-d18d15-d65d15+d8d15;
	J[22][16] =	   d10d16*2.0-d106d16-d142d16-d15d16+d150d16-d170d16-d18d16;
	J[22][17] =	   d10d17*2.0-d106d17-d142d17-d15d17+d150d17-d176d17-d18d17;
	J[22][18] =	   d10d18*2.0-d106d18-d142d18-d15d18+d150d18-d18d18-d65d18;
	J[22][19] =	   d10d19*2.0-d106d19-d142d19-d15d19+d150d19+d16d19-d18d19;
	J[22][20] =	   d10d20*2.0-d106d20-d142d20-d15d20+d150d20+d16d20-d18d20-d24d20;
	J[22][21] =	   d10d21*2.0-d106d21-d142d21-d15d21+d150d21-d18d21+d62d21+d63d21;
	J[22][22] =	   d10d22*2.0-d106d22-d111d22-d119d22-d127d22-d133d22-d141d22-d142d22-d15d22+d150d22-d154d22-d155d22-d158d22+d16d22-d170d22-d176d22-d18d22+d1d22-d24d22-d2d22-d34d22-d35d22-d36d22-d40d22-d51d22+d55d22+d62d22+d63d22-d65d22-d70d22-d71d22-d77d22-d78d22-d7d22-d85d22-d86d22+d8d22-d97d22-d98d22;
	J[22][23] =	   d10d23*2.0-d106d23-d133d23-d142d23-d15d23+d150d23-d18d23;
	J[22][24] =	   d10d24*2.0-d106d24-d142d24-d15d24+d150d24-d18d24+d1d24-d24d24-d2d24-d34d24-d51d24-d85d24;
	J[22][25] =	   d10d25*2.0-d106d25-d142d25-d15d25+d150d25-d18d25-d35d25-d51d25+d55d25-d86d25;
	J[22][26] =	   d10d26*2.0-d106d26-d142d26-d15d26+d150d26-d18d26-d78d26-d85d26-d86d26-d98d26;
	J[22][27] =	   d10d27*2.0-d106d27-d111d27-d119d27-d127d27-d142d27-d15d27+d150d27-d158d27-d18d27+d1d27-d2d27-d35d27-d70d27-d78d27-d7d27+d8d27-d97d27;
	J[22][28] =	   d10d28*2.0-d106d28-d142d28-d15d28+d150d28-d18d28-d7d28;
	J[22][29] =	   d10d29*2.0-d106d29-d142d29-d15d29+d150d29-d18d29-d34d29-d35d29-d36d29-d40d29-d71d29-d77d29;
	J[22][30] =	   d10d30*2.0-d106d30-d127d30-d142d30-d15d30+d150d30-d18d30-d40d30;
	J[22][31] =	   d10d31*2.0-d106d31-d142d31-d15d31+d150d31-d18d31-d70d31-d71d31;
	J[22][32] =	   d10d32*2.0-d106d32-d142d32-d15d32+d150d32-d18d32-d24d32;

	J[23][1] =	   d142d1;
	J[23][2] =	   d142d2;
	J[23][3] =	   d142d3;
	J[23][4] =	   -d136d4-d137d4+d142d4+d145d4+d146d4;
	J[23][5] =	   d142d5;
	J[23][6] =	   -d133d6-d135d6-d136d6-d137d6-d138d6*2.0+d142d6;
	J[23][7] =	   d103d7-d132d7-d133d7-d134d7-d135d7-d136d7-d138d7*2.0+d142d7+d145d7+d146d7;
	J[23][8] =	   d142d8;
	J[23][9] =	   d103d9+d142d9;
	J[23][10] =	   d142d10;
	J[23][11] =	   d142d11;
	J[23][12] =	   d142d12+d162d12;
	J[23][13] =	   d142d13;
	J[23][14] =	   d142d14;
	J[23][15] =	   d142d15+d162d15;
	J[23][16] =	   d142d16;
	J[23][17] =	   d142d17;
	J[23][18] =	   d142d18;
	J[23][19] =	   d142d19;
	J[23][20] =	   d142d20;
	J[23][21] =	   d142d21;
	J[23][22] =	   -d133d22+d142d22;
	J[23][23] =	   d103d23-d132d23-d133d23-d134d23-d135d23-d136d23-d137d23-d138d23*2.0+d142d23+d145d23+d146d23+d162d23;
	J[23][24] =	   -d132d24+d142d24;
	J[23][25] =	   d142d25;
	J[23][26] =	   d103d26+d142d26;
	J[23][27] =	   -d132d27-d134d27-d135d27+d142d27+d162d27;
	J[23][28] =	   -d134d28-d135d28+d142d28;
	J[23][29] =	   d142d29;
	J[23][30] =	   d142d30;
	J[23][31] =	   d142d31;
	J[23][32] =	   d142d32;

	J[24][1] =	   -d109d1-d11d1+d124d1+d23d1+d32d1-d3d1+d68d1+d75d1+d93d1+d9d1*2.0;
	J[24][2] =	   -d109d2-d11d2+d124d2+d23d2+d32d2-d3d2+d68d2+d75d2+d93d2+d9d2*2.0;
	J[24][3] =	   -d109d3-d11d3+d124d3-d151d3-d152d3+d23d3+d32d3-d3d3+d41d3*2.0+d57d3*2.0+d58d3+d68d3+d75d3+d83d3-d84d3+d89d3+d93d3+d9d3*2.0;
	J[24][4] =	   -d109d4-d11d4-d110d4+d124d4+d139d4-d140d4+d153d4+d23d4+d32d4-d3d4+d51d4+d61d4+d68d4+d75d4+d93d4+d9d4*2.0;
	J[24][5] =	   -d109d5-d11d5+d124d5-d151d5-d152d5+d153d5+d23d5+d32d5-d3d5+d68d5+d75d5+d85d5+d93d5-d95d5+d9d5*2.0;
	J[24][6] =	   -d109d6-d11d6-d12d6+d124d6-d1d6+d23d6+d32d6-d3d6-d4d6+d68d6+d75d6+d89d6+d93d6+d9d6*2.0;
	J[24][7] =	   -d109d7-d11d7+d118d7+d124d7-d132d7+d139d7-d140d7-d157d7+d23d7+d32d7-d3d7+d68d7+d75d7+d93d7+d9d7*2.0;
	J[24][8] =	   -d109d8-d11d8-d110d8+d124d8-d125d8+d23d8+d32d8+d34d8-d3d8+d52d8+d68d8+d75d8+d93d8+d9d8*2.0;
	J[24][9] =	   -d109d9-d11d9+d124d9+d23d9+d32d9-d3d9+d68d9+d75d9+d93d9-d95d9-d96d9+d9d9*2.0;
	J[24][10] =	   -d109d10-d11d10+d118d10+d124d10+d23d10+d32d10-d3d10+d68d10+d75d10+d93d10+d9d10*2.0;
	J[24][11] =	   -d109d11-d11d11+d124d11-d157d11+d23d11+d32d11-d3d11+d68d11+d75d11+d93d11+d9d11*2.0;
	J[24][12] =	   -d109d12-d11d12+d124d12+d23d12+d32d12-d3d12+d68d12+d75d12+d93d12+d9d12*2.0;
	J[24][13] =	   -d109d13-d11d13-d110d13+d118d13-d12d13+d124d13-d125d13-d157d13+d164d13+d23d13-d25d13-d29d13+d2d13-d31d13+d32d13-d33d13-d3d13-d50d13-d64d13+d68d13+d75d13-d76d13-d84d13+d93d13-d95d13+d9d13*2.0;
	J[24][14] =	   -d109d14-d11d14+d124d14+d23d14+d32d14-d3d14+d42d14+d68d14+d75d14-d76d14+d93d14+d9d14*2.0;
	J[24][15] =	   -d109d15-d11d15+d124d15+d23d15-d29d15-d30d15+d32d15-d3d15+d68d15+d75d15+d93d15+d9d15*2.0;
	J[24][16] =	   -d109d16-d11d16+d124d16+d23d16-d30d16-d31d16+d32d16-d3d16+d68d16+d75d16+d93d16+d9d16*2.0;
	J[24][17] =	   -d109d17-d11d17+d124d17-d174d17+d23d17+d32d17-d3d17+d68d17+d75d17+d93d17+d9d17*2.0;
	J[24][18] =	   -d109d18-d11d18+d124d18+d23d18+d32d18-d3d18-d64d18+d68d18+d75d18+d93d18+d9d18*2.0;
	J[24][19] =	   -d109d19-d11d19+d124d19+d19d19+d20d19+d23d19-d25d19+d32d19-d3d19+d68d19+d75d19+d93d19+d9d19*2.0;
	J[24][20] =	   -d109d20-d11d20+d124d20+d139d20+d19d20+d20d20+d23d20+d24d20+d32d20-d3d20+d68d20+d75d20+d93d20+d9d20*2.0;
	J[24][21] =	   -d109d21-d11d21+d124d21+d23d21+d32d21-d3d21+d41d21*2.0-d50d21+d58d21+d61d21+d68d21+d75d21+d93d21+d9d21*2.0;
	J[24][22] =	   -d109d22-d11d22+d124d22-d1d22+d23d22+d24d22+d2d22+d32d22+d34d22-d3d22+d51d22+d68d22+d75d22+d85d22+d93d22+d9d22*2.0;
	J[24][23] =	   -d109d23-d11d23+d124d23-d132d23+d23d23+d32d23-d3d23+d68d23+d75d23+d93d23+d9d23*2.0;
	J[24][24] =	   -d109d24-d11d24-d110d24+d118d24-d12d24+d124d24-d125d24-d126d24-d132d24+d139d24-d140d24-d151d24-d152d24+d153d24-d157d24-d174d24+d19d24-d1d24+d20d24+d23d24+d24d24-d25d24-d29d24+d2d24-d30d24-d31d24+d32d24-d33d24+d34d24-d3d24+d41d24*2.0+d42d24-d4d24-d50d24+d51d24+d52d24+d57d24*2.0+d58d24+d61d24-d64d24+d68d24-d69d24-d6d24+d75d24-d76d24+d83d24-d84d24+d85d24+d89d24+d93d24-d95d24-d96d24+d9d24*2.0;
	J[24][25] =	   -d109d25-d11d25+d124d25+d153d25+d23d25+d32d25-d33d25-d3d25+d42d25-d50d25+d51d25+d52d25+d57d25*2.0+d58d25+d68d25+d75d25+d93d25+d9d25*2.0;
	J[24][26] =	   -d109d26-d11d26+d124d26+d23d26+d32d26-d3d26+d68d26+d75d26-d76d26+d83d26-d84d26+d85d26+d89d26+d93d26-d96d26+d9d26*2.0;
	J[24][27] =	   -d109d27-d11d27+d124d27-d126d27-d132d27-d140d27-d151d27-d152d27+d19d27-d1d27+d20d27+d23d27-d29d27+d2d27-d30d27+d32d27-d3d27+d52d27+d61d27+d68d27-d6d27+d75d27+d93d27-d96d27+d9d27*2.0;
	J[24][28] =	   -d109d28-d11d28-d12d28+d124d28+d23d28+d28d28-d31d28+d32d28-d3d28-d4d28+d68d28-d6d28+d75d28+d93d28+d9d28*2.0;
	J[24][29] =	   -d109d29-d11d29+d124d29-d126d29+d164d29+d23d29+d32d29-d33d29+d34d29-d3d29+d41d29*2.0+d42d29+d68d29-d69d29+d75d29+d93d29+d9d29*2.0;
	J[24][30] =	   -d109d30-d11d30+d124d30-d125d30-d126d30+d23d30+d32d30-d3d30+d68d30+d75d30+d93d30+d9d30*2.0;
	J[24][31] =	   -d109d31-d11d31+d124d31+d23d31+d32d31-d3d31-d64d31+d68d31-d69d31+d75d31+d93d31+d9d31*2.0;
	J[24][32] =	   -d109d32-d11d32+d124d32+d139d32+d23d32+d24d32-d25d32+d28d32+d32d32-d3d32+d68d32+d75d32+d93d32+d9d32*2.0;

	J[25][1] =	   0.0;
	J[25][2] =	   0.0;
	J[25][3] =	   d57d3*-2.0-d58d3-d59d3-d90d3;
	J[25][4] =	   -d153d4-d51d4-d54d4-d56d4-d59d4+d86d4;
	J[25][5] =	   -d153d5-d60d5;
	J[25][6] =	   -d55d6-d56d6;
	J[25][7] =	   -d121d7-d60d7;
	J[25][8] =	   -d52d8-d55d8;
	J[25][9] =	   0.0;
	J[25][10] =	   -d121d10;
	J[25][11] =	   0.0;
	J[25][12] =	   0.0;
	J[25][13] =	   d33d13-d50d13-d54d13;
	J[25][14] =	   -d42d14-d74d14-d81d14;
	J[25][15] =	   d38d15-d53d15;
	J[25][16] =	   0.0;
	J[25][17] =	   d44d17;
	J[25][18] =	   0.0;
	J[25][19] =	   0.0;
	J[25][20] =	   0.0;
	J[25][21] =	   -d50d21-d53d21-d58d21;
	J[25][22] =	   d35d22-d51d22-d55d22+d86d22;
	J[25][23] =	   0.0;
	J[25][24] =	   -d153d24+d33d24-d42d24-d50d24-d51d24-d52d24-d57d24*2.0-d58d24;
	J[25][25] =	   -d121d25-d153d25+d33d25+d35d25+d38d25-d42d25+d44d25-d50d25-d51d25-d52d25-d53d25-d54d25-d55d25-d56d25-d57d25*2.0-d58d25-d59d25-d60d25-d74d25-d81d25+d86d25-d90d25;
	J[25][26] =	   -d81d26+d86d26-d90d26;
	J[25][27] =	   d35d27+d38d27-d52d27-d53d27-d54d27-d56d27-d59d27-d60d27;
	J[25][28] =	   0.0;
	J[25][29] =	   -d121d29+d33d29+d35d29+d38d29-d42d29+d44d29-d74d29-d81d29-d90d29;
	J[25][30] =	   0.0;
	J[25][31] =	   -d74d31;
	J[25][32] =	   0.0;

	J[26][1] =	   d75d1;
	J[26][2] =	   d75d2;
	J[26][3] =	   d75d3-d83d3-d84d3-d87d3-d88d3-d89d3-d90d3-d91d3-d92d3;
	J[26][4] =	   d102d4+d46d4+d75d4-d86d4-d92d4;
	J[26][5] =	   d75d5-d85d5;
	J[26][6] =	   d75d6-d88d6-d89d6+d98d6;
	J[26][7] =	   d102d7+d103d7+d75d7;
	J[26][8] =	   d75d8-d92d8;
	J[26][9] =	   d100d9+d102d9+d103d9+d75d9+d96d9+d98d9;
	J[26][10] =	   d75d10;
	J[26][11] =	   d75d11;
	J[26][12] =	   d75d12;
	J[26][13] =	   d75d13+d76d13-d84d13;
	J[26][14] =	   d75d14+d76d14+d78d14+d79d14+d81d14+d82d14;
	J[26][15] =	   d75d15+d79d15-d87d15;
	J[26][16] =	   d75d16;
	J[26][17] =	   d75d17+d82d17-d91d17;
	J[26][18] =	   d75d18;
	J[26][19] =	   d75d19;
	J[26][20] =	   d75d20;
	J[26][21] =	   d75d21;
	J[26][22] =	   d75d22+d78d22-d85d22-d86d22+d98d22;
	J[26][23] =	   d103d23+d75d23;
	J[26][24] =	   d75d24+d76d24-d83d24-d84d24-d85d24-d89d24+d96d24;
	J[26][25] =	   d75d25+d81d25-d86d25-d90d25;
	J[26][26] =	   d100d26+d102d26+d103d26+d46d26+d75d26+d76d26+d78d26+d79d26+d81d26+d82d26-d83d26-d84d26-d85d26-d86d26-d87d26-d88d26-d89d26-d90d26-d91d26-d92d26+d96d26+d98d26;
	J[26][27] =	   d100d27+d46d27+d75d27+d78d27+d79d27-d87d27+d96d27;
	J[26][28] =	   d100d28+d75d28-d88d28;
	J[26][29] =	   d46d29+d75d29+d81d29+d82d29-d90d29-d91d29;
	J[26][30] =	   d75d30;
	J[26][31] =	   d75d31;
	J[26][32] =	   d75d32;

	J[27][1] =	   -d107d1-d11d1-d14d1*2.0-d15d1-d161d1+d94d1;
	J[27][2] =	   -d107d2-d11d2-d14d2*2.0-d15d2-d161d2+d94d2;
	J[27][3] =	   -d107d3-d11d3-d14d3*2.0-d15d3+d151d3+d152d3-d156d3-d161d3+d59d3-d87d3+d94d3;
	J[27][4] =	   d105d4-d107d4-d11d4+d111d4-d112d4+d115d4-d14d4*2.0+d140d4-d15d4-d161d4+d46d4-d54d4+d56d4+d59d4-d61d4+d94d4;
	J[27][5] =	   -d107d5-d11d5+d115d5-d14d5*2.0-d15d5+d151d5+d152d5-d156d5-d161d5+d60d5+d94d5+d97d5-d99d5;
	J[27][6] =	   -d107d6-d11d6+d135d6-d14d6*2.0-d15d6-d161d6+d1d6+d37d6+d56d6-d5d6+d7d6+d94d6;
	J[27][7] =	   d105d7-d107d7-d11d7+d119d7-d120d7+d132d7-d134d7+d135d7-d14d7*2.0+d140d7-d15d7+d158d7-d159d7-d161d7+d60d7+d94d7;
	J[27][8] =	   -d107d8-d11d8+d111d8-d112d8+d115d8+d127d8-d128d8-d14d8*2.0-d15d8-d161d8+d37d8-d52d8-d73d8+d94d8;
	J[27][9] =	   -d100d9-d107d9-d11d9-d14d9*2.0-d15d9-d161d9+d94d9+d96d9+d97d9-d99d9;
	J[27][10] =	   -d107d10-d11d10+d119d10-d120d10-d14d10*2.0-d15d10-d161d10+d94d10;
	J[27][11] =	   -d107d11-d11d11-d14d11*2.0-d15d11+d158d11-d159d11-d161d11+d94d11;
	J[27][12] =	   -d107d12-d11d12-d14d12*2.0-d15d12-d161d12-d162d12+d94d12;
	J[27][13] =	   -d107d13-d11d13-d14d13*2.0-d15d13-d161d13+d29d13+d2d13-d54d13+d94d13;
	J[27][14] =	   -d107d14-d11d14-d14d14*2.0-d15d14-d161d14+d70d14-d72d14+d78d14-d79d14+d94d14;
	J[27][15] =	   -d107d15-d11d15-d112d15-d120d15-d128d15-d14d15*2.0-d15d15-d159d15-d161d15-d162d15+d166d15+d167d15-d26d15+d29d15+d30d15-d38d15-d53d15-d5d15-d66d15-d72d15-d79d15-d87d15-d8d15*2.0+d94d15-d99d15;
	J[27][16] =	   -d107d16-d11d16-d14d16*2.0-d15d16-d161d16-d169d16+d170d16+d30d16+d94d16;
	J[27][17] =	   -d107d17-d11d17-d14d17*2.0-d15d17-d161d17-d175d17+d176d17-d73d17+d94d17;
	J[27][18] =	   -d107d18-d11d18-d14d18*2.0-d15d18-d161d18-d66d18+d94d18;
	J[27][19] =	   -d107d19-d11d19-d14d19*2.0-d15d19-d161d19-d19d19-d20d19+d21d19-d26d19+d94d19;
	J[27][20] =	   -d107d20-d11d20-d14d20*2.0-d15d20-d161d20-d19d20-d20d20+d21d20+d94d20;
	J[27][21] =	   -d107d21-d11d21-d14d21*2.0-d15d21-d161d21-d53d21-d61d21+d94d21;
	J[27][22] =	   -d107d22-d11d22+d111d22+d119d22+d127d22-d14d22*2.0-d15d22+d158d22-d161d22+d170d22+d176d22+d1d22+d2d22+d35d22+d70d22+d78d22+d7d22-d8d22*2.0+d94d22+d97d22;
	J[27][23] =	   -d107d23-d11d23+d132d23-d134d23+d135d23-d14d23*2.0-d15d23-d161d23-d162d23+d94d23;
	J[27][24] =	   -d107d24-d11d24+d126d24+d132d24-d14d24*2.0+d140d24-d15d24+d151d24+d152d24-d161d24-d19d24+d1d24-d20d24+d29d24+d2d24+d30d24-d52d24-d61d24+d6d24*2.0+d94d24+d96d24;
	J[27][25] =	   -d107d25-d11d25-d14d25*2.0-d15d25-d161d25+d35d25-d38d25-d52d25-d53d25-d54d25+d56d25+d59d25+d60d25+d94d25;
	J[27][26] =	   -d100d26-d107d26-d11d26-d14d26*2.0-d15d26-d161d26+d46d26+d78d26-d79d26-d87d26+d94d26+d96d26;
	J[27][27] =	   -d100d27+d105d27-d107d27-d11d27+d111d27-d112d27+d115d27+d119d27-d120d27+d126d27+d127d27-d128d27+d132d27-d134d27+d135d27-d14d27*2.0+d140d27-d15d27+d151d27+d152d27-d156d27+d158d27-d159d27-d161d27-d162d27-d169d27-d175d27-d19d27+d1d27-d20d27+d21d27-d26d27+d29d27+d2d27+d30d27+d35d27+d37d27-d38d27+d39d27+d46d27-d52d27-d53d27-d54d27+d56d27+d59d27-d5d27+d60d27-d61d27-d66d27+d6d27*2.0+d70d27-d72d27-d73d27+d78d27-d79d27+d7d27-d87d27-d8d27*2.0+d94d27+d96d27+d97d27-d99d27;
	J[27][28] =	   -d100d28+d105d28-d107d28-d11d28-d134d28+d135d28-d14d28*2.0-d15d28-d156d28-d161d28+d166d28+d21d28+d28d28+d39d28-d5d28+d6d28*2.0+d7d28+d94d28;
	J[27][29] =	   -d107d29-d11d29+d126d29-d14d29*2.0-d15d29-d161d29+d167d29+d35d29+d37d29-d38d29+d39d29+d46d29+d94d29;
	J[27][30] =	   -d107d30-d11d30+d126d30+d127d30-d128d30-d14d30*2.0-d15d30-d161d30+d39d30+d94d30;
	J[27][31] =	   -d107d31-d11d31-d14d31*2.0-d15d31-d161d31-d66d31+d70d31-d72d31-d73d31+d94d31;
	J[27][32] =	   -d107d32-d11d32-d14d32*2.0-d15d32-d161d32-d26d32+d28d32+d94d32;

	J[28][1] =	   d15d1+d3d1;
	J[28][2] =	   d15d2+d3d2;
	J[28][3] =	   d15d3+d156d3+d3d3+d88d3;
	J[28][4] =	   -d105d4+d113d4+d15d4+d3d4;
	J[28][5] =	   d15d5+d156d5+d3d5;
	J[28][6] =	   d113d6-d12d6-d13d6*2.0-d135d6-d143d6+d15d6+d17d6+d173d6+d3d6+d4d6-d5d6-d7d6+d88d6;
	J[28][7] =	   -d105d7+d134d7-d135d7-d143d7+d15d7+d168d7+d3d7;
	J[28][8] =	   d113d8+d15d8+d3d8;
	J[28][9] =	   d100d9+d15d9+d3d9;
	J[28][10] =	   -d143d10+d15d10+d3d10;
	J[28][11] =	   d15d11+d3d11;
	J[28][12] =	   d15d12+d3d12;
	J[28][13] =	   -d12d13+d15d13+d31d13+d3d13;
	J[28][14] =	   d15d14+d3d14;
	J[28][15] =	   d15d15-d166d15+d3d15-d5d15;
	J[28][16] =	   d13d16*-2.0+d15d16+d168d16+d169d16+d170d16+d172d16-d27d16+d31d16+d3d16;
	J[28][17] =	   d15d17+d173d17-d177d17+d3d17;
	J[28][18] =	   d15d18+d3d18;
	J[28][19] =	   d15d19+d17d19-d21d19-d27d19+d3d19;
	J[28][20] =	   d15d20-d21d20+d3d20;
	J[28][21] =	   d15d21+d3d21;
	J[28][22] =	   d15d22+d170d22+d3d22-d7d22;
	J[28][23] =	   d134d23-d135d23+d15d23+d3d23;
	J[28][24] =	   -d12d24+d15d24+d31d24+d3d24+d4d24-d6d24;
	J[28][25] =	   d15d25+d3d25;
	J[28][26] =	   d100d26+d15d26+d3d26+d88d26;
	J[28][27] =	   d100d27-d105d27+d134d27-d135d27+d15d27+d156d27+d169d27-d21d27-d39d27+d3d27-d5d27-d6d27-d7d27;
	J[28][28] =	   d100d28-d105d28+d113d28-d12d28-d13d28*2.0+d134d28-d135d28-d143d28+d15d28+d156d28-d166d28+d17d28-d177d28-d21d28-d27d28-d28d28+d31d28-d39d28+d3d28+d4d28-d5d28-d6d28-d7d28+d88d28;
	J[28][29] =	   d15d29+d172d29-d39d29+d3d29;
	J[28][30] =	   d15d30-d39d30+d3d30;
	J[28][31] =	   d15d31+d3d31;
	J[28][32] =	   d15d32+d17d32-d27d32-d28d32+d3d32;

	J[29][1] =	   d32d1-d45d1*2.0;
	J[29][2] =	   d32d2-d45d2*2.0;
	J[29][3] =	   d32d3-d41d3-d45d3*2.0-d47d3+d90d3-d91d3;
	J[29][4] =	   -d114d4+d32d4-d36d4-d45d4*2.0-d46d4-d47d4-d49d4+d77d4+d80d4;
	J[29][5] =	   d32d5-d45d5*2.0-d48d5+d80d5;
	J[29][6] =	   d173d6+d32d6-d37d6+d40d6-d45d6*2.0;
	J[29][7] =	   d121d7-d122d7+d32d7-d45d7*2.0-d48d7-d49d7;
	J[29][8] =	   -d114d8-d130d8+d32d8-d34d8-d37d8-d45d8*2.0+d71d8;
	J[29][9] =	   d32d9-d45d9*2.0;
	J[29][10] =	   d121d10-d122d10+d32d10-d45d10*2.0;
	J[29][11] =	   d32d11-d45d11*2.0;
	J[29][12] =	   d32d12-d45d12*2.0;
	J[29][13] =	   -d164d13+d32d13-d33d13-d36d13-d43d13*2.0-d45d13*2.0;
	J[29][14] =	   d32d14-d42d14-d43d14*2.0-d45d14*2.0+d74d14+d77d14+d80d14+d81d14-d82d14;
	J[29][15] =	   -d167d15+d32d15-d38d15-d45d15*2.0-d47d15-d48d15;
	J[29][16] =	   -d172d16+d32d16-d45d16*2.0;
	J[29][17] =	   -d114d17-d122d17-d130d17+d173d17+d174d17+d175d17+d176d17+d177d17+d32d17-d44d17*2.0-d45d17*2.0-d67d17-d82d17-d91d17;
	J[29][18] =	   d32d18-d45d18*2.0-d67d18;
	J[29][19] =	   d32d19-d45d19*2.0;
	J[29][20] =	   d32d20-d45d20*2.0;
	J[29][21] =	   d32d21-d41d21-d45d21*2.0;
	J[29][22] =	   d176d22+d32d22-d34d22-d35d22-d36d22+d40d22-d45d22*2.0+d71d22+d77d22;
	J[29][23] =	   d32d23-d45d23*2.0;
	J[29][24] =	   d126d24+d174d24+d32d24-d33d24-d34d24-d41d24-d42d24-d45d24*2.0+d69d24*2.0;
	J[29][25] =	   d121d25+d32d25-d33d25-d35d25-d38d25-d42d25-d44d25*2.0-d45d25*2.0+d74d25+d81d25+d90d25;
	J[29][26] =	   d32d26-d45d26*2.0-d46d26+d81d26-d82d26+d90d26-d91d26;
	J[29][27] =	   d126d27+d175d27+d32d27-d35d27-d37d27-d38d27-d39d27-d45d27*2.0-d46d27;
	J[29][28] =	   d177d28+d32d28-d39d28-d45d28*2.0;
	J[29][29] =	   -d114d29+d121d29-d122d29+d126d29-d130d29-d164d29-d167d29-d172d29+d32d29-d33d29-d34d29-d35d29-d36d29-d37d29-d38d29-d39d29+d40d29-d41d29-d42d29-d43d29*2.0-d44d29*2.0-d45d29*2.0-d46d29-d47d29-d48d29-d49d29-d67d29+d69d29*2.0+d71d29+d74d29+d77d29+d80d29+d81d29-d82d29+d90d29-d91d29;
	J[29][30] =	   d126d30-d130d30+d32d30-d39d30+d40d30-d45d30*2.0-d49d30;
	J[29][31] =	   d32d31-d45d31*2.0-d67d31+d69d31*2.0+d71d31+d74d31;
	J[29][32] =	   d32d32-d45d32*2.0;

	J[30][1] =	   -d124d1;
	J[30][2] =	   -d124d2;
	J[30][3] =	   -d124d3;
	J[30][4] =	   -d124d4-d129d4+d49d4;
	J[30][5] =	   -d124d5;
	J[30][6] =	   -d124d6-d40d6;
	J[30][7] =	   -d124d7-d131d7+d49d7;
	J[30][8] =	   -d124d8-d125d8-d127d8-d128d8-d129d8-d130d8-d131d8;
	J[30][9] =	   -d124d9;
	J[30][10] =	   -d124d10-d131d10;
	J[30][11] =	   -d124d11;
	J[30][12] =	   -d124d12;
	J[30][13] =	   -d124d13-d125d13;
	J[30][14] =	   -d124d14;
	J[30][15] =	   -d124d15-d128d15;
	J[30][16] =	   -d124d16;
	J[30][17] =	   -d124d17-d130d17;
	J[30][18] =	   -d124d18;
	J[30][19] =	   -d124d19;
	J[30][20] =	   -d124d20;
	J[30][21] =	   -d124d21;
	J[30][22] =	   -d124d22-d127d22-d40d22;
	J[30][23] =	   -d124d23;
	J[30][24] =	   -d124d24-d125d24-d126d24;
	J[30][25] =	   -d124d25;
	J[30][26] =	   -d124d26;
	J[30][27] =	   -d124d27-d126d27-d127d27-d128d27+d39d27;
	J[30][28] =	   -d124d28+d39d28;
	J[30][29] =	   -d124d29-d126d29-d130d29+d39d29-d40d29+d49d29;
	J[30][30] =	   -d124d30-d125d30-d126d30-d127d30-d128d30-d129d30-d130d30-d131d30+d39d30-d40d30+d49d30;
	J[30][31] =	   -d124d31;
	J[30][32] =	   -d124d32;

	J[31][1] =	   -d68d1;
	J[31][2] =	   -d68d2;
	J[31][3] =	   -d68d3;
	J[31][4] =	   -d68d4;
	J[31][5] =	   -d68d5;
	J[31][6] =	   -d68d6;
	J[31][7] =	   -d68d7;
	J[31][8] =	   -d68d8-d71d8-d73d8;
	J[31][9] =	   -d68d9;
	J[31][10] =	   -d68d10;
	J[31][11] =	   -d68d11;
	J[31][12] =	   -d68d12;
	J[31][13] =	   d64d13-d68d13;
	J[31][14] =	   -d68d14-d70d14-d72d14-d74d14;
	J[31][15] =	   d66d15-d68d15-d72d15;
	J[31][16] =	   -d68d16;
	J[31][17] =	   d67d17-d68d17-d73d17;
	J[31][18] =	   d64d18+d66d18+d67d18-d68d18;
	J[31][19] =	   -d68d19;
	J[31][20] =	   -d68d20;
	J[31][21] =	   -d68d21;
	J[31][22] =	   -d68d22-d70d22-d71d22;
	J[31][23] =	   -d68d23;
	J[31][24] =	   d64d24-d68d24-d69d24;
	J[31][25] =	   -d68d25-d74d25;
	J[31][26] =	   -d68d26;
	J[31][27] =	   d66d27-d68d27-d70d27-d72d27-d73d27;
	J[31][28] =	   -d68d28;
	J[31][29] =	   d67d29-d68d29-d69d29-d71d29-d74d29;
	J[31][30] =	   -d68d30;
	J[31][31] =	   d64d31+d66d31+d67d31-d68d31-d69d31-d70d31-d71d31-d72d31-d73d31-d74d31;
	J[31][32] =	   -d68d32;

	J[32][1] =	   -d23d1;
	J[32][2] =	   -d23d2;
	J[32][3] =	   -d23d3;
	J[32][4] =	   -d108d4-d139d4-d23d4;
	J[32][5] =	   -d23d5;
	J[32][6] =	   -d17d6-d23d6;
	J[32][7] =	   -d139d7-d148d7-d23d7;
	J[32][8] =	   -d108d8-d23d8;
	J[32][9] =	   -d23d9;
	J[32][10] =	   -d148d10-d23d10;
	J[32][11] =	   -d23d11;
	J[32][12] =	   -d23d12;
	J[32][13] =	   -d23d13-d25d13;
	J[32][14] =	   -d23d14;
	J[32][15] =	   -d23d15-d26d15;
	J[32][16] =	   -d23d16-d27d16;
	J[32][17] =	   -d23d17;
	J[32][18] =	   -d23d18;
	J[32][19] =	   -d108d19-d148d19-d17d19-d23d19-d25d19-d26d19-d27d19;
	J[32][20] =	   -d139d20-d23d20-d24d20;
	J[32][21] =	   -d23d21;
	J[32][22] =	   -d23d22-d24d22;
	J[32][23] =	   -d23d23;
	J[32][24] =	   -d139d24-d23d24-d24d24-d25d24;
	J[32][25] =	   -d23d25;
	J[32][26] =	   -d23d26;
	J[32][27] =	   -d23d27-d26d27;
	J[32][28] =	   -d17d28-d23d28-d27d28-d28d28;
	J[32][29] =	   -d23d29;
	J[32][30] =	   -d23d30;
	J[32][31] =	   -d23d31;
	J[32][32] =	   -d108d32-d139d32-d148d32-d17d32-d23d32-d24d32-d25d32-d26d32-d27d32-d28d32;
}