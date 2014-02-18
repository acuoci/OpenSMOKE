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

#include "symbolickinetics/gri12/OpenSMOKE_SymbolicKinetics_GRI12.h"

void OpenSMOKE_SymbolicKinetics_GRI12::giveReactionRates(double cTot, BzzVector &c, BzzVector &R) 
{

	// ============================================================ 
	// ===== CORRECTION COEFFICIENTS FOR THIRD BODY REACTIONS ===== 
	// ============================================================ 
	coeffM1 =	   c[1]*1.4+c[14]*1.0+c[15]*7.5E-1+c[16]*2.6+c[27]*2.0-c[32]*1.7E-1+c[6]*1.44E1+cTot;
	coeffM2 =	   c[1]*1.0+c[14]*1.0+c[15]*5.0E-1+c[16]*1.0+c[27]*2.0-c[32]*3.0E-1+c[6]*5.0+cTot;
	coeffM12 =	   c[1]*1.0+c[14]*1.0+c[15]*5.0E-1+c[16]*2.5+c[27]*2.0-c[32]*5.0E-1+c[4]*5.0+c[6]*5.0+cTot;
	coeffM33 =	   c[15]*-2.5E-1+c[16]*5.0E-1+c[27]*5.0E-1-c[31]*1.0-c[32]*1.0-c[4]*1.0-c[6]*1.0+cTot;
	coeffM39 =	   c[1]*-1.0+c[14]*1.0-c[16]*1.0+c[27]*2.0-c[32]*3.7E-1-c[6]*1.0+cTot;
	coeffM43 =	   c[1]*-2.7E-1+c[14]*1.0+c[27]*2.0-c[32]*6.2E-1+c[6]*2.65+cTot;
	coeffM50 =	   c[1]*1.0+c[14]*1.0+c[15]*5.0E-1+c[16]*1.0+c[27]*2.0-c[32]*3.0E-1+c[6]*5.0+cTot;
	coeffM52 =	   c[1]*1.0+c[14]*1.0+c[15]*5.0E-1+c[16]*1.0+c[27]*2.0-c[32]*3.0E-1+c[6]*5.0+cTot;
	coeffM54 =	   c[1]*1.0+c[14]*1.0+c[15]*5.0E-1+c[16]*1.0+c[27]*2.0-c[32]*3.0E-1+c[6]*5.0+cTot;
	coeffM56 =	   c[1]*1.0+c[14]*1.0+c[15]*5.0E-1+c[16]*1.0+c[27]*2.0+c[6]*5.0+cTot;
	coeffM57 =	   c[1]*1.0+c[14]*1.0+c[15]*5.0E-1+c[16]*1.0+c[27]*2.0+c[6]*5.0+cTot;
	coeffM59 =	   c[1]*1.0+c[14]*1.0+c[15]*5.0E-1+c[16]*1.0+c[27]*2.0+c[6]*5.0+cTot;
	coeffM63 =	   c[1]*1.0+c[14]*1.0+c[15]*5.0E-1+c[16]*1.0+c[27]*2.0+c[6]*5.0+cTot;
	coeffM70 =	   c[1]*1.0+c[14]*1.0+c[15]*5.0E-1+c[16]*1.0+c[27]*2.0-c[32]*3.0E-1+c[6]*5.0+cTot;
	coeffM71 =	   c[1]*1.0+c[14]*1.0+c[15]*5.0E-1+c[16]*1.0+c[27]*2.0-c[32]*3.0E-1+c[6]*5.0+cTot;
	coeffM72 =	   c[1]*1.0+c[14]*1.0+c[15]*5.0E-1+c[16]*1.0+c[27]*2.0-c[32]*3.0E-1+c[6]*5.0+cTot;
	coeffM74 =	   c[1]*1.0+c[14]*1.0+c[15]*5.0E-1+c[16]*1.0+c[27]*2.0-c[32]*3.0E-1+c[6]*5.0+cTot;
	coeffM76 =	   c[1]*1.0+c[14]*1.0+c[15]*5.0E-1+c[16]*1.0+c[27]*2.0-c[32]*3.0E-1+c[6]*5.0+cTot;
	coeffM83 =	   c[1]*1.0+c[14]*1.0+c[15]*5.0E-1+c[16]*1.0+c[27]*2.0-c[32]*3.0E-1+c[6]*5.0+cTot;
	coeffM85 =	   c[1]*1.0+c[14]*1.0+c[15]*5.0E-1+c[16]*1.0+c[27]*2.0-c[32]*3.0E-1+c[6]*5.0+cTot;
	coeffM95 =	   c[1]*1.0+c[14]*1.0+c[15]*5.0E-1+c[16]*1.0+c[27]*2.0+c[6]*5.0+cTot;
	coeffM131 =	   c[1]*1.0+c[14]*1.0+c[15]*5.0E-1+c[16]*1.0+c[27]*2.0-c[32]*3.0E-1+c[6]*5.0+cTot;
	coeffM140 =	   c[1]*1.0+c[14]*1.0+c[15]*5.0E-1+c[16]*1.0+c[27]*2.0-c[32]*3.0E-1+c[6]*5.0+cTot;
	coeffM147 =	   c[1]*1.0+c[14]*1.0+c[15]*5.0E-1+c[16]*1.0+c[27]*2.0+c[6]*5.0+cTot;
	coeffM158 =	   c[1]*1.0+c[14]*1.0+c[15]*5.0E-1+c[16]*1.0+c[27]*2.0-c[32]*3.0E-1+c[6]*5.0+cTot;
	coeffM167 =	   c[1]*1.0+c[14]*1.0+c[15]*5.0E-1+c[16]*1.0+c[27]*2.0-c[6]*1.0+cTot;
	coeffM174 =	   c[1]*1.0+c[14]*1.0+c[15]*5.0E-1+c[16]*1.0+c[27]*2.0-c[32]*3.0E-1+c[6]*5.0+cTot;


	// ============================================================ 
	// ===== CORRECTION COEFFICIENTS FOR FALL OFF REACTIONS ======= 
	// ============================================================ 
	CFO50 =	   coeffM50/((coeffM50*k50)/kFallOff50+1.0);
	CFO52 =	   coeffM52/((coeffM52*k52)/kFallOff52+1.0);
	CFO54 =	   coeffM54/((coeffM54*k54)/kFallOff54+1.0);
	CFO56 =	   coeffM56/((coeffM56*k56)/kFallOff56+1.0);
	CFO57 =	   coeffM57/((coeffM57*k57)/kFallOff57+1.0);
	CFO59 =	   coeffM59/((coeffM59*k59)/kFallOff59+1.0);
	CFO63 =	   coeffM63/((coeffM63*k63)/kFallOff63+1.0);
	CFO70 =	   coeffM70/((coeffM70*k70)/kFallOff70+1.0);
	CFO71 =	   coeffM71/((coeffM71*k71)/kFallOff71+1.0);
	CFO72 =	   coeffM72/((coeffM72*k72)/kFallOff72+1.0);
	CFO74 =	   coeffM74/((coeffM74*k74)/kFallOff74+1.0);
	CFO76 =	   coeffM76/((coeffM76*k76)/kFallOff76+1.0);
	CFO83 =	   coeffM83/((coeffM83*k83)/kFallOff83+1.0);
	CFO85 =	   coeffM85/((coeffM85*k85)/kFallOff85+1.0);
	CFO95 =	   coeffM95/((coeffM95*k95)/kFallOff95+1.0);
	CFO131 =	   coeffM131/((coeffM131*k131)/kFallOff131+1.0);
	CFO140 =	   coeffM140/((coeffM140*k140)/kFallOff140+1.0);
	CFO147 =	   coeffM147/((coeffM147*k147)/kFallOff147+1.0);
	CFO158 =	   coeffM158/((coeffM158*k158)/kFallOff158+1.0);
	CFO174 =	   coeffM174/((coeffM174*k174)/kFallOff174+1.0);


	lnPr50 =	   log((coeffM50*k50)/kFallOff50);
	lnPr52 =	   log((coeffM52*k52)/kFallOff52);
	lnPr54 =	   log((coeffM54*k54)/kFallOff54);
	lnPr56 =	   log((coeffM56*k56)/kFallOff56);
	lnPr57 =	   log((coeffM57*k57)/kFallOff57);
	lnPr59 =	   log((coeffM59*k59)/kFallOff59);
	lnPr63 =	   log((coeffM63*k63)/kFallOff63);
	lnPr70 =	   log((coeffM70*k70)/kFallOff70);
	lnPr71 =	   log((coeffM71*k71)/kFallOff71);
	lnPr72 =	   log((coeffM72*k72)/kFallOff72);
	lnPr74 =	   log((coeffM74*k74)/kFallOff74);
	lnPr76 =	   log((coeffM76*k76)/kFallOff76);
	lnPr83 =	   log((coeffM83*k83)/kFallOff83);
	lnPr85 =	   log((coeffM85*k85)/kFallOff85);
	lnPr95 =	   log((coeffM95*k95)/kFallOff95);
	lnPr131 =	   log((coeffM131*k131)/kFallOff131);
	lnPr140 =	   log((coeffM140*k140)/kFallOff140);
	lnPr147 =	   log((coeffM147*k147)/kFallOff147);
	lnPr158 =	   log((coeffM158*k158)/kFallOff158);
	lnPr174 =	   log((coeffM174*k174)/kFallOff174);


	wF50 =	pow(10.0, logFcent50/(1.0+BzzPow2((0.434294481903*lnPr50-0.40-0.67*logFcent50)/(0.806-0.0608012274665*lnPr50-1.1762*logFcent50)))); 
	wF52 =	pow(10.0, logFcent52/(1.0+BzzPow2((0.434294481903*lnPr52-0.40-0.67*logFcent52)/(0.806-0.0608012274665*lnPr52-1.1762*logFcent52)))); 
	wF54 =	pow(10.0, logFcent54/(1.0+BzzPow2((0.434294481903*lnPr54-0.40-0.67*logFcent54)/(0.806-0.0608012274665*lnPr54-1.1762*logFcent54)))); 
	wF56 =	pow(10.0, logFcent56/(1.0+BzzPow2((0.434294481903*lnPr56-0.40-0.67*logFcent56)/(0.806-0.0608012274665*lnPr56-1.1762*logFcent56)))); 
	wF57 =	pow(10.0, logFcent57/(1.0+BzzPow2((0.434294481903*lnPr57-0.40-0.67*logFcent57)/(0.806-0.0608012274665*lnPr57-1.1762*logFcent57)))); 
	wF59 =	pow(10.0, logFcent59/(1.0+BzzPow2((0.434294481903*lnPr59-0.40-0.67*logFcent59)/(0.806-0.0608012274665*lnPr59-1.1762*logFcent59)))); 
	wF63 =	pow(10.0, logFcent63/(1.0+BzzPow2((0.434294481903*lnPr63-0.40-0.67*logFcent63)/(0.806-0.0608012274665*lnPr63-1.1762*logFcent63)))); 
	wF70 =	pow(10.0, logFcent70/(1.0+BzzPow2((0.434294481903*lnPr70-0.40-0.67*logFcent70)/(0.806-0.0608012274665*lnPr70-1.1762*logFcent70)))); 
	wF71 =	pow(10.0, logFcent71/(1.0+BzzPow2((0.434294481903*lnPr71-0.40-0.67*logFcent71)/(0.806-0.0608012274665*lnPr71-1.1762*logFcent71)))); 
	wF72 =	pow(10.0, logFcent72/(1.0+BzzPow2((0.434294481903*lnPr72-0.40-0.67*logFcent72)/(0.806-0.0608012274665*lnPr72-1.1762*logFcent72)))); 
	wF74 =	pow(10.0, logFcent74/(1.0+BzzPow2((0.434294481903*lnPr74-0.40-0.67*logFcent74)/(0.806-0.0608012274665*lnPr74-1.1762*logFcent74)))); 
	wF76 =	pow(10.0, logFcent76/(1.0+BzzPow2((0.434294481903*lnPr76-0.40-0.67*logFcent76)/(0.806-0.0608012274665*lnPr76-1.1762*logFcent76)))); 
	wF83 =	pow(10.0, logFcent83/(1.0+BzzPow2((0.434294481903*lnPr83-0.40-0.67*logFcent83)/(0.806-0.0608012274665*lnPr83-1.1762*logFcent83)))); 
	wF85 =	pow(10.0, logFcent85/(1.0+BzzPow2((0.434294481903*lnPr85-0.40-0.67*logFcent85)/(0.806-0.0608012274665*lnPr85-1.1762*logFcent85)))); 
	wF95 =	pow(10.0, logFcent95/(1.0+BzzPow2((0.434294481903*lnPr95-0.40-0.67*logFcent95)/(0.806-0.0608012274665*lnPr95-1.1762*logFcent95)))); 
	wF131 =	pow(10.0, logFcent131/(1.0+BzzPow2((0.434294481903*lnPr131-0.40-0.67*logFcent131)/(0.806-0.0608012274665*lnPr131-1.1762*logFcent131)))); 
	wF140 =	pow(10.0, logFcent140/(1.0+BzzPow2((0.434294481903*lnPr140-0.40-0.67*logFcent140)/(0.806-0.0608012274665*lnPr140-1.1762*logFcent140)))); 
	wF147 =	pow(10.0, logFcent147/(1.0+BzzPow2((0.434294481903*lnPr147-0.40-0.67*logFcent147)/(0.806-0.0608012274665*lnPr147-1.1762*logFcent147)))); 
	wF158 =	pow(10.0, logFcent158/(1.0+BzzPow2((0.434294481903*lnPr158-0.40-0.67*logFcent158)/(0.806-0.0608012274665*lnPr158-1.1762*logFcent158)))); 
	wF174 =	pow(10.0, logFcent174/(1.0+BzzPow2((0.434294481903*lnPr174-0.40-0.67*logFcent174)/(0.806-0.0608012274665*lnPr174-1.1762*logFcent174)))); 


	coeffFallOff50 =	   CFO50*wF50;
	coeffFallOff52 =	   CFO52*wF52;
	coeffFallOff54 =	   CFO54*wF54;
	coeffFallOff56 =	   CFO56*wF56;
	coeffFallOff57 =	   CFO57*wF57;
	coeffFallOff59 =	   CFO59*wF59;
	coeffFallOff63 =	   CFO63*wF63;
	coeffFallOff70 =	   CFO70*wF70;
	coeffFallOff71 =	   CFO71*wF71;
	coeffFallOff72 =	   CFO72*wF72;
	coeffFallOff74 =	   CFO74*wF74;
	coeffFallOff76 =	   CFO76*wF76;
	coeffFallOff83 =	   CFO83*wF83;
	coeffFallOff85 =	   CFO85*wF85;
	coeffFallOff95 =	   CFO95*wF95;
	coeffFallOff131 =	   CFO131*wF131;
	coeffFallOff140 =	   CFO140*wF140;
	coeffFallOff147 =	   CFO147*wF147;
	coeffFallOff158 =	   CFO158*wF158;
	coeffFallOff174 =	   CFO174*wF174;


	// ============================================================ 
	// ===== REACTION RATES FOR THIRD BODY REACTIONS (Flat) ======= 
	// ============================================================ 
	rFlat1 =	   -k1*(c[4]*uK1-c[3]*c[3]);
	rFlat2 =	   k2*(c[2]*c[3]-c[5]*uK2);
	rFlat12 =	   k12*(c[15]*c[3]-c[16]*uK12);
	rFlat33 =	   k33*(c[2]*c[4]-c[7]*uK33);
	rFlat39 =	   -k39*(c[1]*uK39-c[2]*c[2]);
	rFlat43 =	   k43*(c[2]*c[5]-c[6]*uK43);
	rFlat50 =	   k50*(c[11]*c[2]-c[13]*uK50);
	rFlat52 =	   k52*(c[13]*c[2]-c[14]*uK52);
	rFlat54 =	   k54*(c[17]*c[2]-c[18]*uK54);
	rFlat56 =	   k56*(c[18]*c[2]-c[19]*uK56);
	rFlat57 =	   k57*(c[18]*c[2]-c[20]*uK57);
	rFlat59 =	   k59*(c[19]*c[2]-c[21]*uK59);
	rFlat63 =	   k63*(c[2]*c[20]-c[21]*uK63);
	rFlat70 =	   k70*(c[2]*c[22]-c[23]*uK70);
	rFlat71 =	   k71*(c[2]*c[23]-c[24]*uK71);
	rFlat72 =	   k72*(c[2]*c[24]-c[25]*uK72);
	rFlat74 =	   k74*(c[2]*c[25]-c[26]*uK74);
	rFlat76 =	   k76*(c[2]*c[26]-c[27]*uK76);
	rFlat83 =	   k83*(c[1]*c[15]-c[18]*uK83);
	rFlat85 =	   -k85*(c[8]*uK85-c[5]*c[5]);
	rFlat95 =	   k95*(c[13]*c[5]-c[21]*uK95);
	rFlat131 =	   k131*(c[10]*c[15]-c[28]*uK131);
	rFlat140 =	   k140*(c[11]*c[15]-c[29]*uK140);
	rFlat147 =	   k147*(c[12]*c[6]-c[21]*uK147);
	rFlat158 =	   -k158*(c[27]*uK158-c[13]*c[13]);
	rFlat167 =	   k167*(c[17]-c[15]*c[2]*uK167);
	rFlat174 =	   k174*(c[25]-c[1]*c[23]*uK174);


	// ============================================================ 
	// ===== REACTION RATES FOR EVERY REACTION ==================== 
	// ============================================================ 
	r1 =	   coeffM1*rFlat1;
	r2 =	   coeffM2*rFlat2;
	r3 =	   k3*(c[1]*c[3]-c[2]*c[5]*uK3);
	r4 =	   k4*(c[3]*c[7]-c[4]*c[5]*uK4);
	r5 =	   k5*(c[3]*c[8]-c[5]*c[7]*uK5);
	r6 =	   k6*(c[10]*c[3]-c[15]*c[2]*uK6);
	r7 =	   k7*(c[11]*c[3]-c[17]*c[2]*uK7);
	r8 =	   k8*(c[12]*c[3]-c[1]*c[15]*uK8);
	r9 =	   k9*(c[12]*c[3]-c[17]*c[2]*uK9);
	r10 =	   k10*(c[13]*c[3]-c[18]*c[2]*uK10);
	r11 =	   k11*(c[14]*c[3]-c[13]*c[5]*uK11);
	r12 =	   coeffM12*rFlat12;
	r13 =	   k13*(c[17]*c[3]-c[15]*c[5]*uK13);
	r14 =	   k14*(c[17]*c[3]-c[16]*c[2]*uK14);
	r15 =	   k15*(c[18]*c[3]-c[17]*c[5]*uK15);
	r16 =	   k16*(c[19]*c[3]-c[18]*c[5]*uK16);
	r17 =	   k17*(c[20]*c[3]-c[18]*c[5]*uK17);
	r18 =	   k18*(c[21]*c[3]-c[19]*c[5]*uK18);
	r19 =	   k19*(c[21]*c[3]-c[20]*c[5]*uK19);
	r20 =	   k20*(c[22]*c[3]-c[10]*c[15]*uK20);
	r21 =	   k21*(c[23]*c[3]-c[2]*c[28]*uK21);
	r22 =	   k22*(c[23]*c[3]-c[22]*c[5]*uK22);
	r23 =	   k23*(c[23]*c[3]-c[11]*c[15]*uK23);
	r24 =	   k24*(c[24]*c[3]-c[2]*c[29]*uK24);
	r25 =	   k25*(c[25]*c[3]-c[13]*c[17]*uK25);
	r26 =	   k26*(c[26]*c[3]-c[13]*c[18]*uK26);
	r27 =	   k27*(c[27]*c[3]-c[26]*c[5]*uK27);
	r28 =	   k28*(c[28]*c[3]-(c[15]*c[15])*c[2]*uK28);
	r29 =	   k29*(c[29]*c[3]-c[28]*c[5]*uK29);
	r30 =	   k30*(c[29]*c[3]-c[11]*c[16]*uK30);
	r31 =	   k31*(c[15]*c[4]-c[16]*c[3]*uK31);
	r32 =	   k32*(c[18]*c[4]-c[17]*c[7]*uK32);
	r33 =	   coeffM33*rFlat33;
	r34 =	   k34*(c[2]*(c[4]*c[4])-c[4]*c[7]*uK34);
	r35 =	   k35*(c[2]*c[4]*c[6]-c[6]*c[7]*uK35);
	r36 =	   k36*(c[2]*c[31]*c[4]-c[31]*c[7]*uK36);
	r37 =	   k37*(c[2]*c[32]*c[4]-c[32]*c[7]*uK37);
	r38 =	   k38*(c[2]*c[4]-c[3]*c[5]*uK38);
	r39 =	   coeffM39*rFlat39;
	r40 =	   k40*(c[1]*(c[2]*c[2])-(c[1]*c[1])*uK40);
	r41 =	   k41*((c[2]*c[2])*c[6]-c[1]*c[6]*uK41);
	r42 =	   k42*(c[16]*(c[2]*c[2])-c[1]*c[16]*uK42);
	r43 =	   coeffM43*rFlat43;
	r44 =	   k44*(c[2]*c[7]-c[3]*c[6]*uK44);
	r45 =	   k45*(c[2]*c[7]-c[1]*c[4]*uK45);
	r46 =	   k46*(c[2]*c[7]-(c[5]*c[5])*uK46);
	r47 =	   k47*(c[2]*c[8]-c[1]*c[7]*uK47);
	r48 =	   k48*(c[2]*c[8]-c[5]*c[6]*uK48);
	r49 =	   k49*(c[10]*c[2]-c[1]*c[9]*uK49);
	r50 =	   coeffFallOff50*rFlat50;
	r51 =	   k51*(c[12]*c[2]-c[1]*c[10]*uK51);
	r52 =	   coeffFallOff52*rFlat52;
	r53 =	   k53*(c[14]*c[2]-c[1]*c[13]*uK53);
	r54 =	   coeffFallOff54*rFlat54;
	r55 =	   k55*(c[17]*c[2]-c[1]*c[15]*uK55);
	r56 =	   coeffFallOff56*rFlat56;
	r57 =	   coeffFallOff57*rFlat57;
	r58 =	   k58*(c[18]*c[2]-c[1]*c[17]*uK58);
	r59 =	   coeffFallOff59*rFlat59;
	r60 =	   k60*(c[19]*c[2]-c[1]*c[18]*uK60);
	r61 =	   k61*(c[19]*c[2]-c[13]*c[5]*uK61);
	r62 =	   k62*(c[19]*c[2]-c[12]*c[6]*uK62);
	r63 =	   coeffFallOff63*rFlat63;
	r64 =	   k64*(c[2]*c[20]-c[19]*c[2]*uK64);
	r65 =	   k65*(c[2]*c[20]-c[1]*c[18]*uK65);
	r66 =	   k66*(c[2]*c[20]-c[13]*c[5]*uK66);
	r67 =	   k67*(c[2]*c[20]-c[12]*c[6]*uK67);
	r68 =	   k68*(c[2]*c[21]-c[1]*c[19]*uK68);
	r69 =	   k69*(c[2]*c[21]-c[1]*c[20]*uK69);
	r70 =	   coeffFallOff70*rFlat70;
	r71 =	   coeffFallOff71*rFlat71;
	r72 =	   coeffFallOff72*rFlat72;
	r73 =	   k73*(c[2]*c[24]-c[1]*c[23]*uK73);
	r74 =	   coeffFallOff74*rFlat74;
	r75 =	   k75*(c[2]*c[25]-c[1]*c[24]*uK75);
	r76 =	   coeffFallOff76*rFlat76;
	r77 =	   k77*(c[2]*c[26]-c[1]*c[25]*uK77);
	r78 =	   k78*(c[2]*c[27]-c[1]*c[26]*uK78);
	r79 =	   k79*(c[2]*c[28]-c[12]*c[15]*uK79);
	r80 =	   k80*(c[2]*c[29]-c[1]*c[28]*uK80);
	r81 =	   k81*(c[2]*c[29]-c[13]*c[15]*uK81);
	r82 =	   k82*(c[2]*c[30]-c[2]*c[29]*uK82);
	r83 =	   coeffFallOff83*rFlat83;
	r84 =	   k84*(c[1]*c[5]-c[2]*c[6]*uK84);
	r85 =	   coeffFallOff85*rFlat85;
	r86 =	   k86*(c[5]*c[5]-c[3]*c[6]*uK86);
	r87 =	   k87*(c[5]*c[7]-c[4]*c[6]*uK87);
	r88 =	   k88*(c[5]*c[8]-c[6]*c[7]*uK88);
	r89 =	   k89*(c[5]*c[8]-c[6]*c[7]*uK89);
	r90 =	   k90*(c[5]*c[9]-c[15]*c[2]*uK90);
	r91 =	   k91*(c[10]*c[5]-c[17]*c[2]*uK91);
	r92 =	   k92*(c[11]*c[5]-c[18]*c[2]*uK92);
	r93 =	   k93*(c[11]*c[5]-c[10]*c[6]*uK93);
	r94 =	   k94*(c[12]*c[5]-c[18]*c[2]*uK94);
	r95 =	   coeffFallOff95*rFlat95;
	r96 =	   k96*(c[13]*c[5]-c[11]*c[6]*uK96);
	r97 =	   k97*(c[13]*c[5]-c[12]*c[6]*uK97);
	r98 =	   k98*(c[14]*c[5]-c[13]*c[6]*uK98);
	r99 =	   k99*(c[15]*c[5]-c[16]*c[2]*uK99);
	r100 =	   k100*(c[17]*c[5]-c[15]*c[6]*uK100);
	r101 =	   k101*(c[18]*c[5]-c[17]*c[6]*uK101);
	r102 =	   k102*(c[19]*c[5]-c[18]*c[6]*uK102);
	r103 =	   k103*(c[20]*c[5]-c[18]*c[6]*uK103);
	r104 =	   k104*(c[21]*c[5]-c[19]*c[6]*uK104);
	r105 =	   k105*(c[21]*c[5]-c[20]*c[6]*uK105);
	r106 =	   k106*(c[22]*c[5]-c[2]*c[28]*uK106);
	r107 =	   k107*(c[23]*c[5]-c[2]*c[29]*uK107);
	r108 =	   k108*(c[23]*c[5]-c[2]*c[30]*uK108);
	r109 =	   k109*(c[23]*c[5]-c[22]*c[6]*uK109);
	r110 =	   k110*(c[23]*c[5]-c[13]*c[15]*uK110);
	r111 =	   k111*(c[24]*c[5]-c[23]*c[6]*uK111);
	r112 =	   k112*(c[25]*c[5]-c[24]*c[6]*uK112);
	r113 =	   k113*(c[27]*c[5]-c[26]*c[6]*uK113);
	r114 =	   k114*(c[29]*c[5]-c[28]*c[6]*uK114);
	r115 =	   k115*(c[7]*c[7]-c[4]*c[8]*uK115);
	r116 =	   k116*(c[7]*c[7]-c[4]*c[8]*uK116);
	r117 =	   k117*(c[11]*c[7]-c[18]*c[5]*uK117);
	r118 =	   k118*(c[13]*c[7]-c[14]*c[4]*uK118);
	r119 =	   k119*(c[13]*c[7]-c[20]*c[5]*uK119);
	r120 =	   k120*(c[15]*c[7]-c[16]*c[5]*uK120);
	r121 =	   k121*(c[18]*c[7]-c[17]*c[8]*uK121);
	r122 =	   k122*(c[4]*c[9]-c[15]*c[3]*uK122);
	r123 =	   k123*(c[11]*c[9]-c[2]*c[22]*uK123);
	r124 =	   k124*(c[13]*c[9]-c[2]*c[23]*uK124);
	r125 =	   k125*(c[10]*c[4]-c[17]*c[3]*uK125);
	r126 =	   k126*(c[1]*c[10]-c[11]*c[2]*uK126);
	r127 =	   k127*(c[10]*c[6]-c[18]*c[2]*uK127);
	r128 =	   k128*(c[10]*c[11]-c[2]*c[23]*uK128);
	r129 =	   k129*(c[10]*c[13]-c[2]*c[24]*uK129);
	r130 =	   k130*(c[10]*c[14]-c[2]*c[25]*uK130);
	r131 =	   coeffFallOff131*rFlat131;
	r132 =	   k132*(c[10]*c[16]-c[15]*c[17]*uK132);
	r133 =	   k133*(c[10]*c[18]-c[2]*c[29]*uK133);
	r134 =	   k134*(c[10]*c[28]-c[15]*c[23]*uK134);
	r135 =	   k135*(c[11]*c[4]-c[17]*c[5]*uK135);
	r136 =	   k136*(c[1]*c[11]-c[13]*c[2]*uK136);
	r137 =	   k137*(c[11]*c[11]-c[1]*c[23]*uK137);
	r138 =	   k138*(c[11]*c[13]-c[2]*c[25]*uK138);
	r139 =	   k139*(c[11]*c[14]-(c[13]*c[13])*uK139);
	r140 =	   coeffFallOff140*rFlat140;
	r141 =	   k141*(c[11]*c[28]-c[15]*c[24]*uK141);
	r142 =	   k142*(c[12]*c[31]-c[11]*c[31]*uK142);
	r143 =	   k143*(c[12]*c[32]-c[11]*c[32]*uK143);
	r144 =	   k144*(c[12]*c[4]-c[15]*c[2]*c[5]*uK144);
	r145 =	   k145*(c[12]*c[4]-c[15]*c[6]*uK145);
	r146 =	   k146*(c[1]*c[12]-c[13]*c[2]*uK146);
	r147 =	   coeffFallOff147*rFlat147;
	r148 =	   k148*(c[12]*c[6]-c[11]*c[6]*uK148);
	r149 =	   k149*(c[12]*c[13]-c[2]*c[25]*uK149);
	r150 =	   k150*(c[12]*c[14]-(c[13]*c[13])*uK150);
	r151 =	   k151*(c[12]*c[15]-c[11]*c[15]*uK151);
	r152 =	   k152*(c[12]*c[16]-c[11]*c[16]*uK152);
	r153 =	   k153*(c[12]*c[16]-c[15]*c[18]*uK153);
	r154 =	   k154*(c[12]*c[27]-c[13]*c[26]*uK154);
	r155 =	   k155*(c[13]*c[4]-c[20]*c[3]*uK155);
	r156 =	   k156*(c[13]*c[4]-c[18]*c[5]*uK156);
	r157 =	   k157*(c[13]*c[8]-c[14]*c[7]*uK157);
	r158 =	   coeffFallOff158*rFlat158;
	r159 =	   k159*(c[13]*c[13]-c[2]*c[26]*uK159);
	r160 =	   k160*(c[13]*c[17]-c[14]*c[15]*uK160);
	r161 =	   k161*(c[13]*c[18]-c[14]*c[17]*uK161);
	r162 =	   k162*(c[13]*c[21]-c[14]*c[19]*uK162);
	r163 =	   k163*(c[13]*c[21]-c[14]*c[20]*uK163);
	r164 =	   k164*(c[13]*c[25]-c[14]*c[24]*uK164);
	r165 =	   k165*(c[13]*c[27]-c[14]*c[26]*uK165);
	r166 =	   k166*(c[17]*c[6]-c[15]*c[2]*c[6]*uK166);
	r167 =	   coeffM167*rFlat167;
	r168 =	   k168*(c[17]*c[4]-c[15]*c[7]*uK168);
	r169 =	   k169*(c[19]*c[4]-c[18]*c[7]*uK169);
	r170 =	   k170*(c[20]*c[4]-c[18]*c[7]*uK170);
	r171 =	   k171*(c[22]*c[4]-c[15]*c[17]*uK171);
	r172 =	   k172*(c[1]*c[22]-c[2]*c[23]*uK172);
	r173 =	   k173*(c[24]*c[4]-c[17]*c[18]*uK173);
	r174 =	   coeffFallOff174*rFlat174;
	r175 =	   k175*(c[26]*c[4]-c[25]*c[7]*uK175);
	r176 =	   k176*(c[28]*c[4]-(c[15]*c[15])*c[5]*uK176);
	r177 =	   k177*(c[28]*c[28]-(c[15]*c[15])*c[23]*uK177);


	// ============================================================ 
	// ===== REACTION RATES FOR EVERY SPECIES ===================== 
	// ============================================================ 
	R[1] =	   -r126-r136+r137-r146-r172+r174-r3+r39+r40+r41+r42+r45+r47+r49+r51+r53+r55+r58+r60+r65+r68+r69+r73+r75+r77+r78+r8+r80-r83-r84;
	R[2] =	   r10+r106+r107+r108+r123+r124+r126+r127+r128+r129+r130+r133+r136+r138+r14+r144+r146+r149+r159+r166+r167+r172-r2+r21+r24+r28+r3-r33-r34-r35-r36-r37-r38-r39*2.0-r40*2.0-r41*2.0-r42*2.0-r43-r44-r45-r46-r47-r48-r49-r50-r51-r52-r53-r54-r55-r56-r57-r58-r59+r6-r60-r61-r62-r63-r65-r66-r67-r68-r69+r7-r70-r71-r72-r73-r74-r75-r76-r77-r78-r79-r80-r81+r84+r9+r90+r91+r92+r94+r99;
	R[3] =	   r1*-2.0-r10-r11-r12+r122+r125-r13-r14-r15+r155-r16-r17-r18-r19-r2-r20-r21-r22-r23-r24-r25-r26-r27-r28-r29-r3-r30+r31+r38-r4+r44-r5-r6-r7-r8+r86-r9;
	R[4] =	   r1+r115+r116+r118-r122-r125-r135-r144-r145-r155-r156-r168-r169-r170-r171-r173-r175-r176-r31-r32-r33-r34-r35-r36-r37-r38+r4+r45+r87;
	R[5] =	   -r100-r101-r102-r103-r104-r105-r106-r107-r108-r109+r11-r110-r111-r112-r113-r114+r117+r119+r120+r13+r135+r144+r15+r156+r16+r17+r176+r18+r19+r2+r22+r27+r29+r3+r38+r4-r43+r46*2.0+r48+r5+r61+r66-r84-r85*2.0-r86*2.0-r87-r88-r89-r90-r91-r92-r93-r94-r95-r96-r97-r98-r99;
	R[6] =	   r100+r101+r102+r103+r104+r105+r109+r111+r112+r113+r114-r127+r145-r147+r43+r44+r48+r62+r67+r84+r86+r87+r88+r89+r93+r96+r97+r98;
	R[7] =	   r115*-2.0-r116*2.0-r117-r118-r119-r120-r121+r157+r168+r169+r170+r175+r32+r33+r34+r35+r36+r37-r4-r44-r45-r46+r47+r5-r87+r88+r89;
	R[8] =	   r115+r116+r121-r157-r47-r48-r5+r85-r88-r89;
	R[9] =	   -r122-r123-r124+r49-r90;
	R[10] =	   -r125-r126-r127-r128-r129-r130-r131-r132-r133-r134+r20-r49+r51-r6-r91+r93;
	R[11] =	   -r117-r123+r126-r128-r135-r136-r137*2.0-r138-r139-r140-r141+r142+r143+r148+r151+r152+r23+r30-r50-r7-r92-r93+r96;
	R[12] =	   -r142-r143-r144-r145-r146-r147-r148-r149-r150-r151-r152-r153-r154-r51+r62+r67+r79-r8-r9-r94+r97;
	R[13] =	   -r10+r11+r110-r118-r119-r124-r129+r136-r138+r139*2.0+r146-r149+r150*2.0+r154-r155-r156-r157-r158*2.0-r159*2.0-r160-r161-r162-r163-r164-r165+r25+r26+r50-r52+r53+r61+r66+r81-r95-r96-r97+r98;
	R[14] =	   -r11+r118-r130-r139-r150+r157+r160+r161+r162+r163+r164+r165+r52-r53-r98;
	R[15] =	   r100+r110-r12-r120+r122+r13-r131+r132+r134-r140+r141+r144+r145+r153+r160+r166+r167+r168+r171+r176*2.0+r177*2.0+r20+r23+r28*2.0-r31+r55+r6+r79+r8+r81-r83+r90-r99;
	R[16] =	   r12+r120-r132+r14-r153+r30+r31+r99;
	R[17] =	   -r100+r101+r121+r125-r13+r132+r135-r14+r15-r160+r161-r166-r167-r168+r171+r173+r25+r32-r54-r55+r58+r7+r9+r91;
	R[18] =	   r10-r101+r102+r103+r117-r121+r127-r133-r15+r153+r156+r16-r161+r169+r17+r170+r173+r26-r32+r54-r56-r57-r58+r60+r65+r83+r92+r94;
	R[19] =	   -r102+r104-r16+r162-r169+r18+r56-r59-r60-r61-r62+r64+r68;
	R[20] =	   -r103+r105+r119+r155+r163-r17-r170+r19+r57-r63-r64-r65-r66-r67+r69;
	R[21] =	   -r104-r105+r147-r162-r163-r18-r19+r59+r63-r68-r69+r95;
	R[22] =	   -r106+r109+r123-r171-r172-r20+r22-r70;
	R[23] =	   -r107-r108-r109-r110+r111+r124+r128+r134+r137+r172+r174+r177-r21-r22-r23+r70-r71+r73;
	R[24] =	   -r111+r112+r129+r141+r164-r173-r24+r71-r72-r73+r75;
	R[25] =	   -r112+r130+r138+r149-r164-r174+r175-r25+r72-r74-r75+r77;
	R[26] =	   r113+r154+r159+r165-r175-r26+r27+r74-r76-r77+r78;
	R[27] =	   -r113-r154+r158-r165-r27+r76-r78;
	R[28] =	   r106+r114+r131-r134-r141-r176-r177*2.0+r21-r28+r29-r79+r80;
	R[29] =	   r107-r114+r133+r140+r24-r29-r30-r80-r81+r82;
	R[30] =	   r108-r82;
	R[31] =	   0.0;
	R[32] =	   0.0;

}