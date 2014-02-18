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

#include "symbolickinetics/fluent_drm22_polimi_nox/OpenSMOKE_SymbolicKinetics_Fluent_DRM22_Polimi_NOX.h"

void OpenSMOKE_SymbolicKinetics_Fluent_DRM22_Polimi_NOX::giveReactionRates(double cTot, BzzVector &c, BzzVector &R) 
{
		// ============================================================ 
	// ===== CORRECTION COEFFICIENTS FOR THIRD BODY REACTIONS ===== 
	// ============================================================ 
	coeffM1 =	   c[1]*1.0+c[11]*1.0+c[19]*2.0+c[2]*5.0E-1+c[3]*1.0+c[4]*5.0+cTot;
	coeffM9 =	   c[1]*1.0+c[11]*1.0+c[19]*2.0+c[2]*5.0E-1+c[3]*2.5+c[4]*5.0+c[5]*5.0+cTot;
	coeffM20 =	   c[19]*5.0E-1-c[2]*2.5E-1-c[24]*1.0+c[3]*5.0E-1-c[4]*1.0-c[5]*1.0+cTot;
	coeffM25 =	   c[1]*1.0-c[11]*1.0+c[19]*2.0-c[3]*1.0-c[4]*1.0+cTot;
	coeffM29 =	   c[1]*1.0-c[11]*2.7E-1+c[19]*2.0+c[4]*2.65+cTot;
	coeffM33 =	   c[1]*1.0+c[11]*1.0+c[19]*2.0+c[2]*5.0E-1+c[3]*1.0+c[4]*5.0+cTot;
	coeffM34 =	   c[1]*1.0+c[11]*1.0+c[19]*2.0+c[2]*5.0E-1+c[3]*1.0+c[4]*5.0+cTot;
	coeffM36 =	   c[1]*1.0+c[11]*1.0+c[19]*2.0+c[2]*5.0E-1+c[3]*1.0+c[4]*5.0+cTot;
	coeffM38 =	   c[1]*1.0+c[11]*1.0+c[19]*2.0+c[2]*5.0E-1+c[3]*1.0+c[4]*5.0+cTot;
	coeffM41 =	   c[1]*1.0+c[11]*1.0+c[19]*2.0+c[2]*5.0E-1+c[3]*1.0+c[4]*5.0+cTot;
	coeffM42 =	   c[1]*1.0+c[11]*1.0+c[19]*2.0+c[2]*5.0E-1+c[3]*1.0+c[4]*5.0+cTot;
	coeffM44 =	   c[1]*1.0+c[11]*1.0+c[19]*2.0+c[2]*5.0E-1+c[3]*1.0+c[4]*5.0+cTot;
	coeffM46 =	   c[1]*1.0+c[11]*1.0+c[19]*2.0+c[2]*5.0E-1+c[3]*1.0+c[4]*5.0+cTot;
	coeffM48 =	   c[1]*1.0+c[11]*1.0+c[19]*2.0+c[2]*5.0E-1+c[3]*1.0+c[4]*5.0+cTot;
	coeffM50 =	   c[1]*1.0+c[11]*1.0+c[19]*2.0+c[2]*5.0E-1+c[3]*1.0+c[4]*5.0+cTot;
	coeffM90 =	   c[1]*1.0+c[11]*1.0+c[19]*2.0+c[2]*5.0E-1+c[3]*1.0+c[4]*5.0+cTot;
	coeffM97 =	   c[1]*1.0+c[11]*1.0+c[19]*2.0+c[2]*5.0E-1+c[3]*1.0-c[4]*1.0+cTot;
	coeffM101 =	   c[1]*1.0+c[11]*1.0+c[19]*2.0+c[2]*5.0E-1+c[3]*1.0+c[4]*5.0+cTot;
	coeffM117 =	   c[24]*5.0E-1+c[4]*1.76E1+c[5]*5.0E-1+cTot;
	coeffM121 =	   c[24]*5.0E-1+c[4]*1.76E1+c[5]*5.0E-1+cTot;
	coeffM146 =	   c[11]*1.0+c[24]*1.0+c[4]*9.0+c[5]*1.0+cTot;
	coeffM152 =	   cTot;


	// ============================================================ 
	// ===== CORRECTION COEFFICIENTS FOR FALL OFF REACTIONS ======= 
	// ============================================================ 
	CFO33 =	   coeffM33/((coeffM33*k33)/kFallOff33+1.0);
	CFO34 =	   coeffM34/((coeffM34*k34)/kFallOff34+1.0);
	CFO36 =	   coeffM36/((coeffM36*k36)/kFallOff36+1.0);
	CFO38 =	   coeffM38/((coeffM38*k38)/kFallOff38+1.0);
	CFO41 =	   coeffM41/((coeffM41*k41)/kFallOff41+1.0);
	CFO42 =	   coeffM42/((coeffM42*k42)/kFallOff42+1.0);
	CFO44 =	   coeffM44/((coeffM44*k44)/kFallOff44+1.0);
	CFO46 =	   coeffM46/((coeffM46*k46)/kFallOff46+1.0);
	CFO48 =	   coeffM48/((coeffM48*k48)/kFallOff48+1.0);
	CFO50 =	   coeffM50/((coeffM50*k50)/kFallOff50+1.0);
	CFO90 =	   coeffM90/((coeffM90*k90)/kFallOff90+1.0);
	CFO101 =	   coeffM101/((coeffM101*k101)/kFallOff101+1.0);
	CFO152 =	   coeffM152/((coeffM152*k152)/kFallOff152+1.0);


	lnPr33 =	   log((coeffM33*k33)/kFallOff33);
	lnPr34 =	   log((coeffM34*k34)/kFallOff34);
	lnPr36 =	   log((coeffM36*k36)/kFallOff36);
	lnPr38 =	   log((coeffM38*k38)/kFallOff38);
	lnPr41 =	   log((coeffM41*k41)/kFallOff41);
	lnPr42 =	   log((coeffM42*k42)/kFallOff42);
	lnPr44 =	   log((coeffM44*k44)/kFallOff44);
	lnPr46 =	   log((coeffM46*k46)/kFallOff46);
	lnPr48 =	   log((coeffM48*k48)/kFallOff48);
	lnPr50 =	   log((coeffM50*k50)/kFallOff50);
	lnPr90 =	   log((coeffM90*k90)/kFallOff90);
	lnPr101 =	   log((coeffM101*k101)/kFallOff101);


	wF33 =	pow(10.0, logFcent33/(1.0+BzzPow2((0.434294481903*lnPr33-0.40-0.67*logFcent33)/(0.806-0.0608012274665*lnPr33-1.1762*logFcent33)))); 
	wF34 =	pow(10.0, logFcent34/(1.0+BzzPow2((0.434294481903*lnPr34-0.40-0.67*logFcent34)/(0.806-0.0608012274665*lnPr34-1.1762*logFcent34)))); 
	wF36 =	pow(10.0, logFcent36/(1.0+BzzPow2((0.434294481903*lnPr36-0.40-0.67*logFcent36)/(0.806-0.0608012274665*lnPr36-1.1762*logFcent36)))); 
	wF38 =	pow(10.0, logFcent38/(1.0+BzzPow2((0.434294481903*lnPr38-0.40-0.67*logFcent38)/(0.806-0.0608012274665*lnPr38-1.1762*logFcent38)))); 
	wF41 =	pow(10.0, logFcent41/(1.0+BzzPow2((0.434294481903*lnPr41-0.40-0.67*logFcent41)/(0.806-0.0608012274665*lnPr41-1.1762*logFcent41)))); 
	wF42 =	pow(10.0, logFcent42/(1.0+BzzPow2((0.434294481903*lnPr42-0.40-0.67*logFcent42)/(0.806-0.0608012274665*lnPr42-1.1762*logFcent42)))); 
	wF44 =	pow(10.0, logFcent44/(1.0+BzzPow2((0.434294481903*lnPr44-0.40-0.67*logFcent44)/(0.806-0.0608012274665*lnPr44-1.1762*logFcent44)))); 
	wF46 =	pow(10.0, logFcent46/(1.0+BzzPow2((0.434294481903*lnPr46-0.40-0.67*logFcent46)/(0.806-0.0608012274665*lnPr46-1.1762*logFcent46)))); 
	wF48 =	pow(10.0, logFcent48/(1.0+BzzPow2((0.434294481903*lnPr48-0.40-0.67*logFcent48)/(0.806-0.0608012274665*lnPr48-1.1762*logFcent48)))); 
	wF50 =	pow(10.0, logFcent50/(1.0+BzzPow2((0.434294481903*lnPr50-0.40-0.67*logFcent50)/(0.806-0.0608012274665*lnPr50-1.1762*logFcent50)))); 
	wF90 =	pow(10.0, logFcent90/(1.0+BzzPow2((0.434294481903*lnPr90-0.40-0.67*logFcent90)/(0.806-0.0608012274665*lnPr90-1.1762*logFcent90)))); 
	wF101 =	pow(10.0, logFcent101/(1.0+BzzPow2((0.434294481903*lnPr101-0.40-0.67*logFcent101)/(0.806-0.0608012274665*lnPr101-1.1762*logFcent101)))); 
	wF152 =	1.00; 


	coeffFallOff33 =	   CFO33*wF33;
	coeffFallOff34 =	   CFO34*wF34;
	coeffFallOff36 =	   CFO36*wF36;
	coeffFallOff38 =	   CFO38*wF38;
	coeffFallOff41 =	   CFO41*wF41;
	coeffFallOff42 =	   CFO42*wF42;
	coeffFallOff44 =	   CFO44*wF44;
	coeffFallOff46 =	   CFO46*wF46;
	coeffFallOff48 =	   CFO48*wF48;
	coeffFallOff50 =	   CFO50*wF50;
	coeffFallOff90 =	   CFO90*wF90;
	coeffFallOff101 =	   CFO101*wF101;
	coeffFallOff152 =	   CFO152*wF152;


	// ============================================================ 
	// ===== REACTION RATES FOR THIRD BODY REACTIONS (Flat) ======= 
	// ============================================================ 
	rFlat1 =	   k1*(c[12]*c[13]-c[6]*uK1);
	rFlat9 =	   k9*(c[13]*c[2]-c[3]*uK9);
	rFlat20 =	   k20*(c[12]*c[5]-c[7]*uK20);
	rFlat25 =	   -k25*(c[11]*uK25-c[12]*c[12]);
	rFlat29 =	   k29*(c[12]*c[6]-c[4]*uK29);
	rFlat33 =	   k33*(c[12]*c[8]-c[10]*uK33);
	rFlat34 =	   k34*(c[10]*c[12]-c[1]*uK34);
	rFlat36 =	   k36*(c[12]*c[14]-c[15]*uK36);
	rFlat38 =	   k38*(c[12]*c[15]-c[16]*uK38);
	rFlat41 =	   k41*(c[12]*c[21]-c[22]*uK41);
	rFlat42 =	   k42*(c[12]*c[22]-c[17]*uK42);
	rFlat44 =	   k44*(c[12]*c[17]-c[18]*uK44);
	rFlat46 =	   k46*(c[12]*c[18]-c[19]*uK46);
	rFlat48 =	   k48*(c[11]*c[2]-c[15]*uK48);
	rFlat50 =	   -k50*(c[20]*uK50-c[6]*c[6]);
	rFlat90 =	   -k90*(c[19]*uK90-c[10]*c[10]);
	rFlat97 =	   k97*(c[14]-c[12]*c[2]*uK97);
	rFlat101 =	   k101*(c[17]-c[11]*c[21]*uK101);
	rFlat117 =	   k117*(c[28]-c[2]*c[26]*uK117);
	rFlat121 =	   k121*(c[30]-c[2]*c[29]*uK121);
	rFlat146 =	   k146*(c[33]-c[12]*c[25]*uK146);
	rFlat152 =	   k152*(c[34]-c[13]*c[24]*uK152);


	// ============================================================ 
	// ===== REACTION RATES FOR EVERY REACTION ==================== 
	// ============================================================ 
	r1 =	   coeffM1*rFlat1;
	r2 =	   k2*(c[11]*c[13]-c[12]*c[6]*uK2);
	r3 =	   k3*(c[13]*c[7]-c[5]*c[6]*uK3);
	r4 =	   k4*(c[13]*c[8]-(c[12]*c[12])*c[2]*uK4);
	r5 =	   k5*(c[13]*c[8]-c[11]*c[2]*uK5);
	r6 =	   k6*(c[13]*c[9]-c[12]*c[14]*uK6);
	r7 =	   k7*(c[10]*c[13]-c[12]*c[15]*uK7);
	r8 =	   k8*(c[1]*c[13]-c[10]*c[6]*uK8);
	r9 =	   coeffM9*rFlat9;
	r10 =	   k10*(c[13]*c[14]-c[2]*c[6]*uK10);
	r11 =	   k11*(c[13]*c[14]-c[12]*c[3]*uK11);
	r12 =	   k12*(c[13]*c[15]-c[14]*c[6]*uK12);
	r13 =	   k13*(c[13]*c[21]-c[2]*c[9]*uK13);
	r14 =	   k14*(c[13]*c[21]-c[2]*c[8]*uK14);
	r15 =	   k15*(c[13]*c[17]-c[10]*c[14]*uK15);
	r16 =	   k16*(c[13]*c[18]-c[10]*c[15]*uK16);
	r17 =	   k17*(c[13]*c[19]-c[18]*c[6]*uK17);
	r18 =	   k18*(c[2]*c[5]-c[13]*c[3]*uK18);
	r19 =	   k19*(c[15]*c[5]-c[14]*c[7]*uK19);
	r20 =	   coeffM20*rFlat20;
	r21 =	   k21*(c[12]*(c[5]*c[5])-c[5]*c[7]*uK21);
	r22 =	   k22*(c[12]*c[4]*c[5]-c[4]*c[7]*uK22);
	r23 =	   k23*(c[12]*c[24]*c[5]-c[24]*c[7]*uK23);
	r24 =	   k24*(c[12]*c[5]-c[13]*c[6]*uK24);
	r25 =	   coeffM25*rFlat25;
	r26 =	   k26*(c[11]*(c[12]*c[12])-(c[11]*c[11])*uK26);
	r27 =	   k27*((c[12]*c[12])*c[4]-c[11]*c[4]*uK27);
	r28 =	   k28*((c[12]*c[12])*c[3]-c[11]*c[3]*uK28);
	r29 =	   coeffM29*rFlat29;
	r30 =	   k30*(c[12]*c[7]-c[11]*c[5]*uK30);
	r31 =	   k31*(c[12]*c[7]-(c[6]*c[6])*uK31);
	r32 =	   k32*(c[12]*c[20]-c[11]*c[7]*uK32);
	r33 =	   coeffFallOff33*rFlat33;
	r34 =	   coeffFallOff34*rFlat34;
	r35 =	   k35*(c[1]*c[12]-c[10]*c[11]*uK35);
	r36 =	   coeffFallOff36*rFlat36;
	r37 =	   k37*(c[12]*c[14]-c[11]*c[2]*uK37);
	r38 =	   coeffFallOff38*rFlat38;
	r39 =	   k39*(c[12]*c[15]-c[11]*c[14]*uK39);
	r40 =	   k40*(c[12]*c[16]-c[10]*c[6]*uK40);
	r41 =	   coeffFallOff41*rFlat41;
	r42 =	   coeffFallOff42*rFlat42;
	r43 =	   k43*(c[12]*c[22]-c[11]*c[21]*uK43);
	r44 =	   coeffFallOff44*rFlat44;
	r45 =	   k45*(c[12]*c[17]-c[11]*c[22]*uK45);
	r46 =	   coeffFallOff46*rFlat46;
	r47 =	   k47*(c[12]*c[19]-c[11]*c[18]*uK47);
	r48 =	   coeffFallOff48*rFlat48;
	r49 =	   k49*(c[11]*c[6]-c[12]*c[4]*uK49);
	r50 =	   coeffFallOff50*rFlat50;
	r51 =	   k51*(c[6]*c[6]-c[13]*c[4]*uK51);
	r52 =	   k52*(c[6]*c[7]-c[4]*c[5]*uK52);
	r53 =	   k53*(c[20]*c[6]-c[4]*c[7]*uK53);
	r54 =	   k54*(c[6]*c[8]-c[12]*c[15]*uK54);
	r55 =	   k55*(c[6]*c[9]-c[12]*c[15]*uK55);
	r56 =	   k56*(c[10]*c[6]-c[4]*c[9]*uK56);
	r57 =	   k57*(c[1]*c[6]-c[10]*c[4]*uK57);
	r58 =	   k58*(c[2]*c[6]-c[12]*c[3]*uK58);
	r59 =	   k59*(c[14]*c[6]-c[2]*c[4]*uK59);
	r60 =	   k60*(c[15]*c[6]-c[14]*c[4]*uK60);
	r61 =	   k61*(c[21]*c[6]-c[10]*c[2]*uK61);
	r62 =	   k62*(c[22]*c[6]-c[21]*c[4]*uK62);
	r63 =	   k63*(c[17]*c[6]-c[22]*c[4]*uK63);
	r64 =	   k64*(c[19]*c[6]-c[18]*c[4]*uK64);
	r65 =	   k65*(c[7]*c[7]-c[20]*c[5]*uK65);
	r66 =	   k66*(c[7]*c[7]-c[20]*c[5]*uK66);
	r67 =	   k67*(c[7]*c[8]-c[15]*c[6]*uK67);
	r68 =	   k68*(c[10]*c[7]-c[1]*c[5]*uK68);
	r69 =	   k69*(c[10]*c[7]-c[16]*c[6]*uK69);
	r70 =	   k70*(c[2]*c[7]-c[3]*c[6]*uK70);
	r71 =	   k71*(c[15]*c[7]-c[14]*c[20]*uK71);
	r72 =	   k72*(c[5]*c[8]-c[14]*c[6]*uK72);
	r73 =	   k73*(c[11]*c[8]-c[10]*c[12]*uK73);
	r74 =	   k74*(c[8]*c[8]-c[11]*c[21]*uK74);
	r75 =	   k75*(c[10]*c[8]-c[12]*c[17]*uK75);
	r76 =	   k76*(c[1]*c[8]-(c[10]*c[10])*uK76);
	r77 =	   k77*(c[24]*c[9]-c[24]*c[8]*uK77);
	r78 =	   k78*(c[5]*c[9]-c[12]*c[2]*c[6]*uK78);
	r79 =	   k79*(c[5]*c[9]-c[2]*c[4]*uK79);
	r80 =	   k80*(c[11]*c[9]-c[10]*c[12]*uK80);
	r81 =	   k81*(c[4]*c[9]-c[4]*c[8]*uK81);
	r82 =	   k82*(c[10]*c[9]-c[12]*c[17]*uK82);
	r83 =	   k83*(c[1]*c[9]-(c[10]*c[10])*uK83);
	r84 =	   k84*(c[2]*c[9]-c[2]*c[8]*uK84);
	r85 =	   k85*(c[3]*c[9]-c[3]*c[8]*uK85);
	r86 =	   k86*(c[3]*c[9]-c[15]*c[2]*uK86);
	r87 =	   k87*(c[10]*c[5]-c[13]*c[16]*uK87);
	r88 =	   k88*(c[10]*c[5]-c[15]*c[6]*uK88);
	r89 =	   k89*(c[10]*c[20]-c[1]*c[7]*uK89);
	r90 =	   coeffFallOff90*rFlat90;
	r91 =	   k91*(c[10]*c[10]-c[12]*c[18]*uK91);
	r92 =	   k92*(c[10]*c[14]-c[1]*c[2]*uK92);
	r93 =	   k93*(c[10]*c[15]-c[1]*c[14]*uK93);
	r94 =	   k94*(c[10]*c[17]-c[1]*c[22]*uK94);
	r95 =	   k95*(c[10]*c[19]-c[1]*c[18]*uK95);
	r96 =	   k96*(c[14]*c[4]-c[12]*c[2]*c[4]*uK96);
	r97 =	   coeffM97*rFlat97;
	r98 =	   k98*(c[14]*c[5]-c[2]*c[7]*uK98);
	r99 =	   k99*(c[16]*c[5]-c[15]*c[7]*uK99);
	r100 =	   k100*(c[22]*c[5]-c[14]*c[15]*uK100);
	r101 =	   coeffFallOff101*rFlat101;
	r102 =	   k102*(c[18]*c[5]-c[17]*c[7]*uK102);
	r103 =	   k103*(c[12]*c[8]-c[11]*c[23]*uK103);
	r104 =	   k104*(c[6]*c[8]-c[23]*c[4]*uK104);
	r105 =	   k105*(c[12]*c[9]-c[11]*c[23]*uK105);
	r106 =	   k106*(c[1]*c[23]-c[12]*c[17]*uK106);
	r107 =	   k107*(c[13]*c[23]-c[12]*c[2]*uK107);
	r108 =	   k108*(c[23]*c[6]-c[12]*c[14]*uK108);
	r109 =	   k109*(c[23]*c[5]-c[13]*c[14]*uK109);
	r110 =	   k110*(c[23]*c[4]-c[12]*c[15]*uK110);
	r111 =	   k111*(c[23]*c[3]-c[14]*c[2]*uK111);
	r112 =	   k112*(c[13]*c[24]-c[25]*c[26]*uK112);
	r113 =	   k113*(c[26]*c[5]-c[13]*c[25]*uK113);
	r114 =	   k114*(c[26]*c[6]-c[12]*c[25]*uK114);
	r115 =	   k115*(c[23]*c[24]-c[26]*c[27]*uK115);
	r116 =	   k116*(c[13]*c[27]-c[12]*c[28]*uK116);
	r117 =	   coeffM117*rFlat117;
	r118 =	   k118*(c[12]*c[28]-c[2]*c[29]*uK118);
	r119 =	   k119*(c[13]*c[28]-c[2]*c[25]*uK119);
	r120 =	   k120*(c[11]*c[28]-c[12]*c[30]*uK120);
	r121 =	   coeffM121*rFlat121;
	r122 =	   k122*(c[12]*c[30]-c[2]*c[31]*uK122);
	r123 =	   k123*(c[13]*c[30]-c[28]*c[6]*uK123);
	r124 =	   k124*(c[13]*c[30]-c[29]*c[3]*uK124);
	r125 =	   k125*(c[30]*c[6]-c[28]*c[4]*uK125);
	r126 =	   k126*(c[11]*c[32]-c[12]*c[27]*uK126);
	r127 =	   k127*(c[32]*c[4]-c[27]*c[6]*uK127);
	r128 =	   k128*(c[32]*c[6]-c[12]*c[28]*uK128);
	r129 =	   k129*(c[32]*c[5]-c[13]*c[28]*uK129);
	r130 =	   k130*(c[12]*c[29]-c[11]*c[26]*uK130);
	r131 =	   k131*(c[13]*c[29]-c[12]*c[25]*uK131);
	r132 =	   k132*(c[29]*c[6]-c[12]*c[33]*uK132);
	r133 =	   k133*(c[29]*c[6]-c[26]*c[4]*uK133);
	r134 =	   k134*(c[29]*c[5]-c[13]*c[33]*uK134);
	r135 =	   k135*(c[25]*c[29]-c[12]*c[34]*uK135);
	r136 =	   k136*(c[25]*c[29]-c[24]*c[6]*uK136);
	r137 =	   k137*(c[12]*c[31]-c[11]*c[29]*uK137);
	r138 =	   k138*(c[13]*c[31]-c[12]*c[33]*uK138);
	r139 =	   k139*(c[31]*c[6]-c[29]*c[4]*uK139);
	r140 =	   k140*(c[25]*c[31]-c[24]*c[4]*uK140);
	r141 =	   k141*(c[25]*c[31]-c[35]*c[6]*uK141);
	r142 =	   k142*(c[35]-c[12]*c[24]*uK142);
	r143 =	   k143*(c[12]*c[35]-c[11]*c[24]*uK143);
	r144 =	   k144*(c[13]*c[35]-c[12]*c[34]*uK144);
	r145 =	   k145*(c[35]*c[6]-c[24]*c[4]*uK145);
	r146 =	   coeffM146*rFlat146;
	r147 =	   k147*(c[12]*c[33]-c[11]*c[25]*uK147);
	r148 =	   k148*(c[33]*c[6]-c[25]*c[4]*uK148);
	r149 =	   c[10]*c[25]*k149;
	r150 =	   k150*(c[25]*c[8]-c[12]*c[30]*uK150);
	r151 =	   k151*(c[23]*c[25]-c[13]*c[27]*uK151);
	r152 =	   coeffFallOff152*rFlat152;
	r153 =	   k153*(c[12]*c[34]-c[24]*c[6]*uK153);
	r154 =	   k154*(c[13]*c[34]-(c[25]*c[25])*uK154);
	r155 =	   k155*(c[34]*c[6]-c[24]*c[7]*uK155);


	// ============================================================ 
	// ===== REACTION RATES FOR EVERY SPECIES ===================== 
	// ============================================================ 
	R[1] =	   -r106+r34-r35-r57+r68-r76-r8-r83+r89+r92+r93+r94+r95;
	R[2] =	   r10+r107+r111+r117+r118+r119+r121+r122+r13+r14-r18+r37+r4-r48+r5-r58+r59+r61-r70+r78+r79+r86-r9+r92+r96+r97+r98;
	R[3] =	   r11-r111+r124+r18+r58+r70-r86+r9;
	R[4] =	   r104-r110+r125-r127+r133+r139+r140+r145+r148+r149+r29+r49+r51+r52+r53+r56+r57+r59+r60+r62+r63+r64+r79;
	R[5] =	   -r100-r102-r109-r113-r129-r134-r18-r19-r20-r21-r22-r23-r24+r3+r30+r52+r65+r66+r68-r72-r78-r79-r87-r88-r98-r99;
	R[6] =	   r1+r10-r104-r108-r114+r12+r123-r125+r127-r128-r132-r133+r136-r139+r141-r145-r148+r153-r155+r17+r2+r24-r29+r3+r31*2.0+r40-r49-r50*2.0-r51*2.0-r52-r53-r54-r55-r56-r57-r58-r59-r60-r61-r62-r63-r64+r67+r69+r70+r72+r78+r8+r88;
	R[7] =	   r102+r155+r19+r20+r21+r22+r23-r3-r30-r31+r32-r52+r53-r65*2.0-r66*2.0-r67-r68-r69-r70-r71+r89+r98+r99;
	R[8] =	   -r103-r104+r14-r150-r33-r4-r5-r54-r67-r72-r73-r74*2.0-r75-r76+r77+r81+r84+r85;
	R[9] =	   -r105+r13-r55+r56-r6-r77-r78-r79-r80-r81-r82-r83-r84-r85-r86;
	R[10] =	   -r149+r15+r16+r33-r34+r35+r40-r56+r57+r61-r68-r69-r7+r73-r75+r76*2.0+r8+r80-r82+r83*2.0-r87-r88-r89-r90*2.0-r91*2.0-r92-r93-r94-r95;
	R[11] =	   r101+r103+r105-r120-r126+r130+r137+r143+r147-r2+r25+r26+r27+r28+r30+r32+r35+r37+r39+r43+r45+r47-r48-r49+r5-r73+r74-r80;
	R[12] =	   -r1-r103-r105+r106+r107+r108+r11+r110+r114+r116-r118+r120-r122+r126+r128-r130+r131+r132+r135-r137+r138+r142-r143+r144+r146-r147+r150-r153+r2-r20-r21-r22-r23-r24-r25*2.0-r26*2.0-r27*2.0-r28*2.0-r29-r30-r31-r32-r33-r34-r35-r36-r37-r38-r39+r4*2.0-r40-r41-r42-r43-r44-r45-r46-r47+r49+r54+r55+r58+r6+r7+r73+r75+r78+r80+r82+r91+r96+r97;
	R[13] =	   -r1-r10-r107+r109-r11-r112+r113-r116-r119-r12-r123-r124+r129-r13-r131+r134-r138-r14-r144-r15+r151+r152-r154-r16-r17+r18-r2+r24-r3-r4-r5+r51-r6-r7-r8+r87-r9;
	R[14] =	   -r10+r100+r108+r109-r11+r111+r12+r15+r19-r36-r37+r39-r59+r6+r60+r71+r72-r92+r93-r96-r97-r98;
	R[15] =	   r100+r110-r12+r16-r19+r36-r38-r39+r48+r54+r55-r60+r67+r7-r71+r86+r88-r93+r99;
	R[16] =	   r38-r40+r69+r87-r99;
	R[17] =	   -r101+r102+r106-r15+r42-r44-r45-r63+r75+r82-r94;
	R[18] =	   -r102-r16+r17+r44-r46+r47+r64+r91+r95;
	R[19] =	   -r17+r46-r47-r64+r90-r95;
	R[20] =	   -r32+r50-r53+r65+r66+r71-r89;
	R[21] =	   r101-r13-r14-r41+r43-r61+r62+r74;
	R[22] =	   -r100+r41-r42-r43+r45-r62+r63+r94;
	R[23] =	   r103+r104+r105-r106-r107-r108-r109-r110-r111-r115-r151;
	R[24] =	   -r112-r115+r136+r140+r142+r143+r145+r152+r153+r155;
	R[25] =	   r112+r113+r114+r119+r131-r135-r136-r140-r141+r146+r147+r148-r149-r150-r151+r154*2.0;
	R[26] =	   r112-r113-r114+r115+r117+r130+r133;
	R[27] =	   r115-r116+r126+r127+r149+r151;
	R[28] =	   r116-r117-r118-r119-r120+r123+r125+r128+r129;
	R[29] =	   r118+r121+r124-r130-r131-r132-r133-r134-r135-r136+r137+r139;
	R[30] =	   r120-r121-r122-r123-r124-r125+r150;
	R[31] =	   r122-r137-r138-r139-r140-r141;
	R[32] =	   -r126-r127-r128-r129;
	R[33] =	   r132+r134+r138-r146-r147-r148;
	R[34] =	   r135+r144-r152-r153-r154-r155;
	R[35] =	   r141-r142-r143-r144-r145;
}