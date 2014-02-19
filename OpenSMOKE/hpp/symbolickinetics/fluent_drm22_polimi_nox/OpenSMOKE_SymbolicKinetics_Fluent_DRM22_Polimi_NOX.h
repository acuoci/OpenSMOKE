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

#if !defined(OPENSMOKE_SYMBOLICKINETICS_FLUENT_DRM22_POLIMI_NOX_H)
#define OPENSMOKE_SYMBOLICKINETICS_FLUENT_DRM22_POLIMI_NOX_H

#include "BzzMath.hpp"
#include "symbolickinetics/OpenSMOKE_SymbolicKinetics.h"

class OpenSMOKE_SymbolicKinetics_Fluent_DRM22_Polimi_NOX : public OpenSMOKE_SymbolicKinetics
{
public:
	
	OpenSMOKE_SymbolicKinetics_Fluent_DRM22_Polimi_NOX();

	virtual void assignKineticConstants(BzzVector &xxxk, BzzVector &xxxuK, BzzVector &xxxlogFcent, BzzVector &xxxkFallOff);
	virtual void giveReactionRates(double cTot, BzzVector &c, BzzVector &R);
	virtual void giveJacobian(BzzVector &c, BzzMatrix &J);
	virtual void SetAccurateJacobian();
	virtual void UnsetAccurateJacobian();

private:

	 double k1, uK1, kFallOff1, logFcent1, r1;
 double k2, uK2, kFallOff2, logFcent2, r2;
 double k3, uK3, kFallOff3, logFcent3, r3;
 double k4, uK4, kFallOff4, logFcent4, r4;
 double k5, uK5, kFallOff5, logFcent5, r5;
 double k6, uK6, kFallOff6, logFcent6, r6;
 double k7, uK7, kFallOff7, logFcent7, r7;
 double k8, uK8, kFallOff8, logFcent8, r8;
 double k9, uK9, kFallOff9, logFcent9, r9;
 double k10, uK10, kFallOff10, logFcent10, r10;
 double k11, uK11, kFallOff11, logFcent11, r11;
 double k12, uK12, kFallOff12, logFcent12, r12;
 double k13, uK13, kFallOff13, logFcent13, r13;
 double k14, uK14, kFallOff14, logFcent14, r14;
 double k15, uK15, kFallOff15, logFcent15, r15;
 double k16, uK16, kFallOff16, logFcent16, r16;
 double k17, uK17, kFallOff17, logFcent17, r17;
 double k18, uK18, kFallOff18, logFcent18, r18;
 double k19, uK19, kFallOff19, logFcent19, r19;
 double k20, uK20, kFallOff20, logFcent20, r20;
 double k21, uK21, kFallOff21, logFcent21, r21;
 double k22, uK22, kFallOff22, logFcent22, r22;
 double k23, uK23, kFallOff23, logFcent23, r23;
 double k24, uK24, kFallOff24, logFcent24, r24;
 double k25, uK25, kFallOff25, logFcent25, r25;
 double k26, uK26, kFallOff26, logFcent26, r26;
 double k27, uK27, kFallOff27, logFcent27, r27;
 double k28, uK28, kFallOff28, logFcent28, r28;
 double k29, uK29, kFallOff29, logFcent29, r29;
 double k30, uK30, kFallOff30, logFcent30, r30;
 double k31, uK31, kFallOff31, logFcent31, r31;
 double k32, uK32, kFallOff32, logFcent32, r32;
 double k33, uK33, kFallOff33, logFcent33, r33;
 double k34, uK34, kFallOff34, logFcent34, r34;
 double k35, uK35, kFallOff35, logFcent35, r35;
 double k36, uK36, kFallOff36, logFcent36, r36;
 double k37, uK37, kFallOff37, logFcent37, r37;
 double k38, uK38, kFallOff38, logFcent38, r38;
 double k39, uK39, kFallOff39, logFcent39, r39;
 double k40, uK40, kFallOff40, logFcent40, r40;
 double k41, uK41, kFallOff41, logFcent41, r41;
 double k42, uK42, kFallOff42, logFcent42, r42;
 double k43, uK43, kFallOff43, logFcent43, r43;
 double k44, uK44, kFallOff44, logFcent44, r44;
 double k45, uK45, kFallOff45, logFcent45, r45;
 double k46, uK46, kFallOff46, logFcent46, r46;
 double k47, uK47, kFallOff47, logFcent47, r47;
 double k48, uK48, kFallOff48, logFcent48, r48;
 double k49, uK49, kFallOff49, logFcent49, r49;
 double k50, uK50, kFallOff50, logFcent50, r50;
 double k51, uK51, kFallOff51, logFcent51, r51;
 double k52, uK52, kFallOff52, logFcent52, r52;
 double k53, uK53, kFallOff53, logFcent53, r53;
 double k54, uK54, kFallOff54, logFcent54, r54;
 double k55, uK55, kFallOff55, logFcent55, r55;
 double k56, uK56, kFallOff56, logFcent56, r56;
 double k57, uK57, kFallOff57, logFcent57, r57;
 double k58, uK58, kFallOff58, logFcent58, r58;
 double k59, uK59, kFallOff59, logFcent59, r59;
 double k60, uK60, kFallOff60, logFcent60, r60;
 double k61, uK61, kFallOff61, logFcent61, r61;
 double k62, uK62, kFallOff62, logFcent62, r62;
 double k63, uK63, kFallOff63, logFcent63, r63;
 double k64, uK64, kFallOff64, logFcent64, r64;
 double k65, uK65, kFallOff65, logFcent65, r65;
 double k66, uK66, kFallOff66, logFcent66, r66;
 double k67, uK67, kFallOff67, logFcent67, r67;
 double k68, uK68, kFallOff68, logFcent68, r68;
 double k69, uK69, kFallOff69, logFcent69, r69;
 double k70, uK70, kFallOff70, logFcent70, r70;
 double k71, uK71, kFallOff71, logFcent71, r71;
 double k72, uK72, kFallOff72, logFcent72, r72;
 double k73, uK73, kFallOff73, logFcent73, r73;
 double k74, uK74, kFallOff74, logFcent74, r74;
 double k75, uK75, kFallOff75, logFcent75, r75;
 double k76, uK76, kFallOff76, logFcent76, r76;
 double k77, uK77, kFallOff77, logFcent77, r77;
 double k78, uK78, kFallOff78, logFcent78, r78;
 double k79, uK79, kFallOff79, logFcent79, r79;
 double k80, uK80, kFallOff80, logFcent80, r80;
 double k81, uK81, kFallOff81, logFcent81, r81;
 double k82, uK82, kFallOff82, logFcent82, r82;
 double k83, uK83, kFallOff83, logFcent83, r83;
 double k84, uK84, kFallOff84, logFcent84, r84;
 double k85, uK85, kFallOff85, logFcent85, r85;
 double k86, uK86, kFallOff86, logFcent86, r86;
 double k87, uK87, kFallOff87, logFcent87, r87;
 double k88, uK88, kFallOff88, logFcent88, r88;
 double k89, uK89, kFallOff89, logFcent89, r89;
 double k90, uK90, kFallOff90, logFcent90, r90;
 double k91, uK91, kFallOff91, logFcent91, r91;
 double k92, uK92, kFallOff92, logFcent92, r92;
 double k93, uK93, kFallOff93, logFcent93, r93;
 double k94, uK94, kFallOff94, logFcent94, r94;
 double k95, uK95, kFallOff95, logFcent95, r95;
 double k96, uK96, kFallOff96, logFcent96, r96;
 double k97, uK97, kFallOff97, logFcent97, r97;
 double k98, uK98, kFallOff98, logFcent98, r98;
 double k99, uK99, kFallOff99, logFcent99, r99;
 double k100, uK100, kFallOff100, logFcent100, r100;
 double k101, uK101, kFallOff101, logFcent101, r101;
 double k102, uK102, kFallOff102, logFcent102, r102;
 double k103, uK103, kFallOff103, logFcent103, r103;
 double k104, uK104, kFallOff104, logFcent104, r104;
 double k105, uK105, kFallOff105, logFcent105, r105;
 double k106, uK106, kFallOff106, logFcent106, r106;
 double k107, uK107, kFallOff107, logFcent107, r107;
 double k108, uK108, kFallOff108, logFcent108, r108;
 double k109, uK109, kFallOff109, logFcent109, r109;
 double k110, uK110, kFallOff110, logFcent110, r110;
 double k111, uK111, kFallOff111, logFcent111, r111;
 double k112, uK112, kFallOff112, logFcent112, r112;
 double k113, uK113, kFallOff113, logFcent113, r113;
 double k114, uK114, kFallOff114, logFcent114, r114;
 double k115, uK115, kFallOff115, logFcent115, r115;
 double k116, uK116, kFallOff116, logFcent116, r116;
 double k117, uK117, kFallOff117, logFcent117, r117;
 double k118, uK118, kFallOff118, logFcent118, r118;
 double k119, uK119, kFallOff119, logFcent119, r119;
 double k120, uK120, kFallOff120, logFcent120, r120;
 double k121, uK121, kFallOff121, logFcent121, r121;
 double k122, uK122, kFallOff122, logFcent122, r122;
 double k123, uK123, kFallOff123, logFcent123, r123;
 double k124, uK124, kFallOff124, logFcent124, r124;
 double k125, uK125, kFallOff125, logFcent125, r125;
 double k126, uK126, kFallOff126, logFcent126, r126;
 double k127, uK127, kFallOff127, logFcent127, r127;
 double k128, uK128, kFallOff128, logFcent128, r128;
 double k129, uK129, kFallOff129, logFcent129, r129;
 double k130, uK130, kFallOff130, logFcent130, r130;
 double k131, uK131, kFallOff131, logFcent131, r131;
 double k132, uK132, kFallOff132, logFcent132, r132;
 double k133, uK133, kFallOff133, logFcent133, r133;
 double k134, uK134, kFallOff134, logFcent134, r134;
 double k135, uK135, kFallOff135, logFcent135, r135;
 double k136, uK136, kFallOff136, logFcent136, r136;
 double k137, uK137, kFallOff137, logFcent137, r137;
 double k138, uK138, kFallOff138, logFcent138, r138;
 double k139, uK139, kFallOff139, logFcent139, r139;
 double k140, uK140, kFallOff140, logFcent140, r140;
 double k141, uK141, kFallOff141, logFcent141, r141;
 double k142, uK142, kFallOff142, logFcent142, r142;
 double k143, uK143, kFallOff143, logFcent143, r143;
 double k144, uK144, kFallOff144, logFcent144, r144;
 double k145, uK145, kFallOff145, logFcent145, r145;
 double k146, uK146, kFallOff146, logFcent146, r146;
 double k147, uK147, kFallOff147, logFcent147, r147;
 double k148, uK148, kFallOff148, logFcent148, r148;
 double k149, uK149, kFallOff149, logFcent149, r149;
 double k150, uK150, kFallOff150, logFcent150, r150;
 double k151, uK151, kFallOff151, logFcent151, r151;
 double k152, uK152, kFallOff152, logFcent152, r152;
 double k153, uK153, kFallOff153, logFcent153, r153;
 double k154, uK154, kFallOff154, logFcent154, r154;
 double k155, uK155, kFallOff155, logFcent155, r155;

 double coeffM1;
 double coeffM9;
 double coeffM20;
 double coeffM25;
 double coeffM29;
 double coeffM33;
 double coeffM34;
 double coeffM36;
 double coeffM38;
 double coeffM41;
 double coeffM42;
 double coeffM44;
 double coeffM46;
 double coeffM48;
 double coeffM50;
 double coeffM90;
 double coeffM97;
 double coeffM101;
 double coeffM117;
 double coeffM121;
 double coeffM146;
 double coeffM152;

 double coeffFallOff33; 
 double wF33; 
 double CFO33; 
 double dCFOdM33; 
 double dwFdM33; 
 double sigma33; 
 double lnPr33; 
 double coeffFallOff34; 
 double wF34; 
 double CFO34; 
 double dCFOdM34; 
 double dwFdM34; 
 double sigma34; 
 double lnPr34; 
 double coeffFallOff36; 
 double wF36; 
 double CFO36; 
 double dCFOdM36; 
 double dwFdM36; 
 double sigma36; 
 double lnPr36; 
 double coeffFallOff38; 
 double wF38; 
 double CFO38; 
 double dCFOdM38; 
 double dwFdM38; 
 double sigma38; 
 double lnPr38; 
 double coeffFallOff41; 
 double wF41; 
 double CFO41; 
 double dCFOdM41; 
 double dwFdM41; 
 double sigma41; 
 double lnPr41; 
 double coeffFallOff42; 
 double wF42; 
 double CFO42; 
 double dCFOdM42; 
 double dwFdM42; 
 double sigma42; 
 double lnPr42; 
 double coeffFallOff44; 
 double wF44; 
 double CFO44; 
 double dCFOdM44; 
 double dwFdM44; 
 double sigma44; 
 double lnPr44; 
 double coeffFallOff46; 
 double wF46; 
 double CFO46; 
 double dCFOdM46; 
 double dwFdM46; 
 double sigma46; 
 double lnPr46; 
 double coeffFallOff48; 
 double wF48; 
 double CFO48; 
 double dCFOdM48; 
 double dwFdM48; 
 double sigma48; 
 double lnPr48; 
 double coeffFallOff50; 
 double wF50; 
 double CFO50; 
 double dCFOdM50; 
 double dwFdM50; 
 double sigma50; 
 double lnPr50; 
 double coeffFallOff90; 
 double wF90; 
 double CFO90; 
 double dCFOdM90; 
 double dwFdM90; 
 double sigma90; 
 double lnPr90; 
 double coeffFallOff101; 
 double wF101; 
 double CFO101; 
 double dCFOdM101; 
 double dwFdM101; 
 double sigma101; 
 double lnPr101; 
 double coeffFallOff152; 
 double wF152; 
 double CFO152; 
 double dCFOdM152; 
 double dwFdM152; 
 double sigma152; 

 double rFlat1; 
 double rFlat9; 
 double rFlat20; 
 double rFlat25; 
 double rFlat29; 
 double rFlat33; 
 double rFlat34; 
 double rFlat36; 
 double rFlat38; 
 double rFlat41; 
 double rFlat42; 
 double rFlat44; 
 double rFlat46; 
 double rFlat48; 
 double rFlat50; 
 double rFlat90; 
 double rFlat97; 
 double rFlat101; 
 double rFlat117; 
 double rFlat121; 
 double rFlat146; 
 double rFlat152; 

  double d1d1;
  double d1d2;
  double d1d3;
  double d1d4;
  double d1d5;
  double d1d6;
  double d1d7;
  double d1d8;
  double d1d9;
  double d1d10;
  double d1d11;
  double d1d12;
  double d1d13;
  double d1d14;
  double d1d15;
  double d1d16;
  double d1d17;
  double d1d18;
  double d1d19;
  double d1d20;
  double d1d21;
  double d1d22;
  double d1d23;
  double d1d24;
  double d1d25;
  double d1d26;
  double d1d27;
  double d1d28;
  double d1d29;
  double d1d30;
  double d1d31;
  double d1d32;
  double d1d33;
  double d1d34;
  double d1d35;
  double d2d6;
  double d2d11;
  double d2d12;
  double d2d13;
  double d3d5;
  double d3d6;
  double d3d7;
  double d3d13;
  double d4d2;
  double d4d8;
  double d4d12;
  double d4d13;
  double d5d2;
  double d5d8;
  double d5d11;
  double d5d13;
  double d6d9;
  double d6d12;
  double d6d13;
  double d6d14;
  double d7d10;
  double d7d12;
  double d7d13;
  double d7d15;
  double d8d1;
  double d8d6;
  double d8d10;
  double d8d13;
  double d9d1;
  double d9d2;
  double d9d3;
  double d9d4;
  double d9d5;
  double d9d6;
  double d9d7;
  double d9d8;
  double d9d9;
  double d9d10;
  double d9d11;
  double d9d12;
  double d9d13;
  double d9d14;
  double d9d15;
  double d9d16;
  double d9d17;
  double d9d18;
  double d9d19;
  double d9d20;
  double d9d21;
  double d9d22;
  double d9d23;
  double d9d24;
  double d9d25;
  double d9d26;
  double d9d27;
  double d9d28;
  double d9d29;
  double d9d30;
  double d9d31;
  double d9d32;
  double d9d33;
  double d9d34;
  double d9d35;
  double d10d2;
  double d10d6;
  double d10d13;
  double d10d14;
  double d11d3;
  double d11d12;
  double d11d13;
  double d11d14;
  double d12d6;
  double d12d13;
  double d12d14;
  double d12d15;
  double d13d2;
  double d13d9;
  double d13d13;
  double d13d21;
  double d14d2;
  double d14d8;
  double d14d13;
  double d14d21;
  double d15d10;
  double d15d13;
  double d15d14;
  double d15d17;
  double d16d10;
  double d16d13;
  double d16d15;
  double d16d18;
  double d17d6;
  double d17d13;
  double d17d18;
  double d17d19;
  double d18d2;
  double d18d3;
  double d18d5;
  double d18d13;
  double d19d5;
  double d19d7;
  double d19d14;
  double d19d15;
  double d20d1;
  double d20d2;
  double d20d3;
  double d20d5;
  double d20d6;
  double d20d7;
  double d20d8;
  double d20d9;
  double d20d10;
  double d20d11;
  double d20d12;
  double d20d13;
  double d20d14;
  double d20d15;
  double d20d16;
  double d20d17;
  double d20d18;
  double d20d19;
  double d20d20;
  double d20d21;
  double d20d22;
  double d20d23;
  double d20d25;
  double d20d26;
  double d20d27;
  double d20d28;
  double d20d29;
  double d20d30;
  double d20d31;
  double d20d32;
  double d20d33;
  double d20d34;
  double d20d35;
  double d21d5;
  double d21d7;
  double d21d12;
  double d22d4;
  double d22d5;
  double d22d7;
  double d22d12;
  double d23d5;
  double d23d7;
  double d23d12;
  double d23d24;
  double d24d5;
  double d24d6;
  double d24d12;
  double d24d13;
  double d25d1;
  double d25d2;
  double d25d5;
  double d25d6;
  double d25d7;
  double d25d8;
  double d25d9;
  double d25d10;
  double d25d11;
  double d25d12;
  double d25d13;
  double d25d14;
  double d25d15;
  double d25d16;
  double d25d17;
  double d25d18;
  double d25d19;
  double d25d20;
  double d25d21;
  double d25d22;
  double d25d23;
  double d25d24;
  double d25d25;
  double d25d26;
  double d25d27;
  double d25d28;
  double d25d29;
  double d25d30;
  double d25d31;
  double d25d32;
  double d25d33;
  double d25d34;
  double d25d35;
  double d26d11;
  double d26d12;
  double d27d4;
  double d27d11;
  double d27d12;
  double d28d3;
  double d28d11;
  double d28d12;
  double d29d1;
  double d29d2;
  double d29d3;
  double d29d4;
  double d29d5;
  double d29d6;
  double d29d7;
  double d29d8;
  double d29d9;
  double d29d10;
  double d29d11;
  double d29d12;
  double d29d13;
  double d29d14;
  double d29d15;
  double d29d16;
  double d29d17;
  double d29d18;
  double d29d19;
  double d29d20;
  double d29d21;
  double d29d22;
  double d29d23;
  double d29d24;
  double d29d25;
  double d29d26;
  double d29d27;
  double d29d28;
  double d29d29;
  double d29d30;
  double d29d31;
  double d29d32;
  double d29d33;
  double d29d34;
  double d29d35;
  double d30d5;
  double d30d7;
  double d30d11;
  double d30d12;
  double d31d6;
  double d31d7;
  double d31d12;
  double d32d7;
  double d32d11;
  double d32d12;
  double d32d20;
  double d33d1;
  double d33d2;
  double d33d3;
  double d33d4;
  double d33d5;
  double d33d6;
  double d33d7;
  double d33d8;
  double d33d9;
  double d33d10;
  double d33d11;
  double d33d12;
  double d33d13;
  double d33d14;
  double d33d15;
  double d33d16;
  double d33d17;
  double d33d18;
  double d33d19;
  double d33d20;
  double d33d21;
  double d33d22;
  double d33d23;
  double d33d24;
  double d33d25;
  double d33d26;
  double d33d27;
  double d33d28;
  double d33d29;
  double d33d30;
  double d33d31;
  double d33d32;
  double d33d33;
  double d33d34;
  double d33d35;
  double d34d1;
  double d34d2;
  double d34d3;
  double d34d4;
  double d34d5;
  double d34d6;
  double d34d7;
  double d34d8;
  double d34d9;
  double d34d10;
  double d34d11;
  double d34d12;
  double d34d13;
  double d34d14;
  double d34d15;
  double d34d16;
  double d34d17;
  double d34d18;
  double d34d19;
  double d34d20;
  double d34d21;
  double d34d22;
  double d34d23;
  double d34d24;
  double d34d25;
  double d34d26;
  double d34d27;
  double d34d28;
  double d34d29;
  double d34d30;
  double d34d31;
  double d34d32;
  double d34d33;
  double d34d34;
  double d34d35;
  double d35d1;
  double d35d10;
  double d35d11;
  double d35d12;
  double d36d1;
  double d36d2;
  double d36d3;
  double d36d4;
  double d36d5;
  double d36d6;
  double d36d7;
  double d36d8;
  double d36d9;
  double d36d10;
  double d36d11;
  double d36d12;
  double d36d13;
  double d36d14;
  double d36d15;
  double d36d16;
  double d36d17;
  double d36d18;
  double d36d19;
  double d36d20;
  double d36d21;
  double d36d22;
  double d36d23;
  double d36d24;
  double d36d25;
  double d36d26;
  double d36d27;
  double d36d28;
  double d36d29;
  double d36d30;
  double d36d31;
  double d36d32;
  double d36d33;
  double d36d34;
  double d36d35;
  double d37d2;
  double d37d11;
  double d37d12;
  double d37d14;
  double d38d1;
  double d38d2;
  double d38d3;
  double d38d4;
  double d38d5;
  double d38d6;
  double d38d7;
  double d38d8;
  double d38d9;
  double d38d10;
  double d38d11;
  double d38d12;
  double d38d13;
  double d38d14;
  double d38d15;
  double d38d16;
  double d38d17;
  double d38d18;
  double d38d19;
  double d38d20;
  double d38d21;
  double d38d22;
  double d38d23;
  double d38d24;
  double d38d25;
  double d38d26;
  double d38d27;
  double d38d28;
  double d38d29;
  double d38d30;
  double d38d31;
  double d38d32;
  double d38d33;
  double d38d34;
  double d38d35;
  double d39d11;
  double d39d12;
  double d39d14;
  double d39d15;
  double d40d6;
  double d40d10;
  double d40d12;
  double d40d16;
  double d41d1;
  double d41d2;
  double d41d3;
  double d41d4;
  double d41d5;
  double d41d6;
  double d41d7;
  double d41d8;
  double d41d9;
  double d41d10;
  double d41d11;
  double d41d12;
  double d41d13;
  double d41d14;
  double d41d15;
  double d41d16;
  double d41d17;
  double d41d18;
  double d41d19;
  double d41d20;
  double d41d21;
  double d41d22;
  double d41d23;
  double d41d24;
  double d41d25;
  double d41d26;
  double d41d27;
  double d41d28;
  double d41d29;
  double d41d30;
  double d41d31;
  double d41d32;
  double d41d33;
  double d41d34;
  double d41d35;
  double d42d1;
  double d42d2;
  double d42d3;
  double d42d4;
  double d42d5;
  double d42d6;
  double d42d7;
  double d42d8;
  double d42d9;
  double d42d10;
  double d42d11;
  double d42d12;
  double d42d13;
  double d42d14;
  double d42d15;
  double d42d16;
  double d42d17;
  double d42d18;
  double d42d19;
  double d42d20;
  double d42d21;
  double d42d22;
  double d42d23;
  double d42d24;
  double d42d25;
  double d42d26;
  double d42d27;
  double d42d28;
  double d42d29;
  double d42d30;
  double d42d31;
  double d42d32;
  double d42d33;
  double d42d34;
  double d42d35;
  double d43d11;
  double d43d12;
  double d43d21;
  double d43d22;
  double d44d1;
  double d44d2;
  double d44d3;
  double d44d4;
  double d44d5;
  double d44d6;
  double d44d7;
  double d44d8;
  double d44d9;
  double d44d10;
  double d44d11;
  double d44d12;
  double d44d13;
  double d44d14;
  double d44d15;
  double d44d16;
  double d44d17;
  double d44d18;
  double d44d19;
  double d44d20;
  double d44d21;
  double d44d22;
  double d44d23;
  double d44d24;
  double d44d25;
  double d44d26;
  double d44d27;
  double d44d28;
  double d44d29;
  double d44d30;
  double d44d31;
  double d44d32;
  double d44d33;
  double d44d34;
  double d44d35;
  double d45d11;
  double d45d12;
  double d45d17;
  double d45d22;
  double d46d1;
  double d46d2;
  double d46d3;
  double d46d4;
  double d46d5;
  double d46d6;
  double d46d7;
  double d46d8;
  double d46d9;
  double d46d10;
  double d46d11;
  double d46d12;
  double d46d13;
  double d46d14;
  double d46d15;
  double d46d16;
  double d46d17;
  double d46d18;
  double d46d19;
  double d46d20;
  double d46d21;
  double d46d22;
  double d46d23;
  double d46d24;
  double d46d25;
  double d46d26;
  double d46d27;
  double d46d28;
  double d46d29;
  double d46d30;
  double d46d31;
  double d46d32;
  double d46d33;
  double d46d34;
  double d46d35;
  double d47d11;
  double d47d12;
  double d47d18;
  double d47d19;
  double d48d1;
  double d48d2;
  double d48d3;
  double d48d4;
  double d48d5;
  double d48d6;
  double d48d7;
  double d48d8;
  double d48d9;
  double d48d10;
  double d48d11;
  double d48d12;
  double d48d13;
  double d48d14;
  double d48d15;
  double d48d16;
  double d48d17;
  double d48d18;
  double d48d19;
  double d48d20;
  double d48d21;
  double d48d22;
  double d48d23;
  double d48d24;
  double d48d25;
  double d48d26;
  double d48d27;
  double d48d28;
  double d48d29;
  double d48d30;
  double d48d31;
  double d48d32;
  double d48d33;
  double d48d34;
  double d48d35;
  double d49d4;
  double d49d6;
  double d49d11;
  double d49d12;
  double d50d1;
  double d50d2;
  double d50d3;
  double d50d4;
  double d50d5;
  double d50d6;
  double d50d7;
  double d50d8;
  double d50d9;
  double d50d10;
  double d50d11;
  double d50d12;
  double d50d13;
  double d50d14;
  double d50d15;
  double d50d16;
  double d50d17;
  double d50d18;
  double d50d19;
  double d50d20;
  double d50d21;
  double d50d22;
  double d50d23;
  double d50d24;
  double d50d25;
  double d50d26;
  double d50d27;
  double d50d28;
  double d50d29;
  double d50d30;
  double d50d31;
  double d50d32;
  double d50d33;
  double d50d34;
  double d50d35;
  double d51d4;
  double d51d6;
  double d51d13;
  double d52d4;
  double d52d5;
  double d52d6;
  double d52d7;
  double d53d4;
  double d53d6;
  double d53d7;
  double d53d20;
  double d54d6;
  double d54d8;
  double d54d12;
  double d54d15;
  double d55d6;
  double d55d9;
  double d55d12;
  double d55d15;
  double d56d4;
  double d56d6;
  double d56d9;
  double d56d10;
  double d57d1;
  double d57d4;
  double d57d6;
  double d57d10;
  double d58d2;
  double d58d3;
  double d58d6;
  double d58d12;
  double d59d2;
  double d59d4;
  double d59d6;
  double d59d14;
  double d60d4;
  double d60d6;
  double d60d14;
  double d60d15;
  double d61d2;
  double d61d6;
  double d61d10;
  double d61d21;
  double d62d4;
  double d62d6;
  double d62d21;
  double d62d22;
  double d63d4;
  double d63d6;
  double d63d17;
  double d63d22;
  double d64d4;
  double d64d6;
  double d64d18;
  double d64d19;
  double d65d5;
  double d65d7;
  double d65d20;
  double d66d5;
  double d66d7;
  double d66d20;
  double d67d6;
  double d67d7;
  double d67d8;
  double d67d15;
  double d68d1;
  double d68d5;
  double d68d7;
  double d68d10;
  double d69d6;
  double d69d7;
  double d69d10;
  double d69d16;
  double d70d2;
  double d70d3;
  double d70d6;
  double d70d7;
  double d71d7;
  double d71d14;
  double d71d15;
  double d71d20;
  double d72d5;
  double d72d6;
  double d72d8;
  double d72d14;
  double d73d8;
  double d73d10;
  double d73d11;
  double d73d12;
  double d74d8;
  double d74d11;
  double d74d21;
  double d75d8;
  double d75d10;
  double d75d12;
  double d75d17;
  double d76d1;
  double d76d8;
  double d76d10;
  double d77d8;
  double d77d9;
  double d77d24;
  double d78d2;
  double d78d5;
  double d78d6;
  double d78d9;
  double d78d12;
  double d79d2;
  double d79d4;
  double d79d5;
  double d79d9;
  double d80d9;
  double d80d10;
  double d80d11;
  double d80d12;
  double d81d4;
  double d81d8;
  double d81d9;
  double d82d9;
  double d82d10;
  double d82d12;
  double d82d17;
  double d83d1;
  double d83d9;
  double d83d10;
  double d84d2;
  double d84d8;
  double d84d9;
  double d85d3;
  double d85d8;
  double d85d9;
  double d86d2;
  double d86d3;
  double d86d9;
  double d86d15;
  double d87d5;
  double d87d10;
  double d87d13;
  double d87d16;
  double d88d5;
  double d88d6;
  double d88d10;
  double d88d15;
  double d89d1;
  double d89d7;
  double d89d10;
  double d89d20;
  double d90d1;
  double d90d2;
  double d90d3;
  double d90d4;
  double d90d5;
  double d90d6;
  double d90d7;
  double d90d8;
  double d90d9;
  double d90d10;
  double d90d11;
  double d90d12;
  double d90d13;
  double d90d14;
  double d90d15;
  double d90d16;
  double d90d17;
  double d90d18;
  double d90d19;
  double d90d20;
  double d90d21;
  double d90d22;
  double d90d23;
  double d90d24;
  double d90d25;
  double d90d26;
  double d90d27;
  double d90d28;
  double d90d29;
  double d90d30;
  double d90d31;
  double d90d32;
  double d90d33;
  double d90d34;
  double d90d35;
  double d91d10;
  double d91d12;
  double d91d18;
  double d92d1;
  double d92d2;
  double d92d10;
  double d92d14;
  double d93d1;
  double d93d10;
  double d93d14;
  double d93d15;
  double d94d1;
  double d94d10;
  double d94d17;
  double d94d22;
  double d95d1;
  double d95d10;
  double d95d18;
  double d95d19;
  double d96d2;
  double d96d4;
  double d96d12;
  double d96d14;
  double d97d1;
  double d97d2;
  double d97d3;
  double d97d5;
  double d97d6;
  double d97d7;
  double d97d8;
  double d97d9;
  double d97d10;
  double d97d11;
  double d97d12;
  double d97d13;
  double d97d14;
  double d97d15;
  double d97d16;
  double d97d17;
  double d97d18;
  double d97d19;
  double d97d20;
  double d97d21;
  double d97d22;
  double d97d23;
  double d97d24;
  double d97d25;
  double d97d26;
  double d97d27;
  double d97d28;
  double d97d29;
  double d97d30;
  double d97d31;
  double d97d32;
  double d97d33;
  double d97d34;
  double d97d35;
  double d98d2;
  double d98d5;
  double d98d7;
  double d98d14;
  double d99d5;
  double d99d7;
  double d99d15;
  double d99d16;
  double d100d5;
  double d100d14;
  double d100d15;
  double d100d22;
  double d101d1;
  double d101d2;
  double d101d3;
  double d101d4;
  double d101d5;
  double d101d6;
  double d101d7;
  double d101d8;
  double d101d9;
  double d101d10;
  double d101d11;
  double d101d12;
  double d101d13;
  double d101d14;
  double d101d15;
  double d101d16;
  double d101d17;
  double d101d18;
  double d101d19;
  double d101d20;
  double d101d21;
  double d101d22;
  double d101d23;
  double d101d24;
  double d101d25;
  double d101d26;
  double d101d27;
  double d101d28;
  double d101d29;
  double d101d30;
  double d101d31;
  double d101d32;
  double d101d33;
  double d101d34;
  double d101d35;
  double d102d5;
  double d102d7;
  double d102d17;
  double d102d18;
  double d103d8;
  double d103d11;
  double d103d12;
  double d103d23;
  double d104d4;
  double d104d6;
  double d104d8;
  double d104d23;
  double d105d9;
  double d105d11;
  double d105d12;
  double d105d23;
  double d106d1;
  double d106d12;
  double d106d17;
  double d106d23;
  double d107d2;
  double d107d12;
  double d107d13;
  double d107d23;
  double d108d6;
  double d108d12;
  double d108d14;
  double d108d23;
  double d109d5;
  double d109d13;
  double d109d14;
  double d109d23;
  double d110d4;
  double d110d12;
  double d110d15;
  double d110d23;
  double d111d2;
  double d111d3;
  double d111d14;
  double d111d23;
  double d112d13;
  double d112d24;
  double d112d25;
  double d112d26;
  double d113d5;
  double d113d13;
  double d113d25;
  double d113d26;
  double d114d6;
  double d114d12;
  double d114d25;
  double d114d26;
  double d115d23;
  double d115d24;
  double d115d26;
  double d115d27;
  double d116d12;
  double d116d13;
  double d116d27;
  double d116d28;
  double d117d1;
  double d117d2;
  double d117d3;
  double d117d4;
  double d117d5;
  double d117d6;
  double d117d7;
  double d117d8;
  double d117d9;
  double d117d10;
  double d117d11;
  double d117d12;
  double d117d13;
  double d117d14;
  double d117d15;
  double d117d16;
  double d117d17;
  double d117d18;
  double d117d19;
  double d117d20;
  double d117d21;
  double d117d22;
  double d117d23;
  double d117d24;
  double d117d25;
  double d117d26;
  double d117d27;
  double d117d28;
  double d117d29;
  double d117d30;
  double d117d31;
  double d117d32;
  double d117d33;
  double d117d34;
  double d117d35;
  double d118d2;
  double d118d12;
  double d118d28;
  double d118d29;
  double d119d2;
  double d119d13;
  double d119d25;
  double d119d28;
  double d120d11;
  double d120d12;
  double d120d28;
  double d120d30;
  double d121d1;
  double d121d2;
  double d121d3;
  double d121d4;
  double d121d5;
  double d121d6;
  double d121d7;
  double d121d8;
  double d121d9;
  double d121d10;
  double d121d11;
  double d121d12;
  double d121d13;
  double d121d14;
  double d121d15;
  double d121d16;
  double d121d17;
  double d121d18;
  double d121d19;
  double d121d20;
  double d121d21;
  double d121d22;
  double d121d23;
  double d121d24;
  double d121d25;
  double d121d26;
  double d121d27;
  double d121d28;
  double d121d29;
  double d121d30;
  double d121d31;
  double d121d32;
  double d121d33;
  double d121d34;
  double d121d35;
  double d122d2;
  double d122d12;
  double d122d30;
  double d122d31;
  double d123d6;
  double d123d13;
  double d123d28;
  double d123d30;
  double d124d3;
  double d124d13;
  double d124d29;
  double d124d30;
  double d125d4;
  double d125d6;
  double d125d28;
  double d125d30;
  double d126d11;
  double d126d12;
  double d126d27;
  double d126d32;
  double d127d4;
  double d127d6;
  double d127d27;
  double d127d32;
  double d128d6;
  double d128d12;
  double d128d28;
  double d128d32;
  double d129d5;
  double d129d13;
  double d129d28;
  double d129d32;
  double d130d11;
  double d130d12;
  double d130d26;
  double d130d29;
  double d131d12;
  double d131d13;
  double d131d25;
  double d131d29;
  double d132d6;
  double d132d12;
  double d132d29;
  double d132d33;
  double d133d4;
  double d133d6;
  double d133d26;
  double d133d29;
  double d134d5;
  double d134d13;
  double d134d29;
  double d134d33;
  double d135d12;
  double d135d25;
  double d135d29;
  double d135d34;
  double d136d6;
  double d136d24;
  double d136d25;
  double d136d29;
  double d137d11;
  double d137d12;
  double d137d29;
  double d137d31;
  double d138d12;
  double d138d13;
  double d138d31;
  double d138d33;
  double d139d4;
  double d139d6;
  double d139d29;
  double d139d31;
  double d140d4;
  double d140d24;
  double d140d25;
  double d140d31;
  double d141d6;
  double d141d25;
  double d141d31;
  double d141d35;
  double d142d12;
  double d142d24;
  double d142d35;
  double d143d11;
  double d143d12;
  double d143d24;
  double d143d35;
  double d144d12;
  double d144d13;
  double d144d34;
  double d144d35;
  double d145d4;
  double d145d6;
  double d145d24;
  double d145d35;
  double d146d1;
  double d146d2;
  double d146d3;
  double d146d4;
  double d146d5;
  double d146d6;
  double d146d7;
  double d146d8;
  double d146d9;
  double d146d10;
  double d146d11;
  double d146d12;
  double d146d13;
  double d146d14;
  double d146d15;
  double d146d16;
  double d146d17;
  double d146d18;
  double d146d19;
  double d146d20;
  double d146d21;
  double d146d22;
  double d146d23;
  double d146d24;
  double d146d25;
  double d146d26;
  double d146d27;
  double d146d28;
  double d146d29;
  double d146d30;
  double d146d31;
  double d146d32;
  double d146d33;
  double d146d34;
  double d146d35;
  double d147d11;
  double d147d12;
  double d147d25;
  double d147d33;
  double d148d4;
  double d148d6;
  double d148d25;
  double d148d33;
  double d149d10;
  double d149d25;
  double d150d8;
  double d150d12;
  double d150d25;
  double d150d30;
  double d151d13;
  double d151d23;
  double d151d25;
  double d151d27;
  double d152d1;
  double d152d2;
  double d152d3;
  double d152d4;
  double d152d5;
  double d152d6;
  double d152d7;
  double d152d8;
  double d152d9;
  double d152d10;
  double d152d11;
  double d152d12;
  double d152d13;
  double d152d14;
  double d152d15;
  double d152d16;
  double d152d17;
  double d152d18;
  double d152d19;
  double d152d20;
  double d152d21;
  double d152d22;
  double d152d23;
  double d152d24;
  double d152d25;
  double d152d26;
  double d152d27;
  double d152d28;
  double d152d29;
  double d152d30;
  double d152d31;
  double d152d32;
  double d152d33;
  double d152d34;
  double d152d35;
  double d153d6;
  double d153d12;
  double d153d24;
  double d153d34;
  double d154d13;
  double d154d25;
  double d154d34;
  double d155d6;
  double d155d7;
  double d155d24;
  double d155d34;
};

#endif //(OPENSMOKE_KINPP_SYMBOLICKINETICS_FLUENT_DRM22_POLIMI_NOX_H)
