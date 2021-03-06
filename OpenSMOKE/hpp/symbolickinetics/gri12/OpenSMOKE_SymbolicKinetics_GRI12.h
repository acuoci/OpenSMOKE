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

#if !defined(OPENSMOKE_SYMBOLIC_KINETICS_GRI12)
#define OPENSMOKE_SYMBOLIC_KINETICS_GRI12

#include "BzzMath.hpp"
#include "symbolickinetics/OpenSMOKE_SymbolicKinetics.h"

class OpenSMOKE_SymbolicKinetics_GRI12 : public OpenSMOKE_SymbolicKinetics
{
public:
	
	OpenSMOKE_SymbolicKinetics_GRI12();

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
 double k156, uK156, kFallOff156, logFcent156, r156;
 double k157, uK157, kFallOff157, logFcent157, r157;
 double k158, uK158, kFallOff158, logFcent158, r158;
 double k159, uK159, kFallOff159, logFcent159, r159;
 double k160, uK160, kFallOff160, logFcent160, r160;
 double k161, uK161, kFallOff161, logFcent161, r161;
 double k162, uK162, kFallOff162, logFcent162, r162;
 double k163, uK163, kFallOff163, logFcent163, r163;
 double k164, uK164, kFallOff164, logFcent164, r164;
 double k165, uK165, kFallOff165, logFcent165, r165;
 double k166, uK166, kFallOff166, logFcent166, r166;
 double k167, uK167, kFallOff167, logFcent167, r167;
 double k168, uK168, kFallOff168, logFcent168, r168;
 double k169, uK169, kFallOff169, logFcent169, r169;
 double k170, uK170, kFallOff170, logFcent170, r170;
 double k171, uK171, kFallOff171, logFcent171, r171;
 double k172, uK172, kFallOff172, logFcent172, r172;
 double k173, uK173, kFallOff173, logFcent173, r173;
 double k174, uK174, kFallOff174, logFcent174, r174;
 double k175, uK175, kFallOff175, logFcent175, r175;
 double k176, uK176, kFallOff176, logFcent176, r176;
 double k177, uK177, kFallOff177, logFcent177, r177;

 double coeffM1;
 double coeffM2;
 double coeffM12;
 double coeffM33;
 double coeffM39;
 double coeffM43;
 double coeffM50;
 double coeffM52;
 double coeffM54;
 double coeffM56;
 double coeffM57;
 double coeffM59;
 double coeffM63;
 double coeffM70;
 double coeffM71;
 double coeffM72;
 double coeffM74;
 double coeffM76;
 double coeffM83;
 double coeffM85;
 double coeffM95;
 double coeffM131;
 double coeffM140;
 double coeffM147;
 double coeffM158;
 double coeffM167;
 double coeffM174;

 double coeffFallOff50; 
 double wF50; 
 double CFO50; 
 double dCFOdM50; 
 double dwFdM50; 
 double sigma50; 
 double lnPr50; 
 double coeffFallOff52; 
 double wF52; 
 double CFO52; 
 double dCFOdM52; 
 double dwFdM52; 
 double sigma52; 
 double lnPr52; 
 double coeffFallOff54; 
 double wF54; 
 double CFO54; 
 double dCFOdM54; 
 double dwFdM54; 
 double sigma54; 
 double lnPr54; 
 double coeffFallOff56; 
 double wF56; 
 double CFO56; 
 double dCFOdM56; 
 double dwFdM56; 
 double sigma56; 
 double lnPr56; 
 double coeffFallOff57; 
 double wF57; 
 double CFO57; 
 double dCFOdM57; 
 double dwFdM57; 
 double sigma57; 
 double lnPr57; 
 double coeffFallOff59; 
 double wF59; 
 double CFO59; 
 double dCFOdM59; 
 double dwFdM59; 
 double sigma59; 
 double lnPr59; 
 double coeffFallOff63; 
 double wF63; 
 double CFO63; 
 double dCFOdM63; 
 double dwFdM63; 
 double sigma63; 
 double lnPr63; 
 double coeffFallOff70; 
 double wF70; 
 double CFO70; 
 double dCFOdM70; 
 double dwFdM70; 
 double sigma70; 
 double lnPr70; 
 double coeffFallOff71; 
 double wF71; 
 double CFO71; 
 double dCFOdM71; 
 double dwFdM71; 
 double sigma71; 
 double lnPr71; 
 double coeffFallOff72; 
 double wF72; 
 double CFO72; 
 double dCFOdM72; 
 double dwFdM72; 
 double sigma72; 
 double lnPr72; 
 double coeffFallOff74; 
 double wF74; 
 double CFO74; 
 double dCFOdM74; 
 double dwFdM74; 
 double sigma74; 
 double lnPr74; 
 double coeffFallOff76; 
 double wF76; 
 double CFO76; 
 double dCFOdM76; 
 double dwFdM76; 
 double sigma76; 
 double lnPr76; 
 double coeffFallOff83; 
 double wF83; 
 double CFO83; 
 double dCFOdM83; 
 double dwFdM83; 
 double sigma83; 
 double lnPr83; 
 double coeffFallOff85; 
 double wF85; 
 double CFO85; 
 double dCFOdM85; 
 double dwFdM85; 
 double sigma85; 
 double lnPr85; 
 double coeffFallOff95; 
 double wF95; 
 double CFO95; 
 double dCFOdM95; 
 double dwFdM95; 
 double sigma95; 
 double lnPr95; 
 double coeffFallOff131; 
 double wF131; 
 double CFO131; 
 double dCFOdM131; 
 double dwFdM131; 
 double sigma131; 
 double lnPr131; 
 double coeffFallOff140; 
 double wF140; 
 double CFO140; 
 double dCFOdM140; 
 double dwFdM140; 
 double sigma140; 
 double lnPr140; 
 double coeffFallOff147; 
 double wF147; 
 double CFO147; 
 double dCFOdM147; 
 double dwFdM147; 
 double sigma147; 
 double lnPr147; 
 double coeffFallOff158; 
 double wF158; 
 double CFO158; 
 double dCFOdM158; 
 double dwFdM158; 
 double sigma158; 
 double lnPr158; 
 double coeffFallOff174; 
 double wF174; 
 double CFO174; 
 double dCFOdM174; 
 double dwFdM174; 
 double sigma174; 
 double lnPr174; 

 double rFlat1; 
 double rFlat2; 
 double rFlat12; 
 double rFlat33; 
 double rFlat39; 
 double rFlat43; 
 double rFlat50; 
 double rFlat52; 
 double rFlat54; 
 double rFlat56; 
 double rFlat57; 
 double rFlat59; 
 double rFlat63; 
 double rFlat70; 
 double rFlat71; 
 double rFlat72; 
 double rFlat74; 
 double rFlat76; 
 double rFlat83; 
 double rFlat85; 
 double rFlat95; 
 double rFlat131; 
 double rFlat140; 
 double rFlat147; 
 double rFlat158; 
 double rFlat167; 
 double rFlat174; 

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
  double d2d1;
  double d2d2;
  double d2d3;
  double d2d4;
  double d2d5;
  double d2d6;
  double d2d7;
  double d2d8;
  double d2d9;
  double d2d10;
  double d2d11;
  double d2d12;
  double d2d13;
  double d2d14;
  double d2d15;
  double d2d16;
  double d2d17;
  double d2d18;
  double d2d19;
  double d2d20;
  double d2d21;
  double d2d22;
  double d2d23;
  double d2d24;
  double d2d25;
  double d2d26;
  double d2d27;
  double d2d28;
  double d2d29;
  double d2d30;
  double d2d31;
  double d2d32;
  double d3d1;
  double d3d2;
  double d3d3;
  double d3d5;
  double d4d3;
  double d4d4;
  double d4d5;
  double d4d7;
  double d5d3;
  double d5d5;
  double d5d7;
  double d5d8;
  double d6d2;
  double d6d3;
  double d6d10;
  double d6d15;
  double d7d2;
  double d7d3;
  double d7d11;
  double d7d17;
  double d8d1;
  double d8d3;
  double d8d12;
  double d8d15;
  double d9d2;
  double d9d3;
  double d9d12;
  double d9d17;
  double d10d2;
  double d10d3;
  double d10d13;
  double d10d18;
  double d11d3;
  double d11d5;
  double d11d13;
  double d11d14;
  double d12d1;
  double d12d2;
  double d12d3;
  double d12d4;
  double d12d5;
  double d12d6;
  double d12d7;
  double d12d8;
  double d12d9;
  double d12d10;
  double d12d11;
  double d12d12;
  double d12d13;
  double d12d14;
  double d12d15;
  double d12d16;
  double d12d17;
  double d12d18;
  double d12d19;
  double d12d20;
  double d12d21;
  double d12d22;
  double d12d23;
  double d12d24;
  double d12d25;
  double d12d26;
  double d12d27;
  double d12d28;
  double d12d29;
  double d12d30;
  double d12d31;
  double d12d32;
  double d13d3;
  double d13d5;
  double d13d15;
  double d13d17;
  double d14d2;
  double d14d3;
  double d14d16;
  double d14d17;
  double d15d3;
  double d15d5;
  double d15d17;
  double d15d18;
  double d16d3;
  double d16d5;
  double d16d18;
  double d16d19;
  double d17d3;
  double d17d5;
  double d17d18;
  double d17d20;
  double d18d3;
  double d18d5;
  double d18d19;
  double d18d21;
  double d19d3;
  double d19d5;
  double d19d20;
  double d19d21;
  double d20d3;
  double d20d10;
  double d20d15;
  double d20d22;
  double d21d2;
  double d21d3;
  double d21d23;
  double d21d28;
  double d22d3;
  double d22d5;
  double d22d22;
  double d22d23;
  double d23d3;
  double d23d11;
  double d23d15;
  double d23d23;
  double d24d2;
  double d24d3;
  double d24d24;
  double d24d29;
  double d25d3;
  double d25d13;
  double d25d17;
  double d25d25;
  double d26d3;
  double d26d13;
  double d26d18;
  double d26d26;
  double d27d3;
  double d27d5;
  double d27d26;
  double d27d27;
  double d28d2;
  double d28d3;
  double d28d15;
  double d28d28;
  double d29d3;
  double d29d5;
  double d29d28;
  double d29d29;
  double d30d3;
  double d30d11;
  double d30d16;
  double d30d29;
  double d31d3;
  double d31d4;
  double d31d15;
  double d31d16;
  double d32d4;
  double d32d7;
  double d32d17;
  double d32d18;
  double d33d1;
  double d33d2;
  double d33d3;
  double d33d4;
  double d33d5;
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
  double d34d2;
  double d34d4;
  double d34d7;
  double d35d2;
  double d35d4;
  double d35d6;
  double d35d7;
  double d36d2;
  double d36d4;
  double d36d7;
  double d36d31;
  double d37d2;
  double d37d4;
  double d37d7;
  double d37d32;
  double d38d2;
  double d38d3;
  double d38d4;
  double d38d5;
  double d39d1;
  double d39d2;
  double d39d3;
  double d39d4;
  double d39d5;
  double d39d7;
  double d39d8;
  double d39d9;
  double d39d10;
  double d39d11;
  double d39d12;
  double d39d13;
  double d39d14;
  double d39d15;
  double d39d17;
  double d39d18;
  double d39d19;
  double d39d20;
  double d39d21;
  double d39d22;
  double d39d23;
  double d39d24;
  double d39d25;
  double d39d26;
  double d39d27;
  double d39d28;
  double d39d29;
  double d39d30;
  double d39d31;
  double d39d32;
  double d40d1;
  double d40d2;
  double d41d1;
  double d41d2;
  double d41d6;
  double d42d1;
  double d42d2;
  double d42d16;
  double d43d1;
  double d43d2;
  double d43d3;
  double d43d4;
  double d43d5;
  double d43d6;
  double d43d7;
  double d43d8;
  double d43d9;
  double d43d10;
  double d43d11;
  double d43d12;
  double d43d13;
  double d43d14;
  double d43d15;
  double d43d16;
  double d43d17;
  double d43d18;
  double d43d19;
  double d43d20;
  double d43d21;
  double d43d22;
  double d43d23;
  double d43d24;
  double d43d25;
  double d43d26;
  double d43d27;
  double d43d28;
  double d43d29;
  double d43d30;
  double d43d31;
  double d43d32;
  double d44d2;
  double d44d3;
  double d44d6;
  double d44d7;
  double d45d1;
  double d45d2;
  double d45d4;
  double d45d7;
  double d46d2;
  double d46d5;
  double d46d7;
  double d47d1;
  double d47d2;
  double d47d7;
  double d47d8;
  double d48d2;
  double d48d5;
  double d48d6;
  double d48d8;
  double d49d1;
  double d49d2;
  double d49d9;
  double d49d10;
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
  double d51d1;
  double d51d2;
  double d51d10;
  double d51d12;
  double d52d1;
  double d52d2;
  double d52d3;
  double d52d4;
  double d52d5;
  double d52d6;
  double d52d7;
  double d52d8;
  double d52d9;
  double d52d10;
  double d52d11;
  double d52d12;
  double d52d13;
  double d52d14;
  double d52d15;
  double d52d16;
  double d52d17;
  double d52d18;
  double d52d19;
  double d52d20;
  double d52d21;
  double d52d22;
  double d52d23;
  double d52d24;
  double d52d25;
  double d52d26;
  double d52d27;
  double d52d28;
  double d52d29;
  double d52d30;
  double d52d31;
  double d52d32;
  double d53d1;
  double d53d2;
  double d53d13;
  double d53d14;
  double d54d1;
  double d54d2;
  double d54d3;
  double d54d4;
  double d54d5;
  double d54d6;
  double d54d7;
  double d54d8;
  double d54d9;
  double d54d10;
  double d54d11;
  double d54d12;
  double d54d13;
  double d54d14;
  double d54d15;
  double d54d16;
  double d54d17;
  double d54d18;
  double d54d19;
  double d54d20;
  double d54d21;
  double d54d22;
  double d54d23;
  double d54d24;
  double d54d25;
  double d54d26;
  double d54d27;
  double d54d28;
  double d54d29;
  double d54d30;
  double d54d31;
  double d54d32;
  double d55d1;
  double d55d2;
  double d55d15;
  double d55d17;
  double d56d1;
  double d56d2;
  double d56d3;
  double d56d4;
  double d56d5;
  double d56d6;
  double d56d7;
  double d56d8;
  double d56d9;
  double d56d10;
  double d56d11;
  double d56d12;
  double d56d13;
  double d56d14;
  double d56d15;
  double d56d16;
  double d56d17;
  double d56d18;
  double d56d19;
  double d56d20;
  double d56d21;
  double d56d22;
  double d56d23;
  double d56d24;
  double d56d25;
  double d56d26;
  double d56d27;
  double d56d28;
  double d56d29;
  double d56d30;
  double d56d31;
  double d56d32;
  double d57d1;
  double d57d2;
  double d57d3;
  double d57d4;
  double d57d5;
  double d57d6;
  double d57d7;
  double d57d8;
  double d57d9;
  double d57d10;
  double d57d11;
  double d57d12;
  double d57d13;
  double d57d14;
  double d57d15;
  double d57d16;
  double d57d17;
  double d57d18;
  double d57d19;
  double d57d20;
  double d57d21;
  double d57d22;
  double d57d23;
  double d57d24;
  double d57d25;
  double d57d26;
  double d57d27;
  double d57d28;
  double d57d29;
  double d57d30;
  double d57d31;
  double d57d32;
  double d58d1;
  double d58d2;
  double d58d17;
  double d58d18;
  double d59d1;
  double d59d2;
  double d59d3;
  double d59d4;
  double d59d5;
  double d59d6;
  double d59d7;
  double d59d8;
  double d59d9;
  double d59d10;
  double d59d11;
  double d59d12;
  double d59d13;
  double d59d14;
  double d59d15;
  double d59d16;
  double d59d17;
  double d59d18;
  double d59d19;
  double d59d20;
  double d59d21;
  double d59d22;
  double d59d23;
  double d59d24;
  double d59d25;
  double d59d26;
  double d59d27;
  double d59d28;
  double d59d29;
  double d59d30;
  double d59d31;
  double d59d32;
  double d60d1;
  double d60d2;
  double d60d18;
  double d60d19;
  double d61d2;
  double d61d5;
  double d61d13;
  double d61d19;
  double d62d2;
  double d62d6;
  double d62d12;
  double d62d19;
  double d63d1;
  double d63d2;
  double d63d3;
  double d63d4;
  double d63d5;
  double d63d6;
  double d63d7;
  double d63d8;
  double d63d9;
  double d63d10;
  double d63d11;
  double d63d12;
  double d63d13;
  double d63d14;
  double d63d15;
  double d63d16;
  double d63d17;
  double d63d18;
  double d63d19;
  double d63d20;
  double d63d21;
  double d63d22;
  double d63d23;
  double d63d24;
  double d63d25;
  double d63d26;
  double d63d27;
  double d63d28;
  double d63d29;
  double d63d30;
  double d63d31;
  double d63d32;
  double d64d2;
  double d64d19;
  double d64d20;
  double d65d1;
  double d65d2;
  double d65d18;
  double d65d20;
  double d66d2;
  double d66d5;
  double d66d13;
  double d66d20;
  double d67d2;
  double d67d6;
  double d67d12;
  double d67d20;
  double d68d1;
  double d68d2;
  double d68d19;
  double d68d21;
  double d69d1;
  double d69d2;
  double d69d20;
  double d69d21;
  double d70d1;
  double d70d2;
  double d70d3;
  double d70d4;
  double d70d5;
  double d70d6;
  double d70d7;
  double d70d8;
  double d70d9;
  double d70d10;
  double d70d11;
  double d70d12;
  double d70d13;
  double d70d14;
  double d70d15;
  double d70d16;
  double d70d17;
  double d70d18;
  double d70d19;
  double d70d20;
  double d70d21;
  double d70d22;
  double d70d23;
  double d70d24;
  double d70d25;
  double d70d26;
  double d70d27;
  double d70d28;
  double d70d29;
  double d70d30;
  double d70d31;
  double d70d32;
  double d71d1;
  double d71d2;
  double d71d3;
  double d71d4;
  double d71d5;
  double d71d6;
  double d71d7;
  double d71d8;
  double d71d9;
  double d71d10;
  double d71d11;
  double d71d12;
  double d71d13;
  double d71d14;
  double d71d15;
  double d71d16;
  double d71d17;
  double d71d18;
  double d71d19;
  double d71d20;
  double d71d21;
  double d71d22;
  double d71d23;
  double d71d24;
  double d71d25;
  double d71d26;
  double d71d27;
  double d71d28;
  double d71d29;
  double d71d30;
  double d71d31;
  double d71d32;
  double d72d1;
  double d72d2;
  double d72d3;
  double d72d4;
  double d72d5;
  double d72d6;
  double d72d7;
  double d72d8;
  double d72d9;
  double d72d10;
  double d72d11;
  double d72d12;
  double d72d13;
  double d72d14;
  double d72d15;
  double d72d16;
  double d72d17;
  double d72d18;
  double d72d19;
  double d72d20;
  double d72d21;
  double d72d22;
  double d72d23;
  double d72d24;
  double d72d25;
  double d72d26;
  double d72d27;
  double d72d28;
  double d72d29;
  double d72d30;
  double d72d31;
  double d72d32;
  double d73d1;
  double d73d2;
  double d73d23;
  double d73d24;
  double d74d1;
  double d74d2;
  double d74d3;
  double d74d4;
  double d74d5;
  double d74d6;
  double d74d7;
  double d74d8;
  double d74d9;
  double d74d10;
  double d74d11;
  double d74d12;
  double d74d13;
  double d74d14;
  double d74d15;
  double d74d16;
  double d74d17;
  double d74d18;
  double d74d19;
  double d74d20;
  double d74d21;
  double d74d22;
  double d74d23;
  double d74d24;
  double d74d25;
  double d74d26;
  double d74d27;
  double d74d28;
  double d74d29;
  double d74d30;
  double d74d31;
  double d74d32;
  double d75d1;
  double d75d2;
  double d75d24;
  double d75d25;
  double d76d1;
  double d76d2;
  double d76d3;
  double d76d4;
  double d76d5;
  double d76d6;
  double d76d7;
  double d76d8;
  double d76d9;
  double d76d10;
  double d76d11;
  double d76d12;
  double d76d13;
  double d76d14;
  double d76d15;
  double d76d16;
  double d76d17;
  double d76d18;
  double d76d19;
  double d76d20;
  double d76d21;
  double d76d22;
  double d76d23;
  double d76d24;
  double d76d25;
  double d76d26;
  double d76d27;
  double d76d28;
  double d76d29;
  double d76d30;
  double d76d31;
  double d76d32;
  double d77d1;
  double d77d2;
  double d77d25;
  double d77d26;
  double d78d1;
  double d78d2;
  double d78d26;
  double d78d27;
  double d79d2;
  double d79d12;
  double d79d15;
  double d79d28;
  double d80d1;
  double d80d2;
  double d80d28;
  double d80d29;
  double d81d2;
  double d81d13;
  double d81d15;
  double d81d29;
  double d82d2;
  double d82d29;
  double d82d30;
  double d83d1;
  double d83d2;
  double d83d3;
  double d83d4;
  double d83d5;
  double d83d6;
  double d83d7;
  double d83d8;
  double d83d9;
  double d83d10;
  double d83d11;
  double d83d12;
  double d83d13;
  double d83d14;
  double d83d15;
  double d83d16;
  double d83d17;
  double d83d18;
  double d83d19;
  double d83d20;
  double d83d21;
  double d83d22;
  double d83d23;
  double d83d24;
  double d83d25;
  double d83d26;
  double d83d27;
  double d83d28;
  double d83d29;
  double d83d30;
  double d83d31;
  double d83d32;
  double d84d1;
  double d84d2;
  double d84d5;
  double d84d6;
  double d85d1;
  double d85d2;
  double d85d3;
  double d85d4;
  double d85d5;
  double d85d6;
  double d85d7;
  double d85d8;
  double d85d9;
  double d85d10;
  double d85d11;
  double d85d12;
  double d85d13;
  double d85d14;
  double d85d15;
  double d85d16;
  double d85d17;
  double d85d18;
  double d85d19;
  double d85d20;
  double d85d21;
  double d85d22;
  double d85d23;
  double d85d24;
  double d85d25;
  double d85d26;
  double d85d27;
  double d85d28;
  double d85d29;
  double d85d30;
  double d85d31;
  double d85d32;
  double d86d3;
  double d86d5;
  double d86d6;
  double d87d4;
  double d87d5;
  double d87d6;
  double d87d7;
  double d88d5;
  double d88d6;
  double d88d7;
  double d88d8;
  double d89d5;
  double d89d6;
  double d89d7;
  double d89d8;
  double d90d2;
  double d90d5;
  double d90d9;
  double d90d15;
  double d91d2;
  double d91d5;
  double d91d10;
  double d91d17;
  double d92d2;
  double d92d5;
  double d92d11;
  double d92d18;
  double d93d5;
  double d93d6;
  double d93d10;
  double d93d11;
  double d94d2;
  double d94d5;
  double d94d12;
  double d94d18;
  double d95d1;
  double d95d2;
  double d95d3;
  double d95d4;
  double d95d5;
  double d95d6;
  double d95d7;
  double d95d8;
  double d95d9;
  double d95d10;
  double d95d11;
  double d95d12;
  double d95d13;
  double d95d14;
  double d95d15;
  double d95d16;
  double d95d17;
  double d95d18;
  double d95d19;
  double d95d20;
  double d95d21;
  double d95d22;
  double d95d23;
  double d95d24;
  double d95d25;
  double d95d26;
  double d95d27;
  double d95d28;
  double d95d29;
  double d95d30;
  double d95d31;
  double d95d32;
  double d96d5;
  double d96d6;
  double d96d11;
  double d96d13;
  double d97d5;
  double d97d6;
  double d97d12;
  double d97d13;
  double d98d5;
  double d98d6;
  double d98d13;
  double d98d14;
  double d99d2;
  double d99d5;
  double d99d15;
  double d99d16;
  double d100d5;
  double d100d6;
  double d100d15;
  double d100d17;
  double d101d5;
  double d101d6;
  double d101d17;
  double d101d18;
  double d102d5;
  double d102d6;
  double d102d18;
  double d102d19;
  double d103d5;
  double d103d6;
  double d103d18;
  double d103d20;
  double d104d5;
  double d104d6;
  double d104d19;
  double d104d21;
  double d105d5;
  double d105d6;
  double d105d20;
  double d105d21;
  double d106d2;
  double d106d5;
  double d106d22;
  double d106d28;
  double d107d2;
  double d107d5;
  double d107d23;
  double d107d29;
  double d108d2;
  double d108d5;
  double d108d23;
  double d108d30;
  double d109d5;
  double d109d6;
  double d109d22;
  double d109d23;
  double d110d5;
  double d110d13;
  double d110d15;
  double d110d23;
  double d111d5;
  double d111d6;
  double d111d23;
  double d111d24;
  double d112d5;
  double d112d6;
  double d112d24;
  double d112d25;
  double d113d5;
  double d113d6;
  double d113d26;
  double d113d27;
  double d114d5;
  double d114d6;
  double d114d28;
  double d114d29;
  double d115d4;
  double d115d7;
  double d115d8;
  double d116d4;
  double d116d7;
  double d116d8;
  double d117d5;
  double d117d7;
  double d117d11;
  double d117d18;
  double d118d4;
  double d118d7;
  double d118d13;
  double d118d14;
  double d119d5;
  double d119d7;
  double d119d13;
  double d119d20;
  double d120d5;
  double d120d7;
  double d120d15;
  double d120d16;
  double d121d7;
  double d121d8;
  double d121d17;
  double d121d18;
  double d122d3;
  double d122d4;
  double d122d9;
  double d122d15;
  double d123d2;
  double d123d9;
  double d123d11;
  double d123d22;
  double d124d2;
  double d124d9;
  double d124d13;
  double d124d23;
  double d125d3;
  double d125d4;
  double d125d10;
  double d125d17;
  double d126d1;
  double d126d2;
  double d126d10;
  double d126d11;
  double d127d2;
  double d127d6;
  double d127d10;
  double d127d18;
  double d128d2;
  double d128d10;
  double d128d11;
  double d128d23;
  double d129d2;
  double d129d10;
  double d129d13;
  double d129d24;
  double d130d2;
  double d130d10;
  double d130d14;
  double d130d25;
  double d131d1;
  double d131d2;
  double d131d3;
  double d131d4;
  double d131d5;
  double d131d6;
  double d131d7;
  double d131d8;
  double d131d9;
  double d131d10;
  double d131d11;
  double d131d12;
  double d131d13;
  double d131d14;
  double d131d15;
  double d131d16;
  double d131d17;
  double d131d18;
  double d131d19;
  double d131d20;
  double d131d21;
  double d131d22;
  double d131d23;
  double d131d24;
  double d131d25;
  double d131d26;
  double d131d27;
  double d131d28;
  double d131d29;
  double d131d30;
  double d131d31;
  double d131d32;
  double d132d10;
  double d132d15;
  double d132d16;
  double d132d17;
  double d133d2;
  double d133d10;
  double d133d18;
  double d133d29;
  double d134d10;
  double d134d15;
  double d134d23;
  double d134d28;
  double d135d4;
  double d135d5;
  double d135d11;
  double d135d17;
  double d136d1;
  double d136d2;
  double d136d11;
  double d136d13;
  double d137d1;
  double d137d11;
  double d137d23;
  double d138d2;
  double d138d11;
  double d138d13;
  double d138d25;
  double d139d11;
  double d139d13;
  double d139d14;
  double d140d1;
  double d140d2;
  double d140d3;
  double d140d4;
  double d140d5;
  double d140d6;
  double d140d7;
  double d140d8;
  double d140d9;
  double d140d10;
  double d140d11;
  double d140d12;
  double d140d13;
  double d140d14;
  double d140d15;
  double d140d16;
  double d140d17;
  double d140d18;
  double d140d19;
  double d140d20;
  double d140d21;
  double d140d22;
  double d140d23;
  double d140d24;
  double d140d25;
  double d140d26;
  double d140d27;
  double d140d28;
  double d140d29;
  double d140d30;
  double d140d31;
  double d140d32;
  double d141d11;
  double d141d15;
  double d141d24;
  double d141d28;
  double d142d11;
  double d142d12;
  double d142d31;
  double d143d11;
  double d143d12;
  double d143d32;
  double d144d2;
  double d144d4;
  double d144d5;
  double d144d12;
  double d144d15;
  double d145d4;
  double d145d6;
  double d145d12;
  double d145d15;
  double d146d1;
  double d146d2;
  double d146d12;
  double d146d13;
  double d147d1;
  double d147d2;
  double d147d3;
  double d147d4;
  double d147d5;
  double d147d6;
  double d147d7;
  double d147d8;
  double d147d9;
  double d147d10;
  double d147d11;
  double d147d12;
  double d147d13;
  double d147d14;
  double d147d15;
  double d147d16;
  double d147d17;
  double d147d18;
  double d147d19;
  double d147d20;
  double d147d21;
  double d147d22;
  double d147d23;
  double d147d24;
  double d147d25;
  double d147d26;
  double d147d27;
  double d147d28;
  double d147d29;
  double d147d30;
  double d147d31;
  double d147d32;
  double d148d6;
  double d148d11;
  double d148d12;
  double d149d2;
  double d149d12;
  double d149d13;
  double d149d25;
  double d150d12;
  double d150d13;
  double d150d14;
  double d151d11;
  double d151d12;
  double d151d15;
  double d152d11;
  double d152d12;
  double d152d16;
  double d153d12;
  double d153d15;
  double d153d16;
  double d153d18;
  double d154d12;
  double d154d13;
  double d154d26;
  double d154d27;
  double d155d3;
  double d155d4;
  double d155d13;
  double d155d20;
  double d156d4;
  double d156d5;
  double d156d13;
  double d156d18;
  double d157d7;
  double d157d8;
  double d157d13;
  double d157d14;
  double d158d1;
  double d158d2;
  double d158d3;
  double d158d4;
  double d158d5;
  double d158d6;
  double d158d7;
  double d158d8;
  double d158d9;
  double d158d10;
  double d158d11;
  double d158d12;
  double d158d13;
  double d158d14;
  double d158d15;
  double d158d16;
  double d158d17;
  double d158d18;
  double d158d19;
  double d158d20;
  double d158d21;
  double d158d22;
  double d158d23;
  double d158d24;
  double d158d25;
  double d158d26;
  double d158d27;
  double d158d28;
  double d158d29;
  double d158d30;
  double d158d31;
  double d158d32;
  double d159d2;
  double d159d13;
  double d159d26;
  double d160d13;
  double d160d14;
  double d160d15;
  double d160d17;
  double d161d13;
  double d161d14;
  double d161d17;
  double d161d18;
  double d162d13;
  double d162d14;
  double d162d19;
  double d162d21;
  double d163d13;
  double d163d14;
  double d163d20;
  double d163d21;
  double d164d13;
  double d164d14;
  double d164d24;
  double d164d25;
  double d165d13;
  double d165d14;
  double d165d26;
  double d165d27;
  double d166d2;
  double d166d6;
  double d166d15;
  double d166d17;
  double d167d1;
  double d167d2;
  double d167d3;
  double d167d4;
  double d167d5;
  double d167d7;
  double d167d8;
  double d167d9;
  double d167d10;
  double d167d11;
  double d167d12;
  double d167d13;
  double d167d14;
  double d167d15;
  double d167d16;
  double d167d17;
  double d167d18;
  double d167d19;
  double d167d20;
  double d167d21;
  double d167d22;
  double d167d23;
  double d167d24;
  double d167d25;
  double d167d26;
  double d167d27;
  double d167d28;
  double d167d29;
  double d167d30;
  double d167d31;
  double d167d32;
  double d168d4;
  double d168d7;
  double d168d15;
  double d168d17;
  double d169d4;
  double d169d7;
  double d169d18;
  double d169d19;
  double d170d4;
  double d170d7;
  double d170d18;
  double d170d20;
  double d171d4;
  double d171d15;
  double d171d17;
  double d171d22;
  double d172d1;
  double d172d2;
  double d172d22;
  double d172d23;
  double d173d4;
  double d173d17;
  double d173d18;
  double d173d24;
  double d174d1;
  double d174d2;
  double d174d3;
  double d174d4;
  double d174d5;
  double d174d6;
  double d174d7;
  double d174d8;
  double d174d9;
  double d174d10;
  double d174d11;
  double d174d12;
  double d174d13;
  double d174d14;
  double d174d15;
  double d174d16;
  double d174d17;
  double d174d18;
  double d174d19;
  double d174d20;
  double d174d21;
  double d174d22;
  double d174d23;
  double d174d24;
  double d174d25;
  double d174d26;
  double d174d27;
  double d174d28;
  double d174d29;
  double d174d30;
  double d174d31;
  double d174d32;
  double d175d4;
  double d175d7;
  double d175d25;
  double d175d26;
  double d176d4;
  double d176d5;
  double d176d15;
  double d176d28;
  double d177d15;
  double d177d23;
  double d177d28;

};

#endif // !defined(OPENSMOKE_SYMBOLIC_KINETICS_GRI12)

