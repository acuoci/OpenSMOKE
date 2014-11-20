/***************************************************************************
 *   Copyright (C) 2009 by Alberto Cuoci								   *
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

#include <iomanip>
#include "basic/OpenSMOKE_Utilities.h"
#include "basic/OpenSMOKE_WarningFile.h"
#include "preprocessing_chemkin/OpenSMOKE_CHEMKINInterpreter_ElementsData.h"

void OpenSMOKE_CHEMKINInterpreter_ElementsData::ErrorMessage(const std::string message)
{
    cout << endl;
    cout << "Class:  OpenSMOKE_CHEMKINInterpreter_ElementsData"	<< endl;
    cout << "Object: " << name_object		<< endl;
    cout << "Error:  " << message           << endl;
    cout << "Press a key to continue... "   << endl;
    getchar();
    exit(-1);
}

void OpenSMOKE_CHEMKINInterpreter_ElementsData::WarningMessage(const std::string message)
{
    cout << endl;
    cout << "Class:  OpenSMOKE_CHEMKINInterpreter_ElementsData"	<< endl;
    cout << "Object: "		<< name_object	<< endl;
    cout << "Warning:  "	<< message      << endl;
    cout << "Press a key to continue... "   << endl;
    getchar();
}

void OpenSMOKE_CHEMKINInterpreter_ElementsData::Setup(ofstream *_fLog, OpenSMOKE_WarningFile *_fWarning)
{
	fLog		= _fLog;
	fWarning	= _fWarning; 
}

OpenSMOKE_CHEMKINInterpreter_ElementsData::OpenSMOKE_CHEMKINInterpreter_ElementsData()
{
	name_object = "Element parser";

	elements_name.push_back("Elements name");
	elements_description.push_back("Elements description");
	elements_mw.push_back(0);

	elements_name.push_back("H");
	elements_description.push_back("hydrogen");
	elements_mw.push_back(1.008000016212463);

	elements_name.push_back("C");
	elements_description.push_back("carbon");
	elements_mw.push_back(12.010999679565430);

	elements_name.push_back("O");
	elements_description.push_back("oxygen");
	elements_mw.push_back(15.998999595642090);

	elements_name.push_back("N");
	elements_description.push_back("nitrogen");
	elements_mw.push_back(14.0069999694824);

	elements_name.push_back("AR");
	elements_description.push_back("argon");
	elements_mw.push_back(39.948001861572270);

	elements_name.push_back("HE");
	elements_description.push_back("helium");
	elements_mw.push_back(4.002999782562256);

	elements_name.push_back("CL");
	elements_description.push_back("chloride");
	elements_mw.push_back(35.453);

	elements_name.push_back("LI");
	elements_description.push_back("litium");
	elements_mw.push_back(6.941);

	elements_name.push_back("BE");
	elements_description.push_back("berillium");
	elements_mw.push_back(9.01218);

	elements_name.push_back("B");
	elements_description.push_back("boron");
	elements_mw.push_back(10.81);

	elements_name.push_back("F");
	elements_description.push_back("fluorine");
	elements_mw.push_back(18.998403);
	
	elements_name.push_back("NE");
	elements_description.push_back("neon");
	elements_mw.push_back(20.179);

	elements_name.push_back("NA");
	elements_description.push_back("sodium");
	elements_mw.push_back(22.98977);

	elements_name.push_back("MG");
	elements_description.push_back("magnesium");
	elements_mw.push_back(24.305);

	elements_name.push_back("AL");
	elements_description.push_back("aluminium");
	elements_mw.push_back(26.98154);

	elements_name.push_back("SI");
	elements_description.push_back("silicium");
	elements_mw.push_back(28.0855);

	elements_name.push_back("P");
	elements_description.push_back("phosphorum");
	elements_mw.push_back(30.97376);

	elements_name.push_back("S");
	elements_description.push_back("solfurum");
	elements_mw.push_back(32.06);

	elements_name.push_back("K");
	elements_description.push_back("potassium");
	elements_mw.push_back(39.0983);

	elements_name.push_back("CA");
	elements_description.push_back("calcium");
	elements_mw.push_back(40.08);

	elements_name.push_back("SC");
	elements_description.push_back("scandium");
	elements_mw.push_back(44.9559);

	elements_name.push_back("TI");
	elements_description.push_back("titanium");
	elements_mw.push_back(47.90);

	elements_name.push_back("V");
	elements_description.push_back("vanadium");
	elements_mw.push_back(50.9415);

	elements_name.push_back("CR");
	elements_description.push_back("chromium");
	elements_mw.push_back(51.996);

	elements_name.push_back("MN");
	elements_description.push_back("manganese");
	elements_mw.push_back(54.9380);

	elements_name.push_back("FE");
	elements_description.push_back("ferrum");
	elements_mw.push_back(55.847);

	elements_name.push_back("CO");
	elements_description.push_back("cobaltum");
	elements_mw.push_back(58.9332);

	elements_name.push_back("NI");
	elements_description.push_back("nickel");
	elements_mw.push_back(58.70);

	elements_name.push_back("CU");
	elements_description.push_back("copper");
	elements_mw.push_back(63.546);

	elements_name.push_back("ZN");
	elements_description.push_back("zincum");
	elements_mw.push_back(65.38);

	elements_name.push_back("GA");
	elements_description.push_back("gallium");
	elements_mw.push_back(69.72);

	elements_name.push_back("GE");
	elements_description.push_back("germanium");
	elements_mw.push_back(72.59);

	elements_name.push_back("AS");
	elements_description.push_back("arsenicum");
	elements_mw.push_back(74.9216);

	elements_name.push_back("SE");
	elements_description.push_back("selenium");
	elements_mw.push_back(78.96);

	elements_name.push_back("BR");
	elements_description.push_back("bromum");
	elements_mw.push_back(79.904);

	elements_name.push_back("KR");
	elements_description.push_back("kripton");
	elements_mw.push_back(83.80);
 
	elements_name.push_back("RB");
	elements_description.push_back("rubidium");
	elements_mw.push_back(85.4678);

	elements_name.push_back("SR");
	elements_description.push_back("stronzium");
	elements_mw.push_back(87.62);

	elements_name.push_back("Y");
	elements_description.push_back("ittrium");
	elements_mw.push_back(88.9059);

	elements_name.push_back("ZR");
	elements_description.push_back("zirconium");
	elements_mw.push_back(91.22);

	elements_name.push_back("NB");
	elements_description.push_back("niobium");
	elements_mw.push_back(92.9064);

	elements_name.push_back("MO");
	elements_description.push_back("molibdenum");
	elements_mw.push_back(95.94);

	elements_name.push_back("TC");
	elements_description.push_back("tecnezium");
	elements_mw.push_back(98.);

	elements_name.push_back("RU");
	elements_description.push_back("rutenium");
	elements_mw.push_back(101.07);

	elements_name.push_back("RH");
	elements_description.push_back("rhodium");
	elements_mw.push_back(102.9055);

	elements_name.push_back("PD");
	elements_description.push_back("palladium");
	elements_mw.push_back(106.4);

	elements_name.push_back("AG");
	elements_description.push_back("silver");
	elements_mw.push_back(107.868);

	elements_name.push_back("CD");
	elements_description.push_back("cadmium");
	elements_mw.push_back(112.41);

	elements_name.push_back("IN");
	elements_description.push_back("indium");
	elements_mw.push_back(114.82);

	elements_name.push_back("SN");
	elements_description.push_back("stagnum");
	elements_mw.push_back(118.69);

	elements_name.push_back("SB");
	elements_description.push_back("antimonium");
	elements_mw.push_back(121.75);

	elements_name.push_back("TE");
	elements_description.push_back("tellurium");
	elements_mw.push_back(127.60);

	elements_name.push_back("I");
	elements_description.push_back("iodium");
	elements_mw.push_back(126.9045);

	elements_name.push_back("XE");
	elements_description.push_back("xenon");
	elements_mw.push_back(131.30);

	elements_name.push_back("CS");
	elements_description.push_back("cesium");
	elements_mw.push_back(132.9054);

	elements_name.push_back("BA");
	elements_description.push_back("barium");
	elements_mw.push_back(137.33);

	elements_name.push_back("LA");
	elements_description.push_back("lantanium");
	elements_mw.push_back(138.9055);

	elements_name.push_back("CE");
	elements_description.push_back("cerium");
	elements_mw.push_back(140.12);

	elements_name.push_back("PR");
	elements_description.push_back("praseodimium");
	elements_mw.push_back(140.9077);

	elements_name.push_back("ND");
	elements_description.push_back("neodimium");
	elements_mw.push_back(144.24);

	elements_name.push_back("PM");
	elements_description.push_back("promezium");
	elements_mw.push_back(145.);

	elements_name.push_back("SM");
	elements_description.push_back("samarium");
	elements_mw.push_back(150.4);

	elements_name.push_back("EU");
	elements_description.push_back("europium");
	elements_mw.push_back(151.96);

	elements_name.push_back("GD");
	elements_description.push_back("gadolinium");
	elements_mw.push_back(157.25);

	elements_name.push_back("TB");
	elements_description.push_back("terbium");
	elements_mw.push_back(158.9254);

	elements_name.push_back("DY");
	elements_description.push_back("dysprosium");
	elements_mw.push_back(162.50);

	elements_name.push_back("HO");
	elements_description.push_back("olmium");
	elements_mw.push_back(164.9304);

	elements_name.push_back("ER");
	elements_description.push_back("erbium");
	elements_mw.push_back(167.26);

	elements_name.push_back("TM");
	elements_description.push_back("tulium");
	elements_mw.push_back(168.9342);

	elements_name.push_back("YB");
	elements_description.push_back("itterbium");
	elements_mw.push_back(173.04);

	elements_name.push_back("LU");
	elements_description.push_back("lutezium");
	elements_mw.push_back(174.967);

	elements_name.push_back("HF");
	elements_description.push_back("afnium");
	elements_mw.push_back(178.49);

	elements_name.push_back("W");
	elements_description.push_back("tungstenum");
	elements_mw.push_back(183.85);

	elements_name.push_back("RE");
	elements_description.push_back("renium");
	elements_mw.push_back(186.207);

	elements_name.push_back("OS");
	elements_description.push_back("osmium");
	elements_mw.push_back(190.2);

	elements_name.push_back("IR");
	elements_description.push_back("iridium");
	elements_mw.push_back(192.22);

	elements_name.push_back("PT");
	elements_description.push_back("platinum");
	elements_mw.push_back(195.09);

	elements_name.push_back("AU");
	elements_description.push_back("gold");
	elements_mw.push_back(196.9665);

	elements_name.push_back("HG");
	elements_description.push_back("mercurium");
	elements_mw.push_back(200.59);

	elements_name.push_back("TL");
	elements_description.push_back("tallium");
	elements_mw.push_back(204.37);

	elements_name.push_back("PB");
	elements_description.push_back("lead");
	elements_mw.push_back(207.2);

	elements_name.push_back("BI");
	elements_description.push_back("bismutum");
	elements_mw.push_back(208.9804);

	elements_name.push_back("PO");
	elements_description.push_back("polonium");
	elements_mw.push_back(209.);

	elements_name.push_back("AT");
	elements_description.push_back("astatum");
	elements_mw.push_back(210.);
 
	elements_name.push_back("RN");
	elements_description.push_back("radon");
	elements_mw.push_back(222.);

	elements_name.push_back("FR");
	elements_description.push_back("francium");
	elements_mw.push_back(223.);

	elements_name.push_back("RA");
	elements_description.push_back("radium");
	elements_mw.push_back(226.0254);

	elements_name.push_back("AC");
	elements_description.push_back("attinium");
	elements_mw.push_back(227.0278);

	elements_name.push_back("TH");
	elements_description.push_back("thorium");
	elements_mw.push_back(232.0381);

	elements_name.push_back("PA");
	elements_description.push_back("protoattinium");
	elements_mw.push_back(231.0359);

	elements_name.push_back("U");
	elements_description.push_back("uranium");
	elements_mw.push_back(238.029);

	elements_name.push_back("NP");
	elements_description.push_back("nettunium");
	elements_mw.push_back(237.0482);

	elements_name.push_back("PU");
	elements_description.push_back("plutonium");
	elements_mw.push_back(244.);

	elements_name.push_back("AM");
	elements_description.push_back("americium");
	elements_mw.push_back(243.);

	elements_name.push_back("CM");
	elements_description.push_back("curium");
	elements_mw.push_back(247.);

	elements_name.push_back("BK");
	elements_description.push_back("berkelium");
	elements_mw.push_back(247.);

	elements_name.push_back("CF");
	elements_description.push_back("californium");
	elements_mw.push_back(251.);

	elements_name.push_back("ES");
	elements_description.push_back("einstenium");
	elements_mw.push_back(254.);

	elements_name.push_back("FM");
	elements_description.push_back("fermium");
	elements_mw.push_back(257.);

	elements_name.push_back("MD");
	elements_description.push_back("mendelenium");
	elements_mw.push_back(258.);

	elements_name.push_back("NO");
	elements_description.push_back("nobelium");
	elements_mw.push_back(259.);

	elements_name.push_back("LR");
	elements_description.push_back("laurenzium");
	elements_mw.push_back(260.);

	elements_name.push_back("D");
	elements_description.push_back("deuterium");
	elements_mw.push_back(2.0160000324248);

	elements_name.push_back("E");
	elements_description.push_back("electron");
	elements_mw.push_back(0.0005486579);


	elements_in_database = elements_mw.size()-1;
	elements_in_list.push_back("actual list");
	isotope_in_list.push_back("actual list");
	isotope_mw_in_list.push_back(0);
	elements_in_list_indices.push_back(0);

}

bool OpenSMOKE_CHEMKINInterpreter_ElementsData::Parse_Element_Name(const std::string element_name)
{
	for(int i=1;i<=elements_in_database;i++)
		if ( caseInsCompare(element_name, elements_name[i]) )
		{
			for(int j=1;j<=int(elements_in_list.size())-1;j++)
				if (element_name == elements_in_list[j])
					ErrorMessage("The element " + element_name + " was specified more than one time!");
			elements_in_list.push_back(elements_name[i]);
			elements_in_list_indices.push_back(i);
			return true;
		}

	return false;
}

void OpenSMOKE_CHEMKINInterpreter_ElementsData::Parse_Isotope_Name(const std::string isotope_name, const double mw)
{
	int i;

	for (i=1;i<=int(isotope_in_list.size())-1;i++)
		if (caseInsCompare(isotope_name,isotope_in_list[i])	== true)
			ErrorMessage("The " + isotope_name + " isotope was declared more than once");

	for (i=1;i<=int(elements_name.size())-1;i++)
		if (caseInsCompare(isotope_name,elements_name[i])	== true)
			ErrorMessage("The " + isotope_name + " isotope is already included in the database. Please choose a different name");

	isotope_in_list.push_back(isotope_name);
	isotope_mw_in_list.push_back(mw);

	if (mw <= 0.)	
		ErrorMessage("The molecular weigth of " + isotope_name + " must be positive");
}

void OpenSMOKE_CHEMKINInterpreter_ElementsData::Summary()
{
	int i;

	*fLog << " ----------------------------------------------------------------" << endl;
	*fLog << "                            Elements                             " << endl;
	*fLog << " ----------------------------------------------------------------" << endl;
	for(i=1;i<=int(elements_in_list.size())-1;i++)
		*fLog << elements_name[elements_in_list_indices[i]] << "\t" << elements_mw[elements_in_list_indices[i]] << endl;
	*fLog  << endl;

	*fLog << " ----------------------------------------------------------------" << endl;
	*fLog << "                            Isotopes                             " << endl;
	*fLog << " ----------------------------------------------------------------" << endl;
	for(i=1;i<=int(isotope_in_list.size())-1;i++)
		*fLog  << isotope_in_list[i] << "\t" << isotope_mw_in_list[i] << endl;
	*fLog  << endl;
}

void OpenSMOKE_CHEMKINInterpreter_ElementsData::SummaryOnFile(ofstream &fOutput)
{
	int i;
	int k=1;

	fOutput << "------------------------------"	<< endl;
	fOutput << "  Elements     Atomic Weight  "	<< endl;
	fOutput << "------------------------------"	<< endl;
	for(i=1;i<=int(elements_in_list.size())-1;i++)
	{
		fOutput << "  " << k++ << ".\t";
 		fOutput << setw(9) << left << elements_name[elements_in_list_indices[i]];
		fOutput << fixed   << setw(12) << setprecision(8) << right << elements_mw[elements_in_list_indices[i]];
		fOutput << endl;
	}

	for(i=1;i<=int(isotope_in_list.size())-1;i++)
	{
		fOutput << "  " << k++ << ".\t";
 		fOutput << setw(9) << left << isotope_in_list[i];
		fOutput << fixed   << setw(12) << setprecision(8) << right << isotope_mw_in_list[i];
		fOutput << endl;
	}

	fOutput << "------------------------------"	<< endl;
	fOutput << endl << endl;
}

