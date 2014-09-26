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
#include "preprocessing_chemkin/OpenSMOKE_CHEMKINInterpreter_SpeciesData.h"

void OpenSMOKE_CHEMKINInterpreter_SpeciesData::ErrorMessage(const std::string message)
{
    cout << endl;
    cout << "Class:  OpenSMOKE_CHEMKINInterpreter_ElementsData"	<< endl;
    cout << "Object: " << name_object		<< endl;
    cout << "Error:  " << message           << endl;
    cout << "Press a key to continue... "   << endl;
    getchar();
    exit(-1);
}

void OpenSMOKE_CHEMKINInterpreter_SpeciesData::WarningMessage(const std::string message)
{
    cout << endl;
    cout << "Class:  OpenSMOKE_CHEMKINInterpreter_ElementsData"	<< endl;
    cout << "Object: "		<< name_object	<< endl;
    cout << "Warning:  "	<< message      << endl;
    cout << "Press a key to continue... "   << endl;
    getchar();
}

void OpenSMOKE_CHEMKINInterpreter_SpeciesData::Setup(ofstream *_fLog, OpenSMOKE_WarningFile *_fWarning)
{
	fLog		= _fLog;
	fWarning	= _fWarning; 
}

OpenSMOKE_CHEMKINInterpreter_SpeciesData::OpenSMOKE_CHEMKINInterpreter_SpeciesData()
{
	name_object = "Species parser";
	species_in_list.push_back("actual list");
}

bool OpenSMOKE_CHEMKINInterpreter_SpeciesData::Parse_Species_Name(const std::string species_name)
{
	for(int j=1;j<=int(species_in_list.size())-1;j++)
		if (species_name == species_in_list[j])
			ErrorMessage("The species " + species_name + " was specified more than one time!");

	species_in_list.push_back(species_name);
	return true;
}

void OpenSMOKE_CHEMKINInterpreter_SpeciesData::Summary()
{
	*fLog << " ----------------------------------------------------------------" << endl;
	*fLog << "                     Species List                                " << endl;
	*fLog << " ----------------------------------------------------------------" << endl;
	for(int i=1;i<=int(species_in_list.size())-1;i++)
		*fLog << setw(6) << right << i << setw(20) << left << species_in_list[i] << endl;
	*fLog << endl << endl;
}
