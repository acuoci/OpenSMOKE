/***************************************************************************
 *   Copyright (C) 2006-2009 by Alberto Cuoci						       *
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

#include "OpenSMOKE_LookUp_Table_Flame.h"
#include "basic/OpenSMOKE_Utilities.h"

void OpenSMOKE_LookUp_Table_Flame::ErrorMessage(const string message)
{
    cout << endl;
    cout << "Class:  OpenSMOKE_LookUp_Table_Flame"	<< endl;
    cout << "Object: " << name_object		<< endl;
    cout << "Error:  " << message           << endl;
    cout << "Press a key to continue... "   << endl;
    getchar();
    exit(-1);
}

void OpenSMOKE_LookUp_Table_Flame::WarningMessage(const string message)
{
    cout << endl;
    cout << "Class:  OpenSMOKE_LookUp_Table_Flame"	<< endl;
    cout << "Object: "		<< name_object	<< endl;
    cout << "Warning:  "	<< message      << endl;
    cout << "Press a key to continue... "   << endl;
    getchar();
}

void OpenSMOKE_LookUp_Table_Flame::read_from_file(const string fileName)
{
	int i,j;
	string dummy;

	ifstream fInput;
	openInputFileAndControl(fInput, fileName);
	
	fInput >> dummy;
	if (dummy != "mf")	ErrorMessage("Expected: mf - Found: " + dummy);
	fInput >> nCsi;
	ChangeDimensions(nCsi, &csi);
	for(j=1;j<=nCsi;j++)
		fInput >> csi[j];

	fInput >> dummy;
	if (dummy != "mfv")	ErrorMessage("Expected: mfv - Found: " + dummy);
	fInput >> nCsiV;
	ChangeDimensions(nCsiV, &csiV);
	for(j=1;j<=nCsiV;j++)
		fInput >> csiV[j];

	fInput >> dummy;
	if (dummy != "temperature")	ErrorMessage("Expected: temperature - Found: " + dummy);
	ChangeDimensions(nCsi, nCsiV, &temperature);
	for(j=1;j<=nCsiV;j++)
		for(i=1;i<=nCsi;i++)
			fInput >> temperature[i][j];

	fInput >> dummy;
	if (dummy != "density")	ErrorMessage("Expected: density - Found: " + dummy);
	ChangeDimensions(nCsi, nCsiV, &density);
	for(j=1;j<=nCsiV;j++)
		for(i=1;i<=nCsi;i++)
			fInput >> density[i][j];

	fInput >> dummy;
	if (dummy != "cp")	ErrorMessage("Expected: cp - Found: " + dummy);
	ChangeDimensions(nCsi, nCsiV, &cp);
	for(j=1;j<=nCsiV;j++)
		for(i=1;i<=nCsi;i++)
			fInput >> cp[i][j];

	fInput >> dummy;
	if (dummy != "as")	ErrorMessage("Expected: as - Found: " + dummy);
	ChangeDimensions(nCsi, nCsiV, &as);
	for(j=1;j<=nCsiV;j++)
		for(i=1;i<=nCsi;i++)
			fInput >> as[i][j];

	fInput.close();
}
