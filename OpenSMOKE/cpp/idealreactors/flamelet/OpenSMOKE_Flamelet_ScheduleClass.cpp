/***************************************************************************
 *   Copyright (C) 2006 by Alberto Cuoci								   *
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

#include "basic/OpenSMOKE_Utilities.h"
#include "idealreactors/flamelet/OpenSMOKE_Flamelet_ScheduleClass.h"

void OpenSMOKE_Flamelet_ScheduleClass::ErrorMessage(const std::string message)
{
    cout << endl;
    cout << "Class:  ScheduleClass"			<< endl;
    cout << "Object: " << name_object		<< endl;
    cout << "Error:  " << message           << endl;
    cout << "Press a key to continue... "   << endl;
    getchar();
    exit(-1);
}

void OpenSMOKE_Flamelet_ScheduleClass::WarningMessage(const std::string message)
{
    cout << endl;
    cout << "Class:  ScheduleClass"			<< endl;
    cout << "Object: "		<< name_object	<< endl;
    cout << "Warning:  "	<< message      << endl;
    cout << "Press a key to continue... "   << endl;
    getchar();
}

OpenSMOKE_Flamelet_ScheduleClass::OpenSMOKE_Flamelet_ScheduleClass()
{
	name_object = "Undefined name";
}

void OpenSMOKE_Flamelet_ScheduleClass::SetName(const std::string name)
{
	name_object = name;
}

void OpenSMOKE_Flamelet_ScheduleClass::ReadOperations(const std::string fileName)
{
	int i, k, nCycle, onCycle, startCycle, endCycle;
	std::string stringa;
	std::string stringaOptional;
	double t;
	char charIndex;

	ifstream fInput;
	openInputFileAndControl(fInput, fileName);

	for(;;)
	{	
		fInput >> stringa;

		// Soluzione del sistema differenziale
		// ---------------------------------------------
		if (stringa == "SOLVE")
		{
			iOperation.Append(1);
			fInput >> t;
		
			iOptionA.Append(0);
			iOptionB.Append(t);
		}

		// Aggiunta di punti
		// ---------------------------------------------
		else if (stringa == "NEWPOINTS")
		{
			fInput >> stringaOptional;
			if (stringaOptional == "TEMPERATURE")
				iOperation.Append(111);
			else if (stringaOptional == "QREACTION")
				iOperation.Append(112);
			else
				ErrorMessage("You can adapt the grid only on the temperature or Qreaction profile!");

			fInput >> charIndex;
			if (charIndex == 'D')
			{	
				iOptionA.Append('D'); 
				iOptionB.Append(0.);
			}
			else if (charIndex == 'G')

			{
				iOptionA.Append('G'); 
				iOptionB.Append(0.);
			}

			else ErrorMessage(stringa);
		}

		// Raffinamento griglia
		// ---------------------------------------------
		else if (stringa == "DOUBLEGRID")
		{
			iOperation.Append(12);
			iOptionA.Append(0);
			iOptionB.Append(0.);
		}

		// Raffinamento griglia
		// ---------------------------------------------
		else if (stringa == "REFINEPEAK")
		{
			double double_option;
			fInput >> double_option;
			iOperation.Append(13);
			iOptionA.Append(0);
			iOptionB.Append(double_option);
		}

		// Raffinamento griglia
		// ---------------------------------------------
		else if (stringa == "REFINELEANSIDE")
		{
			double double_option;
			fInput >> double_option;
			iOperation.Append(14);
			iOptionA.Append(0);
			iOptionB.Append(double_option);
		}

		// Dichiarazione di un ciclo
		// ---------------------------------------------
		else if (stringa == "CYCLE")
		{
			fInput >> nCycle;
			onCycle = 1;
			startCycle = iOperation.Size();
		}

		// Dichiarazione di un ciclo
		// ---------------------------------------------
		else if (stringa == "ENDCYCLE")
		{
			if (onCycle == 0) ErrorMessage("Cycle assignement is not correct!");
			onCycle = 0;
			endCycle = iOperation.Size();

			for (i=1;i<=nCycle-1;i++)
			for (k=1;k<=endCycle-startCycle;k++)
			{
				iOperation.Append(iOperation[startCycle+k]);
				iOptionA.Append(iOptionA[startCycle+k]);
				iOptionB.Append(iOptionB[startCycle+k]);
			}
		}

		// Dichiarazione di un ciclo
		// ---------------------------------------------
		else if (stringa == "END") break;

		// Error
		// ---------------------------------------------
		else 
			ErrorMessage("Wrong key word: " + stringa);
	}

	// Summary
	nOperations = iOperation.Size();
}
