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

#include "idealreactors/flame1d/OpenSMOKE_Flame1D_ScheduleClass.h"

const int SparkOn  = 1;
const int SparkOff = 0;

void OpenSMOKE_Flame1D_ScheduleClass::readOperations(const std::string fileName)
{
	const int MAXOPERATIONS=100;
	int    i, k;
	int    nCycle, onCycle, startCycle, endCycle;
	std::string stringa;
	std::string stringaOptional;
	char   charIndex;
	double t;
	species = new std::string[MAXOPERATIONS];
	std::string option;

	ChangeDimensions(0, &iOperation);
	ChangeDimensions(0, &iOptionA);
	ChangeDimensions(0, &iOptionB);

	// LEGENDA
	// ------------------------------------------
	//
	//	0 = INITIALCOLD				---						---
	//	1 = NLS				1 = ONLY_MOMENTUM					---		
	//						2 = ALL
	//	2 = DAESTART				---						time
	//	3 = DAE						---						time
	//	4 = NEWPOINTS		1 = TEMPERATURE					---


	ifstream fInput;
	openInputFileAndControl(fInput, fileName);

	for(;;)
	{	fInput >> stringa;
	
		// Recupero BackUp
		// -----------------------------------------------
		if (stringa == "RECOVER_FROM_BACK_UP")
		{
			iOperation.Append(-1);
			fInput >> stringa;

			option = stringa;
			iOptionA.Append(0);
			iOptionB.Append(0);
		}

		// Soluzione di partenza a freddo
		// -----------------------------------------------
		//else if (stringa == "START_FROM_ZERO")
		//{
		//	iOperation.Append(0);
		//	fInput >> t;
		//	iOptionA.Append(0);
		//	iOptionB.Append(t);
		//}
	
		// Soluzione del sistema non lineare
		// -----------------------------------------------
		else if (stringa == "NLS")
		{
			iOperation.Append(1);
			fInput >> stringa;
			if (stringa == "ONLY_MOMENTUM")
				{ iOptionA.Append(1); iOptionB.Append(0.);}
			else if (stringa == "ALL")
				{ iOptionA.Append(2); iOptionB.Append(0.);}
			else if (stringa == "ONLY_MASS_FRACTIONS")
				{ iOptionA.Append(3); iOptionB.Append(0.);}
			else if (stringa == "FLAMESPEED")
				{ iOptionA.Append(4); iOptionB.Append(0.);}
			else ErrorMessage(stringa);
		}
		
		// Soluzione del sistema non lineare
		// -----------------------------------------------
		else if (stringa == "NLS_ALL")
		{
			iOperation.Append(111);
			fInput >> t;
			iOptionA.Append(0); 
			iOptionB.Append(t);
		}

		// Soluzione del sistema non lineare
		// -----------------------------------------------
		else if (stringa == "SENSITIVITY")
		{
			iOperation.Append(6);
			fInput >> stringa;
			if (stringa == "FREQUENCY_FACTOR")
				{ iOptionA.Append(1); iOptionB.Append(0.);}
			else if (stringa == "TRANSPORT_PROPERTIES")
				{ iOptionA.Append(2); iOptionB.Append(0.);}
			else if (stringa == "FREQUENCY_FACTOR_AND_TRANSPORT_PROPERTIES")
				{ iOptionA.Append(3); iOptionB.Append(0.);}
			else ErrorMessage(stringa);
		}

		// Soluzione del sistema differenziale
		// ---------------------------------------------
		else if (stringa == "ODE")
		{
			iOperation.Append(2);
			fInput >> t;
		
			iOptionA.Append(SparkOff);
			iOptionB.Append(t);
		}

		// Soluzione del sistema differenziale-algebrico
		// ---------------------------------------------
		else if (stringa == "DAE")
		{
			iOperation.Append(3);
			fInput >> t;
		
			iOptionA.Append(SparkOff);
			iOptionB.Append(t);
		}

		// Soluzione del sistema differenziale-algebrico
		// ---------------------------------------------
		else if (stringa == "ODESINGLEREACTOR")
		{
			iOperation.Append(61);
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
				iOperation.Append(41);
			else if (stringaOptional == "QREACTION")
				iOperation.Append(42);
			else if (stringaOptional == "ALL")
				iOperation.Append(43);
			else 
			{
				iOperation.Append(44);
				species[iOperation.Size()] = stringaOptional;
			}


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
			else if (charIndex == 'C')
			{
				iOptionA.Append('C');
				iOptionB.Append(0.);
			}

			else ErrorMessage(stringa);
		}

		// Raffinamento griglia
		// ---------------------------------------------
		else if (stringa == "DOUBLEGRID")
		{
			iOperation.Append(51);
			iOptionA.Append(0);
			iOptionB.Append(0.);
		}

		// Raffinamento griglia
		// ---------------------------------------------
		else if (stringa == "REFINEPEAK")
		{
			iOperation.Append(52);
			iOptionA.Append(0);
			iOptionB.Append(0.);
		}

		else if (stringa == "REFINE")
		{
			double xA;
			double xB;
			fInput >> xA;
			fInput >> xB;
			if (xA>=xB)
			{
				cout << "Error in Operations.inp: xA >= xB!" << endl;
				cout << "Press enter to continue..." << endl;
				getchar();
				exit(-1);
			}	
			iOperation.Append(54);
			iOptionA.Append(xA);
			iOptionB.Append(xB);
		}

		else if (stringa == "REFINE_STAGNATION_PLANE")
		{
			double xA;
			fInput >> xA;
			iOperation.Append(58);
			iOptionA.Append(xA);
			iOptionB.Append(0.);
		}

		else if (stringa == "ADAPT")
		{
			int xA;
			fInput >> xA;
			iOperation.Append(55);
			iOptionA.Append(xA);
			iOptionB.Append(0.);
		}

		// Raffinamento griglia
		// ---------------------------------------------
		else if (stringa == "REFINE_FLAME_BASE")
		{
			iOperation.Append(53);
			iOptionA.Append(0);
			iOptionB.Append(0.);
		}

				// Raffinamento griglia
		// ---------------------------------------------
		else if (stringa == "REFINE_BASE")
		{
			iOperation.Append(56);
			iOptionA.Append(0);
			iOptionB.Append(0.);
		}

		else if (stringa == "REFINE_ATTACH_POINT")
		{
			iOperation.Append(57);
			iOptionA.Append(0);
			iOptionB.Append(0.);
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
			if (onCycle == 0) ErrorMessage("The cycle must is not open");
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

		// Errore
		// ---------------------------------------------
		else 
		{
			ErrorMessage(stringa);
		}
	}

	// Riepilogo
	nOperations = iOperation.Size();
}


void OpenSMOKE_Flame1D_ScheduleClass::ErrorMessage(const std::string message)
{
    cout << endl;
    cout << "Class:  OpenSMOKE_Flame1D_ScheduleClass"	<< endl;
    cout << "Object: " << name_object					<< endl;
    cout << "Error:  " << message						<< endl;
    cout << "Press a key to continue... "				<< endl;
    getchar();
    exit(-1);
}

void OpenSMOKE_Flame1D_ScheduleClass::WarningMessage(const std::string message)
{
    cout << endl;
    cout << "Class:  OpenSMOKE_Flame1D_ScheduleClass"	<< endl;
    cout << "Object: "		<< name_object				<< endl;
    cout << "Warning:  "	<< message					<< endl;
    cout << "Press a key to continue... "				<< endl;
    getchar();
}

OpenSMOKE_Flame1D_ScheduleClass::OpenSMOKE_Flame1D_ScheduleClass()
{
	name_object = "[undefined name]";
}

void OpenSMOKE_Flame1D_ScheduleClass::SetName(const std::string name)
{
	name_object = name;
}