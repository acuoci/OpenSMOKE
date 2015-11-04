/***************************************************************************
 *   Copyright (C) 2011 by Alberto Cuoci								   *
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

#include <string>
#include <iomanip>
#include "kinpp/OpenSMOKE_DirectLinearSolver_Unsymmetric.h"
using namespace std;

void OpenSMOKE_DirectLinearSolver_Unsymmetric::ErrorMessage(const std::string message_)
{
    cout << endl;
    cout << "Class: " << name_solver_			<< endl;
    cout << "Error: " << message_				<< endl;
    cout << "Press enter to continue... "		<< endl;
    getchar();
    exit(-1);
}

void OpenSMOKE_DirectLinearSolver_Unsymmetric::WarningMessage(const std::string message_)
{
    cout << endl;
    cout << "Class:	  " << name_solver_	<< endl;
    cout << "Warning: " << message_		<< endl;
	cout << endl;
}

void OpenSMOKE_DirectLinearSolver_Unsymmetric::MessageOnTheScreen(const std::string message_)
{
	cout << endl;
	cout << "|| ********************************************************************************** ||" << endl;
	cout << "||   ";
	cout << message_ << endl;
	cout << "|| ********************************************************************************** ||" << endl;
	cout << endl;
}

void OpenSMOKE_DirectLinearSolver_Unsymmetric::CleanMemory()
{
	delete[] rows;
	delete[] columns;
	delete[] values;
}