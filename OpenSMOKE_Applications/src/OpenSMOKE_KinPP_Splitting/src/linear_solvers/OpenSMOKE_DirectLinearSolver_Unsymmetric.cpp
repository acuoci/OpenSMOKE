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

#include <iomanip>
#include "OpenSMOKE_DirectLinearSolver_Unsymmetric.h"

void OpenSMOKE_DirectLinearSolver_Unsymmetric::ErrorMessage(const string message_)
{
    cout << endl;
    cout << "Class: " << name_solver_			<< endl;
    cout << "Error: " << message_				<< endl;
    cout << "Press enter to continue... "		<< endl;
    getchar();
    exit(-1);
}

void OpenSMOKE_DirectLinearSolver_Unsymmetric::WarningMessage(const string message_)
{
    cout << endl;
    cout << "Class:	  " << name_solver_	<< endl;
    cout << "Warning: " << message_		<< endl;
	cout << endl;
}

void OpenSMOKE_DirectLinearSolver_Unsymmetric::MessageOnTheScreen(const string message_)
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

void OpenSMOKE_DirectLinearSolver_Unsymmetric::CalculateSparsityPattern(const int index, const int blockDimension, const BzzVector& mConvectionDiffusion, const BzzVectorInt& iConvectionDiffusion, BzzVectorInt &indicesColumns)
{
	// Global columns for row block corresponding to the current reactor)	
	ChangeDimensions (blockDimension * (blockDimension+iConvectionDiffusion.Size()), &indicesColumns);
	{
		int jReactor = (index-1)*blockDimension;

		int count=1;
		for (int i=1;i<=blockDimension;i++)
		{
			for (int k=1;k<=mConvectionDiffusion.Size();k++)
				if (iConvectionDiffusion[k] < int(index))
					indicesColumns[count++] = (iConvectionDiffusion[k]-1)*blockDimension+i;

			for (int k=1;k<=blockDimension;k++)
				indicesColumns[count++] = jReactor+k;

			for (int k=1;k<=mConvectionDiffusion.Size();k++)
				if (iConvectionDiffusion[k] > int(index))
					indicesColumns[count++] = (iConvectionDiffusion[k]-1)*blockDimension+i;

			for (int k=1;k<=mConvectionDiffusion.Size();k++)
				if (iConvectionDiffusion[k] == int(index))
					ErrorMessage("Wrong sparsity structure");
		}
	}
}
