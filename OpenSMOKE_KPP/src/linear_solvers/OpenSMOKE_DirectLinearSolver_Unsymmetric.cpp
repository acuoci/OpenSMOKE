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
#include "../kpp/OpenSMOKE_KPP_ReactorNetwork.h"
#include <vector>

void OpenSMOKE_DirectLinearSolver_Unsymmetric::ErrorMessage(const std::string message_)
{
    std::cout << std::endl;
    std::cout << "Class: " << name_solver_			<< std::endl;
    std::cout << "Error: " << message_				<< std::endl;
    std::cout << "Press enter to continue... "		<< std::endl;
    getchar();
    exit(-1);
}

void OpenSMOKE_DirectLinearSolver_Unsymmetric::WarningMessage(const std::string message_)
{
    std::cout << std::endl;
    std::cout << "Class:	  " << name_solver_	<< std::endl;
    std::cout << "Warning: " << message_		<< std::endl;
    std::cout << std::endl;
}

void OpenSMOKE_DirectLinearSolver_Unsymmetric::MessageOnTheScreen(const std::string message_)
{
	std::cout << std::endl;
	std::cout << "|| ********************************************************************************** ||" << std::endl;
	std::cout << "||   ";
	std::cout << message_ << std::endl;
	std::cout << "|| ********************************************************************************** ||" << std::endl;
	std::cout << std::endl;
}

void OpenSMOKE_DirectLinearSolver_Unsymmetric::CleanMemory()
{
	delete [] rows;
	delete [] columns;
	delete [] values;
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

void OpenSMOKE_DirectLinearSolver_Unsymmetric::CalculateSparsityPattern(const int index, const int blockDimension, const std::vector<double*>& mConvDiff, const std::vector<int*>& iConvDiff, BzzVectorInt &indicesColumns, OpenSMOKE_KPP_ReactorNetwork& network)
{

	int ConvDiffSize = network.CD_Size()[index];
	// Global columns for row block corresponding to the current reactor)	
	ChangeDimensions (blockDimension * (blockDimension+ConvDiffSize), &indicesColumns);

	{
		int jReactor = (index-1)*blockDimension;

		int count=1;
		for (int i=1;i<=blockDimension;i++)
		{
			for (int k=1;k<=ConvDiffSize;k++)
				if (iConvDiff[index][k] < int(index))
					indicesColumns[count++] = (iConvDiff[index][k]-1)*blockDimension+i;

			for (int k=1;k<=blockDimension;k++)
				indicesColumns[count++] = jReactor+k;

			for (int k=1;k<=ConvDiffSize;k++)
				if (iConvDiff[index][k] > int(index))
					indicesColumns[count++] = (iConvDiff[index][k]-1)*blockDimension+i;

			for (int k=1;k<=ConvDiffSize;k++)
				if (iConvDiff[index][k] == int(index))
					ErrorMessage("Wrong sparsity structure");
		}
	}
}

void OpenSMOKE_DirectLinearSolver_Unsymmetric::CalculateSparsityPattern(const int index, const int blockDimension, const std::vector<double*>& mConvDiff, const std::vector<int*>& iConvDiff, LIS_INT* &indicesColumns, OpenSMOKE_KPP_ReactorNetwork& network)
{

}
