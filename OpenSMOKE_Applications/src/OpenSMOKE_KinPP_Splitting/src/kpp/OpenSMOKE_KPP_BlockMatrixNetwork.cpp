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

#include "OpenSMOKE_KPP_DataManager.h"
#include "OpenSMOKE_KPP_Dictionary.h"
#include "OpenSMOKE_KPP_BlockMatrixNetwork.h"
#include <iostream>
#include <iomanip>

void OpenSMOKE_KPP_BlockMatrixNetwork::ErrorMessage(const string message_)
{
    cout << endl;
    cout << "Class: OpenSMOKE_KPP_BlockMatrixNetwork"	<< endl;
    cout << "Error: " << message_					<< endl;
    cout << "Press enter to continue... "			<< endl;
    getchar();
    exit(-1);
}

void OpenSMOKE_KPP_BlockMatrixNetwork::WarningMessage(const string message_)
{
    cout << endl;
    cout << "Class:	  OpenSMOKE_KPP_BlockMatrixNetwork"	<< endl;
    cout << "Warning: " << message_					<< endl;
	cout << endl;
}

OpenSMOKE_KPP_BlockMatrixNetwork::OpenSMOKE_KPP_BlockMatrixNetwork()
{
	iAlreadyInitialized = false;
	absTolerance_  = 1.e-11;
	relTolerance_  = 1.e-8;
	maxIterations_ = 1000;
}

OpenSMOKE_KPP_BlockMatrixNetwork::~OpenSMOKE_KPP_BlockMatrixNetwork()
{
}

void OpenSMOKE_KPP_BlockMatrixNetwork::SetAbsoluteTolerance(const double absTolerance)
{
	absTolerance_  = absTolerance;
}

void OpenSMOKE_KPP_BlockMatrixNetwork::SetRelativeTolerance(const double relTolerance)
{
	relTolerance_  = relTolerance;
}

void OpenSMOKE_KPP_BlockMatrixNetwork::SetSparsityPattern(const int NR, const int NC)
{
	if (iAlreadyInitialized == false)
	{
		NR_ = NR;
		NC_ = NC;
		NE_ = NR_*NC_;

		ChangeDimensions(NC_, &localRow_);

		J_			= new BzzFactorizedGauss[NR_+1];
		Block_		= new BzzMatrix[NR_+1];
		diagonals_	= new BzzVector[NR_+1];
		iDiagonals_ = new BzzVectorInt[NR_+1];
		b_			= new BzzVector[NR_+1];
		rhs_		= new BzzVector[NR_+1];
		x_			= new BzzVector[NR_+1];
		xOld_		= new BzzVector[NR_+1];

		for(int k=1;k<=NR_;k++)
		{
			ChangeDimensions(NC_, &b_[k]);
			ChangeDimensions(NC_, &rhs_[k]);
			ChangeDimensions(NC_, &x_[k]);
			ChangeDimensions(NC_, &xOld_[k]);
		}

		ChangeDimensions(NR_, &maxDx_);
		ChangeDimensions(NR_, &maxX_);

		maxIterations_ = NE_+50; 

		iAlreadyInitialized = true;
	}
	else
	{
		if (NC_!=NC || NR_!=NR)
			ErrorMessage("Already initialized with a different structure!");
	}
}

void OpenSMOKE_KPP_BlockMatrixNetwork::SetSparsityPattern(const int iBlock, const BzzVectorInt& indices)
{
	iDiagonals_[iBlock] = indices;
	ChangeDimensions(indices.Size(), &diagonals_[iBlock]);
}

void OpenSMOKE_KPP_BlockMatrixNetwork::SetBlock(const int iBlock, BzzMatrix &block)
{
	Block_[iBlock] = block;
	J_[iBlock] = block;
}

void OpenSMOKE_KPP_BlockMatrixNetwork::SetDiagonals(const int iBlock, BzzVector &diagonals)
{
	diagonals_[iBlock] = diagonals;
}

void OpenSMOKE_KPP_BlockMatrixNetwork::SetBlockAndDiagonals(const int iBlock, BzzMatrix &block, BzzVector &diagonals)
{
	SetBlock(iBlock, block);
	SetDiagonals(iBlock, diagonals);
}

void OpenSMOKE_KPP_BlockMatrixNetwork::SetRHS(const int iBlock, const BzzVector &rhs_local, const BzzVector &rhs_non_local)
{
	Sum(rhs_local, rhs_non_local, &rhs_[iBlock]);
}

int OpenSMOKE_KPP_BlockMatrixNetwork::GaussSiedel(BzzVector& solution)
{
	int countFailure_;
	double maxDxGlobal_;
	double maxDxOldGlobal_;
	double maxXGlobal_;

	// Copying the solution
	{
		int j=1;
		for(int k=1;k<=NR_;k++)
			for(int i=1;i<=NC_;i++)
				x_[k][i] = solution[j++];
	}


	cout << setw(12) << "Iteration" << setw(16) << "NormInf" << setw(16) << "Norm1" << setw(16) << "Norm1(mean)" << endl;

	double timeStart = BzzGetCpuTime();
	int iteration;
	int iForward=1;
	for(iteration=1;iteration<=maxIterations_;iteration++)
	{
		for(int k=1;k<=NR_;k++)
			xOld_[k] = x_[k];

		if (iForward==1)
		{
			for(int k=1;k<=NR_;k++)
			{
				// Adding outer terms
				b_[k] = rhs_[k];
				for(int kk=1;kk<=iDiagonals_[k].Size();kk++)
				{
					double coeff = diagonals_[k][kk];
					int block = iDiagonals_[k][kk];
					for(int i=1;i<=NC_;i++)
						b_[k][i] -= x_[block][i]*coeff;
				}

				// Solving
				Solve(&J_[k],&b_[k]);
				x_[k] = b_[k];
			}
		}
		else
		{
			for(int k=NR_;k>=1;k--)
			{
				// Adding outer terms
				b_[k] = rhs_[k];
				for(int kk=1;kk<=iDiagonals_[k].Size();kk++)
				{
					double coeff = diagonals_[k][kk];
					int block = iDiagonals_[k][kk];
					for(int i=1;i<=NC_;i++)
						b_[k][i] -= x_[block][i]*coeff;
				}

				// Solving
				Solve(&J_[k],&b_[k]);
				x_[k] = b_[k];
			}
		}
		iForward*=-1;

		for(int k=1;k<=NR_;k++)
		{
			Difference(x_[k],&xOld_[k]);
			maxDx_[k] = xOld_[k].MaxAbs();
			maxX_[k] = x_[k].GetSumAbsElements();
		}

		maxDxGlobal_	= maxDx_.Max();
		maxXGlobal_		= maxX_.GetSumElements();

		// Write residuals on video
		cout << setw(12) << iteration << setw(16) << maxDxGlobal_ << setw(16) << maxXGlobal_ << setw(16) << maxXGlobal_/double(NE_) << endl;

		if( maxDxGlobal_ < relTolerance_*maxXGlobal_ + double(NE_)*absTolerance_)
		{
			cout << "Solution succesfully reached in Block-Gauss-Siedel mode..." << endl;
			break;
		}

		// In case of first iteration
		if(iteration == 1)
		{
			countFailure_ = 0;
			maxDxOldGlobal_ = maxDxGlobal_;
		}
		else
		{
			if(maxDxGlobal_ > maxDxOldGlobal_ || maxDxGlobal_ > 1000.)
				countFailure_++;
			
			maxDxOldGlobal_ = maxDxGlobal_;
			
			if( (countFailure_>100) || (maxDxGlobal_>10.*maxXGlobal_) || (maxDxGlobal_ > 1.e10) )
			{
				if (countFailure_>100)
					cout << "Failure of Block-Gauss-Siedel mode: max number of failures..." << endl;
				if (maxDxGlobal_>10.*maxXGlobal_)
					cout << "Failure of Block-Gauss-Siedel mode: max relative dx..." << endl;
				if (maxDxGlobal_ > 1.e10)
					cout << "Failure of Block-Gauss-Siedel mode: max absolute dx..." << endl;
				return -1;
			}
		}
	}
	double timeEnd = BzzGetCpuTime();

	cout << "Block Gauss-Siedel CPU Time: " << timeEnd-timeStart << " (" << (timeEnd-timeStart)/double(iteration) << ")" << endl;
		
	// Copying solution
	{
		int j=1;
		for(int k=1;k<=NR_;k++)
			for(int i=1;i<=NC_;i++)
				solution[j++] = x_[k][i];
	}

	if(iteration >= maxIterations_)
	{
		cout << "Maximum number of iterations in Block-Gauss-Siedel mode..." << endl;
		return 1;
	}

	return 0;
}
