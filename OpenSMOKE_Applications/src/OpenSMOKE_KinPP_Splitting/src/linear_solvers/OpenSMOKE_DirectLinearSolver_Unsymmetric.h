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

#ifndef OpenSMOKE_DirectLinearSolver_Unsymmetric_H
#define OpenSMOKE_DirectLinearSolver_Unsymmetric_H

#include "BzzMath.hpp"
using namespace std;

enum OpenSMOKE_DirectLinearSolver_Unsymmetric_Kind   { OPENSMOKE_DIRECTSOLVER_SQUAREMATRIX, OPENSMOKE_DIRECTSOLVER_SQUAREMATRIX_ROW_BY_ROW };
enum OpenSMOKE_DirectLinearSolver_Unsymmetric_Status { OPENSMOKE_DIRECTSOLVER_STATUS_OPEN, OPENSMOKE_DIRECTSOLVER_STATUS_TOFACTORIZE, OPENSMOKE_DIRECTSOLVER_STATUS_FACTORIZED };

class OpenSMOKE_DirectLinearSolver_Unsymmetric
{
public:

	OpenSMOKE_DirectLinearSolver_Unsymmetric(void) {};
	OpenSMOKE_DirectLinearSolver_Unsymmetric(const OpenSMOKE_DirectLinearSolver_Unsymmetric_Kind kind) {};
	~OpenSMOKE_DirectLinearSolver_Unsymmetric(void) {};

	// 1a. Set from BzzMatrix Sparse + Symbolic Factorization
	virtual void SetSparsityPattern(BzzMatrixSparse &M) = 0;
	virtual void UpdateMatrix(BzzMatrixSparse &M) = 0;
	
	// 1b. Set from BzzVector + Symbolic Factorization
	virtual void OpenMatrix(const int nRows, const int numberNonZeroElements) = 0;
	virtual void SetSparsityPattern(const int nElementsPerRows, const BzzVectorInt &indicesColumns) = 0;
	virtual void SetSparsityPattern(const int index, const int blockDimension, const BzzVector& mConvectionDiffusion, const BzzVectorInt& iConvectionDiffusion, const int nElementsPerRows) = 0;	
	virtual void CloseMatrix() = 0;
	virtual void UpdateMatrix(const BzzVector &valuesVector) = 0;
	virtual void ResetCounters() = 0;

	// 2. Numerical Factorization
	virtual void NumericalFactorization() = 0;
	virtual void Solve(BzzVector &b, BzzVector &x) = 0;
	virtual void Solve(BzzMatrix &b, BzzMatrix &x, const bool iVerticalStorage) = 0;

	// 3. Delete
	virtual void Delete() = 0;

public:

	virtual void UnsetMessage() = 0;

protected:

	// Kind of matrix
	OpenSMOKE_DirectLinearSolver_Unsymmetric_Kind kind_;
	OpenSMOKE_DirectLinearSolver_Unsymmetric_Status status_;
	
	// Error message
	int error_;		

	// System data
	int n;						// Number of equations
	int	nrhs;					// Number of rhs
	int numberNonZeroElements;	// Number of non zero elements
	
	// Matrix data
	double *values;		//	Matrix elements
	int    *rows;		//	Matrix rows
	int    *columns;	//	Matrix columns

	// Counters
	int countGlobal_;
	int countGlobalRows_;
		
protected:

	virtual void ErrorAnalysis() = 0;
	virtual void CompleteMatrix() = 0;
	virtual void SetUserDefaultOptions(const OpenSMOKE_DirectLinearSolver_Unsymmetric_Kind kind) = 0;
	
protected:

	string name_solver_;
	void ErrorMessage(const string message_);
	void WarningMessage(const string message_);
	void MessageOnTheScreen(const string message_);
	void CleanMemory();
	void CalculateSparsityPattern(const int index, const int blockDimension, const BzzVector& mConvectionDiffusion, const BzzVectorInt& iConvectionDiffusion, BzzVectorInt &indicesColumns);

};

#endif	// OpenSMOKE_DirectLinearSolver_Unsymmetric_H
