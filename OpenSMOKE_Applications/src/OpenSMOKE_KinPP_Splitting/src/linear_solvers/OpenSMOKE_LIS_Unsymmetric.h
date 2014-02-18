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

#ifndef OpenSMOKE_LIS_Unsymmetric_H
#define OpenSMOKE_LIS_Unsymmetric_H

#include "lis.h"
#include "BzzMath.hpp"
#include "OpenSMOKE_KPP_Definitions.h"
#include "linear_solvers/OpenSMOKE_DirectLinearSolver_Unsymmetric.h"

class OpenSMOKE_LIS_Unsymmetric : public OpenSMOKE_DirectLinearSolver_Unsymmetric
{
public:

	// Constructors and destructors
	OpenSMOKE_LIS_Unsymmetric(void);
	OpenSMOKE_LIS_Unsymmetric(OpenSMOKE_DirectLinearSolver_Unsymmetric_Kind kind);
	~OpenSMOKE_LIS_Unsymmetric(void);

	// 1a. Set from BzzMatrix Sparse + Symbolic Factorization
	void SetSparsityPattern(BzzMatrixSparse &M);
	void UpdateMatrix(BzzMatrixSparse &M);
	
	// 1b. Set from BzzVector + Symbolic Factorization
	void OpenMatrix(const int nRows, const int numberNonZeroElements);
	void SetSparsityPattern(const int nElementsPerRows, const BzzVectorInt &indicesColumns);
	void SetSparsityPattern(const int index, const int blockDimension, const BzzVector& mConvectionDiffusion, const BzzVectorInt& iConvectionDiffusion, const int nElementsPerRows);
	void CloseMatrix();
	void UpdateMatrix(const BzzVector &valuesVector);
	void ResetCounters();
	
	// 2. Numerical Factorization and Numerical Solution
	void NumericalFactorization();
	void SetInitialSolution(BzzVector &firstguess);
	void Solve(BzzVector &b, BzzVector &x);
	void Solve(BzzMatrix &b, BzzMatrix &x, const bool iVerticalStorage);

	// 3. Delete
	void Delete();

public:

	void SetLinearEquationSolver(const linearEquationSolver kind);
	void SetPreconditioner(const lisPreconditioner kind);
	void SetMaximumIterations(const int maximumIterations);
	void SetTolerance(const double tolerance);
	void SetMessage();
	void DisableInitialSolution();
	void EnableInitialSolution();

public:

	// Options (Set/Unset)
	void SetDefaultOptions();
	void UnsetMessage();

	// Get data
    LIS_MATRIX        A;
    LIS_VECTOR        lis_b, lis_x;
	LIS_SOLVER        solver;

private:

	void CompleteMatrix();
	void ErrorAnalysis();
	void SetUserDefaultOptions(const OpenSMOKE_DirectLinearSolver_Unsymmetric_Kind kind);



	string linearEquationSolver_;
	string lisPreconditioner_;
	string messageLevel_;
	string initialSolution_;
	string maximumIterations_;
	string tolerance_;

private:
	
};

#endif	// OpenSMOKE_PARDISO_Unsymmetric_H
