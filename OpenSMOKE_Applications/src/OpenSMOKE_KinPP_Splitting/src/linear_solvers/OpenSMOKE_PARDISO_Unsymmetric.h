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

#ifndef OpenSMOKE_PARDISO_Unsymmetric_H
#define OpenSMOKE_PARDISO_Unsymmetric_H

#include "BzzMath.hpp"
#include "linear_solvers/OpenSMOKE_DirectLinearSolver_Unsymmetric.h"

enum PARDISO_FillInReducingOrdering { PARDISO_FILLIN_MDA, PARDISO_FILLIN_METIS, PARDISO_FILLIN_OPENMP };

class OpenSMOKE_PARDISO_Unsymmetric : public OpenSMOKE_DirectLinearSolver_Unsymmetric
{
public:

	// Constructors and destructors
	OpenSMOKE_PARDISO_Unsymmetric(void);
	OpenSMOKE_PARDISO_Unsymmetric(OpenSMOKE_DirectLinearSolver_Unsymmetric_Kind kind);
	~OpenSMOKE_PARDISO_Unsymmetric(void);

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
	void Solve(BzzVector &b, BzzVector &x);
	void Solve(BzzMatrix &b, BzzMatrix &x, const bool iVerticalStorage);

	// 3. Delete
	void Delete();

public:

	// Options (Set/Unset)
	void SetDefaultOptions();
	void SetFillInReducingOrdering(const PARDISO_FillInReducingOrdering value);
	void SetMFlopsFactorization();
	void SetTwoLevelFactorization();
	void SetCStyleIndexing();
	void SetFortranStyleIndexing();
	void UnsetMatrixChecker();
	void UnsetMessage();

	// Get data
	int NumberOfPerformedIterativeRefinementSteps();
	int NumberOfPerturbedPivot();
	int PeakMemorySymbolicFactorization();
	int PermanentMemorySymbolicFactorization();
	int SizeOfFactors();
	int NumberOfNonZeroElementsInFactors();
	int MFlopsFactorization();
	int CGCGSDiagnostics();
	int InertiaPositive();
	int InertiaNegative();
	int ZeroOrNegativePivot();

private:

	void CompleteMatrix();
	void ErrorAnalysis();
	void SetUserDefaultOptions(const OpenSMOKE_DirectLinearSolver_Unsymmetric_Kind kind);

private:
	
	void* pt[64];			// Internal data
	int maxfct;				// Max number of factors with identical sparsity pattern
	int mnum;				// Actual matrix
	int mtype;				// Matrix type
	int phase;				// Phase of calculation
	int msglvl;				// Message level information
	int iparm[64];			// Options
	bool iCStyleIndexing;	// C-Style indices
	
	int    *perm;			// Permutation index
};

#endif	// OpenSMOKE_PARDISO_Unsymmetric_H
