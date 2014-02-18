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

#include "mkl.h"
#include <iomanip>
#include "OpenSMOKE_PARDISO_Unsymmetric.h"

OpenSMOKE_PARDISO_Unsymmetric::OpenSMOKE_PARDISO_Unsymmetric(void)
{
	ErrorMessage("No default constructor is available");
}

OpenSMOKE_PARDISO_Unsymmetric::OpenSMOKE_PARDISO_Unsymmetric(const OpenSMOKE_DirectLinearSolver_Unsymmetric_Kind kind)
{
	SetUserDefaultOptions(kind);
}

void OpenSMOKE_PARDISO_Unsymmetric::SetUserDefaultOptions(const OpenSMOKE_DirectLinearSolver_Unsymmetric_Kind kind)
{
	iCStyleIndexing				= true;
	kind_						= kind;
	status_						= OPENSMOKE_DIRECTSOLVER_STATUS_OPEN;

	// Solver internal data address pointer
	for (int j=0;j<64;j++)
		pt[j] = 0;
	
	// (1) Default options vs User options
	iparm[0] = 1;	// 0=Default values    
					// 1=User options
	
	// (2) Fill-in reducing ordering for input matrix
	iparm[1] = 2;	// 0=MDA  2=METIS  3=OpenMP

	// (4) Preconditioned CGS
	iparm[3] = 0;

	// (5) User permutation
	iparm[4] = 0;	// 0=array perm is not used by PARDISO  1=user supplied fill-in reducing permutation
					// 2=PARDISO returns the permutation vector into the array perm

	// (6) Write solution on x
	iparm[5] = 0;	// 0=the array x contains the solution and the value of b is not changed
					// 1=the solver stores the solution on the right-hand side b

	// (8) Iterative refinement step
	iparm[7] = 0;	// must be set to the maximum number of iterative refinement steps that the solver performs
					// The solver automatically performs two steps of iterative refinements when perturbed pivots are obtained during the numerical factorization and iparm(8) = 0

	// (10) Pivoting perturbation
	iparm[9] = 13;	// means eps 1.e-13

	// (11) Scaling vectors
	iparm[10] = 1;	

	// (13) Improved accuracy using (non-)symmetric weighted matchings
	iparm[12] = 1;

	// (18) Number of non-zero elements info
	iparm[17] = -1;	// <0 means the solver reports the numbers of non-zero elements in the factors 

	// (19) MFlops of factorization
	iparm[18] = 0;	// <0 means the solver reports MFlop (1.0E6) that are necessary to factor the matrix A 

	// (21) Pivoting for symmetric indefinite matrices
	iparm[20] = 1;	// 0=1x1 diagonal pivoting is used
					// 1=1x1 and 2x2 Bunch and Kaufman pivoting is used in the factorization process

	// (24) Parallel factorization control
	iparm[23] = 0;	// 0 = PARDISO uses the previous algorithm for factorization
					// 1 = PARDISO uses new two-level factorization algorithm. This algorithm generally
	                //     improves scalability in case of parallel factorization on many threads (>= 8).

	// (27) Matrix checker
	iparm[26] = 1;	// 0=PARDISO does not check the sparse matrix representation
					// 1=PARDISO check integer arrays ia and ja

	// (28) Single or double precision of PARDISO
	iparm[27] = 0;	// 0=double precision
					// 1=single precision
	
	// (31) Partial solution for sparse right-hand sides and sparse solution
	iparm[30] = 0;

	// (35) C or Fortran style array indexing
	iparm[34] = 1;	// 0=PARDISO uses Fortran style indexing
					// 1=PARDISO uses C style indexing

	// (60) Version of PARDISO
	iparm[59] = 0;	// 0=in-core PARDISO used
					// 1=in-core PARDISO used if total memory (MB) needed for storing the factors is < the value of MKL_PARDISO_OOC_MAX_CORE_SIZE (default=2000), and OOC PARDISO is used otherwise
					// 2=the OOC PARDISO used

	// Message level information
	msglvl = 1;	// 0=PARDISO generates no output
				// 1=the solver prints statistical information to the screen

	// Maximal number of factors with identical nonzero sparsity structure that the user 
	// would like to keep at the same time in memor
	maxfct = 1;
    
	// Actual matrix for the solution phase. 
	// With this scalar you can define the matrix that you would like to factorize. 
	mnum = 1;
	
	// This scalar value defines the matrix type
	mtype = 11;		// Unsymmetric matrix
	// mtype = 1;		// Structurally symmetric matrix 
}


OpenSMOKE_PARDISO_Unsymmetric::~OpenSMOKE_PARDISO_Unsymmetric(void)
{
}

void OpenSMOKE_PARDISO_Unsymmetric::SetSparsityPattern(BzzMatrixSparse &M)
{
	if (kind_ != OPENSMOKE_DIRECTSOLVER_SQUAREMATRIX)
		ErrorMessage("SetSparsityPattern(BzzMatrixSparse &M) can be used only for OPENSMOKE_DIRECTSOLVER_SQUAREMATRIX");

	if (status_	!= OPENSMOKE_DIRECTSOLVER_STATUS_OPEN)
		ErrorMessage("SetSparsityPattern(BzzMatrixSparse &M) cannot be used if the Matrix was already closed");

	n = M.Rows();
	rows = new int[n+1];
	perm = new int[n];
	rows[0] = 1;

	// Counting non zero elements
	{
		double* ptrVal;
		int i, j;
		double val;

		int count = 0;
		int iRows = 1;
		M.BeginScanning();
		while(ptrVal = M.Scanning(&i,&j,&val))
		{
			if (iRows != i)
			{
				rows[iRows] = count+1;
				iRows = i;
			}
			count++;
		}
		rows[n] = count+1;
		numberNonZeroElements = count;
	}

	values  = new double[numberNonZeroElements];
	columns = new int[numberNonZeroElements];
	
	// Filling non zero elements
	{	
		double* ptrVal;
		int i, j;
		double val;
		int count = 0;

		M.BeginScanning();
		while(ptrVal = M.Scanning(&i,&j,&val))
		{
			columns[count] = j;
			values[count]  = val;
			count++;
		}
	}

	// Complete Matrix
	CompleteMatrix();

	// Update Status
	status_	= OPENSMOKE_DIRECTSOLVER_STATUS_TOFACTORIZE;
}

void OpenSMOKE_PARDISO_Unsymmetric::UpdateMatrix(BzzMatrixSparse &M)
{
	if (kind_ != OPENSMOKE_DIRECTSOLVER_SQUAREMATRIX)
		ErrorMessage("UpdateMatrix(BzzMatrixSparse &M) can be used only for OPENSMOKE_DIRECTSOLVER_SQUAREMATRIX");

	if (status_	== OPENSMOKE_DIRECTSOLVER_STATUS_OPEN)
		ErrorMessage("SetSparsityPattern(BzzMatrixSparse &M) cannot be used if the Matrix was not closed");

	// Filling non zero elements
	int count = 0;
	{	
		double* ptrVal;
		int i, j;
		double val;

		M.BeginScanning();
		while(ptrVal = M.Scanning(&i,&j,&val))
		{
			values[count]  = val;
			count++;
		}
	}

	// Update Status
	status_	= OPENSMOKE_DIRECTSOLVER_STATUS_TOFACTORIZE;
}

void OpenSMOKE_PARDISO_Unsymmetric::OpenMatrix(const int nRows, const int numberNonZeroElements_)
{
	if (kind_ != OPENSMOKE_DIRECTSOLVER_SQUAREMATRIX_ROW_BY_ROW)
		ErrorMessage("OpenMatrix(const int nRows, const int numberNonZeroElements) can be used only for OPENSMOKE_DIRECTSOLVER_SQUAREMATRIX_ROW_BY_ROW");

	if (status_	!= OPENSMOKE_DIRECTSOLVER_STATUS_OPEN)	
	{
		Delete();
		SetUserDefaultOptions(kind_);
	}

	// Initialization
	numberNonZeroElements = numberNonZeroElements_;
	n = nRows;
	rows = new int[n+1];
	perm = new int[n];
	
	values  = new double[numberNonZeroElements];
	columns = new int[numberNonZeroElements];

	countGlobal_=0;
	countGlobalRows_=0;
	rows[n] = numberNonZeroElements+1;
}

void OpenSMOKE_PARDISO_Unsymmetric::SetSparsityPattern(const int nElementsPerRows, const BzzVectorInt &indicesColumns)
{
	if (kind_ != OPENSMOKE_DIRECTSOLVER_SQUAREMATRIX_ROW_BY_ROW)
		ErrorMessage("SetSparsityPattern(const int nElementsPerRows, const BzzVectorInt &indicesColumns) can be used only for OPENSMOKE_DIRECTSOLVER_SQUAREMATRIX_ROW_BY_ROW");

	if (status_	!= OPENSMOKE_DIRECTSOLVER_STATUS_OPEN)
		ErrorMessage("SetSparsityPattern(const int nElementsPerRows, const BzzVectorInt &indicesColumns) can be used only for Open Matrices");

	// Non zero elements
	{		
		int nLocalRows = indicesColumns.Size() / nElementsPerRows;
		for(int k=1;k<=nLocalRows;k++)
			rows[countGlobalRows_++] = countGlobal_+1 + nElementsPerRows*(k-1);

		for(int k=1;k<=indicesColumns.Size();k++)
		{
			columns[countGlobal_] = indicesColumns[k];
			values[countGlobal_]  = 0.;
			countGlobal_++;
		}
	}
}

void OpenSMOKE_PARDISO_Unsymmetric::SetSparsityPattern(const int index, const int blockDimension, const BzzVector& mConvectionDiffusion, const BzzVectorInt& iConvectionDiffusion, const int nElementsPerRows)
{
	if (kind_ != OPENSMOKE_DIRECTSOLVER_SQUAREMATRIX_ROW_BY_ROW)
		ErrorMessage("SetSparsityPattern(const int nElementsPerRows, const BzzVectorInt &indicesColumns) can be used only for OPENSMOKE_DIRECTSOLVER_SQUAREMATRIX_ROW_BY_ROW");

	if (status_	!= OPENSMOKE_DIRECTSOLVER_STATUS_OPEN)
		ErrorMessage("SetSparsityPattern(const int nElementsPerRows, const BzzVectorInt &indicesColumns) can be used only for Open Matrices");

	BzzVectorInt indicesColumns;
	CalculateSparsityPattern(index, blockDimension, mConvectionDiffusion, iConvectionDiffusion, indicesColumns);

	// Non zero elements
	{		
		int nLocalRows = indicesColumns.Size() / nElementsPerRows;
		for(int k=1;k<=nLocalRows;k++)
			rows[countGlobalRows_++] = countGlobal_+1 + nElementsPerRows*(k-1);

		for(int k=1;k<=indicesColumns.Size();k++)
		{
			columns[countGlobal_] = indicesColumns[k];
			values[countGlobal_]  = 0.;
			countGlobal_++;
		}
	}
}

void OpenSMOKE_PARDISO_Unsymmetric::CloseMatrix()
{
	if (kind_ != OPENSMOKE_DIRECTSOLVER_SQUAREMATRIX_ROW_BY_ROW)
		ErrorMessage("CloseMatrix() can be used only for OPENSMOKE_DIRECTSOLVER_SQUAREMATRIX_ROW_BY_ROW");

	if (status_	!= OPENSMOKE_DIRECTSOLVER_STATUS_OPEN)
		ErrorMessage("CloseMatrix() can be used only for Open Matrices");

	CompleteMatrix();

	// Update status
	status_ = OPENSMOKE_DIRECTSOLVER_STATUS_TOFACTORIZE;
}

void OpenSMOKE_PARDISO_Unsymmetric::UpdateMatrix(const BzzVector &valuesVector)
{
	if (kind_ != OPENSMOKE_DIRECTSOLVER_SQUAREMATRIX_ROW_BY_ROW)
		ErrorMessage("UpdateMatrix(const BzzVector &valuesVector) can be used only for OPENSMOKE_DIRECTSOLVER_SQUAREMATRIX_ROW_BY_ROW");

	if (status_	== OPENSMOKE_DIRECTSOLVER_STATUS_OPEN)
		ErrorMessage("UpdateMatrix(const BzzVector &valuesVector) can be used only for Open Matrices");

	// Non zero elements
	{		
		for(int k=1;k<=valuesVector.Size();k++)
		{
			values[countGlobal_]  = valuesVector[k];
			countGlobal_++;
		}
	}

	// Update status
	status_ = OPENSMOKE_DIRECTSOLVER_STATUS_TOFACTORIZE;
}

void OpenSMOKE_PARDISO_Unsymmetric::CompleteMatrix()
{
	if (status_	!= OPENSMOKE_DIRECTSOLVER_STATUS_OPEN)
		ErrorMessage("CompleteMatrix() cannot be used if the Matrix was already closed");

	cout << " Number of non-zero elements:    " << numberNonZeroElements << endl;
	cout << " Sparsity fill-in:               " << double(numberNonZeroElements)/double(n*n)*100. << endl;
	cout << " Mean non zero elements per row: " << double(numberNonZeroElements)/double(n) << endl;

	// Indices revision
	if (iCStyleIndexing==true)
	{
		SetCStyleIndexing();
		for(int j=1;j<=n+1;j++)
			rows[j-1] -= 1;
		for(int j=1;j<=numberNonZeroElements;j++)
			columns[j-1] -= 1;
	}

	// Only Analysis in this phase
	nrhs  = 1;  // Number of RHSs
	phase = 11;	// Analysis

	double* b_ = new double[n];	// dummy variables
	double* x_ = new double[n];	// dummy variables	

	for(int j=0;j<n;j++)
		b_[j] = 0.;
	
	if (msglvl==1)	MessageOnTheScreen("PARDISO: Sparse linear system analysis");
	pardiso(pt, &maxfct, &mnum, &mtype, &phase, &n, values, rows, columns, perm, &nrhs, iparm, &msglvl, b_, x_, &error_);
	ErrorAnalysis();

	// Update Status
	status_	= OPENSMOKE_DIRECTSOLVER_STATUS_TOFACTORIZE;
}

void OpenSMOKE_PARDISO_Unsymmetric::NumericalFactorization()
{
	if (status_	== OPENSMOKE_DIRECTSOLVER_STATUS_OPEN)
		ErrorMessage("NumericalFactorization() cannot be used for open matrices");

	if (status_	== OPENSMOKE_DIRECTSOLVER_STATUS_FACTORIZED)
		return;

	double* x_ = new double[n];
	double* b_ = new double[n];

	for(int j=0;j<n;j++)
		b_[j] = 0.;//b[j+1];

	nrhs  = 1;  // Number of RHSs
	phase = 22;	// Numerical factorization
	
	if (msglvl==1)	MessageOnTheScreen("PARDISO: Numerical factorization");
	pardiso(pt, &maxfct, &mnum, &mtype, &phase, &n, values, rows, columns, perm, &nrhs, iparm, &msglvl, b_, x_, &error_);
	ErrorAnalysis();

	// Update Status
	status_	= OPENSMOKE_DIRECTSOLVER_STATUS_FACTORIZED;
}

void OpenSMOKE_PARDISO_Unsymmetric::Solve(BzzVector &b, BzzVector &x)
{	
	if (status_	!= OPENSMOKE_DIRECTSOLVER_STATUS_FACTORIZED)
		ErrorMessage("Solve(BzzVector &b, BzzVector &x) can be used only for Factorized Matrices");

	// Checking dimensions
	if (b.Size() != x.Size())	ErrorMessage("The size of b and x do not fit!");
	if (b.Size() != n)			ErrorMessage("The size of matrix and rhs do not fit!");

	// Linear system data
	nrhs  = 1;	// Number of RHSs
	phase = 33;	// Solve, iterative refinement

	// Memory allocation
	double* x_ = x.GetHandle();
	double* b_ = b.GetHandle();
	
	if (msglvl==1)	MessageOnTheScreen("PARDISO: Solve sparse linear system");
	pardiso(pt, &maxfct, &mnum, &mtype, &phase, &n, values, rows, columns, perm, &nrhs, iparm, &msglvl, b_, x_, &error_);
	ErrorAnalysis();
}

void OpenSMOKE_PARDISO_Unsymmetric::Solve(BzzMatrix &b, BzzMatrix &x, const bool iVerticalStorage)
{	
	if (status_	!= OPENSMOKE_DIRECTSOLVER_STATUS_FACTORIZED)
		ErrorMessage("Solve(BzzMatrix &b, BzzMatrix &x, const bool iVerticalStorage) can be used only for Factorized Matrices");

	// Checking dimensions
	if (b.Rows() != x.Rows())		ErrorMessage("The row size of b and x do not fit!");
	if (b.Columns() != x.Columns())	ErrorMessage("The column size of b and x do not fit!");

	if (iVerticalStorage == true)
		if (b.Rows() != n)		ErrorMessage("The size of matrix and rhs do not fit!");
	if (iVerticalStorage == false)
		if (b.Columns() != n)	ErrorMessage("The size of matrix and rhs do not fit!");

	// Linear system data
	phase = 33;												// Solve, iterative refinement
	if (iVerticalStorage == true)	nrhs  = b.Columns();	// Number of RHSs
	else							nrhs  = b.Rows();		// Number of RHSs
	
	// Memory Allocation
	double* x_ = new double[n*nrhs];
	double* b_ = new double[n*nrhs];

	// From BzzVector to C vector
	if (iVerticalStorage == true)
	{
		int j=0;
		for(int k=1;k<=nrhs;k++)
			for(int i=1;i<=n;i++)
				b_[j++] = b[i][k];
	}
	else
	{
		int j=0;
		for(int k=1;k<=nrhs;k++)
			for(int i=1;i<=n;i++)
				b_[j++] = b[k][i];
	}

	if (msglvl==1)	MessageOnTheScreen("PARDISO: Solve sparse linear system");
	pardiso(pt, &maxfct, &mnum, &mtype, &phase, &n, values, rows, columns, perm, &nrhs, iparm, &msglvl, b_, x_, &error_);
	ErrorAnalysis();

	// From C Vector to BzzVector
	if (iVerticalStorage == true)
	{
		int j=0;
		for(int k=1;k<=nrhs;k++)
			for(int i=1;i<=n;i++)
				x[i][k] = x_[j++];
	}
	else
	{
		int j=0;
		for(int k=1;k<=nrhs;k++)
			for(int i=1;i<=n;i++)
				x[k][i] = x_[j++];
	}
}

void OpenSMOKE_PARDISO_Unsymmetric::ResetCounters()
{
	countGlobal_ = 0;
	status_ = OPENSMOKE_DIRECTSOLVER_STATUS_TOFACTORIZE;
}

void OpenSMOKE_PARDISO_Unsymmetric::Delete()
{
	CleanMemory();

	delete[] perm;
}

void OpenSMOKE_PARDISO_Unsymmetric::SetDefaultOptions()
{
	iparm[0] = 0;	// parameter(1)
}

void OpenSMOKE_PARDISO_Unsymmetric::SetFillInReducingOrdering(const PARDISO_FillInReducingOrdering value)
{
	if (value == PARDISO_FILLIN_MDA)	iparm[1] = 0;	// parameter(2)
	if (value == PARDISO_FILLIN_METIS)	iparm[1] = 2;	// parameter(2)
	if (value == PARDISO_FILLIN_OPENMP)	iparm[1] = 3;	// parameter(2)
}

void OpenSMOKE_PARDISO_Unsymmetric::SetMFlopsFactorization()
{
	iparm[18] = -1;	// parameter(19)
}
	
void OpenSMOKE_PARDISO_Unsymmetric::SetTwoLevelFactorization()
{
	iparm[23] =	1;	// parameter(24)
}

void OpenSMOKE_PARDISO_Unsymmetric::UnsetMatrixChecker()
{
	iparm[26] = 0;	// parameter(27)
}

void OpenSMOKE_PARDISO_Unsymmetric::SetCStyleIndexing()
{
	iparm[34] = 1;	// parameter(35)
	iCStyleIndexing = true;
}

void OpenSMOKE_PARDISO_Unsymmetric::SetFortranStyleIndexing()
{
	iparm[34] = 0;	// parameter(35)
	iCStyleIndexing = false;
}

void OpenSMOKE_PARDISO_Unsymmetric::UnsetMessage()
{
	msglvl = 0;
}

int OpenSMOKE_PARDISO_Unsymmetric::NumberOfPerformedIterativeRefinementSteps()
{
	return iparm[6];	// parameter(7)
}

int OpenSMOKE_PARDISO_Unsymmetric::NumberOfPerturbedPivot()
{
	return iparm[13];	// parameter(14)
}

int OpenSMOKE_PARDISO_Unsymmetric::PeakMemorySymbolicFactorization()
{
	// Reports the total peak memory in KBytes that the solver needs during the analysis and 
	// symbolic factorization phase. This value is only computed in phase 1.

	return iparm[14];	// parameter(15)
}

int OpenSMOKE_PARDISO_Unsymmetric::PermanentMemorySymbolicFactorization()
{
	// Permanent memory in KBytes from the analysis and symbolic factorization phase that the solver 
	// needs in the factorization and solve phases. This value is only computed in phase 1.

	return iparm[15];	// parameter(16)
}

int OpenSMOKE_PARDISO_Unsymmetric::SizeOfFactors()
{
    // This parameter is computed in phase 1 and 2 with different meanings.
    // 1. The parameter computed in phase 1 provides the user with an estimate of size of factors (KBytes)
	//    In OOC mode, this information gives an estimate of disk space for storing data required.
	// 2. The parameter computed in phase 2 provides the user with the total double precision memory 
	//    consumption (KBytes) of the solver for the factorization and solve phases.

	return iparm[16];	// parameter(17)
}

int OpenSMOKE_PARDISO_Unsymmetric::NumberOfNonZeroElementsInFactors()
{
    // Number of non-zero elements in factors

	return iparm[17];	// parameter(18)
}

int OpenSMOKE_PARDISO_Unsymmetric::MFlopsFactorization()
{
    // The solver reports the number of operations in MFlop (1.0E6 operations) that are 
	// necessary to factor the matrix A.

	return iparm[18];	// parameter(19)
}

int OpenSMOKE_PARDISO_Unsymmetric::CGCGSDiagnostics()
{
	// The value is used to give CG/CGS diagnostics:
	// >0: CGS succeeded, and the number of iterations executed are reported in
	// <0, iterations executed, but CG/CGS failed. The error report details in are of the form: iparm= - it_cgs*10 - cgs_error.
	// if phase==23, then the factors L, U are recomputed for the matrix A and the error flag error=0 
	// in case of a successful factorization. If phase =33, then error = -4 signals the failure.

    // Description of cgs_error is given in the table below:
    //   1: fluctuations of the residuum are too large
    //   2: ||dxmax_it_cgs/2|| is too large (slow convergence)
    //   3: stopping criterion is not reached at max_it_cgs
    //   4: perturbed pivots causes iterative refinement
    //   5: factorization is too fast for this matrix. It is better to use the factorization 
	//      method with iparm(4)=0

	return iparm[19];	// parameter(20)
}

int OpenSMOKE_PARDISO_Unsymmetric::InertiaPositive()
{
	// Number of positive eigenvalues for symmetric indefinite matrices
	
	return iparm[21];	// parameter(22)
}   

int OpenSMOKE_PARDISO_Unsymmetric::InertiaNegative()
{
	// Number of negative eigenvalues for symmetric indefinite matrices
	
	return iparm[22];	// parameter(23)
}   

int OpenSMOKE_PARDISO_Unsymmetric::ZeroOrNegativePivot()
{
	// Number of equation where PARDISO detects zero or negative pivot for MTYPE=2 
	// (real positive definite matrix) and MTYPE=4 (complex and Hermitian positive definite matrices). 
	// If the solver detects a zero or negative pivot for these matrix types, the factorization is stopped,
	// PARDISO returns immediately with an error (error=-4) and this parameter contains the number
	// of the equation where the first zero or negative pivot is detected

	return iparm[29];	// parameter(30)
} 

void OpenSMOKE_PARDISO_Unsymmetric::ErrorAnalysis()
{
	if (error_ == 0)	return;
	else if (error_ == -1)	ErrorMessage("Input inconsistent");
	else if (error_ == -2)	ErrorMessage("Not enough memory");
	else if (error_ == -3)	ErrorMessage("Reordering problem");
	else if (error_ == -4)	ErrorMessage("Zero pivot, numerical factorization or iterative refinement problem");
	else if (error_ == -5)	ErrorMessage("Unclassified (internal) error");
	else if (error_ == -6)	ErrorMessage("Preordering failed (matrix types 11, 13 only)");
	else if (error_ == -7)	ErrorMessage("Diagonal matrix is singular");
	else if (error_ == -8)	ErrorMessage("32-bit integer overflow problem");
	else if (error_ == -9)	ErrorMessage("Not enough memory for OOC");
	else if (error_ == -10)	ErrorMessage("Problems with opening OOC temporary files");
	else if (error_ == -11)	ErrorMessage("Read/write problems with the OOC data file");
}



