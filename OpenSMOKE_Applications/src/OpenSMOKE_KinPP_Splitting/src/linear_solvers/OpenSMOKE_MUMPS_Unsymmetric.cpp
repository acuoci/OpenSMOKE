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

#include "OpenSMOKE_MUMPS_Unsymmetric.h"

#define ICNTL(I) icntl[(I)-1]
#define RINFO(I) rinfo[(I)-1]
#define INFO(I) info[(I)-1]
#define USE_COMM_WORLD -987654
#define MUMPS_JOB_INIT			-1
#define MUMPS_JOB_END			-2
#define MUMPS_JOB_ANALYSIS		 1
#define MUMPS_JOB_FACTORIZATION  	2
#define MUMPS_JOB_SOLUTION		 3



OpenSMOKE_MUMPS_Unsymmetric::OpenSMOKE_MUMPS_Unsymmetric(void)
{
	ErrorMessage("No default constructor is available");
}

OpenSMOKE_MUMPS_Unsymmetric::OpenSMOKE_MUMPS_Unsymmetric(const OpenSMOKE_DirectLinearSolver_Unsymmetric_Kind kind)
{
	SetUserDefaultOptions(kind);
}

OpenSMOKE_MUMPS_Unsymmetric::~OpenSMOKE_MUMPS_Unsymmetric(void)
{
}

void OpenSMOKE_MUMPS_Unsymmetric::SetUserDefaultOptions(const OpenSMOKE_DirectLinearSolver_Unsymmetric_Kind kind)
{
	kind_						= kind;
	status_						= OPENSMOKE_DIRECTSOLVER_STATUS_OPEN;

	// Current job
	id.job=MUMPS_JOB_INIT;

	// Host policy
	id.par=1;	// 0 = host is not involved in factorization/solve phases
				// 1 = host is involved in factorization/solve phases
	
	// Kind of matrix
	id.sym=0;	// 0 = unsymmetric matrix
	
	// Communicator
	id.comm_fortran=USE_COMM_WORLD;

	// ICNTL(1) is the output stream for error messages. 
	// If it is negative or zero, these messages will be suppressed. Default value is 6.	
	id.ICNTL(1) = 6;	
	
	// ICNTL(2) is the output stream for diagnostic printing, statistics, and warning messages. 
	// If it is negative or zero, these messages will be suppressed. Default value is 0.
	id.ICNTL(2) = 0;	
	
	// ICNTL(3) is the output stream for global information, collected on the host. 
	// If it is negative or zero, these messages will be suppressed. Default value is 6.
	id.ICNTL(3) = 0;
	
	// ICNTL(4) is the level of printing for error, warning, and diagnostic messages.
	// Maximum value is 4 and default value is 2 (errors and warnings printed). Possible values are
	// 0 : No messages output.
	// 1 : Only error messages printed.
	// 2 : Errors, warnings, and main statistics printed.
	// 3 : Errors and warnings and terse diagnostics (only first ten entries of arrays) printed.
	// 4 : Errors and warnings and information on input and output parameters printed.
	id.ICNTL(4) = 3;

	// ICNTL(5) has default value 0 and is only accessed by the host and only during the analysis phase. 
	// If ICNTL(5) = 0, the input matrix must be given in assembled format
	// If ICNTL(5) = 1, the input matrix must be given in elemental format
//	id.ICNTL(5) =  0;

	// ICNTL(6) has default value 7 (automatic choice done by the package) and is used to control an option for permuting and/or scaling the matrix.
	// 7=Automatic choice	
	id.ICNTL(6) = 7;	

	// ICNTL(7) it determines the pivot order to be used for the factorization
	// 7=Automatic choice	
	id.ICNTL(7) = 7;	
				

	// Initialize
	MessageOnTheScreen("MUMPS: Initialization");
	dmumps_c(&id);
}

void OpenSMOKE_MUMPS_Unsymmetric::SetSparsityPattern(BzzMatrixSparse &M)
{
	if (kind_ != OPENSMOKE_DIRECTSOLVER_SQUAREMATRIX)
		ErrorMessage("SetSparsityPattern(BzzMatrixSparse &M) can be used only for OPENSMOKE_DIRECTSOLVER_SQUAREMATRIX");

	if (status_	!= OPENSMOKE_DIRECTSOLVER_STATUS_OPEN)
		ErrorMessage("SetSparsityPattern(BzzMatrixSparse &M) cannot be used if the Matrix was already closed");

	// Size
	n = M.Rows();

	// Counting non zero elements
	cout << "Counting non zero elements" << endl;
	{
		double* ptrVal;
		int i, j;
		double val;

		numberNonZeroElements = 0;
		M.BeginScanning();
		while(ptrVal = M.Scanning(&i,&j,&val))
			numberNonZeroElements++;
	}

	values  = new double[numberNonZeroElements];
	rows    = new int[numberNonZeroElements];
	columns = new int[numberNonZeroElements];

	// Parsing non zero elements
	cout << "Parsing non zero elements" << endl;
	{
		double* ptrVal;
		int i, j;
		double val;

		numberNonZeroElements = 0;
		M.BeginScanning();
		while(ptrVal = M.Scanning(&i,&j,&val))
		{
			rows[numberNonZeroElements]    = i;
			columns[numberNonZeroElements] = j;
			values[numberNonZeroElements]  = val;

			numberNonZeroElements++;
		}
	}

	// Setting values
	MessageOnTheScreen("MUMPS: Set Sparsity Pattern");
	CompleteMatrix();
}

void OpenSMOKE_MUMPS_Unsymmetric::UpdateMatrix(BzzMatrixSparse &M)
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

void OpenSMOKE_MUMPS_Unsymmetric::OpenMatrix(const int nRows, const int nonZeroElements)
{
	MessageOnTheScreen("MUMPS: OpenMatrix(const int nRows, const int nonZeroElements)");

	if (kind_ != OPENSMOKE_DIRECTSOLVER_SQUAREMATRIX_ROW_BY_ROW)
		ErrorMessage("OpenMatrix(const int nRows, const int numberNonZeroElements) can be used only for OPENSMOKE_DIRECTSOLVER_SQUAREMATRIX_ROW_BY_ROW");

	if (status_	!= OPENSMOKE_DIRECTSOLVER_STATUS_OPEN)
	{
		Delete();
		SetUserDefaultOptions(kind_);
	}

	n = nRows;
	numberNonZeroElements = nonZeroElements;
	rows    = new int[nonZeroElements];
	columns = new int[nonZeroElements];
	values  = new double[nonZeroElements];

	countGlobal_=0;
	countGlobalRows_=0;
}

void OpenSMOKE_MUMPS_Unsymmetric::SetSparsityPattern(const int nElementsPerRows, const BzzVectorInt &indicesColumns)
{
	MessageOnTheScreen("MUMPS: SetSparsityPattern(const int nElementsPerRows, const BzzVectorInt &indicesColumns)");

	if (kind_ != OPENSMOKE_DIRECTSOLVER_SQUAREMATRIX_ROW_BY_ROW)
		ErrorMessage("SetSparsityPattern(const int nElementsPerRows, const BzzVectorInt &indicesColumns) can be used only for OPENSMOKE_DIRECTSOLVER_SQUAREMATRIX_ROW_BY_ROW");

	if (status_	!= OPENSMOKE_DIRECTSOLVER_STATUS_OPEN)
		ErrorMessage("SetSparsityPattern(const int nElementsPerRows, const BzzVectorInt &indicesColumns) can be used only for Open Matrices");

	int nLocalRows = indicesColumns.Size() / nElementsPerRows;
	int count = 1;

	for(int k=1;k<=nLocalRows;k++)
	{
		countGlobalRows_++;
		for(int j=1;j<=nElementsPerRows;j++)
		{
			rows[countGlobal_]    = countGlobalRows_;
			columns[countGlobal_] = indicesColumns[count++];
			values[countGlobal_]  = 1.;

			countGlobal_++;
		}
	}
}

void OpenSMOKE_MUMPS_Unsymmetric::SetSparsityPattern(const int index, const int blockDimension, const BzzVector& mConvectionDiffusion, const BzzVectorInt& iConvectionDiffusion, const int nElementsPerRows)
{
	if (kind_ != OPENSMOKE_DIRECTSOLVER_SQUAREMATRIX_ROW_BY_ROW)
		ErrorMessage("SetSparsityPattern(const int nElementsPerRows, const BzzVectorInt &indicesColumns) can be used only for OPENSMOKE_DIRECTSOLVER_SQUAREMATRIX_ROW_BY_ROW");

	if (status_	!= OPENSMOKE_DIRECTSOLVER_STATUS_OPEN)
		ErrorMessage("SetSparsityPattern(const int nElementsPerRows, const BzzVectorInt &indicesColumns) can be used only for Open Matrices");

	BzzVectorInt indicesColumns;
	CalculateSparsityPattern(index, blockDimension, mConvectionDiffusion, iConvectionDiffusion, indicesColumns);

	int nLocalRows = indicesColumns.Size() / nElementsPerRows;
	int count = 1;

	for(int k=1;k<=nLocalRows;k++)
	{
		countGlobalRows_++;
		for(int j=1;j<=nElementsPerRows;j++)
		{
			rows[countGlobal_]    = countGlobalRows_;
			columns[countGlobal_] = indicesColumns[count++];
			values[countGlobal_]  = 1.;

			countGlobal_++;
		}
	}
}



void OpenSMOKE_MUMPS_Unsymmetric::CloseMatrix()
{
	if (kind_ != OPENSMOKE_DIRECTSOLVER_SQUAREMATRIX_ROW_BY_ROW)
		ErrorMessage("CloseMatrix() can be used only for OPENSMOKE_DIRECTSOLVER_SQUAREMATRIX_ROW_BY_ROW");

	if (status_	!= OPENSMOKE_DIRECTSOLVER_STATUS_OPEN)
		ErrorMessage("CloseMatrix() can be used only for Open Matrices");

	if (numberNonZeroElements != countGlobal_)
		ErrorMessage("The OPENSMOKE_DIRECTSOLVER_SQUAREMATRIX_ROW_BY_ROW was not correctly initialized");
	
	if (n != countGlobalRows_)
		ErrorMessage("The OPENSMOKE_DIRECTSOLVER_SQUAREMATRIX_ROW_BY_ROW was not correctly initialized");

	CompleteMatrix();

	// Update status
	status_ = OPENSMOKE_DIRECTSOLVER_STATUS_TOFACTORIZE;
}

void OpenSMOKE_MUMPS_Unsymmetric::UpdateMatrix(const BzzVector &valuesVector)
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

void OpenSMOKE_MUMPS_Unsymmetric::CompleteMatrix()
{
	if (status_	!= OPENSMOKE_DIRECTSOLVER_STATUS_OPEN)
		ErrorMessage("CompleteMatrix() cannot be used if the Matrix was already closed");
 
	// Setting values
	id.n   = n; 
	id.nz  = numberNonZeroElements; 
	id.irn = rows; 
	id.jcn = columns;
	id.a   = values; 

	id.ICNTL(1) = 6;	// ICNTL(1) is the output stream for error messages. If it is negative or zero, these messages will be suppressed. Default value is 6.
	id.ICNTL(2) = 0;	// ICNTL(2) is the output stream for diagnostic printing, statistics, and warning messages. If it is negative or zero, these messages will be suppressed. Default value is 0.
	id.ICNTL(3) = 0;	// ICNTL(3) is the output stream for global information, collected on the host. If it is negative or zero, these messages will be suppressed. Default value is 6.
	id.ICNTL(4) = 3;	// ICNTL(4) is the level of printing for error, warning, and diagnostic messages. Maximum value is 4 and default value is 2 (errors and warnings printed).
	id.ICNTL(6) = 7;	// ICNTL(6) has default value 7 (automatic choice done by the package) and is used to control an option for permuting and/or scaling the matrix.
				// 7=Automatic choice
	id.ICNTL(7) = 7;	// ICNTL(7) it determines the pivot order to be used for the factorization
				// 7=Automatic choice

	// Statistics
	cout << " Number of non-zero elements:    " << numberNonZeroElements << endl;
	cout << " Sparsity fill-in:               " << double(numberNonZeroElements)/double(n*n)*100. << endl;
	cout << " Mean non zero elements per row: " << double(numberNonZeroElements)/double(n) << endl;

	cout << " Number of rows:                 " << n    << endl; 
	cout << " Number of columns:              " << n << endl; 


	// Only Analysis in this phase
	id.job = MUMPS_JOB_ANALYSIS;	// Analysis
	
	MessageOnTheScreen("MUMPS: Sparse linear system analysis");
	dmumps_c(&id);
	MessageOnTheScreen("MUMPS: Done");
	ErrorAnalysis();

	cout << " After analysis: The estimated number of floating-point operations on the processor for the elimination process: " << id.RINFO(1) << endl; 
	cout << " After analysis:estimated size in Megabytes of all working space to run the numerical phases:                    " << id.INFO(15) << endl; 
 
	// Update Status
	status_	= OPENSMOKE_DIRECTSOLVER_STATUS_TOFACTORIZE;
}


void OpenSMOKE_MUMPS_Unsymmetric::NumericalFactorization()
{
	if (status_	== OPENSMOKE_DIRECTSOLVER_STATUS_OPEN)
		ErrorMessage("NumericalFactorization() cannot be used for open matrices");

	if (status_	== OPENSMOKE_DIRECTSOLVER_STATUS_FACTORIZED)
		return;

	id.job = MUMPS_JOB_FACTORIZATION;	// Numerical factorization
	
	MessageOnTheScreen("MUMPS: Numerical factorization");
	dmumps_c(&id);
	MessageOnTheScreen("MUMPS: Done");
	ErrorAnalysis();

	cout << " After factorization: The number of floating-point operations on the processor for the assembly process:    " << id.RINFO(2) << endl; 
	cout << " After factorization: The number of floating-point operations on the processor for the elimination process: " << id.RINFO(3) << endl; 

	// Update Status
	status_	= OPENSMOKE_DIRECTSOLVER_STATUS_FACTORIZED;
}

void OpenSMOKE_MUMPS_Unsymmetric::Solve(BzzVector &b, BzzVector &x)
{	
	if (status_	!= OPENSMOKE_DIRECTSOLVER_STATUS_FACTORIZED)
		ErrorMessage("Solve(BzzVector &b, BzzVector &x) can be used only for Factorized Matrices");

	// Checking dimensions
	if (b.Size() != x.Size())	ErrorMessage("The size of b and x do not fit!");
	if (b.Size() != n)			ErrorMessage("The size of matrix and rhs do not fit!");

	// Linear system data
	nrhs  = 1;						// Number of RHSs
	id.job = MUMPS_JOB_SOLUTION;	// Solve

	// Memory allocation
	double* b_ = new double[n];

	// From BzzVector to C vector
	for(int j=0;j<n;j++)
		b_[j] = b[j+1];

	id.rhs = b_;
	
	MessageOnTheScreen("MUMPS: Solve sparse linear system");
	dmumps_c(&id);
	MessageOnTheScreen("MUMPS: Done");
	ErrorAnalysis();

	// From C Vector to BzzVector
	for(int j=0;j<n;j++)
		x[j+1] = id.rhs[j];
}

void OpenSMOKE_MUMPS_Unsymmetric::Solve(BzzMatrix &b, BzzMatrix &x, const bool iVerticalStorage)
{	
	ErrorMessage("Solve(BzzMatrix &b, BzzMatrix &x, const bool iVerticalStorage) not yet available");
}

void OpenSMOKE_MUMPS_Unsymmetric::ResetCounters()
{
	countGlobal_ = 0;
	countGlobalRows_ = 0;
}

void OpenSMOKE_MUMPS_Unsymmetric::Delete()
{
	CleanMemory();

	id.job=MUMPS_JOB_END;
	MessageOnTheScreen("MUMPS: Delete");
	dmumps_c(&id);
	MessageOnTheScreen("MUMPS: Done");
}

void OpenSMOKE_MUMPS_Unsymmetric::UnsetMessage()
{
}

void OpenSMOKE_MUMPS_Unsymmetric::ErrorAnalysis()
{
	error_ = id.info[0];

	if (error_ == 0)	return;
	else if (error_ == -1)	ErrorMessage("An error occurred on processor");
	else if (error_ == -2)	ErrorMessage("NZ is out of range");
	else if (error_ == -3)	ErrorMessage("MUMPS was called with an invalid value for JOB");
	else if (error_ == -4)	ErrorMessage("Error in user-provided permutation array PERM IN at position INFO(2)");
	else if (error_ == -5)	ErrorMessage("Problem of REAL or COMPLEX workspace allocation of size INFO(2) during analysis");
	else if (error_ == -6)	ErrorMessage("Matrix is singular in structure. INFO(2) holds the structural rank");
	else if (error_ == -7)	ErrorMessage("Problem of INTEGER workspace allocation of size INFO(2) during analysis");
	else if (error_ == -8)	ErrorMessage("Main internal integer workarray IS too small for factorization");
	else if (error_ == -9)	ErrorMessage("Main internal real/complex workarray S too small");
	else if (error_ == -10)	ErrorMessage("Numerically singular matrix");
	else if (error_ == -11)	ErrorMessage("Internal real/complex workarray S too small for solution");
	else if (error_ == -12)	ErrorMessage("Internal real/complex workarray S too small for iterative refinement");
	else if (error_ == -13)	ErrorMessage("An error occurred in a Fortran ALLOCATE statement");
	else if (error_ == -14)	ErrorMessage("Internal integer workarray IS too small for solution");
	else if (error_ == -15)	ErrorMessage("Integer workarray IS too small for iterative refinement and/or error analysis");
	else if (error_ == -16)	ErrorMessage("N is out of range");
	else if (error_ == -17)	ErrorMessage("The internal send buffer that was allocated dynamically by MUMPS on the processor is too small");
	else if (error_ == -20)	ErrorMessage("The internal reception buffer that was allocated dynamically by MUMPS is too small");
	else if (error_ == -21)	ErrorMessage("Value of PAR=0 is not allowed because only one processor is available");
	else if (error_ == -22)	ErrorMessage("A pointer array is provided by the user that is either");

	else if (error_ == -23)	ErrorMessage("MPI was not initialized by the user prior to a call to MUMPS with JOB = -1");
	else if (error_ == -24)	ErrorMessage("NELT is out of range");
	else if (error_ == -25)	ErrorMessage("A problem has occurred in the initialization of the BLACS");
	else if (error_ == -26)	ErrorMessage("LRHS is out of range");
	else if (error_ == -27)	ErrorMessage("NZ RHS and IRHS PTR(NRHS+1) do not match");
	else if (error_ == -28)	ErrorMessage("IRHS PTR(1) is not equal to 1");
	else if (error_ == -29)	ErrorMessage("LSOL loc is smaller than INFO(23)");
	else if (error_ == -30)	ErrorMessage("CHUR LLD is out of range");
	else					ErrorMessage("TODO");
}
