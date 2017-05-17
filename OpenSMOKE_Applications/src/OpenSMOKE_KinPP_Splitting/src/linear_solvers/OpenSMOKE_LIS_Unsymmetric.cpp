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

#include "lis.h"
#include "mkl.h"
#include <iomanip>
#include <sstream>
#include "OpenSMOKE_LIS_Unsymmetric.h"

OpenSMOKE_LIS_Unsymmetric::OpenSMOKE_LIS_Unsymmetric(void)
{
	ErrorMessage("No default constructor is available");
}

OpenSMOKE_LIS_Unsymmetric::OpenSMOKE_LIS_Unsymmetric(const OpenSMOKE_DirectLinearSolver_Unsymmetric_Kind kind)
{
	SetUserDefaultOptions(kind);
}

void OpenSMOKE_LIS_Unsymmetric::SetUserDefaultOptions(const OpenSMOKE_DirectLinearSolver_Unsymmetric_Kind kind)
{
	kind_						= kind;
	status_						= OPENSMOKE_DIRECTSOLVER_STATUS_OPEN;

	maximumIterations_			= " -maxiter 1000 ";
	tolerance_				= " -tol 1.e-12 ";
	messageLevel_				= " -print out ";
	initialSolution_			= " -initx_zeros true ";
	linearEquationSolver_			= " -i bicg ";
	lisPreconditioner_			= " -p ilu ";
}

OpenSMOKE_LIS_Unsymmetric::~OpenSMOKE_LIS_Unsymmetric(void)
{
}

void OpenSMOKE_LIS_Unsymmetric::SetSparsityPattern(BzzMatrixSparse &M)
{
	if (kind_ != OPENSMOKE_DIRECTSOLVER_SQUAREMATRIX)
		ErrorMessage("SetSparsityPattern(BzzMatrixSparse &M) can be used only for OPENSMOKE_DIRECTSOLVER_SQUAREMATRIX");

	if (status_	!= OPENSMOKE_DIRECTSOLVER_STATUS_OPEN)
		ErrorMessage("SetSparsityPattern(BzzMatrixSparse &M) cannot be used if the Matrix was already closed");
	

	n = M.Rows();

	// Counting non zero elements
	{
		double* ptrVal;
		int i, j;
		double val;

		numberNonZeroElements = 0;
		M.BeginScanning();
		while(ptrVal = M.Scanning(&i,&j,&val))
			numberNonZeroElements++;
	}

	// Memory Allocation
	lis_matrix_malloc_csr(n,numberNonZeroElements,&rows,&columns,&values);

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

void OpenSMOKE_LIS_Unsymmetric::UpdateMatrix(BzzMatrixSparse &M)
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

void OpenSMOKE_LIS_Unsymmetric::OpenMatrix(const int nRows, const int numberNonZeroElements_)
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

	// Memory Allocation
	lis_matrix_malloc_csr(n,numberNonZeroElements,&rows,&columns,&values);

	// Counters
	countGlobal_=0;
	countGlobalRows_=0;
	rows[n] = numberNonZeroElements+1;
}

void OpenSMOKE_LIS_Unsymmetric::SetSparsityPattern(const int nElementsPerRows, const BzzVectorInt &indicesColumns)
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

void OpenSMOKE_LIS_Unsymmetric::SetSparsityPattern(const int index, const int blockDimension, const BzzVector& mConvectionDiffusion, const BzzVectorInt& iConvectionDiffusion, const int nElementsPerRows)
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


void OpenSMOKE_LIS_Unsymmetric::CloseMatrix()
{
	if (kind_ != OPENSMOKE_DIRECTSOLVER_SQUAREMATRIX_ROW_BY_ROW)
		ErrorMessage("CloseMatrix() can be used only for OPENSMOKE_DIRECTSOLVER_SQUAREMATRIX_ROW_BY_ROW");

	if (status_	!= OPENSMOKE_DIRECTSOLVER_STATUS_OPEN)
		ErrorMessage("CloseMatrix() can be used only for Open Matrices");

	CompleteMatrix();

	// Update status
	status_ = OPENSMOKE_DIRECTSOLVER_STATUS_TOFACTORIZE;
}

void OpenSMOKE_LIS_Unsymmetric::UpdateMatrix(const BzzVector &valuesVector)
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

void OpenSMOKE_LIS_Unsymmetric::CompleteMatrix()
{
	if (status_	!= OPENSMOKE_DIRECTSOLVER_STATUS_OPEN)
		ErrorMessage("CompleteMatrix() cannot be used if the Matrix was already closed");

	if (status_	!= OPENSMOKE_DIRECTSOLVER_STATUS_OPEN)
		ErrorMessage("CompleteMatrix() cannot be used if the Matrix was already closed");

	cout << " Number of non-zero elements:    " << numberNonZeroElements << endl;
	cout << " Sparsity fill-in:               " << double(numberNonZeroElements)/double(n*n)*100. << endl;
	cout << " Mean non zero elements per row: " << double(numberNonZeroElements)/double(n) << endl;
	
	// Indices revision (C-Style)
	for(int j=1;j<=n+1;j++)
		rows[j-1] -= 1;
	for(int j=1;j<=numberNonZeroElements;j++)
		columns[j-1] -= 1;
	
	// Assembling
	lis_matrix_create(0,&A);
	lis_matrix_set_size(A,0,n);
	lis_matrix_set_csr(numberNonZeroElements, rows, columns, values, A);
	lis_matrix_assemble(A);

	// Vectors
	lis_vector_create(0, &lis_b);
	lis_vector_create(0, &lis_x);
	lis_vector_set_size(lis_b,n,0);
	lis_vector_set_size(lis_x,n,0);

	// Update Status
	status_	= OPENSMOKE_DIRECTSOLVER_STATUS_TOFACTORIZE;
}

void OpenSMOKE_LIS_Unsymmetric::NumericalFactorization()
{
	if (status_	== OPENSMOKE_DIRECTSOLVER_STATUS_OPEN)
		ErrorMessage("NumericalFactorization() cannot be used for open matrices");

	// Update Status
	status_	= OPENSMOKE_DIRECTSOLVER_STATUS_FACTORIZED;
}

void OpenSMOKE_LIS_Unsymmetric::Solve(BzzVector &b, BzzVector &x)
{	
	if (status_	!= OPENSMOKE_DIRECTSOLVER_STATUS_FACTORIZED)
		ErrorMessage("Solve(BzzVector &b, BzzVector &x) can be used only for Factorized Matrices");

		// Checking dimensions
	if (b.Size() != x.Size())	ErrorMessage("The size of b and x do not fit!");
	if (b.Size() != n)			ErrorMessage("The size of matrix and rhs do not fit!");
	
	// From BzzVector to LIS Vector
	double* b_ = b.GetHandle();
	for(int i=0;i<n;i++)
		lis_vector_set_value(LIS_INS_VALUE, i, b_[i], lis_b);

	// Set options
	char message[200];
	string option_string =	linearEquationSolver_ + lisPreconditioner_ + maximumIterations_ + 
							tolerance_ + messageLevel_ + initialSolution_;

	lis_solver_create(&solver);
	strcpy(message, option_string.c_str());
//	lis_solver_set_option("-i gmres -p ilu -print out -initx_zeros false",solver);
//	lis_solver_set_option("-p ilu -print out -initx_zeros false",solver);
	lis_solver_set_option("-p ilut -i bicgsafe -print out -initx_zeros false -maxiter 5000 -tol 1e-13 ",solver);

	// Solve linear system
	error_ = lis_solve(A,lis_b,lis_x,solver);
	ErrorAnalysis();

	// From LIS Vector to BzzVector
	double* x_ = x.GetHandle();
	for(int i=0;i<n;i++)
		lis_vector_get_value(lis_x, i, &x_[i]);
}

void OpenSMOKE_LIS_Unsymmetric::Solve(BzzMatrix &b, BzzMatrix &x, const bool iVerticalStorage)
{	
	if (status_	!= OPENSMOKE_DIRECTSOLVER_STATUS_FACTORIZED)
		ErrorMessage("Solve(BzzMatrix &b, BzzMatrix &x, const bool iVerticalStorage) can be used only for Factorized Matrices");
}

void OpenSMOKE_LIS_Unsymmetric::SetInitialSolution(BzzVector &firstguess)
{
	double* firstguess_ = firstguess.GetHandle();
	for(int i=0;i<n;i++)
		lis_vector_set_value(LIS_INS_VALUE, i, firstguess_[i], lis_x);
}

void OpenSMOKE_LIS_Unsymmetric::ResetCounters()
{
	countGlobal_ = 0;
	status_ = OPENSMOKE_DIRECTSOLVER_STATUS_TOFACTORIZE;
}

void OpenSMOKE_LIS_Unsymmetric::Delete()
{
	// Memory cleaning is managed by the following functions
	lis_solver_destroy(solver);
	lis_matrix_destroy(A);
	lis_vector_destroy(lis_b);
	lis_vector_destroy(lis_x);
}

void OpenSMOKE_LIS_Unsymmetric::SetDefaultOptions()
{
}

void OpenSMOKE_LIS_Unsymmetric::ErrorAnalysis()
{
}

void OpenSMOKE_LIS_Unsymmetric::UnsetMessage()
{
	messageLevel_ = "-print none";
}

void OpenSMOKE_LIS_Unsymmetric::SetLinearEquationSolver(const linearEquationSolver kind)
{
	if (kind == LINEAR_EQUATION_SOLVER_CG)
		linearEquationSolver_ = "-i cg";
	else if (kind == LINEAR_EQUATION_SOLVER_BiCG)
		linearEquationSolver_ = " -i bicg ";
	else if (kind == LINEAR_EQUATION_SOLVER_CGS)
		linearEquationSolver_ = " -i cgs ";
	else if (kind == LINEAR_EQUATION_SOLVER_BiCGSTAB)
		linearEquationSolver_ = " -i bicgstab ";
	else if (kind == LINEAR_EQUATION_SOLVER_BiCGSTABl)
		linearEquationSolver_ = " -i bicgstabl ";
	else if (kind == LINEAR_EQUATION_SOLVER_GPBiCG)
		linearEquationSolver_ = " -i gpbicg ";
	else if (kind == LINEAR_EQUATION_SOLVER_TFQMR)
		linearEquationSolver_ = " -i tfqmr ";
	else if (kind == LINEAR_EQUATION_SOLVER_Orthomin)
		linearEquationSolver_ = " -i orthomin ";
	else if (kind == LINEAR_EQUATION_SOLVER_GMRES)
		linearEquationSolver_ = " -i gmres ";
	else if (kind == LINEAR_EQUATION_SOLVER_Jacobi)
		linearEquationSolver_ = " -i jacobi ";
	else if (kind == LINEAR_EQUATION_SOLVER_GaussSeidel)
		linearEquationSolver_ = " -i gs ";
	else if (kind == LINEAR_EQUATION_SOLVER_SOR)
		linearEquationSolver_ = " -i sor ";
	else if (kind == LINEAR_EQUATION_SOLVER_BiCGSafe)
		linearEquationSolver_ = " -i bicgsafe ";
	else if (kind == LINEAR_EQUATION_SOLVER_CR)
		linearEquationSolver_ = " -i cr ";
	else if (kind == LINEAR_EQUATION_SOLVER_BiCR)
		linearEquationSolver_ = " -i bicr ";
	else if (kind == LINEAR_EQUATION_SOLVER_CRS)
		linearEquationSolver_ = " -i csr ";
	else if (kind == LINEAR_EQUATION_SOLVER_BiCRSTAB)
		linearEquationSolver_ = "-i bicrstab";
	else if (kind == LINEAR_EQUATION_SOLVER_GPBiCR)
		linearEquationSolver_ = "-i gpbicr";
	else if (kind == LINEAR_EQUATION_SOLVER_BiCRSafe)
		linearEquationSolver_ = "-i bicrsafe";
	else if (kind == LINEAR_EQUATION_SOLVER_FGMRESm)
		linearEquationSolver_ = "-i fgmres";
	else if (kind == LINEAR_EQUATION_SOLVER_IDRs)
		linearEquationSolver_ = "-i idrs";
	else if (kind == LINEAR_EQUATION_SOLVER_MINRES)
		linearEquationSolver_ = "-i minres";
}

void OpenSMOKE_LIS_Unsymmetric::SetPreconditioner(const lisPreconditioner kind)
{
	if (kind == LIS_PRECONDITIONER_NONE)
		lisPreconditioner_ = "-p none";
	else if (kind == LIS_PRECONDITIONER_Jacobi)
		lisPreconditioner_ = "-p jacobi";
	else if (kind == LIS_PRECONDITIONER_ILUk)
		lisPreconditioner_ = " -p ilu ";
	else if (kind == LIS_PRECONDITIONER_SSOR)
		lisPreconditioner_ = " -p ssor ";
	else if (kind == LIS_PRECONDITIONER_Hybrid)
		lisPreconditioner_ = " -p hybrid ";
	else if (kind == LIS_PRECONDITIONER_IS)
		lisPreconditioner_ = " -p is ";
	else if (kind == LIS_PRECONDITIONER_SAINV)
		lisPreconditioner_ = " -p sainv ";
	else if (kind == LIS_PRECONDITIONER_SAAMG)
		lisPreconditioner_ = " -p saamg ";
	else if (kind == LIS_PRECONDITIONER_CROUTILU)
		lisPreconditioner_ = " -p iluc ";
	else if (kind == LIS_PRECONDITIONER_ILUT)
		lisPreconditioner_ = " -p ilut ";
}

void OpenSMOKE_LIS_Unsymmetric::SetMaximumIterations(const int maximumIterations)
{
	stringstream maxiter;
	maxiter << maximumIterations;
	maximumIterations_ = "-maxiter " + maxiter.str();
}

void OpenSMOKE_LIS_Unsymmetric::SetTolerance(const double tolerance)
{
	stringstream tol;
	tol << tolerance;
	tolerance_ = "-tol " + tol.str();
}

void OpenSMOKE_LIS_Unsymmetric::SetMessage()
{
	messageLevel_ = "-print out";
}

void OpenSMOKE_LIS_Unsymmetric::DisableInitialSolution()
{
	initialSolution_ = "-initx_zeros true";
}

void OpenSMOKE_LIS_Unsymmetric::EnableInitialSolution()
{
	initialSolution_ = "-initx_zeros false";
}

