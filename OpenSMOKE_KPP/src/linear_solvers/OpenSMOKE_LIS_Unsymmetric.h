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
#include <petscvec.h>
#include "OpenSMOKE_KPP_Definitions.h"
#include "linear_solvers/OpenSMOKE_DirectLinearSolver_Unsymmetric.h"
#include <vector>

class OpenSMOKE_KPP_DataManager;
class OpenSMOKE_KPP_ReactorNetwork;
class OpenSMOKE_KPP_Communicator;

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
	inline long long int localnonzero()		const { return local_nonzero;}
	inline int* countLocal()			const { return countLocal_; }
	inline int* local_rows()			const { return rows_per_worker;}
	inline int* lis_offset()			const { return offset_;}

	void SetData(OpenSMOKE_KPP_DataManager* data);
	void SetCommunicator(OpenSMOKE_KPP_Communicator* communicator);
	
	// 1b. Set from BzzVector + Symbolic Factorization
	void OpenMatrix(const int nRows, const long long int numberNonZeroElements);
	void SetSparsityPattern(const int nElementsPerRows, const BzzVectorInt &indicesColumns);
	void SetSparsityPattern(const int index, const int blockDimension, const BzzVector& mConvectionDiffusion, const BzzVectorInt& iConvectionDiffusion, const int nElementsPerRows);
	void SetSparsityPattern(const int index, const int blockDimension, const std::vector<double*>& mConvDiff, const std::vector<int*>& iConvDiff, int*& nElementsPerRows, OpenSMOKE_KPP_ReactorNetwork& network_);
	void SetLocalSparsityPattern(const int index, const int blockDimension, const std::vector<double*>& mConvDiff, const std::vector<int*>& iConvDiff, int*& nElementsPerRows, OpenSMOKE_KPP_ReactorNetwork& network_);
        void GetBooleanSparsityPattern(const int numberofrows, const int numberofreactors);
	void CloseMatrix();
	void UpdateMatrix(const BzzVector &valuesVector);
	void UpdateMatrix(const BzzVector &valuesVector, int position);
	void UpdateMatrix(const Vec &valuesVector);
	void UpdateMatrix(const double* values_array);
	void ResetCounters();
	void ResetAllCounters();
	void MatrixDistribution();
	void ValuesDistribution();
	void CompleteMatrix();
	void UpdateVectors(BzzVector &b, BzzVector &x);
	void SetLocalRows(OpenSMOKE_KPP_ReactorNetwork& network_);
	
	// 2. Numerical Factorization and Numerical Solution
	void NumericalFactorization();
	void SetInitialSolution(BzzVector &firstguess);
	void Solve(BzzVector &b, BzzVector &x);
	void Solve(Vec &b, BzzVector &x);
	void Solve(double* &b, BzzVector &x);
	void Solve(BzzMatrix &b, BzzMatrix &x, const bool iVerticalStorage);
	inline bool convergence_index() const {return convergence_index_;}

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
	void MemoryAllocation();	
	void MemoryAllocation(OpenSMOKE_KPP_ReactorNetwork& network_);

public:

	// Options (Set/Unset)
	void SetDefaultOptions();
	void UnsetMessage();

	// Get data
	LIS_MATRIX        A;
	LIS_VECTOR        lis_b, lis_x;
	LIS_SOLVER        solver;
	long long int *border_rows;

private:

	void ErrorAnalysis();
	void SetUserDefaultOptions(const OpenSMOKE_DirectLinearSolver_Unsymmetric_Kind kind);
	void UpdateOptions();
	void MassFlowRateUpdateOptions();
	void ParallelMemoryAllocation();
	void SetLocalRows();
	void InitializePlaceVectors();

	bool convergence_index_;

	//Parallel utilities
	int nprocs_, numworkers_, procrank_, MASTER, FROM_MASTER, FROM_WORKER, mtype;
	int averows, extra;
	int proc_counter;
	int values_proccounter;
	int count_val;
	int k_stop;
	int* offset_;
	int *rows_per_worker, *rows_per_worker_offset;
	long long int *nonzero_offset;
	int nequations;
	int *countLocal_, *countLocalRows_;
	PetscInt petsc_nrows, petsc_nonzero;
	LIS_INT local_nrows;
	LIS_INT local_nonzero;
	LIS_INT local_startindex, local_endindex;
	double* proc_vals;
	int process_rows, process_values;
	bool iAllocation;
	bool iSend;

	LIS_SCALAR *loc_values;		//	Matrix local elements
	LIS_INT    *loc_rows;		//	Matrix local rows
	LIS_INT    *loc_columns;		//	Matrix local columns
	
	double *lis_rows, *lis_columns;
	double *x_;
	LIS_INT solver_status;
	std::vector<int*> proc_rows;
	std::vector<int*> proc_cols;
	std::vector<double*> proc_values;
        std::vector<bool*> sparsity_;

	Vec rows_vec, columns_vec, values_vec;
	PetscInt rows_lower, rows_upper, columns_lower, columns_upper, values_lower, values_upper;
	PetscScalar *rows_ptr, *columns_ptr, *values_ptr;
	Petsc64bitInt *rows_place, *columns_place, *values_place;


	std::string maxIter_;
	std::string tol_;
	std::string linEqSolv_;
	std::string lisPrecond_;

	std::string linearEquationSolver_;
	std::string lisPreconditioner_;
	std::string messageLevel_;
	std::string initialSolution_;
	std::string maximumIterations_;
	std::string tolerance_;

private:

	OpenSMOKE_KPP_DataManager* data_;
	OpenSMOKE_KPP_Communicator* communicator_;
	
};

#endif	// OpenSMOKE_PARDISO_Unsymmetric_H
