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
#include <mpi.h>
#include <petscvec.h>
#include "OpenSMOKE_LIS_Unsymmetric.h"
#include "../kpp/OpenSMOKE_KPP_ReactorNetwork.h"
#include "../kpp/OpenSMOKE_KPP_DataManager.h"
#include "../kpp/OpenSMOKE_KPP_Communicator.h"
#include <vector>

OpenSMOKE_LIS_Unsymmetric::OpenSMOKE_LIS_Unsymmetric(void)
{
	ErrorMessage("No default constructor is available");
}

OpenSMOKE_LIS_Unsymmetric::OpenSMOKE_LIS_Unsymmetric(const OpenSMOKE_DirectLinearSolver_Unsymmetric_Kind kind)
{
	SetUserDefaultOptions(kind);

	nprocs_ = MPI::COMM_WORLD.Get_size();
	numworkers_ = nprocs_ - 1;
	procrank_ = MPI::COMM_WORLD.Get_rank();
	MASTER = 0;
	FROM_MASTER = 0;
	FROM_WORKER = 1;

	countLocal_ = new int[nprocs_];
	countLocalRows_ = new int[nprocs_];

	convergence_index_ = false;
}

void OpenSMOKE_LIS_Unsymmetric::SetData(OpenSMOKE_KPP_DataManager* data)
{
	data_ = data;
}

void OpenSMOKE_LIS_Unsymmetric::SetUserDefaultOptions(const OpenSMOKE_DirectLinearSolver_Unsymmetric_Kind kind)
{
	kind_						= kind;
	status_						= OPENSMOKE_DIRECTSOLVER_STATUS_OPEN;
	
	maxIter_				= " -maxiter ";
	tol_					= " -tol ";
	linEqSolv_				= " -i ";
	lisPrecond_				= " -p ";
	
	maximumIterations_			= "1000";
	linearEquationSolver_			= "bicgsafe";
	lisPreconditioner_			= "ilu";
	tolerance_				= "1e-20";

	messageLevel_				= " -print out ";
	initialSolution_			= " -initx_zeros false ";
	
/*	
	maximumIterations_			= " -maxiter 1000 ";
	tolerance_					= " -tol 1.e-12 ";
	linearEquationSolver_			= " -i ssor ";
	lisPreconditioner_			= " -p ilu ";*/
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

	SetLocalRows();

	// Counting non zero elements
	if(procrank_ == 0)
	{
	    double* ptrVal;
	    int i, j;
	    double val;
	    int proc_counter = 0;

	    numberNonZeroElements = 0;

	    M.BeginScanning();
	    while(ptrVal = M.Scanning(&i,&j,&val))
            {
                if(i == offset_[proc_counter + 1] + 1 && proc_counter != numworkers_)
                {
                    proc_counter++;
                }

		countLocal_[proc_counter]++;
		numberNonZeroElements++;
            }
	}

    	MPI::COMM_WORLD.Bcast(&numberNonZeroElements, 1, MPI::LONG_LONG, MASTER);
    	MPI::COMM_WORLD.Bcast(&countLocal_[0], nprocs_, MPI::INT, MASTER);

	// Memory Allocation
	MemoryAllocation();

	lis_matrix_malloc_csr(rows_per_worker[procrank_],countLocal_[procrank_],&loc_rows,&loc_columns,&loc_values);

	if(procrank_ == 0)
	{
	    lis_rows[0] = 1;
	    // Counting non zero elements
	    {
	        double* ptrVal;
                int i, j;
                double val;
                int count = 0;
	        int iRows = 1;
	
	        ParallelMemoryAllocation();

	        M.BeginScanning();
	        while(ptrVal = M.Scanning(&i,&j,&val))
	        {
		    if (iRows != i)
		    {
		        lis_rows[iRows] = count+1;
		        iRows = i;
		    }
		        count++;
	        }
		    lis_rows[n] = count+1;
		    numberNonZeroElements = count;

		delete [] ptrVal;

		rows_place = new Petsc64bitInt[n+1];

		for(int i = 0; i <= n; i++)
		    rows_place[i] = i;

		VecSetValues(rows_vec, n + 1, rows_place, lis_rows, INSERT_VALUES);

		delete [] rows_place;

		for(int i = 0; i <= numworkers_; i++)
		    border_rows[i] = lis_rows[offset_[i]];
	    }	
	}

	VecAssemblyBegin(rows_vec);
	VecAssemblyEnd(rows_vec);
	VecGetArray(rows_vec, &rows_ptr);
	
	MPI::COMM_WORLD.Bcast(&border_rows[0], nprocs_, MPI::LONG_LONG, MASTER);

	int counter = 0;
	if(procrank_ == 0) counter = 1;
	loc_rows[0] = 1;
	for(int j = 1; j <= rows_per_worker[procrank_]; j++)
	{
	    loc_rows[j] = rows_ptr[counter] - border_rows[procrank_] + 1;
	    counter++;
	}

	VecRestoreArray(rows_vec, &rows_ptr);

	
        if(procrank_ == 0)
	{
	    values = new double[numberNonZeroElements];
	    // Filling non zero elements
	    {	
		double* ptrVal;
		int i, j;
		double val;
		int count = 0;

		M.BeginScanning();
		while(ptrVal = M.Scanning(&i,&j,&val))
		{
			lis_columns[count] = j;
			values[count]  = val;
			count++;
		}
	    }

	    columns_place = new Petsc64bitInt[numberNonZeroElements];
	    values_place = new Petsc64bitInt[numberNonZeroElements];

	    for(int i = 0; i < numberNonZeroElements; i++)
	    {
		columns_place[i] = i;
		values_place[i] = i;
	    }

	    VecSetValues(columns_vec, numberNonZeroElements, columns_place, lis_columns, INSERT_VALUES);
	    VecSetValues(values_vec, numberNonZeroElements, values_place, values, INSERT_VALUES);

	    delete [] columns_place;
	    delete [] values_place;
	    delete [] values;

	}

	VecAssemblyBegin(columns_vec);
	VecAssemblyBegin(values_vec);
	VecAssemblyEnd(columns_vec);
	VecAssemblyEnd(values_vec);
	VecGetArray(columns_vec, &columns_ptr);
	VecGetArray(values_vec, &values_ptr);

	for(int i = 0; i < countLocal_[procrank_]; i++)
	{
	    loc_columns[i] = columns_ptr[i];
	    loc_values[i] = values_ptr[i];
	}

	// Update Status
//	status_	= OPENSMOKE_DIRECTSOLVER_STATUS_TOFACTORIZE;
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

void OpenSMOKE_LIS_Unsymmetric::OpenMatrix(const int nRows, const long long int numberNonZeroElements_)
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

	// Counters
	countGlobal_=0;
	countGlobalRows_=0;

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

void OpenSMOKE_LIS_Unsymmetric::SetSparsityPattern(const int index, const int blockDimension, const std::vector<double*>& mConvDiff, const std::vector<int*>& iConvDiff, int*& nElementsPerRows, OpenSMOKE_KPP_ReactorNetwork& network_)
{
	if (kind_ != OPENSMOKE_DIRECTSOLVER_SQUAREMATRIX_ROW_BY_ROW)
		ErrorMessage("SetSparsityPattern(const int nElementsPerRows, const BzzVectorInt &indicesColumns) can be used only for OPENSMOKE_DIRECTSOLVER_SQUAREMATRIX_ROW_BY_ROW");

	if (status_	!= OPENSMOKE_DIRECTSOLVER_STATUS_OPEN)
		ErrorMessage("SetSparsityPattern(const int nElementsPerRows, const BzzVectorInt &indicesColumns) can be used only for Open Matrices");

	BzzVectorInt indicesColumns;
	CalculateSparsityPattern(index, blockDimension, mConvDiff, iConvDiff, indicesColumns, network_);

	// Non zero elements
	{		
	    int nLocalRows = indicesColumns.Size() / nElementsPerRows[index];
	    for(int k=1;k<=nLocalRows;k++)
	    {
		loc_rows[countGlobalRows_++] = countGlobal_ + 1 + nElementsPerRows[index]*(k-1) + nonzero_offset[procrank_];
	    }

	    for(int k=1;k<=indicesColumns.Size();k++)
	    {
		loc_columns[countGlobal_] = indicesColumns[k];
		countGlobal_++;
	    }
	}

}

void OpenSMOKE_LIS_Unsymmetric::SetLocalSparsityPattern(const int index, const int blockDimension, const std::vector<double*>& mConvDiff, const std::vector<int*>& iConvDiff, int*& nElementsPerRows, OpenSMOKE_KPP_ReactorNetwork& network_)
{
	if (kind_ != OPENSMOKE_DIRECTSOLVER_SQUAREMATRIX_ROW_BY_ROW)
		ErrorMessage("SetSparsityPattern(const int nElementsPerRows, const BzzVectorInt &indicesColumns) can be used only for OPENSMOKE_DIRECTSOLVER_SQUAREMATRIX_ROW_BY_ROW");

	if (status_	!= OPENSMOKE_DIRECTSOLVER_STATUS_OPEN)
		ErrorMessage("SetSparsityPattern(const int nElementsPerRows, const BzzVectorInt &indicesColumns) can be used only for Open Matrices");

	BzzVectorInt indicesColumns;

	CalculateSparsityPattern(index, blockDimension, mConvDiff, iConvDiff, indicesColumns, network_);

	//Local Non zero elements
	{		
		int nLocalRows = indicesColumns.Size() / nElementsPerRows[index];


		for(int k=0;k<indicesColumns.Size();k++)
		{
		    if(k % nElementsPerRows[index] == 0)
			countGlobalRows_++;

		    if(countGlobalRows_ == offset_[proc_counter + 1] + 1 && proc_counter != numworkers_)
			proc_counter++;

			countLocal_[proc_counter]++;
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
	    for(int k = 0 ;k < local_nonzero;k++)
	    {
		loc_values[k]  = valuesVector[k+1];
	    }
	}

	// Update status
	status_ = OPENSMOKE_DIRECTSOLVER_STATUS_TOFACTORIZE;
}

void OpenSMOKE_LIS_Unsymmetric::UpdateMatrix(const BzzVector &valuesVector, int position)
{
	if (kind_ != OPENSMOKE_DIRECTSOLVER_SQUAREMATRIX_ROW_BY_ROW)
		ErrorMessage("UpdateMatrix(const BzzVector &valuesVector) can be used only for OPENSMOKE_DIRECTSOLVER_SQUAREMATRIX_ROW_BY_ROW");

	if (status_	== OPENSMOKE_DIRECTSOLVER_STATUS_OPEN)
		ErrorMessage("UpdateMatrix(const BzzVector &valuesVector) can be used only for Open Matrices");

	// Non zero elements
	{		
		for(int k=1;k<=valuesVector.Size();k++)
		{
			values[k - 1 + position]  = valuesVector[k];
		}
	}

	// Update status
	status_ = OPENSMOKE_DIRECTSOLVER_STATUS_TOFACTORIZE;
}

void OpenSMOKE_LIS_Unsymmetric::UpdateMatrix(const Vec &valuesVector)
{
	if (kind_ != OPENSMOKE_DIRECTSOLVER_SQUAREMATRIX_ROW_BY_ROW)
		ErrorMessage("UpdateMatrix(const BzzVector &valuesVector) can be used only for OPENSMOKE_DIRECTSOLVER_SQUAREMATRIX_ROW_BY_ROW");

	if (status_	== OPENSMOKE_DIRECTSOLVER_STATUS_OPEN)
		ErrorMessage("UpdateMatrix(const BzzVector &valuesVector) can be used only for Open Matrices");

	// Non zero elements

	PetscScalar *value_pointer;
	VecGetArray(valuesVector, &value_pointer);

	{		
	    for(int k = 0 ;k < countLocal_[procrank_];k++)
	    {
		loc_values[k]  = value_pointer[k];
	    }
	}

	VecRestoreArray(valuesVector, &value_pointer);

	// Update status
	status_ = OPENSMOKE_DIRECTSOLVER_STATUS_TOFACTORIZE;
}

void OpenSMOKE_LIS_Unsymmetric::UpdateMatrix(const double* values_array)
{
	if (kind_ != OPENSMOKE_DIRECTSOLVER_SQUAREMATRIX_ROW_BY_ROW)
		ErrorMessage("UpdateMatrix(const BzzVector &valuesVector) can be used only for OPENSMOKE_DIRECTSOLVER_SQUAREMATRIX_ROW_BY_ROW");

	if (status_	== OPENSMOKE_DIRECTSOLVER_STATUS_OPEN)
		ErrorMessage("UpdateMatrix(const BzzVector &valuesVector) can be used only for Open Matrices");

	// Non zero elements

	{		
	    for(int k = 0 ;k < local_nonzero;k++)
	    {
		loc_values[k]  = values_array[k];
	    }
	}

	// Update status
	status_ = OPENSMOKE_DIRECTSOLVER_STATUS_TOFACTORIZE;
}

/*void OpenSMOKE_LIS_Unsymmetric::UpdateAndDistributeMatrix(const BzzVector &valuesVector, int position)
{
	if (kind_ != OPENSMOKE_DIRECTSOLVER_SQUAREMATRIX_ROW_BY_ROW)
		ErrorMessage("UpdateMatrix(const BzzVector &valuesVector) can be used only for OPENSMOKE_DIRECTSOLVER_SQUAREMATRIX_ROW_BY_ROW");

	if (status_	== OPENSMOKE_DIRECTSOLVER_STATUS_OPEN)
		ErrorMessage("UpdateMatrix(const BzzVector &valuesVector) can be used only for Open Matrices");


	if(values_proccounter == 0)	
	{	
	    if(procrank_ == 0)
	    {
	        for(int k=1;k<=valuesVector.Size();k++)
	        {
		    loc_values[countGlobal_]  = valuesVector[k];
		    countGlobal_++;
		    if(countGlobal_ == countLocal_[0])
		    {
		        values_proccounter++;
		        break;
		    }
	        }
	    }
	}

	else
	{
	    if(procrank_ == 0)
	    {		
		for(int k=1;k<=valuesVector.Size();k++)
		{
		    proc_vals[count_val]  = valuesVector[k];
		    count_val++;
		    if(count_val == countLocal_[values_proccounter])
		    {
			iSend = true;
			k_stop = k;
		        break;
		    }
		}
		mtype = FROM_MASTER;

		if(iSend == true)
		{
		    MPI::COMM_WORLD.Send(&proc_vals[0], countLocal_[values_proccounter], MPI::DOUBLE, values_proccounter, mtype);
		    values_proccounter++;	
		    iSend = false;
		    iAllocation = false;

		    delete [] proc_vals;
		    proc_vals = new double[countLocal_[values_proccounter]];
		    count_val = 0;
		    iAllocation = true;

		    for(int k = k_stop; k <= valuesVector.Size(); k++)
		    {
			proc_vals[count_val]  = valuesVector[k];
		        count_val++;
		    }

		    iSend = false;
		}
	    }

	    if(procrank_ == values_proccounter)
	    {
		MPI::COMM_WORLD.Recv(&loc_values[0], countLocal_[values_proccounter], MPI::DOUBLE, MASTER, mtype, status);
		values_proccounter++;	
		iSend = false;
		iAllocation = false;
	    }

	    
	}

	// Update status
	status_ = OPENSMOKE_DIRECTSOLVER_STATUS_TOFACTORIZE;
}*/

void OpenSMOKE_LIS_Unsymmetric::ValuesDistribution()
{
	MPI::Status status;

	//Distributing local values
	mtype = FROM_MASTER;
	int counter = countLocal_[0];
	for(int p = 1; p <= numworkers_; p++)
    	{
            MPI::COMM_WORLD.Barrier();
            if(procrank_ == 0)
            {
                int source = MASTER;
                MPI::COMM_WORLD.Send(&values[counter], countLocal_[p], MPI::DOUBLE, p, mtype);

		counter += countLocal_[p];

		if(p == numworkers_ && counter != numberNonZeroElements)
		{
		    std::cout << countGlobal_ << "   " << counter << "      " << numberNonZeroElements << std::endl;
		    ErrorMessage("Incorrect Distribution");
		}

            }
            if(p == procrank_)
            {
                MPI::COMM_WORLD.Recv(&loc_values[0], countLocal_[p], MPI::DOUBLE, MASTER, mtype, status);
            }
        }

        if(procrank_ == 0)
        {
            for(int i = 0; i < countLocal_[0]; i++)
            {
                loc_values[i] = values[i];
            }
        }
	
	if(procrank_ == 0)
	    delete [] values;
		
}

void OpenSMOKE_LIS_Unsymmetric::CompleteMatrix()
{
	if (status_	!= OPENSMOKE_DIRECTSOLVER_STATUS_OPEN)
		ErrorMessage("CompleteMatrix() cannot be used if the Matrix was already closed");

	if (status_	!= OPENSMOKE_DIRECTSOLVER_STATUS_OPEN)
		ErrorMessage("CompleteMatrix() cannot be used if the Matrix was already closed");

	if(procrank_ == 0)
	{
	    std::cout << " Number of non-zero elements:    " << numberNonZeroElements << std::endl;
	    std::cout << " Sparsity fill-in:               " << double(numberNonZeroElements)/double(n)/double(n)*100. << " % " << std::endl;
	    std::cout << " Mean non zero elements per row: " << double(numberNonZeroElements)/double(n) << std::endl;
	}

	if(kind_ == OPENSMOKE_DIRECTSOLVER_SQUAREMATRIX_ROW_BY_ROW)
	{

	    border_rows[procrank_] = loc_rows[0];
	    
	    for(int p = 0; p <= numworkers_; p++)
	    {
		MPI::COMM_WORLD.Bcast(&border_rows[p], 1, MPI::LONG_LONG, p);
	    }

	    if(procrank_ < numworkers_) loc_rows[rows_per_worker[procrank_]] = border_rows[procrank_ + 1];
	    if(procrank_ == numworkers_) loc_rows[rows_per_worker[procrank_]] = numberNonZeroElements + 1;

    	    for(int j = 0; j <= rows_per_worker[procrank_]; j++)
		loc_rows[j] -= border_rows[procrank_] - 1;


	    // Indices revision (C-Style)
	    for(int j=1;j<=rows_per_worker[procrank_] + 1;j++)
		loc_rows[j-1] -= 1;
	    for(int j=1;j<=countLocal_[procrank_];j++)
		loc_columns[j-1] -= 1;


	    // Assembling
	    lis_matrix_create(LIS_COMM_WORLD,&A);
	    lis_matrix_set_size(A,local_nrows,0);

	    lis_matrix_set_csr(local_nonzero, loc_rows, loc_columns, loc_values, A);	

	    lis_matrix_assemble(A);

	    // Vectors
	    lis_vector_create(LIS_COMM_WORLD, &lis_b);
	    lis_vector_create(LIS_COMM_WORLD, &lis_x);
	    lis_vector_set_size(lis_b,rows_per_worker[procrank_],0);
	    lis_vector_set_size(lis_x,rows_per_worker[procrank_],0);
	    lis_vector_get_range(lis_b, &local_startindex, &local_endindex);
	    lis_vector_get_range(lis_x, &local_startindex, &local_endindex);

	    VecDestroy(&rows_vec);
	    VecDestroy(&columns_vec);
	}

	else if(kind_ == OPENSMOKE_DIRECTSOLVER_SQUAREMATRIX)
	{
	    // Indices revision (C-Style)
	    for(int j=1;j<=rows_per_worker[procrank_] + 1;j++)
		loc_rows[j-1] -= 1;
	    for(int j=1;j<=countLocal_[procrank_];j++)
		loc_columns[j-1] -= 1;

	    // Assembling
	    lis_matrix_create(LIS_COMM_WORLD,&A);

	    lis_matrix_set_size(A,rows_per_worker[procrank_],0);

	    lis_matrix_set_csr(countLocal_[procrank_], loc_rows, loc_columns, loc_values, A);
	    lis_matrix_assemble(A);


	    // Vectors
	    lis_vector_create(LIS_COMM_WORLD, &lis_b);
	    lis_vector_create(LIS_COMM_WORLD, &lis_x);
	    lis_vector_set_size(lis_b,0,n);
	    lis_vector_set_size(lis_x,0,n);
	    lis_vector_get_range(lis_b, &local_startindex, &local_endindex);
	    lis_vector_get_range(lis_x, &local_startindex, &local_endindex);
	}

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

	MPI::Status status;

		// Checking dimensions
	if (b.Size() != x.Size())	ErrorMessage("The size of b and x do not fit!");
	if (b.Size() != n)			ErrorMessage("The size of matrix and rhs do not fit!");
	
	// From BzzVector to LIS Vector
	double* b_ = new double [n + 1];
	
	for(int i = 1; i <= n; i++)
	    b_[i] = b[i];

	double *localb_ = new double [rows_per_worker[procrank_] + 1];

	for(int i = 1; i <= rows_per_worker[procrank_]; i++)
	    localb_[i] = b_[i + offset_[procrank_]];

/*	for(int p = 1; p <= numworkers_; p++)
    	{
            MPI::COMM_WORLD.Barrier();
            if(procrank_ == 0)
            {
                int source = MASTER;
                MPI::COMM_WORLD.Send(&b_[1 + offset_[p]], rows_per_worker[p], MPI::DOUBLE, p, mtype);
            }
            if(p == procrank_)
            {
                MPI::COMM_WORLD.Recv(&localb_[1], rows_per_worker[p], MPI::DOUBLE, MASTER, mtype, status);
            }
        }

	if(procrank_ == 0)
	{
	    for(int i = 1; i <= rows_per_worker[0]; i++)
		localb_[i] = b[i];
	}*/


	for(int i = local_startindex; i < local_endindex; i++)
	    lis_vector_set_value(LIS_INS_VALUE, i, localb_[i + 1 - local_startindex], lis_b);



	//Update Options
	UpdateOptions();

	if(kind_ == OPENSMOKE_DIRECTSOLVER_SQUAREMATRIX)
	    MassFlowRateUpdateOptions();

	// Set options
	char message[300];
	std::string space = " ";
	std::string option_string =	linEqSolv_ + linearEquationSolver_ + space + lisPrecond_ + lisPreconditioner_ + space + maxIter_ + maximumIterations_ + space + tol_ + tolerance_ + messageLevel_ + initialSolution_;
	
	SetInitialSolution(x);

	lis_solver_create(&solver);
	strcpy(message, option_string.c_str());
	lis_solver_set_option(message,solver);

	// Solve linear system
	error_ = lis_solve(A,lis_b,lis_x,solver);


	lis_solver_get_status(solver, &solver_status);
	lis_solver_get_residualnorm(solver, &residual_norm_);

	if(solver_status == 4 && residual_norm_ > 1e-4 && data_->networkStatus() == KPP_NETWORK_STATUS_GLOBALODE)
	{
	    if(procrank_ == 0)
	    {
		std::cout << "Convergence failed with LIS. Residual norm = " << residual_norm_ << std::endl;
		std::cout << "Changing method..." << std::endl;
	    }
	    convergence_index_ = false;
	}

	else if(solver_status != 0 && solver_status != 4 && procrank_ == 0)
	{
	    convergence_index_ = false;
	    std::cout << "There were issues in convergence. Check LIS manual for more" << std::endl;
	    std::cout << "Error code = " << solver_status << std::endl;
	    getchar();
	    exit(-1);
	}

	else convergence_index_ = true;

	// From LIS Vector to BzzVector
	double* localx_ = new double[rows_per_worker[procrank_] + 1];

	for(int i=local_startindex;i<local_endindex;i++)
		lis_vector_get_value(lis_x, i, &localx_[i - local_startindex]);

	double* x_ = new double[n + 1];

	for(int p = 1; p <= numworkers_; p++)
	{
	    MPI::COMM_WORLD.Barrier();
	    int mtype = FROM_WORKER;
	    if(procrank_ == p)
	    {
		MPI::COMM_WORLD.Send(&localx_[0], rows_per_worker[p], MPI::DOUBLE, MASTER, mtype);
	    }
	    if(procrank_ == 0)
	    {
		MPI::COMM_WORLD.Recv(&x_[1 + offset_[p]], rows_per_worker[p], MPI::DOUBLE, p, mtype, status);
	    }
	}

	if(procrank_ == 0)
	{
	    for(int i = 0; i < rows_per_worker[0]; i++)
	    {
		x_[i+1] = localx_[i];
	    }
	}

	communicator_->BroadcastArray(x_, n);

	for(int i = 1; i <= n; i++)
	    x[i] = x_[i];

	delete [] x_;
	delete [] localx_;
	delete [] b_;
	delete [] localb_;

/*	double* x_ = x.GetHandle();
	for(int i=0;i<n;i++)
		lis_vector_get_value(lis_x, i, &x_[i]);*/
}

void OpenSMOKE_LIS_Unsymmetric::Solve(BzzMatrix &b, BzzMatrix &x, const bool iVerticalStorage)
{	
	if (status_	!= OPENSMOKE_DIRECTSOLVER_STATUS_FACTORIZED)
		ErrorMessage("Solve(BzzMatrix &b, BzzMatrix &x, const bool iVerticalStorage) can be used only for Factorized Matrices");
}

void OpenSMOKE_LIS_Unsymmetric::SetInitialSolution(BzzVector &firstguess)
{
//	double* firstguess_ = firstguess.GetHandle();
	MPI::Status status;
	mtype = FROM_MASTER;

	if(kind_ == OPENSMOKE_DIRECTSOLVER_SQUAREMATRIX)
	{
	    double* firstguess_ = new double[n + 1];

	    if(procrank_ == 0)
	    {
	        for(int i = 1; i <= n; i++)
	            firstguess_[i] = firstguess[i];
	    }

	    double *localfirstguess_ = new double [rows_per_worker[procrank_] + 1];

	    for(int p = 1; p <= numworkers_; p++)
    	    {
                MPI::COMM_WORLD.Barrier();
                if(procrank_ == 0)
                {
                    int source = MASTER;
                    MPI::COMM_WORLD.Send(&firstguess_[1 + offset_[p]], rows_per_worker[p], MPI::DOUBLE, p, mtype);
                }
                if(p == procrank_)
                {
                    MPI::COMM_WORLD.Recv(&localfirstguess_[1], rows_per_worker[p], MPI::DOUBLE, MASTER, mtype, status);
                }
            }

	    if(procrank_ == 0)
	    {
	        for(int i = 1; i <= rows_per_worker[0]; i++)
		    localfirstguess_[i] = firstguess[i];
	    }

	    for(int i=local_startindex;i < local_endindex; i++)
	        lis_vector_set_value(LIS_INS_VALUE, i, localfirstguess_[i + 1 - local_startindex], lis_x);


	    delete [] firstguess_;
	    delete [] localfirstguess_;
	}

	if(kind_ == OPENSMOKE_DIRECTSOLVER_SQUAREMATRIX_ROW_BY_ROW)
	{
	    double* firstguess_ = new double[n + 1];
	    if(procrank_ == 0)
	    {
	        for(int i = 1; i <= n; i++)
	            firstguess_[i] = firstguess[i];
	    }

	    MPI::COMM_WORLD.Bcast(&firstguess_[1], n, MPI::DOUBLE, MASTER);

	    for(int i=local_startindex;i < local_endindex; i++)
		lis_vector_set_value(LIS_INS_VALUE, i, firstguess_[i + 1 + offset_[procrank_] - local_startindex], lis_x);

	    delete [] firstguess_;
	}
}

void OpenSMOKE_LIS_Unsymmetric::ResetCounters()
{
	countGlobal_ = 0;
	if(kind_ == OPENSMOKE_DIRECTSOLVER_SQUAREMATRIX_ROW_BY_ROW)	countGlobal_ = nonzero_offset[procrank_];

/*	if(procrank_ == 0)
	    values = new double[numberNonZeroElements];*/

	if(kind_ != OPENSMOKE_DIRECTSOLVER_SQUAREMATRIX_ROW_BY_ROW) 
	    status_ = OPENSMOKE_DIRECTSOLVER_STATUS_TOFACTORIZE;
}

void OpenSMOKE_LIS_Unsymmetric::Delete()
{
	// Memory cleaning is managed by the following functions
	lis_solver_destroy(solver);
	lis_matrix_destroy(A);
	lis_vector_destroy(lis_b);
	lis_vector_destroy(lis_x);

	VecDestroy(&rows_vec);
	VecDestroy(&columns_vec);
	VecDestroy(&values_vec);

	if(procrank_ == 0)
	{
	    delete [] lis_rows;
	    delete [] lis_columns;
	}

	delete [] border_rows;

/*	if(procrank_ == 0)
	{
	    communicator_->DeleteVector(proc_values, nprocs_);
	    delete [] values;
	}*/
/*
	delete [] loc_rows;
	delete [] loc_columns;
	delete [] loc_values;*/
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
	std::stringstream maxiter;
	maxiter << maximumIterations;
	maximumIterations_ = "-maxiter " + maxiter.str();
}

void OpenSMOKE_LIS_Unsymmetric::SetTolerance(const double tolerance)
{
	std::stringstream tol;
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

void OpenSMOKE_LIS_Unsymmetric::UpdateOptions()
{
	if(data_->networkStatus() == KPP_NETWORK_STATUS_GLOBALODE)
	{
	    maximumIterations_			= data_->GlobalODE_LisMaxIterations();
	    linearEquationSolver_		= data_->GlobalODE_LisSolvingMethod();
	    lisPreconditioner_			= data_->GlobalODE_LisPreconditioner();
	    tolerance_				= data_->GlobalODE_LisConvCriteria();	
	}

	if(data_->networkStatus() == KPP_NETWORK_STATUS_GLOBALNLS)
	{
	    maximumIterations_			= data_->GlobalNLS_LisMaxIterations();
	    linearEquationSolver_		= data_->GlobalNLS_LisSolvingMethod();
	    lisPreconditioner_			= data_->GlobalNLS_LisPreconditioner();
	    tolerance_				= data_->GlobalNLS_LisConvCriteria();
	}
}

void OpenSMOKE_LIS_Unsymmetric::MassFlowRateUpdateOptions()
{
	maximumIterations_			= data_->massFlowRate_LisMaxIterations();
	linearEquationSolver_			= data_->massFlowRate_LisSolvingMethod();
	lisPreconditioner_			= data_->massFlowRate_LisPreconditioner();
	tolerance_				= data_->massFlowRate_LisConvCriteria();
}

void OpenSMOKE_LIS_Unsymmetric::SetCommunicator(OpenSMOKE_KPP_Communicator* communicator)
{
	communicator_ = communicator;
}


void OpenSMOKE_LIS_Unsymmetric::SetLocalRows(OpenSMOKE_KPP_ReactorNetwork& network_)
{
	communicator_->InitializeArray(rows_per_worker, numworkers_);
	communicator_->InitializeArray(offset_, numworkers_);


	for(int i = 0; i <= numworkers_; i++)
	{
	    rows_per_worker[i] = network_.LocalNumberOfReactors()[i] * network_.NumberOfSpecies();
            if(i != numworkers_)
                offset_[i+1] = offset_[i] + rows_per_worker[i];
	}

	local_startindex = offset_[procrank_];
	local_endindex = offset_[procrank_] + rows_per_worker[procrank_];
	
}


void OpenSMOKE_LIS_Unsymmetric::ResetAllCounters()
{
	countGlobal_ = 0;
	countGlobalRows_ = 0;
	proc_counter = 0;
	for(int i = 0; i < nprocs_; i++)
	{
	    countLocal_[i] = 0;
	    countLocalRows_[i] = 0;
	}
}

void OpenSMOKE_LIS_Unsymmetric::MatrixDistribution()
{
	
/*
	//Separating process' vectors
	if(procrank_ == 0)
	{
	    long long int counter = 0;
	
	    for(int i = 0; i <= numworkers_; i++)
	    {
	        proc_rows[i][0] = 1;
	        for(int j = 1; j <= rows_per_worker[i]; j++)
	        {
                    counter++;
                    proc_rows[i][j] = rows[counter] - rows[offset_[i]] + 1;
                }
	//	proc_rows[i][rows_per_worker[i]] = countLocal_[i] + 1;
            }
		if(counter != n)
		    ErrorMessage("Incorrect Distribution");

	    counter = 0;
	    for(int i = 0; i <= numworkers_; i++)
	    {
		for(int j = 0; j < countLocal_[i]; j++)
                {
                    proc_cols[i][j] = columns[counter];
                    counter++;
                }
	    }
	    

		if(counter != numberNonZeroElements)
		    ErrorMessage("Incorrect Distribution");
	}

	//Sending local vectors to matrices
	mtype = FROM_MASTER;
	for(int p = 1; p <= numworkers_; p++)
    	{
            MPI::COMM_WORLD.Barrier();
            if(procrank_ == 0)
            {
                int source = MASTER;
                MPI::COMM_WORLD.Send(&proc_rows[p][0], rows_per_worker[p] + 1, MPI::INT, p, mtype);
                MPI::COMM_WORLD.Send(&proc_cols[p][0], countLocal_[p], MPI::INT, p, mtype);
            }
            if(p == procrank_)
            {
                MPI::COMM_WORLD.Recv(&loc_rows[0], rows_per_worker[p] + 1, MPI::INT, MASTER, mtype, status);
                MPI::COMM_WORLD.Recv(&loc_columns[0], countLocal_[p], MPI::INT, MASTER, mtype, status);
            }
        }

        if(procrank_ == 0)
        {
            for(int i = 0; i <= rows_per_worker[0]; i++)
                loc_rows[i] = proc_rows[0][i];

            for(int i = 0; i < countLocal_[0]; i++)
            {
                loc_columns[i] = proc_cols[0][i];
            }

	    communicator_->DeleteVector(proc_rows, nprocs_);
	    communicator_->DeleteVector(proc_cols, nprocs_);
	    communicator_->DeleteVector(proc_values, nprocs_);
	    delete [] rows;
	    delete [] columns;
	}*/
	
}

void OpenSMOKE_LIS_Unsymmetric::ParallelMemoryAllocation()
{

	if(procrank_ == 0 && kind_ == OPENSMOKE_DIRECTSOLVER_SQUAREMATRIX)
	{
	    int* proc_rows_length = new int[nprocs_];
	    for(int k = 0; k < nprocs_; k++)
    		proc_rows_length[k] = rows_per_worker[k] + 1;
	}	
}

void OpenSMOKE_LIS_Unsymmetric::MemoryAllocation()
{
	if(kind_ == OPENSMOKE_DIRECTSOLVER_SQUAREMATRIX)
	{
	    local_nrows = rows_per_worker[procrank_];
	    if(procrank_ == 0)	local_nrows = rows_per_worker[procrank_] + 1;
	    local_nonzero = countLocal_[procrank_];

	    petsc_nrows = local_nrows;
	    petsc_nonzero = local_nonzero;

	    communicator_->InitializePetscVector(rows_vec, petsc_nrows, rows_lower, rows_upper);
	    communicator_->InitializePetscVector(columns_vec, petsc_nonzero, columns_lower, columns_upper, columns_place);
	    communicator_->InitializePetscVector(values_vec, petsc_nonzero, values_lower, values_upper);

	    border_rows = new long long int[nprocs_];

	    if(procrank_ == 0)
	    {
	        lis_rows = new double[n + 1];
	        lis_columns = new double[numberNonZeroElements];

	        lis_rows[n] = numberNonZeroElements+1;
	    }
	}
}

void OpenSMOKE_LIS_Unsymmetric::MemoryAllocation(OpenSMOKE_KPP_ReactorNetwork& network_)
{
	if(kind_ == OPENSMOKE_DIRECTSOLVER_SQUAREMATRIX_ROW_BY_ROW)
	{
	    x_ = new double[n];

	    local_nrows = rows_per_worker[procrank_];
//	    if(procrank_ == 0)	local_nrows = rows_per_worker[procrank_] + 1;
	    local_nonzero = countLocal_[procrank_];

	    petsc_nrows = local_nrows;
	    petsc_nonzero = local_nonzero;

	    communicator_->InitializePetscVector(rows_vec, petsc_nrows, rows_lower, rows_upper);
	    communicator_->InitializePetscVector(columns_vec, petsc_nonzero, columns_lower, columns_upper);

	    lis_rows = new double[rows_per_worker[procrank_] + 1];
	    lis_columns = new double[countLocal_[procrank_]];

	    border_rows = new long long int[nprocs_];
	    communicator_->InitializeArray(nonzero_offset, numworkers_);
	    communicator_->InitializeArray(rows_per_worker_offset, numworkers_);

	    rows_per_worker_offset[0] = 0;
	    for(int i = 0; i < numworkers_; i++)
	    {
		nonzero_offset[i+1] = nonzero_offset[i] + countLocal_[i];
		rows_per_worker_offset[i+1] = rows_per_worker_offset[i] + rows_per_worker[i];
	    }

//	    InitializePlaceVectors();

	    // Memory Allocation
	    lis_matrix_malloc_csr(rows_per_worker[procrank_],countLocal_[procrank_],&loc_rows,&loc_columns,&loc_values);
	    countGlobalRows_ = 0;
	    countGlobal_ = 0;
	}
}

void OpenSMOKE_LIS_Unsymmetric::SetLocalRows()
{
	communicator_->InitializeArray(rows_per_worker, numworkers_);
	communicator_->InitializeArray(offset_, numworkers_);
//	communicator_->InitializeArray(local_nonzero, numworkers_);

	averows = n / nprocs_;
	extra = n % nprocs_;

	for(int i = 0; i <= numworkers_; i++)
	{
	    rows_per_worker[i] = (i+1 <= extra) ? averows + 1 : averows;
            if(i != numworkers_)
                offset_[i+1] = offset_[i] + rows_per_worker[i];
	}

	local_startindex = offset_[procrank_];
	local_endindex = offset_[procrank_] + rows_per_worker[procrank_];	
}

void OpenSMOKE_LIS_Unsymmetric::UpdateVectors(BzzVector &b, BzzVector &x)
{
	if(procrank_ > 0)
	{
	    ChangeDimensions(n, &x);
	    ChangeDimensions(n, &b);
	}

	double *x_ = new double[n + 1];
	double *b_ = new double[n + 1];

	if(procrank_ == 0)
	{
	    for(int i = 1; i <= n; i++)
	    {
		x_[i] = x[i];
		b_[i] = b[i];
	    }
	}

	MPI::COMM_WORLD.Bcast(&x_[1], n, MPI::DOUBLE, MASTER);
	MPI::COMM_WORLD.Bcast(&b_[1], n, MPI::DOUBLE, MASTER);

	if(procrank_ > 0)
	{
	    for(int i = 1; i <= n; i++)
	    {
		x[i] = x_[i];
		b[i] = b_[i];
	    }
	}

	delete [] x_;
	delete [] b_;
}

void OpenSMOKE_LIS_Unsymmetric::Solve(Vec &b, BzzVector &x)
{	

	if (status_	!= OPENSMOKE_DIRECTSOLVER_STATUS_FACTORIZED)
		ErrorMessage("Solve(BzzVector &b, BzzVector &x) can be used only for Factorized Matrices");

	MPI::Status status;

	// Checking dimensions
	PetscInt b_size;
	VecGetSize(b, &b_size);
	if (b_size != x.Size())	ErrorMessage("The size of b and x do not fit!");
	if (b_size != n)			ErrorMessage("The size of matrix and rhs do not fit!");
	
	// From Petsc Vector to LIS Vector
	PetscScalar *RHS_pointer;
	VecGetArray(b, &RHS_pointer);

	double *localb_ = new double[rows_per_worker[procrank_] + 1];

	{
	    for(int i = 1; i <= rows_per_worker[procrank_]; i++)
		localb_[i] = RHS_pointer[i-1];
	}
	
	VecRestoreArray(b, &RHS_pointer);

	for(int i = local_startindex; i < local_endindex; i++)
	    lis_vector_set_value(LIS_INS_VALUE, i, localb_[i + 1 - local_startindex], lis_b);


	//Update Options
	UpdateOptions();

	if(kind_ == OPENSMOKE_DIRECTSOLVER_SQUAREMATRIX)
	    MassFlowRateUpdateOptions();

	// Set options
	char message[300];
	std::string space = " ";
	std::string option_string =	linEqSolv_ + linearEquationSolver_ + space + lisPrecond_ + lisPreconditioner_ + space + maxIter_ + maximumIterations_ + space + tol_ + tolerance_ + messageLevel_ + initialSolution_;

	lis_solver_create(&solver);
	strcpy(message, option_string.c_str());
	lis_solver_set_option(message,solver);

	// Solve linear system
	error_ = lis_solve(A,lis_b,lis_x,solver);

	lis_solver_get_status(solver, &solver_status);

	lis_solver_get_residualnorm(solver, &residual_norm_);

	if(solver_status == 4 && residual_norm_ > 1e-4 && data_->networkStatus() == KPP_NETWORK_STATUS_GLOBALODE)
	{
	    if(procrank_ == 0)
	    {
		std::cout << "Convergence failed with LIS. Residual norm = " << residual_norm_ << std::endl;
		std::cout << "Changing method..." << std::endl;
	    }
	    convergence_index_ = false;
	}

	else if(solver_status != 0 && solver_status != 4 && procrank_ == 0)
	{
	    convergence_index_ = false;
	    std::cout << "There were issues in convergence. Check LIS manual for more" << std::endl;
	    std::cout << "Error code = " << solver_status << std::endl;
	    getchar();
	    exit(-1);
	}

	else convergence_index_ = true;

	double *localx_ = new double[rows_per_worker[procrank_] + 1];

	for(int i=local_startindex;i<local_endindex;i++)
		lis_vector_get_value(lis_x, i, &localx_[i - local_startindex]);

	// From LIS Vector to BzzVector

	for(int p = 0; p <= numworkers_; p++)
	{
	    MPI::COMM_WORLD.Barrier();
	    int mtype = FROM_WORKER;
	    if(procrank_ == p)
	    {
		MPI::COMM_WORLD.Send(&localx_[0], rows_per_worker[p], MPI::DOUBLE, MASTER, mtype);
	    }
	    if(procrank_ == 0)
	    {
		MPI::COMM_WORLD.Recv(&x_[offset_[p]], rows_per_worker[p], MPI::DOUBLE, p, mtype, status);
	    }
	}

	for(int i = 0; i < n; i++)
	{
	    x[i+1] = x_[i];
	}

	lis_solver_destroy(solver);
	delete [] localx_;
	delete [] localb_;
}

void OpenSMOKE_LIS_Unsymmetric::InitializePlaceVectors()
{
	rows_place = new Petsc64bitInt[process_rows];
	columns_place = new Petsc64bitInt[process_values];
	values_place = new Petsc64bitInt[process_values];

	int* rows_vector = new int[nprocs_];
	int* rows_offset = new int[nprocs_];

	for(int i = 0; i < nprocs_; i++)
	{
	    rows_vector[i] = process_rows;
	    MPI::COMM_WORLD.Bcast(&rows_vector[i], 1, MPI::INT, i);
	}

	for(int i = 0; i < nprocs_; i++)
	    rows_offset[i] = 0;

	for(int i = 1; i < nprocs_; i++)
	    rows_offset[i] = rows_offset[i-1] + rows_vector[i-1];


	for(int i = 0; i < process_rows; i++)
	{
	    rows_place[i] = i + rows_offset[procrank_];
	}

	for(int i = 0; i < process_values; i++)
	{
	    columns_place[i] = i + nonzero_offset[procrank_];
	}


	delete [] rows_vector;
	delete [] rows_offset;
}

void OpenSMOKE_LIS_Unsymmetric::Solve(double* &b, BzzVector &x)
{	

	if (status_	!= OPENSMOKE_DIRECTSOLVER_STATUS_FACTORIZED)
		ErrorMessage("Solve(BzzVector &b, BzzVector &x) can be used only for Factorized Matrices");

	MPI::Status status;


	for(int i = local_startindex; i < local_endindex; i++)
	    lis_vector_set_value(LIS_INS_VALUE, i, b[i - local_startindex], lis_b);


	//Update Options
	UpdateOptions();

	if(kind_ == OPENSMOKE_DIRECTSOLVER_SQUAREMATRIX)
	    MassFlowRateUpdateOptions();

	// Set options
	char message[400];
	std::string space = " ";
	std::string option_string =	linEqSolv_ + linearEquationSolver_ + space + lisPrecond_ + lisPreconditioner_ + space + maxIter_ + maximumIterations_ + space + tol_ + tolerance_ + messageLevel_ + initialSolution_;

	lis_solver_create(&solver);
	strcpy(message, option_string.c_str());
	lis_solver_set_option(message,solver);

	// Solve linear system
	if(procrank_ == 0)	std::cout << "   Solving..." << std::endl;
	error_ = lis_solve(A,lis_b,lis_x,solver);

	lis_solver_get_status(solver, &solver_status);

	lis_solver_get_residualnorm(solver, &residual_norm_);

	if(solver_status == 4 && residual_norm_ > 1e-4 && data_->networkStatus() == KPP_NETWORK_STATUS_GLOBALODE)
	{
	    if(procrank_ == 0)
	    {
		std::cout << "Convergence failed with LIS. Residual norm = " << residual_norm_ << std::endl;
		std::cout << "Changing method..." << std::endl;
	    }
	    convergence_index_ = false;
	}

	else if(solver_status != 0 && solver_status != 4 && procrank_ == 0)
	{
	    convergence_index_ = false;
	    std::cout << "There were issues in convergence. Check LIS manual for more" << std::endl;
	    std::cout << "Error code = " << solver_status << std::endl;
	    getchar();
	    exit(-1);
	}

	else convergence_index_ = true;

	double *localx_ = new double[rows_per_worker[procrank_] + 1];


	for(int i=local_startindex;i<local_endindex;i++)
		lis_vector_get_value(lis_x, i, &localx_[i - local_startindex]);


	// From LIS Vector to BzzVector

	for(int p = 1; p <= numworkers_; p++)
	{
	    MPI::COMM_WORLD.Barrier();
	    int mtype = FROM_WORKER;
	    if(procrank_ == p)
	    {
		MPI::COMM_WORLD.Send(&localx_[0], rows_per_worker[p], MPI::DOUBLE, MASTER, mtype);
	    }
	    if(procrank_ == 0)
	    {
		MPI::COMM_WORLD.Recv(&x_[offset_[p]], rows_per_worker[p], MPI::DOUBLE, p, mtype, status);
	    }
	}

	if(procrank_ == 0)
	{
	    for(int i = 0; i < rows_per_worker[0]; i++)
		x_[i] = localx_[i];
	}

	
	if(procrank_ == 0)
	{
	    for(int i = 0; i < n; i++)
	    {
		x[i+1] = x_[i];
	    }
	}


	lis_solver_destroy(solver);
	delete [] localx_;
}


void OpenSMOKE_LIS_Unsymmetric::GetBooleanSparsityPattern(const int numberofrows, const int numberofreactors)
{
    std::ofstream fSparsity("SparsityPattern.out", std::ios::out);
    
    if(nprocs_ > 1)
        ErrorMessage("Boolean Matrix is available only with 1 processor");

    
    sparsity_.resize(numberofrows);
    
    for(int i = 0; i < numberofrows; i++)
    {
        sparsity_[i] = new bool[numberofrows];
    }
    
    for(int i = 0; i < numberofrows; i++)
    {
        for(int j = 0; j < numberofrows; j++)
        {
            sparsity_[i][j] = false;
        }
    }
    
    
    
    int index = 0;
    int global_index = 0;
    int *global_row = new int[numberNonZeroElements];
    
    for(int i = 0; i < numberofrows; i++)
    {
        int ElementsOnTheRow = loc_rows[index+1] - loc_rows[index];
        for(int j = 0; j < ElementsOnTheRow; j++)
        {
            if(loc_values[global_index] != 0)
            {
                fSparsity << index + 1 << std::setw(15) << std::right << loc_columns[global_index] + 1 << std::setw(15) << "1.0" << std::endl;
/*                sparsity_[index][loc_columns[global_index]] = true;
                global_row[global_index] = index;*/
            }
            global_index++;
        }
        
        index++;
    }
     
/*    index = 0;
    global_index = 0;
    
        for(int i = 0; i < numberofrows; i++)
    {
        int ElementsOnTheRow = loc_rows[index+1] - loc_rows[index];
        for(int j = 0; j < ElementsOnTheRow; j++)
        {
            if(loc_values[global_index] != 0)
            {
                global_row[global_index] = index;
            }
            global_index++;
        }
        
        index++;
    }*/

    
/*    for(int i = 0; i < numberofrows; i++)                //Sostituire con numberofrows
    {
        for(int j = 0; j < numberofrows; j++)
        {
            fSparsity << sparsity_[i][j] << " ";
        }
        fSparsity << std::endl;
    }*/
    /*
    for(int i = 0; i < global_index; i++)
    {
        fSparsity << global_row[i] << setw(15) << loc_columns[i] << std::endl;
    }*/
        
        
/*
        // Non zero elements
    {
        int nLocalRows = indicesColumns.Size() / nElementsPerRows[index];
        for (int k = 1; k <= nLocalRows; k++)
        {
            loc_rows[countGlobalRows_++] = countGlobal_ + 1 + nElementsPerRows[index]*(k - 1) + nonzero_offset[procrank_];
        }

        for (int k = 1; k <= indicesColumns.Size(); k++)
        {
            loc_columns[countGlobal_] = indicesColumns[k];
            countGlobal_++;
        }
    }*/
}