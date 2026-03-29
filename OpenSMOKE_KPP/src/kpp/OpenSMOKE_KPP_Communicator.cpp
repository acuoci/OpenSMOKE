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

#include "OpenSMOKE_KPP_Communicator.h"
#include "OpenSMOKE_KPP_ReactorNetwork.h"
#include <mpi.h>
#include <vector>
#include <petscvec.h>
#include "BzzMath.hpp"

OpenSMOKE_KPP_Communicator::OpenSMOKE_KPP_Communicator(int nprocesses, int processrank, OpenSMOKE_KPP_ReactorNetwork& network) : network_(network)
{
	nprocs_ = nprocesses;
	numworkers_ = nprocs_ - 1;
	procrank_ = processrank;
	MASTER = 0;
	FROM_MASTER = 1;
	FROM_WORKER = 2;

	NR_P = new int[nprocs_];
	offset = new int[nprocs_];
}

OpenSMOKE_KPP_Communicator::~OpenSMOKE_KPP_Communicator()
{
}

void OpenSMOKE_KPP_Communicator::InitializeParameters()
{	
	for(int i = 0; i <= numworkers_; i++)
	{
	    NR_P[i] = network_.LocalNumberOfReactors()[i];
	    offset[i] = network_.offs()[i];
	}

	NR_ = network_.NumberOfReactors();
	Nspecies = network_.NumberOfSpecies();
}

void OpenSMOKE_KPP_Communicator::SendVectorToMaster(std::vector<int*>& local_vector, std::vector<int*>& global_vector, int length)
{
        mtype = FROM_WORKER;
	for(int p = 1; p <= numworkers_; p++)
        {
	    MPI::COMM_WORLD.Barrier();
	    if(p == procrank_)
	    {
		for(int k = 1; k <= NR_P[procrank_]; k++)
		{
		    int index = p + numworkers_ * (k-1);
		    MPI::COMM_WORLD.Send(&local_vector[k][1], length, MPI::INT, MASTER, mtype);
		}
	    }
	    if(procrank_ == 0)
	    {
		source = p;
		for(int k = 1; k <= NR_P[source]; k++)
		{
		    int index = source + numworkers_ * (k-1);
		    MPI::COMM_WORLD.Recv(&local_vector[k + offset[source]][1], length, MPI::INT, source, mtype, status);
		}
	    }
	}

	if(procrank_ == 0)
	{
	    for(int i = 1; i <= numworkers_; i++)
	    {
		for(int j = 1; j <= NR_P[i]; j++)
		{
		    int global_index = i + numworkers_ * (j-1);
		    for(int k = 1; k <= length; k++)
			global_vector[global_index][k] = local_vector[j + offset[i]][k];
		}
	    }
	}
}

void OpenSMOKE_KPP_Communicator::SendVectorToMaster(std::vector<double*>& local_vector, std::vector<double*>& global_vector, int length)
{
        mtype = FROM_WORKER;
	for(int p = 1; p <= numworkers_; p++)
        {
	    MPI::COMM_WORLD.Barrier();
	    if(p == procrank_)
	    {
		for(int k = 1; k <= NR_P[procrank_]; k++)
		{
		    int index = p + numworkers_ * (k-1);
		    MPI::COMM_WORLD.Send(&local_vector[k][1], length, MPI::DOUBLE, MASTER, mtype);
		}
	    }
	    if(procrank_ == 0)
	    {
		source = p;
		for(int k = 1; k <= NR_P[source]; k++)
		{
		    int index = source + numworkers_ * (k-1);
		    MPI::COMM_WORLD.Recv(&local_vector[k + offset[source]][1], length, MPI::DOUBLE, source, mtype, status);
		}
	    }
	}

	if(procrank_ == 0)
	{
	    for(int i = 1; i <= numworkers_; i++)
	    {
		for(int j = 1; j <= NR_P[i]; j++)
		{
		    int global_index = i + numworkers_ * (j-1);
		    for(int k = 1; k <= length; k++)
			global_vector[global_index][k] = local_vector[j + offset[i]][k];
		}
	    }
	}
}


void OpenSMOKE_KPP_Communicator::BroadcastArray(int*& array, int length)
{
	MPI::COMM_WORLD.Bcast(&array[0], length + 1, MPI::INT, MASTER);
}


void OpenSMOKE_KPP_Communicator::BroadcastArray(double*& array, int length)
{
	MPI::COMM_WORLD.Bcast(&array[0], length + 1, MPI::DOUBLE, MASTER);
}

void OpenSMOKE_KPP_Communicator::SendVectorToMaster(std::vector<int*>& local_vector, std::vector<int*>& global_vector, int* length)
{
        mtype = FROM_WORKER;
	for(int p = 1; p <= numworkers_; p++)
        {
	    MPI::COMM_WORLD.Barrier();
	    if(p == procrank_)
	    {
		for(int k = 1; k <= NR_P[procrank_]; k++)
		{
		    int index = p + numworkers_ * (k-1);
		    MPI::COMM_WORLD.Send(&local_vector[k][1], length[index], MPI::INT, MASTER, mtype);
		}
	    }
	    if(procrank_ == 0)
	    {
		source = p;
		for(int k = 1; k <= NR_P[source]; k++)
		{
		    int index = source + numworkers_ * (k-1);
		    MPI::COMM_WORLD.Recv(&local_vector[k + offset[source]][1], length[index], MPI::INT, source, mtype, status);
		}
	    }
	}

	if(procrank_ == 0)
	{
	    for(int i = 1; i <= numworkers_; i++)
	    {
		for(int j = 1; j <= NR_P[i]; j++)
		{
		    int global_index = i + numworkers_ * (j-1);
		    for(int k = 1; k <= length[global_index]; k++)
			global_vector[global_index][k] = local_vector[j + offset[i]][k];
		}
	    }
	}
}

void OpenSMOKE_KPP_Communicator::SendVectorToMaster(std::vector<double*>& local_vector, std::vector<double*>& global_vector, int* length)
{

        mtype = FROM_WORKER;

	for(int p = 1; p <= numworkers_; p++)
        {
	    MPI::COMM_WORLD.Barrier();
	    if(p == procrank_)
	    {
		for(int k = 1; k <= NR_P[procrank_]; k++)
		{
		    int index = p + numworkers_ * (k-1);
		    MPI::COMM_WORLD.Send(&local_vector[k][1], length[index], MPI::DOUBLE, MASTER, mtype);
		}
	    }
	    if(procrank_ == 0)
	    {
		source = p;
		for(int k = 1; k <= NR_P[source]; k++)
		{
		    int index = source + numworkers_ * (k-1);
		    MPI::COMM_WORLD.Recv(&local_vector[k + offset[source]][1], length[index], MPI::DOUBLE, source, mtype, status);
		}
	    }
	}


	if(procrank_ == 0)
	{
	    for(int i = 1; i <= numworkers_; i++)
	    {
		for(int j = 1; j <= NR_P[i]; j++)
		{
		    int global_index = i + numworkers_ * (j-1);
		    for(int k = 1; k <= length[global_index]; k++)
		    {
			global_vector[global_index][k] = local_vector[j + offset[i]][k];
		    }
		}
	    }
	}
}

void OpenSMOKE_KPP_Communicator::SendArrayToMaster(int*& local_array, int*& global_array, int* length)
{
	mtype = FROM_WORKER;
	for(int p = 1; p <= numworkers_; p++)
        {
	    MPI::COMM_WORLD.Barrier();
	    if(p == procrank_)
	    {
		MPI::COMM_WORLD.Send(&local_array[1], length[p], MPI::INT, MASTER, mtype);
	    }
	    if(procrank_ == 0)
	    {
		source = p;
		MPI::COMM_WORLD.Recv(&local_array[1 + offset[p]], length[p], MPI::INT, source, mtype, status);
	    }
	}

	if(procrank_ == 0)
	{
	    for(int i = 1; i <= numworkers_; i++)
            {
                for(int j = 1; j <= NR_P[i]; j++)
                {
                    global_array[i + numworkers_*(j-1)] = local_array[j + offset[i]];
                }
            }
	}	
}

void OpenSMOKE_KPP_Communicator::InitializeArray(int*& local_array, int*& global_array)
{
	if(procrank_ == 0)
	{
	    local_array = new int[NR_ + 1];
	    global_array = new int[NR_ + 1];
	    for(int i = 0; i <= NR_; i++)
	    {
		local_array[i] = 0;
		global_array[i] = 0;
	    }
	}

	if(procrank_ > 0)
	{
	    local_array = new int[NR_P[procrank_] + 1];
	    global_array = new int[NR_ + 1];

	    for(int i = 0; i <= NR_; i++)
	    {
		global_array[i] = 0;
	    }

	    for(int i = 0; i <= NR_P[procrank_]; i++)
		local_array[i] = 0;
	}
}

void OpenSMOKE_KPP_Communicator::InitializeArray(long long int*& local_array, long long int*& global_array)
{
	if(procrank_ == 0)
	{
	    local_array = new long long int[NR_ + 1];
	    global_array = new long long int[NR_ + 1];
	    for(int i = 0; i <= NR_; i++)
	    {
		local_array[i] = 0;
		global_array[i] = 0;
	    }
	}

	if(procrank_ > 0)
	{
	    local_array = new long long int[NR_P[procrank_] + 1];
	    global_array = new long long int[NR_ + 1];

	    for(int i = 0; i <= NR_; i++)
	    {
		global_array[i] = 0;
	    }

	    for(int i = 0; i <= NR_P[procrank_]; i++)
		local_array[i] = 0;
	}
}

void OpenSMOKE_KPP_Communicator::SendArrayToMaster(double*& local_array, double*& global_array, int* length)
{
	mtype = FROM_WORKER;
	for(int p = 1; p <= numworkers_; p++)
        {
	    MPI::COMM_WORLD.Barrier();
	    if(p == procrank_)
	    {
		MPI::COMM_WORLD.Send(&local_array[1], length[p], MPI::DOUBLE, MASTER, mtype);
	    }
	    if(procrank_ == 0)
	    {
		source = p;
		MPI::COMM_WORLD.Recv(&local_array[1 + offset[p]], length[p], MPI::DOUBLE, source, mtype, status);
	    }
	}

	if(procrank_ == 0)
	{
	    for(int i = 1; i <= numworkers_; i++)
            {
                for(int j = 1; j <= NR_P[i]; j++)
                {
                    global_array[i + numworkers_*(j-1)] = local_array[j + offset[i]];
                }
            }
	}	
}

void OpenSMOKE_KPP_Communicator::InitializeArray(double*& local_array, double*& global_array)
{
	if(procrank_ == 0)
	{
	    local_array = new double[NR_ + 1];
	    global_array = new double[NR_ + 1];
	    for(int i = 0; i <= NR_; i++)
	    {
		local_array[i] = 0;
		global_array[i] = 0;
	    }
	}

	if(procrank_ > 0)
	{
	    local_array = new double[NR_P[procrank_] + 1];
	    global_array = new double[NR_ + 1];

	    for(int i = 0; i <= NR_; i++)
	    {
		global_array[i] = 0;
	    }

	    for(int i = 0; i <= NR_P[procrank_]; i++)
		local_array[i] = 0;
	}
}

void OpenSMOKE_KPP_Communicator::InitializeVector(std::vector<int*>& local_vector, std::vector<int*>& global_vector, int* length)
{
	if(procrank_ == 0)
	{
	    local_vector.resize(NR_ + 1);
	    for(int i = 1; i <= numworkers_; i++)
            {
                for(int j = 1; j <= NR_P[i]; j++)
                {
                    local_vector[j + offset[i]] = new int[length[i + numworkers_*(j-1)] + 1];
                }
            }
		
	    global_vector.resize(NR_ + 1);
	    for(int k = 1; k <= NR_; k++)
		global_vector[k] = new int[length[k] + 1];
	}

	if(procrank_ > 0)
	{
	    local_vector.resize(NR_P[procrank_] + 1);
	    for(int k = 1; k <= NR_P[procrank_]; k++)
		local_vector[k] = new int[length[procrank_ + numworkers_ * (k-1)] + 1];

	    global_vector.resize(NR_ + 1);
	    for(int k = 1; k <= NR_; k++)
		global_vector[k] = new int[length[k] + 1];
	}
}

void OpenSMOKE_KPP_Communicator::InitializeVector(std::vector<double*>& local_vector, std::vector<double*>& global_vector, int* length)
{
	if(procrank_ == 0)
	{
	    local_vector.resize(NR_ + 1);
	    for(int i = 1; i <= numworkers_; i++)
            {
                for(int j = 1; j <= NR_P[i]; j++)
                {
                    local_vector[j + offset[i]] = new double[length[i + numworkers_*(j-1)] + 1];
                }
            }

	    global_vector.resize(NR_ + 1);
	    for(int k = 1; k <= NR_; k++)
		global_vector[k] = new double[length[k] + 1];
	}

	if(procrank_ > 0)
	{
	    local_vector.resize(NR_P[procrank_] + 1);
	    for(int k = 1; k <= NR_P[procrank_]; k++)
		local_vector[k] = new double[length[procrank_ + numworkers_ * (k-1)] + 1];

	    global_vector.resize(NR_ + 1);
	    for(int k = 1; k <= NR_; k++)
		global_vector[k] = new double[length[k] + 1];
	}
}

void OpenSMOKE_KPP_Communicator::BroadcastVector(std::vector<int*>& vec, int* length)
{
	for(int k = 1; k <= NR_; k++)
        {
	    MPI::COMM_WORLD.Bcast(&vec[k][0], length[k] + 1, MPI::INT, MASTER);
	}
}

void OpenSMOKE_KPP_Communicator::BroadcastVector(std::vector<double*>& vec, int length)
{
	for(int k = 1; k <= NR_; k++)
        {
	    MPI::COMM_WORLD.Bcast(&vec[k][0], length + 1, MPI::DOUBLE, MASTER);
	}
}

void OpenSMOKE_KPP_Communicator::BroadcastVector(std::vector<int*>& vec, int length)
{
	for(int k = 1; k <= NR_; k++)
        {
	    MPI::COMM_WORLD.Bcast(&vec[k][0], length + 1, MPI::INT, MASTER);
	}
}

void OpenSMOKE_KPP_Communicator::BroadcastVector(std::vector<double*>& vec, int* length)
{
	for(int k = 1; k <= NR_; k++)
        {
	    MPI::COMM_WORLD.Bcast(&vec[k][0], length[k] + 1, MPI::DOUBLE, MASTER);
	}
}

void OpenSMOKE_KPP_Communicator::InitializeArray(int*& array)
{
	array = new int[NR_ + 1];

	for(int i = 0; i <= NR_; i++)
	    array[i] = 0;
}

void OpenSMOKE_KPP_Communicator::InitializeArray(int*& array, int length)
{
	array = new int[length + 1];

	for(int i = 0; i <= length; i++)
	    array[i] = 0;
}

void OpenSMOKE_KPP_Communicator::InitializeArray(long long int*& array, int length)
{
	array = new long long int[length + 1];

	for(int i = 0; i <= length; i++)
	    array[i] = 0;
}

void OpenSMOKE_KPP_Communicator::InitializeArray(double*& array)
{
	array = new double[NR_ + 1];

	for(int i = 0; i <= NR_; i++)
	    array[i] = 0;
}

void OpenSMOKE_KPP_Communicator::InitializeArray(double*& array, int length)
{
	array = new double[length + 1];

	for(int i = 0; i <= length; i++)
	    array[i] = 0;
}

void OpenSMOKE_KPP_Communicator::InitializeVector(std::vector<int*>& vec, int* length)
{
	vec.resize(NR_ + 1);
	for(int k = 1; k <= NR_; k++)
	    vec[k] = new int[length[k] + 1];
}

void OpenSMOKE_KPP_Communicator::InitializeVector(std::vector<double*>& vec, int* length)
{
	vec.resize(NR_ + 1);
	for(int k = 1; k <= NR_; k++)
	    vec[k] = new double[length[k] + 1];
}

void OpenSMOKE_KPP_Communicator::InitializeVector(std::vector<double*>& vec, int length)
{
	vec.resize(NR_ + 1);
	for(int k = 1; k <= NR_; k++)
	    vec[k] = new double[length + 1];
}

void OpenSMOKE_KPP_Communicator::InitializeVector(std::vector<int*>& vec, int length)
{
	vec.resize(NR_ + 1);
	for(int k = 1; k <= NR_; k++)
	    vec[k] = new int[length + 1];
}

void OpenSMOKE_KPP_Communicator::InitializeVector(std::vector<double*>& local_vector, std::vector<double*>& global_vector, int length)
{
	if(procrank_ == 0)
	{
	    local_vector.resize(NR_ + 1);
	    for(int k = 1; k <= NR_; k++)
		local_vector[k] = new double[length + 1];

	    global_vector.resize(NR_ + 1);
	    for(int k = 1; k <= NR_; k++)
		global_vector[k] = new double[length + 1];
	}

	if(procrank_ > 0)
	{
	    local_vector.resize(NR_P[procrank_] + 1);
	    for(int k = 1; k <= NR_P[procrank_]; k++)
		local_vector[k] = new double[length + 1];

	    global_vector.resize(NR_ + 1);
	    for(int k = 1; k <= NR_; k++)
		global_vector[k] = new double[length + 1];
	}
}


void OpenSMOKE_KPP_Communicator::InitializeVector(std::vector<int*>& local_vector, std::vector<int*>& global_vector, int length)
{
	if(procrank_ == 0)
	{
	    local_vector.resize(NR_ + 1);
	    for(int k = 1; k <= NR_; k++)
		local_vector[k] = new int[length + 1];

	    global_vector.resize(NR_ + 1);
	    for(int k = 1; k <= NR_; k++)
		global_vector[k] = new int[length + 1];
	}

	if(procrank_ > 0)
	{
	    local_vector.resize(NR_P[procrank_] + 1);
	    for(int k = 1; k <= NR_P[procrank_]; k++)
		local_vector[k] = new int[length + 1];

	    global_vector.resize(NR_ + 1);
	    for(int k = 1; k <= NR_; k++)
		global_vector[k] = new int[length + 1];
	}
}


void OpenSMOKE_KPP_Communicator::DeleteArray(int*& array)
{
	delete [] array;
}

void OpenSMOKE_KPP_Communicator::DeleteArray(double*& array)
{
	delete [] array;
}

void OpenSMOKE_KPP_Communicator::DeleteVector(std::vector<int*>& local_vector, std::vector<int*>& global_vector)
{
	if(procrank_ == 0)
	{
	    for(int k = 1; k <= NR_; k++)
	    {
		delete [] local_vector[k];
		delete [] global_vector[k];
	    }
	}

	if(procrank_ > 0)
	{
	    for(int k = 1; k <= NR_P[procrank_]; k++)
		delete [] local_vector[k];

	    for(int k = 1; k <= NR_; k++)
		delete [] global_vector[k];
	}
}

void OpenSMOKE_KPP_Communicator::DeleteVector(std::vector<double*>& local_vector, std::vector<double*>& global_vector)
{
	if(procrank_ == 0)
	{
	    for(int k = 1; k <= NR_; k++)
	    {
		delete [] local_vector[k];
		delete [] global_vector[k];
	    }
	}

	if(procrank_ > 0)
	{
	    for(int k = 1; k <= NR_P[procrank_]; k++)
		delete [] local_vector[k];

	    for(int k = 1; k <= NR_; k++)
		delete [] global_vector[k];
	}
}

void OpenSMOKE_KPP_Communicator::DeleteVector(std::vector<double*>& vec, int size)
{
	for(int k = 1; k <= size; k++)
	{
	    delete [] vec[k];
	}
}

void OpenSMOKE_KPP_Communicator::DeleteVector(std::vector<int*>& vec, int size)
{
	for(int k = 1; k <= size; k++)
	{
	    delete [] vec[k];
	}
}

void OpenSMOKE_KPP_Communicator::InitializeVector(std::vector<int*>& vec, int size, int* length)
{
	vec.resize(size + 1);
	
	for(int k = 1; k <= size; k++)
    	    vec[k] = new int[length[k]];
}

void OpenSMOKE_KPP_Communicator::InitializeVector(std::vector<double*>& vec, int size, int* length)
{
	vec.resize(size + 1);
	
	for(int k = 1; k <= size; k++)
    	    vec[k] = new double[length[k]];
}

void OpenSMOKE_KPP_Communicator::BroadcastVector(std::vector<int*>& vec, int* length, bool* communication)
{
	if(procrank_ == 0)
	{	
	    for(int k = 1; k <= NR_; k++)
            {
	        MPI::COMM_WORLD.Bcast(&vec[k][0], length[k] + 1, MPI::INT, MASTER);
	    }
	}

	if(procrank_ > 0)
	{
	    for(int k = 1; k <= NR_; k++)
            {
		if(communication[k] == true)
	            MPI::COMM_WORLD.Bcast(&vec[k][0], length[k] + 1, MPI::INT, MASTER);

		else
		{
		    int* dummy = new int[length[k] + 1];
		    MPI::COMM_WORLD.Bcast(&dummy[0], length[k] + 1, MPI::INT, MASTER);
		    delete [] dummy;
		}
	    }
	}
}

void OpenSMOKE_KPP_Communicator::BroadcastVector(std::vector<double*>& vec, int* length, bool* communication)
{
	if(procrank_ == 0)
	{	
	    for(int k = 1; k <= NR_; k++)
            {
	        MPI::COMM_WORLD.Bcast(&vec[k][0], length[k] + 1, MPI::DOUBLE, MASTER);
	    }
	}

	if(procrank_ > 0)
	{
	    for(int k = 1; k <= NR_; k++)
            {
		if(communication[k] == true)
	            MPI::COMM_WORLD.Bcast(&vec[k][0], length[k] + 1, MPI::DOUBLE, MASTER);

		else
		{
		    double *dummy = new double[length[k] + 1];
		    MPI::COMM_WORLD.Bcast(&dummy[0], length[k] + 1, MPI::DOUBLE, MASTER);
		    delete [] dummy;
		}
	    }
	}
}

void OpenSMOKE_KPP_Communicator::BroadcastVector(std::vector<int*>& vec, int length, bool* communication)
{
	if(procrank_ == 0)
	{	
	    for(int k = 1; k <= NR_; k++)
            {
	        MPI::COMM_WORLD.Bcast(&vec[k][0], length + 1, MPI::INT, MASTER);
	    }
	}

	if(procrank_ > 0)
	{
	    for(int k = 1; k <= NR_; k++)
            {
		if(communication[k] == true)
	            MPI::COMM_WORLD.Bcast(&vec[k][0], length + 1, MPI::INT, MASTER);

		else
		{
		    int* dummy = new int[length + 1];
		    MPI::COMM_WORLD.Bcast(&dummy[0], length + 1, MPI::INT, MASTER);
		    delete [] dummy;
		}
	    }
	}
}

void OpenSMOKE_KPP_Communicator::BroadcastVector(std::vector<double*>& vec, int length, bool* communication)
{
	if(procrank_ == 0)
	{	
	    for(int k = 1; k <= NR_; k++)
            {
	        MPI::COMM_WORLD.Bcast(&vec[k][0], length + 1, MPI::DOUBLE, MASTER);
	    }
	}

	if(procrank_ > 0)
	{
	    for(int k = 1; k <= NR_; k++)
            {
		if(communication[k] == true)
	            MPI::COMM_WORLD.Bcast(&vec[k][0], length + 1, MPI::DOUBLE, MASTER);

		else
		{
		    double *dummy = new double[length + 1];
		    MPI::COMM_WORLD.Bcast(&dummy[0], length + 1, MPI::DOUBLE, MASTER);
		    delete [] dummy;
		}
	    }
	}
}

void OpenSMOKE_KPP_Communicator::InitializePetscVector(Vec &vect, PetscInt &loc_size, PetscInt &upper_index, PetscInt &lower_index)
{
	VecCreate(PETSC_COMM_WORLD, &vect);
	VecSetSizes(vect, loc_size, PETSC_DECIDE);
	VecSetFromOptions(vect);
	VecGetOwnershipRange(vect, &lower_index, &upper_index);
}

void OpenSMOKE_KPP_Communicator::InitializePetscVector(Vec &vect, PetscInt &loc_size, PetscInt &upper_index, PetscInt &lower_index, Petsc64bitInt*& place)
{
	VecCreate(PETSC_COMM_WORLD, &vect);
	VecSetSizes(vect, loc_size, PETSC_DECIDE);
	VecSetFromOptions(vect);
	VecGetOwnershipRange(vect, &lower_index, &upper_index);

	place = new Petsc64bitInt[loc_size];
	for(int i = 0; i < loc_size; i++)
    	{
	    place[i] = i + lower_index;
    	}
}

void OpenSMOKE_KPP_Communicator::GatherArray(int*& array, int* size, int* offs)
{

	MPI::COMM_WORLD.Allgatherv(&array[1 + offs[procrank_]], NR_P[procrank_], MPI::INT, &array[1], size, offs, MPI::INT);

}

void OpenSMOKE_KPP_Communicator::GatherArray(double*& array, int* size, int* offs)
{

	MPI::COMM_WORLD.Allgatherv(&array[1 + offs[procrank_]], NR_P[procrank_], MPI::DOUBLE, &array[1], size, offs, MPI::DOUBLE);

}

void OpenSMOKE_KPP_Communicator::GatherVector(std::vector<int*>& vec, int* length)
{

	for(int p = 0; p < nprocs_; p++)
	{
	    for(int k = 1; k <= NR_P[p]; k++)
	    {
		int index = k + offset[p];
	 	MPI::COMM_WORLD.Bcast(&vec[index][1], length[index], MPI::INT, p);
	    }
	}
}

void OpenSMOKE_KPP_Communicator::GatherVector(std::vector<double*>& vec, int* length)
{

	for(int p = 0; p < nprocs_; p++)
	{
	    for(int k = 1; k <= NR_P[p]; k++)
	    {
		int index = k + offset[p];
	 	MPI::COMM_WORLD.Bcast(&vec[index][1], length[index], MPI::DOUBLE, p);
	    }
	}
}

void OpenSMOKE_KPP_Communicator::GatherVector(std::vector<int*>& vec, int length)
{

	for(int p = 0; p < nprocs_; p++)
	{
	    for(int k = 1; k <= NR_P[p]; k++)
	    {
		int index = k + offset[p];
	 	MPI::COMM_WORLD.Bcast(&vec[index][1], length, MPI::INT, p);
	    }
	}
}

void OpenSMOKE_KPP_Communicator::GatherVector(std::vector<double*>& vec, int length)
{

	for(int p = 0; p < nprocs_; p++)
	{
	    for(int k = 1; k <= NR_P[p]; k++)
	    {
		int index = k + offset[p];
	 	MPI::COMM_WORLD.Bcast(&vec[index][1], length, MPI::DOUBLE, p);
	    }
	}
}

void OpenSMOKE_KPP_Communicator::MakeArrayGlobal(int* array)
{
	int *local_array = new int[NR_ + 1];

	for(int i = 0; i <= numworkers_; i++)
	{
	    for(int j = 1; j <= NR_P[i]; j++)
	    {
		int index = i + 1 + nprocs_ * (j-1);
		local_array[index] = array[j + offset[i]];
	    }
	}

	for(int i = 1; i <= NR_; i++)
	    array[i] = local_array[i];

	delete [] local_array;
}

void OpenSMOKE_KPP_Communicator::MakeArrayGlobal(double* array)
{
	double *local_array = new double[NR_ + 1];

	for(int i = 0; i <= numworkers_; i++)
	{
	    for(int j = 1; j <= NR_P[i]; j++)
	    {
		int index = i + 1 + nprocs_ * (j-1);
		local_array[index] = array[j + offset[i]];
	    }
	}

	for(int i = 1; i <= NR_; i++)
	    array[i] = local_array[i];

	delete [] local_array;
}

void OpenSMOKE_KPP_Communicator::MakeArrayLocal(double* array)
{
	double *local_array = new double[NR_ + 1];

	for(int i = 0; i <= numworkers_; i++)
	{
	    for(int j = 1; j <= NR_P[i]; j++)
	    {
		int index = i + 1 + nprocs_ * (j-1);
		local_array[j + offset[i]] = array[index];
	    }
	}

	for(int i = 1; i <= NR_; i++)
	    array[i] = local_array[i];

	delete [] local_array;
}
