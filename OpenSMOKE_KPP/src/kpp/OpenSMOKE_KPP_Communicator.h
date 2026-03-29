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

#ifndef OpenSMOKE_KPP_Communicator_H
#define OpenSMOKE_KPP_Communicator_H

#include "BzzMath.hpp"
#include <mpi.h>
#include <vector>
#include <petscvec.h>

class OpenSMOKE_KPP_ReactorNetwork;

class OpenSMOKE_KPP_Communicator
{
public:

	//Constructor and destructor
	OpenSMOKE_KPP_Communicator(int nprocesses, int processrank, OpenSMOKE_KPP_ReactorNetwork& network);
	~OpenSMOKE_KPP_Communicator(void); 

	inline int nprocs()				const { return nprocs_; }
	inline int procrank()				const { return procrank_; }

	//Functions
	void InitializeParameters();

	void InitializeArray(int*& local_array, int*& global_array);
	void InitializeArray(long long int*& local_array, long long int*& global_array);
	void InitializeArray(double*& local_array, double*& global_array);
	void InitializeArray(int*& array);
	void InitializeArray(double*& array);
	void InitializeArray(double*& array, int length);
	void InitializeArray(int*& array, int length);
	void InitializeArray(long long int*& array, int length);
	void InitializeVector(std::vector<int*>& local_vector, std::vector<int*>& global_vector, int* length);
	void InitializeVector(std::vector<double*>& local_vector, std::vector<double*>& global_vector, int* length);
	void InitializeVector(std::vector<int*>& vec, int* length);
	void InitializeVector(std::vector<double*>& vec, int* length);
	void InitializeVector(std::vector<int*>& vec, int length);
	void InitializeVector(std::vector<double*>& vec, int length);
	void InitializeVector(std::vector<int*>& local_vector, std::vector<int*>& global_vector, int length);
	void InitializeVector(std::vector<double*>& local_vector, std::vector<double*>& global_vector, int length);
	void InitializeVector(std::vector<int*>& vec, int size, int* length);
	void InitializeVector(std::vector<double*>& vec, int size, int* length);
	void InitializePetscVector(Vec &vect, PetscInt &loc_size, PetscInt &upper_index, PetscInt &lower_index);
	void InitializePetscVector(Vec &vect, PetscInt &loc_size, PetscInt &upper_index, PetscInt &lower_index, Petsc64bitInt*& place);
	void BroadcastArray(int*& array, int length);
	void BroadcastArray(double*& array, int length);
	void BroadcastVector(std::vector<int*>& vec, int* length);
	void BroadcastVector(std::vector<double*>& vec, int* length);
	void BroadcastVector(std::vector<int*>& vec, int length);
	void BroadcastVector(std::vector<double*>& vec, int length);
	void BroadcastVector(std::vector<int*>& vec, int length, bool* communication);
	void BroadcastVector(std::vector<double*>& vec, int length, bool* communication);
	void BroadcastVector(std::vector<int*>& vec, int* length, bool* communication);
	void BroadcastVector(std::vector<double*>& vec, int* length, bool* communication);
	void SendArrayToMaster(int*& local_array, int*& global_array, int* length);
	void SendArrayToMaster(double*& local_array, double*& global_array, int* length);
	void SendVectorToMaster(std::vector<int*>& local_vector, std::vector<int*>& global_vector, int* length);
	void SendVectorToMaster(std::vector<double*>& local_vector, std::vector<double*>& global_vector, int* length);
	void SendVectorToMaster(std::vector<int*>& local_vector, std::vector<int*>& global_vector, int length);
	void SendVectorToMaster(std::vector<double*>& local_vector, std::vector<double*>& global_vector, int length);
	void DeleteArray(int*& array);
	void DeleteArray(double*& array);
	void DeleteVector(std::vector<int*>& local_vector, std::vector<int*>& global_vector);
	void DeleteVector(std::vector<double*>& local_vector, std::vector<double*>& global_vector);
	void DeleteVector(std::vector<double*>& vec, int size);
	void DeleteVector(std::vector<int*>& vec, int size);
	void GatherArray(int*& array, int* size, int* offset);
	void GatherArray(double*& array, int* size, int* offs);
	void GatherVector(std::vector<int*>& vec, int* length);
	void GatherVector(std::vector<double*>& vec, int* length);
	void GatherVector(std::vector<int*>& vec, int length);
	void GatherVector(std::vector<double*>& vec, int length);
	void MakeArrayGlobal(int* array);
	void MakeArrayGlobal(double* array);
	void MakeArrayLocal(double* array);

private:

	int nprocs_, procrank_, numworkers_;
	int FROM_MASTER, FROM_WORKER;
	int MASTER;
	int source, mtype;
	int NR_, Nspecies;
	int* NR_P;
	int* offset;
	MPI::Status status;

	//Communication utilities


	//Functions

	OpenSMOKE_KPP_ReactorNetwork& network_;

};


#endif //OpenSMOKE_KPP_Communicator_H
