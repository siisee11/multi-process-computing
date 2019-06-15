#include<stdio.h>
#include<mpi.h>

int main(int argc, char **argv){

	int nprocs, myrank;
	
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
	MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
	
	while(1){

	}
	
	MPI_Finalize();
}
