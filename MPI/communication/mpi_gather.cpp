#include<stdio.h>
#include<mpi.h>

int main(int argc, char **argv){

	int i, nprocs, myrank;
	int isend, irecv[3];
	
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
	MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
	isend = myrank + 1;
	
	MPI_Gather(&isend, 1, MPI_INT, irecv, 1, MPI_INT, 0, MPI_COMM_WORLD);

	if(myrank==0){
		printf("irecv = ");
		for(i=0; i<3; i++) printf("%d, ", irecv[i]);
		printf("\n");
	}
	MPI_Finalize();
}
