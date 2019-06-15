#include<stdio.h>
#include<mpi.h>

int main(int argc, char **argv){
	int i, nprocs, myrank;
	int index[6], edges[12], reorder;
	MPI_Comm oldcomm, newcomm;
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
	MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
	reorder =0;
	if(nprocs == 6){
		oldcomm =MPI_COMM_WORLD;
		index[0]=2;
		for(i=1; i<nprocs; ++i) index[i] = index[i-1]+2;
		edges[0]=5; edges[1]=2;
		edges[2]=2; edges[3]=3;
		edges[4]=0; edges[5]=1;
		edges[6]=1; edges[7]=4;
		edges[8]=3; edges[9]=5;
		edges[10]=4; edges[11]=0;
		MPI_Graph_create(oldcomm, nprocs, index, edges, reorder, &newcomm);

	}else{
		printf("NoP mustbe 6\n");
	}
	MPI_Finalize();
}

