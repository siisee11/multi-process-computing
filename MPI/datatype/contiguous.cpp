#include <stdio.h>
#include<mpi.h>

int main(int argc, char **argv){
	int i, myrank, ibuf[20];
	MPI_Datatype inewtype;
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
	
	if(myrank==0) for(i=0; i<20; i++) ibuf[i]=i+1;
	else for(i=0; i<20; i++) ibuf[i]=0;

	MPI_Type_contiguous(3, MPI_INT, &inewtype);
	MPI_Type_commit(&inewtype);
	MPI_Bcast(ibuf, 3, inewtype, 0, MPI_COMM_WORLD);
//	if(myrank==0){
		printf("%d: ibuf = ", myrank);
		for(i = 0; i < 20 ; i++) printf(" %d", ibuf[i]);
		printf("\n");
//	}
	MPI_Finalize();
}
