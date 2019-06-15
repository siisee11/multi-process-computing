#include<mpi.h>
#include<stdio.h>

int main (int argc, char **argv){
	int nprocs, myrank;
	int coords[2], rank;
	int ndims, newprocs, newrank;
	int i,j;
	MPI_Comm newcomm;
	int dimsize[2], periods[2], reorder;
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
	MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

	ndims =2; 
	dimsize[0]=3; dimsize[1]=2;
	//MPI_Dims_create(nprocs, ndims, dimsize);
	periods[0] = 1; periods[1] = 0; reorder =1;
	MPI_Cart_create(MPI_COMM_WORLD, ndims, dimsize, periods, reorder, &newcomm);

	if(myrank ==0){
		for (i=0; i<dimsize[0]; ++i){
			for (j=0; j<dimsize[1]; ++j){
				coords[0]=i;
				coords[1]=j;
				MPI_Cart_rank(newcomm, coords, &rank);
				printf("[%d] coords = %d, %d, rank = %d\n",myrank, coords[0], coords[1], rank);
			}
		}
	}


	MPI_Comm_size(newcomm, &newprocs);
	MPI_Comm_rank(newcomm, &newrank);

	printf("[%d]", myrank);
	printf(" newprocs=%d", newprocs);
	printf(" newrank=%d", newrank);
	printf("\n");
	MPI_Finalize();

}
