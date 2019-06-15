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
	dimsize[0]=4; dimsize[1]=4;
	periods[0] = 1; periods[1] = 0;
	reorder = 0;
	MPI_Cart_create(MPI_COMM_WORLD, ndims, dimsize, periods, reorder, &newcomm);

	MPI_Comm_size(newcomm, &newprocs);
	MPI_Comm_rank(newcomm, &newrank);

	MPI_Cart_coords(newcomm, newrank, ndims, coords);
	int direction=0;
	int displ=1;
	int source, dest;

	MPI_Cart_shift(newcomm, direction, displ, &source, &dest);
	printf("myrank= %d, coords=[%d][%d]\n", newrank, coords[0], coords[1]);
	printf("source= %d, dest= %d \n", source, dest);

	printf("\n");
	MPI_Finalize();

}
