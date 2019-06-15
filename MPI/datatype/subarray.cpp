#include<stdio.h>
#include<mpi.h>

int main(int argc, char **argv){
	int ndims = 2;
	int ibuf1[6][8];
	int array_of_sizes[ndims],array_of_subsizes[ndims], array_of_starts[ndims];
	int i,j,myrank;

	MPI_Datatype newtype;
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

	if(myrank==0){ 
		for(i=0; i<6; i++){
			for(j=0; j<8; j++)
					ibuf1[i][j]=i+2;
		}
	}else{
		for(i=0; i<6; i++){
			for(j=0; j<8; j++)
					ibuf1[i][j]=0;
		}
	}

	array_of_sizes[0]=6; array_of_sizes[1] =8;
	array_of_subsizes[0]=3; array_of_subsizes[1]=5;
	array_of_starts[0]=2; array_of_starts[1]=1;

	MPI_Type_create_subarray(ndims, array_of_sizes, array_of_subsizes, array_of_starts, MPI_ORDER_C, MPI_INT, &newtype);

	MPI_Type_commit(&newtype);
	MPI_Bcast(ibuf1, 1, newtype, 0, MPI_COMM_WORLD);

	if(myrank !=0) {
		printf("I am : %d processor \n", myrank);
		for(i=0;i<6;i++){
			for(j=0;j<8;j++){
				printf(" %d", ibuf1[i][j]);
			}
		printf("\n");
		}
	}
	MPI_Finalize();
}
