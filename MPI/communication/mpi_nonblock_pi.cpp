#include <stdio.h>
#include <math.h>
#include <mpi.h>

int main(int argc, char** argv) {
	long n, i;
	double sum, step, pi, tsum, x;
	int cont;
	int myrank, nprocs, ierr, is, ie;
	double t_start, t_finish;
	MPI_Status status;
	MPI_Request req;

	cont =1; 
	ierr = MPI_Init(&argc, &argv);
	ierr = MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
	ierr = MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
	
	while (cont) {
		if(myrank==0) {
			printf("Enter the Number of Intervals : (n<1 : quits)\n");
			scanf("%ld", &n);
			t_start = MPI_Wtime();
			for(i = 1; i < nprocs ; ++i){
				ierr=MPI_Isend(&n, 1, MPI_LONG, i, 1, MPI_COMM_WORLD, &req);
				ierr=MPI_Wait(&req, &status);
			}
		}else{
			ierr =MPI_Irecv(&n, 1, MPI_LONG, 0, 1, MPI_COMM_WORLD, &req);
			ierr = MPI_Wait(&req, &status);
		}
		
		if(n <=0){
			cont = 0;
			break;
		}
		step = 1.0/(double)n;
		sum = 0.0;
		tsum = 0.0;

		for(i = myrank; i < n; i = i+nprocs){
			x = ((double)i + 0.5)*step;
			sum = sum+4.0/(1.0 + x*x);
		}
		
		if(myrank == 0){
			tsum = sum;
			for(i=1; i<nprocs; i++){
				ierr = MPI_Irecv(&sum, 1, MPI_DOUBLE, i, 1, MPI_COMM_WORLD, &req);
				ierr = MPI_Wait(&req, &status);
				tsum = tsum + sum;
			}
			pi = step*tsum;
			t_finish = MPI_Wtime();

			printf("-------------------------------------------------------\n");
			printf("[%.5f] PI = %.15f (Error = %E)\n",t_finish-t_start, pi, fabs(acos(-1.0)-pi));
			printf("-------------------------------------------------------\n");
		}else{
			ierr = MPI_Isend(&sum, 1, MPI_DOUBLE, 0, 1, MPI_COMM_WORLD, &req);
			ierr = MPI_Wait(&req, &status);
		}
	}
	ierr=MPI_Finalize();
}
		
