#include <stdio.h>
#include <math.h>
#include <mpi.h>

int main(int argc, char** argv) {
	long n, i;
	double sum, step, pi, x;
	int cont;
	int myrank, nprocs, ierr, is, ie;
	double t_start, t_finish;
	MPI_Status status;
	cont =1; 
	ierr = MPI_Init(&argc, &argv);
	ierr = MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
	ierr = MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

	double tsum[nprocs];
	
	while (cont) {
		if(myrank==0) {
			printf("Enter the Number of Intervals : (n<1 : quits)\n");
			scanf("%ld", &n);
			t_start = MPI_Wtime();
		}
		
		ierr = MPI_Bcast(&n, 1, MPI_LONG, 0, MPI_COMM_WORLD);

		if(n <=0){
			cont = 0;
			break;
		}
		step = 1.0/(double)n;
		sum = 0.0;
		tsum[myrank] = 0.0;

		for(i = myrank; i < n; i = i+nprocs){
			x = ((double)i + 0.5)*step;
			sum = sum+4.0/(1.0 + x*x);
		}
		
		MPI_Gather(&sum,1,MPI_DOUBLE, tsum, 1, MPI_DOUBLE, 0 , MPI_COMM_WORLD);

		if(myrank == 0){
			for(i = 1 ; i < nprocs; i++)
				tsum[0] = tsum[0]+tsum[i];
			pi = step*tsum[0];
			t_finish = MPI_Wtime();

			printf("-------------------------------------------------------\n");
			printf("[%.5f] PI = %.15f (Error = %E)\n",t_finish-t_start, pi, fabs(acos(-1.0)-pi));
			printf("-------------------------------------------------------\n");
		}else{
			ierr = MPI_Send(&sum, 1, MPI_DOUBLE, 0, 1, MPI_COMM_WORLD);
		}
	}
	ierr=MPI_Finalize();
}
		
