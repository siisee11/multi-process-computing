#include <stdio.h>
#include "mpi.h"
#include <math.h>
#define NBRILO 10
#define NBRIHI 11
#define NBRJLO 12
#define NBRJHI 13
#define SHIFT_ROW 0
#define SHIFT_COL 1
#define DISP 1
int main(int argc, char *argv[]) {
	int ndims=2, size, my_rank, reorder, my_cart_rank, ierr, errs, nrows, ncols, source, dest, nbr_i_lo, nbr_i_hi;
	MPI_Comm comm2D;
	int dims[ndims],coord[ndims];
	int wrap_around[ndims];
	/* start up initial MPI environment */
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
	/* process command line arguments*/
	if (argc == 3) {
		nrows = atoi (argv[1]);
		ncols = atoi (argv[2]);
		dims[0] = nrows; /* number of rows */
		dims[1] = ncols; /* number of columns */
		if( (nrows*ncols) != size) {
			if( my_rank ==0) printf("ERROR: nrows*ncols)=%d * %d = %d != size\n", nrows, ncols, nrows*ncols,size);
			MPI_Finalize(); exit(0);
		}
	} else {
		nrows=ncols=(int)sqrt(size);
		dims[0]=dims[1]=0;
	}

	/*************************************************************/
	/* create cartesian topology for processes */
	/*************************************************************/
	MPI_Dims_create(size, ndims, dims);
	if(my_rank==0)
		printf("PW[%d], CommSz[%d%]: PEdims = [%d x %d] \n",my_rank,size,dims[0],dims[1]);
	/* create cartesian mapping */
	wrap_around[0] = wrap_around[1] = 0; /* periodic shift is .false. */
	reorder = 1;
	ierr =0;
	ierr = MPI_Cart_create(MPI_COMM_WORLD, ndims, dims,
			wrap_around, reorder, &comm2D);
	if(ierr != 0) printf("ERROR[%d] creating CART\n",ierr);
	/* find my coordinates in the cartesian communicator group */
	MPI_Cart_coords(comm2D, my_rank, ndims, coord);
	/* use my cartesian coordinates to find my rank in cartesian group*/
	MPI_Cart_rank(comm2D, coord, &my_cart_rank);
	/* get my neighbors; axis is coordinate dimension of shift */
	/* axis=0 ==> shift along the rows: P[my_row-1]: P[me] : P[my_row+1] */
	/* axis=1 ==> shift along the columns P[my_col-1]: P[me] : P[my_col+1] */
	MPI_Cart_shift( comm2D, SHIFT_ROW, DISP, &nbr_i_lo, &nbr_i_hi );
	MPI_Cart_shift( comm2D, SHIFT_COL, DISP, &nbr_j_lo, &nbr_j_hi );
	printf( "PW[%2d] Coord(%d,%d): SHIFT_DIM[%d], Shift=%d: nbr_lo[%2d] P[%2d] nbr_hi[%2d] \n",
			my_rank, coord[0], coord[1],SHIFT_ROW, DISP,nbr_i_lo ,my_rank,nbr_i_hi );
	printf( "PW[%2d] Coord(%d,%d): SHIFT_DIM[%d], Shift=%d: nbr_lo[%2d] P[%2d] nbr_hi[%2d] \n",
			my_rank, coord[0], coord[1],SHIFT_COL, DISP,nbr_j_lo ,my_rank,nbr_j_hi );
	fflush(stdout);
	MPI_Comm_free( &comm2D );
	MPI_Finalize();
}
