/* mpi-cart-2D-row.c --
 * creates a 1D row cartesian communicator
 *
 * Written by Mary Thomas
 * - Updated Mar, 2015
 *
 * Based loosely on code from Pacheco'97,
 * Chap 7, PPMPI
 */
#include <stdio.h>
#include "mpi.h"
#include <math.h>

int main(int argc, char *argv[]) {
	int ndims=2;
	int p;
	int my_rank;
	int dims[ndims],coord[ndims];
	int wrap_around[ndims];
	int reorder;
	int my_cart_rank;
	int ierr;
	int nrows, ncols;
	char name[200], nameout[200], rname[100];
	int cnt, rlen;
	MPI_Comm comm2D;
	/* for row split */
	int my_row_rank;
	int free_coords[ndims];
	int row_coord[ndims];
	int row_val;
	int row_sum=0, row_root;
	MPI_Comm comm1D_row;

	/*************************************************************/
	/* start up initial MPI environment: MPI_COMM_WORLD */
	/*************************************************************/
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &p);
	MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
	/* turn on error handling */
	MPI_Errhandler_set(MPI_COMM_WORLD, MPI_ERRORS_RETURN);
	/* get the name of communicator */
	nameout[0] = 0;
	MPI_Comm_get_name( MPI_COMM_WORLD, nameout, &rlen );
	if(my_rank == 0) printf( "Name of comm world is: %s\n", nameout );
	fflush(stdout);
	/* process command line arguments*/
	if (argc == 3) {
		nrows = atoi (argv[1]);
		ncols = atoi (argv[2]);
		dims[0] = nrows; /* number of rows */
		dims[1] = ncols; /* number of columns */
		if( (nrows*ncols) != p) {
			if( my_rank ==0) printf("ERROR: nrows*ncols)=%d * %d = %d != p\n",
					nrows, ncols, nrows*ncols,p);
			MPI_Finalize();
			exit(0);
		}
	} else {
		nrows=ncols=(int)sqrt(p);
		dims[0]=dims[1]=0;
	}

	/*************************************************************/
	/* create a 2D cartesian communicator: comm2D */
	/*************************************************************/
	/* set dimensions for cartesian topology for processes */
	MPI_Dims_create(p, ndims, dims);
	if( my_rank == 0 )
		printf("PW[%d]/[%d%]: PEdims = [%d x %d] \n",my_rank,p,dims[0],dims[1]);
	/* create cartesian mapping */
	wrap_around[0] = wrap_around[1] = 0; // set periodicity
	reorder = 1;
	ierr =0;
	ierr = MPI_Cart_create(MPI_COMM_WORLD, ndims, dims,
			wrap_around, reorder, &comm2D);
	if(ierr != 0) printf("ERROR[%d] creating CART\n",ierr);
	/* set the name of cartesian communicator */
	strcpy( name, "comm2D" );
	MPI_Comm_set_name( comm2D, name );
	nameout[0] = 0;
	/* get the name of communicator */
	MPI_Comm_get_name( comm2D, nameout, &rlen );
	if(my_rank == 0) printf( "Name of comm world is: %s\n", nameout );
	fflush(stdout);
	/* find my coordinates in the cartesian communicator group */
	MPI_Cart_coords(comm2D, my_rank, ndims, coord);
	/* use my cartesian coordinates to find my rank in cartesian group*/
	MPI_Cart_rank(comm2D, coord, &my_cart_rank);

	/*************************************************************/
	/* split comm2D into a 1D row-based communicator: comm1D_row */
	/*************************************************************/
	free_coords[0] = 0; /* rows */; free_coords[1] = 1; /* cols */
	MPI_Cart_sub(comm2D, free_coords, &comm1D_row);
	/* get my_row_rank in my comm1D_row group */
	MPI_Comm_rank(comm1D_row, &my_row_rank);
	/* use my_row_rank to find my coordinates in my comm1D_row group */
	MPI_Cart_coords(comm1D_row, my_row_rank, 1, &row_coord);
	/* set the name of cartesian communicator */
	row_root=-999;
	if(coord[1] == 0) row_root = my_rank; /* use rows root rank for ID of row */
	MPI_Bcast(&row_root, 1, MPI_INT, 0, comm1D_row);
	/* set the name of cartesian communicator and build a unique name for it */
	sprintf(rname, "%d", row_root);
	strcpy( name, "comm1D_row"); strcat( name, "_"); strcat( name, rname);
	MPI_Comm_set_name( comm1D_row, name );
	/* get the name of communicator */
	MPI_Comm_get_name( comm1D_row, nameout, &rlen );
	/* test row com: row_sum processor ranks across rows */
	if (coord[1] == 0)
		row_val = coord[0];
	else
		row_val = -1;
	MPI_Bcast(&row_val, 1, MPI_INT, 0, comm1D_row);
	MPI_Reduce(&my_rank, &row_sum, 1, MPI_INT, MPI_SUM, 0, comm1D_row);
	MPI_Bcast(&row_sum, 1, MPI_INT, 0, comm1D_row);
	printf("PW[%d]: PCM[%d%]: C(%d,%d), RowComm: Root PW[ %d], Name=%s,RSum=%d \n",
			my_rank, my_row_rank , coord[0], coord[1], row_root, nameout,row_sum);
	MPI_Comm_free(comm2D);(comm1D_row);
	MPI_Finalize();
} /* main */
