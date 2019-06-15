#include "mpi.h"
#include <iostream>
#include <unistd.h>
#include <string>
#include <fstream>
#include <stdlib.h>

#define ALIVE '#'
#define DEAD '.'
#define SHIFT_ROW 0
#define SHIFT_COL 1
#define DISP 1

using namespace std;

int malloc2D(int ***array, int n, int m) {
	int i;
	/* allocate the n*m contiguous items */
	int *p = (int*)malloc(n*m*sizeof(int));
	if (!p) return -1;

	/* allocate the row pointers into the memory */
	(*array) = (int**)malloc(n*sizeof(int*));
	if (!(*array)) {
		free(p);
		return -1;
	}

	/* set up the pointers into the contiguous memory */
	for (i=0; i<n; i++)
		(*array)[i] = &(p[i*m]);

	return 0;
}

int free2D(int ***array) {
	/* free the memory - the first element of the array is at the start */
	free(&((*array)[0][0]));

	/* free the pointers into the memory */
	free(*array);

	return 0;
}
		
int doCalc(int **array,int **newarray, int rs, int rf, int cs, int cf, int move){
	int sum = 0;
	for (int x = rs; x < rf; x++){
		for (int y = cs; y < cf; y++){
			for(int i = -1; i <= 1; i++){
				for(int j = -1; j <=1; j++){
					if(i!=0||j!=0)
						sum+=array[x+move+i][y+move+j];
				}
			}
			if (array[x+move][y+move]==1 && (sum==2 || sum==3)) newarray[x+move][y+move]=1;
			else if (array[x+move][y+move]==1 && sum>3) newarray[x+move][y+move]=0;
			else if (array[x+move][y+move]==1 && sum<1) newarray[x+move][y+move]=0;
			else if (array[x+move][y+move]==0 && sum==3) newarray[x+move][y+move]=1;
			else newarray[x+move][y+move]=0;
			sum = 0;
		}
	}
}

int main(int argc, char* argv[])
{
	int nprocs, myrank;
	int my_cart_rank;
	int ndims=2;
	int dims[ndims], coord[ndims];
	int periods[ndims];
	int reorder, nrows, ncols;
	int tag, rc, N, generations, nghosts;
	int row_per_proc, col_per_proc, real_row, real_col;
	int row_padding, col_padding;
	int data[3];
	int full_0_row_proc, full_0_col_proc;
	int ierr;
	int source, dest;
	int nbr_up, nbr_down, nbr_left, nbr_right;
	int **theBoard;
	int **myslice;
	MPI_Request req[8];
	MPI_Status status;
	MPI_Comm comm2D;

	ierr = MPI_Init(&argc,&argv);
	ierr = MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
	ierr = MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

	// create cartesian topology for processes
	dims[0] = dims[1] = 0;
	MPI_Dims_create(nprocs, ndims, dims);
	nrows = dims[0];
	ncols = dims[1];
	
	if(myrank == 0)
		printf("PW[%d]/[%d]: PEdims = [%d x %d] \n", myrank, nprocs, dims[0], dims[1]);

	// create cartesian mapping
	periods[0] = periods[1] =0;	
	reorder = 1;
	ierr = MPI_Cart_create(MPI_COMM_WORLD, ndims, dims, periods, reorder, &comm2D);
	if(ierr != 0) printf("ERROR[%d] creating Cart\n", ierr);

	// find my coordinates in the cartesian communicator group
	MPI_Cart_coords(comm2D, myrank, ndims, coord);
	// use my coords to find my rank
	MPI_Cart_rank(comm2D, coord, &my_cart_rank);
//	printf("PW[%d]: my_cart_rank PCM[%d], my coords = (%d, %d)\n", myrank, my_cart_rank, coord[0], coord[1]);

	// find my neighbors. axis=0 : row , axis=1 : col
	MPI_Cart_shift(comm2D, SHIFT_ROW, DISP, &nbr_up, &nbr_down);
	MPI_Cart_shift(comm2D, SHIFT_COL, DISP, &nbr_left, &nbr_right);

	if (myrank==0){
		cin >> N >> generations >> nghosts;
		data[0]=N;
		data[1]=generations;
		data[2]=nghosts;
	} /* end of rank 0 */

	MPI_Bcast(data, 3, MPI_INT, 0, MPI_COMM_WORLD);

	N = data[0];
	generations = data[1];
	nghosts = data[2];
	row_padding = (N%nrows == 0 ? 0 :nrows - N%nrows);
	col_padding = (N%ncols == 0 ? 0 :ncols - N%ncols);
	row_per_proc = (N+row_padding)/nrows + 2*(nghosts + 1);
	col_per_proc = (N+col_padding)/ncols + 2*(nghosts + 1);
	real_row = (N+row_padding)/nrows;
	real_col = (N+col_padding)/ncols;
	full_0_row_proc = row_padding/real_row;
	full_0_col_proc = col_padding/real_col;
	int half = (real_row-2)/2;
/*
	if (myrank == 0){
		printf("--------------------------------------------\n");
		printf("N: %d  \tnprocs: %d  nghosts: %d\n", N, nprocs, nghosts);
		printf("nrows: %d\t\tncols: %d\n", nrows, ncols);
		printf("row_padding: %d\t\tcol_padding: %d\n", row_padding, col_padding);
		printf("row_real: %d\t\tcol_real: %d\n", real_row, real_col);
		printf("row_per_proc: %d\tcol_per_proc: %d\n", row_per_proc, col_per_proc);
		printf("--------------------------------------------\n");
	}
*/
	
	if (myrank == 0){

		int boardrow = N+row_padding;
		int boardcol = N+col_padding;
		
		malloc2D(&theBoard, boardrow, boardcol);

		/* Make a board */
		for (int i=0;i<N;i++)
		{
			string temp;
			cin >> temp;
			for (int j=0; j<N; j++){
				temp[j]==DEAD? theBoard[i][j]=0 : theBoard[i][j]=1;
			}
			for (int j=N; j<N+col_padding; j++){
				theBoard[i][j] = 0;
			}
		}
		for (int i=N;i<N+row_padding; i++){
			for (int j=0; j<N+col_padding; j++){
				theBoard[i][j] =0;
			}
		}
	}

	malloc2D(&myslice, row_per_proc, col_per_proc);
	
	/* Create new datatype to describe the subarray of the Board */
	int sizes[2] = 		{N+row_padding, N+col_padding};
	int subsizes[2] = 	{real_row, real_col};
	int starts[2] = 	{0,0};
	MPI_Datatype type1, subarrtype_board, subarrtype_slice;
	MPI_Type_create_subarray(2, sizes, subsizes, starts, MPI_ORDER_C, MPI_INT, &type1);
	MPI_Type_create_resized(type1, 0, real_col*sizeof(int), &subarrtype_board);
	MPI_Type_commit(&subarrtype_board);
	int sizes_[2] = 	{row_per_proc, col_per_proc};	
	int subsizes_[2] = 	{real_row,real_col};
	int starts_[2] = 	{1+nghosts, 1+nghosts};
	MPI_Type_create_subarray(2, sizes_, subsizes_, starts_, MPI_ORDER_C, MPI_INT, &subarrtype_slice);
	MPI_Type_commit(&subarrtype_slice);


	int *theBoardptr=NULL;
	if (myrank == 0)
		theBoardptr = &(theBoard[0][0]);

	/* scatter the array to all processors */
	int sendcounts[nprocs];
	int displs[nprocs];

	if (myrank ==0){
		for (int i = 0; i<nprocs; i++)
			sendcounts[i]=1;	
		int disp =0;
		for (int i = 0; i < nrows; i++){
			for (int j = 0; j < ncols; j++){
				displs[i*ncols+j] = disp;
				disp += 1;
			}
			disp += (real_row-1)*ncols;
		}
	}


	MPI_Scatterv(theBoardptr, sendcounts, displs, subarrtype_board, &(myslice[0][0]), 1, subarrtype_slice, 0, MPI_COMM_WORLD);
	
	
	int toleft[row_per_proc][nghosts+1]; 
	int toright[row_per_proc][nghosts+1]; 
	int fromleft[row_per_proc][nghosts+1];
	int fromright[row_per_proc][nghosts+1];
	int todown[nghosts+1][col_per_proc];
	int toup[nghosts+1][col_per_proc]; 
	int fromdown[nghosts+1][col_per_proc]; 
	int fromup[nghosts+1][col_per_proc];
	int row_buffer_size = (nghosts+1)*(col_per_proc);
	int col_buffer_size = (nghosts+1)*(row_per_proc);

	int **mynewslice;
	malloc2D(&mynewslice, row_per_proc, col_per_proc);

	for (int g=0; g < generations; g++)
	{	
		int sum=0;

		/* outermost process's ghost cell must always be 0. */
		if (nbr_up < 0){
			for (int x = 0; x < 1+nghosts ; x++){
				for (int y = 0 ; y < col_per_proc ; y++)
					myslice[x][y]=0;
			}
		}
		if (nbr_left < 0){
			for (int x = 0; x < row_per_proc ; x++){
				for (int y = 0 ; y < 1+nghosts ; y++)
					myslice[x][y]=0;
			}
		}
		if (nbr_down < 0){
			for (int x = 0; x < 1+nghosts ; x++){
				for (int y = 0 ; y < col_per_proc ; y++)
					myslice[x+row_per_proc-1-nghosts][y]=0;
			}
		}
		if (nbr_right < 0){
			for (int x = 0; x < row_per_proc ; x++){
				for (int y = 0 ; y < 1+nghosts ; y++)
					myslice[x][y+col_per_proc-1-nghosts]=0;
			}
		}

/*	
		if(myrank == 1){
			printf("PW[%2d] coord(%d,%d): NEWSLICE!\n", myrank, coord[0], coord[1]);
			for (int i = 0; i < row_per_proc ; i++){
				for(int j = 0 ; j < col_per_proc ; j++){
					cout << myslice[i][j];
				}
				cout << endl;
			}
			cout << endl;
		}
		//sleep(1);
*/

		/* Let's do communication!! 
		 * -----------------------------------------------------
		 * Row first	: P0 <-> P1 <-> P2 <-> ... <-> P[ncols-1]
		 * Col last		: p0 <-> P[ncols] <-> P[2*ncols] <-> ... <-> P[(nrows-1)ncols]
		 */

		if (g%(nghosts+1)==0){

			/* Send to right!!! */
			if (nbr_right > -1) // If I have right processor.
			{
				for (int j=0; j<row_per_proc; j++){
					for (int k =0; k<1+nghosts ; k++){
						toright[j][k]=myslice[j][col_per_proc -2*(1+nghosts)+k];
					}
				}
				MPI_Isend(toright, col_buffer_size, MPI_INT, nbr_right, 1, comm2D, &req[0]);
			}

			if (nbr_left > -1) // If I have left processor.
			{
				MPI_Irecv(fromleft, col_buffer_size, MPI_INT, nbr_left, 1, comm2D, &req[1]);

			}else{
				for (int j=0; j<row_per_proc; j++){
					for (int k=0; k<1+nghosts; k++)
						fromleft[j][k]=0; 
				}
			}

			/* Send to left!!! */
			if (nbr_left > -1) // If I have left processor.
			{
				for (int j=0; j<row_per_proc; j++){
					for (int k =0; k<1+nghosts ; k++){
						toleft[j][k]=myslice[j][k+1+nghosts];
					}
				}
				MPI_Isend(toleft, col_buffer_size, MPI_INT, nbr_left, 1, comm2D, &req[2]);
			}

			if (nbr_right > -1) // If I have right processor.
			{
				MPI_Irecv(fromright, col_buffer_size, MPI_INT, nbr_right, 1, comm2D, &req[3]);

			}else {
				for (int j=0; j<row_per_proc; j++){
					for (int k=0; k<1+nghosts; k++)
						fromright[j][k]=0; 
				}
			}

			/* Calculate Inner part of myslice that is independent with outer cells. */
			doCalc(myslice, mynewslice, 0, half, 0, real_col-2, 2+nghosts);

			/* Guarantee all colume side communication finished beform row side communication. */
			if(nbr_right > -1){
				MPI_Wait(&req[0], &status);
				MPI_Wait(&req[3], &status);
			}
			if(nbr_left > -1){
				MPI_Wait(&req[1], &status);
				MPI_Wait(&req[2], &status);
			}

			for (int i = 0; i < row_per_proc; i++){
				for (int j = 0; j < 1+nghosts ; j++){
					myslice[i][j] = fromleft[i][j];
				}
			}

			for (int i = 0; i < row_per_proc; i++){
				for (int j = 0; j < 1+nghosts ; j++){
					myslice[i][j+col_per_proc -1 - nghosts] = fromright[i][j];
				}
			}


			/* Send to down!! */
			if (nbr_down >-1)
			{
				for (int j=0; j<1+nghosts; j++){
					for (int k=0; k < col_per_proc; k++){
						todown[j][k]=myslice[row_per_proc-2*(1+nghosts)+j][k];
					}
				}

				MPI_Isend(todown, row_buffer_size, MPI_INT, nbr_down, 1, comm2D, &req[0]);

			}

			if (nbr_up>-1) // if I have up processor.
			{
				MPI_Irecv(fromup, row_buffer_size, MPI_INT, nbr_up, 1, comm2D, &req[1]);
			} else {
				for (int j=0; j<1+nghosts; j++){
					for (int k=0; k<col_per_proc; k++)
						fromup[j][k]=0; 
				}
			}

			/* Send to Up!!! */
			if (nbr_up > -1) // If I have up processor.
			{
				for (int j=0; j<1+nghosts; j++){
					for (int k =0; k<col_per_proc ; k++){
						toup[j][k]=myslice[1+nghosts+j][k];
					}
				}
				MPI_Isend(toup, row_buffer_size, MPI_INT, nbr_up, 1, comm2D, &req[2]);
			}

			if (nbr_down > -1) // If I have down processor.
			{
				MPI_Irecv(fromdown, row_buffer_size, MPI_INT, nbr_down, 1, comm2D, &req[3]);
			}else {
				for (int j=0; j<1+nghosts; j++){
					for (int k=0; k<col_per_proc; k++)
						fromdown[j][k]=0; 
				}
			}

			/* Calculate Inner part of myslice that is independent with outer cells. */
			doCalc(myslice, mynewslice, half, real_row-2, 0, real_col-2, 2+nghosts);
		
			/* Guarantee all colume side communication finished. */
			if(nbr_down> -1){
				MPI_Wait(&req[0], &status);
				MPI_Wait(&req[3], &status);
			}
			if(nbr_up > -1){
				MPI_Wait(&req[1], &status);
				MPI_Wait(&req[2], &status);
			}

			// Copy received slice into my slice
			for (int i = 0; i < 1+nghosts; i++){
				for (int j = 0; j < col_per_proc ; j++){
					myslice[i][j] = fromup[i][j];
				}
			}

			for (int i = 0; i < 1+nghosts; i++){
				for (int j = 0; j < col_per_proc ; j++){
					myslice[i+ row_per_proc - 1 - nghosts][j] = fromdown[i][j];
				}
			}

			/* calculate rest part of myslice */
			doCalc(myslice, mynewslice, 1, 2+nghosts, 1+nghosts, col_per_proc-1-nghosts, 0);
			doCalc(myslice, mynewslice, row_per_proc-2-nghosts, row_per_proc-1, 1+nghosts, col_per_proc-1-nghosts, 0);
			doCalc(myslice, mynewslice, 1, row_per_proc-1, 1, nghosts+2, 0);
			doCalc(myslice, mynewslice, 1, row_per_proc-1, col_per_proc-2-nghosts, col_per_proc-1, 0);
		 	/* end of calculation of rest part. */

		}	// end of nghosts if
		else
		{
			doCalc(myslice, mynewslice, 1, row_per_proc-1, 1, col_per_proc - 1, 0);
		}

		
		// copy new slice onto myslice
		for (int x=1; x< row_per_proc - 1; x++)
			for (int y= 1; y< col_per_proc -1; y++)
				myslice[x][y]=mynewslice[x][y];
	

	} // end of generation loop

	/* Let's collect result!! */

	MPI_Gatherv(&(myslice[0][0]), 1, subarrtype_slice, theBoardptr, sendcounts, displs, subarrtype_board, 0, MPI_COMM_WORLD);

	if(myrank == 0){
		for(int i = 0; i < N; i++){
			for(int j = 0; j < N; j++){
				cout << (theBoard[i][j] == 0 ? DEAD : ALIVE);
			}
			cout << endl;
		}
		free2D(&theBoard);
	}

	free2D(&myslice);
	MPI_Type_free(&type1);
	MPI_Type_free(&subarrtype_board);
	MPI_Type_free(&subarrtype_slice);

	MPI_Finalize();
}
