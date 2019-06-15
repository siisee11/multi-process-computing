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
//	int **myslice;
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
	/*
	if(myrank == 0)
		printf("PW[%d]/[%d]: PEdims = [%d x %d] \n", myrank, nprocs, dims[0], dims[1]);
*/
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
//	printf("PW[%d] Coord(%d,%d): nbr_up[%2d] P[%2d] nbr_down[%2d]\n",myrank, coord[0], coord[1], nbr_up, myrank, nbr_down);
//	printf("PW[%d] Coord(%d,%d): nbr_left[%2d] P[%2d] nbr_right[%2d]\n",myrank, coord[0], coord[1], nbr_left, myrank, nbr_right);

//	sleep(1);
//	MPI_Barrier(MPI_COMM_WORLD);

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
	+2;
	real_row = (N+row_padding)/nrows;
	real_col = (N+col_padding)/ncols;
	full_0_row_proc = row_padding/real_row;
	full_0_col_proc = col_padding/real_col;
	int half = (real_row-2)/2;
/*
	if (myrank == 0){
		printf("--------------------------------------------\n");
		printf("N: %d  \tnprocs: %d  nghosts: %d\n", N, nprocs, nghosts);
		printf("row_padding: %d\t\tcol_padding: %d\n", row_padding, col_padding);
		printf("row_real: %d\t\tcol_real: %d\n", real_row, real_col);
		printf("row_per_proc: %d\t\tcol_per_proc: %d\n", row_per_proc, col_per_proc);
		printf("full_0_row_proc: %d\tfull_0_col_proc: %d\n", full_0_row_proc, full_0_col_proc);
		printf("--------------------------------------------\n");
	}
 */
/*
	// alloc myslice size
	myslice = (int **)malloc(row_per_proc * sizeof(int *));
	for (int i = 0; i < row_per_proc ; i++)
		myslice[i] = (int *)malloc(col_per_proc * sizeof(int));
*/
	int myslice[row_per_proc][col_per_proc];


	if (myrank == 0){

		int boardrow = N+row_padding+2*(nghosts+1);
		int boardcol = N+col_padding+2*(nghosts+1);
		int theBoard[boardrow][boardcol]={0,};
		/*
		int **theBoard = (int **)malloc((N+row_padding)*sizeof(int *));
		for (int i = 0; i < N+row_padding; i++){
			theBoard[i] = (int *)malloc(N * sizeof(int));
		}
		*/

		/* Make a board */
		for (int i=0; i<1+nghosts; i++){
			for(int j=0; j <boardcol ; j++)
				theBoard[i][j]=0;
		}
		for (int i=1+nghosts;i<N+nghosts+1;i++)
		{
			string temp;
			cin >> temp;
			for (int j=0; j<1+nghosts; j++){
				theBoard[i][j]=0;
			}
			for (int j=1+nghosts; j<N+1+nghosts; j++){
				temp[j-1-nghosts]==DEAD? theBoard[i][j]=0 : theBoard[i][j]=1;
			}
			for (int j=N+nghosts+1; j<N+col_padding+2*(1+nghosts); j++){
				theBoard[i][j] = 0;
			}
		}
		for (int i=N+1+nghosts;i<N+row_padding+2*(1+nghosts);i++){
			for (int j=0; j<N+col_padding+2*(1+nghosts); j++){
				theBoard[i][j] =0;
			}
		}

		/* print a board */
		/*
		printf("the Board [%d][%d]:\n", boardrow, boardcol);
		cout<<endl;
		for (int i = 0 ; i < boardrow; i++){
			for(int j = 0 ; j < boardcol; j++){
				cout << theBoard[i][j];
			}
			cout << endl;
		}
		cout << "---------------------------------------" << endl;
*/
		for (int row = nrows - 1; row>=0; row--)
		{
			for (int col = ncols - 1; col>=0; col--){
				for (int x=0; x < row_per_proc ; x++){
					for (int y=0; y < col_per_proc; y++){
						myslice[x][y] = theBoard[row*real_row+x][col*real_col+y];
					}
				}

				if (row*ncols + col != 0){
					MPI_Isend(myslice, row_per_proc*col_per_proc, MPI_INT, row*ncols+col, 1, MPI_COMM_WORLD, &req[0]);
					MPI_Wait(&req[0], &status);
				}
			}
		}

	}else{
		MPI_Irecv(myslice, row_per_proc*col_per_proc, MPI_INT, 0, 1, MPI_COMM_WORLD, &req[0]);
		MPI_Wait(&req[0], &status);
	}

	int todown[nghosts+1][col_per_proc];
	int toup[nghosts+1][col_per_proc]; 
	int toleft[row_per_proc][nghosts+1]; 
	int toright[row_per_proc][nghosts+1]; 
	int fromdown[nghosts+1][col_per_proc]; 
	int fromup[nghosts+1][col_per_proc];
	int fromleft[row_per_proc][nghosts+1];
	int fromright[row_per_proc][nghosts+1];
	int row_buffer_size = (nghosts+1)*(col_per_proc);
	int col_buffer_size = (nghosts+1)*(row_per_proc);

/*	
	MPI_Barrier(MPI_COMM_WORLD);
	sleep(0.5);
	if(myrank==1){
		cout<<endl;
		cout<< "Myslice [" << myrank << "]" << endl;
		for (int i = 0 ; i < row_per_proc; i++){
			for(int j = 0 ; j < col_per_proc; j++){
				cout << myslice[i][j];
			}
			cout << endl;
		}
		cout << "---------------------------------------" << endl;
	
	}
	MPI_Barrier(MPI_COMM_WORLD);
	sleep(0.5);
*/
	for (int g=0; g < generations; g++)
	{	
		int sum=0;
		int mynewslice[row_per_proc][col_per_proc]={0,};
		/* Let's do communication!! 
		 * -----------------------------------------------------
		 * Row first	: P0 <-> P1 <-> P2 <-> ... <-> P[ncols-1]
		 * Col last		: p0 <-> P[ncols] <-> P[2*ncols] <-> ... <-> P[(nrows-1)ncols]
		 */
		if (g%(nghosts+1)==0 && g!=0){

//			if(myrank ==0 ) cout << "-------------------DOCOMM--------------------" << endl;
			/* Send to right!!! */
			if (nbr_right > -1) // If I have right processor.
			{
				for (int j=0; j<row_per_proc; j++){
					for (int k =0; k<1+nghosts ; k++){
						toright[j][k]=myslice[j][col_per_proc -2*(1+nghosts)+k];
					}
				}
/*
				if (myrank == 0) {
					printf("PW[%2d] Coord(%d,%d): Ready to send to right P[%2d]\n", myrank, coord[0], coord[1], nbr_right);
					for (int j=0; j<row_per_proc; j++){
						for (int k=0; k<1+nghosts ; k++){
							cout << toright[j][k];
						}
						cout<<endl;
					}
					cout<<endl<<endl;
				}
*/
				MPI_Isend(toright, col_buffer_size, MPI_INT, nbr_right, 1, comm2D, &req[0]);
//				cout << myrank <<" : Send to up" << endl; 
			}

			if (nbr_left > -1) // If I have left processor.
			{
				MPI_Irecv(fromleft, col_buffer_size, MPI_INT, nbr_left, 1, comm2D, &req[1]);

			}else {
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
/*
				if (myrank == 1) {
					printf("PW[%2d] Coord(%d,%d): Ready to send to left P[%2d]\n", myrank, coord[0], coord[1], nbr_left);
					for (int j=0; j<row_per_proc; j++){
						for (int k=0; k<1+nghosts ; k++){
							cout << toleft[j][k];
						}
						cout<<endl;
					}
					cout<<endl<<endl;
				}
*/
				MPI_Isend(toleft, col_buffer_size, MPI_INT, nbr_left, 1, comm2D, &req[2]);
//				cout << myrank <<" : Send to up" << endl; 
			}

			if (nbr_right > -1) // If I have right processor.
			{
				MPI_Irecv(fromright, col_buffer_size, MPI_INT, nbr_right, 1, comm2D, &req[3]);
//				MPI_Recv(fromright, col_buffer_size, MPI_INT, nbr_right, 1, comm2D, &status);

			}else {
				for (int j=0; j<row_per_proc; j++){
					for (int k=0; k<1+nghosts; k++)
						fromright[j][k]=0; 
				}
			}

			/* Calculate Inner part of myslice that is independent with outer cells. */
			for (int x = 0; x < half; x++){
				for (int y = 0; y < real_col-2; y++){
					for(int i = -1; i <= 1; i++){
						for(int j = -1; j <=1; j++){
							if(i!=0||j!=0)
								sum+=myslice[x+2+nghosts+i][y+2+nghosts+j];
						}
					}
					if (myslice[x+2+nghosts][y+2+nghosts]==1 && (sum==2 || sum==3)) mynewslice[x+2+nghosts][y+2+nghosts]=1;
					else if (myslice[x+2+nghosts][y+2+nghosts]==1 && sum>3) mynewslice[x+2+nghosts][y+2+nghosts]=0;
					else if (myslice[x+2+nghosts][y+2+nghosts]==1 && sum<1) mynewslice[x+2+nghosts][y+2+nghosts]=0;
					else if (myslice[x+2+nghosts][y+2+nghosts]==0 && sum==3) mynewslice[x+2+nghosts][y+2+nghosts]=1;
					else mynewslice[x+2+nghosts][y+2+nghosts]=0;
					sum = 0;
				}
			}
	
			

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
//				printf("PW[%2d] Coord(%d,%d): Finished to send to down P[%2d]\n", myrank, coord[0], coord[1], nbr_down);

			} else {
				for (int j = 0 ; j < 1+nghosts ; j++) 
					for (int k = 0 ; k < col_per_proc ; k++)
						fromdown[j][k]=0; 
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
//				cout << myrank <<" : Send to up" << endl; 
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
			for (int x = half; x < real_row-2; x++){
				for (int y = 0; y < real_col-2; y++){
					for(int i = -1; i <= 1; i++){
						for(int j = -1; j <=1; j++){
							if(i!=0||j!=0)
								sum+=myslice[x+2+nghosts+i][y+2+nghosts+j];
						}
					}
					if (myslice[x+2+nghosts][y+2+nghosts]==1 && (sum==2 || sum==3)) mynewslice[x+2+nghosts][y+2+nghosts]=1;
					else if (myslice[x+2+nghosts][y+2+nghosts]==1 && sum>3) mynewslice[x+2+nghosts][y+2+nghosts]=0;
					else if (myslice[x+2+nghosts][y+2+nghosts]==1 && sum<1) mynewslice[x+2+nghosts][y+2+nghosts]=0;
					else if (myslice[x+2+nghosts][y+2+nghosts]==0 && sum==3) mynewslice[x+2+nghosts][y+2+nghosts]=1;
					else mynewslice[x+2+nghosts][y+2+nghosts]=0;
					sum = 0;
				}
			}
		
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
			for (int x = 1; x <= 1+nghosts; x++){
				for (int y = 1+nghosts; y <= col_per_proc-2-nghosts; y++){
					for(int i = -1; i <= 1; i++){
						for(int j = -1; j <=1; j++){
							if(i!=0||j!=0)
								sum+=myslice[x+i][y+j];
						}
					}
					if (myslice[x][y]==1 && (sum==2 || sum==3)) mynewslice[x][y]=1;
					else if (myslice[x][y]==1 && sum>3) mynewslice[x][y]=0;
					else if (myslice[x][y]==1 && sum<1) mynewslice[x][y]=0;
					else if (myslice[x][y]==0 && sum==3) mynewslice[x][y]=1;
					else mynewslice[x][y]=0;
					sum = 0;
				}
			}
			for (int x = row_per_proc-2-nghosts; x <= row_per_proc-2; x++){
				for (int y = 1+nghosts; y <= col_per_proc-2-nghosts; y++){
					for(int i = -1; i <= 1; i++){
						for(int j = -1; j <=1; j++){
							if(i!=0||j!=0)
								sum+=myslice[x+i][y+j];
						}
					}
					if (myslice[x][y]==1 && (sum==2 || sum==3)) mynewslice[x][y]=1;
					else if (myslice[x][y]==1 && sum>3) mynewslice[x][y]=0;
					else if (myslice[x][y]==1 && sum<1) mynewslice[x][y]=0;
					else if (myslice[x][y]==0 && sum==3) mynewslice[x][y]=1;
					else mynewslice[x][y]=0;
					sum = 0;
				}
			}
			for (int x = 1; x <= row_per_proc-2; x++){
				for (int y = 1; y <= nghosts+1; y++){
					for(int i = -1; i <= 1; i++){
						for(int j = -1; j <=1; j++){
							if(i!=0||j!=0)
								sum+=myslice[x+i][y+j];
						}
					}
					if (myslice[x][y]==1 && (sum==2 || sum==3)) mynewslice[x][y]=1;
					else if (myslice[x][y]==1 && sum>3) mynewslice[x][y]=0;
					else if (myslice[x][y]==1 && sum<1) mynewslice[x][y]=0;
					else if (myslice[x][y]==0 && sum==3) mynewslice[x][y]=1;
					else mynewslice[x][y]=0;
					sum = 0;
				}
			}
			for (int x = 1; x <= row_per_proc-2; x++){
				for (int y = col_per_proc-2-nghosts; y <= col_per_proc-2; y++){
					for(int i = -1; i <= 1; i++){
						for(int j = -1; j <=1; j++){
							if(i!=0||j!=0)
								sum+=myslice[x+i][y+j];
						}
					}
					if (myslice[x][y]==1 && (sum==2 || sum==3)) mynewslice[x][y]=1;
					else if (myslice[x][y]==1 && sum>3) mynewslice[x][y]=0;
					else if (myslice[x][y]==1 && sum<1) mynewslice[x][y]=0;
					else if (myslice[x][y]==0 && sum==3) mynewslice[x][y]=1;
					else mynewslice[x][y]=0;
					sum = 0;
				}
			} /* end of calculation of rest part. */



			/* Print after all communication. */
			/*
			if (myrank == 0){
				cout << "---------------------------------------------------------" << endl;
				cout << "After all commutication" << endl;
				for (int i = 0 ;  i < row_per_proc ; i++){
					for (int j = 0 ; j < col_per_proc ; j++){
						cout << myslice[i][j];
					}
					cout << endl;
				}
				cout << "---------------------------------------------------------" << endl;
				cout << endl;
			}
			//sleep(0.5);
*/
		}	// end of nghosts if
		else
		{

			for (int x = 1; x < row_per_proc - 1; x++){
				for (int y = 1; y < col_per_proc - 1; y++){
					for(int i = -1; i <= 1; i++){
						for(int j = -1; j <=1; j++){
							if(i!=0||j!=0)
								sum+=myslice[x+i][y+j];
						}
					}
					if (myslice[x][y]==1 && (sum==2 || sum==3)) mynewslice[x][y]=1;
					else if (myslice[x][y]==1 && sum>3) mynewslice[x][y]=0;
					else if (myslice[x][y]==1 && sum<1) mynewslice[x][y]=0;
					else if (myslice[x][y]==0 && sum==3) mynewslice[x][y]=1;
					else mynewslice[x][y]=0;
					sum = 0;
				}
			}
		}

		
		// copy new slice onto myslice
		for (int x=1; x< row_per_proc - 1; x++)
			for (int y= 1; y< col_per_proc -1; y++)
				myslice[x][y]=mynewslice[x][y];
/*
		// last process's (who has value) ghost cell and padding must always be 0.
		if (myrank == nprocs -1 -full_0_row_proc){
			for (int x = row_per_proc - 1 -nghosts - row_padding%real_row; x < row_per_proc ; x++){
				for (int y = 0 ; y < col_per_proc; y++)
					myslice[x][y]=0;
			}
		}
*/
		// outermost process's ghost cell must always be 0.
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
		if(myrank == 0){
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
	} // end of generation loop


//	MPI_Barrier(MPI_COMM_WORLD);
	if (myrank==0) 
	{
		char theBoard[N][N];
		int aBoard[row_per_proc][col_per_proc];
		for (int x=0; x< real_row; x++) 
		{
			for (int y=0; y< real_col; y++) {
				theBoard[x][y] = (myslice[x+1+nghosts][y+1+nghosts] == 0? DEAD : ALIVE);
			}
		}

		for (int i=1; i<nprocs - full_0_row_proc; i++)
		{
			int row = i/ncols;
			int col = i%ncols;
			MPI_Recv(&aBoard, row_per_proc*col_per_proc, MPI_INT, i, 1, MPI_COMM_WORLD, &status); //receive all others'
			/*
			if (i == nprocs-1-full_0_row_proc){
				real_row = real_row - row_padding%real_row;
			}
			*/
			for (int x=0; x<real_row; x++)
			{
				for (int y=0; y<real_col; y++)
					theBoard[x+real_col*row][y+real_row*col] = (aBoard[x+1+nghosts][y+1+nghosts] == 0? DEAD : ALIVE);
			}
		}

		for(int i = 0; i < N; i++){
			for(int j =0; j< N; j++){
				cout << theBoard[i][j];
			}
			cout << endl;
		}
	}
	else MPI_Send(&myslice, row_per_proc*col_per_proc, MPI_INT, 0,1, MPI_COMM_WORLD);

	MPI_Finalize();
}
