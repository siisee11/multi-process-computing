#include <getopt.h>
#include <limits.h>
#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <omp.h>

#define SWAP(a, b) do { \
	REAL temp; \
	temp = a;  \
	a = b;     \
	b = temp;  \
} while (0)

#define MIN(a, b) ( (a) < (b) ? (a) : (b) )

#define REAL double
#define BS 32				// block size
#define TOLERANCE 10.0e-0


int n;
int bs;			
int BNAM;		// # of bloks in a row.
int rand_seed;
int d_flag;
int padding;
int NP;			// n with padding.
REAL *A;
REAL *A_save;
REAL *L;
REAL *U;

// options
bool debug_mode = false;
bool check_mode = false;

void lu(REAL* , REAL* , int );
REAL* mult(REAL*, REAL*);

int compare_REAL(REAL x, REAL y)
{
	REAL diff = x - y;
	if (fabs(diff) <= TOLERANCE)
		return 0;

	return (diff > 0) ? 1 : -1;
}

void transpose(REAL *, REAL *, int);
void inverseL(REAL *, REAL *, int);
void inverseU(REAL *, REAL *, int);

/*
 * Prints a matrix to standard output.
 * 
 */
void print_matrix(REAL *mat, int size, int cols)
{
	for (int row = 0; row < size; row++) {
		for (int col = 0; col < size; col++) {
			printf("%.4e ", mat[row*cols + col]);
		}
		printf("\n");
	}
}

void matrix_init(){
	A_save = (REAL*)calloc(NP*NP,sizeof(REAL));
	L = (REAL*)calloc(NP*NP, sizeof(REAL));
	U = (REAL*)calloc(NP*NP, sizeof(REAL));
	A = (REAL*)calloc(NP*NP, sizeof(REAL));
	REAL r;
	
	srand(rand_seed);

	// generate random matrix entries
	
	for (int row = 0; row < n; row++) {
		for (int col = 0; col < n; col++) {
			//r = rand()%10+1;
			r = rand();
			A[row*NP+col] = r;
			A_save[row*NP+col] = r;
		}
		L[row*NP + row] = 1;
	}

	for (int row = n ; row < NP; row++){
		for (int col = 0 ; col < NP; col++) {
			A[row*NP+col] = 1;
			A_save[row*NP+col] = 1;
		}
	}
	
}

void put_L(REAL* L_sub, int row, int col){
	for(int i = row*bs, i_ = 0 ; i_ < bs ; i++, i_++)
		for(int j = col*bs, j_ = 0 ; j_ < bs ; j++, j_++)
			L[i*NP+j] = L_sub[i_*bs+j_];
}

void put_U(REAL* U_sub, int row, int col){
	for(int i = row*bs, i_ = 0 ; i_ < bs ; i++, i_++)
		for(int j = col*bs, j_ = 0 ; j_ < bs ; j++, j_++)
			U[i*NP+j] = U_sub[i_*bs+j_];
}


void block_lu(int p){
	REAL* L_ii = (REAL*)calloc(bs*bs,sizeof(REAL));	
	REAL* U_ii = (REAL*)calloc(bs*bs,sizeof(REAL));	
	REAL* L_inv_ii = (REAL*)calloc(bs*bs,sizeof(REAL));
	REAL* U_inv_ii = (REAL*)calloc(bs*bs,sizeof(REAL));

	// L_ii init
	for(int i = 0 ; i < bs ; i++)
		L_ii[i*bs+i] = 1;

	lu(L_ii,U_ii, p);		// compute LU value of A_pp

	put_L(L_ii, p, p);
	put_U(U_ii, p, p);
	if(p==BNAM-1)
		return;
	inverseL(L_ii, L_inv_ii, bs);
	inverseU(U_ii, U_inv_ii, bs);
	
	REAL *UT = (REAL*)malloc(bs*bs*sizeof(REAL));
	transpose(U_inv_ii,UT,bs);

	// figure out L value
#pragma omp parallel
{
#pragma omp for nowait
	for(int i = p + 1; i < BNAM ; i++){
		int row_s_point = i*bs;
		int col_s_point = p*bs;
		//#pragma omp parallel for schedule(auto)
		for(int row = 0 ; row < bs ; row++){
			for(int col = 0 ; col < bs ; col++){
				double sum=0;
				for(int bias = 0 ; bias < bs ; bias++){
					sum += A[(row+row_s_point)*NP+ bias+col_s_point]*UT[col*bs+bias];
				}
				L[(row+row_s_point)*NP+col+col_s_point] = sum;
			}
		}
	}


#pragma omp for
	for(int i = p + 1; i < BNAM ; i++){
		int row_s_point = p*bs;
		int col_s_point = i*bs;
		for(int row = 0 ; row < bs ; row++){
			for(int col = 0 ; col < bs ; col++){
				double sum=0.0;
				for(int bias = 0 ; bias < bs ; bias++){
					sum += L_inv_ii[row*bs+bias]*A[(row_s_point+bias)*NP+col_s_point+col];
				}
				U[(row+row_s_point)*NP+col+col_s_point] = sum;
			}
		}
	}
} /* omp parallel end */
	
	free(UT);	

#pragma omp parallel for
	for(int i = p + 1 ; i < BNAM ; i++){
		for(int j = p + 1 ; j < BNAM ; j++){
			int A_row_s_point = bs*i;
			int A_col_s_point = bs*j;
			int L_row_s_point = bs*i;
			int L_col_s_point = bs*p;
			int U_row_s_point = bs*p;
			int U_col_s_point = bs*j;

			for(int row = 0 ; row < bs ; row++){
				for(int col = 0 ; col < bs ; col++){
					double sum=0.0;
					for(int bias = 0 ; bias < bs ; bias++){
						sum += L[(row+L_row_s_point)*NP+bias+L_col_s_point]*U[(U_row_s_point+bias)*NP+U_col_s_point+col];
					}
					A[(row+A_row_s_point)*NP+col+A_col_s_point] -= sum;
				}
			}
		}
	}
}

void lu(REAL* L_ii, REAL* U_ii, int p){
	int row_s_point=p*bs;
	int col_s_point=p*bs;

	for(int k = 0 ; k < bs ; k++){
		U_ii[k*bs + k] = A[(k+row_s_point)*NP + (k+col_s_point)];
		for(int i = k + 1 ; i < bs ; i++){
			L_ii[i*bs + k] = A[(i+row_s_point)*NP + (k+col_s_point)]/U_ii[k*bs + k];
			U_ii[k*bs + i] = A[(k+row_s_point)*NP + (col_s_point+i)];
		}

#pragma omp parallel for schedule(auto)
		for(int i = k + 1 ; i < bs ; i++){
			for(int j = k + 1 ; j < bs ; j++){
				A[(i+row_s_point)*NP+(j+col_s_point)] = A[(i+row_s_point)*NP + (j+col_s_point)] - L_ii[i*bs+k] * U_ii[k*bs + j];
			}
		}

	} // outter k for loop end
}

void inverseL(REAL *L_ii, REAL *L_, int n){

	double sum;
	for (int k = 0 ; k < n ; k++){
		for(int i = k ; i < n ; i++){
			sum = 0.0;
			for(int j = k ; j < i ; j++){
				sum+=L_ii[i*n+j]*L_[j*n+k];
			}
			L_[i*n+k] = ((i == k ? 1 : 0) - sum);
		}
	}
}

void inverseU(REAL *U_ii, REAL *U_, int n){
	
	double sum;
	for (int k = 0 ; k < n ; k++){
		for(int i = k ; i >= 0 ; i--){
			sum = 0.0;
			for(int j = i+1 ; j <= k ; j++){
				sum+=U_ii[i*n+j]*U_[j*n+k];
			}
			U_[i*n+k] = ((i == k ? 1.0 : 0.0) - sum) / U_ii[i*n+i];
		}
	}
}

// transpose matrix
void transpose(REAL *X, REAL *Y,int n){
	int nbj=16;
#pragma omp parallel for
	for(int j1 = 0 ; j1 < n ; j1+=nbj)
		for(int i = 0 ; i < n ; i++)
			for(int j2 = 0 ; j2 < MIN(n-j1,nbj); j2++)
				Y[i*n+j1+j2] = X[(j1+j2)*n+i];
}

int check(){
	REAL *LU = (REAL*)calloc(NP*NP,sizeof(REAL));
	REAL *UT = (REAL*)malloc(NP*NP*sizeof(REAL));
	int retval = 0;

	transpose(U,UT,NP);
#pragma omp parallel for schedule(auto)
	for(int row = 0 ; row < n ; row++){
		for(int col = 0 ; col < n ; col++){
			double sum=0;
			for(int bias = 0 ; bias < n ; bias++){
					sum += L[row*NP+bias]*UT[col*NP+bias];
			}
			LU[row*NP+col] = sum;
		}
	}
	free(UT);

#pragma omp for schedule(auto) 
	for(int row = 0 ; row < n; row++){
		for(int col = 0 ; col < n ; col++){
			if(compare_REAL(A_save[row*NP+col], LU[row*NP+col])!=0){
				retval = -1;
				break;
			}
		}
	}

	free(LU);
	return retval;
}

int main(int argc, char *argv[]){

	REAL time_in_start, time_in_finish, time_lu_start, time_lu_finish, time_ch_start, time_ch_finish;

	// check and parse command line options
	int c;
	while ((c = getopt(argc, argv, "nd")) != -1) {
		switch (c) {
			case 'd':
				debug_mode = true;
				break;
			case 'n':
				check_mode = true;
				break;
			default:
				printf("Usage: %s [-n] <size> <seed> <thread number> <print>\n", argv[0]);
				exit(EXIT_FAILURE);
		}
	}
	if (optind != argc-4) {
		printf("Usage: %s [-n] <size> <seed> <thread number> <print>\n", argv[0]);
		exit(EXIT_FAILURE);
	}

	// read or generate linear system
	long int size = strtol(argv[optind], NULL, 10);
	rand_seed = strtol(argv[optind+1], NULL, 10);
	int threads_num = strtol(argv[optind+2], NULL, 10);
	d_flag = strtol(argv[optind+3], NULL, 10);
	omp_set_num_threads(threads_num);

	time_in_start  = omp_get_wtime();
	n = (int)size;
	BNAM = n / BS;
	BNAM += (n % BS == 0 ? 0 : 1);

	if (n%BS==0)
		padding=0;
	else
		padding = BS - (n % BS);
	
	NP = n+padding;
	int p_loop_cnt = (n > BNAM ? BNAM : n);
	bs = BS;

	if(debug_mode)
		printf("<bs:%d>>\b<NP:%d>>\b<pad:%d>>\b<BNAM:%d>\n", bs, NP, padding, BNAM);

	matrix_init();

	// perform gaussian elimination
	time_in_finish = time_lu_start = omp_get_wtime();
	for(int p = 0 ; p < p_loop_cnt ; p++)
		block_lu(p);
	time_lu_finish = omp_get_wtime();

	if (d_flag == 1) {
		printf("L:\n");
		print_matrix(L, n, NP);
		printf("U:\n");
		print_matrix(U, n, NP);
		printf("A:\n");
		print_matrix(A_save, n , NP);
	}


	time_ch_start = omp_get_wtime();
	int check_ok;
	if(!check_mode)
		check_ok = 0;
	else
		check_ok = check();
	time_ch_finish = omp_get_wtime();

	if(debug_mode){
		printf("%.5f\t\t%.5f\t\t%.5f\t\t%d\n",
				(float)(time_in_finish-time_in_start),
				(float)(time_lu_finish-time_lu_start),
				(float)(time_ch_finish-time_ch_start),
				check_ok
			  );
	}

	// clean up and exit
	free(A);
	free(A_save);
	free(L);
	free(U);
	return EXIT_SUCCESS;
}
