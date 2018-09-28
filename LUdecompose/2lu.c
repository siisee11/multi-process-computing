#include <getopt.h>
#include <limits.h>
#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <omp.h>

#define SWAP(a, b, type) do { \
	type temp; \
	temp = a;  \
	a = b;     \
	b = temp;  \
} while (0)

// use 64-bit IEEE arithmetic (change to "float" to use 32-bit arithmetic)
#define REAL double

// linear system: PA = LU   (A is n x n matrix; b and x are n x 1 vectors)
int n;
REAL *A;
REAL *P;
REAL *L;
REAL *U;

// enable/disable debugging output (don't enable for large matrix sizes!)
bool debug_mode = false;

// enable/disable triangular mode (to skip the Gaussian elimination phase)
bool triangular_mode = false;

/*
 * Generate a random linear system of size n.
 */
void rand_system()
{
	// allocate space for matrices
	A = (REAL*)calloc(n*n, sizeof(REAL));

	// initialize pseudorandom number generator
	// (see https://en.wikipedia.org/wiki/Linear_congruential_generator)
	unsigned long seed = 0;

	// generate random matrix entries

//#pragma omp parallel for
	for (int row = 0; row < n; row++) {
		int col = triangular_mode ? row : 0;
		for (; col < n; col++) {
			if (row != col) {
				seed = (1103515245*seed + 12345) % (1<<31);
				A[row*n + col] = (REAL)seed / (REAL)ULONG_MAX;
			} else {
				A[row*n + col] = n/10.0;
			}
		}
	}
}

/*
 * Reads a linear system of equations from a file in the form of an augmented
 * matrix [A][b].
 */
void read_system(const char *fn)
{
	// open file and read matrix dimensions
	FILE* fin = fopen(fn, "r");
	if (fin == NULL) {
		printf("Unable to open file \"%s\"\n", fn);
		exit(EXIT_FAILURE);
	}
	if (fscanf(fin, "%d\n", &n) != 1) {
		printf("Invalid matrix file format\n");
		exit(EXIT_FAILURE);
	}

	// allocate space for matrices
	A = (REAL*)malloc(sizeof(REAL) * n*n);

	// read all values
	for (int row = 0; row < n; row++) {
		for (int col = 0; col < n; col++) {
			if (fscanf(fin, "%lf", &A[row*n + col]) != 1) {
				printf("Invalid matrix file format\n");
				exit(EXIT_FAILURE);
			}
		}
	}
	fclose(fin);
}

/*
 * Prints a matrix to standard output in a fixed-width format.
 */
void print_matrix(REAL *mat, int rows, int cols)
{
	for (int row = 0; row < rows; row++) {
		for (int col = 0; col < cols; col++) {
			printf("%8.1e ", mat[row*cols + col]);
		}
		printf("\n");
	}
}

void matrix_init(){
	L = (REAL*)calloc(n*n, sizeof(REAL));
	U = (REAL*)calloc(n*n, sizeof(REAL));
	P = (REAL*)calloc(n, sizeof(REAL));

	for(int i = 0; i < n ; i++){
		L[i*n + i] = 1;
		P[i] = i;
	}
}

void lu(){
	int tmp;
	int k_;
	int max;

	for(int k = 0 ; k < n ; k++){
		max = 0;
		for(int i = k ; i < n ; i++){
			tmp = A[i*n+k];
			tmp = tmp > 0 ? tmp : -tmp;
			if(max < tmp){
				max = tmp;	
			}
		}
		if (max == 0)
			exit(EXIT_FAILURE);
		SWAP(P[k],P[k_],double);
		for(int j = 1; j <= n ; j++){ 
			SWAP(A[k][j],A[k_][j],double);
		}
		for(int j = 1; j <= k-1 ; j++){
			SWAP(L[k][j],L[k_][j],double);
		}
		U[k][k] = A[k][k];
		for(int i = k+1 ; i <=n ; i++){
			L[i][k] = A[i][k]/U[k][k];
			U[k][i] = A[k][i];
		}
		for(int i = k+1 ; i <= n ; i++){
			for(int j = k+1 ; j <=n ; j++){
				A[i][j] = A[i][j] - L[i][k]*U[k][j];
			}
		}
	}
}

/*
 * Performs Gaussian elimination on the linear system.
 * Assumes the matrix is singular and doesn't require any pivoting.
 */
int main(int argc, char *argv[]){

	REAL time_in_start, time_in_finish, time_ge_start, time_ge_finish;

	// check and parse command line options
	int c;
	while ((c = getopt(argc, argv, "dt")) != -1) {
		switch (c) {
			case 'd':
				debug_mode = true;
				break;
			case 't':
				triangular_mode = true;
				break;
			default:
				printf("Usage: %s [-dt] <file|size>\n", argv[0]);
				exit(EXIT_FAILURE);
		}
	}
	if (optind != argc-1) {
		printf("Usage: %s [-dt] <file|size>\n", argv[0]);
		exit(EXIT_FAILURE);
	}

	// read or generate linear system
	long int size = strtol(argv[optind], NULL, 10);
	time_in_start  = omp_get_wtime();
	if (size == 0) {
		read_system(argv[optind]);
	} else {
		n = (int)size;
		rand_system();
	}

	matrix_init();

	// perform gaussian elimination
	time_in_finish = time_ge_start = omp_get_wtime();
	lu();
	time_ge_finish = omp_get_wtime();

	if (debug_mode) {
		printf("A = \n");
		print_matrix(A, n, n);
		printf("P = \n");
		print_matrix(P, n, 1);
		printf("L = \n");
		print_matrix(L, n, n);
		printf("U = \n");
		print_matrix(U, n, n);
	}

	// print results
	printf("Ntheads=%2d  INIT: %8.4fs  GAUS: %8.4fs\n",
			1, 
			(float)(time_in_finish-time_in_start),
			(float)(time_ge_finish-time_ge_start)
		  );

	// clean up and exit
	free(A);
	free(L);
	free(U);
	free(P);
	return EXIT_SUCCESS;
}
