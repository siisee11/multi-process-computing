#include <stdio.h>
#include <omp.h>

int main(){
	printf("proc::%d\n",omp_get_num_procs());

	return 0;
}

