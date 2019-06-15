#include <stdio.h>

int main(){

	int sum = 0;

#pragma omp parallel reduction(+:sum)
{
	int i = 0;
	#pragma omp parallel for
	for( i = 0 ; i < 1000 ; i++){
		sum += i;
	}
}


