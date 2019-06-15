#include<stdio.h>
#include<unistd.h>
#include<omp.h>

#define NUM 30000

int foo(int a){
	return a + 2;
}

void bar(){
	for(int i = 0 ; i < NUM ; i++){
		int a = foo(i);
		double s = a * 1.000001;
		s *= i;
		s /= 0.482931;
	}
}

int main(){
	double s_time = omp_get_wtime();
	
#pragma omp parallel for num_threads(10)
	for(int i = 0 ; i < NUM ; i++){
		bar();
	}

	printf("%lf\n",omp_get_wtime()-s_time);

	return 0;
}
