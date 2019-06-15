#include <stdio.h>

int arr[4];

void insert(int x){
	for (int i = 0 ; i < 4 ; i++){
		if(arr[i] == 0){
			arr[i] = x;
			break;
		}
	}
}

void print_arr(){
	for (int i = 0 ; i< 4 ; i++){
		printf("%d\t", arr[i]);
	}
	printf("\n");
}

int main(){
	
#pragma omp parallel sections
	{
#pragma omp section
		insert(10);

#pragma omp section
		insert(20);
	}

	print_arr();
}
