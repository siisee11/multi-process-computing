#include <omp.h>
#include <stdio.h>
#include <unistd.h>

int main(){
#pragma omp parallel sections
	{
#pragma omp section
		{
			printf("first section\n");
			usleep(1);
			printf("first section\n");
#pragma omp parallel sections
			{
#pragma omp section
				{
					printf("first section - nested 1\n");
					usleep(1);
					printf("first section - nested 1\n");
				}
#pragma omp section
				{
					printf("first section - nested 2\n");
					usleep(1);
					printf("first section - nested 2\n");
				}
			}
		}
#pragma omp section
		{
			printf("second section\n");
			usleep(1);
			printf("second section\n");
		}
	}
	return 0;
}
