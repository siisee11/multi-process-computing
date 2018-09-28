default: lu

lu: lu.c
	gcc -g -Wall -o lu lu.c -fopenmp -O2

par_lu: par_lu.c
	gcc -g -O2 -fopenmp -Wall -o par_lu par_lu.c

clean:
	rm -f lu par_lu 
