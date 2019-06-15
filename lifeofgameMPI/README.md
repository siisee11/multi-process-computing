# Conway's Life of game (MPI version)

Life of game on MPI

## Directory

All source code are in src directory.
Backup directory is used for backup some files, so you can ignore it.

## Compiling
'''
src\$> mpic++ life.cpp -o life
'''

## Running

(restriction: size of the input grid should be dividable by the number of processors)
example for give input (N=100)

'''
&> mpirun -n 10 --machinefile hosts --map-by node ./life input1.txt
'''

If program terminated with no errors, the output is written into output.txt file.
