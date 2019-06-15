#include <stdio.h>
#include <fstream>
#include <stdlib.h>
#include <string>
#include <iostream>

using namespace std;

int main(int argc, char **argv){

	string file_name = "in";
	int N;
	string nfile;

	if(argc < 5){
		cout << "usage: <program> <N> <stage> <file#>" << endl;
		return -1;
	}

	N = atoi(argv[1]);
	int stage = atoi(argv[2]);
	int ngosts = atoi(argv[3]);
	nfile = argv[4];
	file_name.append(nfile);
	file_name.append(".txt");

	ofstream output(file_name.c_str());

	int random;

	output << N << endl;
	output << stage << endl;
	output << ngosts << endl;

	for (int i = 0; i < N; i++){
		for (int j = 0 ; j < N ; j++){
			random = rand() % 10;
			output << (random == 1 ? '#' : '.');
		}
		output << endl;
	}



	return 0;
}
