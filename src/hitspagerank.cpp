//	ISKANDAR ASKAROV 	cs610	8515	prp

//============================================================================
// Name        : hitspagerank.cpp
// Author      : Iskandar Askarov
// Version     :
// Copyright   : Use it all you want!!!
// Description : Hello World in C++, Ansi-style
//============================================================================

#include <iostream>
#include <cstdlib>
#include "udef.h"
#include "hits.h"

#include<string>

using std::stoi;
using std::cout;
using std::endl;
using std::cerr;
using std::ios;

int main(int argc, char **argv) {

	if(argc != 5){
		HELP();
		return EXIT_FAILURE;
	}

	ifstream in;
	in.open(CMD_INPUTFILE, ios::in);

	if(!in.is_open()){
		cerr << "Unable to open input file" << endl;
		return EXIT_FAILURE;
	}

	if(CMD_ALGO == '-h'){
		hits hts(stoi(CMD_ITERATIONS), stoi(CMD_INITVAL));
		int **pmatrix = hts.generate_matrix(in);
		hts.init_hub_auth();
		hts.run(pmatrix);

		for(int i = 0; i < hts.get_vertex_size(); i++){
			delete [] pmatrix[i];
		}
		delete [] pmatrix;
	} else if(CMD_ALGO == '-p'){

	} else{
		HELP();
		return EXIT_FAILURE;
	}



	return EXIT_SUCCESS;
}

//	ISKANDAR ASKAROV 	cs610	8515	prp
