
#include <iostream>
#include <cstdlib>
#include "udef.h"
#include "hits.h"
#include "pagerank.h"

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

	string algo = CMD_ALGO;
	if(algo == "-h"){
		hits hts(stoi(CMD_ITERATIONS), stoi(CMD_INITVAL));
		int **pmatrix = hts.generate_matrix(in);
		hts.init_hub_auth();
		hts.run_hits(pmatrix);

		for(int i = 0; i < hts.get_vertex_size(); i++){
			delete [] pmatrix[i];
		}
		delete [] pmatrix;
	} else if(algo == "-p"){

		pagerank rank(stoi(CMD_ITERATIONS), stoi(CMD_INITVAL));
		int **pmatrix = rank.generate_matrix(in);
		rank.init_src_out(pmatrix);
		rank.run_pagerank(pmatrix);
		for(int i = 0; i < rank.get_vertex_size(); i++){
			delete [] pmatrix[i];
		}
		delete [] pmatrix;

	} else{
		HELP();
		return EXIT_FAILURE;
	}



	return EXIT_SUCCESS;
}
