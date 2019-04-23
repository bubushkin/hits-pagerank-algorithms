/*
 * pagerank.cpp
 *
 *  Created on: Apr 23, 2019
 *      Author: iskandar
 */

#include "pagerank.h"

#include <iostream>
#include <string>
#include <vector>
#include <string.h>
#include <cstdlib>
#include <math.h>
#include "hits.h"
#include "udef.h"

using std::cerr;
using std::endl;
using std::stoi;


int* pagerank::matsplit(string line){

	char *temp;
	int *pval = new int[2];
	char cline[line.size() + 1];
	line.copy(cline, line.size() + 1);

	temp = strtok(cline, " ");
	int i = 0;
	while(temp != NULL && i < 2){
		pval[i++] = stoi(temp);
		temp = strtok(NULL, " ");
	}

	return pval;

}

bool pagerank::converged(double *p, double *q) const{

    for(int i = 0 ; i < this->vertex_count; i++) {
        if (abs(p[i] - q[i]) > this->errrate)
          return false;
    }
    return true;
}


int **pagerank::generate_matrix(ifstream &rfp){

	string line;
	int *limits;
	getline(rfp, line);

	limits = this->matsplit(line);
	this->vertex_count = limits[0];
	this->edge_count = limits[1];
	delete [] limits;

	int **pmatrix = new int*[this->vertex_count];
	for(int i = 0; i < this->vertex_count; i++){
		pmatrix[i] = new int[this->vertex_count];
		for(int j = 0; j < this->vertex_count; j++){
			pmatrix[i][j] = 0;
		}
	}
	while(getline(rfp, line)){
		limits = this->matsplit(line);
		pmatrix[limits[0]][limits[1]] = 1;
		delete [] limits;

	}

	return pmatrix;
}


int pagerank::get_vertex_size(){
	return this->vertex_count;
}


void pagerank::init_src_out(int **pmatrix){

	this->pOutgoing = new int[this->vertex_count];
    this->pSrc = new double[this->vertex_count];

    for(int i = 0; i < this->vertex_count; i++) {
    	this->pOutgoing[i] = 0;
        for(int j = 0; j < this->vertex_count; j++) {
        	this->pOutgoing[i] += pmatrix[i][j];
        }
    }

    for(int i = 0; i < this->vertex_count; i++) {
        switch(this->initialval) {
        case 0:
        	this->pSrc[i] = 0;
        	break;
    	case 1:
    		this->pSrc[i] = 1;
    		break;
    	case -1:
    		this->pSrc[i] = 1.0/this->vertex_count;
    		break;
    	case -2:
    		this->pSrc[i] = 1.0/sqrt(this->vertex_count);
    		break;
        }
    }
}

void pagerank::run_pagerank(int **pmatrix){


	double D[this->vertex_count];

	bool init = true;
	if(this->vertex_count > 10) {

		this->iteration = 0;

		for(int i =0; i < this->vertex_count; i++) { this->pSrc[i] = 1.0 / this->vertex_count; }

		int iter = 0;
		while (!this->converged(this->pSrc, D)){
			if(init)  init = false;
			else  for(int l = 0; l < this->vertex_count; l++) { this->pSrc[l] = D[l]; }

			for(int l = 0; l < this->vertex_count; l++) { D[l] = 0.0; }

			for(int j = 0; j < this->vertex_count; j++) {
				for(int k = 0; k < this->vertex_count; k++){
					if(pmatrix[k][j] == 1) {
						D[j] += this->pSrc[k]/this->pOutgoing[k];
					}
				}
			}
			for(int l = 0; l < this->vertex_count; l++) {
				D[l] = (d * D[l]) + ((1 - d) / this->vertex_count);
			}
			iter++;
		}

		std::cout << "Iter: " << iter;
		std::cout.flush();
		for(int l = 0 ; l < this->vertex_count; l++) {
			double d = D[l];
			printf("P[%d] = %.6f\n", l, round(d * 1000000.0) / 1000000.0);
			std::cout.flush();
		}
		NEWLINE
		return;
	}

	std::cout << "Base    : 0";
	std::cout.flush();
	for(int j = 0; j < this->vertex_count; j++) {
		double src = this->pSrc[j];
		printf(" :P[%d]=%.6f", j, round(src * 1000000.0) / 1000000.0);
		std::cout.flush();
	}


	if (this->iteration != 0) {
		for(int i = 0; i < this->iteration; i++){
			for(int l = 0; l < this->vertex_count; l++) {
				D[l] = 0.0;
			}
			for(int j = 0; j < this->vertex_count; j++) {
				for(int k = 0; k < this->vertex_count; k++){
					if(pmatrix[k][j] == 1) {
						D[j] += this->pSrc[k]/this->pOutgoing[k];
					}
				}
			}
			NEWLINE
			std::cout << "Iter    : " << (i+1);
			std::cout.flush();
			for(int l = 0; l < this->vertex_count; l++) {
				D[l] = (d * D[l]) + ((1-d) / this->vertex_count);
				printf(" :P[%d]=%.6f", l, round(D[l] * 1000000.0) / 1000000.0);

			}
			for(int l = 0; l < this->vertex_count; l++) {
				this->pSrc[l] = D[l];
			}
		}
	}
	else {
		int i = 0;
		while (!this->converged(this->pSrc, D)){
			if(init) init = false;
			else for(int l = 0; l < this->vertex_count; l++) {  this->pSrc[l] = D[l]; }

			for(int l = 0; l < this->vertex_count; l++) {
				D[l] = 0.0;
			}
			for(int j = 0; j < this->vertex_count; j++) {
				for(int k = 0; k < this->vertex_count; k++)
				{
					if(pmatrix[k][j] == 1) {
						D[j] += this->pSrc[k]/this->pOutgoing[k];
					}
				}
			}

			std::cout << "Iter    : " << (i+1);
			std::cout.flush();
			for(int l = 0; l < this->vertex_count; l++) {
				D[l] = (d * D[l]) + ((1-d)/this->vertex_count);
				printf(" :P[%d]=%.6f", l, round(D[l] * 1000000.0) / 1000000.0);
				std::cout.flush();
			}
			i++;
		}
		NEWLINE
	}
}

pagerank::~pagerank() {
	delete [] this->pSrc;
	delete [] this->pOutgoing;
}

