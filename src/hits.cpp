/*
 * hits.cpp
 *
 *  Created on: Apr 23, 2019
 *      Author: iskandar
 */

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


int* hits::matsplit(string line){

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

int hits::get_vertex_size(){
	return this->vertex_count;
}

int **hits::generate_matrix(ifstream &rfp){

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

void hits::init_hub_auth(){
    this->hub0 = new double[this->vertex_count];
    this->auth0 = new double[this->vertex_count];

    for(int i = 0; i < this->vertex_count; i++) {
        switch(this->initialval) {
        case 0:
        	this->hub0[i] = 0;
        	this->auth0[i] = 0;
        	break;
    	case 1:
    		this->hub0[i] = 1;
    		this->auth0[i] = 1;
    		break;
    	case -1:
    		this->hub0[i] = 1.0/this->vertex_count;
    		this->auth0[i] = 1.0/this->vertex_count;
    		break;
    	case -2:
    		this->hub0[i] = 1.0/sqrt(this->vertex_count);
    		this->auth0[i] = 1.0/sqrt(this->vertex_count);
    		break;
        }
    }
}

bool hits::converged(double *p, double *q) const{

    for(int i = 0 ; i < this->vertex_count; i++) {
        if (abs(p[i] - q[i]) > this->errrate)
          return false;
    }
    return true;
}

void hits::run(int **pmatrix){

	double hub[this->vertex_count];
    double auth[this->vertex_count];
    double auth_prev[this->vertex_count];
    double hub_prev[this->vertex_count];
    double a_scale_factor = 0.0;
    double a_sum_square = 0.0;
    double h_scale_factor = 0.0;
    double h_sum_square = 0.0;

    if(this->vertex_count > 10) {
        this->iteration = 0;
        for(int i = 0; i < this->vertex_count; i++) {
            hub[i] = 1.0 / this->vertex_count;
            auth[i] = 1.0 / this->vertex_count;
            hub_prev[i] = hub[i];
            auth_prev[i] = auth[i];
        }

      int iter = 0;

      while(false == this->converged(auth, auth_prev) || false == this->converged(hub, hub_prev)){
          for(int r = 0; r < this->vertex_count; r++) {
              auth_prev[r] = auth[r];
              hub_prev[r] = hub[r];
          }

           for(int p = 0; p < this->vertex_count; p++) { auth[p] = 0.0; }

           for(int j = 0; j < this->vertex_count; j++) {
               for(int k = 0; k < this->vertex_count; k++) {
                   if(pmatrix[k][j] == 1) {  auth[j] += hub[k]; }
               }
           }

           for(int p = 0; p < this->vertex_count; p++) {  hub[p] = 0.0; }

           for(int j = 0; j < this->vertex_count; j++) {
               for(int k = 0; k < this->vertex_count; k++) {
                   if(pmatrix[j][k] == 1) { hub[j] += auth[k]; }
               }
           }

           a_scale_factor = 0.0;
           a_sum_square = 0.0;
           for(int l = 0; l < this->vertex_count; l++) {
               a_sum_square += auth[l]*auth[l];
           }
           a_scale_factor = sqrt(a_sum_square);
           for(int l = 0; l < this->vertex_count; l++) {
               auth[l] = auth[l]/a_scale_factor;
           }

           h_scale_factor = 0.0;
           h_sum_square = 0.0;
           for(int l = 0; l < this->vertex_count; l++) {
               h_sum_square += hub[l]*hub[l];
           }
           h_scale_factor = sqrt(h_sum_square);
           for(int l = 0; l < this->vertex_count; l++) {
               hub[l] = hub[l]/h_scale_factor;
           }
           iter++;
      }

      std::cout << "Iter:    " << iter ;
      std::cout.flush();
      for(int l = 0; l < this->vertex_count; l++) {
    	  double a = auth[l];
    	  double h = hub[l];
    	  printf("A/H[%d]=%.6f/%.6f\n",l, round(a * 1000000.0) / 1000000.0, round(h * 1000000.0) / 1000000.0);
      }
      NEWLINE;
      return;
    }

    for(int i = 0; i < this->vertex_count; i++) {
        hub[i] = this->hub0[i];
        auth[i] = this->auth0[i];
        hub_prev[i] = hub[i];
        auth_prev[i] = auth[i];
    }

    std::cout << "Base:    0 :";
    std::cout.flush();
    for(int i = 0; i < this->vertex_count; i++) {
    	double a = this->auth0[i];
    	double h = this->hub0[i];
    	printf(" A/H[%d]=%.6f/%.6f", i, round(a * 1000000.0)/1000000.0, round(h * 1000000.0)/1000000.0);
    	std::cout.flush();
    }
    NEWLINE;

    if (this->iteration != 0) {
        for(int i = 0; i < this->iteration; i++) {

            for(int p = 0; p < this->vertex_count; p++) {
                auth[p] = 0.0;
            }

            for(int j = 0; j < this->vertex_count; j++) {
                for(int k = 0; k < this->vertex_count; k++) {
                    if(pmatrix[k][j] == 1) { auth[j] += hub[k]; }
                }
            }

            for(int p = 0; p < this->vertex_count; p++) {  hub[p] = 0.0; }

            for(int j = 0; j < this->vertex_count; j++) {
            	for(int k = 0; k < this->vertex_count; k++) {
            		if(pmatrix[j][k] == 1) { hub[j] += auth[k]; }
            	}
            }

            a_scale_factor = 0.0;
            a_sum_square = 0.0;
            for(int l = 0; l < this->vertex_count; l++) {  a_sum_square += auth[l] * auth[l]; }

            a_scale_factor = sqrt(a_sum_square);

            for(int l = 0; l < this->vertex_count; l++) { auth[l] = auth[l] / a_scale_factor; }

            h_scale_factor = 0.0;
            h_sum_square = 0.0;
            for(int l = 0; l < this->vertex_count; l++) { h_sum_square += hub[l] * hub[l]; }
            h_scale_factor = sqrt(h_sum_square);

            for(int l = 0; l < this->vertex_count; l++) { hub[l] = hub[l] / h_scale_factor; }

            std::cout << "Iter:    " << (i+1) <<" :";
            std::cout.flush();
            for(int l = 0; l < this->vertex_count; l++) {
            	double a = auth[l];
            	double h = hub[l];
                printf(" A/H[%d]=%.6f/%.6f",l, round(a * 1000000.0) / 1000000.0,  round(h * 1000000.0)/1000000.0);
                std::cout.flush();
            }
            NEWLINE;
        }
    }
    else
    {
      int iter = 0;
      while(false == this->converged(auth, auth_prev) || false == this->converged(hub, hub_prev)){
          for(int i = 0; i < this->vertex_count; i++) {
              auth_prev[i] = auth[i];
              hub_prev[i] = hub[i];
          }

          //A step starts
          for(int i = 0; i < this->vertex_count; i++) {
          	auth[i] = 0.0;
          }

          for(int i = 0; i < this->vertex_count; i++) {
              for(int k = 0; k < this->vertex_count; k++) {
                  if(pmatrix[k][i] == 1) {
                  	auth[i] += hub[k];
                  }
              }
          }

          for(int i = 0; i < this->vertex_count; i++) {
          	hub[i] = 0.0;
          }

          for(int j = 0; j < this->vertex_count; j++) {
              for(int k = 0; k < this->vertex_count; k++) {
                  if(pmatrix[j][k] == 1) {
                      hub[j] += auth[k];
                  }
              }
          }

          a_scale_factor = 0.0;
          a_sum_square = 0.0;
          for(int l = 0; l < this->vertex_count; l++) {
              a_sum_square += auth[l] * auth[l];
          }
          a_scale_factor = sqrt(a_sum_square);
          for(int l = 0; l < this->vertex_count; l++) {
          	auth[l] = auth[l]/a_scale_factor;
          }

          h_scale_factor = 0.0;
          h_sum_square = 0.0;
          for(int l = 0; l < this->vertex_count; l++) {
              h_sum_square += hub[l] * hub[l];
          }
          h_scale_factor = sqrt(h_sum_square);
          for(int l = 0; l < this->vertex_count; l++) {
          	hub[l] = hub[l]/h_scale_factor;
          }
          iter++;

          std::cout << "Iter:    " << iter << " :";
          std::cout.flush();
          for(int l = 0; l < this->vertex_count; l++) {
        	  double a = auth[l];
        	  double h = hub[l];
        	  printf(" A/H[%d]=%.6f/%.6f",l, round(a * 1000000.0) / 1000000.0, round(h * 1000000.0) / 1000000.0);
        	  std::cout.flush();
          }
          NEWLINE;
      }
    }
}

hits::~hits(){
	delete [] this->auth0;
	delete [] this->hub0;
}
