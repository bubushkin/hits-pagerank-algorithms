/*
 * hits.h
 *
 *  Created on: Apr 23, 2019
 *      Author: iskandar
 */

#ifndef HITS_H_
#define HITS_H_

#include <fstream>
#include <string>
#include <vector>

using std::ifstream;
using std::string;
using std::vector;


class hits {

private:

	int iteration;
	int initialval;
	int vertex_count;
	int edge_count;
	double *hub0;
	double *auth0;

	const double errrate = 0.00001;



public:
	//inline constructor

	hits(int iteration, int initval): iteration(iteration), initialval(initval), vertex_count(0), edge_count(0){}

	int** generate_matrix(ifstream &rfp);

	int* matsplit(string line);

	void init_hub_auth();

	bool converged(double *p, double *q) const;

	int get_vertex_size();

	void run(int **pmatrix);

	~hits();



};

#endif /* HITS_H_ */
