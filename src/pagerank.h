/*
 * pagerank.h
 *
 *  Created on: Apr 23, 2019
 *      Author: iskandar
 */

#ifndef PAGERANK_H_
#define PAGERANK_H_

#include <fstream>
#include <string>

using std::ifstream;
using std::string;

class pagerank {

private:
	int iteration;
	int initialval;
	int vertex_count;
	int edge_count;

	int *pOutgoing;
	double *pSrc;

	const double errrate = 0.00001;
	const double d = 0.85;

public:

	pagerank(int iteration, int initval): iteration(iteration), initialval(initval), vertex_count(0), edge_count(0){}

	~pagerank();

	int** generate_matrix(ifstream &rfp);

	int* matsplit(string line);

	bool converged(double *p, double *q) const;

	void init_src_out(int **pmatrix);

	void run_pagerank(int **pmatrix);

	int get_vertex_size();



};

#endif /* PAGERANK_H_ */
