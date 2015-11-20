#pragma once
#include "gene.h"

class network
{
	float mse;
	int correct;
	int nodes_input, nodes_output, nodes_total;

	float* vals;
	float* vals_prev;
	int* updated;
	int* update_order;


	float correct_margin;

	int* connection_number;

	int** connection_point;
	float** connection_weight;

public:
	network();
	network(gene* src);
	//network(gene* src, float margin);
	~network();

	void run(int num_input, int num_output, float* input, float* output);
	int get_correct();
	float get_mse();
	void reset_accumulators();
	void run_recursive(int num_input, float* input);
	void run_recursive(int num_input, int num_output, float* input, float* output);
	float* run_recursive(int num_input, int num_output, float* input);
};

