#pragma once
#include "gene.h"

class evolver
{
public:

	int num_population;
	int num_elite;
	int innovation_number;
	int num_input, num_output;

	gene** genome;

	evolver();
	~evolver();
	void instantiate(int num_input, int num_output, int population, int elite);
	void print_best();
	void select_best();
	void mutate();
	void evolve(int num_samples, float** input, float** output);
	void evolve_time_series(int num_samples, float** input, float** output, int window_size);
	float* predict_time_series(int num_samples, float** input, int window_size, float* data_point);

};

