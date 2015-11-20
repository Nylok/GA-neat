#pragma once

#define TYPE_INPUT 0
#define TYPE_OUTPUT 1
#define TYPE_HIDDEN 2

typedef struct
{
	int innovation_number;
	int active;

	int src;
	int dst;
	
	float strenght;

} connection;

class gene
{
public:

	int fitness;
	float mse;

	int* nodes;
	connection* genome;
	int num_nodes;
	int num_connections;

	gene();
	~gene();

	void instantiate(int num_input, int num_output);
	
	void add_connection(int src, int dst, int innovation_number);
	
	void split_connection(int src, int dst, int innovation_number);	
	
	int run_network(int num_samples, int input_number, int output_number,
		float** input, float** output);

	void run_recursive(int num_samples, int input_number, int output_number,
		float** input, float** output, int window_size);

	float* run_recursive_predict(int num_samples, int input_number, int output_number,
		float** input, float* to_predict, int window_size);

	void inherit(int fittness_one, int fittness_two, gene* parent_one,
		gene* parent_two);

	void mutate_weights();

	int disparity(gene* to_compare);

};

