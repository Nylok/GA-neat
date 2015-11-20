#include "network.h"
#include <math.h>

//
#define activation(X) 1.0/(1+exp(-X))

//step
//#define activation(X) (X>1.0)?1.0:((X<-1.0)?-1.0:X);

//experimental
//#define activation(X) 1.0/(0.5+exp(-X))

network::network()
{
	mse=0.0;
	correct=0;
	nodes_input = 0;
	nodes_output = 0;
	nodes_total = 0;

	vals = nullptr;
	vals_prev = nullptr;
	updated = nullptr;

	correct_margin = 0.0f;

	connection_number=nullptr;
	connection_point=nullptr;
};

network::network(gene* src)
{
	mse = 0.0;
	correct = 0;

	correct_margin = 0.5f;

	nodes_total = src->num_nodes;
	nodes_input = 0;
	nodes_output = 0;

	for (int i = 0; i < nodes_total; i++)
	{
		if (src->nodes[i] == TYPE_INPUT)
			nodes_input++;

		if (src->nodes[i] == TYPE_OUTPUT)
			nodes_output++;
	}

	vals = new float[nodes_total];
	vals_prev = new float[nodes_total];
	updated = new int[nodes_total];
	update_order = new int[nodes_total];


	connection_number = new int[nodes_total];

	for (int i = 0; i < nodes_total; i++)
	{
		update_order[i] = i;

		connection_number[i] = 0;
		for (int j = 0; j < src->num_connections; j++)
		{
			if (i == src->genome[j].dst && src->genome[j].active == 1)
				connection_number[i]++;
		}
	}

	//swap outputs nodes for last position
	for (int outs = 0; outs < nodes_output; outs++)
	{
		int tmp = update_order[nodes_input + outs];
		update_order[nodes_input + outs] = update_order[nodes_total - outs - 1];
		update_order[nodes_total - outs - 1] = tmp;
	}

	int max = 0;
	for (int i = 0; i < nodes_total; i++)
	{
		max = (max>connection_number[i]) ? max : connection_number[i];
	}

	connection_point = new int*[nodes_total];
	connection_weight = new float*[nodes_total];
	for (int i = 0; i < nodes_total; i++)
	{
		connection_point[i] = new int[max];
		connection_weight[i] = new float[max];

		int tmp = 0;
		for (int j = 0; j < src->num_connections; j++)
		{
			if (i == src->genome[j].dst && src->genome[j].active == 1)
			{
				connection_point[i][tmp] = src->genome[j].src;
				connection_weight[i][tmp] = src->genome[j].strenght;
				tmp++;
			}
		}
	}

};

network::~network()
{
	delete[] vals;
	delete[] vals_prev;
	delete[] updated;
	delete[] update_order;
	
	delete[] connection_number;

	

	for (int i = 0; i < nodes_total; i++)
	{
		if (connection_point[i] != nullptr)
			delete[] connection_point[i];

		if (connection_weight[i] != nullptr)
			delete[] connection_weight[i];
		
	}

	delete[] connection_point;
	delete[] connection_weight;

};

void network::run(int num_input, int num_output, float* input, float* output)
{
	if (num_input != nodes_input || num_output != nodes_output)
		return;

	for (int i = 0; i < nodes_total; i++)
	{
		updated[i] = 0;
	}

	for (int i = 0; i < nodes_input; i++)
	{
		updated[i] = 1;
		vals[i] = (float)input[i];
	}

	for (int i = nodes_input; i < nodes_total; i++)
	{
		int pos = update_order[i];
		for (int j = 0; j < connection_number[pos]; j++)
		{
#define activation(X) 1.0/(1+exp(-X))

			if (updated[connection_point[pos][j]]==0)
			{
				vals[i] += connection_weight[pos][j] * activation(vals_prev[connection_point[pos][j]]);
			}
			else
			{
				vals[i] += vals[i] += connection_weight[pos][j] * activation(vals[connection_point[pos][j]]);
			}

		}

		updated[i] = 1;
	}

	//check of cprrect
	correct = 1;
	mse = 0.0;
	for (int i = 0; i < nodes_output; i++)
	{
		float delta = abs(output[i] - activation(vals[nodes_total - i - 1]));

		if (delta > correct_margin)
			correct = 0;

		mse += delta*delta;

	}
}

void network::reset_accumulators()
{
	for (int i = 0; i < nodes_total; i++)
	{
		vals[i] = 0.0f;
		vals_prev[i] = 0.0f;
		updated[i] = 0;
	}
}


int network::get_correct(){
	return correct;
};

float network::get_mse()
{
	return mse;
};

void network::run_recursive(int num_input, float* input)
{
	if (num_input != nodes_input)
		return;


	float* tmp = vals;
	vals = vals_prev;
	vals_prev = tmp;

	for (int i = 0; i < nodes_total; i++)
		vals[i] = 0.0;

	for (int i = 0; i < nodes_total; i++)
	{
		updated[i] = 0;
	}

	for (int i = 0; i < nodes_input; i++)
	{
		updated[i] = 1;
		vals[i] = (float)input[i];
	}

	for (int i = nodes_input; i < nodes_total; i++)
	{
		int pos = update_order[i];
		for (int j = 0; j < connection_number[pos]; j++)
		{

			if (updated[connection_point[pos][j]] == 0)
			{
				vals[i] += connection_weight[pos][j] * vals_prev[connection_point[pos][j]];
			}
			else
			{
				vals[i] += vals[i] += connection_weight[pos][j] * vals[connection_point[pos][j]];
			}

		}
		vals[i] = activation(vals[i]);
		updated[i] = 1;
	}
};

void network::run_recursive(int num_input, int num_output, float* input, float* output)
{
	if (num_input != nodes_input || num_output != nodes_output)
		return;


	float* tmp = vals;
	vals = vals_prev;
	vals_prev = tmp;

	for (int i = 0; i < nodes_total; i++)
		vals[i] = 0.0;

	for (int i = 0; i < nodes_total; i++)
	{
		updated[i] = 0;
	}

	for (int i = 0; i < nodes_input; i++)
	{
		updated[i] = 1;
		vals[i] = (float)input[i];
	}

	for (int i = nodes_input; i < nodes_total; i++)
	{
		int pos = update_order[i];
		for (int j = 0; j < connection_number[pos]; j++)
		{
			if (updated[connection_point[pos][j]] == 0)
			{
				vals[i] += connection_weight[pos][j] * vals_prev[connection_point[pos][j]];
			}
			else
			{
				vals[i] += vals[i] += connection_weight[pos][j] * vals[connection_point[pos][j]];
			}

		}
		vals[i] = activation(vals[i]);
		updated[i] = 1;
	}

	//check of cprrect
	correct = 1;
	mse = 0.0;
	for (int i = 0; i < nodes_output; i++)
	{
		float delta = abs(output[i] - (activation(vals[nodes_total - i - 1])));

		if (delta > correct_margin)
			correct = 0;

		mse += delta*delta;

	}
};

float* network::run_recursive(int num_input, int num_output, float* input)
{
	if (num_input != nodes_input || num_output != nodes_output)
		return nullptr;


	float* tmp = vals;
	vals = vals_prev;
	vals_prev = tmp;

	float* res = new float[nodes_output];

	for (int i = 0; i < nodes_total; i++)
		vals[i] = 0.0;

	for (int i = 0; i < nodes_total; i++)
	{
		updated[i] = 0;
	}

	for (int i = 0; i < nodes_input; i++)
	{
		updated[i] = 1;
		vals[i] = (float)input[i];
	}

	for (int i = nodes_input; i < nodes_total; i++)
	{
		int pos = update_order[i];
		for (int j = 0; j < connection_number[pos]; j++)
		{

			if (updated[connection_point[pos][j]] == 0)
			{
				vals[i] += connection_weight[pos][j] * vals_prev[connection_point[pos][j]];
			}
			else
			{
				vals[i] += vals[i] += connection_weight[pos][j] * vals[connection_point[pos][j]];
			}

		}
		vals[i] = activation(vals[i]);
		updated[i] = 1;
	}

	//check of cprrect
	for (int i = 0; i < nodes_output; i++)
	{
		res[i]= (activation(vals[nodes_total - i - 1]));
		if (res[i]>0.5)
			res[i] = 1.0;
		else
			res[i] = 0.0;
	}

	return res;
};
