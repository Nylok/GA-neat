#include "gene.h"
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include "network.h"

gene::gene()
{
	fitness = 0;

	nodes = nullptr;
	genome = nullptr;
	num_nodes = 0;
	num_connections = 0;
	mse = 0;
}


gene::~gene()
{

}

void gene::instantiate(int num_input, int num_output)
{
	num_nodes = num_input + num_output;
	nodes = new int[(num_input + num_output)];
	for (int i = 0; i < num_input; i++)
		nodes[i] = TYPE_INPUT;
	

	for (int i = num_input; i < (num_input + num_output); i++)
		nodes[i] = TYPE_OUTPUT;

	num_connections = num_input*num_output;
	genome = new connection[num_connections];
	int tmp = 0;

	for (int in = 0; in < num_input; in++)
	{
		for (int out = 0; out < num_output; out++)
		{
			genome[tmp].active = 1;
			genome[tmp].dst = out+num_input;
			genome[tmp].src = in;
			genome[tmp].strenght = 1;
			genome[tmp].innovation_number = tmp;
			tmp++;
		}
	}


};

void gene::add_connection(int src, int dst, int innovation_number)
{
	int tmp_connections = num_connections;
	connection* tmp_genome = genome;

	num_connections = num_connections + 1;
	genome = new connection[num_connections];

	for (int i = 0; i < tmp_connections; i++)
	{
		//coppy old to new

		genome[i].active = tmp_genome[i].active;
		genome[i].dst = tmp_genome[i].dst;
		genome[i].src = tmp_genome[i].src;
		genome[i].strenght = tmp_genome[i].strenght;
		genome[i].innovation_number = tmp_genome[i].innovation_number;
	}

	#pragma warning(suppress: 6386)
	genome[tmp_connections].active = 1;
	genome[tmp_connections].dst = dst;
	genome[tmp_connections].src = src;
	genome[tmp_connections].strenght = 1;
	genome[tmp_connections].innovation_number = innovation_number;

	delete[] tmp_genome;
};

void gene::split_connection(int src, int dst, int innovation_number)
{
	int target_connection;

	//coppy genome and expand by 2 slots
	int tmp_connections = num_connections;
	connection* tmp_genome = genome;

	num_connections = num_connections + 2;
	genome = new connection[num_connections];

	for (int i = 0; i < tmp_connections; i++)
	{
		//coppy old to new

		genome[i].active = tmp_genome[i].active;
		genome[i].dst = tmp_genome[i].dst;
		genome[i].src = tmp_genome[i].src;
		genome[i].strenght = tmp_genome[i].strenght;
		genome[i].innovation_number = tmp_genome[i].innovation_number;
	}

	
	target_connection = 0;
	for (int i = 0; i < tmp_connections; i++)
	{
#pragma warning(suppress: 6385)
		if (genome[i].src == src && genome[i].dst == dst)
		{
			genome[i].active = 0;
			target_connection = i;
		}
	}

	//coppy nodes and expand by one
	int tmp_num_nodes = num_nodes;
	int* tmp_nodes = nodes;


	num_nodes += 1;
	nodes = new int[num_nodes];
	nodes[tmp_num_nodes] = TYPE_HIDDEN;
	for (int i = 0; i < tmp_num_nodes; i++)
#pragma warning(suppress: 6386)
		nodes[i] = tmp_nodes[i];


	genome[tmp_connections].active = 1;
	genome[tmp_connections].dst = tmp_num_nodes;
	genome[tmp_connections].src = src;
	genome[tmp_connections].strenght = 1;
	genome[tmp_connections].innovation_number = innovation_number;

#pragma warning(suppress: 6386)
	genome[tmp_connections+1].active = 1;
	genome[tmp_connections+1].dst = dst;
	genome[tmp_connections+1].src = tmp_num_nodes;
	genome[tmp_connections+1].strenght = genome[target_connection].strenght;
	genome[tmp_connections+1].innovation_number = innovation_number+1;

	delete[] tmp_genome;
	delete[] tmp_nodes;
};


int gene::run_network(int num_samples, int input_number, int output_number,
	float** input, float** output)
{
	network* n = new network(this);
	mse = 0.0f;
	fitness = 0;

	for (int i = 0; i < num_samples; i++)
	{
		n->reset_accumulators();
		n->run(input_number, output_number, input[i], output[i]);
		fitness += n->get_correct();
		mse += n->get_mse();
	}

	delete n;
	return fitness;
}

void gene::run_recursive(int num_samples, int input_number, int output_number,
	float** input, float** output, int window_size)
{
	network* n = new network(this);
	mse = 0.0f;
	fitness = 0;


	for (int i = 0; i < num_samples; i++)
	{
		n->reset_accumulators();
		
		int pre_run_min = ((i - window_size)>0) ? i - window_size:0;
		
		for (int j = pre_run_min; j < i; j++)
		{
			n->run_recursive(input_number, input[j]);
		}

		n->run_recursive(input_number, output_number, input[i], output[i]);
		fitness += n->get_correct();
		mse += n->get_mse();
	}

	delete n;
};

float* gene::run_recursive_predict(int num_samples, int input_number, int output_number,
	float** input, float* to_predict, int window_size)
{
	network* n = new network(this);
	

	n->reset_accumulators();

	int pre_run_min = ((num_samples - window_size)>0) ? num_samples - window_size : 0;

	for (int j = pre_run_min; j < num_samples; j++)
	{
		n->run_recursive(input_number, input[j]);
	}

	float* t = n->run_recursive(input_number, output_number, to_predict);
	
	//needs to be reactivated
	delete n;

	return t;

};

void inherit(int fittness_one, int fittness_two, gene* parent_one,
	gene* parent_two);

void gene::inherit(int fittness_one, int fittness_two, gene* parent_one,
	gene* parent_two)
{
	int tmp_num_nodes;
	int* tmp_nodes;
	int num_new_connection = 0;
	connection* new_connection;



	if (fittness_one == fittness_two)
	{
		tmp_num_nodes = (parent_one->num_nodes > parent_two->num_nodes) ?
			parent_one->num_nodes : parent_two->num_nodes;

		tmp_nodes = new int[tmp_num_nodes];

		if (parent_one->num_nodes > parent_two->num_nodes)
		{
			for (int i = 0; i < parent_one->num_nodes; i++)
			{
				tmp_nodes[i] = parent_one->nodes[i];
			}
		}
		else
		{
			for (int i = 0; i < parent_two->num_nodes; i++)
			{
				tmp_nodes[i] = parent_two->nodes[i];
			}
		}

		int index_one = 0;
		int index_two = 0;

		while (index_one != parent_one->num_connections && index_two != parent_two->num_connections)
		{
			if (parent_one->genome[index_one].innovation_number == parent_two->genome[index_two].innovation_number)
			{
				num_new_connection++;
				index_one++;
				index_two++;
			}
			else if (parent_one->genome[index_one].innovation_number > parent_two->genome[index_two].innovation_number)
			{
				index_two++;
			}
			else if (parent_one->genome[index_one].innovation_number < parent_two->genome[index_two].innovation_number)
			{
				index_one++;
			}
		}

		index_one = 0;
		index_two = 0;

		new_connection = new connection[num_new_connection];

		int pos = 0;
		while (index_one != parent_one->num_connections && index_two != parent_two->num_connections)
		{
			if (parent_one->genome[index_one].innovation_number == parent_two->genome[index_two].innovation_number)
			{
				if (rand() % 2)
				{
					new_connection[pos].active = parent_one->genome[index_one].active;
					new_connection[pos].dst = parent_one->genome[index_one].dst;
					new_connection[pos].src = parent_one->genome[index_one].src;
					new_connection[pos].innovation_number = parent_one->genome[index_one].innovation_number;
					new_connection[pos].strenght = parent_one->genome[index_one].strenght;

				}
				else
				{
					new_connection[pos].active = parent_two->genome[index_two].active;
					new_connection[pos].dst = parent_two->genome[index_two].dst;
					new_connection[pos].src = parent_two->genome[index_two].src;
					new_connection[pos].innovation_number = parent_two->genome[index_two].innovation_number;
					new_connection[pos].strenght = parent_two->genome[index_two].strenght;
				}

				pos++;
				index_one++;
				index_two++;
			}
			else if (parent_one->genome[index_one].innovation_number > parent_two->genome[index_two].innovation_number)
			{
				index_two++;
			}
			else if (parent_one->genome[index_one].innovation_number < parent_two->genome[index_two].innovation_number)
			{
				index_one++;
			}
		}

	}
	else
	{
		gene* fit_min = nullptr;
		gene* fit_max = nullptr;
		if (fittness_one < fittness_two)
		{
			fit_min = parent_one;
			fit_max = parent_two;
		}

		if (fittness_one > fittness_two)
		{
			fit_max = parent_one;
			fit_min = parent_two;
		}

		tmp_num_nodes = fit_max->num_nodes;
		tmp_nodes = new int[tmp_num_nodes];
		for (int i = 0; i < tmp_num_nodes; i++)
		{
			tmp_nodes[i] = fit_max->nodes[i];
		}


		num_new_connection = fit_max->num_connections;
		new_connection = new connection[num_new_connection];

		int index_one = 0;
		int index_two = 0;

		new_connection = new connection[num_new_connection];

		int pos = 0;
		for (int i = 0; i<num_new_connection; i++)
		{
			new_connection[i].active = fit_max->genome[i].active;
			new_connection[i].dst = fit_max->genome[i].dst;
			new_connection[i].src = fit_max->genome[i].src;
			new_connection[i].innovation_number = fit_max->genome[i].innovation_number;
			new_connection[i].strenght = fit_max->genome[i].strenght;

			for (; pos < fit_min->num_connections; pos++)
			{
				if (new_connection[i].dst == fit_min->genome[pos].dst && new_connection[i].src == fit_min->genome[pos].src)
				{
					if (rand() % 2)
					{
						new_connection[i].active = fit_min->genome[pos].active;
						new_connection[i].strenght = fit_min->genome[pos].strenght;
					}
				}

				pos++;
			}
		}


		



		/*int pos = 0;
		while (index_one != fit_max->num_connections)
		{
			if (fit_max->genome[index_one].innovation_number == fit_min->genome[index_two].innovation_number)
			{
				if (rand() % 2)
				{
					new_connection[pos].active = fit_max->genome[index_one].active;
					new_connection[pos].dst = fit_max->genome[index_one].dst;
					new_connection[pos].src = fit_max->genome[index_one].src;
					new_connection[pos].innovation_number = fit_max->genome[index_one].innovation_number;
					new_connection[pos].strenght = fit_max->genome[index_one].strenght;

				}
				else
				{
					new_connection[pos].active = fit_min->genome[index_two].active;
					new_connection[pos].dst = fit_min->genome[index_two].dst;
					new_connection[pos].src = fit_min->genome[index_two].src;
					new_connection[pos].innovation_number = fit_min->genome[index_two].innovation_number;
					new_connection[pos].strenght = fit_min->genome[index_two].strenght;
				}

				pos++;
				index_one++;
				index_two++;
			}
			else if (fit_max->genome[index_one].innovation_number > fit_min->genome[index_two].innovation_number)
			{
				index_two++;
			}
			else if (fit_max->genome[index_one].innovation_number < fit_min->genome[index_two].innovation_number)
			{
				new_connection[pos].active = fit_max->genome[index_one].active;
				new_connection[pos].dst = fit_max->genome[index_one].dst;
				new_connection[pos].src = fit_max->genome[index_one].src;
				new_connection[pos].innovation_number = fit_max->genome[index_one].innovation_number;
				new_connection[pos].strenght = fit_max->genome[index_one].strenght;

				pos++;
				index_one++;
			}
		}
		*/


	}
	


	delete[] genome;
	genome = new_connection;
	delete[] nodes;
	nodes = tmp_nodes;
	num_nodes = tmp_num_nodes;
	num_connections = num_new_connection;
};

void gene::mutate_weights()
{
	for (int i = 0; i < num_connections; i++)
	{
		if (genome[i].active == 0)
			continue;


		int rval = rand() % 1000;

		if (rval < 10)//1% flip sign
		{
			genome[i].strenght = -genome[i].strenght;
		}
		else if (rval < 20)//1% random
		{
			genome[i].strenght = (2.0f*((float)rand() / (float)RAND_MAX)) - 1.0f;
		}
		else if (rval < 30)//1% increase from 100% up to 200%
		{
			float mlt = ((float)rand() / (float)RAND_MAX)+1.0f;
			genome[i].strenght = genome[i].strenght*mlt;
		}
		else if (rval < 40)//1% decrese from 100% up to 0%
		{
			float mlt = ((float)rand() / (float)RAND_MAX);
			genome[i].strenght = genome[i].strenght*mlt;
		}

	}
};

int gene::disparity(gene* to_compare)
{
	int d = 0;
	d = abs(this->num_nodes - to_compare->num_nodes);
	d += abs(this->num_connections - to_compare->num_connections);

	int index_one, index_two;
	index_one = index_two = 0;

	while (index_one != this->num_connections && index_two != to_compare->num_connections)
	{
		if (this->genome[index_one].innovation_number == to_compare->genome[index_two].innovation_number)
		{
			if (this->genome[index_one].strenght != to_compare->genome[index_two].strenght)
			{
				d++;
			}

			index_one++;
			index_two++;
		}
		else if (this->genome[index_one].innovation_number > to_compare->genome[index_two].innovation_number)
		{
			d++;
			index_two++;
		}
		else if (this->genome[index_one].innovation_number < to_compare->genome[index_two].innovation_number)
		{
			d++;
			index_one++;
		}
	}
	return d;
};
