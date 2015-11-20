#include "evolver.h"
#include <math.h>
#include <stdlib.h>
#include <stdio.h>

evolver::evolver()
{

};

evolver::~evolver()
{

};

void evolver::instantiate(int num_input, int num_output, int population, int elite)
{
	num_population = population;
	num_elite = elite;
	innovation_number = num_input*num_output;
	this->num_input = num_input;
	this->num_output = num_output;
	genome = new gene*[num_population];
	for (int i = 0; i < num_population; i++)
	{
		genome[i] = new gene;
		genome[i]->instantiate(num_input, num_output);
	}

};

void evolver::print_best()
{
	printf("fittness %d  mse %.2f  num nodes %d  num connects %d \n", genome[0]->fitness, genome[0]->mse, genome[0]->num_nodes, genome[0]->num_connections);
};

void evolver::select_best()
{
	//dst src
	int** disparity_map;

	//buid disparity map
	disparity_map = new int*[num_population];
	for (int i = 0; i < num_population; i++)
	{	
		disparity_map[i] = new int[num_population];
	}

	for (int d = 0; d < num_population; d++)
	{
		for (int s = 0; s < num_population; s++)
		{
			if (d == s)
			{
				disparity_map[d][s] = 0;
			}
			else
			{
				disparity_map[d][s] = genome[d]->disparity(genome[s]);
			}
		}
	}


	//search for the elite with mare difrences
	int* elite_candidate = new int[num_elite];

	for (int i = 0; i < num_elite; i++)
	{
		elite_candidate[i] = -1;
	}

	int bfit;
	float bmse;
	bfit = genome[0]->fitness;
	bmse = genome[0]->mse;
	for (int i = 1; i < num_population; i++)
	{
		if (genome[i]->fitness < bfit)
		{
			bfit = genome[i]->fitness;
			bmse = genome[i]->mse;
			gene* tmp = genome[0];
			genome[0] = genome[i];
			genome[i] = tmp;
		}
		else if (genome[i]->fitness == bfit && genome[i]->mse <= bmse)
		{
			bfit = genome[i]->fitness;
			bmse = genome[i]->mse;
			gene* tmp = genome[0];
			genome[0] = genome[i];
			genome[i] = tmp;
		}
	}

	elite_candidate[0] = 0;


	for (int i = 1; i < num_elite; i++)
	{
		int candidate_sum = 0;
		int candidate_val = -1;
		for (int candidate = 0; candidate < num_population; candidate++)
		{
			int flag = 0;
			for (int t = 0; t < num_elite; t++)
			{
				if (elite_candidate[t] == candidate)
				{
					flag = 1;
				}
			}

			if (flag == 1)
				continue;

			elite_candidate[i] = candidate;

			int tsum = 0;
			for (int j = 0; j < num_population; j++)
			{
				flag = 0;
				for (int t = 0; t < num_elite; t++)
				{
					if (elite_candidate[t] == j)
					{
						flag = 1;
					}
				}

				if (flag == 1)
					continue;

				tsum += disparity_map[candidate][j];
			}

			if (tsum >= candidate_sum)
			{
				candidate_sum = tsum;
				candidate_val = candidate;
			}
		}

		elite_candidate[i] = candidate_val;
	}

	for (int i = 0; i < num_elite; i++)
	{
		gene* tmp = genome[elite_candidate[i]];
		genome[elite_candidate[i]] = genome[i];
		genome[i] = tmp;
	}

	//cleanup
	for (int i = 0; i < num_population; i++)
	{
		delete[] disparity_map[i];
	}
	delete[] disparity_map;
}

void evolver::mutate()
{
	for (int j = num_elite; j < num_population; j++)
	{
		int ra, rb;
		ra = rand() % num_elite;
		rb = rand() % num_elite;
		genome[j]->inherit(genome[ra]->fitness, genome[rb]->fitness, genome[ra], genome[rb]);

		genome[j]->mutate_weights();
		//falta spit e connection

		if ((rand() % 1000) < 4) //0.1% split
		{
			int give_up_counter = 0;
			//find active connection
			int num = rand() % genome[j]->num_connections;
			while (genome[j]->genome[num].active == 0 && give_up_counter < 100)
			{
				num = rand() % genome[j]->num_connections;
				give_up_counter++;
			}

			if (give_up_counter < 100)
			{
				genome[j]->split_connection(genome[j]->genome[num].src, genome[j]->genome[num].dst, innovation_number);
				innovation_number += 2;
			}
		}
		else if ((rand() % 1000) < 100) //0.1% split
		{
			//find active connection
			int node_s, node_d;

			int flag = 1;
			int give_up_counter = 0;
			while (flag == 1 && give_up_counter<100)
			{
				node_s = rand() % genome[j]->num_nodes;
				node_d = rand() % genome[j]->num_nodes;

				flag = 0;
				for (int k = 0; k < genome[j]->num_connections; k++)
				{
					if (genome[j]->genome[k].src == node_s && genome[j]->genome[k].dst == node_d)
						flag = 1;
				}

				give_up_counter++;
			}
			if (give_up_counter < 100)
			{
				genome[j]->add_connection(node_s, node_d, innovation_number);
				innovation_number++;
			}
		}
	}

}

void evolver::evolve(int num_samples, float** input, float** output)
{
	//mutate
	select_best();
	
	mutate();

	//run
#pragma omp parallel for
	for (int i = 0; i < num_population; i++)
		genome[i]->run_network(num_samples, num_input, num_output, input, output);
};



void evolver::evolve_time_series(int num_samples, float** input, float** output, int window_size)
{
	//mutate
	select_best();

	mutate();

	//run
#pragma omp parallel for
	for (int i = 0; i < num_population; i++)
	{
		genome[i]->run_recursive(num_samples, num_input, num_output, input, output,window_size);
	}
};


float* evolver::predict_time_series(int num_samples, float** input, int window_size, float* data_point)
{
	select_best();
	return genome[0]->run_recursive_predict(num_samples, num_input, num_output, input, data_point, window_size);

};
