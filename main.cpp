#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <math.h>

#include "evolver.h"

int main(void) {

	//int innovation_number = 1;

	srand(time(NULL));

	float** in;
	float** out;
	in = new float*[4];
	in[0] = new float[2];
	in[0][0] = 0;
	in[0][1] = 0;
	in[1] = new float[2];
	in[1][0] = 0;
	in[1][1] = 1;
	in[2] = new float[2];
	in[2][0] = 1;
	in[2][1] = 0;
	in[3] = new float[2];
	in[3][0] = 1;
	in[3][1] = 1;

	out = new float*[4];
	out[0] = new float[1];
	out[0][0] = 0;
	out[1] = new float[1];
	out[1][0] = 1;
	out[2] = new float[1];
	out[2][0] = 1;
	out[3] = new float[1];
	out[3][0] = 0;

	evolver* ev = new evolver;
	ev->instantiate(2, 1, 128, 16);
	while ((1==1))
	{
		ev->evolve(4,in,out);
		ev->print_best();
	}

	
	return 0;
}
