#include "header.h"
extern mt19937 rgen;
extern uniform_real_distribution<double> uran;

bool metropolis(double dE){
	bool acc=(dE<0);  //if dE<0, we accept!
	if(!acc){  //if not
		acc=uran(rgen)<exp(-dE);  // accept with probability equal to the Boltzman factor
	}
	return acc;//return that value
}
