#pragma once
#include "model_parameters.h"
#include <vector>

class Field{
public:
	double **U, **U_past, **U_future, **lambda, **rho, **vel, **epsilon, **RHS, **sig, *x, *y;
	double dt, frequency, source_time_delay;
	int dx, dy,  steps, elements_x, elements_y, x_source, y_source, print_every;
	Field(model_parameters&);	// Custom default constructor
	void Propagator(); //wave propagation
	
private:
	void initialize(int size_x, int size_y, double **A, std::vector<int>  x_range, std::vector<double> val);
	// void print_to_file(int size, double *A, double *pos, int itteration);
	

};