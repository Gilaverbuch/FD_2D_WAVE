#pragma once
#include <vector>

class model_parameters{
public:
	int l, nsteps, dx, dy,  elements, layers, x_s, y_s, print_int, time;
	double dt, frequency, source_time_delay, stability; 
	std::vector<int>  x_range;//dynamic vectors to read velocity profile
	std::vector<double>  density_range , velocity_range; 

	model_parameters();	// Custom default constructor
	void print_parameters();

private:
	double extractNumbersWords(std::string str); 
	void read_vel_profile();

};