#include <iostream>
#include <cmath>
#include <fstream>
#include <sstream> 
#include <vector>
#include <filesystem>
#include "model_parameters.h"
#include "Field.h"


// --------------------------------------------------------------------------------------------

int main()
{	

	model_parameters model;
	model.print_parameters();
	Field displacement(model);
	// displacement.Propagator();

	// std::cout << "number of elements" << " " << elements_ << std::endl;
	// int n=10;
	// int a[n] = {};

	// std::cout << a[0] << " " << a[5] << std::endl;

	return 0;
}



