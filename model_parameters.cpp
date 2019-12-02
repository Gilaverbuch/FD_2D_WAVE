#include <iostream>
#include <cmath>
#include <fstream>
#include <sstream> 
#include <vector>
#include "model_parameters.h"


// ------------------------------------------------------------------------------------------------------------------
// ------------------------------------------------------------------------------------------------------------------
//public functions
model_parameters::model_parameters(){

	std::cout << "Cuctum default constructor. Reading parameters from input.txt" << std::endl;

	std::ifstream inFile;
	std::string line;
	double data[9], v_max;
	int i=0;



	inFile.open("input.txt");
	while(std::getline(inFile, line)) 
	  {
	  	data[i] = extractNumbersWords(line);
	  	i = i + 1;
	  	// std::cout << extractNumbersWords(line) << " " << data[i] << std::endl;
	  }

	time = int(data[0]);  
	dx = int(data[1]);
	dy = int(data[2]); 
	stability = data[3];
	frequency = data[4];
	source_time_delay = data[5];
	x_s = data[6];
	y_s = data[7];
	print_int = data[8];

	read_vel_profile();
	l = x_range[layers-1];
	elements_x = int (l/dx);
	elements_y = elements_x;

	v_max = 0;
	for (i=0; i<layers; i++){
		if (velocity_range[i] > v_max){
			v_max = velocity_range[i];
		}
	}

	dt = (stability * dx)/v_max;
	nsteps = time/dt; 

	
}
// ------------------------------------------------------------------------------------------------------------------


void model_parameters::print_parameters(){

	std::cout << "###############################" << std::endl;
	std::cout << "MODEL PARAMETERS:" << std::endl;
	std::cout << "length" << " " << l << " " << "m" << std::endl;
	std::cout << "dx" << " " << dx << " " << "m" << std::endl;
	std::cout << "dy" << " " << dx << " " << "m" << std::endl;
	std::cout << "source position (x,y)" << " " << x_s << " " << "m" << " " << y_s << " " << "m" << std::endl;
	std::cout << "number of elements" << " " << elements_x <<"X" << elements_y << std::endl;
	std::cout << "simulation time" << " " << time << " " << "sec" <<  std::endl;
	std::cout << "dt" << " " << dt << " " << "sec" << std::endl;
	std::cout << "frequency" << " " << frequency << " " << "Hz" << std::endl;
	std::cout << "source_time_delay" << " " << source_time_delay << " " << "sec" << std::endl;
	std::cout << "number of time steps" << " " << nsteps << std::endl;
	std::cout << "printing interval" << " " << print_int << std::endl;
	std::cout << "###############################" << std::endl;

}

// ------------------------------------------------------------------------------------------------------------------
// ------------------------------------------------------------------------------------------------------------------
//private functions
double model_parameters::extractNumbersWords(std::string str){ 
    std::stringstream ss;     
  
    /* Storing the whole string into string stream */
    ss << str; 
  
    /* Running loop till the end of the stream */
    std::string temp; 
    double found, val; 
    while (!ss.eof()) { 
  
        // extracting word by word from stream 
        ss >> temp; 
  
        /* Checking the given word is integer or not */
        if (std::stringstream(temp) >> found) 
            // std::cout << found << std::endl;
  			val=found;
        /* To save from space at the end of string */
        temp = ""; 
    } 
   	return val;
} 
// ------------------------------------------------------------------------------------------------------------------

void model_parameters::read_vel_profile(){

	std::ifstream inFile;

	inFile.open("vel_profile.txt");

	while(!inFile.eof())
	{
	    int  a;
	    double b, c;
	    inFile >> a >> b >> c; // extracts 3 values for layer


	    x_range.push_back(a);
		density_range.push_back(b);
		velocity_range.push_back(c);
	}
	layers = x_range.size();

	std::cout << layers << " " << velocity_range[0] << " " << velocity_range[1] << std::endl;
}




















