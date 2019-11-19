#include <iostream>
#include <cmath>
#include <fstream>
#include <vector>
#include "Field.h"
#include "model_parameters.h"
#include <filesystem>

// ------------------------------------------------------------------------------------------------------------------
// ------------------------------------------------------------------------------------------------------------------
//public functions
Field::Field(model_parameters & M){
	int i,j;
	std::cout << "Cuctum default constructor. Initializing modeling parameters and velocity profile" << std::endl;

	dt = M.dt;
	dx = M.dx;
	dy = M.dy;
	frequency = M.frequency;
	source_time_delay = M.source_time_delay;
	steps = M.nsteps; 
	elements_x = M.elements_x;
	elements_y = M.elements_y;
	x_source = M.x_s;
	y_source = M.y_s;
	print_every = M.print_int;

	U = new double *[elements_y];
	U_past = new double *[elements_y];
	U_future = new double *[elements_y]; 
	lambda = new double *[elements_y];
	rho = new double *[elements_y];	
	vel = new double *[elements_y];	
	epsilon = new double *[elements_y]; 
	RHS = new double *[elements_y];
	sig = new double *[elements_y];

	for (i=0; i<elements_y; i++){

		U[i] = new double [elements_x];
		U_past[i] = new double [elements_x];
		U_future[i] = new double [elements_x]; 
		lambda[i] = new double [elements_x];
		rho[i] = new double [elements_x];	
		vel[i] = new double [elements_x];	
		epsilon[i] = new double [elements_x]; 
		RHS[i] = new double [elements_x];
		sig[i] = new double [elements_x];

	}

	x = new double [elements_x];
	y = new double [elements_y];




	initialize(elements_y, elements_x, rho, M.x_range, M.density_range);
	initialize(elements_y, elements_x, vel, M.x_range, M.velocity_range);



	
	// calculating lambda and x and y
	for (i=0; i<elements_y; i++){
		y[i] = (i-1) * dy;
		for (j=0; j<elements_x; j++){
			x[j] = (j-1) * dx;
			lambda[i][j] = rho[i][j] * pow(vel[i][j],2);
		}	
	}
	

}
// ------------------------------------------------------------------------------------------------------------------

// void Field::Propagator(){
	
// 	double source, A, B;
// 	int i, j, itteration, pos;
// 	std::cout << "wave propagation!!!" << std::endl;
// 	system("rm results/*.txt");

// 	// initial conditions
// 	// double a=5.5e-6;
// 	// for (i=0; i<elements; i++){
// 	// 	U[i] = exp(-a * pow((x[i] - x_source), 2));
// 	// 	U_past[i] = exp(-a * pow((x[i] - (x_source - vel[i]*dt)), 2));
// 	// }

// 	pos = int(x_source/dx);
// 	itteration = 0;
// 	print_to_file(elements, U, x, itteration);


// 	// propagator

// 	for (i=0; i<steps; i++){
// 		source = -2 * (i*dt - source_time_delay) * pow(frequency, 2) * exp(-1 * pow(frequency, 2) * pow((i*dt - source_time_delay), 2));

// 		for (j=1; j<(elements-1); j++){

// 			A = (U[j+1] - U[j])/( dx);
// 			B = (U[j] - U[j-1])/( dx);
// 			if (j==pos){
// 				RHS[j] = ((1/rho[j])*A - (1/rho[j-1])*B)/( dx)  + source * 1e-3;
// 			}
// 			else{
// 				RHS[j] = ((1/rho[j])*A - (1/rho[j-1])*B)/(dx);
// 			}
			
// 			U_future[j] = rho[j] * pow(vel[j],2) * pow(dt,2) * RHS[j] + 2*U[j] - U_past[j];
// 		}

// 		U_future[0] = 0;
// 		U_future[elements-1] = 0;

// 		if (i%print_every==0){
// 			itteration = i/print_every;

// 			print_to_file(elements, U_future, x, itteration);
// 		}

// 		for (j=0; j<elements; j++){
// 			U_past[j] = U[j];
// 			U[j] = U_future[j];
// 			U_future[j] = 0;
// 		}

// 	}

// 	system("mv wave_signal*.txt results/");	
// }



// ------------------------------------------------------------------------------------------------------------------
// ------------------------------------------------------------------------------------------------------------------
//private functions
void Field::initialize(int size_y, int size_x, double **A, std::vector<int>  x_range, std::vector<double> val){

	for (int i=0; i<size_y; i++){

		for (int j=0; j<size_x; j++){

			A[i][j] = val[0];
		}
	}
}
// ------------------------------------------------------------------------------------------------------------------

// void Field::print_to_file(int size, double *A, double *pos, int itteration){

//     std::ofstream ofile;

//     ofile.open("wave_signal"+std::to_string(itteration)+".txt"); 


//     for (int i=0; i<size; i++){
//          ofile << pos[i] << " " << A[i] << std::endl;
//        }
//     ofile.close();
// }
// ------------------------------------------------------------------------------------------------------------------























