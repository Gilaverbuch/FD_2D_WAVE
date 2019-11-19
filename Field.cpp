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
		y[i] = i * dy;
		for (j=0; j<elements_x; j++){
			x[j] = j * dx;
			lambda[i][j] = rho[i][j] * pow(vel[i][j],2);

		}	
	}



}
// ------------------------------------------------------------------------------------------------------------------

void Field::Propagator(){
	
	double source, d2x, d2y;
	int i, j, k, itteration, pos_s_x, pos_s_y;
	std::cout << "wave propagation!!!" << std::endl;
	system("rm results/*.txt");


	pos_s_x = int(x_source/dx);
	pos_s_y = int(y_source/dy);
	itteration = 0;
	print_to_file(elements_x, elements_y, U, x, y, itteration);


	// propagator

	for (i=0; i<steps; i++){
		source = -2 * (i*dt - source_time_delay) * pow(frequency, 2) * exp(-1 * pow(frequency, 2) * pow((i*dt - source_time_delay), 2));
		

		for (j=1; j<(elements_y-1); j++){
			for (k=1; k<(elements_x-1); k++){

				d2x = (U[j][k-1] - 2*U[j][k] + U[j][k+1])/(pow(dx,2));

				d2y =  (U[j-1][k] - 2*U[j][k] + U[j+1][k])/(pow(dy,2));

				U_future[j][k] = 2*U[j][k] - U_past[j][k] + pow(vel[j][k],2) * pow(dt,2) * (d2y + d2x);
			}
		}

		U_future[pos_s_y][pos_s_x] = U_future[pos_s_y][pos_s_x] + pow(dt,2) * source;


		for (j=0; j<(elements_y-1); j++){
			U_future[j][0] = 0;
			U_future[j][elements_x-1] = 0;

		} 

		for (k=0; j<(elements_x-1); k++){
			U_future[0][k] = 0;
			U_future[elements_y-1][k] = 0;

		} 
		

		if (i%print_every==0){
			itteration = i/print_every;
			std::cout << "time is " << i*dt << std::endl;

			print_to_file(elements_x, elements_y, U_future, x, y, itteration);
		}

		for (j=0; j<elements_y; j++){
			for (k=0; k<elements_x; k++){

				U_past[j][k] = U[j][k];
				U[j][k] = U_future[j][k];
				U_future[j][k] = 0;
			}
		}

	}


	system("mv wave_signal*.txt results/");	
}



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

void Field::print_to_file(int size_x, int size_y, double **A, double *pos_x, double *pos_y, int itteration){

	int ii, jj;
    std::ofstream ofile;

    std::cout << " printing results " << std::endl;

    ofile.open("wave_signal"+std::to_string(itteration)+".txt"); 


    for (ii=0; ii<size_y; ii=ii+5){
    	for (jj=0; jj<size_x; jj=jj+5){

         	// std::cout << ii << " " << jj << " " <<  pos_x[jj] << " " << pos_y[ii] << std::endl;
         	ofile << pos_x[jj] << " " << pos_y[ii] << " " << A[ii][jj] << std::endl;
       }
   }
   // j=20;
   // for (i=0; i<size_y; i++){

   //   ofile << pos_y[i] << " " << A[i][j] << std::endl;
   // }
    ofile.close();
}
// ------------------------------------------------------------------------------------------------------------------























