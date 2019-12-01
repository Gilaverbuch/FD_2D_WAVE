# FD_2D_WAVE
2D finite difference wave propagation. In progress...For now it simulates propagation in horizontally layered media

To compile the program run "make" in the terminal. It uses GNU g++ compiler. 
"./wave" will run the program. Results are saved to the results directory. 

The modeling parameters (time, dx_ , dy_, stability, frequency_, s_time_delay_ ,  x_source, y_source, print_interval ) are specified in "input.txt" in this order.

velocity profile, "vel_profile.txt", contains three arguments in the following order: distance[m]	density[kg/m^3]		velocity [m/s].	

Distance defines  the layers. 

Easy plotting with gnuplot. Type the following commends in the terminal:

-cd results/	

-gnuplot		

-do for [i=1:set_number] {plot  "wave_signal".i.".txt" with image ; pause 0.1; set xrange[0:x_max]; ;set yrange[0:y_max]};

