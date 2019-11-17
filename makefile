all: main

main: wave.o model_parameters.o Field.o
		g++ -g -Wall -o wave wave.o model_parameters.o Field.o

model_parameters.o: model_parameters.cpp
		g++ -c -g -Wall model_parameters.cpp

Field.o: Field.cpp
		g++ -c -g -Wall Field.cpp

wave.o: wave.cpp
		g++ -c -g -Wall wave.cpp 

clean:	
	 rm *.o