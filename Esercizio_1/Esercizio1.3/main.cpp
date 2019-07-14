/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/

#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <vector>
#include "random.h"
#include "functions.h"

using namespace std;
 
int main (int argc, char *argv[]){

	//Random generator
   	Random rnd;
   	int seed[4];
   	int p1, p2;
   	ifstream Primes("Primes");
   	if (Primes.is_open()){
		Primes >> p1 >> p2 ;
	} else cerr << "PROBLEM: Unable to open Primes" << endl;
	Primes.close();

	ifstream input("seed.in");
	string property;
	if (input.is_open()){
		while ( !input.eof() ){
			input >> property;
			if( property == "RANDOMSEED" ){
				input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
            		rnd.SetRandom(seed,p1,p2);
         		}
		}
     	input.close();
	} else cerr << "PROBLEM: Unable to open seed.in" << endl;
	
	//estimating PI value

	int M=100000;            // Total number of throws
	int N=100;            	// Number of blocks
	int num=M/N;    			// Number of throws in each block (every estimation of pi is made with num throws)
	
	double MIN=0., MAX=10.;	//Area in which we can generate the needle
	double d=1., L=0.7;		//parametres d-distance lines, L-lenght of the needle
	int success = 0;			//How many intersection we have
	vector<double> Pi;				//vector to store mean and mean of squares of each block
	vector<double> Pi_squares;
	
	double y, theta; 						//temporary variable to store the y-centre of the needle and the angle that form with x-axis
	
	double av = 0., av_sq = 0.,err = 0.;		//variable for average, average of squares and error
	
	ofstream Pi_file("Pi.txt");				//file to write results
   	if (!Pi_file.is_open()) cerr << "PROBLEM: Unable to open Pi.txt" << endl;
	
	for(int i=0; i<N; i++) {				//cycle over blocks
		Pi.push_back(0.);				//initialize at 0.
		Pi_squares.push_back(0.);
		success = 0;
		
		for(int j=0; j<num; j++) {			//I throw num needles and make an estimation of pi
			y = rnd.Rannyu(MIN, MAX);		//generate y-centre and theta
			
			double x_box=rnd.Rannyu(-1.,1.), y_box=rnd.Rannyu(-1., 1.);
			while(sqrt(x_box*x_box + y_box*y_box)>1){
				x_box=rnd.Rannyu(-1.,1.); y_box=rnd.Rannyu(-1., 1.);
			}
			theta = atan(y_box/x_box);
			
			if(cross(y, theta, L, d, MIN, MAX)){		//Does the needle cross a line?
				success++;
			}
		}
		
		
		Pi[i] = 2 * L * num / (success * d); 			//estimate pi
		Pi_squares[i] = Pi[i]*Pi[i];
		
		av = (av * double(i) + Pi[i]) / double(i+1);
		av_sq = (av_sq * double(i) + Pi_squares[i]) / double(i+1);
		
		err = error(av, av_sq, (i+1)-1);		//estimating error
		
		//write of file
		Pi_file << num*(i+1) << " " << av << " " << err << endl;
		
	}

	Pi_file.close();

	cout << "I saved all datas in Pi.txt using this format -> N_groups PI PI_err (separated by blanck space)" << endl;
	
	
   	rnd.SaveSeed();
   	return 0;
}

/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
