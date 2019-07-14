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
	
	//PART1 -> estimating integral and error using a uniform distribution in [0,1)

	int M=1E4;            	// Total number of throws
	int N=100;            	// Number of blocks
	int L=M/N;    			// Number of throws in each block

	vector<double> sum;				//vector to store mean and mean of squares of each block
	vector<double> sum_squares;
	
	double tmp;						//temporary variable to store the random number
	double av = 0., av_sq = 0.,err = 0.;			//variable for average, average of squares and error
	
	ofstream Uniform("Uniform.txt");				//file to write results
   	if (!Uniform.is_open()) cerr << "PROBLEM: Unable to open Uniform.txt" << endl;

	for(int i=0; i<N; i++) {				//cycle over blocks
		sum.push_back(0.);				//initialize at 0.
		sum_squares.push_back(0.);
		
		for(int j=0; j<L; j++) {					//cycle over throws in each block
			tmp = eval_f( rnd.Rannyu() );			//generate number and calculate f(x)
			sum[i] += tmp;
			sum_squares[i] += tmp*tmp;
			
		}
		
		sum[i] = sum[i] / double(L);					//sum and sum_squares now contains estimation of integral -- b-a=1. (not necessary)
		sum_squares[i] = sum_squares[i] / double(L);
		
		av = (av*double(i) + sum[i]) / double(i+1);
		av_sq = (av_sq*double(i) + sum_squares[i]) / double(i+1);
		
		err = error(av, av_sq, L*(i+1)-1);	//estimating error
		
		Uniform << L*(i+1) << " " << av << " " << err << endl;
		
	}

	Uniform.close();

	cout << "I saved all datas in Uniform.txt mean value and error for integral using this format -> #throws Integral err (separated by blanck space)" << endl;
	
	//PART2 -> estimating integral and error using a linear distribution in [0,1) p(x)=2(1-x)

	fill(sum.begin(), sum.end(), 0.);				//vector to store mean and mean of squares of each block
	fill(sum_squares.begin(), sum_squares.end(), 0.);
	
	ofstream Importance("Importance.txt");				//file to write results
   	if (!Importance.is_open()) cerr << "PROBLEM: Unable to open Importance.txt" << endl;

	for(int i=0; i<N; i++) {				//cycle over blocks
		
		for(int j=0; j<L; j++) {					//cycle over throws in each block
			tmp = rnd.Line();		
			
			tmp = eval_f_importance( tmp );			//generate number and calculate f(x)
			sum[i] += tmp;
			sum_squares[i] += tmp*tmp;
			
		}
		
		sum[i] = sum[i] / double(L);					//sum and sum_squares now contains estimation of integral -- b-a=1. (not necessary)
		sum_squares[i] = sum_squares[i] / double(L);
		
		av = 0.;						//computing mean and means of squares with L*(i+1) throws
		av_sq = 0.;
		
		for(int j = 0; j < i+1; j++) {	
			av += sum[j];
			av_sq += sum_squares[j];
		}
		
		av = av / double(i+1);
		av_sq = av_sq / double(i+1);
		
		err = error(av, av_sq, L*(i+1)-1);	//estimating error
		
		Importance << L*(i+1) << " " << av << " " << err << endl;
		
	}

	Importance.close();

	cout << "I saved all datas in Importance.txt mean value and error for integral using this format -> #throws Integral err (separated by blanck space)" << endl;
	
	
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
