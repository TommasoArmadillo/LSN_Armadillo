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
	
	//PART1 -> Put_option and call_option in one step case t=0 --> 1

	int M=1E4;            // Total number of throws
	int N=100;            	// Number of blocks
	int L=M/N;    			// Number of throws in each block



	vector<double> sum_put;				//vector to store mean and mean of squares of each block
	vector<double> sum_squares_put;
	vector<double> sum_call;				
	vector<double> sum_squares_call;
	
	double S_t, put, call;								//temporary variable to store the random number
	double av_put = 0., av_sq_put = 0.,err_put = 0.;		//variable for average, average of squares and error
	double av_call = 0., av_sq_call = 0.,err_call = 0.;	

	//Parametres
	double S_0=100.;
	double T=1.;
	double K=100.;
	double r=0.1;
	double sig=0.25; 


	ofstream Put_dir("Put_direct.txt");				//file to write results
   	if (!Put_dir.is_open()) cerr << "PROBLEM: Unable to open Put_direct.txt" << endl;

	ofstream Call_dir("Call_direct.txt");				//file to write results
   	if (!Call_dir.is_open()) cerr << "PROBLEM: Unable to open Call_direct.txt" << endl;

	for(int i=0; i<N; i++) {				//cycle over blocks
		sum_put.push_back(0.);				//initialize at 0.
		sum_squares_put.push_back(0.);
		sum_call.push_back(0.);				//initialize at 0.
		sum_squares_call.push_back(0.);

		for(int j=0; j<L; j++) {			//cycle over throws in each block
			S_t = S_0 * exp( (r - sig*sig*0.5)*T + sig*rnd.Gauss(0., T) );			//generate number
			put = exp(-r*T) * max(0., K-S_t);			
			call = exp(-r*T) * max(0., S_t-K);

			sum_put[i] += put;
			sum_squares_put[i] += put*put;

			sum_call[i] += call;
			sum_squares_call[i] += call*call;
			
		}
		
		sum_put[i] = sum_put[i]/double(L);					//sum and sum_squares now contains means of the block of throws
		sum_squares_put[i] = sum_squares_put[i]/double(L);
		sum_call[i] = sum_call[i]/double(L);					
		sum_squares_call[i] = sum_squares_call[i]/double(L);

		av_put = 0.;						//computing mean and means of squares with L*(i+1) throws
		av_sq_put = 0.;
		av_call = 0.;						//computing mean and means of squares with L*(i+1) throws
		av_sq_call = 0.;
		
		for(int j = 0; j < i+1; j++) {	
			av_put += sum_put[j];
			av_sq_put += sum_squares_put[j];

			av_call += sum_call[j];
			av_sq_call += sum_squares_call[j];
		}
		
		av_put = av_put / double(i+1);
		av_sq_put = av_sq_put / double(i+1);
		err_put = error(av_put, av_sq_put, L*(i+1)-1);	//estimating error

		av_call = av_call / double(i+1);
		av_sq_call = av_sq_call / double(i+1);
		err_call = error(av_call, av_sq_call, L*(i+1)-1);
		
		Put_dir << L*(i+1) << " " << av_put << " " << err_put << endl;
		Call_dir << L*(i+1) << " " << av_call << " " << err_call << endl;
	}

	Put_dir.close();
	Call_dir.close();

	cout << "I saved all datas in Put_dir.txt mean value and error for Put direct using this format -> N Put Puterr (separated by blanck space)" << endl;
	cout << "I saved all datas in Call_dir.txt mean value and error for Call direct using this format -> N Call Callerr (separated by blanck space)" << endl;
	
	//PART2 -> Put_option and call_option in 100 step case t=0 --> t=0.01 --> .. --> t=1
	int n_step = 100;


	ofstream Put_dis("Put_discrete.txt");				//file to write results
   	if (!Put_dis.is_open()) cerr << "PROBLEM: Unable to open Put_discrete.txt" << endl;

	ofstream Call_dis("Call_discrete.txt");				//file to write results
   	if (!Call_dis.is_open()) cerr << "PROBLEM: Unable to open Call_discrete.txt" << endl;

	for(int i=0; i<N; i++) {				//cycle over blocks
		sum_put.push_back(0.);				//initialize at 0.
		sum_squares_put.push_back(0.);
		sum_call.push_back(0.);				//initialize at 0.
		sum_squares_call.push_back(0.);

		for(int j=0; j<L; j++) {			//cycle over throws in each block

			S_t = S_0;
			for(int l=0; l<n_step; l++) {
				S_t = S_t * exp( (r - sig*sig*0.5)*T/double(n_step) + sig*rnd.Gauss(0., 1.)*sqrt(T/double(n_step)) );			//generate number
			}
			
			put = exp(-r*T) * max(0., K-S_t);			
			call = exp(-r*T) * max(0., S_t-K);

			sum_put[i] += put;
			sum_squares_put[i] += put*put;

			sum_call[i] += call;
			sum_squares_call[i] += call*call;
			
		}
		
		sum_put[i] = sum_put[i]/double(L);					//sum and sum_squares now contains means of the block of throws
		sum_squares_put[i] = sum_squares_put[i]/double(L);
		sum_call[i] = sum_call[i]/double(L);					
		sum_squares_call[i] = sum_squares_call[i]/double(L);

		av_put = 0.;						//computing mean and means of squares with L*(i+1) throws
		av_sq_put = 0.;
		av_call = 0.;						//computing mean and means of squares with L*(i+1) throws
		av_sq_call = 0.;
		
		for(int j = 0; j < i+1; j++) {	
			av_put += sum_put[j];
			av_sq_put += sum_squares_put[j];

			av_call += sum_call[j];
			av_sq_call += sum_squares_call[j];
		}
		
		av_put = av_put / double(i+1);
		av_sq_put = av_sq_put / double(i+1);
		err_put = error(av_put, av_sq_put, L*(i+1)-1);	//estimating error

		av_call = av_call / double(i+1);
		av_sq_call = av_sq_call / double(i+1);
		err_call = error(av_call, av_sq_call, L*(i+1)-1);
		
		Put_dis << L*(i+1) << " " << av_put << " " << err_put << endl;
		Call_dis << L*(i+1) << " " << av_call << " " << err_call << endl;
	}

	Put_dis.close();
	Call_dis.close();

	cout << "I saved all datas in Put_discrete.txt mean value and error for Put discrete using this format -> N Put Puterr (separated by blanck space)" << endl;
	cout << "I saved all datas in Call_discrete.txt mean value and error for Call discrete using this format -> N Call Callerr (separated by blanck space)" << endl;
	
	
	
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
