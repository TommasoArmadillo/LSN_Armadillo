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
	
	
	//PART 1 -> Generate 1E4 mean value for 1,2,10,100 value
	int M = 1E4;
	vector<int> N = {1, 2, 10, 100};
	int count = N.size();
	double sum = 0.;
	
	ofstream Dice("Dice.txt");
	if (!Dice.is_open()){
		cerr << "PROBLEM: Unable to open Dice.txt" << endl;
	}
	
	for(int i =0; i<count; i++){
		
		for(int j=0; j<M; j++){
			sum = 0.;
			for(int k=0; k<N[i]; k++){
				sum += rnd.Dice(6);
			}
			Dice << sum/double(N[i]) << endl;
		}
	}
	
	cout << "I wrote sum of 1, 2, 10, 100 variables originated from the throw of a standard dice in Dice.txt" << endl;
	Dice.close();
	
	
	double lambda = 1.;
	ofstream Exp("Exp.txt");
	if (!Exp.is_open()){
		cerr << "PROBLEM: Unable to open Exp.txt" << endl;
	}
	
	for(int i =0; i<count; i++){
		
		for(int j=0; j<M; j++){
			sum = 0.;
			for(int k=0; k<N[i]; k++){
				sum += rnd.Exp(lambda);
			}
			Exp << sum/double(N[i]) << endl;
		}
	}
	
	cout << "I wrote sum of 1, 2, 10, 100 variables originated from an exponential distribution with lambda = 1. in Exp.txt" << endl;
	Exp.close();
	
	
	double mu = 0., gamma = 1.;
	ofstream Lorentz("Lorentz.txt");
	if (!Lorentz.is_open()){
		cerr << "PROBLEM: Unable to open Lorentz.txt" << endl;
	}
	
	for(int i =0; i<count; i++){
		
		for(int j=0; j<M; j++){
			sum = 0.;
			for(int k=0; k<N[i]; k++){
				sum += rnd.Lorentz(mu, gamma);
			}
			Lorentz << sum/double(N[i]) << endl;
		}
	}
	
	cout << "I wrote sum of 1, 2, 10, 100 variables originated from a Lorentz distribution with mu = 0. and gamma = 1. in Lorentz.txt" << endl;
	Lorentz.close();
	
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
