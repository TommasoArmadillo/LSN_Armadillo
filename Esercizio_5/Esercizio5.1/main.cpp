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

//Global variables
Random rnd;
int seed[4];
int p1, p2;

int M=100000;            								// Total number of throws
int N=100;            									// Number of blocks
int L=M/N;    												// Number of throws in each block

vector<double> sum;										//vector to store mean and mean of squares of each block
vector<double> sum_squares;

double delta = 1.2;										//set delta
double x_old = 0., y_old = 0., z_old = 0.;		//set initial position

double radius;												//temporary variable to store the random number
double x_new, y_new, z_new;
double p_old, p_new;
int acc = 0;
int n_equil = 100;
double av = 0., av_sq = 0.,err = 0.;				//variable for average, average of squares and error

//functions
void Initialize(void);
void Metropolis100(void);
void Metropolis210(void);

/*********************************************************/

int main (int argc, char *argv[]){
	
	Initialize();
	
	ofstream Psi100_file("Psi100.txt");						//file to write results
   if (!Psi100_file.is_open()) cerr << "PROBLEM: Unable to open Psi100.txt" << endl;

	ofstream Radius100_file("Radius100.txt");						//file to write results
   if (!Radius100_file.is_open()) cerr << "PROBLEM: Unable to open Radius100.txt" << endl;

	for(int j=0; j<n_equil; j++) {							//EQUILIBRATION
		Metropolis100();
	}
	
	
	for(int i=0; i<N; i++) {								//cycle over blocks
		sum.push_back(0.);									//initialize at 0.
		sum_squares.push_back(0.);
		
		for(int j=0; j<L; j++) {							//cycle over throws in each block
			Metropolis100();
			
			sum[i] += radius;
			sum_squares[i] += radius*radius;
			
			Psi100_file << x_new << " " << y_new << " " << z_new << endl;
		}
		
		sum[i] = sum[i]/double(L);					//sum and sum_squares now contains means of the block of throws
		sum_squares[i] = sum_squares[i]/double(L);
		
		av = (av * double(i) + sum[i]) / double(i+1);
		av_sq = (av_sq * double(i) + sum_squares[i]) / double(i+1);
		
		err = error(av, av_sq, i);	//estimating error
		
		Radius100_file << L*(i+1) << " " << av << " " << err << endl;
		
	}

	Psi100_file.close();
	Radius100_file.close();
	
	cout << "Acceptation rate: " << double(acc)/M << endl;
	cout << "I saved all datas in Psi100.txt mean value and error for r using this format -> X Y Z (separated by blanck space)" << endl;
	cout << "I saved all datas in Radius100.txt mean value and error for r using this format -> X Y Yerr (separated by blanck space)" << endl;
	
	
	//PART2 ->Psi n=2 l=1 m=0
		
	x_old = 0.;
	y_old = 0.;
	z_old = 2.5;
	delta = 2.8;
	acc = 0;

	ofstream Psi210_file("Psi210.txt");						//file to write results
   if (!Psi210_file.is_open()) cerr << "PROBLEM: Unable to open Psi210.txt" << endl;

	ofstream Radius210_file("Radius210.txt");				//file to write results
   if (!Radius210_file.is_open()) cerr << "PROBLEM: Unable to open Radius210.txt" << endl;

	for(int j=0; j<n_equil; j++){		//EQUILIBRATION
		Metropolis210();
	}
	
	for(int i=0; i<N; i++) {									//cycle over blocks
		sum[i] = 0.;												//initialize at 0.
		sum_squares[i] = 0.;
		
		for(int j=0; j<L; j++) {								//cycle over throws in each block
			Metropolis210();
			
			radius = sqrt(x_new*x_new + y_new*y_new + z_new*z_new);	//calculate r
			
			sum[i] += radius;
			sum_squares[i] += radius*radius;
			
			Psi210_file << x_new << " " << y_new << " " << z_new << endl;
		}
		
		sum[i] = sum[i]/double(L);					//sum and sum_squares now contains means of the block of throws
		sum_squares[i] = sum_squares[i]/double(L);
		
		av = (av * double(i) + sum[i]) / double(i+1);
		av_sq = (av_sq * double(i) + sum_squares[i]) / double(i+1);
		
		err = error(av, av_sq, i);	//estimating error
		
		Radius210_file << L*(i+1) << " " << av << " " << err << endl;
		
	}

	Psi210_file.close();
	Radius210_file.close();
	
	cout << "Acceptation rate: " << double(acc)/M << endl;
	cout << "I saved all datas in Psi210.txt mean value and error for r using this format -> X Y Z (separated by blanck space)" << endl;
	cout << "I saved all datas in Radius210.txt mean value and error for r using this format -> X Y Yerr (separated by blanck space)" << endl;
	
	
	
   	rnd.SaveSeed();
  	return 0;
}

void Initialize(void){
	//Random generator
	
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
	
	return;
}

void Metropolis100(void){
	x_new = x_old + rnd.Rannyu(-delta, delta);	//propose new positions
	y_new = y_old + rnd.Rannyu(-delta, delta);
	z_new = z_old + rnd.Rannyu(-delta, delta);
	
	p_old = Psi100(x_old, y_old, z_old);			//calculate probabilities
	p_new = Psi100(x_new, y_new, z_new);
	
	if(rnd.Rannyu() < min(p_new/p_old, 1.)) {						//accept? If no, come back to previous position
		acc++;
	} else {
		x_new = x_old;
		y_new = y_old;
		z_new = z_old;
	}
	
	radius = sqrt(x_new*x_new + y_new*y_new + z_new*z_new);	//calculate r
	
	x_old = x_new;											//prepare for next step
	y_old = y_new;
	z_old = z_new;
	return;
}

void Metropolis210(void){
	x_new = x_old + rnd.Rannyu(-delta, delta);	//propose new positions
	y_new = y_old + rnd.Rannyu(-delta, delta);
	z_new = z_old + rnd.Rannyu(-delta, delta);
	
	p_old = Psi210(x_old, y_old, z_old);			//calculate probabilities
	p_new = Psi210(x_new, y_new, z_new);
	
	if(rnd.Rannyu() < min(p_new/p_old, 1.)) {						//accept? If no, come back to previous position
		acc++;
	} else {
		x_new = x_old;
		y_new = y_old;
		z_new = z_old;
	}
	
	x_old = x_new;											//prepare for next step
	y_old = y_new;
	z_old = z_new;
	return;
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
