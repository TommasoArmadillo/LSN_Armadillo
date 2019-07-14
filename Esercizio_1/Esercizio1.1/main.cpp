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

	//Initialize Random generator
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
	
	//PART1 -> estimating mean and error for a uniform distribution in [0,1)

	int M=100000;            // Total number of throws
	int N=100;            	// Number of blocks
	int L=M/N;    			// Number of throws in each block

	vector<double> sum;				//vector to store mean and mean of squares of each block
	vector<double> sum_squares;
	
	double tmp;								//temporary variable to store the random number
	double av = 0., av_sq = 0.,err = 0.;		//variable for average, average of squares and error
	
	ofstream Mean("Mean.txt");				//file to write results
   	if (!Mean.is_open()) cerr << "PROBLEM: Unable to open Mean.txt" << endl;

	for(int i=0; i<N; i++) {				//cycle over blocks
		sum.push_back(0.);				//initialize at 0.
		sum_squares.push_back(0.);
		for(int j=0; j<L; j++) {			//cycle over throws in each block
			tmp = rnd.Rannyu();			//generate number
			sum[i] += tmp;
			sum_squares[i] += tmp*tmp;
			
		}
		
		sum[i] = sum[i]/double(L);					//sum and sum_squares now contains means of the block
		sum_squares[i] = sum_squares[i]/double(L);
		
		av = (av * double(i) + sum[i]) / double(i+1);
		av_sq = (av_sq * double(i) + sum_squares[i]) / double(i+1);
		
		err = error(av, av_sq, L*(i+1)-1);	//estimating error
		
		Mean << L*(i+1) << " " << av << " " << err << endl;
		
	}

	Mean.close();

	cout << "I saved all datas in Mean.txt mean value and error for r using this format -> X Y Yerr (separated by blanck space)" << endl;
	
	//PART2 -> Mean value and error for variance of a uniform distribution in [0,1)
	
	ofstream Variance("Variance.txt");		//file to write results
   	if (!Variance.is_open()) cerr << "PROBLEM: Unable to open Variance.txt" << endl;

	fill(sum.begin(), sum.end(), 0.);		//reinitialize the vectors at 0.
	fill(sum_squares.begin(), sum_squares.end(), 0.);
	av = 0.;
	av_sq = 0.;
	
	for(int i=0; i<N; i++) {			//cycle over blocks
		
		for(int j=0; j<L; j++) {		//cycle over throws in each block
			tmp = rnd.Rannyu();
			tmp = pow(tmp - 0.5, 2.);
			sum[i] += tmp;
			sum_squares[i] += tmp*tmp;
			
		}
		
		sum[i] = sum[i]/double(L);					//sum and sum_squares now contains means of the block of throws
		sum_squares[i] = sum_squares[i]/double(L);
		
		av = (av * double(i) + sum[i]) / double(i+1);
		av_sq = (av_sq * double(i) + sum_squares[i]) / double(i+1);
		
		err = error(av, av_sq, L*(i+1)-1);		//estimating error
		
		Variance << L*(i+1) << " " << av << " " << err << endl;
		
	}

	Variance.close();

	cout << "I saved all datas in Variance.txt mean value and error for variance using the format -> X Y Yerr (separated by blanck space)" << endl;
	
	//PART 3 -> ESTIMATING CHI^2
	M = 100;			//numbers of interval in which I divide [0,1]
	double width = 1./double(M);
	N = 1E4;			//number of data for each CHI2 test
	int n_test = 100;	//number of tests
	vector<int> frequency (M);		//numbers in the n-th interval
	vector<double> CHI2 (n_test);		//vector to collect all the CHI2
	
	ofstream CHI("CHI2.txt");		//file to write results
	if (!CHI.is_open()) cerr << "PROBLEM: Unable to open CHI2" << endl;
	
	for(int i = 0; i<n_test; i++){		//cycle for n_test CHI2 test
		
		fill(frequency.begin(), frequency.end(), 0);	//inizialize to 0.
		for(int j = 0; j<N; j++){		//cycle for each test I generate N numbers
			
			tmp = rnd.Rannyu();			//generate a number
			
			for(int k = 0; k<M; k++){
				if(is_in(double(k)*width, (double(k)+1)*width, tmp)) {frequency[k]++; break;}		//check in which interval it is and increment frequency
			}
			
		}
		
		for(int j = 0; j<M; j++){		//compute CHI2
			CHI2[i] += pow(frequency[j] - double(N)/double(M), 2);
		}
		
		CHI2[i] = CHI2[i] / (double(N)/double(M));
		
		CHI << i << " " << CHI2[i] << endl;
	}
	
	cout << "I saved all datas in CHI2.txt for the 100 test of CHI2 using the format -> X Y (separated by blanck space)" << endl;
	
	
	
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
