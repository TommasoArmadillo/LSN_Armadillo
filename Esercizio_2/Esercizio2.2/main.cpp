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
	
	int M=1E4;            	// Total number of random walks
	int N=1E2;			//Number of steps in each random walk
	int dir;
	int tmp;
	
	ofstream Discrete("Discrete.txt");				//file to write results
   	if (!Discrete.is_open()) cerr << "PROBLEM: Unable to open Discrete.txt" << endl;
	
	ofstream Example("Example.txt");				//file to write results
	if (!Example.is_open()) cerr << "PROBLEM: Unable to open Example.txt" << endl;
	
	int dim=3;
	vector<int> position_d (dim, 0);	//vector that contains three dimensional discrete position
	vector<double> dist (N, 0.);		//vector in which I can store the distance
	vector<double> dist_sq (N, 0.);
	
	for(int i=0; i<M; i++) {				//cycle over M random walk
		position_d.assign(dim, 0);			//set the position back to the origin
		
		for(int j=0; j<N; j++) {					//cycle over N steps in each random walk
			dir = rnd.Direction(dim);			//generate a random direction 0,1,2 corresponding to x,y,z
			if(rnd.Rannyu() > 0.5)			 	//move in that direction with a 50% chance forward or backward
				position_d[dir]++;	
			else 
				position_d[dir]--;

			tmp = distance(position_d, dim);		//calculate distance from the origin
			
			dist[j] += tmp;					//add to previous walk
			dist_sq[j] += tmp*tmp;
			
			if(i<4){
				Example << position_d[0] << " " << position_d[1] << " " << position_d[2] << endl;
			}
		}
	}
	
	double av=0., av_sq=0., err=0.;

	
	
	for(int i=0; i<N; i++) {	
		av = dist[i] / double(M);			//I calculate the mean value of the distance squared and take sqrt
		av_sq = dist_sq[i] / double(M);
		err = sqrt(av_sq - av*av);
		Discrete << i+1 << " " << sqrt(av) << " " << err / (2. * sqrt(av) ) << endl;
	}

	Discrete.close();
	Example.close();
	
	cout << "I saved all datas in Discrete.txt " << endl;
	
	//PART 2 Random walk in continuum case
	ofstream Continuum("Continuum.txt");				//file to write results
   	if (!Continuum.is_open()) cerr << "PROBLEM: Unable to open Continuum.txt" << endl;
	
	double theta, phi;
	vector<double> position_c (dim, 0);	//vector that contains three dimensional continuos position
	dist.assign(N, 0.);
	dist_sq.assign(N, 0.);

	for(int i=0; i<M; i++) {				//cycle over M random walk
		position_c.assign(dim, 0);			//set the position back to the origin
		
		for(int j=0; j<N; j++) {					//cycle over N steps in each random walk
			theta = rnd.Theta();				//generate a random direction
			phi = rnd.Rannyu(0., 2.*M_PI); 
			
			position_c[0] += sin(theta) * cos(phi);
			position_c[1] += sin(theta) * sin(phi);
			position_c[2] += cos(theta);

			tmp = distance(position_c, dim);		//calculate distance from the origin
			
			dist[j] += tmp;					//add to previous walk
			dist_sq[j] += tmp*tmp;
		}
	}
	
	for(int i=0; i<N; i++) {	
		av = dist[i] / double(M);			//I calculate the mean value of the distance squared and take sqrt
		av_sq = dist_sq[i] / double(M);
		err = sqrt(av_sq - av*av);
		Continuum << i+1 << " " << sqrt(av) << " " << err / (2. * sqrt(av) ) << endl;
	}

	Continuum.close();

	cout << "I saved all datas in Continuum.txt " << endl;

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
