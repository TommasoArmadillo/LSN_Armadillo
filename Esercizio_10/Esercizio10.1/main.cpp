#include <iostream>
#include <fstream>
#include <ostream>
#include <cmath>
#include <iomanip>
#include <string>
#include <vector>
#include "city.h"
#include "random.h"
#include "cromosome.h"
#include "main.h"

using namespace std;

int main(){ 


 	Initialize();
	GenerateStartingPath();

	cout << "STARTING SIMULATED ANNEALING " << endl;

	for(int i=0; i<temp_vec.size(); i++){
		
		temp = temp_vec[i];
		beta = 1./temp;
		accepted = 0;
		attempted =0;
		
		for(int j=0; j<n_step; j++){
			Metropolis ();
		}
		cout << "Temp = " <<setw(5)<< temp << " -----> " << " BestPath Lenght= " <<Best.GetLen() << endl;
			
		cout << "Accettazione = " << double(accepted) / double(attempted) << endl << endl;
		Measure(i);
	}

	Output();

	return 0;
}

void Initialize() {
	//Initializing random generator
	system("bash clean.sh");
	
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
	
	//reading the cities
	ifstream ReadInput;
	ReadInput.open("cities.0");
	
	if (ReadInput.is_open()){
		double x, y;
		int index;
		
		for(int i=0; i<n_city; i++){
			ReadInput >> index >> x >> y;
			Cities_vec[i].SetX(x);
			Cities_vec[i].SetY(y);
			Cities_vec[i].SetIndex(index);
			//cout << Cities_vec[i].GetX() << " " << Cities_vec[i].GetY() << endl;
		}
		
		cout << "I have correctly read cities position from cities.0" << endl;

	} else {
		
		ofstream WriteCities;
		WriteCities.open("cities.0");

		for(int i=0; i<n_city; i++) {
			double theta = rnd.Rannyu(0.,2*M_PI);
			
			//Cities in a circle
			//Cities_vec[i] = City(r*cos(theta), r*sin(theta), i+1);
			
			//Cities in asquare
			Cities_vec[i] = City(rnd.Rannyu(), rnd.Rannyu(), i+1);
			
			WriteCities << Cities_vec[i].GetIndex() << " " << Cities_vec[i].GetX() << " " << Cities_vec[i].GetY() << endl;
		} 
		WriteCities.close();

		cout << "File 'cities.0' didn't exist, I created a new one" << endl;
	}
	ReadInput.close();

	double t=1.;
	temp_vec.push_back( t );
	
	for(int i=0; i<999; i++){
		t -= 0.001;
		temp_vec.push_back( t );
	}
	
	return;
}

void GenerateStartingPath() {
	
	
	Path.Shuffle(rnd, Cities_vec);
	Best = Path;

	cout << "I created a random paths" << endl<<endl;
	return;
}

void Metropolis(){
	
	Cromosome OldPath = Path;
	
	Mutate(Path);

	double p, len_old, len_new;

	len_old = OldPath.GetLen();
	len_new = Path.GetLen();
	
	if(len_new < Best.GetLen() ){
		Best = Path;
	}

	p = exp( beta * (len_old - len_new));
	
	if(p < rnd.Rannyu()) {
    	Path = OldPath;
    } else {
		accepted++;
	}

	attempted++;    

	return;
}

void Mutate(Cromosome& crom){
	int which_mutation = (int)rnd.Rannyu(0., 4.);
	
	switch(which_mutation){
		//inversion of two cities
		case 0:
			crom.Inversion(rnd);
			break;
			
		//finite shift of cities
		case 1:
			crom.FiniteShift(rnd);
			break;
		
		//permutate m cities with other m different cities
		case 2:
			crom.Permutation(rnd);
			break;
			
		//Reverse the order of m cities
		case 3:
			crom.Reverse(rnd);
			break;
			
		default:
			cout << "ERROR: Mutation unexpected ----> which_mutation = " << which_mutation << endl;
			break;
	}

	crom.Cost(Cities_vec);
	
	return;
}



void Measure(int index){
	ofstream SmallestPath;
	
	SmallestPath.open("SmallestPath.dat", ofstream::app);
	
	SmallestPath << index << " " << Best.GetLen() << endl;
	
	SmallestPath.close();
	
	return;
}



void Output(){
	ofstream CitiesOutput;
	CitiesOutput.open("cities.final");
	
	for(int i=0; i<n_city; i++){
		CitiesOutput << Cities_vec[ Best.Get(i) ].GetX() << " " << Cities_vec[ Best.Get(i) ].GetY() << endl;
	}
	
	CitiesOutput << Cities_vec[ Best.Get(0) ].GetX() << " " << Cities_vec[ Best.Get(0) ].GetY() << endl;
	
	cout << endl << "I wrote on file!" << endl;
	CitiesOutput.close();
	return;
}


