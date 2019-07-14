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
	GeneratePopulation();
	
	Order(0, n_path-1);
	Measure(0);

	
	for(int gen=1; gen<=n_generation; gen++) {

		NextGeneration();
		
		if(gen%100 == 0) cout << "Generation number: " << gen << endl;
		
		Order(0, n_path-1);
		Measure(gen);

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
			Cities_vec[i] = City(r*cos(theta), r*sin(theta), i+1);
			
			//Cities in asquare
			//Cities_vec[i] = City(rnd.Rannyu(), rnd.Rannyu(), i+1);
			
			WriteCities << Cities_vec[i].GetIndex() << " " << Cities_vec[i].GetX() << " " << Cities_vec[i].GetY() << endl;
		} 
		WriteCities.close();

		cout << "File 'cities.0' didn't exist, I created a new one" << endl;
	}
	ReadInput.close();


	return;
}

void GeneratePopulation() {
	
	for(int i=0; i<n_path; i++){
		Population[i].Shuffle(rnd, i, Cities_vec);
	}
	cout << "I created 900 random paths" << endl;
	return;
}

void Order(int left, int right){
	int i = left, j = right;
	Cromosome tmp;
	
	double pivot = Population[(left + right) / 2].GetLen();
	
	/* partition */
	while (i <= j) {
		while (Population[i].GetLen() < pivot)
			i++;
		while (Population[j].GetLen() > pivot)
			j--;
		if (i <= j) {
			tmp = Population[i];
			Population[i] = Population[j];
			Population[j] = tmp;
			i++;
			j--;
		}
	}
	/* recursion */
	if (left < j)
		Order(left, j);
	if (i < right)
		Order(i, right);

	return;
}

int Choose(){
	double r=rnd.Rannyu();
	double p=1.3;			//p>1 privileges small numbers
	
	int number = int( double(n_path) * pow(r, p));
	while( number == 0){
		r=rnd.Rannyu();
		number = int( double(n_path) * pow(r, p));
	}
	
	return number;
}

void Measure(int generation){
	ofstream SmallestPath, MeanPath;
	
	SmallestPath.open("SmallestPath.dat", ofstream::app);
	MeanPath.open("MeanPath.dat", ofstream::app);
	
	int counter=0;
	double avg=0.;
	
	for(int i=0; i<(n_path/2); i++){
		avg += Population[i].GetLen();
		counter++;
	}
	
	SmallestPath << generation << " " << Population[0].GetLen() << endl;
	MeanPath << generation << " " << avg/(double)counter << endl;
	
	SmallestPath.close();
	MeanPath.close();
	
	return;
}

void Mutate(Cromosome& crom){
	int which_mutation = (int)rnd.Rannyu(0., 5.);
	
	switch(which_mutation){
		//inversion of two cities
		case 0:
			if(rnd.Rannyu() < prob_inversion)
				crom.Inversion(rnd);
			break;
		
		//shift of cities
		case 1:
			if(rnd.Rannyu() < prob_shift)
				crom.Shift(rnd);
			break;
			
		//finite shift of cities
		case 2:
			if(rnd.Rannyu() < prob_finiteshift)
				crom.FiniteShift(rnd);
			break;
		
		//permutate m cities with other m different cities
		case 3:
			if(rnd.Rannyu() < prob_permutation)
				crom.Permutation(rnd);
			break;
			
		//Reverse the order of m cities
		case 4:
			if(rnd.Rannyu() < prob_reverse)
				crom.Reverse(rnd);
			break;
			
		default:
			cout << "ERROR: Mutation unexpected ----> which_mutation = " << which_mutation << endl;
			break;
	}

	crom.Cost(Cities_vec);
	
	return;
}

void Output(){
	ofstream CitiesOutput;
	CitiesOutput.open("cities.final");
	
	for(int i=0; i<n_city; i++){
		CitiesOutput << Cities_vec[ Population[0].Get(i) ].GetX() << " " << Cities_vec[ Population[0].Get(i) ].GetY() << endl;
	}
	
	CitiesOutput << Cities_vec[ Population[0].Get(0) ].GetX() << " " << Cities_vec[ Population[0].Get(0) ].GetY() << endl;
	
	cout << endl << "I wrote on file!" << endl;
	CitiesOutput.close();
	return;
}

void NextGeneration(){
	vector<Cromosome> NewPopulation (n_path);
	
	//preserve the first two
	NewPopulation[0] = Population[0];
	NewPopulation[1] = Population[1];

	for(int i=2; i<n_path; i+=2){
		//choose 2 Cromosome from the previous generation
		int father = Choose();
		int mother = Choose();
		while (father == mother){
			mother = Choose();
		}
		
		NewPopulation[i]   = Population[father];
		NewPopulation[i+1] = Population[mother];
		
		if(rnd.Rannyu() < prob_crossover){
			NewPopulation[i].Crossover( rnd, NewPopulation[i+1] );
		}
		
		Mutate(NewPopulation[i]);
		Mutate(NewPopulation[i+1]);
		
		NewPopulation[i].Cost(Cities_vec);
		NewPopulation[i+1].Cost(Cities_vec);
		
	}
	
	Population = NewPopulation;
	
	return;
}
