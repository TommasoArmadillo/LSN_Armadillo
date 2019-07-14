#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include <vector>
#include "functions.h"

using namespace std;

double distance (vector<int> position, int dim){
	double dist = 0.;
	for (int i=0; i<dim; i++)
		dist += position[i]*position[i];
	return dist;
}

double distance (vector<double> position, int dim){
	double dist = 0.;
	for (int i=0; i<dim; i++)
		dist += position[i]*position[i];
	return dist;
}
