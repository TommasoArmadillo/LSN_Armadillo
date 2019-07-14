#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include "functions.h"

using namespace std;

double error (double AV, double AV2, int n) {
	if (n==0)
        return 0;
    else
        return sqrt((AV2 - AV*AV)/n);
}

double Psi100 (double x, double y, double z) {
	double r = sqrt(x*x + y*y + z*z);	
	return exp(- 2. * r) / M_PI;
}

double Psi210 (double x, double y, double z) {
	double r = sqrt(x*x + y*y + z*z);	
	double cos_theta = z/r;
	return r * r * exp(-r) * cos_theta * cos_theta / (32 * M_PI);
}

bool accept (double rand, double p_old, double p_new){
	if (p_old == 0.)
		return false;
	
	if( rand > min(p_new/p_old, 1.) )
		return false;

	return true;
}
