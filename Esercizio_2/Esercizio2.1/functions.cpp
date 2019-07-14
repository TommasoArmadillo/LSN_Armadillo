#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include "functions.h"

using namespace std;

double eval_f (double x) {
	return (M_PI/2.) * cos( (M_PI * x) / 2.);
}

double eval_f_importance (double x) {
	return (M_PI/2.) * cos( (M_PI * x) / 2.) / (2 * (1. - x));
}

double error (double AV, double AV2, int n) {
	if (n==0)
        return 0;
    else
        return sqrt((AV2 - AV*AV)/n);
}

