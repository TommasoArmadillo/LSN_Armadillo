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

bool is_in (double a, double b, double x){
	if(a>b) cerr<<"Inserisci a<b!"<<endl;
	if(x>a && x<b) return true;
	return false;
}

