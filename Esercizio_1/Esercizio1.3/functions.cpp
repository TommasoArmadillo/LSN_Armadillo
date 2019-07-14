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

bool cross(double y, double theta, double L, double d, double MIN, double MAX) {
	double y_min, y_max, y_line = MIN;
	y_max = y+(L/2.)*sin(theta);
	y_min = y-(L/2.)*sin(theta);
	if(y_min > y_max){
		double tmp=y_min;
		y_min=y_max;
		y_max=tmp;
	}
	
	while(y_line < MAX+L/2.){
		if(y_min < y_line && y_max > y_line)
			return true;
		y_line += d;
	}
	return false;
}
