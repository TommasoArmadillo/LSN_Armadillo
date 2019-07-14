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
#include "main.h"

using namespace std;
 
int main (int argc, char *argv[]){
	Initialize();
	Present();
	
	Energy(0);
	
	int i_opt = 0;
	while(i_opt < passi_max && delta_H > precision){
		Gradient_Descent();
		i_opt++;
		if(i_opt%10 == 0) cout << "Step " << i_opt << endl << endl;
	}
	
	//print value and error for <H>
	Energy(1);
	
	cout << "Optimazed parameters " << endl;
	cout << "mu = " << mu << endl;
	cout << "sigma = " << sigma << endl;
	cout << "number of step with gradient descent = " << i_opt << endl;
	cout << "Value of <H> = " << ene << " pm " << err_ene << endl;
	
	Histogram();
	
	Output();
	
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
	
	//ReadInput from file input.dat
	ifstream ReadInput("input.dat");
	if (!ReadInput.is_open()){
		cerr << "PROBLEM: Unable to open Primes" << endl;
		exit(-1);
	}
	
	ReadInput >> mu;
	ReadInput >> sigma;
	ReadInput >> delta;
	ReadInput >> n_block;
	ReadInput >> n_step;
	ReadInput >> extremum;
	
	ReadInput.close();
	
	bin_size = 2. * extremum / n_bin;
}

void Present(void){
	//Print the parameters
	cout << "Parameters:" << endl;
	cout << "- mu:                 " << mu << endl;
	cout << "- sigma:              " << sigma << endl;
	cout << "- number of blocks:   " << n_block << endl;
	cout << "- step in each block: " << n_step << endl << endl;
	
}

double Error(double sum, double sum2, int nblock){
	return sqrt( (sum2/double(nblock) - pow((sum/double(nblock)), 2.) )/double(nblock) );
}

double Potential(double x){
	return pow(x, 4.) - (5./2.)*pow(x, 2.);
}

void Metropolis(void){
	x_old = x;
	x = x + rnd.Rannyu(-delta, delta);			//propose new positions
	
	p_old = PsiSquared(x_old, mu, sigma);			//calculate probabilities
	p_new = PsiSquared(x, mu, sigma);
	
	if(rnd.Rannyu() < min(p_new/p_old, 1.)) {		//accept? If no, come back to previous position
		accepted++;
	} else {
		x = x_old;
	}
	attempted++;
	return;
}

double PsiSquared(double x, double mu, double sigma){
	double exp1 = exp(-pow(x-mu, 2.)/(2.*sigma*sigma));
	double exp2 = exp(-pow(x+mu, 2.)/(2.*sigma*sigma));
	return pow(exp1+exp2, 2.);
}

void Energy(int print){
	ofstream H_out;
	if(print ==1){
		H_out.open("Energy.dat");
	}
	Reset();
	
	for(int iblk=0; iblk<n_block; iblk++) {				//cycle over blocks
		
		for(int jstp=0; jstp<n_step; jstp++) {			//cycle over throws in each block
			Metropolis();								//generates a value of x
			Update(iblk);
		}
		Average(iblk);
		if(print == 1){
			H_out << iblk+1 << " " << ene << " " << err_ene << endl;
		}
	}
	
	if(print ==1){
		H_out.close();
	}
	return;
}

void Reset (void){
	sum.assign(n_block, 0.);
	sum2.assign(n_block, 0.);
	
	return;
}

void Update(int iblk){
	double tmp = HPsiOverPsi(x);
	sum[iblk] += tmp;
	sum2[iblk] += tmp*tmp;
	return;
}

double HPsiOverPsi(double x){
	double exp1 = exp(-pow(x-mu, 2.)/(2.*sigma*sigma));
	double Der1 = (pow(x-mu, 2.) - sigma*sigma ) * exp1 / pow(sigma, 4.);
	double exp2 = exp(-pow(x+mu, 2.)/(2.*sigma*sigma));
	double Der2 = (pow(x+mu, 2.) - sigma*sigma ) * exp2 / pow(sigma, 4.);
	double Psi = exp1 + exp2;
	double HPsi = -(1./2.)*(Der1+Der2)+Potential(x)*Psi;
	
	return HPsi/Psi;
}

void Average(int iblk){
	sum[iblk]  /= double(n_step);
	sum2[iblk] /= double(n_step);
	
	double stima_ene=0., stima_ene2=0.;
	for(int i=0; i<=iblk; i++){
		stima_ene += sum[i];
		stima_ene2 += sum2[i];
	}
	
	err_ene = Error(stima_ene/(double)(iblk+1), stima_ene2/(double)(iblk+1), iblk+1);
	ene = stima_ene/(double)(iblk+1);
	//cout << "Stima energia = " << ene << " pm " << err_ene << endl;
	//cout << "Acceptance rate = " << double(accepted)/double(attempted) << endl << endl;
}

void Gradient_Descent(void){
	//variate mu
	mu += d_mu;
	double ene_old = ene;
	
	//calculate new energy
	Energy(0);
	
	//compute gradient
	grad[0] = (ene - ene_old)/d_mu;
	
	//variate mu
	mu -= d_mu;
	sigma += d_sigma;
	
	//calculate new energy
	Energy(0);
	
	//compute gradient
	grad[1] = (ene - ene_old)/d_sigma;
	
	mu = mu - alpha * grad[0];
	sigma -= d_sigma;
	sigma = sigma - alpha * grad[1];
	
	Energy(0);
	
	delta_H = abs(ene - ene_old);
}

void Histogram(void){
	
	ofstream Hist_out;
	Hist_out.open("histogram.dat");
	
	for(int i=0; i<n_points; i++){
		Metropolis();
		Hist_out << x << endl;
	}
	
	Hist_out.close();
}

void Output(void){
	ofstream Out_par;
	Out_par.open("parametres.dat");
	Out_par << mu << endl;
	Out_par << sigma << endl;
	Out_par << ene << " pm " << err_ene << endl;
	Out_par.close();
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
