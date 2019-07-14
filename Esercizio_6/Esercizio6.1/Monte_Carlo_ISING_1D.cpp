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
#include <ostream>
#include <cmath>
#include <vector>
#include <iomanip>
#include "Monte_Carlo_ISING_1D.h"

using namespace std;

int main()
{
	
  	Input(); //Inizialization

	
	for(int i=0; i<=ntemp; i++)
	{
		temp = temp_ini + double(i)*temp_step;
		beta = 1./temp;
		cout << "Temperature = " << temp << endl;
		
		for(int iblk=1; iblk <= nblk; ++iblk) //Simulation
		{
			Reset(iblk);   //Reset block averages
			for(int istep=1; istep <= nstep; ++istep)
			{
				Move(metro, i);
				Measure(i);
				Accumulate(); //Update block averages
			}
			Averages(iblk);   //Print results for current block
		}
	}  
	ConfFinal(); //Write final configuration

  return 0;
}


void Input(void)
{
  ifstream ReadInput;

  cout << "Classic 1D Ising model             " << endl;
  cout << "Monte Carlo simulation             " << endl << endl;
  cout << "Nearest neighbour interaction      " << endl << endl;
  cout << "Boltzmann weight exp(- beta * H ), beta = 1/T " << endl << endl;
  cout << "The program uses k_B=1 and mu_B=1 units " << endl;

//Read seed for random numbers
   int p1, p2;
   ifstream Primes("Primes");
   Primes >> p1 >> p2 ;
   Primes.close();

   ifstream input("seed.in");
   input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
   rnd.SetRandom(seed,p1,p2);
   input.close();
  
//Read input informations
  ReadInput.open("input.dat");

  ReadInput >> temp_ini;
  cout << "Initial Temperature = " << temp_ini << endl;
	
	ReadInput >> temp_fin;
  cout << "Final Temperature = " << temp_fin << endl;

	ReadInput >> temp_step;
  cout << "Temperature Step= " << temp_step << endl;

  ReadInput >> nspin;
  cout << "Number of spins = " << nspin << endl;

  ReadInput >> J;
  cout << "Exchange interaction = " << J << endl;

  ReadInput >> h;
  cout << "External field = " << h << endl;
    
  ReadInput >> metro; // if=1 Metropolis else Gibbs

  ReadInput >> nblk;

  ReadInput >> nstep;
	
  ReadInput >> restart;
  if(restart) cout << "Starting simulation" << endl << endl;
  else cout << "Continuing simulation" << endl << endl;

  if(metro==1) cout << "The program perform Metropolis moves" << endl;
  else cout << "The program perform Gibbs moves" << endl;
  cout << "Number of blocks = " << nblk << endl;
  cout << "Number of steps in one block = " << nstep << endl << endl;
  ReadInput.close();

	ntemp = (temp_fin - temp_ini)/temp_step;

//Prepare arrays for measurements
  iu = 0; //Energy
  ic = 1; //Heat capacity
  im = 2; //Magnetization
  ix = 3; //Magnetic susceptibility
 
  n_props = 4; //Number of observables

	if(restart)
	{
		system("bash clean.sh");
		//initial random configuration
		for(int j=0; j<=ntemp; j++)
		{
			s.push_back( vector<double> () );
			for (int i=0; i<nspin; ++i)
			{
				s[j].push_back(0.);
				if(rnd.Rannyu() >= 0.5) s[j][i] = 1;
				else s[j][i] = -1;
			}
		}
	} else {
		system("rm -rf output.*");
		ifstream InputRead;
		InputRead.open("config.final");
		if (!InputRead.is_open()) cerr << "PROBLEM: Unable to open config.final" << endl;
		for(int j=0; j<=ntemp; j++)
		{
			s.push_back( vector<double> () );
			for (int i=0; i<nspin; ++i)
			{
				s[j].push_back(0.);
				InputRead >> s[j][i];
			}
		}
	}
	//Evaluate energy etc. of the initial configuration
	for(int i=0; i<=ntemp; i++)
	{
		Measure(i);
	}
	
	
}


void Move(int metro, int itemp)
{
  int o;
  double p, energy_old, energy_new;

  for(int i=0; i<nspin; ++i)
  {
	  //Select randomly a particle (for C++ syntax, 0 <= o <= nspin-1)
	  o = (int)(rnd.Rannyu()*nspin);

	  if(metro==1) //Metropolis
	  {
		  energy_old = Boltzmann(s[itemp][o],o,itemp);
		  //energy_old = Energy(itemp);
		  if(s[itemp][o] == 1) s[itemp][o] = -1;				//flip o-th spin
		  else s[itemp][o] = 1;
		
		  energy_new = Boltzmann(s[itemp][o],o, itemp);
		  //energy_new = Energy(itemp);
		  p = exp(- beta*(energy_new - energy_old));
		
		  if(rnd.Rannyu() < min(1., p)) accepted++;
		  else
		  {
			  if(s[itemp][o] == 1) s[itemp][o] = -1;				//re-flip o-th spin
			  else s[itemp][o] = 1;
		  }
		  attempted++;
	  }
	  else //Gibbs sampling
	  {
		  double p_to_be_plus, en_up, en_down;
		
		  if(h!=0.)
		  {
			  s[itemp][o] = +1.;
			  en_up  	= Energy(itemp);
			
			  s[itemp][o] = -1.;
			  en_down = Energy(itemp);
			
			  p_to_be_plus = exp(-beta*en_up)/( exp(-beta*en_up) + exp(-beta*en_down) );
			
		  }else
		  {
			  p_to_be_plus = 1./(1.+exp( -2.*beta*J*(s[itemp][Pbc(o-1)] + s[itemp][Pbc(o+1)])));
		  }
		
		  if(rnd.Rannyu() < p_to_be_plus) s[itemp][o] = +1.;
		  else s[itemp][o] = -1.;
	  }
  }
}

double Boltzmann(int sm, int ip, int itemp)
{
	double ene = -J * sm * ( s[itemp][Pbc(ip-1)] + s[itemp][Pbc(ip+1)] ) - h * sm;
	return ene;
}

double Energy(int itemp)
{
	double ene=0.;
	for(int i=0; i<nspin; i++)
	{
		ene += -J * s[itemp][i] * s[itemp][Pbc(i+1)] - (h/2.) * (s[itemp][i]+s[itemp][Pbc(i+1)]);
	}
	return ene;
}

void Measure(int itemp)
{

  	double u = 0.0, sm = 0.0;

	//cycle over spins
  	for (int i=0; i<nspin; ++i)
  	{
     	u += -J * s[itemp][i] * s[itemp][Pbc(i+1)] - 0.5 * h * (s[itemp][i] + s[itemp][Pbc(i+1)]);
		sm += s[itemp][i];
		
  	}
  	walker[iu] = u;
	walker[ic] = u*u;
	walker[im] = sm;
	walker[ix] = sm*sm;
	
}


void Reset(int iblk) //Reset block averages
{
   
   if(iblk == 1)
   {
       for(int i=0; i<n_props; ++i)
       {
           glob_av[i] = 0;
           glob_av2[i] = 0;
       }
   }

   for(int i=0; i<n_props; ++i)
   {
     blk_av[i] = 0;
   }
   blk_norm = 0;
   attempted = 0;
   accepted = 0;
}



void Accumulate(void) //Update block averages
{

   for(int i=0; i<n_props; ++i)
   {
     blk_av[i] = blk_av[i] + walker[i];
   }
   blk_norm = blk_norm + 1.0;
}


void Averages(int iblk) //Print results for current block
{
    
   	ofstream Ene, Heat, Mag, Chi;
	
    stima_u = blk_av[iu]/blk_norm; //Energy
    glob_av[iu]  += stima_u;
    glob_av2[iu] += stima_u*stima_u;
	
	stima_u_sq = blk_av[ic]/blk_norm;	//Heat Capacity
	stima_c = beta*beta* (stima_u_sq - stima_u*stima_u);
	glob_av[ic]  += stima_c;
	glob_av2[ic] += stima_c*stima_c;
	
	stima_m = blk_av[im]/blk_norm; //Magnetization
	glob_av[im]  += stima_m;
	glob_av2[im] += stima_m*stima_m;
	
	stima_s_sq = blk_av[ix]/blk_norm;	//Susceptibility
	stima_x = beta* stima_s_sq;
	glob_av[ix]  += stima_x;
	glob_av2[ix] += stima_x*stima_x;

	if(iblk == nblk)
	{
		Ene.open("output.ene.0", ofstream::app);
		Heat.open("output.heatcap.0", ofstream::app);
		Chi.open("output.susceptibility.0", ofstream::app);
		Mag.open("output.magnetization.0", ofstream::app);
		
		err_u=Error(glob_av[iu],glob_av2[iu],iblk);
		err_c=Error(glob_av[ic],glob_av2[ic],iblk);
		err_x=Error(glob_av[ix],glob_av2[ix],iblk);
		err_m=Error(glob_av[im],glob_av2[im],iblk);
		
		Ene << temp << " " << glob_av[iu]/(double)iblk << " " << err_u << endl;
		Heat << temp << " " << glob_av[ic]/(double)iblk << " " << err_c << endl;
		Chi << temp << " " << glob_av[ix]/(double)iblk << " " << err_x << endl;
		Mag << temp << " " << glob_av[im]/(double)iblk << " " << err_m << endl;
		
		Ene.close();
		Heat.close();
		Chi.close();
		Mag.close();
	}


    //cout << "----------------------------" << endl;
}


void ConfFinal(void)
{
  	ofstream WriteConf;

  	cout << "Print final configuration to file config.final " << endl << endl;
  	WriteConf.open("config.final");
	for(int j=0; j<=ntemp; j++)
	{
		for (int i=0; i<nspin; ++i)
		{
			WriteConf << s[j][i] << endl;
		}
	}
  	WriteConf.close();

  	rnd.SaveSeed();
}

int Pbc(int i)  //Algorithm for periodic boundary conditions
{
    if(i >= nspin) i = i - nspin;
    else if(i < 0) i = i + nspin;
    return i;
}

double Error(double sum, double sum2, int iblk)
{
    return sqrt((sum2/(double)iblk - pow(sum/(double)iblk,2))/(double)iblk);
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
