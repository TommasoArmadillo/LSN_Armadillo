/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/

#ifndef __fluid_
#define __fluid_

//Random numbers
#include "random.h"
int seed[4];
Random rnd;

//parameters, observables
const int m_props=1000;
int n_props,iu,ic,im,ix,ig;
double nbins;
double walker[m_props];

// averages
double blk_av[m_props],blk_norm,accepted,attempted;
double glob_av[m_props],glob_av2[m_props];
double stima_u,stima_u_sq, stima_c,stima_m,stima_x,stima_g, stima_s_sq;
double err_u,err_c,err_m,err_x,err_g;

//configuration
const int m_spin=50;
std::vector< std::vector<double> > s;

// thermodynamical state
int nspin;
double beta,temp_ini, temp_fin, temp_step, ntemp, temp,J,h;

// simulation
int nstep, nblk, metro;
bool restart;

//functions
void Input(void);
void Reset(int);
void Accumulate(void);
void Averages(int);
void Move(int, int);
void ConfFinal(void);
void Measure(int);
double Boltzmann(int, int, int);
int Pbc(int);
double Error(double,double,int);
double Energy(int);

#endif

/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
