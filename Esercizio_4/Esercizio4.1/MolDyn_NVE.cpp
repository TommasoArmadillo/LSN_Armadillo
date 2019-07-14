/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
#include <stdlib.h>     // srand, rand: to generate random number
#include <iostream>     // cin, cout: Standard Input/Output Streams Library
#include <fstream>      // Stream class to both read and write from/to files.
#include <cmath>        // rint, pow
#include <vector>
#include "random.h"
#include "MolDyn_NVE.h"

using namespace std;

int main(){
	
	//Random generator
	Random rnd;
	int seed[4];
	int p1, p2;
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
	
  	Input(rnd);             //Inizialization
  	int nconf = 1;
	
	
  	for(int iblock=0; iblock < nblock; iblock++){
		
		for(int jblock=0; jblock<lblock; jblock++) {
			int istep = iblock * lblock + jblock + 1;
			Move();           //Move particles with Verlet algorithm
			if(istep%iprint == 0) cout << "Number of time-steps: " << istep << endl;
			if(istep%10 == 0){
				Measure();     //Properties measurement
				Accumulate(iblock);
				
				
				//ConfXYZ(nconf);//Write actual configuration in XYZ format //Commented to avoid "filesystem full"!
				nconf += 1;
			}
			
		}
		Averages(iblock);
  	}
  	ConfFinal();         //Write final configuration to restart
	rnd.SaveSeed();
  	return 0;
}


void Input(Random rnd){ //Prepare all stuff for the simulation
  	ifstream ReadInput,ReadConf;

  	cout << "Classic Lennard-Jones fluid        " << endl;
  	cout << "Molecular dynamics simulation in NVE ensemble  " << endl << endl;
  	cout << "Interatomic potential v(r) = 4 * [(1/r)^12 - (1/r)^6]" << endl << endl;
  	cout << "The program uses Lennard-Jones units " << endl;

  	seed = 1;    //Set seed for random numbers
  
  	ReadInput.open("input.dat"); //Read input

  	ReadInput >> temp;

  	ReadInput >> npart;
  	cout << "Number of particles = " << npart << endl;

  	ReadInput >> rho;
  	cout << "Density of particles = " << rho << endl;
  	vol = (double)npart/rho;
  	cout << "Volume of the simulation box = " << vol << endl;
  	box = pow(vol,1.0/3.0);
  	cout << "Edge of the simulation box = " << box << endl;

  	ReadInput >> rcut;
  	ReadInput >> delta;
  	ReadInput >> nstep;
	ReadInput >> nblock;
  	ReadInput >> iprint;
	ReadInput >> cont;
	ReadInput >> rescale;

	lblock = nstep / nblock;
	
  	cout << "The program integrates Newton equations with the Verlet method " << endl;
  	cout << "Time step = " << delta << endl;
	cout << "Number of steps = " << nstep << endl;
	if(cont) cout << "Continuing simulation" << endl;
	else cout << "Restart simulation" << endl;
	if(rescale) cout << "Rescaling velocities" << endl;
	
	cout << endl;
  	ReadInput.close();

	//Prepare array for measurements
  	iv = 0; //Potential energy
  	ik = 1; //Kinetic energy
  	ie = 2; //Total energy
  	it = 3; //Temperature
	
	//measurement of g(r)
	igofr = 4;
	nbins = 100;
	n_props = n_props + nbins;
	bin_size = (box/2.0)/(double)nbins;
	vec_gofr.assign(nbins, 0.);
	vec_gofr_2.assign(nbins, 0.);
	
	double sumv2 = 0.0, fs;
	
	if(cont && rescale) {
		
		double xmdt[m_part],ymdt[m_part],zmdt[m_part];    //vectors r(t-dt)
		
		ifstream ReadConfOld;
		ReadConfOld.open("old.0");
		
		if (!ReadConfOld.is_open()){
			cerr << "File old.0 not found!" << endl;
			exit(1);
		}
		
		for (int i=0; i<npart; ++i){
			ReadConfOld >> xmdt[i] >> ymdt[i] >> zmdt[i];
			xmdt[i] = xmdt[i] * box;
			ymdt[i] = ymdt[i] * box;
			zmdt[i] = zmdt[i] * box;
			
			xold[i] = xmdt[i];
			yold[i] = ymdt[i];
			zold[i] = zmdt[i];
		}
		ReadConfOld.close();
		
		ReadConf.open("old.final");
		if (!ReadConf.is_open()){
			cerr << "File old.final not found!" << endl;
			exit(1);
		}
		for (int i=0; i<npart; ++i){
			ReadConf >> x[i] >> y[i] >> z[i];
			x[i] = x[i] * box;
			y[i] = y[i] * box;
			z[i] = z[i] * box;
		}
		ReadConf.close();
		
		Move();
		
		for(int i=0; i<npart; i++) {
			vx[i] = (x[i] - xmdt[i])/ (2 * delta);
			vy[i] = (y[i] - ymdt[i])/ (2 * delta);
			vz[i] = (z[i] - zmdt[i])/ (2 * delta);
		}
		
		for (int i=0; i<npart; ++i){
			sumv2 += vx[i]*vx[i] + vy[i]*vy[i] + vz[i]*vz[i];
		}
		sumv2 /= (double)npart;
		
		fs = sqrt(3 * temp / sumv2);   // fs = velocity scale factor
		for (int i=0; i<npart; ++i){
			vx[i] *= fs;
			vy[i] *= fs;
			vz[i] *= fs;
			
			xold[i] = x[i] - vx[i] * delta;
			yold[i] = y[i] - vy[i] * delta;
			zold[i] = z[i] - vz[i] * delta;
		}
		
		for (int i=0; i<npart; ++i){
			xmdt[i] = x[i] - 2 * delta * vx[i];
			ymdt[i] = y[i] - 2 * delta * vy[i];
			zmdt[i] = z[i] - 2 * delta * vz[i];
		}
		
		for (int i=0; i<npart; ++i){
			x[i] = xold[i];
			xold[i] = xmdt[i];
			y[i] = yold[i];
			yold[i] = ymdt[i];
			z[i] = zold[i];
			zold[i] = zmdt[i];
		}
		
	} else if (cont && !rescale) {
		
		system("bash clean.sh");
		
		//Read initial configuration
  		cout << "Read initial configuration from file config.final " << endl << endl;
  		ReadConf.open("config.final");
  		for (int i=0; i<npart; ++i){
    		ReadConf >> x[i] >> y[i] >> z[i];
    		x[i] = x[i] * box;
    		y[i] = y[i] * box;
    		z[i] = z[i] * box;
  		}
  		ReadConf.close();
	
		//Prepare initial velocities
   		cout << "Prepare random velocities with center of mass velocity equal to zero " << endl << endl;
   		double sumv[3] = {0.0, 0.0, 0.0};
   		for (int i=0; i<npart; ++i){
     		vx[i] = rnd.Rannyu() - 0.5;
     		vy[i] = rnd.Rannyu() - 0.5;
     		vz[i] = rnd.Rannyu() - 0.5;
			sumv[0] += vx[i];
			sumv[1] += vy[i];
			sumv[2] += vz[i];
		}
		for (int idim=0; idim<3; ++idim) sumv[idim] /= (double)npart;
		for (int i=0; i<npart; ++i){
			vx[i] = vx[i] - sumv[0];
			vy[i] = vy[i] - sumv[1];
			vz[i] = vz[i] - sumv[2];
		}
		
		for (int i=0; i<npart; ++i){
			sumv2 += vx[i]*vx[i] + vy[i]*vy[i] + vz[i]*vz[i];
		}
		sumv2 /= (double)npart;
		
		fs = sqrt(3 * temp / sumv2);   // fs = velocity scale factor
		for (int i=0; i<npart; ++i){
			vx[i] *= fs;
			vy[i] *= fs;
			vz[i] *= fs;
			
			xold[i] = x[i] - vx[i] * delta;
			yold[i] = y[i] - vy[i] * delta;
			zold[i] = z[i] - vz[i] * delta;
		}
		
	} else if(!cont) {
		system("bash clean.sh");
		
		//Read initial configuration
		cout << "Read initial configuration from file config.0 " << endl << endl;
		ReadConf.open("config.0");
		for (int i=0; i<npart; ++i){
			ReadConf >> x[i] >> y[i] >> z[i];
			x[i] = x[i] * box;
			y[i] = y[i] * box;
			z[i] = z[i] * box;
		}
		ReadConf.close();
		
		//Prepare initial velocities
		cout << "Prepare random velocities with center of mass velocity equal to zero " << endl << endl;
		double sumv[3] = {0.0, 0.0, 0.0};
		for (int i=0; i<npart; ++i){
			vx[i] = rnd.Rannyu() - 0.5;
			vy[i] = rnd.Rannyu() - 0.5;
			vz[i] = rnd.Rannyu() - 0.5;
			sumv[0] += vx[i];
			sumv[1] += vy[i];
			sumv[2] += vz[i];
		}
		for (int idim=0; idim<3; ++idim) sumv[idim] /= (double)npart;
		for (int i=0; i<npart; ++i){
			vx[i] = vx[i] - sumv[0];
			vy[i] = vy[i] - sumv[1];
			vz[i] = vz[i] - sumv[2];
		}
		
		for (int i=0; i<npart; ++i){
			sumv2 += vx[i]*vx[i] + vy[i]*vy[i] + vz[i]*vz[i];
		}
		sumv2 /= (double)npart;
		
		fs = sqrt(3 * temp / sumv2);   // fs = velocity scale factor
		for (int i=0; i<npart; ++i){
			vx[i] *= fs;
			vy[i] *= fs;
			vz[i] *= fs;
			
			xold[i] = Pbc(x[i] - vx[i] * delta);
			yold[i] = Pbc(y[i] - vy[i] * delta);
			zold[i] = Pbc(z[i] - vz[i] * delta);
		}
		
	}
	
	//average block
	vec_epot.assign(nblock, 0.0);				//vector to store mean and mean of squares of each block
	vec_epot_sq.assign(nblock, 0.0);
	vec_etot.assign(nblock, 0.0);
	vec_etot_sq.assign(nblock, 0.0);
	vec_ekin.assign(nblock, 0.0);
	vec_ekin_sq.assign(nblock, 0.0);
	vec_temp.assign(nblock, 0.0);
	vec_temp_sq.assign(nblock, 0.0);
	vec_press.assign(nblock, 0.0);
	vec_press_sq.assign(nblock, 0.0);
	
	
	
	return;
}


void Move(void){ //Move particles with Verlet algorithm
  double xnew, ynew, znew, fx[m_part], fy[m_part], fz[m_part];

  for(int i=0; i<npart; ++i){ //Force acting on particle i
    fx[i] = Force(i,0);
    fy[i] = Force(i,1);
    fz[i] = Force(i,2);
  }

  for(int i=0; i<npart; ++i){ //Verlet integration scheme

    xnew = Pbc( 2.0 * x[i] - xold[i] + fx[i] * pow(delta,2) );
    ynew = Pbc( 2.0 * y[i] - yold[i] + fy[i] * pow(delta,2) );
    znew = Pbc( 2.0 * z[i] - zold[i] + fz[i] * pow(delta,2) );

    vx[i] = Pbc(xnew - xold[i])/(2.0 * delta);
    vy[i] = Pbc(ynew - yold[i])/(2.0 * delta);
    vz[i] = Pbc(znew - zold[i])/(2.0 * delta);

    xold[i] = x[i];
    yold[i] = y[i];
    zold[i] = z[i];

    x[i] = xnew;
    y[i] = ynew;
    z[i] = znew;
  }
  return;
}

double Force(int ip, int idir){ //Compute forces as -Grad_ip V(r)
  double f=0.0;
  double dvec[3], dr;

  for (int i=0; i<npart; ++i){
    if(i != ip){
      dvec[0] = Pbc( x[ip] - x[i] );  // distance ip-i in pbc
      dvec[1] = Pbc( y[ip] - y[i] );
      dvec[2] = Pbc( z[ip] - z[i] );

      dr = dvec[0]*dvec[0] + dvec[1]*dvec[1] + dvec[2]*dvec[2];
      dr = sqrt(dr);

      if(dr < rcut){
        f += dvec[idir] * (48.0/pow(dr,14) - 24.0/pow(dr,8)); // -Grad_ip V(r)
      }
    }
  }
  
  return f;
}

void Measure(){ //Properties measurement
  	int count=0;
	double v, t, vij, p;
  	double dx, dy, dz, dr;
  	ofstream Epot, Ekin, Etot, Temp, Press;
	
	//reset the hystogram of g(r)
	for (int k=igofr; k<igofr+nbins; ++k) walker[k]=0.0;
	
  	Epot.open("output_epot.dat",ios::app);
  	Ekin.open("output_ekin.dat",ios::app);
  	Temp.open("output_temp.dat",ios::app);
  	Etot.open("output_etot.dat",ios::app);
  	Press.open("output_press.dat",ios::app);

  	v = 0.0; //reset observables
  	t = 0.0;
	p = 0.0;
//cycle over pairs of particles
  	for (int i=0; i<npart-1; ++i){
    	for (int j=i+1; j<npart; ++j){
     		count++;
			dx = Pbc( x[i] - x[j] );
			dy = Pbc( y[i] - y[j] );
			dz = Pbc( z[i] - z[j] );

     		dr = dx*dx + dy*dy + dz*dz;
     		dr = sqrt(dr);

			for (int k=igofr; k<igofr+nbins; ++k)
			{
				double rplusdr = (k-igofr+1) * bin_size;
				if(dr < rplusdr)
				{
					walker[k] += 1.;
					break;
				}
			}

     		if(dr < rcut){
       			vij = 4.0/pow(dr,12) - 4.0/pow(dr,6);

				//Potential energy
       			v += vij;
     		}
    	}
  	}

//Kinetic energy
  	for (int i=0; i<npart; ++i) t += 0.5 * (vx[i]*vx[i] + vy[i]*vy[i] + vz[i]*vz[i]);
	
	
    stima_pot = v/(double)npart; //Potential energy
    stima_kin = t/(double)npart; //Kinetic energy
    stima_temp = (2.0 / 3.0) * t/(double)npart; //Temperature
    stima_etot = (t+v)/(double)npart; //Total enery
	stima_press = rho*stima_temp + (16.*v)/(vol*count);

    Epot << stima_pot  << endl;
    Ekin << stima_kin  << endl;
    Temp << stima_temp << endl;
    Etot << stima_etot << endl;
	Press << stima_press << endl;
    Epot.close();
    Ekin.close();
    Temp.close();
    Etot.close();
	Press.close();
	
    return;
}


void ConfFinal(void){ //Write final configuration
	ofstream WriteOld, WriteOldFinal;
	
	cout << "Print final configuration to file old.0 and old.final to restart simulation" << endl << endl;
	
	WriteOld.open("old.0");
	for (int i=0; i<npart; ++i){
		WriteOld << x[i]/box << "   " <<  y[i]/box << "   " << z[i]/box << endl;
	}
	WriteOld.close();
	
	Move();
	
	WriteOldFinal.open("old.final");
	for (int i=0; i<npart; ++i){
		WriteOldFinal << x[i]/box << "   " <<  y[i]/box << "   " << z[i]/box << endl;
	}
	WriteOldFinal.close();
	
	ofstream WriteConf;
  	cout << "Print final configuration to file config.final " << endl << endl;
  	WriteConf.open("config.final");

  	for (int i=0; i<npart; ++i){
    	WriteConf << x[i]/box << "   " <<  y[i]/box << "   " << z[i]/box << endl;
  	}
  	WriteConf.close();
  	return;
}

void ConfXYZ(int nconf){ //Write configuration in .xyz format
  ofstream WriteXYZ;

  WriteXYZ.open("frames/config_" + to_string(nconf) + ".xyz");
  WriteXYZ << npart << endl;
  WriteXYZ << "This is only a comment!" << endl;
  for (int i=0; i<npart; ++i){
    WriteXYZ << "LJ  " << Pbc(x[i]) << "   " <<  Pbc(y[i]) << "   " << Pbc(z[i]) << endl;
  }
  WriteXYZ.close();
}

double Pbc(double r){  //Algorithm for periodic boundary conditions with side L=box
    return r - box * rint(r/box);
}

void Accumulate(int iblock){
	vec_epot[iblock] += stima_pot;
	vec_epot_sq[iblock] += stima_pot * stima_pot;
	vec_etot[iblock] += stima_etot;
	vec_etot_sq[iblock] += stima_etot * stima_etot;
	vec_ekin[iblock] += stima_kin;
	vec_ekin_sq[iblock] += stima_kin * stima_kin;
	vec_temp[iblock] += stima_temp;
	vec_temp_sq[iblock] += stima_temp * stima_temp;
	vec_press[iblock] += stima_press;
	vec_press_sq[iblock] += stima_press * stima_press;
	
}

void Averages(int iblock) {
	double nmeasure = double(lblock)/10.;
	
	vec_epot[iblock] /= nmeasure;
	vec_epot_sq[iblock] /= nmeasure;
	vec_etot[iblock] /= nmeasure;
	vec_etot_sq[iblock] /= nmeasure;
	vec_ekin[iblock] /= nmeasure;
	vec_ekin_sq[iblock] /= nmeasure;
	vec_temp[iblock] /= nmeasure;
	vec_temp_sq[iblock] /= nmeasure;
	vec_press[iblock] /= nmeasure;
	vec_press_sq[iblock] /= nmeasure;
	
	
	
	ofstream Epot, Ekin, Etot, Temp, Press, Gave;
	double epot_ave=0., ekin_ave=0., etot_ave=0., temp_ave=0., press_ave=0.;
	double epot_ave_sq=0., ekin_ave_sq=0., etot_ave_sq=0., temp_ave_sq=0., press_ave_sq=0.;
	cout << "Block number " << iblock+1 << endl;
	
	Epot.open("output_epot_ave.dat",ios::app);
	Ekin.open("output_ekin_ave.dat",ios::app);
	Etot.open("output_etot_ave.dat",ios::app);
	Temp.open("output_temp_ave.dat",ios::app);
	Press.open("output_press_ave.dat",ios::app);
	Gave.open("output.gave.0",ios::app);
	
	for(int i=0; i<=iblock; i++){
		epot_ave += vec_epot[i];
		epot_ave_sq += vec_epot_sq[i];
		etot_ave += vec_etot[i];
		etot_ave_sq += vec_etot_sq[i];
		ekin_ave += vec_ekin[i];
		ekin_ave_sq += vec_ekin_sq[i];
		temp_ave += vec_temp[i];
		temp_ave_sq += vec_temp_sq[i];
		press_ave += vec_press[i];
		press_ave_sq += vec_press_sq[i];
	}
	
	//g(r)
	for (int k=igofr; k<igofr+nbins; ++k)
	{
		double r = bin_size * double(k-igofr+1);
		double deltaV = ((4.*M_PI)/3.) * ( pow(r+bin_size, 3.) - pow(r,3.) );
		walker[k] /= (rho * npart * deltaV);
		
		vec_gofr[k] += walker[k];
		vec_gofr_2[k] += walker[k]*walker[k];
		
		if(iblock == (nblock-1))
		{
			double err_g=Error(vec_gofr[k],vec_gofr_2[k],iblock);
			Gave << r+0.5*bin_size <<" "<< vec_gofr[k]/(double)nblock <<" "<< err_g << endl;
		}
	}
	
	Epot << iblock << " " << epot_ave/(iblock+1) << " " << Error(epot_ave, epot_ave_sq, iblock+1) << endl;
	Etot << iblock << " " << etot_ave/(iblock+1) << " " << Error(etot_ave, etot_ave_sq, iblock+1) << endl;
	Ekin << iblock << " " << ekin_ave/(iblock+1) << " " << Error(ekin_ave, ekin_ave_sq, iblock+1) << endl;
	Temp << iblock << " " << temp_ave/(iblock+1) << " " << Error(temp_ave, temp_ave_sq, iblock+1) << endl;
	Press << iblock << " " << press_ave/(iblock+1) << " " << Error(press_ave, press_ave_sq, iblock+1) << endl;
	
	cout << "----------------------------" << endl << endl;
	
	Epot.close();
	Ekin.close();
	Etot.close();
	Temp.close();
	Press.close();
	Gave.close();
}


double Error(double sum, double sum2, int iblk){
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
