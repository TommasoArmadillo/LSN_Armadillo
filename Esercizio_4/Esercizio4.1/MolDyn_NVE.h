/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
//parameters, observables
const int m_props=1000;
int n_props;
int iv,ik,it,ie;
double stima_pot, stima_kin, stima_etot, stima_temp, stima_press;
double igofr, bin_size, nbins;
// averages
double acc,att;

//configuration
const int m_part=108;
double x[m_part],y[m_part],z[m_part],xold[m_part],yold[m_part],zold[m_part];
double vx[m_part],vy[m_part],vz[m_part];

// thermodynamical state
int npart;
double energy,temp,vol,rho,box,rcut;

// simulation
int nstep, nblock, lblock, iprint, seed;
double delta;
bool cont, rescale;

//functions
void Input(Random);
void Move(void);
void ConfFinal(void);
void ConfXYZ(int);
void Measure(void);
double Force(int, int);
double Pbc(double);
void Accumulate(int);
void Averages(int);
double Error(double, double, int);

//average block
std::vector<double> vec_epot;				//vector to store mean and mean of squares of each block
std::vector<double> vec_epot_sq;
std::vector<double> vec_etot;
std::vector<double> vec_etot_sq;
std::vector<double> vec_ekin;
std::vector<double> vec_ekin_sq;
std::vector<double> vec_temp;
std::vector<double> vec_temp_sq;
std::vector<double> vec_press;
std::vector<double> vec_press_sq;

std::vector<double> walker (m_props, 0.);
std::vector<double> vec_gofr;
std::vector<double> vec_gofr_2;
/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
