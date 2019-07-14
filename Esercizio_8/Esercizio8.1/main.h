//global variables
Random rnd;
int seed[4];
int p1, p2;
int accepted, attempted;
std::vector<double> sum;
std::vector<double> sum2;

//parametres
double mu, sigma, delta;
double norm;
double x_old, x;
double p_old, p_new;
double ene, err_ene;
int n_block, n_step;
int passi_max = 1000;
double d_mu = 0.005, d_sigma=0.005, delta_H = 1.0, precision = 0.001, alpha = 0.01;
std::vector<double> grad (2, 0.);
int n_bin=100, n_points=100000;
double bin_size;
double extremum;

//functions
void Initialize(void);
void Present(void);
double Error(double, double, int);
double Potential(double);
void Metropolis(void);
double PsiSquared(double, double, double);
void Energy(int);
double HPsiOverPsi(double x);
void Average(int);
void Reset(void);
void Update(int);
void Gradient_Descent(void);
void Histogram(void);
void Output(void);
