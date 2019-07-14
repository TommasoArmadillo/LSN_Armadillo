//data
int n_city=30;
int ind_mut;
double r=1.;
std::vector<double> temp_vec;
int Temp_num_step;
double beta;
double temp;
int n_step = 15000;
int attempted, accepted;

Random rnd;
int seed[4];
int p1, p2;

//Settings
std::vector<City> Cities_vec (n_city);
Cromosome Path;
Cromosome Best;

//functions
void GenerateStartingPath();
void Initialize();
void Metropolis();
void Mutate(Cromosome&);
void Output();
void Measure(int);
