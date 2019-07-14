//data
int n_city=30;
int n_path=900;
int n_generation = 900;
int ind_mut;
double r=1.;
Random rnd;
int seed[4];
int p1, p2;

//Settings
std::vector<City> Cities_vec (n_city);
std::vector<Cromosome> Population (n_path);

//probabilities of mutations
double prob_inversion = 0.1;
double prob_shift = 0.1;
double prob_finiteshift = 0.1;
double prob_permutation = 0.1;
double prob_reverse = 0.1;
double prob_crossover = 0.55;

//functions
void GeneratePopulation();
void Initialize();
void Order(int, int);
int Choose();
void Measure(int);
void Mutate(Cromosome&);
void Output();
void NextGeneration();
