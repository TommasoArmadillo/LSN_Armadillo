#ifndef __Cromosome_h__
#define __Cromosome_h__

#include <vector>
#include "city.h"
#include "random.h"

//classe Cromosome
class Cromosome {
public:
	//costruttori e distruttori
	Cromosome ();
	~Cromosome ();
	
	//metodi
	void Shuffle(Random&, int, std::vector<City>);
	void Exchange(int&, int&);
	bool Check();
	void Print(void);
	double GetLen(void);
	std::vector<int> GetGene(void);
	Cromosome operator=(Cromosome);
	void Set(int, int);
	int Get(int);
	void Cost(std::vector<City>);
	void Inversion(Random&);
	void Shift(Random&);
	void FiniteShift(Random&);
	void Permutation(Random&);
	void Reverse(Random&);
	void Crossover(Random&, Cromosome&);
	
protected:
	double m_lenght=0.;
	int m_lpath=30;
	std::vector<int> m_gene;
};

#endif /* __Cromosome_h__ */
