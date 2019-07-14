#include <iostream>
#include <cmath>
#include "cromosome.h"
#include "random.h"

Cromosome::Cromosome(){
	for(int i=0; i<m_lpath; i++){
		m_gene.push_back(i);
	}
}

Cromosome::~Cromosome(){
	
}

void Cromosome::Shuffle(Random& rnd, std::vector<City> vec_city){
	
	for(int i=0; i<100; i++){
		Inversion(rnd);
	}
	
	Check();
	
	Cost(vec_city);
	
	return;
}

void Cromosome::Exchange (int &a, int &b){
	int c;
	c=a;
	a=b;
	b=c;
	return;
}

bool Cromosome::Check(){
	for(int i=0; i<m_lpath-1; i++){
		for(int j=i+1; j<m_lpath; j++){
			if(m_gene[i]==m_gene[j]){
				std::cout << "Path contains a repetition" << std::endl;
				return false;
			}
		}
	}
	return true;
}

void Cromosome::Print(void){
	for(int i=0; i<m_lpath; i++){
		std::cout << m_gene[i] << " ";
	}
	std::cout << " ---> lenght: " << m_lenght << std::endl;
	return;
}

double Cromosome::GetLen() {
	return m_lenght;
}

std::vector<int> Cromosome::GetGene(void){
	return m_gene;
}

void Cromosome::Set(int index, int value){
	m_gene[index] = value;
	return;
}

int Cromosome::Get(int index){
	return m_gene[index];
}

Cromosome Cromosome::operator=(Cromosome crom2){
	m_lenght = crom2.GetLen();
	
	for(int i=0; i<m_lpath; i++){
		m_gene[i] = crom2.Get(i);
	}
	return *this;
}

void Cromosome::Cost(std::vector<City> vec_city){
	int city1, city2;
	m_lenght = 0.;
	
	for(unsigned int i=0; i<vec_city.size(); i++){
		
		
		city1 = m_gene[i];
		
		if(i==vec_city.size()-1){
			city2 = m_gene[0];
		} else{
			city2 = m_gene[i+1];
		}
		//std::cout<< i+1 <<std::endl;
		//std::cout<<city1 << " " << city2 <<std::endl;
		m_lenght += sqrt( pow(vec_city[city1].GetX() - vec_city[city2].GetX(), 2) + pow(vec_city[city1].GetY() - vec_city[city2].GetY(), 2) );

	}
	
	return;
}

void Cromosome::Inversion(Random& rnd){
	int pos1, pos2;
	
	pos1 = (int)rnd.Rannyu(0., m_lpath);
	pos2 = (int)rnd.Rannyu(0., m_lpath);
	while(pos1==pos2){
		pos2 = (int)rnd.Rannyu(0., m_lpath);
	}
		
	Exchange(m_gene[pos1], m_gene[pos2]);
	
	if(!Check())
		std::cout << "ERROR : INVERSION" << std::endl;
	
	return;
}

void Cromosome::Shift(Random& rnd){
	std::vector<int> tmp (m_gene);
	int index = (int)rnd.Rannyu(1., m_lpath);
	
	for(int i=0; i<m_lpath; i++){
		m_gene[index] = tmp[i];
		index++;
		if(index==m_lpath) index=0;
	}
	
	if(!Check())
		std::cout << "ERROR : SHIFT" << std::endl;

	return;
}

void Cromosome::FiniteShift(Random& rnd){

	std::vector<int> tmp (m_gene);
	int start = (int)rnd.Rannyu(0., m_lpath/2);
	int index = (int)rnd.Rannyu(0., m_lpath/2);
	
	for(int i=0; i<start; i++){
		m_gene[i] = tmp[i];
	}
	
	for(int i=start; i<m_lpath; i++){
		m_gene[index+start] = tmp[i];
		index++;
		if(index+start==m_lpath) index=0;
	}
	
	if(!Check())
		std::cout << "ERROR : FINITE SHIFT" << std::endl;
	
	return;
}

void Cromosome::Permutation(Random& rnd){
	
	std::vector<int> tmp (m_gene);
	
	int lenght = (int)rnd.Rannyu(0., m_lpath/2);
	int first = (int)rnd.Rannyu(0., m_lpath/2 - lenght);
	int second = (int)rnd.Rannyu(m_lpath/2, m_lpath - lenght);
	
	for(int i=0; i<first; i++)
		m_gene[i] = tmp[i];
	
	for(int i=0; i<lenght; i++)
		m_gene[first+i] = tmp[second+i];
	
	for(int i=0; i<second-first; i++)
		m_gene[first+lenght+i] = tmp[first+lenght+i];
	
	for(int i=0; i<lenght; i++)
		m_gene[second+i] = tmp[first+i];
	
	for(int i=0; i<m_lpath-second-lenght; i++)
		m_gene[second+lenght+i] = tmp[second+lenght+i];
	
	if(!Check())
		std::cout << "ERROR : PERMUTATION" << std::endl;
	
	return;
}

void Cromosome::Reverse(Random& rnd){
	
	std::vector<int> tmp (m_gene);
	
	int begin = (int)rnd.Rannyu(0., m_lpath-2);
	int lenght = (int)rnd.Rannyu(0., m_lpath-begin+1);
	
	for(int i=0; i<begin; i++)
		m_gene[i] = tmp[i];
	
	for(int i=0; i<lenght; i++)
		m_gene[begin+i] = tmp[begin+lenght-i-1];
	
	for(int i=0; i<m_lpath-begin-lenght; i++)
		m_gene[begin+lenght+i] = tmp[begin+lenght+i];
	
	if(!Check())
		std::cout << "ERROR : REVERSE" << std::endl;
	
	return;
}

void Cromosome::Crossover(Random& rnd, Cromosome& crom2){
	std::vector<int> tmp1 (m_gene);
	std::vector<int> tmp2 (crom2.GetGene());
	
	std::vector<int> new1 (m_gene);
	std::vector<int> new2 (crom2.GetGene());
	
	std::vector<int> missingfrom1;
	std::vector<int> missingfrom2;
	
	int divide = (int)rnd.Rannyu(0., m_lpath);
	
	for(int i=divide; i<m_lpath; i++){
		missingfrom1.push_back( m_gene[i] );
		missingfrom2.push_back( crom2.Get(i) );
		new1.pop_back();
		new2.pop_back();
	}
	
	for(int i=0; i<m_lpath; i++){
		for(unsigned int j=0; j<missingfrom1.size(); j++){
			if(tmp2[i] == missingfrom1[j]){
				new1.push_back(tmp2[i]);
				break;
			}
		}
	}
	
	for(int i=0; i<m_lpath; i++){
		for(unsigned int j=0; j<missingfrom2.size(); j++){
			if(tmp1[i] == missingfrom2[j]){
				new2.push_back(tmp1[i]);
				break;
			}
		}
	}
	
	for(int i=0; i<m_lpath; i++){
		m_gene[i] = new1[i];
		crom2.Set(i, new2[i]);
	}
	
	if(!Check() || !crom2.Check()){
		std::cout << "ERROR: CROSSOVER" << std::endl;
	}
	
	return;
}



