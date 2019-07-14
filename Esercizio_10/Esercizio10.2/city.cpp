#include "city.h"


//costruttori e distruttore
City::City(){
	m_x = 0.;
	m_y = 0.;
	m_index = 0;
}

City::City(double x, double y, int index){
	m_x = x;
	m_y = y;
	m_index = index;
}

City::~City(){
	
}

//Set methods
void City::SetX(double x){
	m_x = x;
	return;
}

void City::SetY(double y){
	m_y = y;
	return;
}

void City::SetIndex(int index){
	m_index = index;
	return;
}

//Get methods
double City::GetX(){
	return m_x;
}

double City::GetY(){
	return m_y;
}

int City::GetIndex(){
	return m_index;
}








