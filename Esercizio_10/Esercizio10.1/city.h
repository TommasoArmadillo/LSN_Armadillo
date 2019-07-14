#ifndef __City_h__
#define __City_h__



//classe City
class City {
public:
	//costruttori e distruttori
	City ();
	City (double, double ,int);
	~City ();
	
	//metodi
	void SetX(double);
	void SetY(double);
	void SetIndex(int);
	double GetX();
	double GetY();
	int GetIndex();
	
protected:
	double m_x, m_y;
	int m_index;
};





#endif //__City_h__


