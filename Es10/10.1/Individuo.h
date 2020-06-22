#ifndef __Individuo__
#define __Individuo__

class Individuo {
private:

protected:
	int* m_nome;

public:
	//constructor
	Individuo();
	Individuo(int l);
	Individuo(int l,int* nomi);
	//destructor
	~Individuo();
	//metodhs
	int* GetSequence();
	int GetCity(int);
	void SetSequence(int*,int);
	void SetCity(int, int);
};
#endif
#pragma once
