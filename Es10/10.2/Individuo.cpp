#define _USE_MATH_DEFINES

#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include "Individuo.h"
#include"random.h"

using namespace std;

Individuo::Individuo() {
	m_nome = new int[32];
	for (int i = 0;i < 32;i++)
		m_nome[i] = i + 1;
}

Individuo::Individuo(int l) {
	m_nome = new int[l];
	for (int i = 0;i < l;i++)
		m_nome[i] = i + 1;
}

Individuo::Individuo(int l,int* nomi) {
	m_nome = new int[l];
	for (int i = 0;i < l;i++)
		m_nome[i] = nomi[i];
}

Individuo :: ~Individuo() {}

int* Individuo::GetSequence() {
	return m_nome;
}

int Individuo::GetCity(int o) {
	return m_nome[o];
}

void Individuo::SetSequence(int* nomi,int Lunghezza_array) {
	for (int i = 0;i < Lunghezza_array;i++)
		m_nome[i] = nomi[i];
	return;
}

void Individuo::SetCity(int posizione, int città) {
	m_nome[posizione] = città;
	return;
}