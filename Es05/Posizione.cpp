#define _USE_MATH_DEFINES

#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include "Posizione.h"

using namespace std;

Posizione::Posizione() {
	m_x = 0;
	m_y = 0;
	m_z = 0;
}

Posizione::Posizione(double x,double y,double z) {
	m_x = x;
	m_y = y;
	m_z = z;
}

Posizione :: ~Posizione() {}

void Posizione::SetX(double x) {
	m_x = x;
	return;
}

void Posizione::SetY(double x) {
	m_y = x;
	return;
}

void Posizione::SetZ(double x) {
	m_z = x;
	return;
}

double Posizione::GetX() {
	return m_x;
}

double Posizione::GetY() {
	return m_y;
}

double Posizione::GetZ() {
	return m_z;
}

double Posizione::GetDist() {
	return sqrt(pow(m_x, 2) + pow(m_y, 2) + pow(m_z, 2));
}