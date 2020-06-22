Random rnd;
int seed[4];
int p1;
int p2;

const int Ncities = 32;
const int population = 200;
double positions[Ncities][2];
int controllo;
Individuo padre;
Individuo madre;
Individuo best;
Individuo g1;
Individuo g2;
double Lmin = 200;
double TIn = 50;
double betaIn = 1 / TIn;
double beta;
double A1, A2;
double prob;
int indice;

int metodo = 1;
const int Ntimes = 10000;

double r_min[Ntimes + 1];

void RandomInitializer();
void PositionGenerator(int);
void RandomPopulationGenerator(Individuo*);
int CheckFunction(Individuo);
double L(Individuo);
Individuo ParentSelector(Individuo*);
//mutazioni
void SingleMutation(Individuo);
void ShiftMutation(Individuo);
void InversionMutation(Individuo);
void PermutationMutation(Individuo);	
void CrossingOver(Individuo, Individuo);
//ordinamento
void Quicksort(Individuo*, int, int);
int PBC(int);
int PBC0(int);

double Acceptance(Individuo, Individuo);
#pragma once
