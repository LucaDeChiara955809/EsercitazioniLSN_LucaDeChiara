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
double Lmin = 200;

int metodo = 1;
const int Ntimes = 10000;
double prob;

double r_mean[Ntimes + 1];
double r_min[Ntimes + 1];
double r_devSt [Ntimes + 1];

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
//medie
double mean(Individuo*);
double DevSt(Individuo*);
#pragma once
