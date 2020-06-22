#include <cstdlib>
#include <iostream>
#include <fstream>
#include <cmath>
#include<string>
#include"random.h"

using namespace std;

double error(double* AV, double* AV2, int n) {
	if (n == 0)
		return 0;
	else
		return sqrt((AV2[n]-pow(AV[n],2))/n);
}


int main() {
	int M = 10000;
	int N = 100;
	int L = M/N;



    //Generatore numeri random
    Random rnd;
    int seed[4];
    int p1, p2;
    ifstream Primes("Primes");
    if (Primes.is_open()) {
        Primes >> p1 >> p2;
    }
    else cerr << "PROBLEM: Unable to open Primes" << endl;
    Primes.close();

    ifstream input("seed.in");
    string property;
    if (input.is_open()) {
        while (!input.eof()) {
            input >> property;
            if (property == "RANDOMSEED") {
                input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
                rnd.SetRandom(seed, p1, p2);
            }
        }
        input.close();
    }
    //cout << seed[0] << seed[1] << seed[2] << seed[3] << endl;
    else cerr << "PROBLEM: Unable to open seed.in" << endl;

    rnd.SaveSeed();



    // creo il vettore di numeri casuali tra 0 e 1
    double* r=new double[M];
    for (int i = 0;i < M;i++)
        r[i] = rnd.Rannyu();

    double* x = new double[N];  //vettore di interi 
    for (int j = 0;j < N;j++)
        x[j] = j;
    double* ave = new double[N] {0};
    double* av2 = new double[N] {0};
    double* sum_prog = new double[N] {0};
    double* su2_prog = new double[N] {0};
    double* err_prog = new double[N] {0};


    //implementazione blocking method
    for (int i = 0;i < N; i++) {
        double sum = 0;
        for (int j = 0;j < L;j++) {
            int k = j + i * L;
            sum = sum+r[k];     //implemento la somma su cui effettuo il blocking method con i numeri casuali
        }
        ave[i] = sum/(double)L;
        av2[i] = pow(ave[i], 2);
    }



    for (int i = 0;i < N;i++) {       
        for (int j = 0;j < i+1;j++) {
            sum_prog[i] = sum_prog[i]+ave[j];
            su2_prog[i] = su2_prog[i]+av2[j];
        }
        sum_prog[i] =sum_prog[i]/((double)(i)+1);
        su2_prog[i] =su2_prog[i]/ ((double)(i)+1);
        err_prog[i] =error(sum_prog, su2_prog, i);
    }
    
    ofstream data1("data01_1.dat");
    for (int i = 0;i < N;i++)
        data1 << x[i] * L <<"\t"<< sum_prog[i] << " \t" << err_prog[i]<<endl;
    data1.close();
    
    
    //01.b
    
    //resetto le medie così da fare nuovamente il blocking method
    for (int i = 0;i < N;i++) {
        ave[i] = 0;
        av2[i] = 0;
        sum_prog[i] = 0;
        su2_prog[i] = 0;
        err_prog[i] = 0;
    }

    //blocking method
    for (int i = 0;i < N;i++) {
        double sum2 = 0;
        for (int j = 0;j < L;j++) {
            int k = j + i * L;
            double funz = pow((r[k] - 0.5),2);  //cambio la funzione con cui implemento il blocking method
            sum2 = sum2+funz;
        }
        ave[i] = sum2 / (double)L;
        av2[i] = pow(ave[i], 2);
    }



    for (int i = 0;i < N;i++) {
        for (int j = 0;j < i + 1;j++) {
            sum_prog[i] = sum_prog[i] + ave[j];
            su2_prog[i] = su2_prog[i] + av2[j];
        }
        sum_prog[i] = sum_prog[i] / ((double)(i)+1);
        su2_prog[i] = su2_prog[i] / ((double)(i)+1);
        err_prog[i] = error(sum_prog, su2_prog, i);
    }

    ofstream data2("data01_2.dat");
    for (int i = 0;i < N;i++)
        data2 << x[i] * L << "\t" << sum_prog[i] << " \t" << err_prog[i] << endl;
    data2.close();

    //01.c
    int step;
    int k = 0;
    int Ncicli = 100;
    int Nstep = 100;
    double* chi = new double[Nstep]{ 0 };
    double* conteggi = new double[Nstep]{ 0 };
    do {
        for (int y = 0;y < Nstep;y++)
            conteggi[y] = 0;        //definisco un vettore di conteggi che azzero all'inizio di ogni ciclo
        for (int i = 0;i < M;i++) {
            r[i] = rnd.Rannyu();    //effettuo l'esperimento
            step = 0;
            //implemento il vettore conteggi
            while ((double)step / 100. < 1) {   
                if (r[i] > (double)step / 100. && r[i] < ((double)step + 1.) / 100.) {
                    conteggi[step]++;
                    break;
                }
                else
                    step++;
            }
        }
        for (int i = 0;i < 100;i++)
            chi[k] = chi[k] + pow((conteggi[i] - 100), 2) / 100;    //implemento il vettore di chi2 con la formula
        k++;
    } while (k < Ncicli);

    ofstream data3("data01_3.dat");
    for (int i = 0;i < N;i++)
        data3 << x[i]+1 << "\t" << chi[i] << endl;
    data3.close();
	
	return 0;
}