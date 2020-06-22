#define _USE_MATH_DEFINES

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
        return sqrt((AV2[n] - pow(AV[n], 2)) / n);
}

double f1(double x) {
    return M_PI / 2 * cos(x * M_PI / 2);    //funzione usata per il calcolo dell'integrale con distribuzione uniforme
}

double f2(double x) {
    return M_PI / 2 * cos(x * M_PI / 2)/(2*(1-x)); //funzione usata per il calcolo dell'integrale con importance sampling (la distribuzione usata è 2*(1-x))
}


int main() {
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


    //blocking method
    int M = 100000;
    int N = 100;
    int L = M / N;
    double* x = new double[N];
    for (int j = 0;j < N;j++)
        x[j] = j;
    double* ave = new double[N] {0};
    double* av2 = new double[N] {0};
    double* sum_prog = new double[N] {0};
    double* su2_prog = new double[N] {0};
    double* err_prog = new double[N] {0};

    for (int i = 0;i < N; i++) {
        double sum = 0;
        for (int j = 0;j < L;j++)
            sum = sum + f1(rnd.Rannyu());//implemento le medie a blocchi con la funzione integranda (distribuzione uniforme)
        ave[i] = sum / (double)L;
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

    ofstream data1("data_unif.dat");
    for (int i = 0;i < N;i++)
        data1 << x[i] * L << "\t" << sum_prog[i] << " \t" << err_prog[i] << endl;
    data1.close();


    //resetto
    for (int i = 0;i < N; i++) {
        ave[i] = 0;
        av2[i] = 0;
        sum_prog[i] = 0;
        su2_prog[i] = 0;
        err_prog[i] = 0;
    }

    //importance sampling
    for (int i = 0;i < N; i++) {
        double sum = 0;
        for (int j = 0;j < L;j++)
            sum = sum + f2(rnd.funz());//implemento le medie a blocchi con la funzione ottenuta dall'importance sampling [2(1-x)]
        ave[i] = sum / (double)L;
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

    ofstream data2("data_importance.dat");
    for (int i = 0;i < N;i++)
        data2 << x[i] * L << "\t" << sum_prog[i] << " \t" << err_prog[i] << endl;
    data2.close();
    return 0;
}