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


    //parametri
    double S0 = 100.;
    double T = 1.;
    double K = 100.;
    double r = 0.1;
    double sig = 0.25;
    int Nint=100;
    int M = 10000;
    int N = 100;
    int L = M / N;

    double ST;  //prezzo dell'opzione al tempo richiesto
    double* C1 = new double [M] {0};
    double* P1 = new double [M] {0};
    //metodo diretto
    for (int i = 0;i < M;i++) {
        ST = S0 * exp((r - 0.5 * pow(sig, 2)) * T + sig * rnd.Gauss(0., T));    //implemento ST direttamente al tempo finale
        //implemento le opzioni call o put
        if (0 > (ST-K)) {
            C1[i] = 0;
            P1[i] = exp(-r * T) * (K - ST);
        }
        else {
            C1[i] = exp(-r * T) * (ST - K);
            P1[i] = 0;
        }       
    }

    //blocking method
    double* aveC = new double[N] {0};
    double* av2C = new double[N] {0};
    double* aveP = new double[N] {0};
    double* av2P = new double[N] {0};
    int* x = new int[N] {0};
    for (int i = 0;i < N;i++)
        x[i] = i + 1;
    for (int i = 0;i < N; i++) {
        double sumC = 0;
        double sumP = 0;
        for (int j = 0;j < L;j++) {
            int k = j + i * L;
            sumC = sumC + C1[k];
            sumP = sumP + P1[k];
        }
        aveC[i] = sumC / (double)L;
        av2C[i] = pow(aveC[i], 2);
        aveP[i] = sumP / (double)L;
        av2P[i] = pow(aveP[i], 2);
    }
    double* sum_progC = new double [N] {0};
    double* su2_progC = new double [N] {0};
    double* err_progC = new double [N] {0};
    double* sum_progP = new double [N] {0};
    double* su2_progP = new double [N] {0};
    double* err_progP = new double [N] {0};


    for (int i = 0;i < N;i++) {
        for (int j = 0;j < i + 1;j++) {
            sum_progC[i] = sum_progC[i] + aveC[j];
            su2_progC[i] = su2_progC[i] + av2C[j];
            sum_progP[i] = sum_progP[i] + aveP[j];
            su2_progP[i] = su2_progP[i] + av2P[j];
        }
        sum_progC[i] = sum_progC[i] / ((double)(i)+1);
        su2_progC[i] = su2_progC[i] / ((double)(i)+1);
        err_progC[i] = error(sum_progC, su2_progC, i);
        sum_progP[i] = sum_progP[i] / ((double)(i)+1);
        su2_progP[i] = su2_progP[i] / ((double)(i)+1);
        err_progP[i] = error(sum_progP, su2_progP, i);
    }

    ofstream data1("data_Direct.dat");
    for (int i = 0;i < N;i++)
        data1 << x[i] * L << "\t" << sum_progC[i] << " \t" << err_progC[i] << "\t" << sum_progP[i] << " \t" << err_progP[i] << endl;
    data1.close();


    double* C2 = new double [M] {0};
    double* P2 = new double [M] {0};

    //metodo standard
    for (int i = 0;i < M;i++) {
        for (int j = 0;j < Nint;j++) {
            //implemtento ST al tempo successivo per ogni tempo fino ad avere un valore finale per ST al tempo finale
            if (j == 0)
                ST = S0 * exp((r - 0.5 * pow(sig, 2)) / Nint + sig * rnd.Gauss(0., 1.) / sqrt(Nint));
            else
                ST = ST * exp((r - 0.5 * pow(sig, 2)) / Nint + sig * rnd.Gauss(0., 1.) / sqrt(Nint));
        }
        //implemento le opzioni call o put
        if (0 > (ST - K)) {
            C2[i] = 0;
            P2[i] = exp(-r * T) * (K - ST);
        }
        else {
            C2[i] = exp(-r * T) * (ST - K);
            P2[i] = 0;
        }
    }

    //resetto
    for (int i = 0;i < N; i++) {
        aveC[i] = 0;
        aveP[i] = 0;
        av2C[i] = 0;
        av2P[i] = 0;
        sum_progC[i] = 0;
        su2_progC[i] = 0;
        err_progC[i] = 0;
        sum_progP[i] = 0;
        su2_progP[i] = 0;
        err_progP[i] = 0;
    }
    //blocking method
    for (int i = 0;i < N; i++) {
        double sumC = 0;
        double sumP = 0;
        for (int j = 0;j < L;j++) {
            int k = j + i * L;
            sumC = sumC + C2[k];
            sumP = sumP + P2[k];
        }
        aveC[i] = sumC / (double)L;
        av2C[i] = pow(aveC[i], 2);
        aveP[i] = sumP / (double)L;
        av2P[i] = pow(aveP[i], 2);
    }

    for (int i = 0;i < N;i++) {
        for (int j = 0;j < i + 1;j++) {
            sum_progC[i] = sum_progC[i] + aveC[j];
            su2_progC[i] = su2_progC[i] + av2C[j];
            sum_progP[i] = sum_progP[i] + aveP[j];
            su2_progP[i] = su2_progP[i] + av2P[j];
        }
        sum_progC[i] = sum_progC[i] / ((double)(i)+1);
        su2_progC[i] = su2_progC[i] / ((double)(i)+1);
        err_progC[i] = error(sum_progC, su2_progC, i);
        sum_progP[i] = sum_progP[i] / ((double)(i)+1);
        su2_progP[i] = su2_progP[i] / ((double)(i)+1);
        err_progP[i] = error(sum_progP, su2_progP, i);
    }

    ofstream data2("data_Std.dat");
    for (int i = 0;i < N;i++)
        data2 << x[i] * L << "\t" << sum_progC[i] << " \t" << err_progC[i] << "\t" << sum_progP[i] << " \t" << err_progP[i] << endl;
    data2.close();
}    