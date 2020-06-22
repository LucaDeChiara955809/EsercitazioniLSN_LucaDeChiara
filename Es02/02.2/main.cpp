#define _USE_MATH_DEFINES

#include <cstdlib>
#include <iostream>
#include <fstream>
#include <cmath>
#include<string>
#include"random.h"
#include"Posizione.h"

using namespace std;

double error(double AV, double AV2, int n) {
    if (n == 0)
        return 0;
    else
        return sqrt((AV2 - pow(AV, 2)) / n);
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
    //fine generatore random
 

    int M = 10000;
    int Nstep = 0;
    int N = 100;
    double* ave=new double [N];
    ofstream data1("dataRWreticolo.dat");
    ofstream data2("dataRWsfera.dat");
    double* av2 = new double[N] {0};
    double sum_prog = 0;
    double su2_prog = 0;
    double err_prog = 0;
    int L = M / N;
    int dado;
    double* posFin = new double [M] {0};
    //poichè devo registrare la posizione ad ogni singolo step effettuo 10^4 esperimenti per ogni step di modo da avere una posizione finale mediata per ogni step
    do {
        for (int j = 0;j < M;j++) {
            Posizione X;
            //per il reticolo cubico con a=1 mi sposto con probabilità dettate da un lancio di un dado
            for (int i = 0;i < Nstep;i++) {
                dado = floor(rnd.Rannyu(1, 7));
                if (dado == 1)
                    X.SetX(X.GetX() + 1);
                else if (dado == 2)
                    X.SetX(X.GetX() - 1);
                else if (dado == 3)
                    X.SetY(X.GetY() + 1);
                else if (dado == 4)
                    X.SetY(X.GetY() - 1);
                else if (dado == 5)
                    X.SetZ(X.GetZ() + 1);
                else if (dado == 6)
                    X.SetZ(X.GetZ() - 1);
            }
            posFin[j] = X.GetDist(); //registro la posizione finale
        }
        //resetto i parametri per il blocking method
        for (int i = 0;i < N; i++) {
            ave[i] = 0;
            av2[i] = 0;
            sum_prog = 0;
            su2_prog = 0;
            err_prog = 0;
        }
        for (int i = 0;i < N; i++) {
            double sum = 0;
            for (int j = 0;j < L;j++) {
                int k = j + i * L;
                sum = sum + pow(posFin[k], 2);  //carico il blocking method con il quadrato del raggio ottenuto
            }
            ave[i] = sum / (double)L;
            ave[i] = sqrt(ave[i]);
            av2[i] = pow(ave[i], 2);
        }
        for (int i = 0;i < N;i++) {
                sum_prog = sum_prog + ave[i];
                su2_prog = su2_prog + av2[i];
        }
        sum_prog = sum_prog / (double)N;
        su2_prog = su2_prog / (double)N;
        err_prog = error(sum_prog, su2_prog, N);
        data1 <<Nstep<<"\t"<<sum_prog<<"\t"<<err_prog << endl;
        Nstep++;    //aumento di uno step
    } while (Nstep <= 100);
    data1.close();

    Nstep = 0;  //resetto gli step
    double* posFin2 = new double [M] {0}; 
    do {
        for (int j = 0;j < M;j++) {
            Posizione X;
            for (int i = 0;i < Nstep;i++) {
                double th = rnd.Theta();    //implemento un theta casuale implementato col metodo accept reject nella classe Random
                double ph = rnd.Phi();      //implemento uno phi casuale implementato col metodo accept reject nella classe Random
                X.SetX(X.GetX() + sin(ph) * cos(th));
                X.SetY(X.GetY() + sin(ph) * sin(th));
                X.SetZ(X.GetZ() + cos(ph));
            }
            posFin2[j] = X.GetDist();
        }
        //resetto i parametri per il blocking method
        for (int i = 0;i < N; i++) {
            ave[i] = 0;
            av2[i] = 0;
            sum_prog = 0;
            su2_prog = 0;
            err_prog = 0;
        }
        for (int i = 0;i < N; i++) {
            double sum = 0;
            for (int j = 0;j < L;j++) {
                int k = i + j * L;
                sum = sum + pow(posFin2[k], 2); //carico il blocking method con il quadrato del raggio ottenuto
            }
            ave[i] = sum / (double)L;
            ave[i] = sqrt(ave[i]);
            av2[i] = pow(ave[i], 2);
        }
        for (int i = 0;i < N;i++) {
            sum_prog = sum_prog + ave[i];
            su2_prog = su2_prog + av2[i];
        }
        sum_prog = sum_prog / (double)N;
        su2_prog = su2_prog / (double)N;
        err_prog = error(sum_prog, su2_prog, N);
        data2 << Nstep << "\t" << sum_prog << "\t" << err_prog << endl;
        Nstep++;

    } while (Nstep <= 100);
    data2.close();

    return 0;
}
