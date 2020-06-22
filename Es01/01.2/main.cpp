#include <iostream>
#include<cstdlib>
#include <fstream>
#include <cmath>
#include<string>
#include"random.h"

using namespace std;


int main(int argc, char* argv[]) {

    int M = 10000;
    double lambda = 1;
    double gamma = 1;
    double mean = 0;
    double* s1 = new double[M] {0};
    double* s2 = new double[M] {0};
    double* s10 = new double[M] {0};
    double* s100 = new double[M] {0};

    //generatore numeri random
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
    else cerr << "PROBLEM: Unable to open seed.in" << endl;
    rnd.SaveSeed();

    //implemento l'istogramma del dado utilizzando una variante di Rannyu in cui 
    //genererò un numero casuale da 1 a 7 approssimato all'intero più basso
    for (int i = 0; i < M; i++) {
        s1[i] = floor(rnd.Rannyu(1,7));
        for (int j = 0;j < 2;j++) 
            s2[i] = s2[i] + floor(rnd.Rannyu(1, 7));
        for (int j = 0;j < 10;j++)
            s10[i] = s10[i] + floor(rnd.Rannyu(1, 7));
        for (int j = 0;j < 100;j++)
            s100[i] = s100[i] + floor(rnd.Rannyu(1, 7));
        s2[i] = s2[i] / 2;
        s10[i] = s10[i] / 10;
        s100[i] = s100[i] / 100;
    }

    ofstream data1("dataDado.dat");
    for (int i = 0;i < M;i++)
        data1 <<s1[i]  <<"\t"<<s2[i]<<"\t"<<s10[i] << "\t" << s100[i]<< endl;
    data1.close();

    //resetto i vettori
    for (int i = 0; i < M; i++) {
        s1[i] = 0;
        s2[i] = 0;
        s10[i] = 0;
        s100[i] = 0;
    }

    //implemento l'istogramma della laurentziana con numeri casuali ottenuti con il metodo implementato nel generatore casuale
    for (int i = 0; i < M; i++) {
        s1[i] = rnd.Lorentz(gamma,mean);
        for (int j = 0;j < 2;j++)
            s2[i] = s2[i] + rnd.Lorentz(gamma, mean);
        for (int j = 0;j < 10;j++)
            s10[i] = s10[i] + rnd.Lorentz(gamma, mean);
        for (int j = 0;j < 100;j++)
            s100[i] = s100[i] + rnd.Lorentz(gamma, mean);
        s2[i] = s2[i] / 2;
        s10[i] = s10[i] / 10;
        s100[i] = s100[i] / 100;
    }
    ofstream data3("dataLorentz.dat");
    for (int i = 0;i < M;i++)
        data3 << s1[i] << "\t" << s2[i] << "\t" << s10[i] << "\t" << s100[i] << endl;
    data3.close();
     
    //resetto
    for (int i = 0; i < M; i++) {
        s1[i] = 0;
        s2[i] = 0;
        s10[i] = 0;
        s100[i] = 0;
    }

    //implemento l'istogramma dell'esponenziale con numeri casuali ottenuti con il metodo implementato nel generatore casuale

    for (int i = 0; i < M; i++) {
         s1[i] = rnd.Exp(lambda);
         for (int j = 0;j < 2;j++)
             s2[i] = s2[i] + rnd.Exp(lambda);
         for (int j = 0;j < 10;j++)
             s10[i] = s10[i] + rnd.Exp(lambda);
         for (int j = 0;j < 100;j++)
             s100[i] = s100[i] + rnd.Exp(lambda);
         s2[i] = s2[i] / 2;
         s10[i] = s10[i] / 10;
         s100[i] = s100[i] / 100;
    }
  
    ofstream data2("dataExp.dat");
    for (int i = 0;i < M;i++)
        data2 << s1[i] << "\t" << s2[i] << "\t" << s10[i] << "\t" << s100[i] << endl;
    data2.close();

    return 0;
}