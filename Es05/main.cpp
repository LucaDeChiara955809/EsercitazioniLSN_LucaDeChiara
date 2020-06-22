#define _USE_MATH_DEFINES

#include <cstdlib>
#include <iostream>
#include <fstream>
#include <cmath>
#include<string>
#include"random.h"
#include"Posizione.h"

using namespace std;

double error(double* AV, double* AV2, int n) {
    if (n == 0)
        return 0;
    else
        return sqrt((AV2[n] - pow(AV[n], 2)) / n);
}

double onda1(Posizione x) {
    return exp(-2.*x.GetDist());
}

double onda2(Posizione x) {
    return pow(x.GetDist(),2) * exp(-x.GetDist())*pow(cos(atan(x.GetY()/x.GetX())),2);
}

int main() {
    //Inizializzazione random
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


    int M = 1000000;
    int N = 1000;
    int L = M / N;
    Posizione O(5.,5.,5.);//settiamo il punto di partenza di campionamento per la prima funzione (n=1 ,l=0 , m=0)
    int start=1;    //0 se campionamento uniforme   1 se campionamento gaussiano


    double RMax;
    if (start == 1)
        RMax = 0.75;    //parametri calcolati a seguito di ottimizzazione
    else
        RMax = 1.23;
    int accept = 0;
    int reject = 0;
    Posizione P;
    double* r = new double[M];
    double sum=0;
    double A,rr;
    for (int i = 0;i < M;i++){
        //campiono un punto preso casualmente da col metodo richiesto
        if (start == 1) {
            P.SetX(rnd.Gauss(O.GetX(), RMax));
            P.SetY(rnd.Gauss(O.GetY(), RMax));
            P.SetZ(rnd.Gauss(O.GetZ(), RMax));
        }
        else {
            P.SetX(O.GetX() + rnd.Rannyu(-RMax, RMax));
            P.SetY(O.GetY() + rnd.Rannyu(-RMax, RMax));
            P.SetZ(O.GetZ() + rnd.Rannyu(-RMax, RMax));
        }
        A = onda1(P) / onda1(O);    //definisco la funzione di accettazione da usare nell'algoritmo di Metropolis
        if (A>1)
            A = 1;
        rr = rnd.Rannyu();
        if (rr < A) {
            O = P;  //se accetto la mossa sostituisco il punto iniziale sennò no
            accept++;
        }
        else
            reject++;
        r[i] = O.GetDist(); //carico il vettore delle distanze
    }

    //blocking method
    double* ave = new double[N] {0};
    double* av2 = new double[N] {0};
    double* sum_prog = new double[N] {0};
    double* su2_prog = new double[N] {0};
    double* err_prog = new double[N] {0};
    for (int i = 0;i < N; i++) {
        sum = 0;
        for (int j = 0;j < L;j++) {
            int k = j + i * L;
            sum = sum + r[k];
        }
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

    double rapporto;
    rapporto = (double)accept / ((double)reject+(double)accept);
    if (start == 0)
        cout << "Campionamento effettuato con funzione uniforme:" << endl;
    else
        cout << "Campionamento effettuato con funzione gaussiana:" << endl;
    cout << endl;
    cout << "Funzione d'onda n=1, l=0, m=0" << endl;
    cout << "Il rate di accettazione è "<<rapporto << endl;
    cout << endl;
    double* x = new double[N];
    for (int i = 0;i < N;i++)
        x[i] = ((double)i + 1) * (double)L; //definisco un vettore x del numero di esperimenti
    ofstream data1;
    if (start == 0)
        data1.open("data_ground_uniforme.dat");
    else
        data1.open("data_ground_Gaussiano.dat");
    for (int i = 0;i < N;i++)
        data1 << x[i]<<"\t"<<sum_prog[i] << " \t" << err_prog[i]<< endl;
    data1.close();





    Posizione O2(10.,10.,10.);//settiamo il punto di partenza di campionamento per la seconda funzione (n=2 ,l=1 , m=0)
    double RMax2;
    if (start == 1)
        RMax2 = 2.;     //parametri calcolati a seguito di ottimizzazione
    else
        RMax2 = 3.15;
    double accept2 = 0;
    double reject2 = 0;
    double* r2 = new double[M];
    sum = 0;
    for (int i = 0;i < M;i++) {
        //campiono un punto preso casualmente da col metodo richiesto
        if (start == 1) {
            P.SetX(rnd.Gauss(O2.GetX(), RMax2));
            P.SetY(rnd.Gauss(O2.GetY(), RMax2));
            P.SetZ(rnd.Gauss(O2.GetZ(), RMax2));
        }
        else {
            P.SetX(O2.GetX() + rnd.Rannyu(-RMax2, RMax2));
            P.SetY(O2.GetY() + rnd.Rannyu(-RMax2, RMax2));
            P.SetZ(O2.GetZ() + rnd.Rannyu(-RMax2, RMax2));
        }
        A = onda2(P) / onda2(O2);      //definisco la funzione di accettazione da usare nell'algoritmo di Metropolis
        if (A > 1)
            A = 1;
        rr = rnd.Rannyu();
        if (rr < A) {
            O2 = P;     //se accetto la mossa sostituisco il punto iniziale sennò no
            accept2++;
        }
        else
            reject2++;
        r2[i] = O2.GetDist();   //carico il vettore delle distanze
    }

    //resetto parametri blocking method
    for (int i = 0;i < N; i++) {
        ave[i] = 0;
        av2[i] = 0;
        sum_prog[i] = 0;
        su2_prog[i] = 0;
        err_prog[i] = 0;
    }
    //blocking method
    for (int i = 0;i < N; i++) {
        sum = 0;
        for (int j = 0;j < L;j++) {
            int k = j + i * L;
            sum = sum + r2[k];
        }
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

    rapporto = (double)accept2 / ((double)reject2+(double)accept2);
    cout << "Funzione d'onda n=2, l=1, m=0" << endl;
    cout << "Il rate di accettazione è " << rapporto << endl;
    cout << endl;
    ofstream data2;
    if (start == 0)
        data2.open("data_ecc_uniforme.dat");
    else
        data2.open("data_ecc_Gaussiano.dat");
    for (int i = 0;i < N;i++)
        data2 << x[i] << "\t" << sum_prog[i] << " \t" << err_prog[i] << endl;
    data2.close();
}