#define _USE_MATH_DEFINES
#include <iostream>
#include <fstream>
#include <string>
#include "random.h"

using namespace std;

double error(double* AV, double* AV2, int n) {
    if (n == 0)
        return 0;
    else
        return sqrt((AV2[n] - pow(AV[n], 2)) / n);
}

int main(int argc, char* argv[]) {

    //generatore numeri casuali
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

    //simulazione esperimento di Buffon

    double Lu = 1. / 15.;   //definisco la lunghezza del segmento (NOTA   Lu<1)
    double d = 1. / 10.;    //definisco la spaziatura tra le linee verticali (NOTA  1>d>Lu)
    int M = 10000;      //Numero ripetizioni dell'esperimento
    int Nthrow = 10000; //Numero tiri
    int Hit;    
    double X1,theta,DX,X2;
    int j;
    double* pi = new double [M];
    for (int k = 0;k < M;k++){
        Hit = 0;
        for (int i = 0;i < Nthrow;i++) {
            X1 = rnd.Rannyu();      //Definisco un punto di partenza casuale all'interno del mio quadrato 1x1
            theta = rnd.Rannyu(0, 2*M_PI);  //definisco un angolo
            DX = Lu*cos(theta);     //effettuo la proiezione del segmento sull'asse x
            X2 = X1 + DX;
            j = 0;
            while (j * d <= 1) {
                if (X1 <= j * d && X2>=j * d)       //effettuo il test di verifica se colpisce o meno le righe verticali
                    Hit++;
                else if (X1>= j * d && X2<= j * d)
                    Hit++;
                j++;
            }
        }
        pi[k] = 2 *Lu * Nthrow / (Hit * d);     //carico tutto su un vettore pi che userò per il data blocking
    }

    //Data blocking
    int N = 100;
    int L = M / N;
    double* x = new double[N];
    for (int i = 0;i < N;i++)
        x[i] = L * i;
    double* ave = new double[N] {0};
    double* av2 = new double[N] {0};
    double* sum_prog = new double[N] {0};
    double* su2_prog = new double[N] {0};
    double* err_prog = new double[N] {0};
    double sum;
    for (int i = 0;i < N; i++) {
        sum = 0;
        for (int j = 0;j < L;j++) {
            int k = j + i * L;
            sum = sum + pi[k];
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

    ofstream data("dataBuffon.dat");
    for (int i = 0;i < N;i++)
        data  << x[i]<<"\t"<<sum_prog[i] << " \t" << err_prog[i] << endl;
    data.close();

    return 0;
}