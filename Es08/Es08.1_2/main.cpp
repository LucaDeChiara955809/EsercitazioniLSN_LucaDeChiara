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

double onda(double x, double mu, double sigma) {
    return exp(-pow(x-mu,2)/(2*pow(sigma,2)))+ exp(-pow(x + mu, 2) / (2 * pow(sigma, 2)));
}

double V(double x) {
    return pow(x, 4) - pow(x, 2) * 5 / 2;
}

double HOnda(double x, double mu, double sigma) {
    double e1 = exp(-pow(x - mu, 2) / (2 * pow(sigma,2)));
    double e2 = exp(-pow(x + mu, 2) / (2 * pow(sigma, 2)));
    return 0.5*(e1 / pow(sigma,2) - pow(x - mu, 2) / pow(sigma, 4)*e1 + e2 / pow(sigma,2) - pow(x + mu, 2) / pow(sigma, 4)*e2) + V(x) * onda(x, mu, sigma);
}


// Es 08.1 (codice in cui mi stampo) e 08.2 (ottimizzazione con mu e sigma)
int main() {
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
	int M = 100000;		//# steps in metropolis
	int N = 100;		//# blocks
	int L = int(M / N);	//# elements in each block

	double x_old;
	double x_new;
	double xmax = 2.8;
	double A, rr;
	double accept = 0;
	double reject = 0;
	double conf = 0;

	double Emin=100;
	double mubest,sigbest;

	double mu_i = 0.45;
	double sigma_i = 0.5;

	ofstream data("data2.dat");

	for (int i = 0; i < 20; i++) {
		double mu = mu_i + i * 0.0275;
		for (int j = 0; j < 20; j++) {
			double sigma = sigma_i + j * 0.025;
			x_old = 2.5;
			double* r = new double[M] {0};
			double* ave = new double[N] {0};
			double* av2 = new double[N] {0};
			double* sum_prog = new double[N] {0};
			double* su2_prog = new double[N] {0};
			double* err_prog = new double[N] {0};
			for (int j = 0; j < M; j++) {
				x_new = x_old + rnd.Rannyu(-xmax, xmax);
				A = pow(onda(x_new,mu,sigma), 2) / pow(onda(x_old,mu,sigma), 2);
				if (A > 1)
					A = 1;
				rr = rnd.Rannyu();
				if (rr <= A) {
					x_old = x_new;
					accept++;
				}
				else {
					reject++;
				}
				r[j] = x_old;
			}

			conf = accept / (reject+accept);
			cout << "The ratio between accepted and rejected jumps is " << conf << endl;

			for (int i = 0; i < N; i++) {
				double sum = 0;
				for (int j = 0; j < L; j++) {
					int k = j + i * L;
					sum = sum + HOnda(r[k], mu, sigma) /onda(r[k], mu, sigma);
				}
				ave[i] = sum / L;
				av2[i] = (ave[i]) * (ave[i]);
			}

			for (int i = 0; i < N; i++) {
				for (int j = 0; j < i + 1; j++) {
					sum_prog[i] = sum_prog[i] + ave[j];
					su2_prog[i] = su2_prog[i] + av2[j];
				}
				sum_prog[i] = sum_prog[i] / (i + 1); // Cumulative average
				su2_prog[i] = su2_prog[i] / (i + 1); // Cumulative square average
				err_prog[i] = error(sum_prog, su2_prog, i); // Statistical uncertainty
			}
			if (Emin > sum_prog[N - 1]) {
				Emin = sum_prog[N - 1];
				mubest = mu;
				sigbest = sigma;
			}
			data << mu << "\t" << sigma << "\t" << sum_prog[N - 1] << "\t" << err_prog[N - 1] << "\t"<<conf<< endl;
			accept = 0;
			reject = 0;
		}
	}
	cout << "Best mu:" << mubest << " best sigma:" << sigbest << endl;
	data.close();

    double xold = 2.5;
    double mu = mubest;
    double sigma = sigbest;
    double xMax = 2.8;
    accept = 0;
    reject = 0;
    double xnew;
    double* r = new double[M];
    double sum = 0;
    
    ofstream data2;
    data2.open("config.dat");
    for (int i = 0;i < M;i++) {
        xnew = xold + rnd.Rannyu(-xMax, xMax);
        A = pow(abs(onda(xnew, mu, sigma)), 2) / pow(onda(xold, mu, sigma), 2);
        if (A > 1)
            A = 1;
        rr = rnd.Rannyu();
        if (rr < A) {
            xold = xnew;
            accept++;
        }
        else
            reject++;
        r[i] = xold;
        data2 << r[i] << endl;
    }
    data2.close();
    double rapporto;
    rapporto = (double)accept / ((double)reject + (double)accept);
    cout <<"Rate di accettazione: "<< rapporto << endl << endl;


    double* ave = new double[N] {0};
    double* av2 = new double[N] {0};
    double* sum_prog = new double[N] {0};
    double* su2_prog = new double[N] {0};
    double* err_prog = new double[N] {0};
    for (int i = 0;i < N; i++) {
        sum = 0;
        for (int j = 0;j < L;j++) {
            int k = j + i * L;
            sum = sum + HOnda(r[k], mu, sigma) / onda(r[k], mu, sigma);
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
    int* x = new int[N];
    for (int i = 0;i < N;i++)
        x[i] = (i + 1) * L;

    ofstream data1("data1.dat");
    for (int i = 0;i < N;i++)
        data1 << x[i] << "\t" << sum_prog[i] << " \t" << err_prog[i] << endl;
    data1.close();

	rnd.SaveSeed();
	return 0;
}