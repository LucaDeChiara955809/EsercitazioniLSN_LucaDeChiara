#define _USE_MATH_DEFINES

#include<iostream>
#include<fstream>
#include"random.h"
#include"Individuo.h"
#include"function.h"
#include"mpi.h"
#include<cmath>

using namespace std;

int main(int argc, char* argv[]) {
    int size, rank; 
    MPI_Init(&argc, &argv); 
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    RandomInitializer(rank);
    Individuo* i = new Individuo[population];
    RandomPopulationGenerator(i);
    if (rank==0)
        PositionGenerator(metodo);
    MPI_Bcast(&X, Ncities, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD);
    MPI_Bcast(&Y, Ncities, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD);
    Individuo* figli = new Individuo[population];
    for (int y = 0;y < Ntimes;y++) {
        Quicksort(i, 0, population - 1);
        r_mean[y] = mean(i);
        r_min[y] = L(i[0]);
        r_devSt[y] = DevSt(i);
        if (y != 0 && y % Nmigr == 0) {
            MPI_Status stat0, stat1, stat2, stat3;
            int recv[Ncities]{ 0 };
            int send[Ncities];
            for (int h = 0;h < Ncities;h++)
                send[h] = i[0].GetCity(h);
            int l1;
            if (rank==0)
                l1 = (int)floor(rnd.Rannyu(0, 3));
            MPI_Bcast(&l1, 1, MPI_INTEGER, 0, MPI_COMM_WORLD);
            if (l1 == 0) {
                if (rank == 0) {
                    MPI_Send(send, Ncities, MPI_INTEGER, 1, 1, MPI_COMM_WORLD);
                    MPI_Recv(recv, Ncities, MPI_INTEGER, 1, 2, MPI_COMM_WORLD, &stat1);
                }
                else if (rank == 1) {
                    MPI_Send(send, Ncities, MPI_INTEGER, 0, 2, MPI_COMM_WORLD);
                    MPI_Recv(recv, Ncities, MPI_INTEGER, 0, 1, MPI_COMM_WORLD, &stat0);
                }
                if (rank == 2) {
                    MPI_Send(send, Ncities, MPI_INTEGER, 3, 3, MPI_COMM_WORLD);
                    MPI_Recv(recv, Ncities, MPI_INTEGER, 3, 4, MPI_COMM_WORLD, &stat2);
                }
                else if (rank == 3) {
                    MPI_Send(send, Ncities, MPI_INTEGER, 2, 4, MPI_COMM_WORLD);
                    MPI_Recv(recv, Ncities, MPI_INTEGER, 2, 3, MPI_COMM_WORLD, &stat3);
                }
            }
            else if (l1 == 1) {
                if (rank == 0) {
                    MPI_Send(send, Ncities, MPI_INTEGER, 2, 1, MPI_COMM_WORLD);
                    MPI_Recv(recv, Ncities, MPI_INTEGER, 2, 2, MPI_COMM_WORLD, &stat1);
                }
                else if (rank == 2) {
                    MPI_Send(send, Ncities, MPI_INTEGER, 0, 2, MPI_COMM_WORLD);
                    MPI_Recv(recv, Ncities, MPI_INTEGER, 0, 1, MPI_COMM_WORLD, &stat0);
                }
                if (rank == 1) {
                    MPI_Send(send, Ncities, MPI_INTEGER, 3, 3, MPI_COMM_WORLD);
                    MPI_Recv(recv, Ncities, MPI_INTEGER, 3, 4, MPI_COMM_WORLD, &stat2);
                }
                else if (rank == 3) {
                    MPI_Send(send, Ncities, MPI_INTEGER, 1, 4, MPI_COMM_WORLD);
                    MPI_Recv(recv, Ncities, MPI_INTEGER, 1, 3, MPI_COMM_WORLD, &stat3);
                }
            }
            else {
                if (rank == 0) {
                    MPI_Send(send, Ncities, MPI_INTEGER, 3, 1, MPI_COMM_WORLD);
                    MPI_Recv(recv, Ncities, MPI_INTEGER, 3, 2, MPI_COMM_WORLD, &stat1);
                }
                else if (rank == 3) {
                    MPI_Send(send, Ncities, MPI_INTEGER, 0, 2, MPI_COMM_WORLD);
                    MPI_Recv(recv, Ncities, MPI_INTEGER, 0, 1, MPI_COMM_WORLD, &stat0);
                }
                if (rank == 1) {
                    MPI_Send(send, Ncities, MPI_INTEGER, 2, 3, MPI_COMM_WORLD);
                    MPI_Recv(recv, Ncities, MPI_INTEGER, 2, 4, MPI_COMM_WORLD, &stat2);
                }
                else if (rank == 2) {
                    MPI_Send(send, Ncities, MPI_INTEGER, 1, 4, MPI_COMM_WORLD);
                    MPI_Recv(recv, Ncities, MPI_INTEGER, 1, 3, MPI_COMM_WORLD, &stat3);
                }
            }
            i[0].SetSequence(recv, Ncities);
        }
        if (Lmin > L(i[0])) {
            Lmin = L(i[0]);
            best.SetSequence(i[0].GetSequence(), 32);
        }

        if (y % 500 == 0) {
            cout << "Rank:" << rank << endl;
            cout << L(i[0]) << " " << L(i[9]) << endl;
            cout << "Vivo fino a qui " << y << endl;
        }
        for (int w = 0;w < population;w = w + 2) {
            madre = ParentSelector(i);
            padre = ParentSelector(i);
            double rr = rnd.Rannyu();
            if (y > 20)
                prob = 0.03;
            else
                prob = 0.10;
            if (rr <= prob) {
                SingleMutation(padre);
                SingleMutation(madre);
            }
            else if (rr > prob && rr <= 2 * prob) {
                ShiftMutation(padre);
                ShiftMutation(madre);
            }
            else if (rr > 2 * prob && rr <= 3 * prob) {
                PermutationMutation(padre);
                PermutationMutation(madre);
            }
            else if (rr > 3 * prob && rr <= 4 * prob) {
                InversionMutation(padre);
                InversionMutation(madre);
            }
            else
                CrossingOver(padre, madre);
            CheckFunction(padre);
            if (controllo == 1) {
                cout << "errore padre" << endl;
                cout << y << "  " << w << endl;
            }
            CheckFunction(madre);
            if (controllo == 1) {
                cout << "errore padre" << endl;
                cout << y << "  " << w << endl;
            }
            figli[w].SetSequence(padre.GetSequence(), Ncities);
            figli[w + 1].SetSequence(madre.GetSequence(), Ncities);
        }
        for (int k = 0;k < population;k++)
            i[k].SetSequence(figli[k].GetSequence(), Ncities);
    }
    Quicksort(i, 0, population - 1);
    r_mean[Ntimes] = mean(i);
    r_min[Ntimes] = L(i[0]);
    r_devSt[Ntimes] = DevSt(i);
    //cout << "ho tutti i dati ora li carico" << endl;
    ofstream dataL;
    if (rank == 0)
        dataL.open("../dataL0.dat");
    else if (rank == 1)
        dataL.open("../dataL1.dat");
    else if (rank == 2)
        dataL.open("../dataL2.dat");
    else if (rank == 3)
        dataL.open("../dataL3.dat");
    //cout << "ho aperto l'ofstream desiderato?" << endl;
    for (int k = 0;k < Ntimes+1;k++)
        dataL << r_min[k]<<"\t"<<r_mean[k]<<"\t"<<r_devSt[k] << endl;

    //cout << "ho fatto, sono solo lento" << endl;
    dataL.close();


    ofstream data;
    if (rank == 0)
        data.open("../data0.dat");
    else if (rank == 1)
        data.open("../data1.dat");
    else if (rank == 2)
        data.open("../data2.dat");
    else if (rank == 3)
        data.open("../data3.dat");
    for (int k = 0;k < Ncities;k++)
        data << X[best.GetCity(k)-1] << "\t" << Y[best.GetCity(k)-1] << endl;
    data.close();
    rnd.SaveSeed();
    MPI_Finalize();
    return 0;
}

void RandomInitializer(int rank) {
    ifstream Primes;
    Primes.open("../Primes");      
    if (Primes.is_open()) {
        for(int k=0;k<=rank;k++)
            Primes >> p1 >> p2;
    }
    else cerr << "PROBLEM: Unable to open Primes" << endl;    
    Primes.close();
    ifstream input("../seed.in");
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
    return;
}

void PositionGenerator(int metodo) {
    if (metodo == 0) {
        double theta;
        //disponiamo le città SU una crf
        for (int i = 0;i < Ncities;i++) {
            theta = rnd.Theta();
            X[i]= cos(theta);
            Y[i]= sin(theta);
        }
    }
    else {
        //disponiamo le città IN un quadrato
        for (int i = 0;i < Ncities;i++) {
            X[i] = rnd.Rannyu();
            Y[i] = rnd.Rannyu();
        }
    }
    return;
}

double L(Individuo y) {
    double L = 0;
    int mm;
    for (int j = 0; j < Ncities; j++) {
        if (j == Ncities - 1)
            mm = 0;
        else
            mm = j + 1;
        L += sqrt(pow(X[y.GetCity(j)-1] - X[y.GetCity(mm)-1], 2) + pow(Y[y.GetCity(j)-1] - Y[y.GetCity(mm)-1], 2));
    }
    return L;
}

void RandomPopulationGenerator(Individuo* i) {
    int* random_seq = new int[Ncities];
    int i_s = 0;
    int j_s = 0;
    for (int n = 0; n < population; n++) { //cicliamo sui differenti individui della popolazione
        i_s = 1;
        for (int k = 0;k < Ncities;k++)
            random_seq[k] = 0;      //resetto random seq
        random_seq[0] = 1;//la città di partenza è sempre la stessa
        while (i_s < Ncities) {
            random_seq[i_s] = (int)floor(rnd.Rannyu(2., 33.)); //genero un etichetta random
            int j_s = 1;
            while (j_s < i_s) { //controllo se l'indice è già nella sequenza
                if (random_seq[j_s] == random_seq[i_s]) { //se è già nella sequenza
                    j_s = i_s; //esco dal while
                    i_s = i_s - 1; //genero un altro indice casuale
                    if (i_s < 0)
                        cout << "Error: i_s < 0" << endl;
                }
                else {
                    j_s++;
                }
            }
            i_s++;
        }
        i[n].SetSequence(random_seq, Ncities);
        //cout << endl << endl;
    }
    return;
}

int CheckFunction(Individuo y) {
    controllo = 0;
    for (int i = 1;i < Ncities;i++) {
        for (int j = 1;j < i;j++) {
            if (y.GetCity(j) == y.GetCity(i))
                controllo = 1;
        }
    }
    return controllo;
}

Individuo ParentSelector(Individuo* i) {
    double p=5;
    int j = (int)floor(((double)population) * pow(rnd.Rannyu(), p));
    return i[j];
}

void SingleMutation(Individuo y) {
    int o1 = (int)floor(rnd.Rannyu(1., 32.));
    int o2 = (int)floor(rnd.Rannyu(1., 32.));
    while (o2==o1)
        o2 = (int)floor(rnd.Rannyu(1., 32.));
    int appo;
    appo = y.GetCity(o1);
    y.SetCity(o1,y.GetCity(o2));
    y.SetCity(o2,appo);
    return;
}

void ShiftMutation(Individuo y) { //shift of +n positions for m contiguous cities
    int o1 = (int)floor(rnd.Rannyu(1., 32.));
    int n = (int)floor(rnd.Rannyu(1., 16.));    //Sono più piccoli perchè serve lo "spazio di manovra" (non posso spostare avanti di 9 30 città perchè mi manca lo spazio)
    int m = (int)floor(rnd.Rannyu(1., 16.));
    //cout << o1 << "  " << n << "  " << m << endl;
    int* appoM = new int[m];
    int* appoN = new int[n];
    for (int k = 0; k < m; k++)
        appoM[k] = y.GetCity(PBC(o1 + k));
    for (int k = 0; k < n; k++)
        appoN[k] = y.GetCity(PBC(o1 + m + k));
    for (int k = 0; k < n; k++)
        y.SetCity(PBC(o1 + k), appoN[k]);
    for (int k = 0; k < m; k++)
        y.SetCity(PBC(o1 + k + n), appoM[k]);   
    return;
}

void InversionMutation(Individuo y) {
    int o1 = (int)floor(rnd.Rannyu(1., 32.));
    int n = (int)floor(rnd.Rannyu(1., 32.));
    int* appo = new int[n];
    for (int k = 0;k < n;k++)
        appo[n-k-1] = y.GetCity(PBC(o1 + k));
    for (int k = o1;k < o1 + n;k++)
        y.SetCity(PBC(k), appo[k-o1]);
    return;
}

void PermutationMutation(Individuo y) {
    int o1 = (int)floor(rnd.Rannyu(1., 32.));
    int n = (int)floor(rnd.Rannyu(1., 16.));
    int o2=PBC(o1 + n + (int)floor(rnd.Rannyu(0., 31. - 2 * (double)n)));   //spazio di manovra
    //cout << o1 << "   " << o2 << "   " << n << endl;
    int* appo1 = new int[n];
    int* appo2 = new int[n];
    for (int k = 0;k < n;k++) {
        appo1[k] = y.GetCity(PBC(o1 + k));
        appo2[k] = y.GetCity(PBC(o2 + k));
    }
    for (int k = 0;k < n;k++) {
        y.SetCity(PBC(o1 + k), appo2[k]);
        y.SetCity(PBC(o2 + k), appo1[k]);
    }
    return;
}

void CrossingOver(Individuo p1, Individuo p2) {
    int cut = (int)floor(rnd.Rannyu(2., 31.));  //scelgo un sito casuale in cui tagliare (COMPRESO) da 2 a 31
    int* missing1 = new int[Ncities - cut];
    int* missing2 = new int[Ncities - cut];
    for (int k = 0;k < Ncities - cut;k++) {
        missing1[k] = p1.GetCity(PBC(cut + k));
        missing2[k] = p2.GetCity(PBC(cut + k));
    }
    int* neworder1 = new int[Ncities - cut];
    int* neworder2 = new int[Ncities - cut];
    int k = 0;
    while (k<Ncities-cut){
        for (int j = 1;j < Ncities;j++) {
            for (int l = 0;l < Ncities - cut;l++) {
                if (p2.GetCity(j) == missing1[l]) {
                    neworder1[k] = p2.GetCity(j);
                    k++;
                }
            }
        }
    }
    int k2 = 0;
    while (k2 < Ncities - cut) {
        for (int j = 1;j < Ncities;j++) {
            for (int l = 0;l < Ncities - cut;l++) {
                if (p1.GetCity(j) == missing2[l]) {
                    neworder2[k2] = p1.GetCity(j);
                    k2++;
                }
            }
        }
    }
    for (int k = 0;k < Ncities - cut;k++) {
        p1.SetCity(cut + k, neworder1[k]);
        p2.SetCity(cut + k, neworder2[k]);
    }
    return;
}

void Quicksort(Individuo* y, int first, int last) {
    int i, j, pivot;
    Individuo temp;
    if (first < last) {
        pivot = first;
        i = first;
        j = last;
        while (i < j) {
            while (L(y[i]) <= L(y[pivot]) && i < last)
                i++;
            while (L(y[j]) > L(y[pivot]))
                j--;
            if (i < j) {
                temp = y[i];
                y[i] = y[j];
                y[j] = temp;
            }
        }
        temp = y[pivot];
        y[pivot] = y[j];
        y[j] = temp;
        Quicksort(y, first, j - 1);
        Quicksort(y, j + 1, last);
    }
}

int PBC(int i) {
        return i - 31 * rint(i / 32);
}

int PBC0(int i) {
    return i - 32 * rint(i / 32);
}

double mean(Individuo* y) {
    double mean = 0;
    for (int i = 0;i < population / 2;i++)
        mean += L(y[i]);
    return mean / (population / 2);
}

double DevSt(Individuo* y) {
    double m = mean(y);
    double mean2 = 0;
    for (int i = 0;i < population / 2;i++)
        mean2 += pow(L(y[i]),2);
    mean2 = mean2 / (population / 2);
    return sqrt(mean2 - m * m);
}