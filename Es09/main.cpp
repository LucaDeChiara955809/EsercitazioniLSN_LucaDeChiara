#define _USE_MATH_DEFINES

#include<iostream>
#include<fstream>
#include"random.h"
#include"Individuo.h"
#include"function.h"
#include<cmath>

using namespace std;

int main() {
    RandomInitializer();
    Individuo* i = new Individuo[population];
    RandomPopulationGenerator(i);
    PositionGenerator(metodo);
    Individuo* figli = new Individuo[population];
    for (int y = 0;y < Ntimes;y++) {
        Quicksort(i, 0, population - 1);
        r_mean[y] = mean(i);
        r_min[y] = L(i[0]);
        r_devSt[y] = DevSt(i);
        if (Lmin > L(i[0])) {   //registro il percorso migliore in assoluto
            Lmin = L(i[0]);
            best.SetSequence(i[0].GetSequence(),32);    
        }
        if (y % 1000 == 0)
            cout << L(i[0])<<" "<<L(i[99]) << endl; //stampa di controllo
        for (int w = 0;w < population;w=w+2) {
            madre = ParentSelector(i);
            padre = ParentSelector(i);
            double rr = rnd.Rannyu();
            if (y > 20)
                prob = 0.03;    //faccio in modo che all'inizio abbia una probabilità di poco maggiore di avere mutazioni
            else
                prob = 0.10;
            if (rr <= prob) {
                SingleMutation(padre);
                SingleMutation(madre);
            }
            else if (rr > prob && rr <= 2*prob) {
                ShiftMutation(padre);
                ShiftMutation(madre);
            }
            else if (rr > 2*prob && rr <= 3*prob) {
                PermutationMutation(padre);
                PermutationMutation(madre);
            }
            else if (rr > 3*prob && rr <= 4*prob) {
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

    ofstream dataL;
    dataL.open("dataL.dat");
    for (int k = 0;k < Ntimes+1;k++)
        dataL << r_min[k]<<"\t"<<r_mean[k]<<"\t"<<r_devSt[k] << endl;
    dataL.close();

    ofstream data;
    data.open("data.dat");
    for (int k = 0;k < Ncities;k++)
        data << positions[best.GetCity(k)-1][0] << "\t" << positions[best.GetCity(k)-1][1] << endl;
    data.close();
    rnd.SaveSeed();
    return 0;
}

void RandomInitializer() {
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
    return;
}

void PositionGenerator(int metodo) {
    if (metodo == 0) {
        double theta;
        //disponiamo le città SU una crf centrata in (0;0) di raggio 1
        for (int i = 0;i < Ncities;i++) {
            theta = rnd.Theta();
            positions[i][0] = cos(theta);
            positions[i][1] = sin(theta);
        }
    }
    else {
        //disponiamo le città IN un quadrato di lato 1 con vertici in (0;0),(1;0),(1;1) e (0;1)
        for (int i = 0;i < Ncities;i++) {
            positions[i][0] = rnd.Rannyu();
            positions[i][1] = rnd.Rannyu();
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
        L += sqrt(pow(positions[y.GetCity(j)-1][0] - positions[y.GetCity(mm)-1][0], 2) + pow(positions[y.GetCity(j)-1][1] - positions[y.GetCity(mm)-1][1], 2));
    }
    //calcolo di L(1) per un individuo
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
                if (random_seq[j_s] == random_seq[i_s]) {   //se è già nella sequenza
                    j_s = i_s;  //esco dal while
                    i_s = i_s - 1;  //genero un altro indice casuale
                    if (i_s < 0)
                        cout << "Error: i_s < 0" << endl;   //parametro di controllo
                }
                else {
                    j_s++;
                }
            }
            i_s++;
        }
        i[n].SetSequence(random_seq, Ncities);
    }
    return;
}

int CheckFunction(Individuo y) {    //verifica solo se non ho ripetizioni nella città in quanto le etichette vanno da 1 a 32
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
    double p=5; //potenza scelta a seguito di ottimizzazione
    int j = (int)floor(((double)population) * pow(rnd.Rannyu(), p));
    return i[j];
}

void SingleMutation(Individuo y) {  //seleziono due posizioni e scambio le etichette
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

void ShiftMutation(Individuo y) { //shift di +n posizioni per m città contigue
    int o1 = (int)floor(rnd.Rannyu(1., 32.));   //posizione scelta casualmente
    int n = (int)floor(rnd.Rannyu(1., 16.));    //Sono più piccoli perchè serve lo "spazio di manovra" (e.g. non posso spostare avanti di 9 30 città perchè mi manca lo spazio)
    int m = (int)floor(rnd.Rannyu(1., 16.));
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

void InversionMutation(Individuo y) {   //inverto l'ordine di n etichette 
    int o1 = (int)floor(rnd.Rannyu(1., 32.));
    int n = (int)floor(rnd.Rannyu(1., 32.));
    int* appo = new int[n];
    for (int k = 0;k < n;k++)
        appo[n-k-1] = y.GetCity(PBC(o1 + k));
    for (int k = o1;k < o1 + n;k++)
        y.SetCity(PBC(k), appo[k-o1]);
    return;
}

void PermutationMutation(Individuo y) { //scambio di m città contigue con altrettante contigue
    int o1 = (int)floor(rnd.Rannyu(1., 32.));
    int n = (int)floor(rnd.Rannyu(1., 16.));    //ridotto perchè mi serve "spazio di manovra"
    int o2=PBC(o1 + n + (int)floor(rnd.Rannyu(0., 31. - 2 * (double)n)));   //prendo o2 di modo che non appartenga alle città da spostare
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
    int cut = (int)floor(rnd.Rannyu(2., 31.));  //scelgo un sito casuale in cui tagliare (COMPRESO) da 2 a 30 (non ha senso tagliare in 1 o 31 in quanto avrei lo stesso ordine)
    int* missing1 = new int[Ncities - cut];
    int* missing2 = new int[Ncities - cut];
    for (int k = 0;k < Ncities - cut;k++) {
        missing1[k] = p1.GetCity(PBC(cut + k)); //registro le etichette perse dopo il taglio
        missing2[k] = p2.GetCity(PBC(cut + k));
    }
    int* neworder1 = new int[Ncities - cut];
    int* neworder2 = new int[Ncities - cut];
    int k = 0;
    while (k<Ncities-cut){
        for (int j = 1;j < Ncities;j++) {
            for (int l = 0;l < Ncities - cut;l++) {
                if (p2.GetCity(j) == missing1[l]) {
                    neworder1[k] = p2.GetCity(j);   //registro due vettori new_order che avranno solo le etichette perse in un nuovo ordine (dipendente dall'altro individuo)
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
        p1.SetCity(cut + k, neworder1[k]);  //registro le etichette perse nel nuovo ordine senza cambiare l'ordine precedente
        p2.SetCity(cut + k, neworder2[k]);
    }
    return;
}

void Quicksort(Individuo* y, int first, int last) { //algoritmo di ordinamento
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
        return i - 31 * rint(i / 32);   //funzione periodica usata per le posizioni dell'array di modo da non cambiare mai la posizione 0
}

int PBC0(int i) {
    return i - 32 * rint(i / 32);   //funzione poeriodica
}

double mean(Individuo* y) { //effettuo la media sulla migliore metà di popolazione
    double mean = 0;
    for (int i = 0;i < population / 2;i++)
        mean += L(y[i]);
    return mean / (population / 2);
}

double DevSt(Individuo* y) {    //errore della media ottenuto come deviazione standard su metà della popolazione
    double m = mean(y);
    double mean2 = 0;
    for (int i = 0;i < population / 2;i++)
        mean2 += pow(L(y[i]),2);
    mean2 = mean2 / (population / 2);
    return sqrt(abs(mean2 - m * m));
}