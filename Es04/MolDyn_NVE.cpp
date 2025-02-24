/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
#include <stdlib.h>     // srand, rand: to generate random number
#include <iostream>     // cin, cout: Standard Input/Output Streams Library
#include <fstream>      // Stream class to both read and write from/to files.
#include <cmath>        // rint, pow
#include <string>
#include "MolDyn_NVE.h"

using namespace std;

int main(){ 
  Input();             //Inizialization
  int nconf = 1;
  double* U = new double[nstep];
  double* K = new double[nstep];
  double* E = new double[nstep];
  double* T = new double[nstep];
  for(int istep=1; istep <= nstep; ++istep){
     Move();           //Move particles with Verlet algorithm
     if(istep%iprint == 0) cout << "Number of time-steps: " << istep << endl;
     Measure();     //Properties measurement
     U[istep-1] = stima_pot;
     E[istep-1] = stima_etot;
     K[istep-1] = stima_kin;
     T[istep-1] = stima_temp;
     if (istep == (nstep - 1))
         OldFinal();
  }
  ConfFinal();         //Write final configuration to restart

  //blocking method
  int N = nblocks;
  int L = nstep / N;

  double* ave = new double[N] {0};
  double* av2 = new double[N] {0};
  double* sum_prog = new double[N] {0};
  double* su2_prog = new double[N] {0};
  double* err_prog = new double[N] {0};
  double sums;
  for (int i = 0;i < N; i++) {
      sums = 0;
      for (int j = 0;j < L;j++) {
          int k = j + i * L;
          sums = sums + U[k];     
      }
      ave[i] = sums / (double)L;
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

  ofstream data1("ave_epot.out");
  for (int i = 0;i < N;i++)
      data1 << sum_prog[i] << " \t" << err_prog[i] << endl;
  data1.close();

  //resetto
  for (int i = 0;i < N;i++) {
      ave[i] = 0;
      av2[i] = 0;
      sum_prog[i] = 0;
      su2_prog[i] = 0;
      err_prog[i] = 0;
  }
  for (int i = 0;i < N; i++) {
      sums = 0;
      for (int j = 0;j < L;j++) {
          int k = j + i * L;
          sums = sums + K[k];     
      }
      ave[i] = sums / (double)L;
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
  ofstream data2("ave_ekin.out");
  for (int i = 0;i < N;i++)
      data2 << sum_prog[i] << " \t" << err_prog[i] << endl;
  data2.close();


  //resetto
  for (int i = 0;i < N;i++) {
      ave[i] = 0;
      av2[i] = 0;
      sum_prog[i] = 0;
      su2_prog[i] = 0;
      err_prog[i] = 0;
  }
  for (int i = 0;i < N; i++) {
      sums = 0;
      for (int j = 0;j < L;j++) {
          int k = j + i * L;
          sums = sums + E[k];     
      }
      ave[i] = sums / (double)L;
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
  ofstream data3("ave_etot.out");
  for (int i = 0;i < N;i++)
      data3 << sum_prog[i] << " \t" << err_prog[i] << endl;
  data3.close();

  //resetto
  for (int i = 0;i < N;i++) {
      ave[i] = 0;
      av2[i] = 0;
      sum_prog[i] = 0;
      su2_prog[i] = 0;
      err_prog[i] = 0;
  }
  for (int i = 0;i < N; i++) {
      sums = 0;
      for (int j = 0;j < L;j++) {
          int k = j + i * L;
          sums = sums + T[k];     
      }
      ave[i] = sums / (double)L;
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
  ofstream data4("ave_temp.out");
  for (int i = 0;i < N;i++)
      data4 << sum_prog[i] << " \t" << err_prog[i] << endl;
  data4.close();

  return 0;
}

double error(double* AV, double* AV2, int n) {
    if (n == 0)
        return 0;
    else
        return sqrt((AV2[n] - pow(AV[n], 2)) / n);
}

void Input(void){ //Prepare all stuff for the simulation
  ifstream ReadInput,ReadConf;
  double ep, ek, pr, et, vir;
  int method;

  cout << "Classic Lennard-Jones fluid        " << endl;
  cout << "Molecular dynamics simulation in NVE ensemble  " << endl << endl;
  cout << "Interatomic potential v(r) = 4 * [(1/r)^12 - (1/r)^6]" << endl << endl;
  cout << "The program uses Lennard-Jones units " << endl;

  seed = 1;    //Set seed for random numbers
  srand(seed); //Initialize random number generator
  
  ReadInput.open("input.dat"); //Read input

  ReadInput >> temp;

  ReadInput >> npart;
  cout << "Number of particles = " << npart << endl;

  ReadInput >> rho;
  cout << "Density of particles = " << rho << endl;
  vol = (double)npart/rho;
  cout << "Volume of the simulation box = " << vol << endl;
  box = pow(vol,1.0/3.0);
  cout << "Edge of the simulation box = " << box << endl;

  ReadInput >> rcut;
  ReadInput >> delta;
  ReadInput >> nstep;
  ReadInput >> iprint;

  cout << "The program integrates Newton equations with the Verlet method " << endl;
  cout << "Time step = " << delta << endl;
  cout << "Number of steps = " << nstep << endl << endl;
  ReadInput >> method;
  ReadInput >> nblocks;
  ReadInput.close();

//Prepare array for measurements
  iv = 0; //Potential energy
  ik = 1; //Kinetic energy
  ie = 2; //Total energy
  it = 3; //Temperature
  n_props = 4; //Number of observables

  if (method == 0) {
      //Read initial configuration
      cout << "Read initial configuration from file config.0 " << endl << endl;
      ReadConf.open("config.0");
      for (int i = 0; i < npart; ++i) {
          ReadConf >> x[i] >> y[i] >> z[i];
          x[i] = x[i] * box;
          y[i] = y[i] * box;
          z[i] = z[i] * box;
      }
      ReadConf.close();

      //Prepare initial velocities
      cout << "Prepare random velocities with center of mass velocity equal to zero " << endl << endl;
      double sumv[3] = { 0.0, 0.0, 0.0 };
      for (int i = 0; i < npart; ++i) {
          vx[i] = rand() / double(RAND_MAX) - 0.5;
          vy[i] = rand() / double(RAND_MAX) - 0.5;
          vz[i] = rand() / double(RAND_MAX) - 0.5;

          sumv[0] += vx[i];
          sumv[1] += vy[i];
          sumv[2] += vz[i];
      }
      for (int idim = 0; idim < 3; ++idim) sumv[idim] /= (double)npart;
      double sumv2 = 0.0, fs;
      for (int i = 0; i < npart; ++i) {
          vx[i] = vx[i] - sumv[0];
          vy[i] = vy[i] - sumv[1];
          vz[i] = vz[i] - sumv[2];

          sumv2 += vx[i] * vx[i] + vy[i] * vy[i] + vz[i] * vz[i];
      }
      sumv2 /= (double)npart;

      fs = sqrt(3 * temp / sumv2);   // fs = velocity scale factor 
      for (int i = 0; i < npart; ++i) {
          vx[i] *= fs;
          vy[i] *= fs;
          vz[i] *= fs;

          xold[i] = Pbc(x[i] - vx[i] * delta);
          yold[i] = Pbc(y[i] - vy[i] * delta);
          zold[i] = Pbc(z[i] - vz[i] * delta);
      }
  }
  else {
      //Read the two initial configurations
      cout << "Read pre-initial configuration from file old.0 " << endl << endl;
      ReadConf.open("old.0");
      for (int i = 0; i < npart; ++i) {
          ReadConf >> xold[i] >> yold[i] >> zold[i];
          xold[i] = xold[i] * box;
          yold[i] = yold[i] * box;
          zold[i] = zold[i] * box;
      }
      ReadConf.close();
      ifstream ReadConf2;
      ReadConf2.open("config.0");
      cout << "Read initial configuration from file config.0 " << endl << endl;
      for (int i = 0; i < npart; ++i) {
          ReadConf2 >> x[i] >> y[i] >> z[i];
          x[i] = x[i] * box;
          y[i] = y[i] * box;
          z[i] = z[i] * box;
      }
      ReadConf2.close();
      Move();
      double sumv2 = 0.0,fs;
      for (int i = 0;i < npart;i++) {
          vx[i] = Pbc(x[i] - xold[i])/delta;
          vy[i] = Pbc(y[i] - yold[i])/delta;
          vz[i] = Pbc(z[i] - zold[i])/delta;
          sumv2 += vx[i] * vx[i] + vy[i] * vy[i] + vz[i] * vz[i];
      }
      sumv2 /= (double)npart;
      double T;
      T = sumv2 / 3.; //Questo � kB*T in udM
      fs = sqrt(temp/T);    //da verificare
      for (int i = 0; i < npart; ++i) {
          vx[i] *= fs;
          vy[i] *= fs;
          vz[i] *= fs;
          xold[i] = Pbc(x[i] - vx[i] * delta);
          yold[i] = Pbc(y[i] - vy[i] * delta);
          zold[i] = Pbc(z[i] - vz[i] * delta);
      }
  }
   return;
}

void Move(void){ //Move particles with Verlet algorithm
  double xnew, ynew, znew, fx[m_part], fy[m_part], fz[m_part];

  for(int i=0; i<npart; ++i){ //Force acting on particle i
    fx[i] = Force(i,0);
    fy[i] = Force(i,1);
    fz[i] = Force(i,2);
  }

  for(int i=0; i<npart; ++i){ //Verlet integration scheme

    xnew = Pbc( 2.0 * x[i] - xold[i] + fx[i] * pow(delta,2) );
    ynew = Pbc( 2.0 * y[i] - yold[i] + fy[i] * pow(delta,2) );
    znew = Pbc( 2.0 * z[i] - zold[i] + fz[i] * pow(delta,2) );

    vx[i] = Pbc(xnew - xold[i])/(2.0 * delta);
    vy[i] = Pbc(ynew - yold[i])/(2.0 * delta);
    vz[i] = Pbc(znew - zold[i])/(2.0 * delta);

    xold[i] = x[i];
    yold[i] = y[i];
    zold[i] = z[i];

    x[i] = xnew;
    y[i] = ynew;
    z[i] = znew;
  }
  return;
}

double Force(int ip, int idir){ //Compute forces as -Grad_ip V(r)
  double f=0.0;
  double dvec[3], dr;

  for (int i=0; i<npart; ++i){
    if(i != ip){
      dvec[0] = Pbc( x[ip] - x[i] );  // distance ip-i in pbc
      dvec[1] = Pbc( y[ip] - y[i] );
      dvec[2] = Pbc( z[ip] - z[i] );

      dr = dvec[0]*dvec[0] + dvec[1]*dvec[1] + dvec[2]*dvec[2];
      dr = sqrt(dr);

      if(dr < rcut){
        f += dvec[idir] * (48.0/pow(dr,14) - 24.0/pow(dr,8)); // -Grad_ip V(r)
      }
    }
  }
  
  return f;
}

void Measure(){ //Properties measurement
  int bin;
  double v, t, vij;
  double dx, dy, dz, dr;
  ofstream Epot, Ekin, Etot, Temp;

  Epot.open("output_epot.dat",ios::app);
  Ekin.open("output_ekin.dat",ios::app);
  Temp.open("output_temp.dat",ios::app);
  Etot.open("output_etot.dat",ios::app);

  v = 0.0; //reset observables
  t = 0.0;

//cycle over pairs of particles
  for (int i=0; i<npart-1; ++i){
    for (int j=i+1; j<npart; ++j){

     dx = Pbc( x[i] - x[j] );
     dy = Pbc( y[i] - y[j] );
     dz = Pbc( z[i] - z[j] );

     dr = dx*dx + dy*dy + dz*dz;
     dr = sqrt(dr);

     if(dr < rcut){
       vij = 4.0/pow(dr,12) - 4.0/pow(dr,6);

//Potential energy
       v += vij;
     }
    }          
  }

//Kinetic energy
  for (int i=0; i<npart; ++i) t += 0.5 * (vx[i]*vx[i] + vy[i]*vy[i] + vz[i]*vz[i]);
   
    stima_pot = v/(double)npart; //Potential energy per particle
    stima_kin = t/(double)npart; //Kinetic energy per particle
    stima_temp = (2.0 / 3.0) * t/(double)npart; //Temperature
    stima_etot = (t+v)/(double)npart; //Total energy per particle

    Epot << stima_pot  << endl;
    Ekin << stima_kin  << endl;
    Temp << stima_temp << endl;
    Etot << stima_etot << endl;

    Epot.close();
    Ekin.close();
    Temp.close();
    Etot.close();

    return;
}

void ConfFinal(void){ //Write final configuration
  ofstream WriteConf;

  cout << "Print final configuration to file config.final " << endl << endl;
  WriteConf.open("config.final");

  for (int i=0; i<npart; ++i){
    WriteConf << x[i]/box << "   " <<  y[i]/box << "   " << z[i]/box << endl;
  }
  WriteConf.close();
  return;
}

void OldFinal(void) { //Write final configuration
    ofstream WriteConf;

    cout << "Print pre-final configuration to file old.final " << endl << endl;
    WriteConf.open("old.final");

    for (int i = 0; i < npart; ++i) {
        WriteConf << x[i] / box << "   " << y[i] / box << "   " << z[i] / box << endl;
    }
    WriteConf.close();
    return;
}

void ConfXYZ(int nconf){ //Write configuration in .xyz format
  ofstream WriteXYZ;

  WriteXYZ.open("frames/config_" + to_string(nconf) + ".xyz");
  WriteXYZ << npart << endl;
  WriteXYZ << "This is only a comment!" << endl;
  for (int i=0; i<npart; ++i){
    WriteXYZ << "LJ  " << Pbc(x[i]) << "   " <<  Pbc(y[i]) << "   " << Pbc(z[i]) << endl;
  }
  WriteXYZ.close();
}

double Pbc(double r){  //Algorithm for periodic boundary conditions with side L=box
    return r - box * rint(r/box);
}
/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
