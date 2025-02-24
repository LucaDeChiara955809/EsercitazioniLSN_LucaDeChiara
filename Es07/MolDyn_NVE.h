/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
const int m_props = 1000;
int n_props, iv, iw, igofr;
double vtail, ptail, bin_size, nbins, sd;
double walker[m_props];

// averages
double blk_av[m_props], blk_norm, accepted, attempted;
double glob_av[m_props], glob_av2[m_props];
double stima_pot, stima_pres, err_pot, err_press, err_gdir;

//configuration
const int m_part = 108;
double x[m_part], y[m_part], z[m_part], xold[m_part], yold[m_part], zold[m_part];
double vx[m_part], vy[m_part], vz[m_part];

// boolean input
bool method;

// thermodynamical state
int npart;
double energy, temp, vol, rho, box, rcut;

// simulation
int nstep, iprint, seed;
double delta;

//blocking method
int nblk = 50;

//functions
void Input(void);
void Move(void);
void ConfFinal(void);
void OldFinal(void);
void ConfXYZ(int);
void Measure(void);
void Reset(int);
void Accumulate(void);
void Averages(int);
double Force(int, int);
double Pbc(double);
double Error(double, double, int);
/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
