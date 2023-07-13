#ifndef __fluid__
#define __fluid__

//Random numbers
#include "random.h"
int seed[4];
Random rnd;

//parameters, observables
const int m_props=1000;
int n_props, iv, ik, it, ie, iw, ip;
double vtail,wtail, bin_size, nbins, sd;
double walker[m_props];

const int ng = 100;
double gdr[ng];
double gdr_ave[ng];

double min_dist=0; // risoluzione dell'istogramma della gdr
int bin_index=0;

// averages
double blk_av[m_props], blk_norm, accepted, attempted;
double glob_av[m_props], glob_av2[m_props];
double stima_pot, stima_pres, stima_kin, stima_etot, stima_temp;
double err_pot, err_pres, err_kin, err_etot, err_temp, err_gdir;

//configuration
const int m_part=108;
double x[m_part],    y[m_part],    z[m_part];
double xold[m_part], yold[m_part], zold[m_part];
double vx[m_part],  vy[m_part],   vz[m_part];

// thermodynamical state
int npart;
double beta,temp,energy,vol,rho,box,rcut;

// simulation
int iNVET, nstep, nblk, restart;
double delta;

//pigreco
const double pi=3.1415927;
//costante Boltzmann
const double Kb=1.380649E-23;

//functions
void Input(std::string);
void Reset(int);
void Accumulate(void);
void Averages(int,std::string);
void Move(void);
void ConfFinal(void);
void ConfXYZ(int);
void Measure(void);
void Equilibra(void);
double Boltzmann(double, double, double, int);
double Pbc(double);
double Error(double,double,int);
double Force(int, int);

#endif
