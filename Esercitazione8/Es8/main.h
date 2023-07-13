#ifndef __EX8__
#define __EX8__

#include <iostream>
#include <fstream>
#include <cmath>
#include <string>
#include <iomanip>
#include "random.h"


Random rnd;

int nstep = 1e6;                      // numero totale di ripetizioni
int nblk = 100;                   // numero di blocchi
int steps = (int)nstep/nblk;          // numero di step in ogni blocco
double L = 3;                   // lato del cubo visitabile / sigma per caso gaussiano
double span;
bool gauss;
int np_rend = 10000;
int mod = (int)nstep/np_rend;
int contastep = 0;

double sgm = 1;
double mu = 1;

bool SA = 0;

int histofill_blk = 100;

// variabili di appoggio per SA
//double T = 10;          // regolabile da input.dat
double beta = 1.0;    // 1/T, avendo dichiarato una T iniziale prima
double db = 1;          // delta beta per ogni passaggio successivo di annealing. Regolabile da input.dat
int step = 0;           // contatore di step
int step_in_beta = 10;  // regolabile da input.dat

double deltaH = 1;
double errH = 0; // basta che sia < deltaH, per entrare nel ciclo
double oldH = 0;
double newH = 0;
double d_mu = 0;
double d_sgm = 0;
double Lmu = 1; // ampiezza iniziale dei salti di mu. Da impostare da input.dat
double Lsgm = 1; // idem per sgm
double errMu = 0.01;
double errSgm = 0.01;

double mean_sgm = 0;
double var_sgm = 0;
double mean_mu = 0;
double var_mu = 0;

bool manualseed = 0;

int amp = 1000;

// variabili di appoggio per metropolis
double H = 0;
double var_prog_H = 0;
double mean_prog_H = 0;

double x0; // posizione iniziale
double y;  // posizione prima del passo
double x;  // posizione dopo il passo

double dx;

double q, A;
int attempted = 0, accepted = 0;

//functions
void MeasureAverageH(int);
void MetroMove(void);
double EvalH(double);
double Psi2(double);
void Input(void);
double min(double, double);
double error(double av, double av2, int n);
#endif
