#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cmath>
#include "random.h"

using namespace std;
 
double errore(double AV1, double AV2, int n);
double rand_expo(double unif, double lambda); //funzione che genera numeri esponenzialmente distribuiti
double rand_cauchy_lorentz(double unif, double mu, double tau);// funzione che genera numeri distribuiti come una lorentziana

int main (int argc, char *argv[]){

    ofstream outfile_Iunif;
    ofstream outfile_Iimp;
    ofstream outfile_errunif;
    ofstream outfile_errimp;
    outfile_Iunif.open("int_unif.txt");
    outfile_Iimp.open("int_imp.txt");
    outfile_errunif.open("err_unif.txt");
    outfile_errimp.open("err_imp.txt");
    
    // Generatore di random
    Random rnd;
    int seed[4];
    int p1, p2;
    ifstream Primes("Primes");
    if (Primes.is_open()){
        Primes >> p1 >> p2 ;
    } else cerr << "PROBLEM: Unable to open Primes" << endl;
    Primes.close();
    
    ifstream input("seed.in");
    string property;
    if (input.is_open()){
        while ( !input.eof() ){
            input >> property;
            if( property == "RANDOMSEED" ){
                input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
                rnd.SetRandom(seed,p1,p2);
            }
        }
        input.close();
    } else cerr << "PROBLEM: Unable to open seed.in" << endl;
    
    //
    
    //integrale con metodo monte carlo estraendo uniformemente
    
    int N=100;
    int M=100000;
    double r=0;
    vector <double> I;
    
    for(int i=0; i<N;i++){
        
        double In=0;
        for(int j=0;j<M;j++){
            
            r=rnd.Rannyu();
            In+=(M_PI/2)*cos(M_PI*r/2);
        }
        In/=M;
        I.push_back(In);
        
    }
    
    vector <double> SumProg;
    vector <double> SumProg2;
    vector <double> err;
    
    double sum=0;
    double sum2=0;
    for(int i=0;i<N;i++){
        sum+=I[i];
        sum2+=pow(I[i],2);
        SumProg.push_back(sum);
        SumProg2.push_back(sum2);
    }
    
    for(int i=0;i<N;i++){
        SumProg[i]/=(i+1);
        SumProg2[i]/=(i+1);
    }

    for(int i=0;i<N;i++){
        err.push_back( errore(SumProg[i],SumProg2[i],i) );
    }
    
    for(int i=0;i<N;i++){
        outfile_Iunif<<SumProg[i]<<endl;
        outfile_errunif<<err[i]<<endl;
    }
    
    outfile_Iunif.close();
    outfile_errunif.close();
    
    I.clear();
    SumProg.clear();
    SumProg2.clear();
    err.clear();
    
    
    //cout<<endl<<"******   importance sampling y=-2x+2    ********"<<endl<<endl;
    //integrale con metodo monte carlo estraendo con importance sampling

    
    vector <double> II;
    
    
    for(int i=0; i<N; i++){
        
        double IIn=0;
        
        for(int w=0; w<M; w++){
        
            r=1-pow(1-rnd.Rannyu(),0.5);
            IIn+=(M_PI*0.5)*cos(M_PI*r*0.5)/(-2*r+2);

        }
        IIn/=M;
        II.push_back(IIn);
    }
 


    
    sum=0;
    sum2=0;
    for(int i=0;i<N;i++){
        sum+=II[i];
        sum2+=pow(II[i],2);
        SumProg.push_back(sum);
        SumProg2.push_back(sum2);
    }
    
    for(int i=0;i<N;i++){
        SumProg[i]/=(i+1);
        SumProg2[i]/=(i+1);
    }

    for(int i=0;i<N;i++){
        err.push_back( errore(SumProg[i],SumProg2[i],i) );
    }
    
    for(int i=0;i<N;i++){
        outfile_Iimp<<SumProg[i]<<endl;
        outfile_errimp<<err[i]<<endl;
    }
    
    outfile_Iimp.close();
    outfile_errimp.close();
    
    II.clear();
    SumProg.clear();
    SumProg2.clear();
    err.clear();
    
    
    



    return 0;
}

double errore(double AV1, double AV2, int n)
{
    if (n==0) return 0;
        else
    
    return pow((AV2 - pow(AV1,2))/n,0.5);
}

double rand_expo(double unif, double lambda)
{
    return -1/lambda * log(1-unif); //distribuzione di probabilità invertendo la cumulativa
}

double rand_cauchy_lorentz(double unif, double mu, double tau){
    
    return tau*tan(M_PI*(unif-0.5)) + mu; //distribuzione di probabilità invertendo la cumulativa
    
}
