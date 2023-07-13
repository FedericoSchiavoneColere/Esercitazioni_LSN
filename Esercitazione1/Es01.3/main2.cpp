#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cmath>
#include <math.h>
#include "random.h"

using namespace std;

Random rnd;
double errore(double AV1, double AV2, int n);
double rand_coss1();
double rand_coss2();
double rand_coss3();
int main (int argc, char *argv[]){
    
    
    // Generatore di random
   
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
    
    cout<<endl;
    
    ofstream outfile_ave;
    ofstream outfile_err;
    outfile_ave.open("medie.txt");
    outfile_err.open("errori.txt");
    
    int N=100000;
    int NN=1000;
    double L=0.8;
    
    vector <double> Pi;
    vector <double> Pi2;
    vector <double> SumProg;
    vector <double> SumProg2;
    vector <double> err;
    
    
    for(int j=0;j<NN;j++){
        
        vector <double> w;
        
        int m=0;
        
        for(int i=0;i<N;i++){
            w.push_back(rnd.Rannyu());
        }
        
        for(int i=0;i<N;i++){
            if(w[i]<=(L/2)){
                double t=0;
                //t=cos(rnd.Rannyu()*M_PI/2);
                t=rand_coss1();
                //t=rand_coss2();
                //t=rand_coss3();
                if(w[i]<=((L/2)*t))m++;
            }
            if(w[i]>=(1-(L/2))){
                double t=0;
                //t=cos(rnd.Rannyu()*M_PI/2);
                t=rand_coss1();
                //t=rand_coss2();
                //t=rand_coss3();
                if(w[i]+((L/2)*t)>=1)m++;
            }
        }
        
        Pi.push_back((2*L*N)/m);
        Pi2.push_back(pow(Pi[j],2));
        cout<<j<<endl;
        }
    
    double prog=0;
    double prog2=0;
    
    for(int i=0;i<NN;i++){
        prog+=Pi[i];
        prog2+=Pi2[i];
        SumProg.push_back(prog);
        SumProg2.push_back(prog2);
    }

    for(int i=0;i<NN;i++){
        SumProg[i]/=(i+1.);
        SumProg2[i]/=(i+1.);
        
        cout<<SumProg[i]<<endl;
    }
    for(int i=0;i<NN;i++){
        err.push_back( errore(SumProg[i],SumProg2[i],i) );
    }
    
    for(int i=0;i<NN;i++){
        outfile_ave<<SumProg[i]<<endl;
        outfile_err<<err[i]<<endl;
    }
    
    outfile_ave.close();
    outfile_err.close();
    
    
// istogramma distribuzione coseno con vari metodi
    
    ofstream outfile_cos;
    ofstream outfile_coss1;
    ofstream outfile_coss2;
    ofstream outfile_coss3;
    outfile_cos.open("cos.txt");
    outfile_coss1.open("coss1.txt");
    outfile_coss2.open("coss2.txt");
    outfile_coss3.open("coss3.txt");
    
    for(int i=0;i<10000;i++){
        
        outfile_cos<<cos(rnd.Rannyu()*M_PI/2)<<endl;
        outfile_coss1<<rand_coss1()<<endl;
        outfile_coss2<<rand_coss2()<<endl;
        outfile_coss3<<rand_coss3()<<endl;
    }
    outfile_cos.close();
    outfile_coss1.close();
    outfile_coss2.close();
    outfile_coss3.close();
//
     
     
    
    return 0;
}

double errore(double AV1, double AV2, int n)
{
    if (n==0) return 0;
        else
    
    return pow((AV2 - pow(AV1,2))/n,0.5);
}

double rand_coss1()
{
    
    double a=0;
    double b=0;
    double c=0;
    
    do{
        a=rnd.Rannyu();
        b=rnd.Rannyu();
        c=pow(pow(a,2)+pow(b,2),0.5);
    }while(c!=0&&c>1);
    
    if(c>0&&c<1) return a/c;
}
 
double rand_coss2()
{
    double a=rnd.Rannyu();
    a=a*(2)-1;
    a=a*M_PI/2;

    return cos(a);
}

double rand_coss3()
{
    double a=rnd.Rannyu();
 
    return asin(1-2*a);
}
