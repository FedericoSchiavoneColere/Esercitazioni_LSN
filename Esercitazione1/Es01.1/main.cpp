#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cmath>
#include "random.h"

using namespace std;
 
double errore(double AV1, double AV2, int n);

int main (int argc, char *argv[]){
    
    ofstream outfile_ave;
    ofstream outfile_err;
    ofstream outfile_avedev;
    ofstream outfile_errdev;
    outfile_ave.open("medie.txt");
    outfile_err.open("errori.txt");
    outfile_avedev.open("medie_dev.txt");
    outfile_errdev.open("errori_dev.txt");
    

    int M=10000;
    int N=100;
    int L=M/N;
    
    vector <double> ave;
    vector <double> ave2;
    vector <double> avedev;
    vector <double> ave2dev;
    vector <double> sum_prog;
    vector <double> sum_prog2;
    vector <double> sum_progdev;
    vector <double> sum_prog2dev;
    vector <double> err;
    vector <double> err_prog;
    vector <double> errdev;
    vector <double> err_progdev;

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
    
    cout<<endl;
    
    //calcolo della media a blocchi
    //per ogni iterazione creo un vettore di random vv, è il "blocco"
    //del vettore di random mi interessa solo il valor medio, che scrivo in un altro vettore: ave
    
    for(int i=0;i<L;i++){
        
        vector <double> vv;
        double sum=0;
        double sumdev=0;
        
        for(int j=0; j<N; j++){
            vv.push_back( rnd.Rannyu() );
            sum+=vv[j];
            sumdev+=pow(vv[j]-0.5,2);
        }
        sum=sum/N;
        sumdev=sumdev/N;
        
        ave.push_back(sum);
        ave2.push_back(pow(sum,2));
        avedev.push_back(sumdev);
        ave2dev.push_back(pow(sumdev,2));
    }
    
    rnd.SaveSeed();
    
    //faccio la somma progressiva delle medie dei blocchi, le scrivo in sum_prog
    double prog=0;
    double prog2=0;
    double progdev=0;
    double prog2dev=0;
    for(int i=0; i<L;i++){
        
        prog+=ave[i];
        sum_prog.push_back(prog);
        prog2+=ave2[i];
        sum_prog2.push_back(prog2);
        
        progdev+=avedev[i];
        sum_progdev.push_back(progdev);
        prog2dev+=ave2dev[i];
        sum_prog2dev.push_back(prog2dev);
    }

    //divido le somme progressive per il numero di blocchi sommati, diventano medie
    for(int i=0; i<L;i++){
        sum_prog[i]=sum_prog[i]/(1+i);
        sum_prog2[i]=sum_prog2[i]/(1+i);
        sum_progdev[i]=sum_progdev[i]/(1+i);
        sum_prog2dev[i]=sum_prog2dev[i]/(1+i);
    }
       
    //utilizzo la funzione che calcola gli errori
    for(int i=0;i<L;i++){
        err.push_back( errore(sum_prog[i],sum_prog2[i],i) );
        errdev.push_back( errore(sum_progdev[i],sum_prog2dev[i],i) );
    }
    
  
    for(int i=0;i<L;i++){
        outfile_ave<<sum_prog[i]<<endl;
        outfile_avedev<<sum_progdev[i]<<endl;
        outfile_err<<err[i]<<endl;
        outfile_errdev<<errdev[i]<<endl;
        
    }

    outfile_ave.close();
    outfile_err.close();
    outfile_avedev.close();
    outfile_errdev.close();
    

    ave.clear();
    ave2.clear();
    avedev.clear();
    ave2dev.clear();
    sum_prog.clear();
    sum_prog2.clear();
    sum_progdev.clear();
    sum_prog2dev.clear();
    err.clear();
    err_prog.clear();
    errdev.clear();
    err_progdev.clear();


    return 0;
}

double errore(double AV1, double AV2, int n)
{
    if (n==0) return 0;
        else
    
    return pow((AV2 - pow(AV1,2))/n,0.5);
}


 
