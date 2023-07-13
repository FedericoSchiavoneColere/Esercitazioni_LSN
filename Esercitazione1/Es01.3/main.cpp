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
double rand_cos1();
double rand_cos2();
int main (int argc, char *argv[]){
    
    
    // Generatore di random uniformi in (0;1)
   
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
    
    int N=10000000; //lanci per ogni ciclo
    int NN=100;  //numero di clicli
    double L=0.8; //lunghezza barretta
                  //la distanza tra le "righe" è normalizzata a 1
    
    vector <double> Pi;
    vector <double> Pi2;
    vector <double> SumProg;
    vector <double> SumProg2;
    vector <double> err;
    
//per ogni iterazione di NN c'è un ciclo N (blocco) in cui si calcola pi N volte
//il valore medio di ogni ciclo N viene immagazzinato in un vettore nel ciclo NN
    
    for(int j=0;j<NN;j++){
        
        vector <double> w;
        
        int m=0;
        
        for(int i=0;i<N;i++){
            w.push_back(rnd.Rannyu());
        }
        
        for(int i=0;i<N;i++){
            if(w[i]<=(L/2)){
                double t=0;
                t=rand_cos1();
                //t=rand_cos2();
                if(w[i]<=((L/2)*t))m++;
            }
            if(w[i]>=(1-(L/2))){
                double t=0;
                t=rand_cos1();
                //t=rand_cos2();
                if(w[i]+((L/2)*t)>=1)m++;
            }
        }
        
        Pi.push_back((2*L*N)/m);
        Pi2.push_back(pow(Pi[j],2));
        cout<<j<<endl;
        }

//passo dal vettore con le medie a vettore con le medie progressive a blocchi
    
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
    
//calcolo l'errore progressivo sui blocchi
    
    for(int i=0;i<NN;i++){
        err.push_back( errore(SumProg[i],SumProg2[i],i) );
    }

// scrivo i dati in due file
    
    for(int i=0;i<NN;i++){
        outfile_ave<<SumProg[i]<<endl;
        outfile_err<<err[i]<<endl;
    }
    
    outfile_ave.close();
    outfile_err.close();
        
    Pi.clear();
    Pi2.clear();
    SumProg.clear();
    SumProg2.clear();
    err.clear();
    
    
// istogramma distribuzione coseno con i due metodi
    
    ofstream outfile_cos1;
    ofstream outfile_cos2;
    outfile_cos1.open("cos1.txt");
    outfile_cos2.open("cos2.txt");
    
    for(int i=0;i<10000;i++){
        
        outfile_cos1<<rand_cos1()<<endl;
        outfile_cos2<<rand_cos2()<<endl;
        
    }
    outfile_cos1.close();
    outfile_cos2.close();
     
//
    
    return 0;
}

double errore(double AV1, double AV2, int n)
{
    if (n==0) return 0;
        else
    
    return pow((AV2 - pow(AV1,2))/n,0.5);
}

double rand_cos1()
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
 
double rand_cos2()
{
    double a=rnd.Rannyu();
    a=a*(2)-1;
    a=a*M_PI/2;

    return cos(a);
}

