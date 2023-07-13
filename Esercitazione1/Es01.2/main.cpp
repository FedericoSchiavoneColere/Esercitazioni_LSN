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
    
    if(argc==1){
        cout<<endl<<"Errore: Inserire anche il numero di punti N"<<endl<<endl;
        return 0;
    }
    
    int N = atoi(argv[1]);
    cout<<"N = "<<N<<endl;
    
    
    ofstream outfile_unif1;
    ofstream outfile_unif2;
    ofstream outfile_unif10;
    ofstream outfile_unif100;
    ofstream outfile_expo1;
    ofstream outfile_expo2;
    ofstream outfile_expo10;
    ofstream outfile_expo100;
    ofstream outfile_cauchy_lorentz1;
    ofstream outfile_cauchy_lorentz2;
    ofstream outfile_cauchy_lorentz10;
    ofstream outfile_cauchy_lorentz100;
    outfile_unif1.open("unif1.txt");
    outfile_unif2.open("unif2.txt");
    outfile_unif10.open("unif10.txt");
    outfile_unif100.open("unif100.txt");
    outfile_expo1.open("expo1.txt");
    outfile_expo2.open("expo2.txt");
    outfile_expo10.open("expo10.txt");
    outfile_expo100.open("expo100.txt");
    outfile_cauchy_lorentz1.open("cauchy_lorentz1.txt");
    outfile_cauchy_lorentz2.open("cauchy_lorentz2.txt");
    outfile_cauchy_lorentz10.open("cauchy_lorentz10.txt");
    outfile_cauchy_lorentz100.open("cauchy_lorentz100.txt");
    
    
    
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
    
    vector <double> s_n1;
    vector <double> s_n2;
    vector <double> s_n10;
    vector <double> s_n100;
    vector <double> s_n1expo;
    vector <double> s_n2expo;
    vector <double> s_n10expo;
    vector <double> s_n100expo;
    vector <double> s_n1lor;
    vector <double> s_n2lor;
    vector <double> s_n10lor;
    vector <double> s_n100lor;
    
    //faccio N set di estrazioni, con j diversi.
    
    for(int i=0;i<N;i++){
    
        //j=1
        
        s_n1.push_back(rnd.Rannyu()*6);
        s_n1expo.push_back(rand_expo(rnd.Rannyu(),1));
        s_n1lor.push_back(rand_cauchy_lorentz(rnd.Rannyu(),0,1));
    
        //j=2
        double s_n22=0;
        double s_n22expo=0;
        double s_n22lor=0;
        
        for(int j=0;j<2;j++){
            s_n22+=rnd.Rannyu()*6;
            s_n22expo+=rand_expo(rnd.Rannyu(),1);
            s_n22lor+=rand_cauchy_lorentz(rnd.Rannyu(),0,1);
        }
        s_n2.push_back(s_n22/=2);
        s_n2expo.push_back(s_n22expo/=2);
        s_n2lor.push_back(s_n22lor/=2);
        
        //j=10
        double s_n33=0;
        double s_n33expo=0;
        double s_n33lor=0;
        
        for(int j=0;j<10;j++){
            s_n33+=rnd.Rannyu()*6;
            s_n33expo+=rand_expo(rnd.Rannyu(),1);
            s_n33lor+=rand_cauchy_lorentz(rnd.Rannyu(),0,1);
        }
        s_n10.push_back(s_n33/=10);
        s_n10expo.push_back(s_n33expo/=10);
        s_n10lor.push_back(s_n33lor/=10);
        
        //j=100
        double s_n44=0;
        double s_n44expo=0;
        double s_n44lor=0;
        
        for(int j=0;j<100;j++){
            s_n44+=rnd.Rannyu()*6;
            s_n44expo+=rand_expo(rnd.Rannyu(),1);
            s_n44lor+=rand_cauchy_lorentz(rnd.Rannyu(),0,1);
        }
        s_n100.push_back(s_n44/=100);
        s_n100expo.push_back(s_n44expo/=100);
        s_n100lor.push_back(s_n44lor/=100);
        
        cout<<i<<endl;
    }
    
    //scrivo i file
    
    for(int i=0;i<N;i++){
      
        outfile_unif1<<s_n1[i]<<endl;
        outfile_unif2<<s_n2[i]<<endl;
        outfile_unif10<<s_n10[i]<<endl;
        outfile_unif100<<s_n100[i]<<endl;
        outfile_expo1<<s_n1expo[i]<<endl;
        outfile_expo2<<s_n2expo[i]<<endl;
        outfile_expo10<<s_n10expo[i]<<endl;
        outfile_expo100<<s_n100expo[i]<<endl;
        outfile_cauchy_lorentz1<<s_n1lor[i]<<endl;
        outfile_cauchy_lorentz2<<s_n2lor[i]<<endl;
        outfile_cauchy_lorentz10<<s_n10lor[i]<<endl;
        outfile_cauchy_lorentz100<<s_n100lor[i]<<endl;
        
        cout<<i<<endl;
        
    }
    
    outfile_unif1.close();
    outfile_unif2.close();
    outfile_unif10.close();
    outfile_unif100.close();
    outfile_expo1.close();
    outfile_expo2.close();
    outfile_expo10.close();
    outfile_expo100.close();
    outfile_cauchy_lorentz1.close();
    outfile_cauchy_lorentz2.close();
    outfile_cauchy_lorentz10.close();
    outfile_cauchy_lorentz100.close();
  
    //svuoto i vettori per essere sicuro di liberare la memoria
    
    s_n1.clear();
    s_n2.clear();
    s_n10.clear();
    s_n100.clear();
    s_n1expo.clear();
    s_n2expo.clear();
    s_n10expo.clear();
    s_n100expo.clear();
    s_n1lor.clear();
    s_n2lor.clear();
    s_n10lor.clear();
    s_n100lor.clear();
    
    
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
