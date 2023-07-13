#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cmath>
#include "random.h"

using namespace std;
 

int B=100; //numero blocchi
int E=5000; //passi di equilibratura
int N=50000+E; //passi per ogni blocco;

int distrib=1; //0 per unifiorme, 1 per gauss...

double errore(double AV1, double AV2, int n);
double psi_mod2_100(double x, double y, double z);
double psi_mod2_210(double x, double y, double z);
void AC_100(double x, double y, double z,int ii,int &aacc);
void AC_210(double x, double y, double z,int ii,int &aacc);

Random rnd;

vector <double> coord(3,0);
vector < vector <double> > pos_100(N,coord);
vector < vector <double> > pos_210(N,coord);

int main (int argc, char *argv[]){
    
    double delta_100=1.2;
    double delta_210=3.0;
    double sig_100=0.7;
    double sig_210=1.8;
    
    double xx0_100=1;
    double yy0_100=1;
    double zz0_100=1;
    
    double xx0_210=0;
    double yy0_210=0;
    double zz0_210=3;
    
    vector <double> raggio_100;
    vector <double> raggio2_100;
    vector <double> raggio_210;
    vector <double> raggio2_210;
    
    vector <int> accettati_100;
    vector <int> accettati_210;
    
    
    string Name;
    if(distrib==0)Name="Unif";
    else if(distrib==1) Name="Gauss";
    
    
    ofstream outfile_100;
    ofstream outfile_210;
    ofstream outfile_err100;
    ofstream outfile_err210;
    outfile_100.open("./"+Name+"/100.txt");
    outfile_210.open("./"+Name+"/210.txt");
    outfile_err100.open("./"+Name+"/errori_100.txt");
    outfile_err210.open("./"+Name+"/errori_210.txt");

    
    
    
    // Generatore di random
    //Random rnd;
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
    
    
    cout<<endl<<"posizione iniziale100: ("<<xx0_100<<", "<<yy0_100<<", "<<zz0_100<<")"<<endl;
    cout<<"posizione iniziale210: ("<<xx0_210<<", "<<yy0_210<<", "<<zz0_210<<")"<<endl;
    cout<<B<<" blocchi da "<<N<<" step ciascuno, di cui i primi "<<E<<" passi sono per equilibrare"<<endl;
    
    
    if(distrib<0.5){
        cout<<"passi da distribuzione uniforme, delta_100 = ± "<<delta_100<<", delta_210 = ± "<<delta_210<<endl<<endl;
    }else{
        cout<<"passi da distribuzione gaussiana, sigma_100 = "<<sig_100<<", sigma_210 = "<<sig_210<<endl<<endl;
    }
    
    
    for(int j=0;j<B;j++){
        
        int acc_100=0;
        int acc_210=0;
        double R_100=0;
        double R_210=0;
        
        pos_100[0][0]=xx0_100;
        pos_100[0][1]=yy0_100;
        pos_100[0][2]=zz0_100;
        
        pos_210[0][0]=xx0_210;
        pos_210[0][1]=yy0_210;
        pos_210[0][2]=zz0_210;
        
        cout<<"eseguo ciclo: "<<j<<endl;
        
        for(int i=1; i<N; i++){
        
            double x_100;
            double y_100;
            double z_100;
            
            double x_210;
            double y_210;
            double z_210;
            
            if(distrib<0.5){
                x_100=pos_100[i-1][0]+rnd.Rannyu(-delta_100,delta_100);
                y_100=pos_100[i-1][1]+rnd.Rannyu(-delta_100,delta_100);
                z_100=pos_100[i-1][2]+rnd.Rannyu(-delta_100,delta_100);
                
                x_210=pos_210[i-1][0]+rnd.Rannyu(-delta_210,delta_210);
                y_210=pos_210[i-1][1]+rnd.Rannyu(-delta_210,delta_210);
                z_210=pos_210[i-1][2]+rnd.Rannyu(-delta_210,delta_210);
            }else{
                x_100=pos_100[i-1][0]+rnd.Gauss(0,sig_100);
                y_100=pos_100[i-1][1]+rnd.Gauss(0,sig_100);
                z_100=pos_100[i-1][2]+rnd.Gauss(0,sig_100);
                
                x_210=pos_210[i-1][0]+rnd.Gauss(0,sig_210);
                y_210=pos_210[i-1][1]+rnd.Gauss(0,sig_210);
                z_210=pos_210[i-1][2]+rnd.Gauss(0,sig_210);
            }
            
             
            if(i==E){       //aspetto che sia equilibrato
                acc_100=0;
                acc_210=0;
            }
            
            AC_100(x_100,y_100,z_100,i,acc_100);
            AC_210(x_210,y_210,z_210,i,acc_210);
            
            if(i>=E){  //aspetto che sia equilibrato
                double ri_100=0;
                double ri_210=0;
                ri_100=pow( pow(pos_100[i][0],2)+pow(pos_100[i][1],2)+pow(pos_100[i][2],2) ,0.5);
                ri_210=pow( pow(pos_210[i][0],2)+pow(pos_210[i][1],2)+pow(pos_210[i][2],2) ,0.5);
                
                R_100+=ri_100;
                R_210+=ri_210;
            }
        }
        accettati_100.push_back(acc_100);
        accettati_210.push_back(acc_210);
        R_100/=(N-E);
        R_210/=(N-E);
        raggio_100.push_back(R_100);
        raggio_210.push_back(R_210);
        raggio2_100.push_back(pow(R_100,2));
        raggio2_210.push_back(pow(R_210,2));
    }
    cout<<endl;
    
    vector <double> Rbl_100;
    vector <double> Rbl2_100;
    vector <double> Rbl_210;
    vector <double> Rbl2_210;
    vector <double> err_100;
    vector <double> err_210;
    
    Rbl_100.push_back(raggio_100[0]);
    Rbl_210.push_back(raggio_210[0]);
    Rbl2_100.push_back(raggio2_100[0]);
    Rbl2_210.push_back(raggio2_210[0]);
    
    for(int i=1; i<B; i++){
        Rbl_100.push_back(Rbl_100[i-1]+raggio_100[i]);
        Rbl_210.push_back(Rbl_210[i-1]+raggio_210[i]);
        Rbl2_100.push_back(Rbl2_100[i-1]+raggio2_100[i]);
        Rbl2_210.push_back(Rbl2_210[i-1]+raggio2_210[i]);
    }
    double ac_100=0;
    double ac_210=0;
    for(int i=0;i<B;i++){
        Rbl_100[i]/=i+1;
        Rbl_210[i]/=i+1;
        Rbl2_100[i]/=i+1;
        Rbl2_210[i]/=i+1;
        ac_100+=accettati_100[i];
        ac_210+=accettati_210[i];
    }
    
    cout<<endl<<"media a blocchi"<<endl;
    cout<<"stato 100       stato 210 "<<endl<<endl;
    for(int i=0;i<B;i++){
        cout<<fixed<<setprecision(6)<<Rbl_100[i]<<"         "<<Rbl_210[i]<<endl;
        outfile_100<<Rbl_100[i]<<endl;
        outfile_210<<Rbl_210[i]<<endl;
    }
    
    for(int i=0; i<B;i++){
        err_100.push_back(errore(Rbl_100[i],Rbl2_100[i],i));
        err_210.push_back(errore(Rbl_210[i],Rbl2_210[i],i));
        outfile_err100<<err_100[i]<<endl;
        outfile_err210<<err_210[i]<<endl;
    }
    
    ac_100/=B;
    ac_210/=B;
    cout<<endl<<"tasso accettazione medio_100: "<<ac_100/(N-E)<<endl;
    cout<<"tasso accettazione medio_210: "<<ac_210/(N-E)<<endl<<endl;
    
    outfile_100.close();
    outfile_210.close();
    outfile_err100.close();
    outfile_err210.close();
    
    
    return 0;
}



double errore(double AV1, double AV2, int n)
{
    if (n==0) return 0;
        else
    
    return pow((AV2 - pow(AV1,2))/n,0.5);
}

double psi_mod2_100(double x, double y, double z)
{
    return 1/M_PI * exp(-2*pow( pow(x,2)+pow(y,2)+pow(z,2) ,0.5));
}

double psi_mod2_210(double x, double y, double z)
{
    return 1/(32*M_PI) * pow(z,2) * exp(-pow( pow(x,2)+pow(y,2)+pow(z,2) ,0.5));
}

void AC_100(double x, double y, double z,int ii,int &aacc)
{
    double p=0;
    vector <double> coord_prov(3,0);
    coord_prov[0]=x;
    coord_prov[1]=y;
    coord_prov[2]=z;
    
    if(ii>0){
        p=psi_mod2_100(x,y,z)/psi_mod2_100(pos_100[ii-1][0],pos_100[ii-1][1],pos_100[ii-1][2]);
    if(p>=1){
        aacc++;
        pos_100[ii]=coord_prov;
    }else{
        double aa=rnd.Rannyu();
        
        if(aa<=p){
            aacc++;
            pos_100[ii]=coord_prov;
        }else{
            coord_prov[0]=pos_100[ii-1][0];
            coord_prov[1]=pos_100[ii-1][1];
            coord_prov[2]=pos_100[ii-1][2];

            pos_100[ii]=coord_prov;
            }
        }
    }
}

void AC_210(double x, double y, double z,int ii,int &aacc)
{
    double p=0;
    vector <double> coord_prov(3,0);
    coord_prov[0]=x;
    coord_prov[1]=y;
    coord_prov[2]=z;
    
    if(ii>0){
        p=psi_mod2_210(x,y,z)/psi_mod2_210(pos_210[ii-1][0],pos_210[ii-1][1],pos_210[ii-1][2]);
    if(p>=1){
        aacc++;
        pos_210[ii]=coord_prov;
    }else{
        double aa=rnd.Rannyu();
        
        if(aa<=p){
            aacc++;
            pos_210[ii]=coord_prov;
        }else{
            coord_prov[0]=pos_210[ii-1][0];
            coord_prov[1]=pos_210[ii-1][1];
            coord_prov[2]=pos_210[ii-1][2];

            pos_210[ii]=coord_prov;
            }
        }
    }
}
