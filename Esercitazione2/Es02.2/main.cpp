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
double distanza(vector <double> part);
double st_dev(vector<double> dat, double av,int N);


int main (int argc, char *argv[]){
    
    ofstream outfile_RW;
    ofstream outfile_err;
    ofstream outfile_RWsfer;
    ofstream outfile_errsfer;
    outfile_RW.open("RW.txt");
    outfile_err.open("err.txt");
    outfile_RWsfer.open("RWsfer.txt");
    outfile_errsfer.open("errsfer.txt");
    
    
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
    
    
    //random walk sugli assi
    
    int N=10000;
    int T=400;
    
    cout<<endl<<"Eseguo il programma con "<<N<<" particelle."<<endl;
    cout<<"Si eseguono "<<T<<" passi."<<endl<<endl;
    
    
    
    //inizializzo i vettori cord alla posizione iniziale (0;0;0)
    
    vector <float> cord;
    vector <vector<float>> pos;
    
    for(int i=0; i<N;i++){
        
        for(int x=0; x<3;x++){
            cord.push_back(0);
        }
        pos.push_back(cord);
    }
    
    
    cout<<"Simulazione con passi a=1 lungo la direzione di un asse"<<endl <<endl;
    
    
    //faccio i T passi e scrivo la posizione per ogni RW ogni 10 passi
    
    
    vector <vector<double>> pos_10(N);
    vector <vector<double>> pos_20(N);
    vector <vector<double>> pos_30(N);
    vector <vector<double>> pos_40(N);
    vector <vector<double>> pos_50(N);
    vector <vector<double>> pos_60(N);
    vector <vector<double>> pos_100(N);
    vector <vector<double>> pos_200(N);
    vector <vector<double>> pos_300(N);
    vector <vector<double>> pos_400(N);
    
    for(int p=0;p<N;p++){
        
        
        for(int t=0; t<T/10;t++){
            
            for(int i=1; i<=10;i++){
                
                int dir=0;
                dir=int(rnd.Rannyu()*3+1);
                
                if(dir==1){
                    int verso=0;
                    verso=int(rnd.Rannyu()*2);
                    if(verso==0){
                        pos[p][0]--;
                    }
                    if(verso==1){
                        pos[p][0]++;
                    }
                }
                if(dir==2){
                    int verso=0;
                    verso=int(rnd.Rannyu()*2);
                    if(verso==0){
                        pos[p][1]--;
                    }
                    if(verso==1){
                        pos[p][1]++;
                    }
                }
                if(dir==3){
                    int verso=0;
                    verso=int(rnd.Rannyu()*2);
                    if(verso==0){
                        pos[p][2]--;
                    }
                    if(verso==1){
                        pos[p][2]++;
                    }
                }
                //  cout<<"RW numero: "<<p+1<<" posizione dopo "<<i+t*10<<" passi ";
                // cout<<pos[p][0]<<"  "<<pos[p][1]<<"  "<<pos[p][2]<<endl;
                
            }
            
            if(t==0){
                pos_10[p].push_back(pos[p][0]);
                pos_10[p].push_back(pos[p][1]);
                pos_10[p].push_back(pos[p][2]);
                //cout<<endl<<"particella "<<p+1<<": posizione dopo 10 passi: ";
                //cout<<pos_10[p][0]<<"  "<<pos_10[p][1]<<"  "<<pos_10[p][2]<<endl<<endl;
            }
            if(t==1){
                pos_20[p].push_back(pos[p][0]);
                pos_20[p].push_back(pos[p][1]);
                pos_20[p].push_back(pos[p][2]);
                //cout<<endl<<"particella "<<p+1<<": posizione dopo 20 passi: ";
                //cout<<pos_20[p][0]<<"  "<<pos_20[p][1]<<"  "<<pos_20[p][2]<<endl<<endl;
            }
            if(t==2){
                pos_30[p].push_back(pos[p][0]);
                pos_30[p].push_back(pos[p][1]);
                pos_30[p].push_back(pos[p][2]);
                //cout<<endl<<"particella "<<p+1<<": posizione dopo 30 passi: ";
                //cout<<pos_30[p][0]<<"  "<<pos_30[p][1]<<"  "<<pos_30[p][2]<<endl<<endl;
            }
            if(t==3){
                pos_40[p].push_back(pos[p][0]);
                pos_40[p].push_back(pos[p][1]);
                pos_40[p].push_back(pos[p][2]);
                //cout<<endl<<"particella "<<p+1<<": posizione dopo 40 passi: ";
                //cout<<pos_40[p][0]<<"  "<<pos_40[p][1]<<"  "<<pos_40[p][2]<<endl<<endl;
            }
            if(t==4){
                pos_50[p].push_back(pos[p][0]);
                pos_50[p].push_back(pos[p][1]);
                pos_50[p].push_back(pos[p][2]);
                //cout<<endl<<"particella "<<p+1<<": posizione dopo 50 passi: ";
                //cout<<pos_50[p][0]<<"  "<<pos_50[p][1]<<"  "<<pos_50[p][2]<<endl<<endl;
            }
            if(t==5){
                pos_60[p].push_back(pos[p][0]);
                pos_60[p].push_back(pos[p][1]);
                pos_60[p].push_back(pos[p][2]);
                //cout<<endl<<"particella "<<p+1<<": posizione dopo 60 passi: ";
                //cout<<pos_60[p][0]<<"  "<<pos_60[p][1]<<"  "<<pos_60[p][2]<<endl<<endl;
            }
            if(t==9){
                pos_100[p].push_back(pos[p][0]);
                pos_100[p].push_back(pos[p][1]);
                pos_100[p].push_back(pos[p][2]);
            }
            if(t==19){
                pos_200[p].push_back(pos[p][0]);
                pos_200[p].push_back(pos[p][1]);
                pos_200[p].push_back(pos[p][2]);
            }
            if(t==29){
                pos_300[p].push_back(pos[p][0]);
                pos_300[p].push_back(pos[p][1]);
                pos_300[p].push_back(pos[p][2]);
            }
            if(t==39){
                pos_400[p].push_back(pos[p][0]);
                pos_400[p].push_back(pos[p][1]);
                pos_400[p].push_back(pos[p][2]);
            }
        }
    }
    
    /*
     for(int i=0;i<N;i++){
     
     cout<<pos_10[i][0]<<"  "<<pos_10[i][1]<<"  "<<pos_10[i][2]<<endl;
     cout<<"distanza: "<<distanza(pos_10[i])<<endl;
     
     }
     */
    double dist_media10=0;
    double dist_media20=0;
    double dist_media30=0;
    double dist_media40=0;
    double dist_media50=0;
    double dist_media60=0;
    double dist_media100=0;
    double dist_media200=0;
    double dist_media300=0;
    double dist_media400=0;
    
    double dist2_media10=0;
    double dist2_media20=0;
    double dist2_media30=0;
    double dist2_media40=0;
    double dist2_media50=0;
    double dist2_media60=0;
    double dist2_media100=0;
    double dist2_media200=0;
    double dist2_media300=0;
    double dist2_media400=0;
    
    for(int i=0;i<N;i++){
        dist_media10+=distanza(pos_10[i]);
        dist_media20+=distanza(pos_20[i]);
        dist_media30+=distanza(pos_30[i]);
        dist_media40+=distanza(pos_40[i]);
        dist_media50+=distanza(pos_50[i]);
        dist_media60+=distanza(pos_60[i]);
        dist_media100+=distanza(pos_100[i]);
        dist_media200+=distanza(pos_200[i]);
        dist_media300+=distanza(pos_300[i]);
        dist_media400+=distanza(pos_400[i]);
        
        dist2_media10+=pow(distanza(pos_10[i]),2);
        dist2_media20+=pow(distanza(pos_20[i]),2);
        dist2_media30+=pow(distanza(pos_30[i]),2);
        dist2_media40+=pow(distanza(pos_40[i]),2);
        dist2_media50+=pow(distanza(pos_50[i]),2);
        dist2_media60+=pow(distanza(pos_60[i]),2);
        dist2_media100+=pow(distanza(pos_100[i]),2);
        dist2_media200+=pow(distanza(pos_200[i]),2);
        dist2_media300+=pow(distanza(pos_300[i]),2);
        dist2_media400+=pow(distanza(pos_400[i]),2);
    }
    dist_media10/=N;
    dist_media20/=N;
    dist_media30/=N;
    dist_media40/=N;
    dist_media50/=N;
    dist_media60/=N;
    dist_media100/=N;
    dist_media200/=N;
    dist_media300/=N;
    dist_media400/=N;
    
    
    
    dist2_media10/=N;
    dist2_media20/=N;
    dist2_media30/=N;
    dist2_media40/=N;
    dist2_media50/=N;
    dist2_media60/=N;
    dist2_media100/=N;
    dist2_media200/=N;
    dist2_media300/=N;
    dist2_media400/=N;
    
    dist2_media10=pow(dist2_media10,0.5);
    dist2_media20=pow(dist2_media20,0.5);
    dist2_media30=pow(dist2_media30,0.5);
    dist2_media40=pow(dist2_media40,0.5);
    dist2_media50=pow(dist2_media50,0.5);
    dist2_media60=pow(dist2_media60,0.5);
    dist2_media100=pow(dist2_media100,0.5);
    dist2_media200=pow(dist2_media200,0.5);
    dist2_media300=pow(dist2_media300,0.5);
    dist2_media400=pow(dist2_media400,0.5);
    
    
    
    vector <double> dist2a_media={dist2_media10,dist2_media20,dist2_media30,dist2_media40,dist2_media50,dist2_media60,dist2_media100,dist2_media200,dist2_media300,dist2_media400};
    
    for(int i=0;i<10;i++){
        outfile_RW<<dist2a_media[i]<<endl;
    }
    
    cout<<"Distanza media dopo 10 passi: "<<dist_media10<<endl;
    cout<<"Distanza media dopo 20 passi: "<<dist_media20<<endl;
    cout<<"Distanza media dopo 30 passi: "<<dist_media30<<endl;
    cout<<"Distanza media dopo 40 passi: "<<dist_media40<<endl;
    cout<<"Distanza media dopo 50 passi: "<<dist_media50<<endl;
    cout<<"Distanza media dopo 60 passi: "<<dist_media60<<endl;
    cout<<"Distanza media dopo 100 passi: "<<dist_media100<<endl;
    cout<<"Distanza media dopo 200 passi: "<<dist_media200<<endl;
    cout<<"Distanza media dopo 300 passi: "<<dist_media300<<endl;
    cout<<"Distanza media dopo 400 passi: "<<dist_media400<<endl;
    cout<<endl;
    cout<<"Distanza media2 dopo 10 passi: "<<dist2_media10<<endl;
    cout<<"Distanza media2 dopo 20 passi: "<<dist2_media20<<endl;
    cout<<"Distanza media2 dopo 30 passi: "<<dist2_media30<<endl;
    cout<<"Distanza media2 dopo 40 passi: "<<dist2_media40<<endl;
    cout<<"Distanza media2 dopo 50 passi: "<<dist2_media50<<endl;
    cout<<"Distanza media2 dopo 60 passi: "<<dist2_media60<<endl;
    cout<<"Distanza media2 dopo 100 passi: "<<dist2_media100<<endl;
    cout<<"Distanza media2 dopo 200 passi: "<<dist2_media200<<endl;
    cout<<"Distanza media2 dopo 300 passi: "<<dist2_media300<<endl;
    cout<<"Distanza media2 dopo 400 passi: "<<dist2_media400<<endl;
    
    
    //random walk direzione random sferica
    //&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
    
    
    vector <float> cord_sfer;
    vector <vector<float>> pos_sfer;
    
    //inizializzo i vettori cord alla posizione iniziale (0;0;0)
    
    
    for(int i=0; i<N;i++){
        
        for(int x=0; x<3;x++){
            cord_sfer.push_back(0);
        }
        pos_sfer.push_back(cord_sfer);
    }
    
    cout<<endl<<endl<<"Simulazione con passi a=1 lungo una direzione random"<<endl<<endl;
    
    vector <vector<double>> pos_10sfer(N);
    vector <vector<double>> pos_20sfer(N);
    vector <vector<double>> pos_30sfer(N);
    vector <vector<double>> pos_40sfer(N);
    vector <vector<double>> pos_50sfer(N);
    vector <vector<double>> pos_60sfer(N);
    
    for(int p=0;p<N;p++){
        
        double x=0;
        double y=0;
        double z=0;
        
        for(int t=0; t<T/10;t++){
            
            for(int i=0; i<10;i++){
                
                double th=rnd.Rannyu()*M_PI;
                double phi=rnd.Rannyu()*2*M_PI;
                
                x+=cos(phi)*sin(th);
                y+=sin(phi)*sin(th);
                z+=cos(th);
                
                
                pos_sfer[p][0]=x;
                pos_sfer[p][1]=y;
                pos_sfer[p][2]=z;
                
                //cout<<"RW numero: "<<p+1<<" posizione dopo "<<i+t*10<<" passi ";
                //cout<<pos_sfer[p][0]<<"  "<<pos_sfer[p][1]<<"  "<<pos_sfer[p][2]<<endl;
                
            }
            
            if(t==0){
                pos_10sfer[p].push_back(pos_sfer[p][0]);
                pos_10sfer[p].push_back(pos_sfer[p][1]);
                pos_10sfer[p].push_back(pos_sfer[p][2]);
                //cout<<endl<<"particella "<<p+1<<": posizione dopo 10 passi: ";
                //cout<<pos_10sfer[p][0]<<"  "<<pos_10sfer[p][1]<<"  "<<pos_10sfer[p][2]<<endl<<endl;
            }
            if(t==1){
                pos_20sfer[p].push_back(pos_sfer[p][0]);
                pos_20sfer[p].push_back(pos_sfer[p][1]);
                pos_20sfer[p].push_back(pos_sfer[p][2]);
                //cout<<endl<<"particella "<<p+1<<": posizione dopo 20 passi: ";
                //cout<<pos_20sfer[p][0]<<"  "<<pos_20sfer[p][1]<<"  "<<pos_20sfer[p][2]<<endl<<endl;
            }
            if(t==2){
                pos_30sfer[p].push_back(pos_sfer[p][0]);
                pos_30sfer[p].push_back(pos_sfer[p][1]);
                pos_30sfer[p].push_back(pos_sfer[p][2]);
                //cout<<endl<<"particella "<<p+1<<": posizione dopo 30 passi: ";
                //cout<<pos_30sfer[p][0]<<"  "<<pos_30sfer[p][1]<<"  "<<pos_30sfer[p][2]<<endl<<endl;
            }
            if(t==3){
                pos_40sfer[p].push_back(pos_sfer[p][0]);
                pos_40sfer[p].push_back(pos_sfer[p][1]);
                pos_40sfer[p].push_back(pos_sfer[p][2]);
                //cout<<endl<<"particella "<<p+1<<": posizione dopo 40 passi: ";
                //cout<<pos_40sfer[p][0]<<"  "<<pos_40sfer[p][1]<<"  "<<pos_40sfer[p][2]<<endl<<endl;
            }
            if(t==4){
                pos_50sfer[p].push_back(pos_sfer[p][0]);
                pos_50sfer[p].push_back(pos_sfer[p][1]);
                pos_50sfer[p].push_back(pos_sfer[p][2]);
                //cout<<endl<<"particella "<<p+1<<": posizione dopo 50 passi: ";
                //cout<<pos_50sfer[p][0]<<"  "<<pos_50sfer[p][1]<<"  "<<pos_50sfer[p][2]<<endl<<endl;
            }
            if(t==5){
                pos_60sfer[p].push_back(pos_sfer[p][0]);
                pos_60sfer[p].push_back(pos_sfer[p][1]);
                pos_60sfer[p].push_back(pos_sfer[p][2]);
                //cout<<endl<<"particella "<<p+1<<": posizione dopo 60 passi: ";
                //cout<<pos_60sfer[p][0]<<"  "<<pos_60sfer[p][1]<<"  "<<pos_60sfer[p][2]<<endl<<endl;
            }
            
        }
    }
    
    
    
    double dist_media10s=0;
    double dist_media20s=0;
    double dist_media30s=0;
    double dist_media40s=0;
    double dist_media50s=0;
    double dist_media60s=0;
    
    double dist2_media10s=0;
    double dist2_media20s=0;
    double dist2_media30s=0;
    double dist2_media40s=0;
    double dist2_media50s=0;
    double dist2_media60s=0;
    
    for(int i=0;i<N;i++){
        dist_media10s+=distanza(pos_10sfer[i]);
        dist_media20s+=distanza(pos_20sfer[i]);
        dist_media30s+=distanza(pos_30sfer[i]);
        dist_media40s+=distanza(pos_40sfer[i]);
        dist_media50s+=distanza(pos_50sfer[i]);
        dist_media60s+=distanza(pos_60sfer[i]);
        
        dist2_media10s+=pow(distanza(pos_10sfer[i]),2);
        dist2_media20s+=pow(distanza(pos_20sfer[i]),2);
        dist2_media30s+=pow(distanza(pos_30sfer[i]),2);
        dist2_media40s+=pow(distanza(pos_40sfer[i]),2);
        dist2_media50s+=pow(distanza(pos_50sfer[i]),2);
        dist2_media60s+=pow(distanza(pos_60sfer[i]),2);
    }
    dist_media10s/=N;
    dist_media20s/=N;
    dist_media30s/=N;
    dist_media40s/=N;
    dist_media50s/=N;
    dist_media60s/=N;
    
    
    dist2_media10s/=N;
    dist2_media20s/=N;
    dist2_media30s/=N;
    dist2_media40s/=N;
    dist2_media50s/=N;
    dist2_media60s/=N;
    
    dist2_media10s=pow(dist2_media10s,0.5);
    dist2_media20s=pow(dist2_media20s,0.5);
    dist2_media30s=pow(dist2_media30s,0.5);
    dist2_media40s=pow(dist2_media40s,0.5);
    dist2_media50s=pow(dist2_media50s,0.5);
    dist2_media60s=pow(dist2_media60s,0.5);
    
    vector <double> dist2s_media={dist2_media10s,dist2_media20s,dist2_media30s,dist2_media40s,dist2_media50s,dist2_media60s};
    
    for(int i=0;i<6;i++){
        outfile_RWsfer<<dist2s_media[i]<<endl;
    }
    
    
    cout<<"Distanza media dopo 10 passi: "<<dist_media10s<<endl;
    cout<<"Distanza media dopo 20 passi: "<<dist_media20s<<endl;
    cout<<"Distanza media dopo 30 passi: "<<dist_media30s<<endl;
    cout<<"Distanza media dopo 40 passi: "<<dist_media40s<<endl;
    cout<<"Distanza media dopo 50 passi: "<<dist_media50s<<endl;
    cout<<"Distanza media dopo 60 passi: "<<dist_media60s<<endl;
    cout<<endl;
    cout<<"Distanza media2 dopo 10 passi: "<<dist2_media10s<<endl;
    cout<<"Distanza media2 dopo 20 passi: "<<dist2_media20s<<endl;
    cout<<"Distanza media2 dopo 30 passi: "<<dist2_media30s<<endl;
    cout<<"Distanza media2 dopo 40 passi: "<<dist2_media40s<<endl;
    cout<<"Distanza media2 dopo 50 passi: "<<dist2_media50s<<endl;
    cout<<"Distanza media2 dopo 60 passi: "<<dist2_media60s<<endl;
    cout<<endl;
    
    
    //calcolo degli errori
    //&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
    
    
    int B=100; //quante particelle in ogni blocco
    int L=N/B; //numero di blocchi
    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    //direzione lungo gli assi
    
    
    
    vector <double> pos10aB;  //distanze medie su un blocco di 100 particelle
    vector <double> pos20aB;
    vector <double> pos30aB;
    vector <double> pos40aB;
    vector <double> pos50aB;
    vector <double> pos60aB;
    vector <double> pos100aB;
    vector <double> pos200aB;
    vector <double> pos300aB;
    vector <double> pos400aB;
    
    vector <double> pos10aB2; //distanze al quadrato medie sotto radice
    vector <double> pos20aB2; //su un blocco di 100 particelle
    vector <double> pos30aB2;
    vector <double> pos40aB2;
    vector <double> pos50aB2;
    vector <double> pos60aB2;
    vector <double> pos100aB2;
    vector <double> pos200aB2;
    vector <double> pos300aB2;
    vector <double> pos400aB2;
    
    for(int i=0; i<L;i++){
        
        double pos10assiB=0;
        double pos20assiB=0;
        double pos30assiB=0;
        double pos40assiB=0;
        double pos50assiB=0;
        double pos60assiB=0;
        double pos100assiB=0;
        double pos200assiB=0;
        double pos300assiB=0;
        double pos400assiB=0;
        
        double pos10assiB2=0;
        double pos20assiB2=0;
        double pos30assiB2=0;
        double pos40assiB2=0;
        double pos50assiB2=0;
        double pos60assiB2=0;
        double pos100assiB2=0;
        double pos200assiB2=0;
        double pos300assiB2=0;
        double pos400assiB2=0;
        
        for(int j=0;j<B;j++){
            
            pos10assiB+=distanza(pos_10[i]);
            pos20assiB+=distanza(pos_20[i]);
            pos30assiB+=distanza(pos_30[i]);
            pos40assiB+=distanza(pos_40[i]);
            pos50assiB+=distanza(pos_50[i]);
            pos60assiB+=distanza(pos_60[i]);
            pos100assiB+=distanza(pos_100[i]);
            pos200assiB+=distanza(pos_200[i]);
            pos300assiB+=distanza(pos_300[i]);
            pos400assiB+=distanza(pos_400[i]);
            
            pos10assiB2+=pow(distanza(pos_10[i]),2);
            pos20assiB2+=pow(distanza(pos_20[i]),2);
            pos30assiB2+=pow(distanza(pos_30[i]),2);
            pos40assiB2+=pow(distanza(pos_40[i]),2);
            pos50assiB2+=pow(distanza(pos_50[i]),2);
            pos60assiB2+=pow(distanza(pos_60[i]),2);
            pos100assiB2+=pow(distanza(pos_100[i]),2);
            pos200assiB2+=pow(distanza(pos_200[i]),2);
            pos300assiB2+=pow(distanza(pos_300[i]),2);
            pos400assiB2+=pow(distanza(pos_400[i]),2);
            
        }
        
        pos10assiB/=B;
        pos20assiB/=B;
        pos30assiB/=B;
        pos40assiB/=B;
        pos50assiB/=B;
        pos60assiB/=B;
        pos100assiB/=B;
        pos200assiB/=B;
        pos300assiB/=B;
        pos400assiB/=B;
        
        pos10assiB2/=B;
        pos20assiB2/=B;
        pos30assiB2/=B;
        pos40assiB2/=B;
        pos50assiB2/=B;
        pos60assiB2/=B;
        pos100assiB2/=B;
        pos200assiB2/=B;
        pos300assiB2/=B;
        pos400assiB2/=B;
        
        pos10aB.push_back(pos10assiB);
        pos20aB.push_back(pos20assiB);
        pos30aB.push_back(pos30assiB);
        pos40aB.push_back(pos40assiB);
        pos50aB.push_back(pos50assiB);
        pos60aB.push_back(pos60assiB);
        pos100aB.push_back(pos100assiB);
        pos200aB.push_back(pos200assiB);
        pos300aB.push_back(pos300assiB);
        pos400aB.push_back(pos400assiB);
        
        pos10aB2.push_back(pow(pos10assiB2,0.5));
        pos20aB2.push_back(pow(pos20assiB2,0.5));
        pos30aB2.push_back(pow(pos30assiB2,0.5));
        pos40aB2.push_back(pow(pos40assiB2,0.5));
        pos50aB2.push_back(pow(pos50assiB2,0.5));
        pos60aB2.push_back(pow(pos60assiB2,0.5));
        pos100aB2.push_back(pow(pos100assiB2,0.5));
        pos200aB2.push_back(pow(pos200assiB2,0.5));
        pos300aB2.push_back(pow(pos300assiB2,0.5));
        pos400aB2.push_back(pow(pos400assiB2,0.5));
        
    }
    
    
    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    //direzione random
    
    vector <double> pos10sB;  //distanze medie su un blocco di 100 particelle
    vector <double> pos20sB;
    vector <double> pos30sB;
    vector <double> pos40sB;
    vector <double> pos50sB;
    vector <double> pos60sB;
    
    vector <double> pos10sB2; //distanze al quadrato medie sotto radice
    vector <double> pos20sB2; //su un blocco di 100 particelle
    vector <double> pos30sB2;
    vector <double> pos40sB2;
    vector <double> pos50sB2;
    vector <double> pos60sB2;
    
    for(int i=0; i<L;i++){
        
        double pos10sferB=0;
        double pos20sferB=0;
        double pos30sferB=0;
        double pos40sferB=0;
        double pos50sferB=0;
        double pos60sferB=0;
        
        double pos10sferB2=0;
        double pos20sferB2=0;
        double pos30sferB2=0;
        double pos40sferB2=0;
        double pos50sferB2=0;
        double pos60sferB2=0;
        
        for(int j=0;j<B;j++){
            
            pos10sferB+=distanza(pos_10sfer[i]);
            pos20sferB+=distanza(pos_20sfer[i]);
            pos30sferB+=distanza(pos_30sfer[i]);
            pos40sferB+=distanza(pos_40sfer[i]);
            pos50sferB+=distanza(pos_50sfer[i]);
            pos60sferB+=distanza(pos_60sfer[i]);
            
            pos10sferB2+=pow(distanza(pos_10sfer[i]),2);
            pos20sferB2+=pow(distanza(pos_20sfer[i]),2);
            pos30sferB2+=pow(distanza(pos_30sfer[i]),2);
            pos40sferB2+=pow(distanza(pos_40sfer[i]),2);
            pos50sferB2+=pow(distanza(pos_50sfer[i]),2);
            pos60sferB2+=pow(distanza(pos_60sfer[i]),2);
            
        }
        
        pos10sferB/=B;
        pos20sferB/=B;
        pos30sferB/=B;
        pos40sferB/=B;
        pos50sferB/=B;
        pos60sferB/=B;
        
        pos10sferB2/=B;
        pos20sferB2/=B;
        pos30sferB2/=B;
        pos40sferB2/=B;
        pos50sferB2/=B;
        pos60sferB2/=B;
        
        pos10sB.push_back(pos10sferB);
        pos20sB.push_back(pos20sferB);
        pos30sB.push_back(pos30sferB);
        pos40sB.push_back(pos40sferB);
        pos50sB.push_back(pos50sferB);
        pos60sB.push_back(pos60sferB);
        
        pos10sB2.push_back(pow(pos10sferB2,0.5));
        pos20sB2.push_back(pow(pos20sferB2,0.5));
        pos30sB2.push_back(pow(pos30sferB2,0.5));
        pos40sB2.push_back(pow(pos40sferB2,0.5));
        pos50sB2.push_back(pow(pos50sferB2,0.5));
        pos60sB2.push_back(pow(pos60sferB2,0.5));
        
    }
    
    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    //utilizzando questi vettori di medie su blocchi posso calcolare la deviazione standard
    
    
    cout<<"Deviazioen standard:"<<endl;
    cout<<"calcolata con la distanza media su blocchi da 100 particelle"<<endl<<endl;
    
    
    double dev10a2=0; //coordinate lungo gli assi
    double dev20a2=0;
    double dev30a2=0;
    double dev40a2=0;
    double dev50a2=0;
    double dev60a2=0;
    double dev100a2=0;
    double dev200a2=0;
    double dev300a2=0;
    double dev400a2=0;
    
    double dev10s2=0; //cordinate sulla sfera
    double dev20s2=0;
    double dev30s2=0;
    double dev40s2=0;
    double dev50s2=0;
    double dev60s2=0;
    
    
    
    dev10a2=st_dev(pos10aB2, dist2_media10,L);
    dev20a2=st_dev(pos20aB2, dist2_media20,L);
    dev30a2=st_dev(pos30aB2, dist2_media30,L);
    dev40a2=st_dev(pos40aB2, dist2_media40,L);
    dev50a2=st_dev(pos50aB2, dist2_media50,L);
    dev60a2=st_dev(pos60aB2, dist2_media60,L);
    dev100a2=st_dev(pos100aB2, dist2_media100,L);
    dev200a2=st_dev(pos200aB2, dist2_media200,L);
    dev300a2=st_dev(pos300aB2, dist2_media300,L);
    dev400a2=st_dev(pos400aB2, dist2_media400,L);
    
    dev10s2=st_dev(pos10sB2, dist2_media10s,L);
    dev20s2=st_dev(pos20sB2, dist2_media20s,L);
    dev30s2=st_dev(pos30sB2, dist2_media30s,L);
    dev40s2=st_dev(pos40sB2, dist2_media40s,L);
    dev50s2=st_dev(pos50sB2, dist2_media50s,L);
    dev60s2=st_dev(pos60sB2, dist2_media60s,L);
    
    
    
    vector <double> err2a_media={dev10a2,dev20a2,dev30a2,dev40a2,dev50a2,dev60a2,dev100a2,dev200a2,dev300a2,dev400a2};
    vector <double> err2s_media={dev10s2,dev20s2,dev30s2,dev40s2,dev50s2,dev60s2};
    
    for(int i=0; i<10;i++){
        outfile_err<<err2a_media[i]<<endl;
    }
    for(int i=0;i<6;i++){
        outfile_errsfer<<err2s_media[i]<<endl;
    }
    
    
    cout<<"passi a=1 lungo un asse"<<endl;
    
    cout<<"10 passi: "<<dev10a2<<endl;
    cout<<"20 passi: "<<dev20a2<<endl;
    cout<<"30 passi: "<<dev30a2<<endl;
    cout<<"40 passi: "<<dev40a2<<endl;
    cout<<"50 passi: "<<dev50a2<<endl;
    cout<<"60 passi: "<<dev60a2<<endl;
    cout<<"100 passi: "<<dev100a2<<endl;
    cout<<"200 passi: "<<dev200a2<<endl;
    cout<<"300 passi: "<<dev300a2<<endl;
    cout<<"400 passi: "<<dev400a2<<endl;
    
    
    cout<<endl<<"passi a=1 lungo una direzione random"<<endl;
    
    cout<<"10 passi: "<<dev10s2<<endl;
    cout<<"20 passi: "<<dev20s2<<endl;
    cout<<"30 passi: "<<dev30s2<<endl;
    cout<<"40 passi: "<<dev40s2<<endl;
    cout<<"50 passi: "<<dev50s2<<endl;
    cout<<"60 passi: "<<dev60s2<<endl;
    cout<<endl;
    
    

    outfile_RW.close();
    outfile_err.close();
    outfile_RWsfer.close();
    outfile_errsfer.close();

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

double distanza(vector <double> part)
{
    
    return pow( pow(part[0],2)+ pow(part[1],2)+pow(part[2],2)  ,0.5);

}


double st_dev(vector<double> dat, double av,int N)
{
    double st_dev=0;
    
    for(int i=0;i<N;i++){
        
        st_dev+=pow(dat[i]-av,2);
        
    }
    
    st_dev/=N;
    
    return pow(st_dev,0.5);
    
}

