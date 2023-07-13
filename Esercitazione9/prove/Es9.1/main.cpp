#include "lib.h"

using namespace std;

int main(){
    
    Random rnd;
    int n=9;
    int B=pow(n,3);
    int S=10;
    
    vector <ListaCittà> regione;
    vector <double> sum_dist;
    
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
    
    int d=0;
    
    for(int j=0; j<B;j++){
        ListaCittà C(n);
        regione.push_back(C);
        regione[j].SetCittà(0,1);
        for(int i=1;i<n;i++){
            d=rnd.Rannyu(1,n+1);
            regione[j].SetCittà(i,d);
            
            if(regione[j].Check(i)==1){
                do{
                    d=rnd.Rannyu(1,n+1);
                    regione[j].SetCittà(i,d);
                }while(regione[j].Check(i)==1);
            }
        }
    }
    
    // setto posizioni e le estraggo
    cout<<endl;
    
    for(int j=0; j<B;j++){
        
        for(int i=0; i<n;i++){
            
           // double r=rnd.Rannyu();
            
            //regione[j].SetPosizione_R(i,r);
            regione[j].SetPosizione(i);
            
        }
        cout<<"Costo percorso "<<j<<": "<<regione[j].Costo()<<endl<<endl;
        
    }
    
    cout<<endl;
  
    vector<GA> simulazione;
    GA gen(regione);
    simulazione.push_back(gen);
    cout<<"configurazione originale"<<endl<<endl;
    for(int i=0; i<simulazione.size(); i++){
        simulazione[i].Show();
        cout<<endl;
    }
  
    for(int i=1; i<1;i++){
        
        ListaCittà città2;
        vector<ListaCittà> reg2;
        reg2=regione;
        città2=simulazione[0].GetList_i(1);
        città2.Swap(2,3);
        reg2[1]=città2;
        GA gen2(reg2);
        simulazione.push_back(gen2);
        
        cout<<"generazione "<<i<<endl;
        simulazione[i].Show();
        cout<<endl;
    }
    
    
    cout<<endl<<"////////////////////simulazione///////////////////////"<<endl<<endl;
    
    for(int g=1;g<S;g++){
        
        cout<<"generazione: "<<g<<endl;
        
        vector<ListaCittà> BH;
        
        BH=simulazione[g-1].BestPart();
        
        vector<ListaCittà> L2;
        L2=BH;
        
        
        //cout<<"BH.size()= "<<BH.size()<<", B-BH.size()= "<<B-BH.size()<<endl;
        
        const int bs=BH.size()/3;
        cout<<"bs: "<<bs<<endl;
        
        if(bs==0){
            cout<<"!! convergenza !!"<<endl;
            for(int i=0;i<g;i++){
                
                cout<<"generazione "<<i<<endl;
                simulazione[i].Show();
                cout<<endl;
            }
            for(int i=0;i<g;i++){
                cout<<"sim["<<i<<"] media costo: "<<simulazione[i].Average()<<endl;
            }
            cout<<"vorrei costo: "<<2*2*M_PI<<endl;
            
            cout<<endl<<"!! convergenza !!"<<endl<<endl;
            simulazione[0].GetList_i(0).PrintPosCart("pos_iniz.txt");
            simulazione[g-1].GetList_i(pow(n,2)-10).PrintPosCart("pos_fin.txt");
                
            
            return 0;
            
        }
        
        for(int i=bs; i<B;i++){
            
            int r1=rnd.Rannyu(1,bs);
            int r2=rnd.Rannyu(1,bs);
            double r=rnd.Rannyu();
            if(r1==r2){
                double rr=rnd.Rannyu();
                int rr2=rnd.Rannyu(0,bs);
                if(rr>0.5){
                    
                    L2.push_back(BH[r1].Figlia(BH[rr2],r));
                }else{
                    L2.push_back(BH[rr2].Figlia(BH[r2],r));
                }
                
            }else{
                L2.push_back(BH[r1].Figlia(BH[r2],r));
            }
        }
        
        
      
        GA gen2(L2);
        simulazione.push_back(gen2);
       
        
       // cout<<endl<<endl;
    }
     
  /*
    for(int i=0;i<S;i++){
        
        cout<<"generazione "<<i<<endl;
        simulazione[i].Show();
    }
    */
    for(int i=0;i<S;i++){
        cout<<"sim["<<i<<"] media costo: "<<simulazione[i].Average()<<endl;
    }
    cout<<"vorrei costo: "<<2*2*M_PI<<endl;
   
    
    
    
    
    
        
    simulazione[0].GetList_i(0).PrintPosCart("pos_iniz.txt");
    simulazione[9].GetList_i(pow(n,3)-1).PrintPosCart("pos_fin.txt");
        
    
    
    
    
    
    return 0;
    
}
