#include "lib.h"

using namespace std;

int main(){
    
    Random rnd;
    int n=5;
    //int B=pow(n,3);
    int B=8;
    int S=1;
    
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
    //&&&&&&&& creo prima le posizioni, poi le carico nelle città sempre uguali
    
    
    map<int,vector<double>> pos_iniz;
  
    for(int i=0; i<n;i++){
        double rr=rnd.Rannyu();
        vector<double> pos;
        double r=2.;
        double theta=2*M_PI*rr;
        
        pos.push_back(r);
        pos.push_back(theta);
        pos_iniz[i]=pos;
    }
    
    
    for (map<int, vector<double>>::const_iterator it = pos_iniz.begin(); it != pos_iniz.end(); ++it) {
        cout<<"** "<<it->first<<" ";
        vector<double> pos;
        pos=it->second;
        cout<<pos[0]<<" "<<pos[1]<<endl;
    }
    
    
 
    
    // ordino le citta random
    
    int d=0;
    
    for(int j=0; j<B;j++){
        ListaCittà C(n);
        regione.push_back(C);
        regione[j].SetCittà(0,1);
        for(int i=1;i<n;i++){
            //d=rnd.Rannyu(1,n+1);
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
            cout<<regione[j].GetCittà(i)<<" ";
        }
        cout<<endl;
        for(int i=0; i<n;i++){
            
            regione[j].SetPosMap(i, pos_iniz);
            cout<<"pos:  "<<regione[j].GetPosizione(i)[0]<<" "<<regione[j].GetPosizione(i)[1]<<endl;
            
        
        }
        cout<<"Costo percorso "<<j<<": "<<regione[j].Costo()<<endl<<endl;
        
    }
    
    cout<<endl;
    
    /*
    map<int, vector<double>>::iterator it = regione[0].GetMap().begin();

    while(it!=regione[0].GetMap().end()){
        cout<<it->first<<" "<<it->second[0]<<" "<<it->second[1]<<endl;
        ++it;
    }
*/
    
    
    

   /*
    for (map<int, vector<double>>::const_iterator it2 = regione[0].GetMap().begin(); it2 != regione[0].GetMap().end(); ++it2) {
        cout<<it2->first<<" "<<it2->second[0]<<" "<<it2->second[1]<<endl;
    }
    
    */

  
    vector<GA> simulazione;
    GA gen(regione);
    simulazione.push_back(gen);
    
    
    
    // Check if GetMap() returns a valid map
    const std::map<int, std::vector<double>>& myMap = regione[0].GetMap();
    if (!myMap.empty()) {
        // Iterate over the map and access its elements
        for (const auto& pair : myMap) {
            int key = pair.first;
            const std::vector<double>& values = pair.second;
            // Access the key and values as needed
            // ...
            
            cout<<key<<" "<<values[0]<<" "<<values[1]<<endl;
            
        }
    } else {
        // Handle the case when the map is empty or nullptr
        std::cout << "Error: GetMap() returned an empty map or nullptr." << std::endl;
    }
    
    

    
    
    cout<<"configurazione originale"<<endl<<endl;
    for(int i=0; i<simulazione.size(); i++){
        cout<<"&&"<<endl;
        simulazione[i].Show();
        cout<<"&&"<<endl;
        cout<<endl;
    }

    
    gen.SortByCost();
    gen.Show();
    
    
    
    
    for(int j=0; j<B;j++){
        
        cout<<"Costo percorso "<<j<<": "<<gen.GetList_i(j).Costo()<<endl<<endl;
        
    }
    /*
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
        
        const int bs=BH.size()/2;
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
            int b=simulazione[g-1].Best();
            simulazione[g-1].GetList_i(1).PrintPosCart("pos_fin.txt");
            cout<<endl<<"costo del best generazione "<<g-1<<" : " <<simulazione[g-1].GetList_i(0).Costo()<<endl;
            
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
     
  
    for(int i=0;i<S;i++){
        
        cout<<"generazione "<<i<<endl;
        simulazione[i].Show();
    }

    for(int i=0;i<S;i++){
        cout<<"sim["<<i<<"] media costo: "<<simulazione[i].Average()<<endl;
    }
    cout<<"vorrei costo: "<<2*2*M_PI<<endl;
   
    
    
    
    
    
        
    simulazione[0].GetList_i(0).PrintPosCart("pos_iniz.txt");
    int b=simulazione[9].Best();
    
    simulazione[9].GetList_i(b).PrintPosCart("pos_fin.txt");
    cout<<endl<<"costo del best generazione 9: "<<simulazione[9].GetList_i(b).Costo()<<endl;

 */
    return 0;
    
}
