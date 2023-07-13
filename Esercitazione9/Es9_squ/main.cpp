#include "lib.h"



int main(){
    
    Random rnd;
    
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
    
    for(int i=0; i<100000;i++){
        rnd.Rannyu();
    }
    
    
    
    int N=34; //quante città
    int B=pow(N,2); //grandezza del bacino genetico
    int S=250; //numero delle generazioni
    
    scrivi(N, "dati/dati.txt");
    
    Lista l1;
    
    ifstream file;
    file.open("dati/dati.txt");
    
    l1.inizializza(N ,file);
    //l1.Shuffle(); // se si vuole che vengano rimescolati anche i dati presi da file
    vector<Lista> vec_list;
    vec_list.push_back(l1);
    
    for(int i=0; i<B;i++){
        l1.Shuffle();
        vec_list.push_back(l1);
    }
    
    cout<<endl;
    vector<GA> simulazione;
    GA gen(vec_list);
    
    simulazione.push_back(gen);
    vector<Lista> Migliori;
    
    
    for(int k=1; k<S;k++){
        
        cout<<"Simulazione "<<k<<endl;
        
        vector<Lista> sim_prec= simulazione[k-1].GetList_2();
        GA genk(sim_prec);

        simulazione.push_back(genk);
        
        
        
        
        if (k > 1) {
            // rimescolo città peggiori
            for (int i = 0; i < 400; i++) {
                int ll=genk.Worst();
                Lista* prova = new Lista(genk.GetList_i(ll));
                prova->Shuffle();
                genk.SetList_i(ll, *prova);
                genk.SetDist_i(ll);
                delete prova;
            }
        }
        

        if(k>1){
            //probabilità mutazione (swap di due città)
            for (int i = 0; i < 400; i++) {
                int rl = static_cast<int>(rnd.Rannyu(1, B));
                Lista* prova = new Lista(genk.GetList_i(genk.Best()));
                prova->Swap(rnd.Rannyu(1,N),rnd.Rannyu(1,N));
                genk.SetList_i(rl, *prova);
                genk.SetDist_i(rl);
                delete prova;
            }
        }
        
        if(k>1){
            //duplico la città migliore e le sostituisco alle peggiori
            vector<int> peggiori=genk.Worst_n_index(100);
            for (int i = 0; i < 100; i++) {
                Lista* prova = new Lista(genk.GetList_i(genk.Best()));
                genk.SetList_i(peggiori[i], *prova);
                genk.SetDist_i(peggiori[i]);
                delete prova;
            }
        }
        
        
        double jj=100;
        double j=40;
        
        bool attivazione_shuffle = false;
        vector<Lista> bn;
        vector<Lista> vec_list2;
        
        
        cout<<genk.GetList_i(genk.Best()).Dist_N()<<"  "<<genk.GetList_i(genk.Average()).Dist_N()<<endl;
        
        if (genk.GetList_i(genk.Best()).Dist_N()>genk.GetList_i(genk.Average()).Dist_N()-7 && k<S-30) {
            attivazione_shuffle = true;
            cout << "attivazione shuffle" << endl;
            vector<int> peggiori=genk.Worst_n_index(500);
            for (int i = 0; i < 500 ; i++) {
                Lista* prova = new Lista(genk.GetList_i(genk.Best()));
                genk.SetList_i(peggiori[i], *prova);
                genk.SetDist_i(peggiori[i]);
                delete prova;
                
            }
            
            cout << "metodo 1" << endl;
            bn = genk.Best_n(10);
            
            
            for(int i=0; i<B; i++){
                
                int rl1 = static_cast<int>(rnd.Rannyu(0, bn.size()));
                int rl2 = static_cast<int>(rnd.Rannyu(0, bn.size()));
                
                Lista prova= Figlia(genk.GetList_i(rl1),genk.GetList_i(rl2),rnd.Rannyu(),rnd.Rannyu());
                vec_list2.push_back(prova);
            }
            simulazione[k] = vec_list2;
        }
        
        
        if (attivazione_shuffle == false) {
            if (k % 2 == 0) {
                cout << "metodo 1" << endl;
                if(jj>60)jj-=0.01;
                else jj=60;
                bn = genk.Best_n(jj);
            }
            else {
                cout << "metodo 2" << endl;
                if(j>10)j-=0.01;
                else j=10;
                bn = genk.Best_n(j);
                
            }
            
            for(int i=0; i<B; i++){
                
                int rl1 = static_cast<int>(rnd.Rannyu(0, bn.size()));
                int rl2 = static_cast<int>(rnd.Rannyu(0, bn.size()));
                
                Lista prova= Figlia(genk.GetList_i(rl1),genk.GetList_i(rl2),rnd.Rannyu(),rnd.Rannyu());
                vec_list2.push_back(prova);
            }
            simulazione[k] = vec_list2;
            
        }
            
        Migliori.push_back(simulazione[k].GetList_i(simulazione[k].Best()));
        
        cout<<"migliore gen "<<k<<endl;
        simulazione[k].GetList_i(simulazione[k].Best()).ShowTag();
        cout<<endl<<simulazione[k].GetList_i(simulazione[k].Best()).Dist_N()<<endl;
        
        string filename = "dati/dati_" + to_string(k) + ".txt";

        ofstream outputFile;
        outputFile.open(filename);
        
        out_file(simulazione[k].GetList_i(simulazione[k].Best()),outputFile );
        outputFile.close();
        
    }
    
    
    GA Mig(Migliori);
    int i=0;
    
    i=Mig.Best();
    
    cout<<"il migliore percorso: ";
    Mig.GetList_i(i).ShowTag();
    cout<<endl<<"Distanza: "<<Mig.GetList_i(i).Dist_N();
    cout<<endl;
    
    
    ofstream outputFile;
    outputFile.open("dati/Dati_mig.txt");
    
    out_file(Mig.GetList_i(i),outputFile );
    outputFile.close();
    
    
    return 0;
}

