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
    
    int N=34; //quante città
    int B=pow(N,3); //grandezza del bacino genetico
    int S=30; //numero delle generazioni
    
    scrivi(N, "dati.txt");
    
    Lista l1;
    
    ifstream file;
    file.open("dati.txt");
    
    l1.inizializza(N ,file);
    
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
    
    
    double i=1; //fine tuning dei migliori percorsi da tenere.
    //media delle distanze - i = elite;
    
    for(int k=1; k<S;k++){
        
        cout<<"Simulazione "<<k<<endl;
        
        
        if(k>3){
            double mm=1;
            
            if(rnd.Rannyu()<mm){
                for(int i=0; i<N;i++){
                    int rl=static_cast<int>(rnd.Rannyu(1, N));
                    simulazione[k-1].GetList_ii(rl).Shuffle();
                }
            }
        }
        
        if(k>5){
            double m=1; //probabilità mutazione (swap di due città)
            
            if(rnd.Rannyu()<m){
                cout<<"$$"<<endl;
                simulazione[k-1].GetList_i(simulazione[k-1].Best()).ShowTag();
                cout<<endl;
                simulazione[k-1].GetList_i(simulazione[k-1].Best()).Swap();
                simulazione[k-1].GetList_i(simulazione[k-1].Best()).ShowTag();
                cout<<endl;
                cout<<"$$"<<endl;
                simulazione[k-1].GetList_i(simulazione[k-1].Best()).Swap();
            }
        }
        
       
        bool attivazione_shuffle = false;
        vector<Lista> bn;
        Lista vec;
        vector<Lista> vec_list2;
                
        if (simulazione[k - 1].GetList_i(simulazione[k - 1].Best()).Dist_N() == simulazione[k - 1].GetList_i(simulazione[k - 1].Worst()).Dist_N()) {
            attivazione_shuffle = true;
            cout << "attivazione shuffle" << endl;
            
            for (int i = 0; i < B / 2; i++) {
                int rl = static_cast<int>(rnd.Rannyu(1, N));
                vec = simulazione[k - 1].GetList_i(rl);
                vec.Shuffle();
                simulazione[k - 1].SetList_i(rl, vec);
            }
            
            if (k % 2 == 0) {
                cout << "metodo 1" << endl;
                bn = simulazione[k - 1].Best_n(2 * N);
            }
            else {
                cout << "metodo 2" << endl;
                bn = simulazione[k - 1].BestPart(i);
                i -= k * 0.2;
            }
            cout<<"&&&&&"<<endl;
            //vector<Lista> vec_list2;
            for(int i=0; i<B; i++){
                Lista f0;
                f0=Figlia(bn[rnd.Rannyu(0,bn.size())],bn[rnd.Rannyu(0,bn.size())],rnd.Rannyu());
                vec_list2.push_back(f0);
            }
            cout<<"&&&&&"<<endl;
            GA genk(vec_list2);
            simulazione.push_back(genk);
            cout<<"&&&&&"<<endl;
        }
        
        
        if (attivazione_shuffle == false) {
            if (k % 2 == 0) {
                cout << "metodo 1" << endl;
                bn = simulazione[k - 1].Best_n(2 * N);
            }
            else {
                cout << "metodo 2" << endl;
                bn = simulazione[k - 1].BestPart(i);
                i -= k * 0.2;
            }
        
        
            vector<Lista> vec_list2;
            for(int i=0; i<B; i++){
                Lista f0;
                f0=Figlia(bn[rnd.Rannyu(0,bn.size())],bn[rnd.Rannyu(0,bn.size())],rnd.Rannyu());
                vec_list2.push_back(f0);
            }
        
            GA genk(vec_list2);
            simulazione.push_back(genk);
        }
        
        cout<<"migliore gen "<<k<<endl;
        simulazione[k].GetList_i(simulazione[k].Best()).ShowTag();
        cout<<endl<<simulazione[k].GetList_i(simulazione[k].Best()).Dist_N()<<endl;
        cout<<"peggiore gen "<<k<<endl;
        simulazione[k].GetList_i(simulazione[k].Worst()).ShowTag();
        cout<<endl<<simulazione[k].GetList_i(simulazione[k].Worst()).Dist_N()<<endl<<endl;
        
        string filename = "dati_" + to_string(k) + ".txt";

        ofstream outputFile;
        outputFile.open(filename);
        
        out_file(simulazione[k].GetList_i(simulazione[k].Best()),outputFile );
        outputFile.close();
    }
    
    
    return 0;
}

