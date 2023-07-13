#include "lib.h"



int main(int argc, char** argv){
    
    Random rnd;
    
    int rank, size;
    
    
    int N=50; //quante città
    int B=2000; //grandezza del bacino genetico
    int S=199; //numero delle generazioni
    
    //scrivi(N, "list.txt");
    
    Lista l1;
    
    ifstream file;
    file.open("list.txt");
    
    l1.inizializza(N ,file);
    //l1.Shuffle(); // se si vuole che vengano rimescolati anche i dati presi da file
    vector<Lista> vec_list;
    vec_list.push_back(l1);
    
    for(int i=0; i<B;i++){
        l1.Shuffle();
        vec_list.push_back(l1);
    }


    
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Status stat;
    
    // Generatore di random
    int seed[4];
    int p1, p2;
    ifstream Primes("Primes");
    if (Primes.is_open()){
        Primes >> p1>> p2;
        p1+=rank;
        p2+=rank;
        cout<<p1<<" "<<p2<<endl;
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
    
    for(int i=0; i<200000;i++){
        rnd.Rannyu();
    }
    
    
    cout<<rank<<" "<<size<<endl;
    
    
    vector<GA> simulazione;
    GA gen(vec_list);
    
    simulazione.push_back(gen);
    vector<Lista> Migliori;
    

    
    for (int k=1;k<S;k++) {
        
        cout<<endl<<"rank: "<<rank<<", Simulazione: "<<k<<endl;
        
        vector<Lista> sim_prec= simulazione[k-1].GetList_2();
        GA genk(sim_prec);
        
        simulazione.push_back(genk);
    
        //voglio sostituire il migliore del rank n con il peggiore del rank n+1
        
        MPI_Barrier(MPI_COMM_WORLD);
        int *best_path=new int[N];
        int itag1= 1;
        int itag2= 2;
        int itag3= 3;
        
        if(k%10==0){
        
            if(rank==0){
                for(int i=0; i<N; i++){
                    vector<City> vc = genk.GetList_i(genk.Best()).GetList();
                    best_path[i]=vc[i].GetN();
                }
                MPI_Send(best_path,N,MPI_INTEGER,1,itag1,MPI_COMM_WORLD);
                
            }else if(rank==1){
                MPI_Recv(best_path,N,MPI_INTEGER,0,itag1,MPI_COMM_WORLD, &stat);
                
                vector<int> bp;
                for(int i=0; i<N; i++){
                    cout<<"*"<<best_path[i]<<" ";
                    bp.push_back(best_path[i]);
                }
                cout<<endl;
                genk.SwapBxW(bp); //scambio da un continente all'altro
            }
        }
        
        
        if(k%10==0){
        
            if(rank==1){
                for(int i=0; i<N; i++){
                    vector<City> vc = genk.GetList_i(genk.Best()).GetList();
                    best_path[i]=vc[i].GetN();
                }
                MPI_Send(best_path,N,MPI_INTEGER,2,itag2,MPI_COMM_WORLD);
                
            }else if(rank==2){
                MPI_Recv(best_path,N,MPI_INTEGER,1,itag2,MPI_COMM_WORLD, &stat);
                
                vector<int> bp;
                for(int i=0; i<N; i++){
                    cout<<"*"<<best_path[i]<<" ";
                    bp.push_back(best_path[i]);
                }
                cout<<endl;
                genk.SwapBxW(bp); //scambio da un continente all'altro
            }
        }
        
        if(k%10==0){
        
            if(rank==2){
                for(int i=0; i<N; i++){
                    vector<City> vc = genk.GetList_i(genk.Best()).GetList();
                    best_path[i]=vc[i].GetN();
                }
                MPI_Send(best_path,N,MPI_INTEGER,0,itag3,MPI_COMM_WORLD);
                
            }else if(rank==0){
                MPI_Recv(best_path,N,MPI_INTEGER,2,itag3,MPI_COMM_WORLD, &stat);
                
                vector<int> bp;
                for(int i=0; i<N; i++){
                    cout<<"*"<<best_path[i]<<" ";
                    bp.push_back(best_path[i]);
                }
                cout<<endl;
                genk.SwapBxW(bp); //scambio da un continente all'altro
            }
        }
        
        
        
        
        MPI_Barrier(MPI_COMM_WORLD);
        
        //
        
        if (k > 1) {
            // rimescolo città peggiori
            for (int i = 0; i < 450; i++) {
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
            vector<int> peggiori=genk.Worst_n_index(150);
            for (int i = 0; i < 150; i++) {
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
        
        
        cout<<"rank: "<<rank<<"  "<<genk.GetList_i(genk.Best()).Dist_N()<<"  "<<genk.GetList_i(genk.Average()).Dist_N()<<endl;
        
        if (genk.GetList_i(genk.Best()).Dist_N()>genk.GetList_i(genk.Average()).Dist_N()-5 && k<S-30) {
            attivazione_shuffle = true;
            cout<<"rank: "<<rank<< " attivazione shuffle" << endl;
            vector<int> peggiori=genk.Worst_n_index(500);
            for (int i = 0; i < 500 ; i++) {
                Lista* prova = new Lista(genk.GetList_i(genk.Best()));
                genk.SetList_i(peggiori[i], *prova);
                genk.SetDist_i(peggiori[i]);
                delete prova;
                
            }
            
            cout<<"rank: "<<rank<< " metodo 1" << endl;
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
                cout<<"rank: "<<rank<< " metodo 1" << endl;
                if(jj>60)jj-=0.01;
                else jj=60;
                bn = genk.Best_n(jj);
            }
            else {
                cout<<"rank: "<<rank<< " metodo 2" << endl;
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
        
        cout<<"rank: "<<rank<<" migliore gen "<<k<<endl;
        simulazione[k].GetList_i(simulazione[k].Best()).ShowTag(rank);
        cout<<endl<<"rank: "<<rank<< " "<<", gen "<<k<<" "<<simulazione[k].GetList_i(simulazione[k].Best()).Dist_N()<<endl;
        
        string filename = "dati/dati_" +to_string(rank)+"_"+to_string(k) + ".txt";
        
        ofstream outputFile;
        outputFile.open(filename);
        
        out_file(simulazione[k].GetList_i(simulazione[k].Best()),outputFile );
        outputFile.close();
        
    }
    MPI_Barrier(MPI_COMM_WORLD);
    cout<<endl<<"Fine simulazione rank: "<<rank<<endl<<endl;
    
    GA Mig(Migliori); //per ogni core è il container dei migliori di ogni generazione
    
    //cout<<"**"<<endl;
    //Mig.Show(rank);
    //cout<<"**"<<endl;
    
    //Mig.Sorted();
    
    MPI_Barrier(MPI_COMM_WORLD);
    //GA Mig_dei_Mig; //vorrei raccogliesse i migliori tra i core
    int appo_pos=Mig.Best();
    double appo_dist=Mig.GetList_i(Mig.Best()).Dist_N();
    cout<<"rank: "<<rank<<"  "<<appo_dist<<endl;
    int list_pos[size];
    double list_dist[size];
    
    MPI_Gather(&appo_pos, 1, MPI_INTEGER, list_pos, 1, MPI_INTEGER, 0, MPI_COMM_WORLD);
    MPI_Gather(&appo_dist, 1, MPI_DOUBLE, list_dist, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    
    
    int rank_min_index = 0;
    int itag=0;
    //MPI_Status stat;
    
    double local_min_value = list_dist[0];
    
    
    
    if (rank==0) {
        cout<<endl<<"%$$$% primo rank 0"<<endl<<endl;
        
        for(int i=0; i<size; i++){
            cout<<list_pos[i]<<endl;
        }
        
        for(int i=0; i<size; i++){
            cout<<list_dist[i]<<endl;
        }
        
        //sort(list_dist, list_dist + size);
        
        // Find the minimum value and its position
        //double min_value = list_dist[0];
        
        for (int i = 1; i < size; i++) {
            if (list_dist[i] < local_min_value) {
                local_min_value = list_dist[i];
                //rank_min_index = i;
                rank_min_index=i;
            }
        }
        
        cout<<"rank che contiene il minimo: "<<rank_min_index<<endl;
        
    }
    
    //if(rank==0)MPI_Send(&rank_min_index,1,MPI_INTEGER,,itag,MPI_COMM_WORLD);
    //else if (rank!=0)MPI_Recv(&rank_min_index,1,MPI_INTEGER,0,itag,MPI_COMM_WORLD, &stat);
    MPI_Barrier(MPI_COMM_WORLD);
    if (rank == 0) {
        cout<<endl<<"send "<<rank<<endl;
        for (int dest_rank = 1; dest_rank < size; dest_rank++) {
            MPI_Send(&rank_min_index, 1, MPI_INT, dest_rank, itag, MPI_COMM_WORLD);
        }
    } else {
        MPI_Barrier(MPI_COMM_WORLD);
        cout<<rank<<" rec"<<endl;
        MPI_Recv(&rank_min_index, 1, MPI_INT, 0, itag, MPI_COMM_WORLD, &stat);
    }
    MPI_Barrier(MPI_COMM_WORLD);
    
    cout<<"rank: "<<rank<<" itag: "<<itag<<" rank_min_index: "<<rank_min_index<<endl;
    
    
    MPI_Barrier(MPI_COMM_WORLD);
    
    if(rank==0){
        cout<<endl<<"%$$$% "<<rank<<endl<<endl;
        Mig.GetList_i(Mig.Best()).ShowTag(rank);
        cout<<endl<<"Distanza: "<<appo_dist<<endl;
    }
    MPI_Barrier(MPI_COMM_WORLD);
    if(rank==1){
        cout<<endl<<"%$$$% "<<rank<<endl<<endl;
        Mig.GetList_i(Mig.Best()).ShowTag(rank);
        cout<<endl<<"Distanza: "<<appo_dist<<endl;
    }
    MPI_Barrier(MPI_COMM_WORLD);
    if(rank==2){
        cout<<endl<<"%$$$% "<<rank<<endl<<endl;
        Mig.GetList_i(Mig.Best()).ShowTag(rank);
        cout<<endl<<"Distanza: "<<appo_dist<<endl;
    }
    
    MPI_Barrier(MPI_COMM_WORLD);
    
    
     if(rank==rank_min_index){
        cout<<endl<<"%%%%%%"<<endl;
        cout<<"Rank percorso migliore: "<<rank<<endl;
        Mig.GetList_i(Mig.Best()).ShowTag(rank_min_index);
        cout<<endl<<"Distanza: "<<appo_dist<<endl;
         
         string filename = "dati/dati_Best.txt";
         
         ofstream outputFile;
         outputFile.open(filename);
         
         out_file(Mig.GetList_i(Mig.Best()),outputFile);
         outputFile.close();
        
         cout<<"%%%%%%"<<endl;
    }
    
    
    
    MPI_Barrier(MPI_COMM_WORLD);
    if(rank==0){
        
        cout<<"Fine del programma"<<endl;
    }
    
 
    int MPI_Abort(MPI_Comm comm, int errorcode);
   
    //MPI_Finalize();
    
    return 0;
}
