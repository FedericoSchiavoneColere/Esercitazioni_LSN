#ifndef __Lib__
#define __Lib__

#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include <vector>
#include <map>

using namespace std;

// This class contains functions for generating random numbers using the RANNYU algorithm
class Random {

private:
        int m1,m2,m3,m4,l1,l2,l3,l4,n1,n2,n3,n4;

protected:

public:
    // Default constructor
    Random();
    // Destructor
    ~Random();
    // Method to set the seed for the RNG
    void SetRandom(int * , int, int);
    // Method to save the seed to a file
    void SaveSeed();
    // Method to generate a random number in the range [0,1)
    double Rannyu(void);
    // Method to generate a random number in the range [min,max)
    double Rannyu(double min, double max);
    // Method to generate a random number with a Gaussian distribution
    double Gauss(double mean, double sigma);
};


class ListaCittà {
    
private:
    
protected:
    int M_N;
    vector <int> cromosoma;
    vector<vector<double>> vec_pos;
    map<int, vector<double>> M_mappa;
    
public:
    // Default constructor
    ListaCittà();
    //Costruttore
    ListaCittà(int n);
    //costruttore
    ListaCittà(vector <int> v,vector<vector<double>>vp);
    //costruttore con map
    ListaCittà(map<int, vector<double>> mappa);
    // Destructor
    ~ListaCittà();
    //Set N
    void SetN(int n);
    //get N
    int GetN();
    //Set Città
    void SetCittà(int pos, int val);
    //controlla
    int Check(int i);
    //restituisce copia del cromosoma
    vector<int> GetList();
    //get
    int GetCittà(int pos);
    // setta la posizione delle città angolo random
    void SetPosizione_R(int j, double r);
    
    void SetPosMap(int i, map<int, vector<double>> & pos);
    //posizioni regolari
    void SetPosizione(int j);
    // restituisce le posizioni in coordinate polari
    vector <double> GetPosizione(int j);
    //restituisce le posizioni in coordinate cartesiane
    vector <double> GetPosizioneCart(int j);
    //restituisce mappa
    map<int, vector<double>> GetMap();
    //calcola la distanza tra la citta i e i-1
    double Dist_i(int i);
    //calcola la distanza tra due città definite dall'utente
    //(utile per il ritonro all'origine finale)
    double Dist_n(int i,int n);
    //costo, distanza totale
    double Costo();
    ListaCittà Figlia(ListaCittà c1, double r);
    ListaCittà Figlia2(ListaCittà c1,double r);
    
    //swappa due città
    void Swap(int i1,int i2);
    //vector<int> Swap(int i1,int i2);
    //esporta le componenti dei vettori
    void PrintPosCart(string nome);
};


class GA : public ListaCittà{
    
private:
    
protected:
    int M_NN;
    vector<ListaCittà> M_vc;
    vector<double> M_Dist;
public:
    GA();
    ~GA();
    GA(vector<ListaCittà> &vc);
    int Best();
    map<int, vector<double>> Best(int q);
    double Average();
    ListaCittà GetList_i(int i);
    vector<ListaCittà> GetList_2();
    vector<ListaCittà> Swap_2(int i1,int i2);
    static bool sort_Cost(ListaCittà& c1, ListaCittà& c2);
    void SortByCost();
    vector<ListaCittà> BestPart();
    void Show();
};


#endif // __Lib__

