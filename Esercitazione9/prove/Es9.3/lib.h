#ifndef __Lib__
#define __Lib__

#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include <vector>
#include <algorithm>
#include <string>



using namespace std;

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


class City
{
private:
    int M_N; //nome della città
    vector <double> M_posPol; //polari: (raggio, angolo)
    vector <double> M_posCart; //cartesiane: (x,y)
public:
    City(int n);
    City(int n, double raggio, double angolo); //costruttore polare
    City(int n, vector<double> p_cart); //costruttore cartesiano
    ~City();
    int GetN(){return M_N;}
    void SetPosPol(double raggio, double angolo);
    void SetPosCart(vector<double> p_cart);
    vector <double> GetPosPol(){return M_posPol;}
    vector <double> GetPosCart(){return M_posCart;}
    double Dist(City c);
    
};


class Lista
{
private:
    int M_NN;
    vector<City> M_VC;
public:
    Lista();
    Lista(int N);
    Lista(vector<City> vc);
    Lista(Lista *l);
    void inizializza(int n, ifstream &a);
    int GetN();
    vector<City> GetList();
    vector<int> Tag();
    void ShowTag();
    void Shuffle();
    vector<double> PosPol_i(int i);
    vector<double> PosCart_i(int i);
    vector<vector<double>> PosPol_N();
    vector<vector<double>> PosCart_N();
    double Dist_N();
    void Swap();
};


class GA : public Lista
{
private:
    int M_NL;
    vector<Lista> M_VL;
    vector<double> M_Dist;
public:
    GA();
    ~GA();
    GA(vector<Lista> &vc);
    int Best();
    int Worst();
    double Average();
    Lista GetList_i(int i);
    Lista GetList_ii(int &i);
    vector<Lista> GetList_2();
    void SetList_i(int i, Lista &l);
    static bool sort_Cost(Lista &c1, Lista &c2);
    void SortByCost();
    vector<Lista> Sorted(); //restituisce le liste in ordine di costo crescente
    vector<Lista> BestPart(double b); //le liste al di sotto della media - b
    vector<Lista> Best_n(int n); //le miglior n liste
    void Show();
    
};


void scrivi(int n, string a); //restituisce coordinate cartesiane
//in [0;5)x[0;5)

Lista Figlia(Lista l1, Lista l2, double r);

void out_file(Lista a, ofstream &output);


#endif // __Lib__


