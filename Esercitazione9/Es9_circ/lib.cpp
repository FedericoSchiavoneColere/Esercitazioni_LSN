#include "lib.h"

using namespace std;

Random :: Random(){}
// Default constructor, does not perform any action

Random :: ~Random(){}
// Default destructor, does not perform any action

void Random :: SaveSeed(){
   // This function saves the current state of the random number generator to a file "seed.out"
   ofstream WriteSeed;
   WriteSeed.open("seed.out");
   if (WriteSeed.is_open()){
      WriteSeed << l1 << " " << l2 << " " << l3 << " " << l4 << endl;;
   } else cerr << "PROBLEM: Unable to open random.out" << endl;
  WriteSeed.close();
  return;
}

double Random :: Gauss(double mean, double sigma) {
   // This function generates a random number from a Gaussian distribution with given mean and sigma
   double s=Rannyu();
   double t=Rannyu();
   double x=sqrt(-2.*log(1.-s))*cos(2.*M_PI*t);
   return mean + x * sigma;
}

double Random :: Rannyu(double min, double max){
   // This function generates a random number in the range [min, max]
   return min+(max-min)*Rannyu();
}

double Random :: Rannyu(void){
  // This function generates a random number in the range [0,1)
  const double twom12=0.000244140625;
  int i1,i2,i3,i4;
  double r;

  i1 = l1*m4 + l2*m3 + l3*m2 + l4*m1 + n1;
  i2 = l2*m4 + l3*m3 + l4*m2 + n2;
  i3 = l3*m4 + l4*m3 + n3;
  i4 = l4*m4 + n4;
  l4 = i4%4096;
  i3 = i3 + i4/4096;
  l3 = i3%4096;
  i2 = i2 + i3/4096;
  l2 = i2%4096;
  l1 = (i1 + i2/4096)%4096;
  r=twom12*(l1+twom12*(l2+twom12*(l3+twom12*(l4))));

  return r;
}

void Random :: SetRandom(int * s, int p1, int p2){
  // This function sets the seed and parameters of the random number generator
  m1 = 502;
  m2 = 1521;
  m3 = 4071;
  m4 = 2107;
  l1 = s[0];
  l2 = s[1];
  l3 = s[2];
  l4 = s[3];
  n1 = 0;
  n2 = 0;
  n3 = p1;
  n4 = p2;

  return;
}


City::City(int n)
{
    M_N=n;
    vector <double> posPol(2,0);
    vector <double> posCart(2,0);
    M_posPol=posPol;
    M_posCart=posCart;}

City::City(int n, double raggio, double angolo)
{
    M_N=n;
    vector <double> posPol(2,0);
    vector <double> posCart(2,0);
    M_posPol=posPol;
    M_posCart=posCart;
    
    M_posPol[0]=raggio;
    M_posPol[1]=angolo;
    
    M_posCart[0]= raggio*cos(angolo);
    M_posCart[1]= raggio*sin(angolo);
}

City::City(int n, vector<double> p_cart)
{
    M_N=n;
    vector <double> posPol(2,0);
    vector <double> posCart(2,0);
    M_posPol=posPol;
    M_posCart=posCart;
    
    M_posCart[0]=p_cart[0];
    M_posCart[1]=p_cart[1];
    
    M_posPol[0]=pow(pow(M_posCart[0],2)+pow(M_posCart[1],2),0.5);
    M_posPol[1]=atan(M_posCart[1]/M_posCart[0]);
}
City::~City()
{
}

void City::SetPosPol(double raggio, double angolo)
{
    vector <double> posPol(2,0);
    vector <double> posCart(2,0);
    M_posPol=posPol;
    M_posCart=posCart;
    
    M_posPol[0]=raggio;
    M_posPol[1]=angolo;
    
    M_posCart[0]= raggio*cos(angolo);
    M_posCart[1]= raggio*sin(angolo);
}

void City::SetPosCart(vector<double> p_cart)
{
    vector <double> posPol(2,0);
    vector <double> posCart(2,0);
    M_posPol=posPol;
    M_posCart=posCart;
    
    M_posCart[0]=p_cart[0];
    M_posCart[1]=p_cart[1];
    
    M_posPol[0]=pow(pow(M_posCart[0],2)+pow(M_posCart[1],2),0.5);
    M_posPol[1]=atan(M_posCart[1]/M_posCart[0]);
}

double City::Dist(City c)
{
    return pow(pow(c.GetPosCart()[0] - M_posCart[0],2)+pow(c.GetPosCart()[1] - M_posCart[1],2),0.5);
}


Lista::Lista()
{
    M_NN=0;
}

Lista::Lista(int n)
{
    M_NN=n;
    vector<City> nn(n,NULL);
    M_VC=nn;
}

Lista::Lista(vector<City> vc)
{
    M_NN=vc.size();
    M_VC=vc;
}

Lista::Lista(Lista *l){
    M_NN=l->GetN();
    M_VC=l->GetList();
}

void Lista::inizializza(int n, ifstream &a)
{
    M_NN=n;
    
    for(int i=0; i<n;i++){
        City c(i+1);
        vector <double> p(2,0);
        for(int j=0; j<2;j++){
            a>>p[j];
        }
        c.SetPosCart(p);
        M_VC.push_back(c);
    }
}

int Lista::GetN()
{
    return M_NN;
}
vector<City> Lista::GetList()
{
    return M_VC;
}

vector<int> Lista::Tag()
{
    vector<int> t;
    
    for(int i=0; i<M_NN;i++){
        t.push_back(M_VC[i].GetN());
    }
    return t;
}

void Lista::ShowTag()
{
    for(int i=0; i<M_NN;i++){
        cout<<M_VC[i].GetN()<<" ";
    }
}

void Lista::Shuffle()
{
    random_shuffle(M_VC.begin()+1, M_VC.end());
}

vector<double> Lista::PosPol_i(int i)
{
    return M_VC[i].GetPosPol();
}

vector<double> Lista::PosCart_i(int i)
{
    return M_VC[i].GetPosCart();
}

vector<vector<double>> Lista::PosPol_N()
{
    vector<vector<double>> lp;
    for(int i=0; i<M_NN;i++){
        lp.push_back(M_VC[i].GetPosPol());
    }
    return lp;
}

vector<vector<double>> Lista::PosCart_N()
{
    vector<vector<double>> lp;
    for(int i=0; i<M_NN;i++){
        lp.push_back(M_VC[i].GetPosCart());
    }
    return lp;
}

double Lista::Dist_N()
{
    double d=0;
    
    for(int i=1; i<M_NN;i++){
        d+=M_VC[i].Dist(M_VC[i-1]);
    }
    
    d+=M_VC[M_NN-1].Dist(M_VC[0]);
    return d;
}

void Lista::Swap(double r1,double r2)
{
    
    
    int index1 = static_cast<int>(r1) ;
    int index2 = static_cast<int>(r2);
    swap(M_VC[index1], M_VC[index2]);
}


GA::GA(){}

GA::~GA(){}

GA::GA(vector<Lista> &vc)
{
    M_VL=vc;
    M_NL=M_VL.size();
    for(int i=0;i<M_NL;i++){
        M_Dist.push_back(M_VL[i].Dist_N());
    }
}

int GA::Best()
{
    int min_element_index = min_element(M_Dist.begin(), M_Dist.end())-M_Dist.begin();
    
    return min_element_index;
}

int GA::Worst()
{
        auto max_element_iterator = max_element(M_Dist.begin(), M_Dist.end());
        int max_element_index = distance(M_Dist.begin(), max_element_iterator);
        return max_element_index;
}

double GA::Average()
{
    double av=0;
    
    for(int i=0;i<M_NL;i++){
        av+=M_Dist[i];
    }
    av/=M_NL;
    return av;
}

vector<Lista> GA::GetList_2()
{
    return M_VL;
}

Lista GA::GetList_i(int i)
{
    return M_VL[i];
}

Lista GA::GetList_ii(int &i)
{
    return &M_VL[i];
}

void GA::SetList_i(int i, Lista &l)
{
    M_VL[i]=l;
}

void GA::SetDist_i(int i)
{
    M_Dist[i]=M_VL[i].Dist_N();
}

bool GA::sort_Cost(Lista &c1, Lista &c2) {
    
    return c1.Dist_N() < c2.Dist_N();
}

void GA::SortByCost()
{
    sort(M_VL.begin(), M_VL.end(), sort_Cost);

    for (Lista &CC : M_VL) {
        cout << CC.Dist_N()<<" ";
    }

    cout<<endl<<endl;
}


vector<Lista> GA::Sorted()
{
    vector<Lista> sorted;
    sorted=M_VL;
    sort(sorted.begin(), sorted.end(), sort_Cost);
    
    return sorted;
}

vector<Lista> GA::BestPart(double b)
{
    double av=Average();
    vector<Lista> BH;
    for(int i=0; i<M_NL;i++){
        if(M_VL[i].Dist_N()<av-b) BH.push_back(M_VL[i]);
    }
    
    return BH;
}

vector<Lista> GA::Best_n(int n)
{
    vector<Lista> mig;
    vector<Lista> mig_n;
    mig=Sorted();
    for(int i=0;i<n;i++){
        mig_n.push_back(mig[i]);
    }
    mig.clear();
    return mig_n;
}


vector<int> GA::Worst_n_index(int n)
{
    vector<double> dist_sort = M_Dist;
    sort(dist_sort.begin(), dist_sort.end());

    vector<int> indices;

    for (int i = 1; i <= n; ++i) {
        double largestNumber = dist_sort[dist_sort.size() - i];  // Get the ith largest number
        auto it = upper_bound(M_Dist.begin(), M_Dist.end(), largestNumber);  // Find the first element greater than largestNumber
        int index = distance(M_Dist.begin(), it) - 1;  // Calculate the index
        indices.push_back(index);
    }

    return indices;
}


void GA::Show()
{
    for(int i=0;i<M_NL;i++){
        M_VL[i].ShowTag();
        cout<<endl<<"Costo percorso "<<i<<": "<<M_VL[i].Dist_N()<<endl<<endl;
    }
}


GA& GA::operator=(const GA& other) {
    if (this != &other) {
    
        M_NL = other.M_NL;
        M_VL = other.M_VL;
        M_Dist = other.M_Dist;
    }
    return *this;
}


void scrivi(int n, string a)
{
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
    ofstream file_out;
    file_out.open(a);
    
    double R=2.5;
    
    for(int i=0;i<n;i++){
        
        double theta=rnd.Rannyu(0,2*M_PI);
        
        file_out<<R*cos(theta)<<" "<<R*sin(theta)<<endl;
        
    }
    file_out.close();
}


Lista Figlia(Lista l1, Lista l2,double r,double rr)
{
    vector<City> v_prova;
    
    int taglia=2+r*(l1.GetN()-2);

    if(taglia>1){
        for(int i=0;i<taglia;i++){
            v_prova.push_back(l1.GetList()[i]);
        }
    }else{
        cout<<endl<<"Errore nella funzione Figlia"<<endl<<endl;
        return 0;
        
    }
    
    for(int i=taglia;i<l1.GetN();i++){
        v_prova.push_back(l2.GetList()[i]);
    }
      
    Lista prova(v_prova);
 
    bool has_duplicates = false;

    for (size_t i = 0; i < l1.GetN() - 1; ++i) {
        for (size_t j = i + 1; j < l1.GetN(); ++j) {
            if (prova.Tag()[i] == prova.Tag()[j]) {
                has_duplicates = true;
                break;
            }
        }
        if (has_duplicates) {
            break;
        }
    }
        
    if (has_duplicates) {
        if(l1.Dist_N()>l2.Dist_N()){
            
            if(rr<0.75){
                l2.Shuffle();
                return l2;
            }else { return l2;
            }
        }else{
            if(rr<0.75){
                l1.Shuffle();
                return l1;
            } else { return l1;
            }
        }
    }else{
        return prova;
    }
}


void out_file(Lista a,ofstream &output)
{
    for(int i=0; i<a.GetN();i++){
        output<<a.GetList()[i].GetPosCart()[0]<<"  "<<a.GetList()[i].GetPosCart()[1]<<endl;
    }
}


