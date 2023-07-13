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



ListaCittà::ListaCittà(){
    M_N=0;
}

ListaCittà::ListaCittà(int n){
    M_N=n;
    for(int i=0; i<M_N;i++){
        cromosoma.push_back(0);
        vector<double> p(2,0);
        vec_pos.push_back(p);
        M_mappa[i]=p;
    }
    
}

ListaCittà::ListaCittà(vector <int> v,vector<vector<double>> vp)
{
    M_N=v.size();
    for(int i=0; i<M_N;i++){
        cromosoma.push_back(v[i]);
        vector<double> p(2,0);
        p[0]=vp[i][0];
        p[1]=vp[i][1];
        vec_pos.push_back(p);
        
        M_mappa[i]=p;
        
        for (map<int, vector<double>>::const_iterator it = M_mappa.begin(); it != M_mappa.end(); ++it) {
            cout << it->first << endl;
        }
    }
}


ListaCittà::ListaCittà(map<int, vector<double>> mappa)
{
    M_mappa=mappa;
    M_N=mappa.size();
    
    map<int, vector<double>>::iterator it = mappa.begin();

    while (it != mappa.end())
      {
        cromosoma.push_back(it->first);
        vec_pos.push_back(it->second);
    }
}

ListaCittà::~ListaCittà(){}

void ListaCittà::SetN(int n)
{
    M_N=n;
}

int ListaCittà::GetN()
{
    return M_N;
}

void ListaCittà::SetCittà(int pos, int val)
{
    cromosoma[pos]=val;
}

int ListaCittà::Check(int i)
{
    int test=0;
    
    for(int j=0;j<M_N;j++){
        if(j!=i){
            if(cromosoma[j]==cromosoma[i]){
                test=1;
                return test;
            }
        }
    }
    return test;
}

vector<int>ListaCittà::GetList()
{
    return cromosoma;
}

int ListaCittà::GetCittà(int pos)
{
    return cromosoma[pos];
}

void ListaCittà:: SetPosizione_R(int j, double rr)
{
    vector<double> pos;
    //int cc=cromosoma[j];
    double r=2.;
    //double theta=2*cc/static_cast<double>(M_N);
    double theta=2*M_PI*rr;
    
    pos.push_back(r);
    pos.push_back(theta);
    vec_pos[j]=pos;
    M_mappa[j]=pos;

}

void ListaCittà::SetPosMap(int i, map<int, vector<double>>& pos)
{
    int k = cromosoma[i] - 1;  // Adjust the index by subtracting 1

    auto it = pos.find(k);
    if (it != pos.end()) {
        cout << "k: " << k + 1 << " pos: " << it->second[0] << " " << it->second[1] << endl;

        vec_pos[i] = it->second;
        M_mappa[i] = it->second;
    }
    else {
        cout << "Error: Position not found for city " << k + 1 << endl;
    }
}


void ListaCittà:: SetPosizione(int j)
{
    vector<double> pos;
    int cc=cromosoma[j];
    double r=2.;
    double theta=2*cc/static_cast<double>(M_N);
        
    pos.push_back(r);
    pos.push_back(theta);
    vec_pos[j]=pos;
    M_mappa[j]=pos;
    
}


vector<double> ListaCittà:: GetPosizione(int j)
{
    vector<double> pos(2,0);
    pos[0]=vec_pos[j][0];
    pos[1]=vec_pos[j][1];
    return pos;
    
}

vector<double> ListaCittà:: GetPosizioneCart(int j)
{
    vector<double> pos_cart;

    double x=vec_pos[j][0]*cos(vec_pos[j][1]*M_PI);
    double y=vec_pos[j][0]*sin(vec_pos[j][1]*M_PI);
    
    pos_cart.push_back(x);
    pos_cart.push_back(y);
    
    return pos_cart;
}
map<int, vector<double>> ListaCittà::GetMap()
{
    return M_mappa;
}
double ListaCittà::Dist_i(int i)
{
    double l=0;
    
    vector<double> pos_cart_i;
    vector<double> pos_cart_i_1;
    pos_cart_i=GetPosizioneCart(i);
    pos_cart_i_1=GetPosizioneCart(i-1);
   
    l=pow(pow(pos_cart_i[0]-pos_cart_i_1[0],2)+pow(pos_cart_i[1]-pos_cart_i_1[1],2), 0.5);
    
    return l;
}

double ListaCittà::Dist_n(int i,int n)
{
    double l=0;
    
    vector<double> pos_cart_i;
    vector<double> pos_cart_i_n;
    pos_cart_i=GetPosizioneCart(i);
    pos_cart_i_n=GetPosizioneCart(n);
    l= pow(pow(pos_cart_i[0]-pos_cart_i_n[0],2)+pow(pos_cart_i[1]-pos_cart_i_n[1],2),0.5);
    
    return l;
}

double ListaCittà::Costo()
{
    double L=0;
    for(int i=0;i<M_N;i++){
        if(i>0){
            L+=Dist_i(i);
        }
        if(i==M_N-1){
        L+=Dist_n(i,0);
        }
    }
    return L;
}

void ListaCittà::Swap(int i1, int i2)
{
    int appo = cromosoma[i1];
    cromosoma[i1] = cromosoma[i2];
    cromosoma[i2] = appo;
}

/*
vector<int> ListaCittà::Swap(int i1, int i2)
{
    vector<int> appo = cromosoma;
    
    swap(appo[i1], appo[i2]);
    
    return appo;
}
*/


ListaCittà ListaCittà::Figlia(ListaCittà c1,double r)
{
    Random rnd;
    //
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
    
    ListaCittà prova(M_N);
    int taglia=2+r*(M_N-2);

    //cout<<"taglio: "<<taglia<<" "<<" su un totoale di "<<M_N<<endl;
    if(taglia>1){
        for(int i=0;i<taglia;i++){
            
            prova.SetCittà(i,cromosoma[i]);
            prova.SetPosizione(i);
        }
    }else{
        cout<<"sbarella"<<endl;
        return 0;
        
    }
    
    //cout<<"vorrei inserire: ";
    for(int i=taglia; i<M_N;i++){
        //cout<<c1.GetCittà(i)<<" ";
    }
    
    //cout<<endl;
    for(int i=taglia;i<M_N;i++){
        //cout<<" prima "<<cromosoma[i]<<endl;
        prova.SetCittà(i,c1.GetCittà(i));
        prova.SetPosizione(i);
        //cout<<"i:"<<i<<" "<<prova.GetCittà(i)<<" dopo "<<c1.GetCittà(i)<<endl;
        if(prova.GetCittà(i)==0)cout<<"qui!!"<<endl;
        
        
    }
      
    bool has_duplicates = false;

    for (size_t i = 0; i < prova.GetList().size() - 1; ++i) {
        for (size_t j = i + 1; j < prova.GetList().size(); ++j) {
            if (prova.GetList()[i] == prova.GetList()[j]) {
                has_duplicates = true;
                break;
            }
        }
        if (has_duplicates) {
            break;
            }
    }

    if (has_duplicates) {
        //std::cout << "There are duplicates in the vector.\n";
            //cout<<endl;
            
            
        if(r>0.5){
                //cout<<"c1"<<endl;
            return c1;
        }else{
                //cout<<"c2"<<endl;
            ListaCittà c2(cromosoma,vec_pos);
            return c2;
        }

    } else {
        //std::cout << "There are no duplicates in the vector.\n";
        //cout<<endl;
        for(int k=0;k<M_N;k++){
            //cout<<prova.GetCittà(k)<<" ";
            vector<double> pp(2,0);
            pp=prova.GetPosizioneCart(k);
            //cout<<" £ "<<pp[0]<<" "<<pp[1]<<endl;
        }
            
            //cout<<"(((()))"<<prova.Costo()<<endl;
        return prova;
    }
}
 
void ListaCittà::PrintPosCart(string nome)
{
    ofstream out_file;
    out_file.open(nome);
    
    
    for(int i=0;i<M_N;i++){
        vector<double> pos(2,0);
        pos=GetPosizioneCart(i);
        out_file<<pos[0]<<" "<<pos[1]<<endl;
    }
    
    out_file.close();
    
    
}


//&&&&&&&&&&&&&&&&&&&&&&&&&&

///*


ListaCittà ListaCittà::Figlia2(ListaCittà c1,double r)
{
    Random rnd;
    //
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
    
    ListaCittà prova(M_N);
    int taglia=2+r*(M_N-2);

    //cout<<"taglio: "<<taglia<<" "<<" su un totoale di "<<M_N<<endl;
    if(taglia>1){
        for(int i=0;i<taglia;i++){
            
            prova.SetCittà(i,cromosoma[i]);
            prova.SetPosizione(i);
        }
    }else{
        cout<<"sbarella"<<endl;
        return 0;
        
    }
    
    //cout<<"vorrei inserire: ";
    for(int i=taglia; i<M_N;i++){
        //cout<<c1.GetCittà(i)<<" ";
    }
    
    //cout<<endl;
    for(int i=taglia;i<M_N;i++){
        //cout<<" prima "<<cromosoma[i]<<endl;
        prova.SetCittà(i,c1.GetCittà(i));
        prova.SetPosizione(i);
        //cout<<"i:"<<i<<" "<<prova.GetCittà(i)<<" dopo "<<c1.GetCittà(i)<<endl;
        if(prova.GetCittà(i)==0)cout<<"qui!!"<<endl;
        
        
    }
      
    bool has_duplicates = false;

    for (size_t i = 0; i < prova.GetList().size() - 1; ++i) {
        for (size_t j = i + 1; j < prova.GetList().size(); ++j) {
            if (prova.GetList()[i] == prova.GetList()[j]) {
                has_duplicates = true;
                break;
            }
        }
        if (has_duplicates) {
            break;
            }
    }

    if (has_duplicates) {
        //std::cout << "There are duplicates in the vector.\n";
            //cout<<endl;
            
            
        if(r>0.5){
                //cout<<"c1"<<endl;
            return c1;
        }else{
                //cout<<"c2"<<endl;
            ListaCittà c2(cromosoma,vec_pos);
            return c2;
        }

    } else {
        //std::cout << "There are no duplicates in the vector.\n";
        //cout<<endl;
        for(int k=0;k<M_N;k++){
            //cout<<prova.GetCittà(k)<<" ";
            vector<double> pp(2,0);
            pp=prova.GetPosizioneCart(k);
            //cout<<" £ "<<pp[0]<<" "<<pp[1]<<endl;
        }
            
            //cout<<"(((()))"<<prova.Costo()<<endl;
        return prova;
    }
}

 //*/

//&&&&&&&&&&&&&&&&&&&&&&&&&&

GA::GA(){}

GA::~GA(){}

GA::GA(vector<ListaCittà> &vc)
{
    M_vc=vc;
    M_NN=M_vc.size();
    for(int i=0;i<M_NN;i++){
        M_Dist.push_back(M_vc[i].Costo());
    }
}
int GA::Best()
{
    int min_element_index = min_element(M_Dist.begin(), M_Dist.end())-M_Dist.begin();
    
    return min_element_index;
}

/*
map<int, vector<double>> GA::Best(int q) //q= quanti migliori voglio selezionare
{
    map<int, vector<double>> migliori;
    map<int, vector<double>> migliori_q;
    vector<double> costi;
    
    // Populate the vector with index and cost information
     for (int i = 0; i < M_vc.size(); i++) {
         costi.push_back(make_pair(i, M_vc[i].Costo()));
     }

     // Sort the vector based on cost in ascending order
     sort(costi.begin(), costi.end(), [](const auto& a, const auto& b) {
         return a.second < b.second;
     });

     // Select the top q best solutions
     for (int i = 0; i < q; i++) {
         int index = costi[i].first;
         migliori_q[index] = M_vc[index].GetPosizione(); // Assuming GetPosizione() returns a vector<double> for the position
     }
    
    return migliori_q;
}

*/


double GA::Average()
{
    double av=0;
    
    for(int i=0;i<M_NN;i++){
        av+=M_Dist[i];
    }
    av/=M_NN;
    return av;
}

vector<ListaCittà> GA::GetList_2()
{
    return M_vc;

}

ListaCittà GA::GetList_i(int i)
{
    return M_vc[i];

}

vector<ListaCittà> GA::Swap_2(int i1,int i2)
{
    /*
    ListaCittà appo = M_vc[i1];
    M_vc[i1] = M_vc[i2];
    M_vc[i2] = appo;
    */
    vector<ListaCittà> appo;
    appo=M_vc;
    
    swap(appo[i1], appo[i2]);
    
    return appo;
  
}



bool GA::sort_Cost(ListaCittà& c1, ListaCittà& c2) {
    
    return c1.Costo() < c2.Costo();
}

void GA::SortByCost()
{
    sort(M_vc.begin(), M_vc.end(), sort_Cost);

    for (ListaCittà& CC : M_vc) {
        cout << CC.Costo() << " ";
    }

    cout<<endl<<endl;
}

    
vector<ListaCittà> GA::BestPart()
{
    double av=Average();
    vector<ListaCittà> BH;
    for(int i=0; i<M_NN;i++){
        if(M_vc[i].Costo()<av-0.5) BH.push_back(M_vc[i]);
    }
    
    return BH;
}

void GA::Show()
{
    for(int i=0;i<M_NN;i++){
    
        const std::map<int, std::vector<double>>& myMap = M_vc[i].GetMap();
        if (!myMap.empty()) {
            // Iterate over the map and access its elements
            for (const auto& pair : myMap) {
                
                const std::vector<double>& values = pair.second;
         
                cout<<values[0]<<" "<<values[1]<<endl;
                
            }
        } else {
            // Handle the case when the map is empty or nullptr
            std::cout << "Error: GetMap() returned an empty map or nullptr." << std::endl;
        }
 
        cout<<"Costo percorso "<<i<<": "<<M_vc[i].Costo()<<endl<<endl;
        cout<<endl;
    }
    
}



