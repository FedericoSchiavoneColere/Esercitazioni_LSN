
#include "Monte_Carlo_ISING_1D.h"

int main()
{
    
    //Equilibro();
    Input(); //Inizialization
    
    
    
    int N_step_temp=50;
    for(int i=0;i<N_step_temp;i++){
        temp=i*0.1 +0.1;
        beta=1/temp;
        
        for(int i=0; i<10000;i++){
            Move(metro);
        }
        
        cout<<endl<<"Simulazione "<<i+1<<"/"<<N_step_temp<<endl;
        cout<<"T = "<<temp<<endl<<endl;
        
        
        stringstream stream;
        stream << fixed << std::setprecision(2) << temp;
        string filename = "output_" + stream.str();
        for(int iblk=1; iblk <= nblk; ++iblk) //Simulation
        {
                Reset(iblk);   //Reset block averages
            
                for(int istep=1; istep <= nstep; ++istep)
                {
                        Move(metro);
                        Measure();
                        Accumulate(); //Update block averages
                }
                Averages(iblk, filename);   //Print results for current block
        }
        ConfFinal(filename); //Write final configuration
    }

return 0;
}



void Input(void)
{
    ifstream ReadInput;

    cout << "Classic 1D Ising model             " << endl;
    cout << "Monte Carlo simulation             " << endl << endl;
    cout << "Nearest neighbour interaction      " << endl << endl;
    cout << "Boltzmann weight exp(- beta * H ), beta = 1/T " << endl << endl;
    cout << "The program uses k_B=1 and mu_B=1 units " << endl;

    //Read seed for random numbers
    int p1, p2;
    ifstream Primes("Primes");
    Primes >> p1 >> p2 ;
    Primes.close();

    ifstream input("seed.in");
    input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
    rnd.SetRandom(seed,p1,p2);
    input.close();
  
    //Read input informations
    ReadInput.open("input.txt");

    ReadInput >> temp;
    //beta = 1.0/temp;
    cout << "Temperature = " << temp << endl;

    ReadInput >> nspin;
    cout << "Number of spins = " << nspin << endl;

    ReadInput >> J;
    cout << "Exchange interaction = " << J << endl;

    ReadInput >> h;
    cout << "External field = " << h << endl << endl;
    
    ReadInput >> metro; // if=1 Metropolis else Gibbs

    ReadInput >> nblk;

    ReadInput >> nstep;
    
    if(metro==1) cout << "The program perform Metropolis moves" << endl;
    else cout << "The program perform Gibbs moves" << endl;
    cout << "Number of blocks = " << nblk << endl;
    cout << "Number of steps in one block = " << nstep << endl << endl;
    ReadInput.close();


    //Prepare arrays for measurements
    iu = 0; //Energy
    ic = 1; //Heat capacity
    im = 2; //Magnetization
    ix = 3; //Magnetic susceptibility
 
    n_props = 4; //Number of observables

    //initial configuration
    for (int i=0; i<nspin; ++i)
    {
        if(rnd.Rannyu() >= 0.5) s[i] = 1;
        else s[i] = -1;
    }
  
    //Evaluate energy etc. of the initial configuration
    Measure();

    //Print initial values for the potential energy and virial
    cout << "Initial energy = " << walker[iu]/(double)nspin << endl;

}


void Move(int metro)
{
   
    for(int i=0; i<nspin; ++i)
    {
            int o=0;
            
            for(int i=0; i<nspin; ++i) // ciclo sul numero di spin totali
            {
                //Select randomly a particle (for C++ syntax, 0 <= o <= nspin-1)
                o = (int)(rnd.Rannyu()*nspin);
                
                if(metro==1)
                {
                    attempted++;
                    // flip
                    int s_try = -s[o];
                    
                    // calcolo variazione energia
                    double DeltaE = 2*Boltzmann(s_try, o);
                    
                    double q = exp(-beta*DeltaE);
                    double P = min(1.,q);
                    
                    if (P==1) // accetto
                    {
                        s[o] = s_try;
                        accepted++;
                    }
                    else // accetto con probabilità P
                    {
                        if(rnd.Rannyu() < P)
                        {
                            s[o] = s_try;
                            accepted++;
                        }
                    }
                }
                
                else //Gibbs
                {
                    attempted++;
                    int snew = ((int)rnd.Rannyu(0,2))*2-1;
                    
                    // variazione energia DeltaE e probabilità
                    double DeltaE = -2*Boltzmann(snew, o);
                    double p = 1.0/(1.0+exp(-beta*DeltaE));
                    
                    //flip con probabilità p
                    if(rnd.Rannyu() < p)
                    {
                        s[o] = snew;
                    } else s[o] = -snew;
                    accepted++;
            }
        }
    }
}
double Boltzmann(int sm, int ip)
{
    double ene = -J * sm * ( s[Pbc(ip-1)] + s[Pbc(ip+1)] ) - h * sm;
    return ene;
}

void Measure()
{
    // int bin;
    double H = 0.0, u = 0.0, m = 0.0; // x = 0.0, c = 0.0;

    //cycle over spins
    for (int i=0; i<nspin; ++i)
    {
        H = -J * s[i] * s[Pbc(i+1)] - 0.5 * h * (s[i] + s[Pbc(i+1)]);

        u += H;
        //c += H*H; // prima era sbagliato!
        m += s[i];
        // x +=  s[i]*s[i]; // non serve, perchè h = 0
    }
    // energia
    walker[iu] = u;
    // capacità termica
    walker[ic] = pow(u,2); //il resto calcolato in averages
    // magnetizzazione
    walker[im] = m;
    // suscettività
    walker[ix] = beta*pow(m,2);
}


void Reset(int iblk) //Reset block averages
{
   
    if(iblk == 1){
        
        for(int i=0; i<n_props; ++i)
        {
            glob_av[i] = 0;
            glob_av2[i] = 0;
        }
    }

    for(int i=0; i<n_props; ++i){
        blk_av[i] = 0;
    }
    blk_norm = 0;
    attempted = 0;
    accepted = 0;
}


void Accumulate(void) //Update block averages
{

    for(int i=0; i<n_props; ++i){
        blk_av[i] = blk_av[i] + walker[i];
    }
    blk_norm = blk_norm + 1.0;
}


void Averages(int iblk, string name) //Print results for current block
{
    ofstream Ene, Heat, Mag, Chi;
    
    cout << "Block number " << iblk << "  |  ";
    cout << "Acceptance rate " << accepted/attempted << endl;
    
    // ENERGIA
    Ene.open("output/"+name+"_ene.txt",ios::app);
    stima_u = blk_av[iu]/blk_norm/(double)nspin;
    glob_av[iu]  += stima_u; // accu globale
    glob_av2[iu] += stima_u*stima_u; // accu quadratico

    err_u=Error(glob_av[iu],glob_av2[iu],iblk);
    Ene << iblk <<  " " << stima_u << " " << glob_av[iu]/(double)iblk << " " << err_u << endl;
    Ene.close();

    // CAPACITÀ TERMICA
    Heat.open("output/"+name+"_heat.txt",ios::app);
    stima_c = beta*beta * (blk_av[ic]/blk_norm - pow(blk_av[iu]/blk_norm,2))/(double)nspin;
    glob_av[ic]  += stima_c; // accu globale
    glob_av2[ic] += stima_c*stima_c; // accu quadratico

    err_c=Error(glob_av[ic],glob_av2[ic],iblk);
    Heat << iblk  <<  " " << stima_c << " " <<  glob_av[ic]/(double)iblk << " " << err_c << endl;
    Heat.close();

    // MAGNETIZZAZIONE
    Mag.open("output/"+name+"_mag.txt",ios::app);
    stima_m = blk_av[im]/blk_norm/(double)nspin;
    glob_av[im]  += stima_m; // accu globale
    glob_av2[im] += stima_m*stima_m; // accu quadratico

    err_m=Error(glob_av[im],glob_av2[im],iblk);
    Mag << iblk <<  " " << stima_m << " " << glob_av[im]/(double)iblk << " " << err_m << endl;
    Mag.close();

    // SUSCETTIVITÀ
    Chi.open("output/"+name+"_chi.txt",ios::app);
    stima_x = blk_av[ix]/blk_norm/(double)nspin;
    glob_av[ix]  += stima_x; // accu globale
    glob_av2[ix] += stima_x*stima_x; // accu quadratico

    err_x=Error(glob_av[ix],glob_av2[ix],iblk);
    Chi << iblk <<  " " << stima_x << " " << glob_av[ix]/(double)iblk << " " << err_x << endl;
    Chi.close();
}


void ConfFinal(string name)
{
    ofstream WriteConf;
    
    cout << "Print final configuration to file config.final " << endl << endl;
    WriteConf.open("./output/"+name+"conf_fin.txt");
    for (int i=0; i<nspin; ++i){

        WriteConf << s[i] << endl;
    }
    WriteConf.close();

    rnd.SaveSeed();
}


int Pbc(int i)  //Algorithm for periodic boundary conditions
{
    if(i >= nspin) i = i - nspin;
    else if(i < 0) i = i + nspin;
    return i;
}

double Error(double sum, double sum2, int iblk)
{
    if(iblk==1) return 0.0;
    else return sqrt((sum2/(double)iblk - pow(sum/(double)iblk,2))/(double)(iblk-1));
}

void Equilibro()
{
    for(int i=0; i<100000;i++){
        Move(metro);
    }
}
