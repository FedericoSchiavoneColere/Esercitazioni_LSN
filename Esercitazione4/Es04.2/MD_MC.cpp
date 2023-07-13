#include <iostream>
#include <fstream>
#include <ostream>
#include <cmath>
#include <iomanip>
#include <string>
#include "MD_MC.h"

using namespace std;

int main(){ //Esercizio 7

    string phase="solid"; //da cambiare tra {gas; liquid; solid}
    
    
    Input(phase); //Inizialization
    
    if(restart==0){
        Equilibra();
        cout<<endl<<"Ho equilibrato, fare ripartire dopo aver modificato file di input"<<endl<<endl;
        return 0;
    }
    if(restart!=0){
        int nconf = 1;
        for(int iblk=1; iblk <= nblk; iblk++) //Simulation //data blocking
        {
            Reset(iblk);   //Reset block averages
            for(int istep=1; istep <= nstep; istep++) //numero di ste per blocco
            {
                Move(); //muovo le particelle
                Measure();
                Accumulate(); //Update block averages
                if(istep%10 == 0){ //ogni 10 step scrivo la configurazione delle particelle
                    //ConfXYZ(nconf);//Write actual configuration in XYZ format //Commented to avoid "filesystem full"!
                    nconf += 1;
                }
            }
            Averages(iblk,phase);   //Print results for current block //l'ultima volta fara anche il calcolo dei valori finali
        }
        ConfFinal(); //Write final configuration
        Measure();
        
        cout << "Final potential energy = " << walker[iv]/(double)npart << endl;
        cout << "Final temperature      = " << walker[it] << endl;
        cout << "Final kinetic energy   = " << walker[ik]/(double)npart << endl;
        cout << "Final total energy     = " << walker[ie]/(double)npart << endl;
        cout << "Final pressure     = " << walker[ip] << endl;
    }
  return 0;
}


void Input(string p)
{
  ifstream ReadInput, ReadConf, ReadVelocity, Primes, Seed;

  cout << "Classic Lennard-Jones fluid        " << endl;
  cout << "MD(NVE) / MC(NVT) simulation       " << endl << endl;
  cout << "Interatomic potential v(r) = 4 * [(1/r)^12 - (1/r)^6]" << endl << endl;
  cout << "Boltzmann weight exp(- beta * sum_{i<j} v(r_ij) ), beta = 1/T " << endl << endl;
  cout << "The program uses Lennard-Jones units " << endl;

//Read seed for random numbers
  int p1, p2;
  Primes.open("Primes");
  Primes >> p1 >> p2 ;
  Primes.close();

//Read input informations
  ReadInput.open("input."+p);

  ReadInput >> iNVET; //0 = MD (NVE); 1 = MC (NVT) dinamica molecolare o monte carlo
  ReadInput >> restart; //0 No restart; 1 = restart

  if(restart) Seed.open("seed.out");
  else Seed.open("seed.in");
  Seed >> seed[0] >> seed[1] >> seed[2] >> seed[3];
  rnd.SetRandom(seed,p1,p2);
  Seed.close();

  ReadInput >> temp;
  beta = 1.0/temp;
  cout << "Temperature = " << temp << endl;

  ReadInput >> npart;
  cout << "Number of particles = " << npart << endl;

  ReadInput >> rho;
  cout << "Density of particles = " << rho << endl;
  vol = (double)npart/rho;
  box = pow(vol,1.0/3.0);
  cout << "Volume of the simulation box = " << vol << endl;
  cout << "Edge of the simulation box = " << box << endl;

  ReadInput >> rcut;
  cout << "Cutoff of the interatomic potential = " << rcut << endl << endl;
    
  ReadInput >> delta;

  ReadInput >> nblk;

  ReadInput >> nstep;

  cout << "The program perform Metropolis moves with uniform translations" << endl;
  cout << "Moves parameter = " << delta << endl;
  cout << "Number of blocks = " << nblk << endl;
  cout << "Number of steps in one block = " << nstep << endl << endl;
  cout << "File: input."<<p<<endl;
  ReadInput.close();

//Prepare arrays for measurements
  iv = 0; //Potential energy
  it = 1; //Temperature
  ik = 2; //Kinetic energy
  ie = 3; //Total energy
  ip = 4; //Pressione
  n_props = 5; //Number of observables. dovrà diventare 5 dopo introduzione della pressione

//Read initial configuration
  cout << "Read initial configuration" << endl << endl;
  if(restart)
  {
    ReadConf.open("config.out");
    ReadVelocity.open("velocity.out");
    for (int i=0; i<npart; ++i) ReadVelocity >> vx[i] >> vy[i] >> vz[i];
  }
  else 
  {
    ReadConf.open("config.in");
    cout << "Prepare velocities with center of mass velocity equal to zero " << endl;
    double sumv[3] = {0.0, 0.0, 0.0};
    for (int i=0; i<npart; ++i)
    {
      vx[i] = rnd.Gauss(0.,sqrt(temp)); //maxwell-boltzmann velocity distribution
      vy[i] = rnd.Gauss(0.,sqrt(temp));
      vz[i] = rnd.Gauss(0.,sqrt(temp));
      sumv[0] += vx[i];                     //compute drift velocity
      sumv[1] += vy[i];
      sumv[2] += vz[i];
    }
    for (int idim=0; idim<3; ++idim) sumv[idim] /= (double)npart;
    double sumv2 = 0.0, fs;
    for (int i=0; i<npart; ++i)
    {
      vx[i] = vx[i] - sumv[0];              //subtract drift velocity per particle
      vy[i] = vy[i] - sumv[1];              //siccome uso maxwell boltzmann v_tot diverso da zero
      vz[i] = vz[i] - sumv[2];              //stato sta traslando e io non voglio
      sumv2 += vx[i]*vx[i] + vy[i]*vy[i] + vz[i]*vz[i];
    }
    sumv2 /= (double)npart;
    fs = sqrt(3 * temp / sumv2);   // fs = velocity scale factor  //v*=pow(3T*,0.5) da termodinamica
    cout << "velocity scale factor: " << fs << endl << endl;
    for (int i=0; i<npart; ++i)
    {
      vx[i] *= fs;
      vy[i] *= fs;
      vz[i] *= fs;
    }
  }

  for (int i=0; i<npart; ++i)
  {
    ReadConf >> x[i] >> y[i] >> z[i];       //read configuration (side units) //posso usare questo sia per gas liquido e solido
    x[i] = Pbc( x[i] * box );               //Scale configuration (reduced units) //cambia il fattore di scala tra gli stati
    y[i] = Pbc( y[i] * box );
    z[i] = Pbc( z[i] * box );               //Pbc( ) per prudenza, se ci sono particelle fuori dalla scatola
                                            //Pbc le rimette nella cella elementare
  }
  ReadConf.close();

  for (int i=0; i<npart; ++i)
  {
    if(iNVET)
    {
      xold[i] = x[i];
      yold[i] = y[i];
      zold[i] = z[i];
    }
    else            //MD simulation
    {
      xold[i] = Pbc(x[i] - vx[i] * delta);          //compuite old configuration from velocities
      yold[i] = Pbc(y[i] - vy[i] * delta);
      zold[i] = Pbc(z[i] - vz[i] * delta);
    }
  }
  
//
    vtail = 8*pi*rho*(1.0/(9*pow(rcut,9))-1.0/(3*pow(rcut,3)))*npart;
    wtail = 32*pi*rho*(1.0/(9*pow(rcut,9))-1.0/(6*pow(rcut,3)))*3*npart;
//
    
  
//Evaluate properties of the initial configuration
  Measure();

//Print initial values for measured properties
  cout << "Initial potential energy = " << walker[iv]/(double)npart << endl;
  cout << "Initial temperature      = " << walker[it] << endl;
  cout << "Initial kinetic energy   = " << walker[ik]/(double)npart << endl;
  cout << "Initial total energy     = " << walker[ie]/(double)npart << endl;
  cout << "Initial pressure     = " << walker[ip] << endl;

  return;
}


void Move()
{
  int o;
  double p, energy_old, energy_new;
  double xnew, ynew, znew;

  if(iNVET) // Monte Carlo (NVT) move
  {
    for(int i=0; i<npart; ++i)
    {
    //Select randomly a particle (for C++ syntax, 0 <= o <= npart-1)
      o = (int)(rnd.Rannyu()*npart);

    //Old
      energy_old = Boltzmann(x[o],y[o],z[o],o);

    //New
      x[o] = Pbc( x[o] + delta*(rnd.Rannyu() - 0.5) );
      y[o] = Pbc( y[o] + delta*(rnd.Rannyu() - 0.5) );
      z[o] = Pbc( z[o] + delta*(rnd.Rannyu() - 0.5) );

      energy_new = Boltzmann(x[o],y[o],z[o],o);

    //Metropolis test
      p = exp(beta*(energy_old-energy_new));
      if(p >= rnd.Rannyu())  
      {
      //Update
        xold[o] = x[o];
        yold[o] = y[o];
        zold[o] = z[o];
          accepted++;
      } else {
        x[o] = xold[o];
        y[o] = yold[o];
        z[o] = zold[o];
      }
        attempted++;
    }
  } else // Molecular Dynamics (NVE) move
  {
    double fx[m_part], fy[m_part], fz[m_part];

    for(int i=0; i<npart; ++i){ //Force acting on particle i
      fx[i] = Force(i,0);       //devo calcolare le forze che agiscono su particelle. ciclo su particelle
      fy[i] = Force(i,1);       //ho calcolato Force su paticelella i-esima su 3 coordinate
      fz[i] = Force(i,2);
    }

    for(int i=0; i<npart; ++i){ //Verlet integration scheme

      xnew = Pbc( 2.0 * x[i] - xold[i] + fx[i] * pow(delta,2) );    //nuove posizioni
      ynew = Pbc( 2.0 * y[i] - yold[i] + fy[i] * pow(delta,2) );
      znew = Pbc( 2.0 * z[i] - zold[i] + fz[i] * pow(delta,2) );

      vx[i] = Pbc(xnew - xold[i])/(2.0 * delta);                //nuove velocità
      vy[i] = Pbc(ynew - yold[i])/(2.0 * delta);
      vz[i] = Pbc(znew - zold[i])/(2.0 * delta);

      xold[i] = x[i];                        //registro perche siamo ad un time-step successivo
      yold[i] = y[i];
      zold[i] = z[i];

      x[i] = xnew;                           //nuove posizioni
      y[i] = ynew;
      z[i] = znew;

      accepted = accepted + 1.0;
      attempted = attempted + 1.0;
    }
  }
  return;
}

double Boltzmann(double xx, double yy, double zz, int ip)
{
  double ene=0.0;
  double dx, dy, dz, dr;

  for (int i=0; i<npart; ++i)
  {
    if(i != ip)
    {
// distance ip-i in pbc
      dx = Pbc(xx - x[i]);
      dy = Pbc(yy - y[i]);
      dz = Pbc(zz - z[i]);

      dr = dx*dx + dy*dy + dz*dz;
      dr = sqrt(dr);

      if(dr < rcut)
      {
        ene += 1.0/pow(dr,12) - 1.0/pow(dr,6);
      }
    }
  }

  return 4.0*ene;
}

double Force(int ip, int idir){ //Compute forces as -Grad_ip V(r)
  double f=0.0;
  double dvec[3], dr;
                                    
  for (int i=0; i<npart; ++i){
    if(i != ip){                    //ciclo su tutte le particelle tranne sè stessa
      dvec[0] = Pbc( x[ip] - x[i] );  // distance ip-i in pbc
      dvec[1] = Pbc( y[ip] - y[i] );
      dvec[2] = Pbc( z[ip] - z[i] );

      dr = dvec[0]*dvec[0] + dvec[1]*dvec[1] + dvec[2]*dvec[2];
      dr = sqrt(dr);

      if(dr < rcut){ //incremento la forza in una direzione
        f += dvec[idir] * (48.0/pow(dr,14) - 24.0/pow(dr,8)); // -Grad_ip V(r)
      }
    }
  }
 
  return f;
}

void Measure() //Properties measurement //misura le 4 propiretà a cui andrà aggiuna la pressione
{
  double v = 0.0, w=0.0 ,kin=0.0;
  double vij,wij;
  double dx, dy, dz, dr;
    for(int i=0; i<ng; i++) gdr[i] = 0; // azzero l'istogramma
    
    
//cycle over pairs of particles //tutte le coppie non ripetute
  for (int i=0; i<npart-1; ++i)
  {
      for (int j=i+1; j<npart; ++j)
      {
          // distance i-j in pbc
          dx = Pbc(x[i] - x[j]);
          dy = Pbc(y[i] - y[j]);
          dz = Pbc(z[i] - z[j]);
          
          dr = dx*dx + dy*dy + dz*dz;
          dr = sqrt(dr);
          
          if(dr < rcut) //se è minore della distanza di cutoff faccio calcolo
          {
              vij = 1.0/pow(dr,12) - 1.0/pow(dr,6); //il 4* lo metto dopo
              wij = 48.0/pow(dr,12) - 24.0/pow(dr,6); // viriale
              v += vij;
              w += wij;
          }
          
          // riempio l'istogramma della gdr
          min_dist = box/2.0; // risoluzione dell'istogramma della gdr
          bin_index = ng*(dr/min_dist); // posizione nell'istogramma
          
          if(dr<min_dist){
              gdr[ bin_index ] += 2; // gdr
          }
      }
   }
    
    //misuro energia cinetica
  for (int i=0; i<npart; ++i) kin += 0.5 * (vx[i]*vx[i] + vy[i]*vy[i] + vz[i]*vz[i]);

   // walker[iv] = (4.0 *v + (8*M_PI*rho)/(9*pow(rcut,9))-(8*M_PI*rho)/(3*pow(rcut,3))/(double)npart); // Potential energy //qui c'è il 4* e il v_tail
  walker[iv]=4.0*v+vtail;
  walker[ik] = kin/(double)npart; // Kinetic energy
  walker[it] = (2.0 / 3.0) * kin/(double)npart; // Temperature
  walker[ie] = 4.0 * v + kin;  // Total energy;
  //walker[ip] = rho * walker[it] + (w + 32*M_PI*rho *(1/(9*pow(rcut,9))-1/(6*pow(rcut,3))))/(3.0*vol);  //Pressione con tail correction del viriale
  walker[ip]= rho*temp+(w+wtail)/(3.0*vol);
  
    return;
}


void Reset(int iblk) //Reset block averages
{
   
    if(iblk == 1){
        for(int i=0; i<n_props; ++i){ //per ogni propiretà
        
            glob_av[i] = 0;          //accumulatore valori medi
            glob_av2[i] = 0;         //accumulatore valori quadratici
        }
        for(int i=0; i<ng; i++){ //per la gdr
            gdr[i]=0;
            gdr_ave[i]=0;
        }
    }

    for(int i=0; i<n_props; ++i){
        blk_av[i] = 0;         //azzero le medie di blocco
    }
    for(int i=0; i<ng; i++){ //<<<<<<< gdr
        gdr[i]=0;
        gdr_ave[i]=0;
    }
    blk_norm = 0;
    attempted = 0;
    accepted = 0;
}


void Accumulate(void) //Update block averages
{
    for (int i=0; i<ng; i++)
       {
        gdr_ave[i]=gdr_ave[i]+gdr[i];
        gdr[i]=0;
       }

   for(int i=0; i<n_props; ++i)
   {
     blk_av[i] = blk_av[i] + walker[i];
   }
   blk_norm = blk_norm + 1.0;
}


void Averages(int iblk, string p) //Print results for current block
{
    //utilizza quello calcolato nel block avareges per fare valore medio e quadratico
    ofstream Epot, Pres, Gdr;
    ofstream Ekin, Etot, Temp;
   const int wd=12;
    
    cout << "Block number " << iblk << endl;
    cout << "Acceptance rate " << accepted/attempted << endl << endl;
    
    Epot.open("./"+p+"/output_epot.dat",ios::app);
    Ekin.open("./"+p+"/output_ekin.dat",ios::app);
    Temp.open("./"+p+"/output_temp.dat",ios::app);
    Etot.open("./"+p+"/output_etot.dat",ios::app);
    Pres.open("./"+p+"/output_pres.dat",ios::app);
    Gdr.open("./"+p+"/output_gdr.dat", ios::app);
    
    stima_pot = blk_av[iv]/blk_norm/(double)npart; //Potential energy //per particella
    glob_av[iv] += stima_pot;
    glob_av2[iv] += stima_pot*stima_pot;
    err_pot=Error(glob_av[iv],glob_av2[iv],iblk); //incertezza sapendo numero di blocchi
    
    
    stima_kin = blk_av[ik]/blk_norm/(double)npart; //Kinetic energy
    glob_av[ik] += stima_kin;
    glob_av2[ik] += stima_kin*stima_kin;
    err_kin=Error(glob_av[ik],glob_av2[ik],iblk);

    stima_etot = blk_av[ie]/blk_norm/(double)npart; //Total energy
    glob_av[ie] += stima_etot;
    glob_av2[ie] += stima_etot*stima_etot;
    err_etot=Error(glob_av[ie],glob_av2[ie],iblk);

    stima_temp = blk_av[it]/blk_norm; //Temperature
    glob_av[it] += stima_temp;
    glob_av2[it] += stima_temp*stima_temp;
    err_temp=Error(glob_av[it],glob_av2[it],iblk);

    stima_pres = blk_av[ip]/blk_norm; //Pressione
    glob_av[ip] += stima_pres;
    glob_av2[ip] += stima_pres*stima_pres;
    err_pres=Error(glob_av[ip],glob_av2[ip],iblk);

    
    double gdr_norm=0;              // GDR
       for(int i=0; i<ng; i++)         // normalizzo le medie della gdr
       {
         double r = i*(min_dist)/ng; // distanza tra due particelle
         double dr =  (min_dist)/ng;  // risoluzione della gdr
         gdr_norm = rho*npart*4.0*M_PI/3.0*(pow(r+dr, 3)-pow(r, 3));
         gdr_ave[i] = double(gdr_ave[i]/(double)blk_norm)/gdr_norm; // medie su gdr.
       }
    
    //Gdr
      for(int i=0; i<ng; i++){
        Gdr << gdr_ave[i] << " "; // stampo le medie di blocco su ogni bin
      }
      Gdr << endl;

    
//Potential energy per particle //a che blocco, stima nel blocco, stimo fino ad ora ed errore
    Epot  << setw(wd)<< iblk <<  setw(wd) << stima_pot << setw(wd) << glob_av[iv]/(double)iblk << setw(wd) << err_pot << endl;
//Kinetic energy
    Ekin << setw(wd) << iblk <<  setw(wd) << stima_kin << setw(wd) << glob_av[ik]/(double)iblk << setw(wd) << err_kin << endl;
//Total energy
    Etot << setw(wd) << iblk <<  setw(wd) << stima_etot << setw(wd) << glob_av[ie]/(double)iblk << setw(wd) << err_etot << endl;
//Temperature
   Temp << setw(wd) << iblk <<  setw(wd) << stima_temp << setw(wd) << glob_av[it]/(double)iblk << setw(wd) << err_temp << endl;
//Pressione
    Pres  << setw(wd)<< iblk <<  setw(wd) << stima_pres << setw(wd) << glob_av[ip]/(double)iblk << setw(wd) << err_pres << endl;

    cout << "----------------------------" << endl << endl;

    Epot.close();
    Ekin.close();
    Etot.close();
    Temp.close();
    Pres.close();
    Gdr.close();
}


void ConfFinal(void) //scrive vel e posizone nel formato x y z
{
  ofstream WriteConf, WriteVelocity, WriteSeed;

  cout << "Print final configuration to file config.out" << endl << endl;
  WriteConf.open("config.out");
  WriteVelocity.open("velocity.out");
  for (int i=0; i<npart; ++i)
  {
    WriteConf << x[i]/box << "   " <<  y[i]/box << "   " << z[i]/box << endl;
    WriteVelocity << vx[i] << "   " <<  vy[i] << "   " << vz[i] << endl;
  }
  WriteConf.close();
  WriteVelocity.close();

  rnd.SaveSeed();
}

void ConfXYZ(int nconf){ //Write configuration in .xyz format
  ofstream WriteXYZ;

  WriteXYZ.open("frames/config_" + to_string(nconf) + ".xyz");
  WriteXYZ << npart << endl;
  WriteXYZ << "This is only a comment!" << endl;
  for (int i=0; i<npart; ++i){
    WriteXYZ << "LJ  " << Pbc(x[i]) << "   " <<  Pbc(y[i]) << "   " << Pbc(z[i]) << endl;
  }
  WriteXYZ.close();
}

double Pbc(double r)  //Algorithm for periodic boundary conditions with side L=box
{
    return r - box * rint(r/box);
}

double Error(double sum, double sum2, int iblk)
{
    return sqrt(fabs(sum2/(double)iblk - pow(sum/(double)iblk,2))/(double)iblk);
}

void Equilibra(void)
{
    int N_step=4000;
    int I=100;
    for(int iblk=1; iblk <= N_step; iblk++) //Simulation //data blocking
    {
        Reset(iblk);   //Reset block averages
        cout<<"Passo di equilibrautra: "<<iblk<<"/"<<N_step<<endl;
        for(int istep=1; istep <= I; istep++) //numero di step per blocco
        {
            Move(); //muovo le particelle
            
        }
    }
    ConfFinal(); //Write configuration
}
