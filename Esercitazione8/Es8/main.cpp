#include "main.h"

using namespace std;

ofstream WriteResult;
ofstream WritePos;

ofstream WriteTraj;


int main (){
    
    // Generatore di random
    //Random rnd;
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
    
    
    Input();
    
    
    if(SA==0)
    {
        MeasureAverageH(SA);
        
    }

    if(SA==1)
    {
        WriteTraj.open("./SA/traj.out");


        int eq = 10;
        for(int i = 0; i < eq; i++)
        {
            MetroMove();
            cout << "equilibro: " << i+1 << "/" << eq << endl << endl;
        }

        H = 0;
        attempted = 0;
        accepted = 0;


        while (2*Lmu/beta > errMu or 2*Lsgm/beta > errSgm)
        {
            
            for(int i = 0; i < step_in_beta; i++)
            {
                oldH = mean_prog_H/nblk;

                d_mu = rnd.Rannyu(-Lmu,Lmu)/beta;
                d_sgm = rnd.Rannyu(-Lsgm,Lsgm)/beta;

                mu += d_mu;
                sgm += d_sgm;

                MeasureAverageH(SA);
                attempted++;
                
                newH = mean_prog_H/nblk; //nuova probabilità
                errH = error(mean_prog_H/nblk,var_prog_H/nblk,(nblk-1));
                deltaH = oldH - newH;

                q = exp(beta*(deltaH));
                A = min(1,q);

                //metropolis
                if (A==1) // accetto
                {
                    oldH = newH;
                    accepted++;
                }
                else // accetto con probabilità A
                {
                    if(rnd.Rannyu() < A)
                    {
                        oldH = newH;
                        accepted++;
                    } else {
                        mu -= d_mu;
                        sgm -= d_sgm;
                    }
                }
            
                step++;

            }

            cout<<"step:    "<<step<<endl;
            cout<<"beta:    "<<beta<<endl;
            cout<<"H:       "<<oldH<<endl;
            cout<<"errH:    "<<errH<<endl;
            cout<<"mu:      "<<mu<<endl;
            cout<<"errMu:   "<<2*Lmu/beta<<endl;
            cout<<"sgm:     "<<sgm<<endl;
            cout<<"errSgm:  "<<2*Lsgm/beta<<endl;
            cout<<endl;
            
            WriteTraj<<beta<<" "<<mu<<" "<<sgm<<" "<<oldH<<" "<<errH<<endl;

            // abbasso temperatura
            beta += db;
            

            if (2*Lmu/beta <= errMu and 2*Lsgm/beta <= errSgm){

                cout<<"Mi muovo un po' fissata la temperatura più bassa, per valutare la STD di mu e sigma."<<endl;
            
                for(int i = 0; i < step_in_beta*amp; i++)
                {
                    if (i%1000 == 0) cout << i << "/" << step_in_beta*amp << endl;

                    oldH = mean_prog_H/nblk;

                    //esploro mu e sgm in funzione di beta
                    d_mu = rnd.Rannyu(-Lmu,Lmu)/beta; // a basse T >> alte beta, l'esplorazione è meno ampia.
                    d_sgm = rnd.Rannyu(-Lsgm,Lsgm)/beta;

                    mu += d_mu;
                    sgm += d_sgm;

                    MeasureAverageH(SA);
                    attempted++;
                
                    // calcolo probabilità
                    newH = mean_prog_H/nblk;
                    errH = error(mean_prog_H/nblk,
                             var_prog_H/nblk,
                             (nblk-1));
                    deltaH = oldH - newH;

                    q = exp(beta*(deltaH)); // con +, perchè devo andare verso il minimo di H
                    A = min(1,q);

                    // giudico con metropolis, confrontantdo newH con oldH
                    if (A==1) // accetto direttamente
                    {
                        oldH = newH;
                        accepted++;
                    }
                    else // accetto con probabilità A
                    {
                        if(rnd.Rannyu() < A)
                        {
                            oldH = newH;
                            accepted++;
                        } else {mu -= d_mu; sgm -= d_sgm;} // se non va, ripristino i parametri
                    }
                    step++;

                    // calcolo errore su sgm: <<< attenzione ragionare sulla correlazione
                    mean_sgm += sgm;
                    var_sgm += sgm*sgm;
                    mean_mu += mu;
                    var_mu += mu*mu;
                }

                mean_sgm /= step_in_beta*amp;
                var_sgm /= step_in_beta*amp;
                mean_mu /= step_in_beta*amp;
                var_mu /= step_in_beta*amp;
            }

            
        }
        cout<<endl;
        cout<<"Minimo di H: "<<oldH<< " ± "<<errH<<endl;
        cout<<"mu:          "<<mean_mu<<" ± "<<error(mean_mu, var_mu, step_in_beta*amp)<<endl;
        cout<<"sigma:       "<<mean_sgm<<" ± "<<error(mean_sgm,var_sgm,step_in_beta*amp)<<endl<<endl;
    }
    WriteTraj.close();

    return 0;
}



void MeasureAverageH(int q)
{
    string Name;
    
    if(q==0)Name="H";
    if(q==1)Name="SA";
    
    WriteResult.open("./"+Name+"/result.out");
    WritePos.open("./"+Name+"/pos.out");

    mean_prog_H = 0;
      var_prog_H = 0;

      // ciclo sui blocchi
      for(int i = 0; i < nblk; i++)
      {
          H = 0; // azzero H prima di iniziare il blocco
          attempted = 0;
          accepted = 0;

          // cliclo nel singolo blocco
          for (int j = 0; j < steps; j++)
          {
              MetroMove();
              
              // accumulo H nel blocco
              H += EvalH(y);
              
              // salvo dove sono
              if(SA==0 && j%(steps/histofill_blk) == 0) WritePos << y << endl;
              
          }

          // stampo acceptance rate
          if(SA==0)
          {
              cout << "Block # " << i+1 << endl;
              cout << "Acceptance rate:   " << (double)accepted/attempted << endl;
              cout << "-----------------------------------" << endl;
          }
                      

          // medie di blocco
          H /= steps;
          mean_prog_H += H;
          var_prog_H += H*H;

          // salvo su file
          if(SA==0)
          {
          if (WriteResult.is_open())
          {
              if(i == 0) WriteResult << H << " " << mean_prog_H/(i+1) << " " << 0 << endl;
              else WriteResult << H << " " << mean_prog_H/(i+1) << " " << error(mean_prog_H/(i+1), var_prog_H/(i+1), i) << " " << endl;
          } else {
              if (!WriteResult.is_open()) cerr << "PROBLEM: Unable to open result.out" << endl;
          }
          }
      }

      WriteResult.close();
      WritePos.close();
  }


void MetroMove()
{
    dx = {rnd.Rannyu(-span,span)}; // uniform transition probability in un cubo di lato L
            
    // provo a spostarmi (x = posizione nuova)
    x = y + dx;
    attempted++;

    // calcolo prob e accetto secondo metropolis
    q = Psi2(x)/Psi2(y);
    A = min(1,q);

    if (A==1) // accetto direttamente
    {
        y = x;
        accepted++;
    }
    else // accetto con probabilità A
    {
        if(rnd.Rannyu() < A)
        {
            y = x;
            accepted++;
        }
    }
}

double EvalH(double x)
{
    double exp1 = exp(-(x-mu)*(x-mu)/(2*sgm*sgm));
    double exp2 = exp(-(x+mu)*(x+mu)/(2*sgm*sgm));

    double psi = exp1+exp2;
    double V = pow(x,4) - 2.5*pow(x,2);

    double Hpsi_pot = V*psi;
    double Hpsi_kin = -0.5*( exp1*(pow(x-mu,2)/pow(sgm,4)-(pow(sgm,-2)))
                           + exp2*(pow(x+mu,2)/pow(sgm,4)-(pow(sgm,-2)))
                           );


    double Hpsi = Hpsi_kin + Hpsi_pot;
    return Hpsi/psi;
}

double Psi2(double x)
{

    return pow(exp(-(x-mu)*(x-mu)/(2*sgm*sgm))+exp(-(x+mu)*(x+mu)/(2*sgm*sgm)),2);
    
}

void Input(void)
{
  ifstream ReadInput;

    
//Read input informations
  ReadInput.open("dati.txt");
    
  ReadInput >> nblk;
  ReadInput >> nstep;
  ReadInput >> L;
  ReadInput >> x0;
  ReadInput >> mu;
  ReadInput >> sgm;
  ReadInput >> SA;
  ReadInput >> histofill_blk;

  ReadInput >> beta;
  ReadInput >> db;
  ReadInput >> step_in_beta;

  ReadInput >> Lmu;     // larghezza iniziale passi mu
  ReadInput >> Lsgm;    // idem per sigma
  ReadInput >> errMu; // errore con cui voglio determinare mu
  ReadInput >> errSgm; // errore con cui voglio determinare mu

  ReadInput >> manualseed;
  int pr;
  ReadInput >> pr;

  
  x = x0;
  steps = nstep/nblk;
  span = L/2.0;

  cout<<"Total number of steps = "<<nstep<<endl;
  cout<<"Number of blocks = "<<nblk<<endl;
  cout<<"Number of steps in one block = "<<steps<<endl<<endl;
  cout<<"Step lenght = "<<L<<endl;
  cout<<"Initial position = "<<x0<<endl;
  cout<<"Num punti per riempire histo = "<<histofill_blk*nblk<<endl;

  ReadInput.close();

}

double min(double s, double t)
{
    if (s <= t) return s;
    else return t;
}

double error(double av, double av2, int n)
{
   if(n==0){
      return 0;
   }
   else{
      return sqrt((av2-av*av)/n);
   }
}
