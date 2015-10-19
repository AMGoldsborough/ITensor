// Disordered Bose-Hubbard model DMRG using ITensor
// 5 bosons per site (maximum) and periodic boundaries
// Filling n/L = 1
// Andrew Goldsborough 16/05/2014
//
#include "core.h"
#include "mpo.h"
#include "model/bosehubbard5.h"
#include "hams/DisorderedBH5pbc.h"
#include "boost/random.hpp"
#include <fstream>
#include <sstream>
using boost::format;
using boost::mt19937;
using boost::uniform_01;
using namespace std;

//create random number generator function
double random01(mt19937 generator)
    {
    static uniform_01<mt19937> dist(generator);
    return dist();
    }

int main(int argc, char* argv[])
{
    //
    // Get system variables
    //
    if(argc != 10)
        {
        cout << "Requires 9 inputs L, tstr, tdis, Ustr, Udis, mustr, mudis, seed, chimax" << endl;
        exit(EXIT_FAILURE);
        }
    
    int N = atoi(argv[1]);              //system size
    double tstr = atof(argv[2]);        //hopping strength
    double tdis = atof(argv[3]);        //hopping disorder
    double Ustr = atof(argv[4]);        //On-site interaction
    double Udis = atof(argv[5]);        //On-site disorder
    double mustr = atof(argv[6]);       //chemical potential
    double mudis = atof(argv[7]);       //disorder strength on mu
    int seed = atoi(argv[8]);           //seed for random number generator  
    int chimax = atoi(argv[9]);         //maximum bond dimension
    clock_t tstart = clock();           //start time
    clock_t texe;                       //execution time

    //divide by 100 so that bash can input integer variables
    tstr = tstr/100;
    tdis = tdis/100;
    Ustr = Ustr/100;
    Udis = Udis/100;
    mustr = mustr/100;
    mudis = mudis/100;

    cout << "L = " << N << ", tstr = " << tstr << ", tdis = " << tdis << ", Ustr = " << Ustr << ", Udis = " << Udis << ", mustr = " << mustr << ", mudis = " << mudis << ", chimax = " << chimax << endl;
    
    //
    // Generate interaction strengths
    //
    Vector t(N);
    Vector U(N);
    Vector mu(N);

    for(int i = 1; i<= N; i++)
    {
      mt19937 generator(seed);
      
      //disorder 8 (box)
      t(i) = tstr + tdis*(random01(generator)-0.5);

//      cout << t(i) << endl;
    }
    
    for(int i = 1; i<= N; i++)
    {
      mt19937 generator(seed);
      
      //disorder 8 (box)
      U(i) = Ustr + Udis*(random01(generator)-0.5);

//      cout << U(i) << endl;
    }
    
    for(int i = 1; i<= N; i++)
    {
      mt19937 generator(seed);
      
      //disorder 8 (box) 
      mu(i) = mustr + mudis*(random01(generator)-0.5);

//      cout << mu(i) << endl;
    }

    //
    // Initialize the site degrees of freedom.
    //
    cout << "Initializing sites" << endl;
    int Npart = N; //number of particles, default is N (half filling)
    
    BoseHubbard5 model(N);
    
    //
    // Create the Hamiltonian matrix product operator.
    // Here we use the IQMPO class which is an MPO of 
    // IQTensors, tensors whose indices are sorted
    // with respect to quantum numbers
    //
    cout << "creating MPO" << endl;
    IQMPO H = DisorderedBH5pbc(model,
                               t, 
                               U,
                               mu);
  
    //
    // Set the initial wavefunction matrix product state
    // to be a Neel state.
    //
    cout << "occupying MPS " << endl;
    InitState initState(model);
    int p = Npart;
    for(int i = N; i >= 1; --i) 
        {
        if(p > i)
            {
            cout << "Doubly occupying site " << i << endl;
            initState.set(i,&BoseHubbard5::Dou);
            p -= 2;
            }
        else
        if(p > 0)
            {
            cout << "Singly occupying site " << i << endl;
            initState.set(i, &BoseHubbard5::Occ);
            p -= 1;
            }
        else
            {
            initState.set(i,&BoseHubbard5::Em);
            }
        }
    
    IQMPS psi(model,initState);
    
    cout << totalQN(psi) << endl;
    
    //
    // Set the parameters controlling the accuracy of the DMRG
    // calculation for each DMRG sweep. 
    // Here less than 5 cutoff values are provided, for example,
    // so all remaining sweeps will use the last one given (= 1E-10).
    //
    Sweeps sweeps(20);
//    sweeps.maxm() = 10,20,50,100,200;
    sweeps.maxm() = 10,20,chimax;
    sweeps.cutoff() = 1E-10,1E-12,1E-15,1E-20;
    sweeps.niter() = 2;
    sweeps.noise() = 1E-7,1E-8,1E-8,0.0;
    cout << sweeps;
    
    Real En = dmrg(psi,H,sweeps,Quiet());
    
    //
    // Print the final energy reported by DMRG
    //
    cout << format("\nGround State Energy = %.10f\n")%En;
    cout << format("\nUsing psiHphi = %.10f\n") % psiHphi(psi,H,psi);
    
    cout << "\nTotal QN of Ground State = " << totalQN(psi) << "\n";
   
    //
    // calculate ee
    //
    cout << "Printing entanglement entropy and spectrum" << endl;
    
    //open files to write to
    ofstream eefile;
    std::stringstream fnamestreamee;
    fnamestreamee << "./ee/" << N << "_" << tstr << "_" << tdis << "_" << Ustr << "_" << Udis << "_" << mustr << "_" << mudis << "_" << seed << "_" << chimax << "_ee_5_pbc.txt";
    eefile.open(fnamestreamee.str().c_str());

    ofstream eespfile;
    std::stringstream fnamestreameesp;
    fnamestreameesp << "./eespec/" << N << "_" << tstr << "_" << tdis << "_" << Ustr << "_" << Udis << "_" << mustr << "_" << mudis << "_" << seed << "_" << chimax << "_eespec_5_pbc.txt";
    eespfile.open(fnamestreameesp.str().c_str());

    for(int i = 1; i<= N-1; ++i)
    {
        Spectrum spec = psi.spectrum(i);

        Vector D = spec.eigsKept();
        
        double ee = 0.;

        eespfile << i << " ";

        for(int j = 1; j <= D.Length(); ++j)
        {
            ee -= D(j)*log(fabs(D(j)));
            //ee = ee - D(i)*D(i) * log2(D(i)*D(i));
            eespfile << D(j) << " ";
        }  
    
        eespfile << endl;

        eefile << i << " " << ee << endl;
    }

    //close files
    eefile.close();
    eespfile.close();

    //check execution time
    texe = clock() - tstart;

    cout << "Elapsed time is " << ((float)texe)/CLOCKS_PER_SEC << " seconds" << endl;

    return 0;
    
    }                           
 
