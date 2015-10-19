// Disordered Bose-Hubbard model DMRG using ITensor
// 5 bosons per site (maximum)
// Filling n/L = 1
// Andrew Goldsborough 07/02/2014
//
#include "core.h"
#include "mpo.h"
#include "model/bosehubbard5.h"
#include "hams/DisorderedBH5.h"
#include "boost/random.hpp"
#include <fstream>
#include <sstream>
using boost::format;
using boost::mt19937;
using boost::uniform_01;
using boost::array;
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

    cout << "L = " << N << ", tstr = " << tstr << ", tdis = " << tdis << ", Ustr = " << Ustr << ", Udis = " << Udis << ", mustr = " << mustr << ", mudis = " << mudis << ", seed = " << seed << ", chimax = " << chimax << endl;
    
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
    IQMPO H = DisorderedBH5(model,
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
    // calculate n_i expectation value, b_i^d b_j, dn_i dn_j correlation functions and string order parameter
    //

    cout << "Printing correlation functions (n, b^d b, dn dn, string order) " << endl;

    //vector of n expectation values
    double nexpect;
    nexpect = 0.;
    
    //open file to write to
    ofstream nfile;
    std::stringstream fnamestreamn;
    fnamestreamn << "./nexpect/" << N << "_" << tstr << "_" << tdis << "_" << Ustr << "_" << Udis << "_" << mustr << "_" << mudis << "_" << seed << "_" << chimax << "_nexpect_5.txt";

    nfile.open(fnamestreamn.str().c_str());

    //now for the bb^d correlations
    double bdbcorr; 
    bdbcorr = 0.0;

    //file to write bdbcorr information to
    ofstream bdbfile;
    std::stringstream fnamestreambdb;
    fnamestreambdb << "./bdbcorr/" << N << "_" << tstr << "_" << tdis << "_" << Ustr << "_" << Udis << "_" << mustr << "_" << mudis << "_" << seed << "_" << chimax << "_bdbcorr_5.txt";
    bdbfile.open(fnamestreambdb.str().c_str());
    
    //now for the nn correlations
    double nncorr; 
    nncorr = 0.0;

    //file to write nncorr information to
    ofstream nnfile;
    std::stringstream fnamestreamnn;
    fnamestreamnn << "./nncorr/" << N << "_" << tstr << "_" << tdis << "_" << Ustr << "_" << Udis << "_" << mustr << "_" << mudis << "_" << seed << "_" << chimax << "_nncorr_5.txt";
    nnfile.open(fnamestreamnn.str().c_str());

    //now for the sring order
    double strorder; 
    strorder = 0.0;
    IQMPO strop(model);

    //make operators
    for(int k = 1; k < N; k++)
    {
        //by Taylor expansion
        // exp(x) = 1 + x +  x^2/2! + x^3/3! ..
        // = 1 + x * (1 + x/2 *(1 + x/3 * (...
        // ~ ((x/3 + 1) * x/2 + 1) * x + 1
        IQTensor expdn = model.op("Id",k);
        IQTensor dn = Complex_i * Pi * (model.op("N",k) - (Npart/N)*model.op("Id",k));
        IQTensor term = dn * Complex_i * Pi;
    
        dn.mapprime(1,2);
        dn.mapprime(0,1);
    
        for(int ord = 100; ord >= 1; --ord)
        {
            term /= ord;
            expdn = model.op("Id",k) + term;
            term = expdn * dn;
            term.mapprime(2,1);
        }

        strop.Anc(k) = expdn;
    }

    //file to write strorder information to
    ofstream sofile;
    std::stringstream fnamestreamso;
    fnamestreamso << "./strorder/" << N << "_" << tstr << "_" << tdis << "_" << Ustr << "_" << Udis << "_" << mustr << "_" << mudis << "_" << seed << "_" << chimax << "_strorder_5.txt";
    sofile.open(fnamestreamso.str().c_str());


    cout.precision(15);

    for(int j=1; j<=N; j++) 
        {
        //move psi to get ready to measure at position j
        psi.position(j);

        //after calling psi.position(j), psi.A(j) returns a 
        //local representation of the wavefunction,
        //which is the proper way to make measurements
        //i.e. take expectation values with local operators.

        //Dirac "ket" for wavefunction
        ITensor ket = psi.A(j);

        //Dirac "bra" for wavefunction
        ITensor bra = conj(primed(ket,Site));

        //operator for n at site j
        ITensor nj_op = model.op("N",j);

        //take an inner product 
        nexpect = (bra*nj_op*ket).toReal();

        nfile << std::scientific << nexpect << endl; //print to file

        for(int i=j+1; i<=N; i++)
            {
            // b^d b correlation function
            //make a bond ket/bra based on wavefunction 
            IQMPS bdbket = psi;
            bdbket.position(j);
            
            //attach b^d and b
            bdbket.Anc(j) = bdbket.Anc(j) * model.op("Bdag",j);
            bdbket.Anc(i) = bdbket.Anc(i) * model.op("B",i);

            //remove primes
            bdbket.Anc(j).noprime();
            bdbket.Anc(i).noprime();

            //contract
//            bdbcorr = psiphi(psi,bdbket);

            //contract making use of normalization
            IQTensor L;
            if(j==1)
            {
                L = bdbket.A(1) * conj(primed(psi.A(1),psi.LinkInd(1)));
            }
            else
            {
                L = bdbket.A(j) * conj(primed(psi.A(j),psi.LinkInd(j)));
            }
            for(int k = j+1; k < i; ++k) 
                { 
                L = L * bdbket.A(k) * conj(primed(psi.A(k),Link)); 
                }
            L = L * bdbket.A(i);

            Complex z = BraKet(primed(psi.A(i),psi.LinkInd(i-1)),L);
            bdbcorr = z.real();

            bdbfile << j << " " << i << " " << std::scientific << bdbcorr << endl; //print to file

            //
            // dn correlation function
            //

            //make a bond ket/bra based on wavefunction 
            IQMPS nnket = psi;
            nnket.position(j);
            
            //attach n's
            nnket.Anc(j) = nnket.Anc(j) * (model.op("N",j) - (Npart/N)*model.op("Id",j));
            nnket.Anc(i) = nnket.Anc(i) * (model.op("N",i) - (Npart/N)*model.op("Id",i));

            //remove primes
            nnket.Anc(j).noprime();
            nnket.Anc(i).noprime();

            //contract
//            nncorr = psiphi(psi,nnket);

            //contract making use of normalization
//            IQTensor L;
            if(j==1)
            {
                L = nnket.A(1) * conj(primed(psi.A(1),psi.LinkInd(1)));
            }
            else
            {
                L = nnket.A(j) * conj(primed(psi.A(j),psi.LinkInd(j)));
            }
            for(int k = j+1; k < i; ++k) 
                { 
                L = L * nnket.A(k) * conj(primed(psi.A(k),Link)); 
                }
            L = L * nnket.A(i);

            z = BraKet(primed(psi.A(i),psi.LinkInd(i-1)),L);
            nncorr = z.real();

            nnfile << j << " " << i << " " << std::scientific << nncorr << endl; //print to file

            //
            // string order
            //

            //make a bond ket/bra based on wavefunction 
            IQMPS soket = psi;
            soket.position(j);
            
            //attach n's
            soket.Anc(j) = soket.Anc(j) * (model.op("N",j) - (Npart/N)*model.op("Id",j));
            soket.Anc(i) = soket.Anc(i) * (model.op("N",i) - (Npart/N)*model.op("Id",i));

            //remove primes
            soket.Anc(j).noprime();
            soket.Anc(i).noprime();

            //attach exp(i*pi*dn)
            for(int k = j; k < i; k++)
            {
                soket.Anc(k) = soket.Anc(k) * strop.A(k);
            
                soket.Anc(k).noprime();
            }
            
            //contract making use of normalization
//            IQTensor L;
            if(j==1)
            {
                L = soket.A(1) * conj(primed(psi.A(1),psi.LinkInd(1)));
            }
            else
            {
                L = soket.A(j) * conj(primed(psi.A(j),psi.LinkInd(j)));
            }
            for(int k = j+1; k < i; ++k) 
                { 
                L = L * soket.A(k) * conj(primed(psi.A(k),Link)); 
                }
            L = L * soket.A(i);

            z = BraKet(primed(psi.A(i),psi.LinkInd(i-1)),L);
            strorder = z.real();

            sofile << j << " " << i << " " << std::scientific << strorder << endl; //print to file
            }
        }
        


    //close files
    nfile.close();
    bdbfile.close();
    nnfile.close();
    sofile.close();
    

    //
    // calculate ee
    //
    cout << "Printing entanglement entropy and spectrum" << endl;
    
    //open files to write to
    ofstream eefile;
    std::stringstream fnamestreamee;
    fnamestreamee << "./ee/" << N << "_" << tstr << "_" << tdis << "_" << Ustr << "_" << Udis << "_" << mustr << "_" << mudis << "_" << seed << "_" << chimax << "_ee_5.txt";
    eefile.open(fnamestreamee.str().c_str());

    ofstream eespfile;
    std::stringstream fnamestreameesp;
    fnamestreameesp << "./eespec/" << N << "_" << tstr << "_" << tdis << "_" << Ustr << "_" << Udis << "_" << mustr << "_" << mudis << "_" << seed << "_" << chimax << "_eespec_5.txt";
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
 
