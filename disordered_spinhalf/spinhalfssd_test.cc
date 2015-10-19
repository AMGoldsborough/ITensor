//
// Testing the sine squared deformation see Hikihara & Nishino
// http://dx.doi.org/10.1103/PhysRevB.83.060414
//
// Andrew Goldsborough 13/05/2014
//
#include "core.h"
#include "mpo.h"
#include "model/spinhalf.h"
#include "hams/DisorderedHeisHalf.h"
#include "boost/random.hpp"
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

typedef SpinHalf
Spin;           //use S=1/2 degrees of freedom

int 
main(int argc, char* argv[])
    {

    //
    // Get system variables
    //
    if(argc != 10)
        {
        cout << "Requires 9 inputs L, Jstr, Jdis, Jzstr, Jzdis, hstr, hdis, seed, chimax" << endl;
        exit(EXIT_FAILURE);
        }
    
    int N = atoi(argv[1]);              //system size
    double Jstr = atof(argv[2]);        //coupling strength
    double Jdis = atof(argv[3]);        //coupling disorder
    double Jzstr = atof(argv[4]);       //anisotropy strength
    double Jzdis = atof(argv[5]);       //anisotropy disorder
    double hstr = atof(argv[6]);        //magnetic field strength
    double hdis = atof(argv[7]);        //magnetic field disorder
    int seed = atoi(argv[8]);           //seed for random number generator  
    int chimax = atoi(argv[9]);         //maximum bond dimension
    clock_t tstart = clock();           //start time
    clock_t texe;                       //execution time
    const long double PI = 3.141592653589793238462643383279502884L;
    double ssdnorm = 0.;                    //ssd normalisation

    cout << "L = " << N << ", Jstr = " << Jstr << ", Jdis = " << Jdis << ", Jzstr = " << Jzstr << ", Jzdis = " << Jzdis << ", hstr = " << hstr << ", hdis = " << hdis << ", seed = " << seed << ", chimax = " << chimax << endl;

    //
    // Generate interaction strengths
    //
    Vector J(N);
    Vector Jz(N);
    Vector h(N);

    for(int i = 1; i<= N; i++)
        {
        mt19937 generator(seed);
          
        //disorder 8 (box)
        J(i) = Jstr + Jdis*(random01(generator)-0.5);

//        cout << t(i) << endl;
        }
    
    for(int i = 1; i<= N; i++)
        {
        mt19937 generator(seed);
      
        //disorder 8 (box)
        Jz(i) = Jzstr + Jzdis*(random01(generator)-0.5);

//        cout << U(i) << endl;
        }
    
    for(int i = 1; i<= N; i++)
        {
        mt19937 generator(seed);
      
        //disorder 8 (box) 
        h(i) = hstr + hdis*(random01(generator)-0.5);

//        cout << mu(i) << endl;
        }

    //introduce sine squared deformation
    for(int i=1; i<=N; i++)
    {
        J(i) = J(i) * sin(PI*i/N)*sin(PI*i/N);
        h(i) = h(i) * sin(PI*(i-0.5)/N)*sin(PI*(i-0.5)/N);

        //calculate normalization
        ssdnorm = ssdnorm + sin(PI*i/N)*sin(PI*i/N);
    }


    //
    // Initialize the site degrees of freedom.
    //
    Spin model(N);    // make a chain of N spins

    //
    // Create the Hamiltonian matrix product operator (MPO)
    //
    cout << "creating MPO" << endl;
    IQMPO H = DisorderedHeisHalf(model,
                                  J,
                                  Jz,
                                  h);

    //
    // Set the initial wavefunction matrix product state (MPS)
    // to be a Neel state.
    //
//    int Npart = N/2; //half-fill

    cout << "occupying MPS " << endl;
    InitState initState(model);
//    int p = Npart;
    for(int i = 1; i <= N; ++i) 
        {
        if(i%2 == 1)
            initState.set(i,"Up");
        else
            initState.set(i,"Dn");
        }

    IQMPS psi(model,initState);

    cout << totalQN(psi) << endl;

    //
    // psiHphi calculates matrix elements of MPO's with respect to MPS's
    // psiHphi(psi,H,psi) = <psi|H|psi>
    //
    cout << format("Initial energy = %.5f") % psiHphi(psi,H,psi) << endl;

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
    cout << format("\nGround State Energy = %.10f")%(En*N/ssdnorm) << endl;

    cout << "\nTotal QN of Ground State = " << totalQN(psi) << "\n";

    //
    //MEASURING SPIN
    //

    //vector of z-components of spin for each site
    Vector Sz(N);
    Sz = 0.0;
    
    //open file to write to
    ofstream szfile;
    std::stringstream fnamestreamsz;
    fnamestreamsz << "./szcorr/" << N << "_" << Jstr << "_" << Jdis << "_" << Jzstr << "_" << seed << "_" << chimax << "_szcorr_ssd.txt";
    szfile.open(fnamestreamsz.str().c_str());

    //Sj dot Sj+1 by means of Splus and Sminus
    double SdotS; 
    SdotS = 0.0;

    //file to write SdotS information to
    ofstream spfile;
    std::stringstream fnamestreamsp;
    fnamestreamsp << "./spcorr/" << N << "_" << Jstr << "_" << Jdis << "_" << Jzstr << "_" << seed << "_" << chimax << "_spcorr_ssd.txt";
    spfile.open(fnamestreamsp.str().c_str());
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

        //operator for sz at site j
        ITensor szj_op = model.op("Sz",j);

        //take an inner product 
        Sz(j) = (bra*szj_op*ket).toReal();

        szfile << std::scientific << Sz(j) << endl; //print to file

        for(int i=j+1; i<=N; i++)
            { 
            //make a bond ket/bra based on wavefunction 
            IQMPS ssket = psi;
            
            //start with z components
            ssket.Anc(j) = ssket.Anc(j) * model.op("Sz",j);
            ssket.Anc(i) = ssket.Anc(i) * model.op("Sz",i);

            //remove primes
            ssket.Anc(j).noprime();
            ssket.Anc(i).noprime();

            //contract
            SdotS = psiphi(psi,ssket);

            //add in S+ and S- components:
            ssket = psi;
            ssket.Anc(j) = ssket.Anc(j) * 0.5 * model.op("Sp",j);
            ssket.Anc(i) = ssket.Anc(i) * model.op("Sm",i);
            ssket.Anc(j).noprime();
            ssket.Anc(i).noprime();
            SdotS += psiphi(psi,ssket);

            ssket = psi;
            ssket.Anc(j) = ssket.Anc(j) * 0.5 * model.op("Sm",j);
            ssket.Anc(i) = ssket.Anc(i) * model.op("Sp",i);
            ssket.Anc(j).noprime();
            ssket.Anc(i).noprime();
            SdotS += psiphi(psi,ssket);

            spfile << j << " " << i << " " << std::scientific << SdotS << endl; //print to file
            }
        }


    //close files
    szfile.close();
    spfile.close();
    
    //
    // calculate ee
    //
    cout << "Printing entanglement entropy" << endl;
    
    //open file to write to
    ofstream eefile;
    std::stringstream fnamestreamee;
    fnamestreamee << "./ee/" << N << "_" << Jstr << "_" << Jdis << "_" << seed << "_" << chimax << "_ee_ssd.txt";
    eefile.open(fnamestreamee.str().c_str());

    for(int i = 1; i<= N-1; ++i)
    {
        Spectrum spec = psi.spectrum(i);

        Vector D = spec.eigsKept();
        
        double ee = 0.;

        for(int j = 1; j <= D.Length(); ++j)
        {
            ee -= D(j)*log2(fabs(D(j)));
        }  
    
        eefile << i << " " << ee << endl;
    }

    //close file
    eefile.close();

    //check execution time
    texe = clock() - tstart;

    cout << "Elapsed time is " << ((float)texe)/CLOCKS_PER_SEC << " seconds" << endl;

    return 0;
    }
