//
// Full disordered DMRG using DisorderedHeisHalf
// outputting entanglement and correlation
//
// Andrew Goldsborough 28/01/2015
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

//        cout << h(i) << endl;
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
    cout << format("\nGround State Energy = %.10f")%En << endl;

    cout << "\nTotal QN of Ground State = " << totalQN(psi) << "\n";

    //open file to write to
    ofstream energyfile;
    std::stringstream fnamestreamenergy;
    fnamestreamenergy << "./energy/" << N << "_" << Jstr << "_" << Jdis << "_" << Jzstr << "_" << Jzdis << "_" << hstr << "_" << hdis << "_" << seed << "_" << chimax << "_energy_spinhalf.txt";
    
    energyfile.open(fnamestreamenergy.str().c_str());

    energyfile << format("%.15e\n")%En; //print to file

    energyfile.close();
    
    
    //
    // Calculate spin correlation functions
    //
    cout << "Printing correlation functions" << endl;
    double Spcorr = 0.0;
    
    MPO SpcorrMPO = MPO(model); 
    
    //open file to write to
    ofstream corrfile;
    std::stringstream fnamestream;
    fnamestream << "./spcorr/" << N << "_" << Jstr << "_" << Jdis << "_" << Jzstr << "_" << Jzdis << "_" << hstr << "_" << hdis << "_" << seed << "_" << chimax << "_spcorr_spinhalf.txt";
    corrfile.open(fnamestream.str().c_str());
    
    for(int j=1; j<=N; j++) 
        {

        //move psi to get ready to measure at position j
        psi.position(j);

        for(int i=j+1; i<=N; i++)
            {
            // Sz.Sz
            //make a bond ket/bra based on wavefunction 
            IQMPS Spket = psi;
            Spket.position(j);
            
            //attach Sz at j and i
            Spket.Anc(j) = Spket.Anc(j) * model.op("Sz",j);
            Spket.Anc(i) = Spket.Anc(i) * model.op("Sz",i);

            //remove primes
            Spket.Anc(j).noprime();
            Spket.Anc(i).noprime();

            //contract making use of normalization
            IQTensor L;
            if(j==1)
            {
                L = Spket.A(1) * conj(primed(psi.A(1),psi.LinkInd(1)));
            }
            else
            {
                L = Spket.A(j) * conj(primed(psi.A(j),psi.LinkInd(j)));
            }
            for(int k = j+1; k < i; ++k) 
                { 
                L = L * Spket.A(k) * conj(primed(psi.A(k),Link)); 
                }
            L = L * Spket.A(i);

            Complex z = BraKet(primed(psi.A(i),psi.LinkInd(i-1)),L);
            Spcorr = z.real();
	    
	    //Add 0.5*Sp.Sm
	    Spket = psi;
            Spket.position(j);
            
            //attach Sp at j and Sm at i
            Spket.Anc(j) = Spket.Anc(j) * model.op("Sp",j);
            Spket.Anc(i) = Spket.Anc(i) * model.op("Sm",i);

            //remove primes
            Spket.Anc(j).noprime();
            Spket.Anc(i).noprime();

            //contract making use of normalization
            if(j==1)
            {
                L = Spket.A(1) * conj(primed(psi.A(1),psi.LinkInd(1)));
            }
            else
            {
                L = Spket.A(j) * conj(primed(psi.A(j),psi.LinkInd(j)));
            }
            for(int k = j+1; k < i; ++k) 
                { 
                L = L * Spket.A(k) * conj(primed(psi.A(k),Link)); 
                }
            L = L * Spket.A(i);

            z = BraKet(primed(psi.A(i),psi.LinkInd(i-1)),L);
            Spcorr += 0.5 * z.real();
	    
	    //Add 0.5*Sm.Sp
	    Spket = psi;
            Spket.position(j);
            
            //attach Sm at j and Sp at i
            Spket.Anc(j) = Spket.Anc(j) * model.op("Sm",j);
            Spket.Anc(i) = Spket.Anc(i) * model.op("Sp",i);

            //remove primes
            Spket.Anc(j).noprime();
            Spket.Anc(i).noprime();

            //contract making use of normalization
            //IQTensor L;
            if(j==1)
            {
                L = Spket.A(1) * conj(primed(psi.A(1),psi.LinkInd(1)));
            }
            else
            {
                L = Spket.A(j) * conj(primed(psi.A(j),psi.LinkInd(j)));
            }
            for(int k = j+1; k < i; ++k) 
                { 
                L = L * Spket.A(k) * conj(primed(psi.A(k),Link)); 
                }
            L = L * Spket.A(i);

            z = BraKet(primed(psi.A(i),psi.LinkInd(i-1)),L);
            Spcorr += 0.5 * z.real();

	    //print to file
            corrfile << j << " " << i << " " << format("%.15e\n")%Spcorr; 

            }
        }
        
    //close file
    corrfile.close();
    
    //
    // calculate ee
    //
    cout << "Printing entanglement entropy" << endl;
    
    //open file to write to
    ofstream eefile;
    std::stringstream fnamestreamee;
    fnamestreamee << "./ee/" << N << "_" << Jstr << "_" << Jdis << "_" << Jzstr << "_" << Jzdis << "_" << hstr << "_" << hdis << "_" << seed << "_" << chimax << "_ee_spinhalf.txt";
    eefile.open(fnamestreamee.str().c_str());

    ofstream eespfile;
    std::stringstream fnamestreameesp;
    fnamestreameesp << "./eespec/" << N << "_" << Jstr << "_" << Jdis << "_" << Jzstr << "_" << Jzdis << "_" << hstr << "_" << hdis << "_" << seed << "_" << chimax << "_eespec_spinhalf.txt";
    eefile.open(fnamestreameesp.str().c_str());
    
    for(int i = 1; i<= N-1; ++i)
    {
        Spectrum spec = psi.spectrum(i);

        Vector D = spec.eigsKept();
        
        double ee = 0.;
	
	eespfile << i << " ";

        for(int j = 1; j <= D.Length(); ++j)
        {
            ee -= D(j)*log2(fabs(D(j)));
            //ee = ee - D(i)*D(i) * log2(D(i)*D(i));
	    eespfile << format("%.15e\n")%D(j) << " ";
        }  
        
        eespfile << endl;
    
        eefile << i << " " << format("%.15e\n")%ee << endl;
    }
    
    //close files
    eefile.close();
    eespfile.close();
    
    //check execution time
    texe = clock() - tstart;

    cout << "Elapsed time is " << ((float)texe)/CLOCKS_PER_SEC << " seconds" << endl;

    return 0;
    }
