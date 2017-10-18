//
// Disordered DMRG using DisorderedHeisHalfpbc
// outputting entanglement and correlation
// for use with Js from dMERA
//
// Andrew Goldsborough 29/12/2016
//
#include "core.h"
#include "mpo.h"
#include "model/spinhalf.h"
#include "hams/DisorderedHeisHalfpbc.h"
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
    if(argc != 9)
    {
        cout << "Requires 8 inputs L, Jstr, Jdis, Jzstr, chimax, Pdist, seed, sweepmax" << endl;
        exit(EXIT_FAILURE);
    }
    
    int N = atoi(argv[1]);              //system size
    double Jstr = atof(argv[2]);        //coupling strength
    double Jdis = atof(argv[3]);        //coupling disorder
    double Jzstr = atof(argv[4]);       //anisotropy strength
    int chimax = atoi(argv[5]);         //maximum bond dimension
    int Pdist = atoi(argv[6]);          //Probability distribution
    int seed = atoi(argv[7]);           //seed for random number generator  
    int sweepmax = atoi(argv[8]);       //max number of sweeps 
    clock_t tstart = clock();           //start time
    clock_t texe;                       //execution time

    cout << "L = " << N << ", Jstr = " << Jstr << ", Jdis = " << Jdis << ", Jzstr = " << Jzstr << ", chimax = " << chimax << ", Pdist = " << Pdist << ", seed = " << seed << ", sweepmax = " << sweepmax << endl;

    //
    // Get interaction strengths from file
    //
    Vector J(N);
    Vector Jz(N);
    Vector h(N);
   
    //open file to read from 
    ifstream Jfile;
    std::stringstream fnamestreamJ;
    fnamestreamJ << "./J/" << N << "_" << Jstr << "_" << Jdis << "_" << Pdist << "_" << seed << "_J.txt";
    Jfile.open(fnamestreamJ.str().c_str());

    if ( Jfile.is_open() )
    {
        for(int i=1; i<=N; i++)
        {
            Jfile >> J(i);
        }
        Jfile.close();
    }
    else 
    {
        cout << "Unable to open " << fnamestreamJ.str() << endl;
        exit (EXIT_FAILURE);
    }

    for(int i = 1; i<= N; i++)
    {
        Jz(i) = Jzstr;
    }
    
    for(int i = 1; i<= N; i++)
    {
        h(i) = 0.0;
    }

    //
    // Initialize the site degrees of freedom.
    //
    Spin model(N); 

    //
    // Create the Hamiltonian matrix product operator (MPO)
    //
    cout << "creating MPO" << endl;
    MPO H = DisorderedHeisHalfpbc(model,
                                  J,
                                  Jz,
                                  h);

    //Random starting point
    MPS psi(model); 

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
    Sweeps sweeps(sweepmax);
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
    
    //open file to write to
    ofstream energyfile;
    std::stringstream fnamestreamenergy;
    fnamestreamenergy << "./energy/" << N << "_" << Jstr << "_" << Jdis << "_" << Jzstr << "_" << chimax << "_" << Pdist << "_" << seed << "_" << sweepmax << "_energy_DMRG_dMERA.txt";
    
    energyfile.open(fnamestreamenergy.str().c_str());

    energyfile << format("%.15e\n")%En; //print to file

    energyfile.close();
    
    //
    // Calculate spin correlation functions
    //
    cout << "Printing correlation functions" << endl;
    double Spcorr = 0.0;
    
    //open file to write Sz corr
    ofstream szcorrfile;
    std::stringstream fnamestreamsz;
    fnamestreamsz << "./szcorr/" << N << "_" << Jstr << "_" << Jdis << "_" << Jzstr << "_" << chimax << "_" << Pdist << "_" << seed << "_" << sweepmax << "_szcorr_DMRG_dMERA.txt";
    szcorrfile.open(fnamestreamsz.str().c_str());
    
    //open file to write SpSm corr
    ofstream spsmcorrfile;
    std::stringstream fnamestreamspsm;
    fnamestreamspsm << "./spsmcorr/" << N << "_" << Jstr << "_" << Jdis << "_" << Jzstr << "_" << chimax << "_" << Pdist << "_" << seed << "_" << sweepmax << "_spsmcorr_DMRG_dMERA.txt";
    spsmcorrfile.open(fnamestreamspsm.str().c_str());
    
    //open file to write SmSp corr
    ofstream smspcorrfile;
    std::stringstream fnamestreamsmsp;
    fnamestreamsmsp << "./smspcorr/" << N << "_" << Jstr << "_" << Jdis << "_" << Jzstr << "_" << chimax << "_" << Pdist << "_" << seed << "_" << sweepmax << "_smspcorr_DMRG_dMERA.txt";
    smspcorrfile.open(fnamestreamsmsp.str().c_str());
    
    for(int j=1; j<=N; j++) 
    {

        //move psi to get ready to measure at position j
        psi.position(j);

        for(int i=j+1; i<=N; i++)
        {
            // Sz.Sz
            //make a bond ket/bra based on wavefunction 
            MPS Spket = psi;
            Spket.position(j);
            
            //attach Sz at j and i
            Spket.Anc(j) = Spket.Anc(j) * model.op("Sz",j);
            Spket.Anc(i) = Spket.Anc(i) * model.op("Sz",i);

            //remove primes
            Spket.Anc(j).noprime();
            Spket.Anc(i).noprime();

            //contract making use of normalization
            ITensor L;
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
	    
	        //output Sz corr
            szcorrfile << j << " " << i << " " << format("%.15e\n")%Spcorr; 
            
            //Sp.Sm
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
            Spcorr = z.real();
	    
	        //output SpSm corr
            spsmcorrfile << j << " " << i << " " << format("%.15e\n")%Spcorr; 
            
	        //Sm.Sp
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
            Spcorr = z.real();

	        //print to file
            smspcorrfile << j << " " << i << " " << format("%.15e\n")%Spcorr; 
        }
    }
        
    //close files
    szcorrfile.close();
    spsmcorrfile.close();
    smspcorrfile.close();
    
    //
    // Calculate spin expectation values
    //
    cout << "Printing expectation values" << endl;
    double Spexp = 0.0;
    
    //open file to write to
    ofstream szexpfile;
    std::stringstream fnamestreamszexp;
    fnamestreamszexp << "./szexp/" << N << "_" << Jstr << "_" << Jdis << "_" << Jzstr << "_" << chimax << "_" << Pdist << "_" << seed << "_" << sweepmax << "_szexp_DMRG_dMERA.txt";
    szexpfile.open(fnamestreamszexp.str().c_str());
    
    //open file to write spexp
    ofstream spexpfile;
    std::stringstream fnamestreamspexp;
    fnamestreamspexp << "./spexp/" << N << "_" << Jstr << "_" << Jdis << "_" << Jzstr << "_" << chimax << "_" << Pdist << "_" << seed << "_" << sweepmax << "_spexp_DMRG_dMERA.txt";
    spexpfile.open(fnamestreamspexp.str().c_str());
    
    //open file to write smexp
    ofstream smexpfile;
    std::stringstream fnamestreamsmexp;
    fnamestreamsmexp << "./smexp/" << N << "_" << Jstr << "_" << Jdis << "_" << Jzstr << "_" << chimax << "_" << Pdist << "_" << seed << "_" << sweepmax << "_smexp_DMRG_dMERA.txt";
    smexpfile.open(fnamestreamsmexp.str().c_str());
    
    for(int j=1; j<=N; j++) 
    {

        //move psi to get ready to measure at position j
        psi.position(j);

        // Sz exp
        ITensor ketz = psi.A(j)*model.op("Sz",j);
        ITensor bra = conj(psi.A(j));
        bra.prime(Site);
        Spexp = Dot(bra,ketz);
        
	    //output Sz exp
        szexpfile << j << " " << format("%.15e\n")%Spexp; 
        
        // Sp exp
        ITensor ketp = psi.A(j)*model.op("Sp",j);
        Spexp = Dot(bra,ketp);
        
	    //output Sp exp
        spexpfile << j << " " << format("%.15e\n")%Spexp; 
        
        // Sm exp
        ITensor ketm = psi.A(j)*model.op("Sm",j);
        Spexp = Dot(bra,ketm);
        
	    //print to file
        smexpfile << j << " " << format("%.15e\n")%Spexp; 
    }
        
    //close files
    szexpfile.close();
    spexpfile.close();
    smexpfile.close();

    //
    //calculate ee for blocks
    //

    cout << "Printing block entanglement entropy" << endl;
    
    //open file to write to
    ofstream beefile;
    std::stringstream fnamestreambee;
    fnamestreambee << "./ee/" << N << "_" << Jstr << "_" << Jdis << "_" << Jzstr << "_" << chimax << "_" << Pdist << "_" << seed << "_" << sweepmax << "_ee_DMRG_dMERA.txt";
    beefile.open(fnamestreambee.str().c_str());

    Matrix Lm;                      //matrix for ee calculation
    Vector nrg;                     //vector of eigenvalues
    Matrix Veig;                    //matrix of eigenvectors

    for(int i=2; i<=N-1; i++) 
    {
        //move psi to get ready to measure at position j
        psi.position(i);

        for(int j=i; j<=N-1; j++)
        {

            //contract transfer matrix
            ITensor L;
            L = psi.A(i) * conj(primed(psi.A(i),Link));

            for(int k = i+1; k <= j; ++k) 
            { 
                L = L * psi.A(k) * conj(primed(psi.A(k),Link)); 
            }
            //combine to make a matrix
            L.toMatrix22(psi.LinkInd(i-1),psi.LinkInd(j),primed(psi.LinkInd(i-1)),primed(psi.LinkInd(j)),Lm);

            //find eigenvalues and vectors
            Lm = 0.5*(Lm + Lm.t());
            EigenValues(Lm,nrg,Veig);
            
            double ee = 0.;
            for(int k = 1; k <= nrg.Length(); ++k)
            {
                ee -= nrg(k)*log2(fabs(nrg(k)));
            }
            
            beefile << i-1 << " " << j << " " <<  format("%.15e")%ee << endl;
        }

        Spectrum spec = psi.spectrum(i-1);
        Vector D = spec.eigsKept();
        double ee = 0.;

        for(int k = 1; k <= D.Length(); ++k)
        {
            ee -= D(k)*log2(fabs(D(k)));
        }  
    
        beefile << i-1 << " " << N << " " << format("%.15e")%ee << endl;
    }

    Spectrum spec = psi.spectrum(N-1);
    Vector D = spec.eigsKept();
    double ee = 0.;

    for(int k = 1; k <= D.Length(); ++k)
    {
        ee -= D(k)*log2(fabs(D(k)));
    }  

    beefile << N-1 << " " << N << " " << format("%.15e")%ee << endl;

    //close files
    beefile.close();
    
    //check execution time
    texe = clock() - tstart;

    cout << "Elapsed time is " << ((float)texe)/CLOCKS_PER_SEC << " seconds" << endl;

    return 0;
}

