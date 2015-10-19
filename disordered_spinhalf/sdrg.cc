//
// recreating my SDRG matlab code in ITensor
// Andrew Goldsborough 20/01/2013
//
// note that C/++ indexes vectors from zero, the indexing scheme is:
//
// Spins:      0    1    2  
// links/J:  0 | 1  | 2  | 3
// MPO:      -|0|--|1|--|2|-- ...
//             |    |    |
//
///////////////////////////////////////////////////////////////////

#include "core.h"
#include "model/spinhalf.h"
#include "model/spinone.h"
#include "hams/Heisenberg.h"
#include "mpo.h"
#include "combiner.h"
#include "boost/random.hpp"
#include <fstream>
#include <sstream>
using boost::format;
using boost::mt19937;
using boost::uniform_01;
using boost::array;
using namespace std;

//use S=1/2 degrees of freedom
typedef SpinHalf
Spin;           

//create random number generator function
double random01(mt19937 generator)
{
      static uniform_01<mt19937> dist(generator);
      return dist();
}

int 
main(int argc, char* argv[])
{
    //
    // Get system variables
    //
    if(argc != 6)
    {
        cout << "Requires 5 inputs L, Jstr, Jdis, Jseed, chimax" << endl;
        exit(EXIT_FAILURE);
    }

    int N = atoi(argv[1]);                //system size
    const double Jstr = atof(argv[2]);    //Interaction strength
    const double Jdis = atof(argv[3]);    //Disorder on interaction
    const int Jseed = atoi(argv[4]);      //seed for random number generator
    const int chimax = atoi(argv[5]);     //maximum bond dimension
    const int Norig = N;                  //original system size
 
    cout << "N =  " << N << " Jstr = " << Jstr << " Jdis = " << Jdis << " Jseed = " << Jseed << " chimax = " << chimax << endl;

    //
    // Initialize the site degrees of freedom.
    //
    Spin model(N);   

    //
    // Make some definitions for later
    //
    const int k = 5;                //MPO bond dimension for Heisenberg Hamiltonian
    std::vector<Index> links(N+1);  //index vector for links
    std::vector<Index> spins(N+1);  //index vector for spins
    std::vector<Index> spinsorig(N+1);  //index vector for spins for use later
    std::vector<ITensor*> MPOvec;   //vector of ITensors as an MPO
    std::vector<ITensor*> isovec;   //vector of isometries
    std::vector<Index> w1(N);       //index vector for iso top leg
    std::vector<Index> w2(N);       //index vector for iso bottom left
    std::vector<Index> w3(N);       //index vector for iso bottom right
    clock_t tstart = clock();       //start time
    clock_t texe;                   //execution time
    std::vector<double> J(N+1);     //vector of interaction strengths
    double Jmaxval = 0.;            //current maximum J
    int Jmax = 0;                   //index of max J
    std::vector<int> Jorder(N-1);   //order for TTN creation
    ITensor block;                  //block for renormalization
    Matrix blockM;                  //block made into a matrix
    Vector nrg;                     //vector of eigenvalues
    Matrix Veig;                    //matrix of eigenvectors
    double diff=0.;                 //difference between eigenvalues
    int diffloc=0;                  //element where diff occurs
    ITensor blockfull;              //ITensor for block + couplings
    ITensor blockcomb;              //combined two-site MPO tensor
    int iter;                       //counts the iteration number
    std::vector<int> tL(N-1,-1);    //tensor below left
    std::vector<int> tR(N-1,-1);    //tensor below right
    std::vector<int> tU(N-1,-1);    //tensor up
    std::vector<int> curt(N,-1);    //current tensor to find tL and tR
    std::vector<int> nofsites(N,1); //number if sites represented by each MPO tensor
    int nocount;                    //counts the number of sites to the left
    std::vector<int> groundL(N,-1); //need to know which tensors connect to the ground
    std::vector<int> groundR(N,-1); //need to know which tensors connect to the ground
    std::vector<int> ground(N,-1);  //need to know which tensors connect to the ground

    //
    // Generate interaction strengths
    //
    for(int i = 1; i< N; i++)
    {
      mt19937 generator(Jseed);
      
      //disorder 7
      J[i] = Jstr + Jstr*Jdis*(random01(generator)-0.5);

      cout << J[i] << endl;
    }

    cout << endl;

    //
    // Create the Hamiltonian matrix product operator (MPO)
    //
    //MPO built as a vector of pointers 
    for(int l = 0; l <= N; ++l) links.at(l) = Index(nameint("hl",l),k);
    for(int l = 0; l <= N-1; ++l) spins.at(l) = Index(nameint("spin",l),2,Site);

    spinsorig = spins;

    for(int i = 1; i <= N; i++)
     {
        //Build MPO tensors
        Index &row = links[i-1], &col = links[i];
        Index &spin = spins[i-1];
        ITensor* W = new ITensor(spin,primed(spin),row,col);
        (*W)(spin(1),primed(spin(1)),row(1),col(1)) = 1.;
        (*W)(spin(2),primed(spin(2)),row(1),col(1)) = 1.;
        (*W)(spin(1),primed(spin(2)),row(1),col(2)) = J[i] * 0.5;
        (*W)(spin(2),primed(spin(1)),row(1),col(3)) = J[i] * 0.5;
        (*W)(spin(1),primed(spin(1)),row(1),col(4)) = J[i] * 0.5;
        (*W)(spin(2),primed(spin(2)),row(1),col(4)) = J[i] * -0.5;
        (*W)(spin(2),primed(spin(1)),row(2),col(5)) = 1.;
        (*W)(spin(1),primed(spin(2)),row(3),col(5)) = 1.;
        (*W)(spin(1),primed(spin(1)),row(4),col(5)) = 0.5;
        (*W)(spin(2),primed(spin(2)),row(4),col(5)) = -0.5;
        (*W)(spin(1),primed(spin(1)),row(5),col(5)) = 1.;
        (*W)(spin(2),primed(spin(2)),row(5),col(5)) = 1.;
        MPOvec.push_back(W);
     }

    //
    // Start algorithm
    //
    for(iter=0; iter <= Norig-3; ++iter)
//    for(int iter=1; iter <= 1; ++iter)
    {
        //find max J
        Jmaxval = 0.;
        Jmax = 0;
    
        for(int i=1;i<N;i++)
        {

            if(J[i]>Jmaxval)
            {
                Jmaxval=J[i];
                Jmax = i;
            }
        }

        //save Jmax in Jorder for future calculations
        Jorder[iter] = Jmax;

        //save tL and tR
        tL[iter] = curt[Jmax-1];
        tR[iter] = curt[Jmax];

        //make block of the two most coupled MPO tensors remember C++ counts from zero, hence [Jmax-1]
        block = (ITensor(links.at(Jmax-1)(1)) * (*MPOvec[Jmax-1])) * ((*MPOvec[Jmax]) * ITensor(links.at(Jmax+1)(k)));

        //combine to make a matrix
        block.toMatrix22(spins.at(Jmax-1),spins.at(Jmax),primed(spins.at(Jmax-1)),primed(spins.at(Jmax)),blockM);
    
        //find eigenvalues and vectors
        blockM = 0.5*(blockM + blockM.t());
        EigenValues(blockM,nrg,Veig);
   
        if(chimax < nrg.Length())
        {
            //take only full blocks
            
            for (int i = chimax; i >= 1; --i)
            {
                diff = nrg(i+1) - nrg(i);
                diffloc = i;
                
                if(diff >= 1e-10)
                    break;
                else
                    continue;
            }
            
            Matrix Veigtrunc(nrg.Length(),diffloc);
    
            for(int i=1; i<=nrg.Length(); ++i)
            {
                for(int j=1; j<=diffloc; ++j)
                {
                    Veigtrunc(i,j) = Veig(i,j);
                }
            }
    
            Veig = Veigtrunc;
        }
    
        //make this an isometry
        w1.at(iter) = Index(nameint("Snew",iter),Veig.Ncols(),Site);
        w2.at(iter) = spins.at(Jmax-1);
        w3.at(iter) = spins.at(Jmax);
        ITensor* wtemp = new ITensor;
        wtemp->fromMatrix12(w3.at(iter),w2.at(iter),w1.at(iter),Veig);
        isovec.push_back(wtemp);

        //contract isometry to perform renormalization
        *MPOvec[Jmax-1] = ((conj(*isovec[iter]) * (*MPOvec[Jmax-1])) * (*MPOvec[Jmax])) * primed(*isovec[iter]);

        //update the spin indices
        spins.at(Jmax-1) = w1.at(iter);
 
        //update curt
        curt[Jmax-1] = iter;

        //store ground
        if(tL[iter] == -1)
        {
            nocount = 0;

            for(int i = 0; i <= Jmax-1; ++i)
                nocount = nocount + nofsites[i];

            groundL[nocount-1] = iter;
        }
        else
        {
            tU[tL[iter]] = iter;
        }

        if(tR[iter] == -1)
        {
            nocount = 0;

            for(int i = 0; i <= Jmax; ++i)
                nocount = nocount + nofsites[i];

            groundR[nocount-1] = iter;
        }
        else
        {
            tU[tR[iter]] = iter;
        }

        //update nofsites
        nofsites[Jmax-1] = nofsites[Jmax-1]+nofsites[Jmax];

        //remove the extra site
        J.erase(J.begin()+Jmax);
        spins.erase(spins.begin()+Jmax);
        links.erase(links.begin()+Jmax);
        delete *(MPOvec.begin()+Jmax);
        MPOvec.erase(MPOvec.begin()+Jmax); 
        curt.erase(curt.begin()+Jmax);
        nofsites.erase(nofsites.begin()+Jmax);

        //find renormalised couplings
        if(Jmax != 1)
        {
            //find left
            block = (ITensor(links.at(Jmax-2)(1)) * (*MPOvec[Jmax-2])) * ((*MPOvec[Jmax-1]) * ITensor(links.at(Jmax)(k)));

            //make into a matrix
            block.toMatrix22(spins.at(Jmax-2),spins.at(Jmax-1),primed(spins.at(Jmax-2)),primed(spins.at(Jmax-1)),blockM);

            //find eigenvalues and vectors
            blockM = 0.5*(blockM + blockM.t());
            EigenValues(blockM,nrg,Veig);
        
            diff=0.;
            diffloc=0; //location of heigest gap
     
            if(chimax < nrg.Length())
            {
                //find the heighest gap
                for (int i = chimax; i >= 1; --i)
                {
                    diff = nrg(i+1) - nrg(i);
                    
                    if(diff >= 1e-10)
                    {
                        break;
                    }
                    else
                        continue;
                }
            }
            else
            {
                //find the heighest gap
                for (int i = nrg.Length()-1; i >= 1; --i)
                {
                    diff = nrg(i+1) - nrg(i);
                   
                    if(diff >= 1e-10)
                        break;
                    else
                        continue;
                }
            }

            J[Jmax-1] = diff;
        }
        
        if(Jmax != N-1)
        {
            //find right
            block = (ITensor(links.at(Jmax-1)(1)) * (*MPOvec[Jmax-1])) * ((*MPOvec[Jmax]) * ITensor(links.at(Jmax+1)(k)));
            
            //make into a matrix
            block.toMatrix22(spins.at(Jmax-1),spins.at(Jmax),primed(spins.at(Jmax-1)),primed(spins.at(Jmax)),blockM);
        
            //find eigenvalues and vectors
            blockM = 0.5*(blockM + blockM.t());
            EigenValues(blockM,nrg,Veig);
        
            diff=0.;
            diffloc=0; //location of heigest gap
     
            if(chimax < nrg.Length())
            {
                //find the heighest gap
                for (int i = chimax; i >= 1; --i)
                {
                    diff = nrg(i+1) - nrg(i);
                    
                    if(diff >= 1e-10)
                    {
                        break;
                    }
                    else
                        continue;
                }
            }
            else
            {
                //find the heighest gap
                for (int i = nrg.Length()-1; i >= 1; --i)
                {
                    diff = nrg(i+1) - nrg(i);
                    
                    if(diff >= 1e-10)
                        break;
                    else
                        continue;
                }
            }
   
            J[Jmax] = diff;
        }
        
        N = N-1;
    }

    //final step
    iter = Norig-2;
    Jmax = 1;
    Jorder[iter] = Jmax;

    //save tL and tR
    tL[iter] = curt[Jmax-1];
    tR[iter] = curt[Jmax];

    //make block of the two most coupled MPO tensors
    block = (ITensor(links.at(Jmax-1)(1)) * (*MPOvec[Jmax-1])) * ((*MPOvec[Jmax]) * ITensor(links.at(Jmax+1)(k)));

    //combine to make a matrix
    block.toMatrix22(spins.at(Jmax-1),spins.at(Jmax),primed(spins.at(Jmax-1)),primed(spins.at(Jmax)),blockM);

    //find eigenvalues and vectors
    blockM = 0.5*(blockM + blockM.t());
    EigenValues(blockM,nrg,Veig);

    Matrix Veigtrunc(nrg.Length(),1);
    
    for(int i=1; i<=nrg.Length(); ++i)
    {
        Veigtrunc(i,1) = Veig(i,1);
    }

    Veig = Veigtrunc;
    
    //make this an isometry
    w1.at(iter) = Index(nameint("Snew",iter),Veig.Ncols(),Site);
    w2.at(iter) = spins.at(Jmax-1);
    w3.at(iter) = spins.at(Jmax);
    ITensor* wtemp = new ITensor;
    wtemp->fromMatrix12(w3.at(iter),w2.at(iter),w1.at(iter),Veig);
    isovec.push_back(wtemp);

    cout << nrg(1) << endl;

    //store ground
    if(tL[iter] == -1)
    {
        nocount = 0;

        for(int i = 0; i <= Jmax-1; ++i)
            nocount = nocount + nofsites[i];

        groundL[nocount-1] = iter;
    }
    else
    {
        tU[tL[iter]] = iter;
    }

    if(tR[iter] == -1)
    {
        nocount = 0;

        for(int i = 0; i <= Jmax; ++i)
            nocount = nocount + nofsites[i];

        groundR[nocount-1] = iter;
    }
    else
    {
        tU[tR[iter]] = iter;
    }

    //combine grounds
    for (int i=0; i<=Norig-1; ++i)
        ground[i] = groundL[i] + groundR[i] + 1;

    //update nofsites
    nofsites[Jmax-1] = nofsites[Jmax-1]+nofsites[Jmax];

//    for (int i=0; i<=Norig-2; ++i)
//        cout << tL[i] << " " ;
//
//    cout << endl;
//
//    for (int i=0; i<=Norig-2; ++i)
//        cout << tR[i] << " " ;
//
//    cout << endl;
//
//    for (int i=0; i<=Norig-2; ++i)
//        cout << tU[i] << " " ;
//
//    cout << endl;
//    
//    for (int i=0; i<=Norig-1; ++i)
//        cout << groundL[i] << " " ;
//    
//    cout << endl;
//
//    for (int i=0; i<=Norig-1; ++i)
//        cout << groundR[i] << " " ;
//
//    cout << endl;
//
//    for (int i=0; i<=Norig-1; ++i)
//        cout << ground[i] << " " ;
//
//    cout << endl;

//    for(int i=0; i<=iter; ++i)
//            cout << Jorder[i] << endl;

    //check execution time
    texe = clock() - tstart;
    
    //
    // Calculate spin correlation functions
    //
    cout << "Printing correlation functions" << endl;

    //variable for correlation function
    int stepprev = 0;                       //previous step of left corr contraction
    int numL = 0;                           //number of tensors in left corr contraction
    int numR = 0;                           //number of tensors in right corr contraction
    Index corrlink;                         //link for spin corr
    ITensor corr1;                          //left spin corr operator
    ITensor corr2;                          //right spin corr operator
    int nextt;                              //next tensor for path to root
    int path_element;                       //element in path vector
    std::vector<int> path1;                 //path to root for left site
    std::vector<int> path2;                 //path to root for right site
    std::vector<int> path_inter(Norig-1);   //intersection of paths
    std::vector<int>::iterator inter_it;    //iterator for path_inter

    //open file to write to
    ofstream corrfile;
    std::stringstream fnamestream;
    fnamestream << "./spcorr/" << Norig << "_" << Jstr << "_" << Jdis << "_" << Jseed << "_" << chimax << "_spcorrSDRG.txt";
    corrfile.open(fnamestream.str().c_str());

    for(int site1=0; site1<=Norig-2; ++site1)
    {
        //
        // Create the Hamiltonian matrix product operator (MPO)
        //
    
        //Build MPO tensor for left spin
        corrlink = Index(nameint("hl",site1),3);
        Index &corrspin1 = spinsorig[site1];
        corr1 = ITensor(corrspin1,primed(corrspin1),corrlink);
        corr1(corrspin1(1),primed(corrspin1(2)),corrlink(1)) = 0.5;
        corr1(corrspin1(2),primed(corrspin1(1)),corrlink(2)) = 0.5;
        corr1(corrspin1(1),primed(corrspin1(1)),corrlink(3)) = 0.5;
        corr1(corrspin1(2),primed(corrspin1(2)),corrlink(3)) = -0.5;

        //initialize stepprev
        stepprev = 0;

        //find path to the root node for site 1
        nextt = ground[site1];

        path1.push_back(nextt);

        while (nextt < Norig-2)
        {
            nextt = tU[nextt];

            path1.push_back(nextt);    
        }

        for(int site2 = site1+1; site2<=Norig-1; ++site2)
        {

        //Build MPO tensor for right spin
        Index &corrspin2 = spinsorig[site2];
        corr2 = ITensor(corrspin2,primed(corrspin2),corrlink);
        corr2(corrspin2(2),primed(corrspin2(1)),corrlink(1)) = 1;
        corr2(corrspin2(1),primed(corrspin2(2)),corrlink(2)) = 1;
        corr2(corrspin2(1),primed(corrspin2(1)),corrlink(3)) = 0.5;
        corr2(corrspin2(2),primed(corrspin2(2)),corrlink(3)) = -0.5;

        //find path to the root node for site 1
        nextt = ground[site2];

        path2.push_back(nextt);

        while (nextt < Norig-2)
        {
            nextt = tU[nextt];

            path2.push_back(nextt);    
        }
 
        //find the intersection between the paths
        inter_it = std::set_intersection(path1.begin(),path1.end(),path2.begin(),path2.end(),path_inter.begin());

        //contract left operator up to the common node
        path_element = stepprev;
        nextt = path1[path_element];
        
        while (nextt < path_inter[0])
        {
            corr1 = conj(*isovec[nextt]) * corr1; 

            corr1.noprime();
            corr1.prime(w2.at(nextt));
            corr1.prime(w3.at(nextt));
                
            corr1 = corr1 * primed(*isovec[nextt]);

            path_element = path_element + 1;

            nextt = path1[path_element];

            numL = numL + 1;
        }
        
        //store how far the left contraction got to save time
        stepprev = path_element;

        //contract right operator up to the common node
        path_element = 0;
        nextt = path2[path_element];

        while (nextt < path_inter[0])
        {
            corr2 = conj(*isovec[nextt]) * corr2;

            corr2.noprime();
            corr2.prime(w2.at(nextt));
            corr2.prime(w3.at(nextt));

            corr2 = corr2 * primed(*isovec[nextt]);

            path_element = path_element + 1;

            numR = numR + 1;

            nextt = path2[path_element];
        }

        //contract the common node
        corr2 = (((corr1 * conj(*isovec[nextt])) * corr2) * primed(*isovec[nextt]));

        path_element = path_element + 1;

        nextt = path2[path_element];

        //contract up to the root node
        for(int i = path_element; i<path2.size(); ++i)
        {
            corr2 = conj(*isovec[nextt]) * corr2;

            corr2.noprime();
            corr2.prime(w2.at(nextt));
            corr2.prime(w3.at(nextt));

            corr2 = corr2 * primed(*isovec[nextt]);

            path_element = path_element + 1;

            nextt = path2[path_element];
        }

        //print to file
        corrfile << site1 + 1 << " " << site2 + 1 << " " << corr2.toReal() << " " << numL + numR + 1 << endl;

        //clean up for next iteration
        path2.clear();
        numR = 0;

        for(int i = 0; i<Norig-1; ++i)
            path_inter[i] = 0;

        }

        path1.clear();
        numL = 0;


     }


    
    
   
//    //close file
    corrfile.close();
 
    cout << "Elapsed time is " << ((float)texe)/CLOCKS_PER_SEC << " seconds" << endl;

    
    return 0;
}

