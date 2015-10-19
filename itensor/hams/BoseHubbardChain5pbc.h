//
// Bose-Hubbard Hamiltonian with up to 5 bosons per site
// Periodic boundaries introduced via a long bond
//
// Andrew Goldsborough - 13/05/2014
// written for ITensor V1.0
//
#ifndef __ITENSOR_HAMS_BOSEHUBBARDCHAINPBC_H
#define __ITENSOR_HAMS_BOSEHUBBARDCHAINPBC_H
#include "../mpo.h"
#include "../model/bosehubbard5.h"

class BoseHubbardChain5pbc
    {
    public:

    BoseHubbardChain5pbc(const BoseHubbard5& model,
                 const OptSet& opts = Global::opts());

    Real
    t() const { return t_; }
    void
    t(Real val) { initted_ = false; t_ = val; }

    Real
    U() const { return U_; }
    void
    U(Real val) { initted_ = false; U_ = val; }

    Real
    mu() const { return mu_; }
    void
    mu(Real val) { initted_ = false; mu_ = val; }
    
    operator MPO() { init_(); return H; }

    operator IQMPO() { init_(); return H; }

    private:

    ///////////////////
    //
    // Data Members

    const BoseHubbard5& model_;
    Real t_,U_,mu_;
    bool initted_;
    MPO H;

    //
    //////////////////

    void 
    init_();

    }; //class BoseHubbardChain

inline BoseHubbardChain5pbc::
BoseHubbardChain5pbc(const BoseHubbard5& model, 
             const OptSet& opts)
    : 
    model_(model), 
    initted_(false)
    { 
    t_ = opts.getReal("t",1);
    U_ = opts.getReal("U",0);
    mu_ = opts.getReal("mu",0);
    }

void inline BoseHubbardChain5pbc::
init_()
    {
    if(initted_) return;

    H = MPO(model_);

    const int Ns = model_.N();
    const int k = 6;

    std::vector<Index> links(Ns+1);
    for(int l = 0; l <= Ns; ++l) links.at(l) = Index(nameint("hl",l),k);

    ITensor W;

    for(int n = 1; n <= Ns; ++n)
        {
        
        if(n==1)
            {
            ITensor& W = H.Anc(n);
            Index &row = links[n-1], &col = links[n];

            W = ITensor(model_.si(n),model_.siP(n),row,col);

            //Identity strings
            W += model_.op("Id",n) * row(1) * col(1);
            
            //PBC term
            W += model_.op("Bdag",n) * row(1) * col(4);
            W += model_.op("B",n) * row(1) * col(5);

            //Hopping terms -t*(b^d_i b_{i+1} + b_i b^d_{i+1})
            W += model_.op("Bdag",n) * row(1) * col(2) * (-0.5*t_);
            W += model_.op("B",n) * row(1) * col(3) * (-0.5*t_);

    	    //on-site terms U/2 * n_i(n_i-1) + mu * n_i
            W += model_.op("N",n) * row(1) * col(k) * (mu_);
            W += model_.op("OS",n) * row(1) * col(k) * (0.5*U_);//OS = N(N-1)
            }
 
        else if(n==Ns)
            {
            ITensor& W = H.Anc(n);
            Index &row = links[n-1], &col = links[n];

            W = ITensor(model_.si(n),model_.siP(n),row,col);

            //Identity strings
            W += model_.op("Id",n) * row(k) * col(k);
            
            //PBC term
            W += model_.op("B",n) * row(4) * col(k) * (-0.5*t_);
            W += model_.op("Bdag",n) * row(5) * col(k) * (-0.5*t_);

            //Hopping terms -t*(b^d_i b_{i+1} + b_i b^d_{i+1})
            W += model_.op("B",n) * row(2) * col(k);
            W += model_.op("Bdag",n) * row(3) * col(k);

    	    //on-site terms U/2 * n_i(n_i-1) + mu * n_i
            W += model_.op("N",n) * row(1) * col(k) * (mu_);
            W += model_.op("OS",n) * row(1) * col(k) * (0.5*U_);//OS = N(N-1)
            }

        else
            {
            ITensor& W = H.Anc(n);
            Index &row = links[n-1], &col = links[n];

            W = ITensor(model_.si(n),model_.siP(n),row,col);

            //Identity strings
            W += model_.op("Id",n) * row(1) * col(1);
            W += model_.op("Id",n) * row(k) * col(k);
            
            //PBC identities
            W += model_.op("Id",n) * row(4) * col(4);
            W += model_.op("Id",n) * row(5) * col(5);

            //Hopping terms -t*(b^d_i b_{i+1} + b_i b^d_{i+1})
            W += model_.op("Bdag",n) * row(1) * col(2) * (-0.5*t_);
            W += model_.op("B",n) * row(1) * col(3) * (-0.5*t_);
            W += model_.op("B",n) * row(2) * col(k);
            W += model_.op("Bdag",n) * row(3) * col(k);

    	    //on-site terms U/2 * n_i(n_i-1) + mu * n_i
            W += model_.op("N",n) * row(1) * col(k) * (mu_);
            W += model_.op("OS",n) * row(1) * col(k) * (0.5*U_);//OS = N(N-1)
            } 
       }

    H.Anc(1) *= ITensor(links.at(0)(1));
    H.Anc(Ns) *= ITensor(links.at(Ns)(k));
    
    initted_ = true;
    }

#endif
