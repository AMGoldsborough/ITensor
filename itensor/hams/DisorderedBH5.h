//
// Bose-Hubbard Hamiltonian with up to 5 bosons per site 
// with vector parameters to enable disorder
//
// Andrew Goldsborough - 06/03/2014
// written for ITensor V1.0
//
#ifndef __ITENSOR_HAMS_DISORDEREDBH_H
#define __ITENSOR_HAMS_DISORDEREDBH_H
#include "../mpo.h"
#include "../model/bosehubbard5.h"

class DisorderedBH5
    {
    public:

    DisorderedBH5(const BoseHubbard5& model,
                 const Vector& t,
                 const Vector& U,
                 const Vector& mu);
    
    operator MPO() { init_(); return H; }

    operator IQMPO() { init_(); return H; }

    private:

    ///////////////////
    //
    // Data Members

    const BoseHubbard5& model_;
    const Vector& t_;
    const Vector& U_;
    const Vector& mu_;
    bool initted_;
    MPO H;

    //
    //////////////////

    void 
    init_();

    }; //class BoseHubbardChain

inline DisorderedBH5::
DisorderedBH5(const BoseHubbard5& model,
             const Vector& t,
             const Vector& U,
             const Vector& mu)
    : 
    model_(model), 
    t_(t),
    U_(U),
    mu_(mu),
    initted_(false)
    { 
    Vector t_(model_.N(),1.0);
    Vector U_(model_.N(),0.0);
    Vector mu_(model_.N(),0.0);
    }

void inline DisorderedBH5::
init_()
    {
    if(initted_) return;

    H = MPO(model_);

    const int Ns = model_.N();
    const int k = 4;

    std::vector<Index> links(Ns+1);
    for(int l = 0; l <= Ns; ++l) links.at(l) = Index(nameint("hl",l),k);

    ITensor W;
    for(int n = 1; n <= Ns; ++n)
        {
        ITensor& W = H.Anc(n);
        Index &row = links[n-1], &col = links[n];

        W = ITensor(model_.si(n),model_.siP(n),row,col);

        //Identity strings
        W += model_.op("Id",n) * row(1) * col(1);
        W += model_.op("Id",n) * row(k) * col(k);

        //Hopping terms -t/2 * (b^d_i b_{i+1} + b_i b^d_{i+1})
        W += model_.op("Bdag",n) * row(1) * col(2) * (-0.5*t_(n));
        W += model_.op("B",n) * row(1) * col(3) * (-0.5*t_(n));
        W += model_.op("B",n) * row(2) * col(k);
        W += model_.op("Bdag",n) * row(3) * col(k);

    	//on-site terms U/2 * n_i(n_i-1) + mu * n_i
        W += model_.op("N",n) * row(1) * col(k) * (mu_(n));
        W += model_.op("OS",n) * row(1) * col(k) * (0.5*U_(n));//OS = N(N-1)
        }

    H.Anc(1) *= ITensor(links.at(0)(1));
    H.Anc(Ns) *= ITensor(links.at(Ns)(k));
    
    initted_ = true;
    }

#endif
