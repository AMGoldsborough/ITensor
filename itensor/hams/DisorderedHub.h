//
// Hubbard Model Hamiltonian  
// with vector parameters to enable disorder
//
// Andrew Goldsborough - 12/02/2014
// written for ITensor V1.0
//
#ifndef __ITENSOR_HAMS_DISORDEREDHUB_H
#define __ITENSOR_HAMS_DISORDEREDHUB_H
#include "../mpo.h"
#include "../model/hubbard.h"

class DisorderedHub
    {
    public:

    DisorderedHub(const Hubbard& model,
                  const Vector& t,
                  const Vector& U,
                  const Vector& mu);

    operator MPO() { init_(); return H; }

    operator IQMPO() { init_(); return H; }

    private:

    ///////////////////
    //
    // Data Members

    const Hubbard& model_;
    const Vector& t_;
    const Vector& U_;
    const Vector& mu_;
    bool initted_;
    MPO H;

    //
    //////////////////

    void 
    init_();

    }; //class DisorderedHub

inline DisorderedHub::
DisorderedHub(const Hubbard& model, 
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

void inline DisorderedHub::
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
        ITensor& W = H.Anc(n);
        Index &row = links[n-1], &col = links[n];

        W = ITensor(model_.si(n),model_.siP(n),row,col);

        //Identity strings
        W += model_.op("Id",n) * row(1) * col(1);
        W += model_.op("Id",n) * row(k) * col(k);

        //Kinetic energy/hopping terms, defined as -t_*(c^d_i c_{i+1} + h.c.)
        W += model_.op("F*Cup",n) * row(k) * col(2) * t_(n);
        W += model_.op("F*Cdn",n) * row(k) * col(3) * t_(n);
        W += model_.op("Cdagup*F",n) * row(k) * col(4) * t_(n);
        W += model_.op("Cdagdn*F",n) * row(k) * col(5) * t_(n);
        W += model_.op("Cdagup",n) * row(2) * col(1) * (-1.0);
        W += model_.op("Cdagdn",n) * row(3) * col(1) * (-1.0);
        W += model_.op("Cup",n) * row(4) * col(1) * (-1.0);
        W += model_.op("Cdn",n) * row(5) * col(1) * (-1.0);
        
        // on-site U * Nup*Ndn + mu * Nup + mu * Ndn
        W += model_.op("Nupdn",n) * row(k) * col(1) * U_(n);
        W += model_.op("Nup",n) * row(k) * col(1) * mu_(n);
        W += model_.op("Ndn",n) * row(k) * col(1) * mu_(n);
        }

    H.Anc(1) *= ITensor(links.at(0)(k));
    H.Anc(Ns) *= ITensor(links.at(Ns)(1));

    initted_ = true;
    }

#endif
