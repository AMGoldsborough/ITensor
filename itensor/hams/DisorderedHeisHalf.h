//
// Spin Half Heisenberg Hamiltonian  
// with vector parameters to enable disorder
//
// Andrew Goldsborough - 18/03/2014
// written for ITensor V1.0
//
#ifndef __ITENSOR_HAMS_DISORDEREDHEISHALF_H
#define __ITENSOR_HAMS_DISORDEREDHEISHALF_H
#include "../mpo.h"
#include "../model/spinhalf.h"

class DisorderedHeisHalf
    {
    public:

    DisorderedHeisHalf(const SpinHalf& model,
                       const Vector& J,
                       const Vector& Jz,
                       const Vector& h);

    operator MPO() { init_(); return H; }

    operator IQMPO() { init_(); return H; }

    private:

    ///////////////////
    //
    // Data Members

    const SpinHalf& model_;
    const Vector& J_;
    const Vector& Jz_;
    const Vector& h_;
    bool initted_;
    MPO H;

    //
    //////////////////

    void 
    init_();

    }; //class DisorderedHeisHalf

inline DisorderedHeisHalf::
DisorderedHeisHalf(const SpinHalf& model, 
                   const Vector& J,
                   const Vector& Jz,
                   const Vector& h)
    : 
    model_(model), 
    J_(J),
    Jz_(Jz),
    h_(h),
    initted_(false)
    { 
    Vector J_(model_.N(),1.0);
    Vector Jz_(model_.N(),1.0);
    Vector h_(model_.N(),0.0);
    }

void inline DisorderedHeisHalf::
init_()
    {
    if(initted_) return;

    H = MPO(model_);

    const int Ns = model_.N();
    const int k = 5;

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

        //Spin.Spin terms (note anisotropy term on Sz.Sz)
        W += model_.op("Sp",n) * row(1) * col(2) * (0.5*J_(n));
        W += model_.op("Sm",n) * row(1) * col(3) * (0.5*J_(n));
        W += model_.op("Sz",n) * row(1) * col(4) * J_(n) * Jz_(n);
        W += model_.op("Sm",n) * row(2) * col(k);
        W += model_.op("Sp",n) * row(3) * col(k);
        W += model_.op("Sz",n) * row(4) * col(k);
        
        // on-site (magnetic field)
        W += model_.op("Sz",n) * row(1) * col(k) * h_(n);
        }

    H.Anc(1) *= ITensor(links.at(0)(1));
    H.Anc(Ns) *= ITensor(links.at(Ns)(k));

    initted_ = true;
    }

#endif
