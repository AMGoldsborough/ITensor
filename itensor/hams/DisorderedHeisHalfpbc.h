//
// Spin Half Heisenberg Hamiltonian  
// with vector parameters to enable disorder
// Periodic boundaries introduced via a long bond
//
// Andrew Goldsborough - 13/05/2014
// written for ITensor V1.0
//
#ifndef __ITENSOR_HAMS_DISORDEREDHEISHALFPBC_H
#define __ITENSOR_HAMS_DISORDEREDHEISHALFPBC_H
#include "../mpo.h"
#include "../model/spinhalf.h"

class DisorderedHeisHalfpbc
    {
    public:

    DisorderedHeisHalfpbc(const SpinHalf& model,
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

    }; //class DisorderedHeisHalfpbc

inline DisorderedHeisHalfpbc::
DisorderedHeisHalfpbc(const SpinHalf& model, 
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

void inline DisorderedHeisHalfpbc::
init_()
    {
    if(initted_) return;

    H = MPO(model_);

    const int Ns = model_.N();
    const int k = 8;

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

            //PBC terms
            W += model_.op("Sp",n) * row(1) * col(5);
            W += model_.op("Sm",n) * row(1) * col(6);
            W += model_.op("Sz",n) * row(1) * col(7);

            //Spin.Spin terms (note anisotropy term on Sz.Sz)
            W += model_.op("Sp",n) * row(1) * col(2) * (0.5*J_(n));
            W += model_.op("Sm",n) * row(1) * col(3) * (0.5*J_(n));
            W += model_.op("Sz",n) * row(1) * col(4) * J_(n) * Jz_(n);
            
            // on-site (magnetic field)
            W += model_.op("Sz",n) * row(1) * col(k) * h_(n);
            }
        else if(n==Ns)
            {
            ITensor& W = H.Anc(n);
            Index &row = links[n-1], &col = links[n];

            W = ITensor(model_.si(n),model_.siP(n),row,col);

            //Identity strings
            W += model_.op("Id",n) * row(k) * col(k);

            //PBC identities
            W += model_.op("Sm",n) * row(5) * col(k) * (0.5*J_(n));
            W += model_.op("Sp",n) * row(6) * col(k) * (0.5*J_(n));
            W += model_.op("Sz",n) * row(7) * col(k) * J_(n) * Jz_(n);
 
            //Spin.Spin terms (note anisotropy term on Sz.Sz)
            W += model_.op("Sm",n) * row(2) * col(k);
            W += model_.op("Sp",n) * row(3) * col(k);
            W += model_.op("Sz",n) * row(4) * col(k);
            
            // on-site (magnetic field)
            W += model_.op("Sz",n) * row(1) * col(k) * h_(n);
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
            W += model_.op("Id",n) * row(5) * col(5);
            W += model_.op("Id",n) * row(6) * col(6);
            W += model_.op("Id",n) * row(7) * col(7);

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
        }

    H.Anc(1) *= ITensor(links.at(0)(1));
    H.Anc(Ns) *= ITensor(links.at(Ns)(k));

    initted_ = true;
    }

#endif
