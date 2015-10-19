//
// Distributed under the ITensor Library License, Version 1.0.
//    (See accompanying LICENSE file.)
//
// Bose Hubbard model with maximum 5 bosons per site
#ifndef __ITENSOR_HUBBARD_H
#define __ITENSOR_HUBBARD_H
#include "../model.h"

class BoseHubbard5 : public Model
    {
    public:

    BoseHubbard5();

    BoseHubbard5(int N, 
            const OptSet& opts = Global::opts());

    bool
    conserveNf() const { return conserveNf_; }

    IQIndexVal
    Em(int i) const;

    IQIndexVal
    Occ(int i) const;

    IQIndexVal
    Dou(int i) const;

    IQIndexVal
    Tri(int i) const;

    IQIndexVal
    Qua(int i) const;

    IQIndexVal
    Qui(int i) const;

    IQIndexVal
    EmP(int i) const;

    IQIndexVal
    OccP(int i) const;

    IQIndexVal
    DouP(int i) const;

    IQIndexVal
    TriP(int i) const;

    IQIndexVal
    QuaP(int i) const;

    IQIndexVal
    QuiP(int i) const;

    private:

    virtual int
    getN() const;

    virtual const IQIndex&
    getSi(int i) const;

    virtual IQIndexVal
    getState(int i, const String& state) const;

    virtual IQTensor
    getOp(int i, const String& opname, const OptSet& opts = Global::opts()) const;

    virtual void
    doRead(std::istream& s);

    virtual void
    doWrite(std::ostream& s) const;

    virtual void
    constructSites();

    //IQTensor
    IQTensor
    makeN(int i) const;
    IQTensor
    makeOS(int i) const;
    IQTensor
    makeB(int i) const;
    IQTensor
    makeBdag(int i) const;
    IQTensor
    makeProjEm(int i) const;
    IQTensor
    makeProjOcc(int i) const;
    IQTensor
    makeProjDou(int i) const;
    IQTensor
    makeProjTri(int i) const;
    IQTensor
    makeProjQua(int i) const;
    IQTensor
    makeProjQui(int i) const;
       
    //Data members -----------------

    int N_;
    bool conserveNf_;

    std::vector<IQIndex> site_;

    };

inline BoseHubbard5::
BoseHubbard5()
    : N_(-1),
    conserveNf_(true)
    { }

inline BoseHubbard5::
BoseHubbard5(int N, const OptSet& opts)
    : N_(N),
      site_(N_+1)
    { 
    conserveNf_ = opts.getBool("ConserveNf",true);
    constructSites();
    }

void inline BoseHubbard5::
constructSites()
    {
    const int One = (conserveNf_ ? 1 : 0);
    const int Two = 2*One;
    const int Three = 3*One;
    const int Four = 4*One;
    const int Five = 5*One;

    for(int j = 1; j <= N_; ++j)
        {
        site_.at(j) = IQIndex(nameint("BoseHubbard5 site=",j),
            Index(nameint("Em for site ",j),1,Site), QN(0,0,0),
            Index(nameint("Occ for site ",j),1,Site), QN(0,One,1),
            Index(nameint("Dou for site ",j),1,Site), QN(0,Two,2),
            Index(nameint("Tri for site ",j),1,Site), QN(0,Three,3),
            Index(nameint("Qua for site ",j),1,Site), QN(0,Four,4),
            Index(nameint("Qui for site ",j),1,Site), QN(0,Five,5));
        }
    }

void inline BoseHubbard5::
doRead(std::istream& s)
    {
    s.read((char*) &N_,sizeof(N_));
    site_.resize(N_+1);
    for(int j = 1; j <= N_; ++j) 
        site_.at(j).read(s);
    if(site_.at(1).qn(2).Nf() == 1)
        conserveNf_ = true;
    else
        conserveNf_ = false;
    }

void inline BoseHubbard5::
doWrite(std::ostream& s) const
    {
    s.write((char*) &N_,sizeof(N_));
    for(int j = 1; j <= N_; ++j) 
        site_.at(j).write(s);
    }

int inline BoseHubbard5::
getN() const
    { return N_; }

inline const IQIndex& BoseHubbard5::
getSi(int i) const
    { return site_.at(i); }

inline IQIndexVal BoseHubbard5::
getState(int i, const String& state) const
    {
    if(state == "0" || state == "Em") 
        {
        return getSi(i)(1);
        }
    else 
    if(state == "1" || state == "Occ") 
        {
        return getSi(i)(2);
        }
    else 
    if(state == "2" || state == "Dou") 
        {
        return getSi(i)(3);
        }
    else
    if(state == "3" || state == "Tri") 
        {
        return getSi(i)(4);
        }
    else
    if(state == "4" || state == "Qua") 
        {
        return getSi(i)(5);
        }
    else
    if(state == "5" || state == "Qui") 
        {
        return getSi(i)(6);
        }
    else
       {
        Error("State " + state + " not recognized");
        return IQIndexVal();
        }
    }

IQIndexVal inline BoseHubbard5::
Em(int i) const
    {
    return getState(i,"0");
    }

IQIndexVal inline BoseHubbard5::
Occ(int i) const
    {
    return getState(i,"1");
    }

IQIndexVal inline BoseHubbard5::
Dou(int i) const
    {
    return getState(i,"2");
    }

IQIndexVal inline BoseHubbard5::
Tri(int i) const
    {
    return getState(i,"3");
    }

IQIndexVal inline BoseHubbard5::
Qua(int i) const
    {
    return getState(i,"4");
    }

IQIndexVal inline BoseHubbard5::
Qui(int i) const
    {
    return getState(i,"5");
    }

IQIndexVal inline BoseHubbard5::
EmP(int i) const
    {
    return primed(getState(i,"0"));
    }

IQIndexVal inline BoseHubbard5::
OccP(int i) const
    {
    return primed(getState(i,"1"));
    }

IQIndexVal inline BoseHubbard5::
DouP(int i) const
    {
    return primed(getState(i,"2"));
    }

IQIndexVal inline BoseHubbard5::
TriP(int i) const
    {
    return primed(getState(i,"3"));
    }

IQIndexVal inline BoseHubbard5::
QuaP(int i) const
    {
    return primed(getState(i,"4"));
    }

IQIndexVal inline BoseHubbard5::
QuiP(int i) const
    {
    return primed(getState(i,"5"));
    }

inline IQTensor BoseHubbard5::
getOp(int i, const String& opname, const OptSet& opts) const
    {
    const
    IQIndex s(si(i));
    const
    IQIndex sP = primed(s);

    IQIndexVal Em(s(1)),
               EmP(sP(1)),
               Occ(s(2)),
               OccP(sP(2)),
               Dou(s(3)),
               DouP(sP(3)),
               Tri(s(4)),
               TriP(sP(4)),
               Qua(s(5)),
               QuaP(sP(5)),
               Qui(s(6)),
               QuiP(sP(6));

    IQTensor Op(conj(s),sP);

    if(opname == "N")
        {
        Op(Occ,OccP) = 1;
        Op(Dou,DouP) = 2;
        Op(Tri,TriP) = 3;
        Op(Qua,QuaP) = 4;
        Op(Qui,QuiP) = 5;
        }
    else
    if(opname == "OS")
        {
        Op(Dou,DouP) = 2;
        Op(Tri,TriP) = 6;
        Op(Qua,QuaP) = 12;
        Op(Qui,QuiP) = 20;
        }//on-site N(N-1)
    else
    if(opname == "B")
        {
        Op(Occ,EmP) = 1; 
        Op(Dou,OccP) = sqrt(2); 
        Op(Tri,DouP) = sqrt(3); 
        Op(Qua,TriP) = sqrt(4); 
        Op(Qui,QuaP) = sqrt(5); 
        }
    else
    if(opname == "Bdag")
        {
        Op(Em,OccP) = 1; 
        Op(Occ,DouP) = sqrt(2);
        Op(Dou,TriP) = sqrt(3);
        Op(Tri,QuaP) = sqrt(4);
        Op(Qua,QuiP) = sqrt(5);
        }
    else
        {
        Error("Operator " + opname + " name not recognized");
        }

    return Op;
    }


//
// The following methods are deprecated but included for 
// backwards compatibility with older Model interface.
//

IQTensor inline BoseHubbard5::
makeN(int i) const
    {
    IQTensor N(conj(si(i)),siP(i));
    N(Occ(i),OccP(i)) = 1;
    N(Dou(i),DouP(i)) = 2;
    N(Tri(i),TriP(i)) = 3;
    N(Qua(i),QuaP(i)) = 4;
    N(Qui(i),QuiP(i)) = 5;
    return N;
    }

IQTensor inline BoseHubbard5::
makeOS(int i) const
    {
    IQTensor OS(conj(si(i)),siP(i));
    OS(Dou(i),DouP(i)) = 2;
    OS(Tri(i),TriP(i)) = 6;
    OS(Qua(i),QuaP(i)) = 12;
    OS(Qui(i),QuiP(i)) = 20;
    return OS;
    }

IQTensor inline BoseHubbard5::
makeB(int i) const
    {
    IQTensor B(conj(si(i)),siP(i));
    B(Occ(i),EmP(i)) = 1;
    B(Dou(i),OccP(i)) = sqrt(2);
    B(Tri(i),DouP(i)) = sqrt(3);
    B(Qua(i),TriP(i)) = sqrt(4);
    B(Qui(i),QuaP(i)) = sqrt(5);
    return B;
    }

IQTensor inline BoseHubbard5::
makeBdag(int i) const
    {
    IQTensor Bdag(conj(si(i)),siP(i));
    Bdag(Em(i),OccP(i)) = 1;
    Bdag(Occ(i),DouP(i)) = sqrt(2);
    Bdag(Dou(i),TriP(i)) = sqrt(3);
    Bdag(Tri(i),QuaP(i)) = sqrt(4);
    Bdag(Qua(i),QuiP(i)) = sqrt(5);
    return Bdag;
    }

#endif
