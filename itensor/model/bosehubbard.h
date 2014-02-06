//
// Distributed under the ITensor Library License, Version 1.0.
//    (See accompanying LICENSE file.)
//
#ifndef __ITENSOR_HUBBARD_H
#define __ITENSOR_HUBBARD_H
#include "../model.h"

class BoseHubbard : public Model
    {
    public:

    BoseHubbard();

    BoseHubbard(int N, 
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
    EmP(int i) const;

    IQIndexVal
    OccP(int i) const;

    IQIndexVal
    DouP(int i) const;

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

        
    //Data members -----------------

    int N_;
    bool conserveNf_;

    std::vector<IQIndex> site_;

    };

inline BoseHubbard::
BoseHubbard()
    : N_(-1),
    conserveNf_(true)
    { }

inline BoseHubbard::
BoseHubbard(int N, const OptSet& opts)
    : N_(N),
      site_(N_+1)
    { 
    conserveNf_ = opts.getBool("ConserveNf",true);
    constructSites();
    }

void inline BoseHubbard::
constructSites()
    {
    const int One = (conserveNf_ ? 1 : 0);
    const int Two = 2*One;
    for(int j = 1; j <= N_; ++j)
        {
        site_.at(j) = IQIndex(nameint("BoseHubbard site=",j),
            Index(nameint("Em for site ",j),1,Site), QN(0,0,0),
            Index(nameint("Occ for site ",j),1,Site), QN(0,One,1),
            Index(nameint("Dou for site ",j),1,Site), QN(0,Two,2));
        }
    }

void inline BoseHubbard::
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

void inline BoseHubbard::
doWrite(std::ostream& s) const
    {
    s.write((char*) &N_,sizeof(N_));
    for(int j = 1; j <= N_; ++j) 
        site_.at(j).write(s);
    }

int inline BoseHubbard::
getN() const
    { return N_; }

inline const IQIndex& BoseHubbard::
getSi(int i) const
    { return site_.at(i); }

inline IQIndexVal BoseHubbard::
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
        {
        Error("State " + state + " not recognized");
        return IQIndexVal();
        }
    }

IQIndexVal inline BoseHubbard::
Em(int i) const
    {
    return getState(i,"0");
    }

IQIndexVal inline BoseHubbard::
Occ(int i) const
    {
    return getState(i,"1");
    }

IQIndexVal inline BoseHubbard::
Dou(int i) const
    {
    return getState(i,"2");
    }

IQIndexVal inline BoseHubbard::
EmP(int i) const
    {
    return primed(getState(i,"0"));
    }

IQIndexVal inline BoseHubbard::
OccP(int i) const
    {
    return primed(getState(i,"1"));
    }

IQIndexVal inline BoseHubbard::
DouP(int i) const
    {
    return primed(getState(i,"2"));
    }

inline IQTensor BoseHubbard::
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
               DouP(sP(3));

    IQTensor Op(conj(s),sP);

    if(opname == "N")
        {
        Op(Occ,OccP) = 1;
        Op(Dou,DouP) = 2;
        }
    else
    if(opname == "OS")
        {
        Op(Dou,DouP) = 2;
        }//on-site N(N-1)
    else
    if(opname == "B")
        {
        Op(Occ,EmP) = 1; 
        Op(Dou,OccP) = sqrt(2); 
        }
    else
    if(opname == "Bdag")
        {
        Op(Em,OccP) = 1; 
        Op(Occ,DouP) = sqrt(2);
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

IQTensor inline BoseHubbard::
makeN(int i) const
    {
    IQTensor N(conj(si(i)),siP(i));
    N(Occ(i),OccP(i)) = 1;
    N(Dou(i),DouP(i)) = 2;
    return N;
    }

IQTensor inline BoseHubbard::
makeOS(int i) const
    {
    IQTensor OS(conj(si(i)),siP(i));
    OS(Dou(i),DouP(i)) = 2;
    return OS;
    }

IQTensor inline BoseHubbard::
makeB(int i) const
    {
    IQTensor B(conj(si(i)),siP(i));
    B(Occ(i),EmP(i)) = 1;
    B(Dou(i),OccP(i)) = sqrt(2);
    return B;
    }

IQTensor inline BoseHubbard::
makeBdag(int i) const
    {
    IQTensor Bdag(conj(si(i)),siP(i));
    Bdag(Em(i),OccP(i)) = 1;
    Bdag(Occ(i),DouP(i)) = sqrt(2);
    return Bdag;
    }

#endif
