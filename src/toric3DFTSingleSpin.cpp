#include "toric3DFTSingleSpin.h"
#include <exception>
#include <iostream>
#include <fstream>
#include <cmath>



using namespace boost;


toricFTSP::toricFTSP(const alps::ProcessList& where,const alps::Parameters& p,int node) : alps::scheduler::LatticeMCRun<graph_type>(where,p,node), 
    L(static_cast<alps::uint32_t>(p["L"])),      // Linear lattice size
    Nb_Steps(static_cast<alps::uint64_t>(p["SWEEPS"])),    // # of simulation steps
    Nb_Therm_Steps(static_cast<alps::uint64_t>(p["THERMALIZATION"])),
    beta(static_cast<double>(p["beta"])),
    B(static_cast<double>(p.value_or_default("B",1.0))),
    n(static_cast<alps::uint32_t>(p.value_or_default("n",2))),      // Renyi index
    d(static_cast<alps::uint32_t>(p.value_or_default("d",3))),      // dimension 
    Total_Steps(0),
    IncStep(static_cast<string>(p["IncStep"]))
{

    numsites=num_sites();
    if (d==2) {
        N=numsites/4;
        numspins=n*2*N;
    }
    else {
        N=pow(L,d);
        numspins=n*3*N;
    }

    spins_.resize(n*numsites,0);
    vertex_defects.resize(n);
    for (int i=0; i<n; ++i) 
        vertex_defects[i].resize(numsites,false); 
    geom.resize(numsites,false);
    edge.resize(numsites,false);
    
    std::cout << "# L: " << L << " Steps: " << Nb_Steps  << " Spins: " <<numspins<<" Sites: "<<numsites<< std::endl;

    expmB.resize(2*n+1);
    expmB[0]=1.0;
    for (int i=1; i<=n; ++i) 
        expmB[2*i]=std::exp(-2.*i*beta*B);

    NofD=-n*N;

}

void toricFTSP::save(alps::ODump& dump) const
{
    dump << Total_Steps<<spins_<<vertex_defects;
}

void toricFTSP::load(alps::IDump& dump)
{
    vertex_defects.clear();
    vertex_defects.resize(n);
    dump >> Total_Steps>>spins_>>vertex_defects;
//    if(!where.empty())
//        dump >> spins_;
}

void toricFTSP::print_copyright(std::ostream & out)
{
    out << " copyright (c) by Johannes Helmes\n";
}

bool toricFTSP::is_thermalized() const {
    return (Total_Steps>=Nb_Therm_Steps);
}

double toricFTSP::work_done() const {
    return (is_thermalized() ? (Total_Steps-Nb_Therm_Steps)/double(Nb_Steps) :0.);
}


void toricFTSP::flip(int replica, int spin) {
    if (geom[spin]) {
        for (int i=0; i<n; ++i) 
            spins_[numsites*i+spin]=!spins_[numsites*i+spin];
    }
    else 
        spins_[numsites*replica+spin]=!spins_[numsites*replica+spin];
}

void toricFTSP::flip_defect(int replica, int vertex) {
    if (vertex_defects[replica][vertex]) {
        vertex_defects[replica][vertex]=false;
        NofD-=2;
    }
    else {
        vertex_defects[replica][vertex]=true;
        NofD+=2;
    }
    /*
    if (geom[vertex])
        vertex_defects[replica^1][vertex]=!vertex_defects[replica^1][vertex];
        */
}

void toricFTSP::dostep() {

    //insert or remove electric defects
    //this code is for zero-field only!!

    for (int j=0; j<N/2; ++j) {
        replica=random_int(n);

        if (d==3)
            start=7*random_int(N)+4+random_int(3); //choose the spin to flip 
        if (d==2)
            start=4*random_int(N)+2+random_int(2); //choose the spin to flip
        
        weight=1.0;

        vneighs.clear();
        for (nit=neighbors(start).first; nit!=neighbors(start).second; ++nit) {
            if (site_type(*nit)==2) {
                vneighs.push_back(*nit);
                if (IsInA(start)) {
                    for (int i=0; i<n; ++i) {
                        weight*=vertex_defects[i][*nit]? 1 : expmB[2];
                    }
                }
                else {
                    weight*=vertex_defects[replica][*nit]? 1 : expmB[2];
                }
            }
        }

        if ((weight>=1)||(random_01()<weight)) {
            flip(replica,start);
            if (IsInA(start)) {
                for (int i=0; i<n; ++i) {
                    flip_defect(i,vneighs[0]);
                    flip_defect(i,vneighs[1]);
                }
            }
            else {
                flip_defect(replica,vneighs[0]);
                flip_defect(replica,vneighs[1]);
            }
        }
    }

    
    if (is_thermalized()) 
        do_measurements();
    ++Total_Steps;
}


bool toricFTSP::IsInA(site_descriptor where) const {
    return (geom[where]);
}

