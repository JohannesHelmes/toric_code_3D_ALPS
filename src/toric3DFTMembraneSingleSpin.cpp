#include "toric3DFTMembraneSingleSpin.h"
#include <exception>
#include <iostream>
#include <fstream>
#include <cmath>



using namespace boost;


toricFTSPMembrane::toricFTSPMembrane(const alps::ProcessList& where,const alps::Parameters& p,int node) : alps::scheduler::LatticeMCRun<graph_type>(where,p,node), 
    L(static_cast<alps::uint32_t>(p["L"])),      // Linear lattice size
    Nb_Steps(static_cast<alps::uint64_t>(p["SWEEPS"])),    // # of simulation steps
    Nb_Therm_Steps(static_cast<alps::uint64_t>(p["THERMALIZATION"])),
    beta(static_cast<double>(p["beta"])),
    B(static_cast<double>(p.value_or_default("B",1.0))),
    n(static_cast<alps::uint32_t>(p.value_or_default("n",2))),      // Renyi index
    d(static_cast<alps::uint32_t>(p.value_or_default("d",3))),      // dimension 
    exc(static_cast<alps::uint32_t>(p.value_or_default("ExcType",2))),      // Type of excitation: 1(plaquettes) 2(vertices) 
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

    geom.resize(numsites,false);

    sit=sites().first;
    for (int i=0; i<IncStep.length(); ++i,++sit) {
        if (sit==sites().second)
            break;
        while (site_type(*sit)!=0)
            ++sit;
        if (sit==sites().second)
            break;
        if (IncStep[i]!='0') {
            geom[*sit]=true;
        }
         
    }

    map_lat_to_spin.resize(n*numsites);
    map_lat_to_plaq.resize(n*numsites);
    cout<<"numsites "<<numsites<<" total "<<n*numsites<<endl;
    for (int i=0; i<n; ++i) {
        for (sit=sites().first; sit!=sites().second; ++sit) {
            if (site_type(*sit)==0) {
                if ((!geom[*sit])||(i==0)) {
                    spin_ptr nspin = std::make_shared<spin>((n-1)*geom[*sit]+1);
                    spins.push_back(nspin);
                }
                else 
                    spins.push_back(spins[map_lat_to_spin[*sit]]);

                map_lat_to_spin[*sit + i*numsites]=spins.size()-1;
            }
            else if (site_type(*sit)==exc) {
                plaq_ptr nplaq = std::make_shared<plaquette>();
                plaqs.push_back(nplaq);
                map_lat_to_plaq[*sit + i*numsites]=plaqs.size()-1;
            }
        }
    }

    //create all neighbor pairs
    for (int i=0; i<n; ++i) {
        for (sit=sites().first; sit!=sites().second; ++sit) {
            if (site_type(*sit)==0) {
                for (nit=neighbors(*sit).first; nit!=neighbors(*sit).second; ++nit) {
                    if (site_type(*nit)==exc) { //plaqs is not a good variable name - can be vertices as well
                        spins[map_lat_to_spin[*sit + i*numsites]]->add_neighbor(plaqs[map_lat_to_plaq[*nit + i*numsites]]);
                        plaqs[map_lat_to_plaq[*nit + i*numsites]]->add_neighbor(spins[map_lat_to_spin[*sit + i*numsites]]);
                    }
                }
            }
        }
    }

    for (sit=sites().first; sit!=sites().second; ++sit) {
        if (site_type(*sit)==0) { 
            cout<<*sit<<" has spin index "<<map_lat_to_spin[*sit]<<endl;
            cout<<"IN A = "<<geom[*sit]<<" and the weight is "<<spins[map_lat_to_spin[*sit]]->get_weight()<<" and has num neighbors"<<spins[map_lat_to_spin[*sit]]->num_neighbors()<<endl;
        }
    }
    std::cout << "# L: " << L << " Steps: " << Nb_Steps  << " Spins: " <<numspins<<" Sites: "<<numsites<< std::endl;

    expmB.resize(8*n+1);
    expmB[0]=1.0;
    for (int i=1; i<=8*n; ++i) 
        expmB[i]=std::exp(-1.*i*beta*B);

    //NofD=-n*N;
    //cout<<"initial NofD "<<NofD<<" vs "<<plaqs.size()<<endl;
    NofD=-plaqs.size();

}

void toricFTSPMembrane::save(alps::ODump& dump) const
{
    dump << Total_Steps<<spins_<<plaquette_defects;
}

void toricFTSPMembrane::load(alps::IDump& dump)
{
    plaquette_defects.clear();
    plaquette_defects.resize(n);
    dump >> Total_Steps>>spins_>>plaquette_defects;
//    if(!where.empty())
//        dump >> spins_;
}

void toricFTSPMembrane::print_copyright(std::ostream & out)
{
    out << " copyright (c) by Johannes Helmes\n";
}

bool toricFTSPMembrane::is_thermalized() const {
    return (Total_Steps>=Nb_Therm_Steps);
}

double toricFTSPMembrane::work_done() const {
    return (is_thermalized() ? (Total_Steps-Nb_Therm_Steps)/double(Nb_Steps) :0.);
}


void toricFTSPMembrane::flip(int replica, int spin) {
    if (geom[spin]) {
        for (int i=0; i<n; ++i) 
            spins_[numsites*i+spin]=!spins_[numsites*i+spin];
    }
    else 
        spins_[numsites*replica+spin]=!spins_[numsites*replica+spin];
}

void toricFTSPMembrane::flip_defect(int replica, int plaquette) {
    if (plaquette_defects[replica][plaquette]) {
        plaquette_defects[replica][plaquette]=false;
        NofD-=2;
    }
    else {
        plaquette_defects[replica][plaquette]=true;
        NofD+=2;
    }
}

void toricFTSPMembrane::dostep() {

    //insert or remove electric defects
    //this code is for zero-field only!!

    for (int j=0; j<N/2; ++j) {

        candidate=spins[random_int(spins.size())];
        int cand_weight = candidate->get_weight();
        //cout<<Total_Steps<<": Try to flip "<<cand<<" with ediff"<<candidate->get_weight()<<" and weight "<<expmB[-2*candidate->get_weight()]<<endl;

        if ((cand_weight>=0)||(random_01()<expmB[-2*cand_weight])) {
            NofD -= 2*cand_weight; //is the old weight
            candidate->flip();
        }
        //cout<<NofD<<endl;


        /*
        replica=random_int(n);

        if (d==3)
            start=7*random_int(N)+4+random_int(3); //choose the spin to flip 
        if (d==2)
            start=4*random_int(N)+2+random_int(2); //choose the spin to flip
        
        weight=1.0;

        pneighs.clear();
        for (nit=neighbors(start).first; nit!=neighbors(start).second; ++nit) {
            if (site_type(*nit)==1) {
                pneighs.push_back(*nit);
                if (IsInA(start)) {
                    for (int i=0; i<n; ++i) {
                        weight*=plaquette_defects[i][*nit]? 1 : expmB[2];
                    }
                }
                else {
                    weight*=plaquette_defects[replica][*nit]? 1 : expmB[2];
                }
            }
        }

        if ((weight>=1)||(random_01()<weight)) {
            flip(replica,start);
            if (IsInA(start)) {
                for (int i=0; i<n; ++i) {
                    flip_defect(i,pneighs[0]);
                    flip_defect(i,pneighs[1]);
                    flip_defect(i,pneighs[2]);
                    flip_defect(i,pneighs[3]);
                }
            }
            else {
                flip_defect(replica,pneighs[0]);
                flip_defect(replica,pneighs[1]);
                flip_defect(replica,pneighs[2]);
                flip_defect(replica,pneighs[3]);
            }
        }
        */
    }

    
    if (is_thermalized()) 
        do_measurements();
    ++Total_Steps;
}


bool toricFTSPMembrane::IsInA(site_descriptor where) {
    return (geom[where]);
}

