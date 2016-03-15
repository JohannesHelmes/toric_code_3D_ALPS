#include "toric3DFTMembraneSingleSpin.h"
#include <exception>
#include <iostream>
#include <fstream>
#include <cmath>



using namespace boost;


/************* MEMBERS ********************************/
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
    map_lat_to_vert.resize(n*numsites);
    cout<<"numsites "<<numsites<<" total "<<n*numsites<<endl;
    for (int i=0; i<n; ++i) {
        for (sit=sites().first; sit!=sites().second; ++sit) {
            if (site_type(*sit)==0) {
                if ((!geom[*sit])||(i==0)) {
                    spin_ptr nspin = std::make_shared<spin>();
                    spins.push_back(nspin);
                }
                else 
                    spins.push_back(spins[map_lat_to_spin[*sit]]);

                map_lat_to_spin[*sit + i*numsites]=spins.size()-1;
            }
            else if (site_type(*sit)==1) {
                plaq_ptr nplaq = std::make_shared<plaquette>();
                plaqs.push_back(nplaq);
                map_lat_to_plaq[*sit + i*numsites]=plaqs.size()-1;
            }
            else if (site_type(*sit)==2) {
                vert_ptr nver = std::make_shared<vertexx>();
                verts.push_back(nver);
                map_lat_to_vert[*sit + i*numsites]=verts.size()-1;
            }
        }
    }

    //create all neighbor pairs
    for (int i=0; i<n; ++i) {
        for (sit=sites().first; sit!=sites().second; ++sit) {
            if (site_type(*sit)==0) {
                for (nit=neighbors(*sit).first; nit!=neighbors(*sit).second; ++nit) {
                    if (site_type(*nit)==1) { //plaqs is not a good variable name - can be vertices as well
                        spins[map_lat_to_spin[*sit + i*numsites]]->add_neighbor(plaqs[map_lat_to_plaq[*nit + i*numsites]]);
                        plaqs[map_lat_to_plaq[*nit + i*numsites]]->add_neighbor(spins[map_lat_to_spin[*sit + i*numsites]]);
                    }
                    if (site_type(*nit)==2) { //plaqs is not a good variable name - can be vertices as well
                        spins[map_lat_to_spin[*sit + i*numsites]]->add_neighbor(verts[map_lat_to_vert[*nit + i*numsites]]);
                        verts[map_lat_to_vert[*nit + i*numsites]]->add_neighbor(spins[map_lat_to_spin[*sit + i*numsites]]);
                    }
                }
            }
        }
    }

    /*
    for (sit=sites().first; sit!=sites().second; ++sit) {
        if (site_type(*sit)==0) { 
            cout<<*sit<<" has spin index "<<map_lat_to_spin[*sit]<<endl;
            cout<<"IN A = "<<geom[*sit]<<" and the weight is "<<spins[map_lat_to_spin[*sit]]->get_weight()<<" and has num neighbors"<<spins[map_lat_to_spin[*sit]]->num_neighbors()<<endl;
        }
    }
    */

    expmB.resize(8*n+1);
    expmB[0]=1.0;
    for (int i=1; i<=8*n; ++i) 
        expmB[i]=std::exp(-1.*i*beta*B);

    //NofD=-n*N;
    NofD=(exc==1)? -plaqs.size() : -verts.size();
    //cout<<"initial NofD "<<NofD<<" vs "<<plaqs.size()<<endl;
    numspins=spins.size();
    std::cout << "# L: " << L << " Steps: " << Nb_Steps  << " Spins: " <<numspins<<" Sites: "<<numsites<< std::endl;

}


using ::operator<<;
using ::operator>>;

void toricFTSPMembrane::save(alps::ODump& dump) const
{
    dump << Total_Steps <<spins<< plaqs<< verts ; //IMPLEMENT OutStream of spins, plaqs, verts

}

void toricFTSPMembrane::load(alps::IDump& dump)
{
    dump >> Total_Steps >>spins>> plaqs >>verts; //IMPLEMENT InStreal of spins, plaqs, verts
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


void toricFTSPMembrane::dostep() {

    //insert or remove electric defects
    //this code is for zero-field only!!

    for (int j=0; j<numspins/7; ++j) {

        candidate=spins[random_int(spins.size())];
        int cand_weight = (exc==1)? candidate->get_weight_from_plaqs() : candidate->get_weight_from_verts();
        //cout<<Total_Steps<<": Try to flip  with ediff"<<cand_weight<<" and weight "<<expmB[-2*cand_weight]<<endl;

        if ((cand_weight>=0)||(random_01()<expmB[-2*cand_weight])) {
            NofD -= 2*cand_weight; //is the old weight
            if (exc==1)
                candidate->flip_and_flip_plaqs();
            else
                candidate->flip_and_flip_verts();
        }
    }
    //cout<<NofD<<endl;

    
    if (is_thermalized()) 
        do_measurements();
    ++Total_Steps;
}



