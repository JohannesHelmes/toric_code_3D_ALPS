#include "haahcode.h"
#include <exception>
#include <iostream>
#include <fstream>
#include <cmath>
#include <alps/osiris/comm.h>
#include <cassert>



using namespace boost;


/************* MEMBERS ********************************/
haahcode::haahcode(const alps::ProcessList& where,const alps::Parameters& p,int node) : alps::scheduler::LatticeMCRun<graph_type>(where,p,node), 
    seed(static_cast<alps::uint32_t>(p["SEED"])),      
    L(static_cast<alps::uint32_t>(p["L"])),      // Linear lattice size
    Nb_Steps(static_cast<alps::uint64_t>(p["SWEEPS"])),    // # of simulation steps
    Nb_Therm_Steps(static_cast<alps::uint64_t>(p["THERMALIZATION"])),
    beta(static_cast<double>(p["beta"])),
    h(static_cast<double>(p.value_or_default("h",0.0))),
    exc(static_cast<alps::uint32_t>(p.value_or_default("ExcType",4))),      // Type of underlying groundstate: 3(plaquettes) 4(vertices) 
    Total_Steps(0)
{

    numsites=num_sites();

    geom.resize(numsites,0);
    sit=sites().first;

    map_lat_to_spin.resize(numsites);
    map_lat_to_plaq.resize(numsites);
    map_lat_to_vert.resize(numsites);
    for (sit=sites().first; sit!=sites().second; ++sit) {
        if (site_type(*sit)<=1) {
            spin_ptr nspin;
            nspin = std::make_shared<spin>(*sit , 0, site_type(*sit));
            spins.push_back(nspin);

            map_lat_to_spin[*sit]=spins.size()-1;
        }
        else if (site_type(*sit)==3) {
            plaq_ptr nplaq = std::make_shared<plaquette>(*sit, site_type(*sit) );
            plaqs.push_back(nplaq);
            map_lat_to_plaq[*sit ]=plaqs.size()-1;
        }
        else if (site_type(*sit)==4) {
            vert_ptr nver = std::make_shared<vertexx>(*sit );
            verts.push_back(nver);
            map_lat_to_vert[*sit ]=verts.size()-1;
        }
    }

    cout<<"numsites "<<numsites<<endl;

    //create all neighbor pairs
    for (sit=sites().first; sit!=sites().second; ++sit) {
        if (site_type(*sit)<=1) {
            for (nit=neighbors(*sit).first; nit!=neighbors(*sit).second; ++nit) {
                if (site_type(*nit)==3) { 
                    spins[map_lat_to_spin[*sit ]]->add_neighbor(static_pointer_cast<plaquette>(plaqs[map_lat_to_plaq[*nit]]));
                    plaqs[map_lat_to_plaq[*nit ]]->add_neighbor(spins[map_lat_to_spin[*sit ]]);
                }
                if (site_type(*nit)==4) { 
                    spins[map_lat_to_spin[*sit ]]->add_neighbor(static_pointer_cast<vertexx>(verts[map_lat_to_vert[*nit ]]));
                    verts[map_lat_to_vert[*nit ]]->add_neighbor(spins[map_lat_to_spin[*sit ]]);
                }
            }
        }
    }
    cout<<"created all neighbors "<<numsites<<endl;

    for (auto s : spins)
        s->copy_neighbors_internally(exc) ;


    measurement_object = std::make_shared<full_energy>(measurements, NofD ); 

    if (exc==3)
        update_object = std::make_shared<single_spin_vert>(seed, 1, beta, spins, plaqs, NofD);
    else if (exc==4)
        update_object = std::make_shared<single_spin_vert>(seed, 1, beta, spins, verts, NofD);

    std::cout << "# L: " << L << " Steps: " << Nb_Steps  << " Spins: " <<spins.size()<<" Sites: "<<numsites<< std::endl;
    //print_information();

}


using ::operator<<;
using ::operator>>;

void haahcode::save(alps::ODump& dump) const
{
    dump << Total_Steps <<NofD <<spins<< plaqs<< verts ; //IMPLEMENT OutStream of spins, plaqs, verts

}

void haahcode::load(alps::IDump& dump)
{
    dump >> Total_Steps >> NofD >>spins>> plaqs >>verts; //IMPLEMENT InStreal of spins, plaqs, verts
}

void haahcode::print_copyright(std::ostream & out)
{
    out << " copyright (c) by Johannes Helmes\n";
}

bool haahcode::is_thermalized() const {
    return (Total_Steps>=Nb_Therm_Steps);
}

double haahcode::work_done() const {
    return (is_thermalized() ? (Total_Steps-Nb_Therm_Steps)/double(Nb_Steps) :0.);
}


void haahcode::dostep() {

    //cout<<"Sweep "<<Total_Steps<<endl;

    update_object->update();
    //cout<<NofD<<endl;

    //TEST for VERTEX METROPOLIS
    /**
    int count;
    for (iit_t plit = plaqs.begin(); plit != plaqs.end(); ++plit) {
        count=0;
        for (spit_t nbit = (*plit)->get_neighbors_begin(); nbit != (*plit)->get_neighbors_end(); ++nbit) {
            count += (*nbit)->get_value();
        }
        if (count %4 == 2)
            cout<<Total_Steps<<": Error at a plaquette "<<endl;
    }
    **/
    
    if (is_thermalized()) 
        measurement_object->measure();
    ++Total_Steps;

}

void haahcode::print_information() {
    for (sit = sites().first; sit!=sites().second; ++sit) {
        cout<<"On the level of ALPS lattice graph ------------------------------------------------------------- "<<endl;
        cout<<"Site "<<*sit<<" has type "<<site_type(*sit)<<", and neighbors ";
        for (nit = neighbors(*sit).first; nit!= neighbors(*sit).second; ++nit) {
            cout<<*nit<<"  ";
        }
        cout<<endl;
        cout<<"On the level of the site object  ------------------------------------------------------------- "<<endl;
        if (site_type(*sit) == 4) {
            int siteindex = map_lat_to_vert[*sit ];
            cout<<"VERTEX, maps to No. "<<siteindex<<" which in turn has name "<<verts[siteindex]->get_name()<<endl;
            cout<<"Neighbors have names ";
            for (spit_t spit = verts[siteindex]->get_neighbors_begin(); spit != verts[siteindex]->get_neighbors_end(); ++spit) {
                cout<<(*spit)->get_name()<<" ";
            }
            cout<<endl;
        }
        else if (site_type(*sit) == 3) {
            int siteindex = map_lat_to_plaq[*sit ];
            cout<<"PLAQUETTE with orientation "<<plaqs[siteindex]->get_orientation()<<", maps to No. "<<siteindex<<" which in turn has name "<<plaqs[siteindex]->get_name()<<endl;
            cout<<"Neighbors have names ";
            for (spit_t spit = plaqs[siteindex]->get_neighbors_begin(); spit != plaqs[siteindex]->get_neighbors_end(); ++spit) {
                cout<<(*spit)->get_name()<<" ";
            }
            cout<<endl;
        }
        else  {
            int siteindex = map_lat_to_spin[*sit ];
            cout<<"SPIN with orientation "<<spins[siteindex]->get_orientation()<<" with geometry "<<spins[siteindex]->get_geometry()<<", maps to No. "<<siteindex<<" which in turn has name "<<spins[siteindex]->get_name()<<endl;
            cout<<"Neighbors have names ";
            for (const_iit_t iit = spins[siteindex]->get_interaction_neighbors_begin(); iit != spins[siteindex]->get_interaction_neighbors_end(); ++iit) {
                cout<<(*iit)->get_name()<<" ";
            }
            cout<<endl;
            cout<<"Initial weight is "<<update_object->get_weight_from_spin(spins[siteindex])<<endl;
        }

        cout<<endl;
    }
}

/****************************************************          MAIN           ***************************************************/
/****************************************************          MAIN           ***************************************************/
/****************************************************          MAIN           ***************************************************/
/****************************************************          MAIN           ***************************************************/
/****************************************************          MAIN           ***************************************************/
/****************************************************          MAIN           ***************************************************/
/****************************************************          MAIN           ***************************************************/

int main(int argc, char** argv)
{
#ifndef BOOST_NO_EXCEPTIONS
try {
#endif

  return alps::scheduler::start(argc,argv,HaahFactory());

#ifndef BOOST_NO_EXCEPTIONS
}
catch (std::exception& exc) {
  std::cerr << exc.what() << "\n";
  alps::comm_exit(true);
  return -1;
}
catch (...) {
  std::cerr << "Fatal Error: Unknown Exception!\n";
  return -2;
}
#endif
}

