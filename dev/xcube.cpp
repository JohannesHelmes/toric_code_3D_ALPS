#include "xcube.h"
#include <exception>
#include <iostream>
#include <fstream>
#include <cmath>
#include <alps/osiris/comm.h>
#include <cassert>



using namespace boost;


/************* MEMBERS ********************************/
xcube::xcube(const alps::ProcessList& where,const alps::Parameters& p,int node) : alps::scheduler::LatticeMCRun<graph_type>(where,p,node), 
    seed(static_cast<alps::uint32_t>(p["SEED"])),      
    L(static_cast<alps::uint32_t>(p["L"])),      // Linear lattice size
    Nb_Steps(static_cast<alps::uint64_t>(p["SWEEPS"])),    // # of simulation steps
    Nb_Therm_Steps(static_cast<alps::uint64_t>(p["THERMALIZATION"])),
    beta(static_cast<double>(p["beta"])),
    h(static_cast<double>(p.value_or_default("h",0.0))),
    Total_Steps(0)
{

    numsites=num_sites();

    sit=sites().first;

    map_lat_to_spin.resize(numsites);
    map_lat_to_cube.resize(numsites);
    for (sit=sites().first; sit!=sites().second; ++sit) {
        if (site_type(*sit)<=2) {
            spin_ptr nspin;
            nspin = std::make_shared<spin>(*sit , 0, site_type(*sit));
            spins.push_back(nspin);

            map_lat_to_spin[*sit ]=spins.size()-1;
        }
        else if (site_type(*sit)==7) {
            vert_ptr nver = std::make_shared<vertexx>(*sit );
            cubes.push_back(nver);
            map_lat_to_cube[*sit ]=cubes.size()-1;
        }
    }
    cout<<"numsites "<<numsites<<endl;

    //create all neighbor pairs
    for (sit=sites().first; sit!=sites().second; ++sit) {
        if (site_type(*sit)<=2) {
            for (nit=neighbors(*sit).first; nit!=neighbors(*sit).second; ++nit) {
                if (site_type(*nit)==7) { 
                    spins[map_lat_to_spin[*sit ]]->add_neighbor(static_pointer_cast<vertexx>(cubes[map_lat_to_cube[*nit ]]));
                    cubes[map_lat_to_cube[*nit ]]->add_neighbor(spins[map_lat_to_spin[*sit ]]);
                }
            }
        }
    }

    for (auto s : spins)
        s->copy_neighbors_internally(4) ;

    measurement_object = std::make_shared<full_energy>(measurements, NofD ); 
    //We only consider the energy from the unreplicated model

    update_object = std::make_shared<single_spin_vert>(seed, 1, beta, spins, cubes, NofD);
    //We only perform single spin flips, cubes can be regarded as vertices


    std::cout << "# L: " << L << " Steps: " << Nb_Steps  << " Spins: " <<spins.size()<<" Sites: "<<numsites<< std::endl;
    //print_information();

}


using ::operator<<;
using ::operator>>;

void xcube::save(alps::ODump& dump) const
{
    dump << Total_Steps <<NofD <<spins ; //IMPLEMENT OutStream of spins, plaqs, verts

}

void xcube::load(alps::IDump& dump)
{
    dump >> Total_Steps >> NofD >>spins; //IMPLEMENT InStreal of spins, plaqs, verts
}

void xcube::print_copyright(std::ostream & out)
{
    out << " copyright (c) by Johannes Helmes\n";
}

bool xcube::is_thermalized() const {
    return (Total_Steps>=Nb_Therm_Steps);
}

double xcube::work_done() const {
    return (is_thermalized() ? (Total_Steps-Nb_Therm_Steps)/double(Nb_Steps) :0.);
}


void xcube::dostep() {

    //cout<<"Sweep "<<Total_Steps<<endl;

    update_object->update();
    //cout<<Total_Steps<<": "<<NofD<<endl;

    
    if (is_thermalized()) 
        measurement_object->measure();
    ++Total_Steps;

}

void xcube::print_information() {
    for (sit = sites().first; sit!=sites().second; ++sit) {
        cout<<"On the level of ALPS lattice graph ------------------------------------------------------------- "<<endl;
        cout<<"Site "<<*sit<<" has type "<<site_type(*sit)<<", and neighbors ";
        for (nit = neighbors(*sit).first; nit!= neighbors(*sit).second; ++nit) {
            cout<<*nit<<"  ";
        }
        cout<<endl;
        cout<<"On the level of the site object  ------------------------------------------------------------- "<<endl;
        if (site_type(*sit) == 7) {
            int siteindex = map_lat_to_cube[*sit ];
            cout<<"Cube, maps to No. "<<siteindex<<" which in turn has name "<<cubes[siteindex]->get_name()<<endl;
            cout<<"Neighbors have names ";
            for (spit_t spit = cubes[siteindex]->get_neighbors_begin(); spit != cubes[siteindex]->get_neighbors_end(); ++spit) {
                cout<<(*spit)->get_name()<<" ";
            }
            cout<<endl;
        }
        else if (site_type(*sit) <= 2)  {
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

  return alps::scheduler::start(argc,argv,XCubeFactory());

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

