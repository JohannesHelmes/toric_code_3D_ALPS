#include "toriccode.h"
#include <exception>
#include <iostream>
#include <fstream>
#include <cmath>
#include <alps/osiris/comm.h>



using namespace boost;


/************* MEMBERS ********************************/
toriccode::toriccode(const alps::ProcessList& where,const alps::Parameters& p,int node) : alps::scheduler::LatticeMCRun<graph_type>(where,p,node), 
    seed(static_cast<alps::uint32_t>(p["SEED"])),      
    L(static_cast<alps::uint32_t>(p["L"])),      // Linear lattice size
    Nb_Steps(static_cast<alps::uint64_t>(p["SWEEPS"])),    // # of simulation steps
    Nb_Therm_Steps(static_cast<alps::uint64_t>(p["THERMALIZATION"])),
    beta(static_cast<double>(p["beta"])),
    ratio(static_cast<double>(p.value_or_default("ratio",1.0))),    // not useful for thermodynamic integration
    n(static_cast<alps::uint32_t>(p.value_or_default("n",2))),      // Renyi index
    exc(static_cast<alps::uint32_t>(p.value_or_default("ExcType",2))),      // Type of excitation: 1(plaquettes) 2(vertices) 
    algo(static_cast<alps::uint32_t>(p.value_or_default("Algorithm",1))),      // local updates (1),  deconfined updates (2) 
    measure(static_cast<alps::uint32_t>(p.value_or_default("Measurement",1))),      // thermodynamic int (1),  ensemble switching (2), specific heat (3)
    Total_Steps(0),
    IncStep(static_cast<string>(p["IncStep"]))
{

    numsites=num_sites();

    geom.resize(numsites,0);
    sit=sites().first;
    if (n > 1) {  //replica trick only senseful for n>=2 
        for (int i=0; i<IncStep.length(); ++i,++sit) {
            if (sit==sites().second)
                break;
            while (site_type(*sit)!=0)
                ++sit;
            if (sit==sites().second)
                break;
            geom[*sit]=(int)(IncStep[i]-'0');
        }
    }

    map_lat_to_spin.resize(n*numsites);
    map_lat_to_plaq.resize(n*numsites);
    map_lat_to_vert.resize(n*numsites);
    cout<<"numsites "<<numsites<<" total "<<n*numsites<<endl;
    for (int i=0; i<n; ++i) {
        for (sit=sites().first; sit!=sites().second; ++sit) {
            if (site_type(*sit)==0) {
                if ((geom[*sit]!=1)||(i==0)) {
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
                    if (site_type(*nit)==1) { 
                        spins[map_lat_to_spin[*sit + i*numsites]]->add_neighbor(plaqs[map_lat_to_plaq[*nit + i*numsites]]);
                        plaqs[map_lat_to_plaq[*nit + i*numsites]]->add_neighbor(spins[map_lat_to_spin[*sit + i*numsites]]);
                    }
                    if (site_type(*nit)==2) { 
                        spins[map_lat_to_spin[*sit + i*numsites]]->add_neighbor(verts[map_lat_to_vert[*nit + i*numsites]]);
                        verts[map_lat_to_vert[*nit + i*numsites]]->add_neighbor(spins[map_lat_to_spin[*sit + i*numsites]]);
                    }
                }
            }
        }
    }

    //label all connected regions and boundaries of vertices/plaquettes
    int counter;
    for (int i=0; i<n; ++i) {
        for (sit=sites().first; sit!=sites().second; ++sit) {
            if (site_type(*sit)==2) {
                counter=0;
                for (nit=neighbors(*sit).first; nit!=neighbors(*sit).second; ++nit) {
                    if (geom[*nit]!=1) {
                        verts[map_lat_to_vert[*sit + i*numsites]]->add_label(geom[*nit]);
                        verts[map_lat_to_vert[*sit + i*numsites]]->set_boundary(true);
                        ++counter;
                    }
                }
                if (counter==0) {
                        verts[map_lat_to_vert[*sit + i*numsites]]->add_label(1); //means completely in subsystem A
                        verts[map_lat_to_vert[*sit + i*numsites]]->set_boundary(false); 
                }
                else if (counter==6) {
                        verts[map_lat_to_vert[*sit + i*numsites]]->set_boundary(false); //means completely in subsystem B
                }
            }
        }
    }

    /*
    cout<<"Check labels and boundaries of all vertices"<<endl;
    counter=0;
    for (vit_t vit=verts.begin(); vit!=verts.end(); ++vit, ++counter)
        cout<<counter<<": has label "<<(*vit)->get_label()<<" and boundary "<<(*vit)->get_boundary()<<endl;
        */


    NofD=(exc==1)? -plaqs.size() : -verts.size();
    if (measure == 1) {
        measurement_object = std::make_shared<thermo_int>(measurements, NofD, spins.size() );
    }

    if (algo==1) {
        if (exc==1)
            update_object = std::make_shared<single_spin_plaq>(seed, n, beta, spins, plaqs, NofD); //spins, plaqs and NofD are referenced
        else if (exc==2) {
            if (ratio==1.0)
                update_object = std::make_shared<single_spin_vert>(seed, n, beta, spins, verts, NofD);
            else  //not useful for thermodynamic int
                update_object = std::make_shared<mix_spin_plaq_for_vert>(seed, n, beta, spins, plaqs, verts, NofD, ratio);
        }
    }
    else if (algo==2) {
        update_object = std::make_shared<deconfined_vert>(seed, n, beta, spins, verts, NofD); 
    }

    numspins=spins.size(); //get rid of this later
    std::cout << "# L: " << L << " Steps: " << Nb_Steps  << " Spins: " <<spins.size()<<" Sites: "<<numsites<< std::endl;

}


using ::operator<<;
using ::operator>>;

void toriccode::save(alps::ODump& dump) const
{
    dump << Total_Steps <<NofD <<spins<< plaqs<< verts ; //IMPLEMENT OutStream of spins, plaqs, verts

}

void toriccode::load(alps::IDump& dump)
{
    dump >> Total_Steps >> NofD >>spins>> plaqs >>verts; //IMPLEMENT InStreal of spins, plaqs, verts
}

void toriccode::print_copyright(std::ostream & out)
{
    out << " copyright (c) by Johannes Helmes\n";
}

bool toriccode::is_thermalized() const {
    return (Total_Steps>=Nb_Therm_Steps);
}

double toriccode::work_done() const {
    return (is_thermalized() ? (Total_Steps-Nb_Therm_Steps)/double(Nb_Steps) :0.);
}


void toriccode::dostep() {

    //cout<<"Sweep "<<Total_Steps<<endl;

    update_object->update();
    //cout<<NofD<<endl;
    
    if (is_thermalized()) 
        measurement_object->measure();
    ++Total_Steps;
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

  return alps::scheduler::start(argc,argv,ToricFactory());

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

