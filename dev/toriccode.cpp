#include "toriccode.h"
#include <exception>
#include <iostream>
#include <fstream>
#include <cmath>
#include <alps/osiris/comm.h>
#include <cassert>



using namespace boost;


/************* MEMBERS ********************************/
toriccode::toriccode(const alps::ProcessList& where,const alps::Parameters& p,int node) : alps::scheduler::LatticeMCRun<graph_type>(where,p,node), 
    seed(static_cast<alps::uint32_t>(p["SEED"])),      
    L(static_cast<alps::uint32_t>(p["L"])),      // Linear lattice size
    W(static_cast<alps::uint32_t>(p["W"])),      // Linear lattice size in z direction
    Nb_Steps(static_cast<alps::uint64_t>(p["SWEEPS"])),    // # of simulation steps
    Nb_Therm_Steps(static_cast<alps::uint64_t>(p["THERMALIZATION"])),
    beta(static_cast<double>(p["beta"])),
    h(static_cast<double>(p.value_or_default("h",0.0))),
    hz(static_cast<double>(p.value_or_default("hz",h))),
    n(static_cast<alps::uint32_t>(p.value_or_default("n",2))),      // Renyi index
    exc(static_cast<alps::uint32_t>(p.value_or_default("ExcType",4))),      // Type of underlying groundstate: 3(plaquettes) 4(vertices) 
    algo(static_cast<alps::uint32_t>(p.value_or_default("Algorithm",1))),      
        // local updates (1),  deconfined updates (2), single-vertex-flips (3), non-isotropic single vertex flips (4), vertex_wolff (5), non-isotropic vertex wolff (6)
    measure(static_cast<alps::uint32_t>(p.value_or_default("Measurement",1))),      
        // thermodynamic int (1),  ensemble switching (2), thermodynamic int with h (3), full_energy for specific heat (4)
    Total_Steps(0),
    IncStep(static_cast<string>(p.value_or_default("IncStep","")))
{

    numsites=num_sites();

    geom.resize(numsites,0);
    sit=sites().first;
    if (n > 1) {  //replica trick only senseful for n>=2 
        if ( (algo == 4) || (algo == 6) ) {
            for (int w=0; w<W; ++w) {
                for (int i=0; i<IncStep.length(); ++i,++sit) {
                    if (sit==sites().second)
                        break;
                    while (site_type(*sit)>2)
                        ++sit;
                    if (sit==sites().second) {
                        cout<<"ERROR: Geometry information and lattice dont fit !!"<<endl;
                        break;
                    }
                    geom[*sit]=(int)(IncStep[i]-'0');
                }
            }
        }
        else {
            for (int i=0; i<IncStep.length(); ++i,++sit) {
                if (sit==sites().second)
                    break;
                while (site_type(*sit)>2)
                    ++sit;
                if (sit==sites().second)
                    break;
                geom[*sit]=(int)(IncStep[i]-'0');
            }
        }
    }
    //cout<<"Achieved site no "<<*sit<<endl;

    map_lat_to_spin.resize(n*numsites);
    map_lat_to_plaq.resize(n*numsites);
    map_lat_to_vert.resize(n*numsites);
    for (int i=0; i<n; ++i) {
        for (sit=sites().first; sit!=sites().second; ++sit) {
            if (site_type(*sit)<=2) {
                if ((geom[*sit]!=1)||(i==0)) {
                    spin_ptr nspin;
                    nspin = std::make_shared<spin>(*sit + i*numsites, geom[*sit], site_type(*sit));
                    spins.push_back(nspin);
                }
                else 
                    spins.push_back(spins[map_lat_to_spin[*sit]]);

                map_lat_to_spin[*sit + i*numsites]=spins.size()-1;
            }
            else if (site_type(*sit)<=5) {
                plaq_ptr nplaq = std::make_shared<plaquette>(*sit + i*numsites, site_type(*sit) );
                plaqs.push_back(nplaq);
                map_lat_to_plaq[*sit + i*numsites]=plaqs.size()-1;
            }
            else if (site_type(*sit)==6) {
                vert_ptr nver = std::make_shared<vertexx>(*sit + i*numsites);
                verts.push_back(nver);
                map_lat_to_vert[*sit + i*numsites]=verts.size()-1;
            }
        }
    }
    cout<<"numsites "<<numsites<<" total "<<n*numsites<<endl;

    //create all neighbor pairs
    for (int i=0; i<n; ++i) {
        for (sit=sites().first; sit!=sites().second; ++sit) {
            if (site_type(*sit)<=2) {
                spins[map_lat_to_spin[*sit + i*numsites]]->set_ninr(spins[map_lat_to_spin[*sit + ((i+1)%n)*numsites]] );
                
                //spins[map_lat_to_spin[(*sit + 1)%L + i*numsites]]->set_ninr(spins[map_lat_to_spin[*sit + ((i+1)%n)*numsites]] );

                //cout<<"Geom "<<geom[*sit]<<", will be a pointer to itself "<< ( spins[map_lat_to_spin[*sit + i*numsites]] == spins[map_lat_to_spin[*sit + ((i+1)%n)*numsites]] ) <<endl;

                for (nit=neighbors(*sit).first; nit!=neighbors(*sit).second; ++nit) {
                    if (site_type(*nit)<=5) { 
                        spins[map_lat_to_spin[*sit + i*numsites]]->add_neighbor(static_pointer_cast<plaquette>(plaqs[map_lat_to_plaq[*nit + i*numsites]]));
                        plaqs[map_lat_to_plaq[*nit + i*numsites]]->add_neighbor(spins[map_lat_to_spin[*sit + i*numsites]]);
                        plaqs[map_lat_to_plaq[*nit + i*numsites]]->set_ninr(plaqs[map_lat_to_plaq[*nit + ((i+1)%n)*numsites]] );
                    }
                    if (site_type(*nit)==6) { 
                        spins[map_lat_to_spin[*sit + i*numsites]]->add_neighbor(static_pointer_cast<vertexx>(verts[map_lat_to_vert[*nit + i*numsites]]));
                        verts[map_lat_to_vert[*nit + i*numsites]]->add_neighbor(spins[map_lat_to_spin[*sit + i*numsites]]);
                        verts[map_lat_to_vert[*nit + i*numsites]]->set_ninr(verts[map_lat_to_vert[*nit + ((i+1)%n)*numsites]] );

                    }
                }
            }
            else if (site_type(*sit) == 6) {
                //set all lattice neighbors
                //THIS PART IS ONLY COMPATIBLE WITH THE 3D TORIC CODE !!!!!!!!!!!!!!!!!!!!!!!!!!!
                int x_off = ((*sit + 8)%(8*L) == 0 )? *sit - 8*(L-1) : *sit+8 ;
                x_off = x_off%numsites;
                verts[map_lat_to_vert[*sit + i*numsites]]->set_nindir(verts[map_lat_to_vert[x_off + i*numsites]] ,0 );
                int y_off = ((*sit + 8*L)%(8*L*L) == 0 )? *sit - 8*L*(L-1) : *sit+8*L ;
                y_off = y_off%numsites;
                verts[map_lat_to_vert[*sit + i*numsites]]->set_nindir(verts[map_lat_to_vert[y_off + i*numsites]] ,1 );
                int z_off = ((*sit + 8*L*L)%(8*L*L*W) == 0 )? *sit - 8*L*L*(W-1) : *sit+8*L*L ;
                z_off = z_off%numsites;
                verts[map_lat_to_vert[*sit + i*numsites]]->set_nindir(verts[map_lat_to_vert[z_off + i*numsites]] ,2 );

            }

        }
    }
    cout<<"created all neighbors "<<n*numsites<<endl;


    switch (measure) {
        case 1: measurement_object = std::make_shared<thermo_int>(measurements, NofD, spins.size() ); break;
        case 2: measurement_object = std::make_shared<switching>(measurements, spins, spins.size()/n ); break;
        case 3: measurement_object = std::make_shared<h_int>(measurements, NofD, spins.size() ); break;
        case 4: measurement_object = std::make_shared<full_energy>(measurements, NofD ); break;
        case 5: measurement_object = std::make_shared<h_int_full>(measurements, NofD, trans_NofD ); break; //can be used for h_plane only
    }

    switch (algo) {
        case 1: {
                    if (exc==3) {
                        for (auto s : spins)
                            s->copy_neighbors_internally(exc) ;
                        update_object = std::make_shared<single_spin_plaq>(seed, n, beta, spins, plaqs, NofD); //spins, plaqs and NofD are referenced
                    }
                    else if (exc==4) {
                        update_object = std::make_shared<single_spin_vert>(seed, n, beta, spins, verts, NofD);
                    }
                }
                break;
        case 2: update_object = std::make_shared<deconfined_vert>(seed, n, beta, spins, verts, NofD);  break;
        case 3: {
                    assert (measure == 3);
                    for (auto s : spins)
                        s->copy_neighbors_internally(exc) ;
                    if (exc==3)
                        update_object = std::make_shared<interaction_metropolis>(seed, n, h, spins, verts, NofD);  //metropolis on vertices = plaquette groundstate
                    else if (exc==4) {
                        update_object = std::make_shared<interaction_metropolis>(seed, n, h, spins, plaqs, NofD);  //metropolis on plaquettes  = vertex groundstate
                    }
                }
                break;
        case 4: {
                    assert (measure == 5);
                    for (auto s : spins)
                        s->copy_neighbors_internally(exc) ;
                    if (exc==3)
                        update_object = std::make_shared<interaction_metropolis_aniso>(seed, n, h, spins, verts, NofD, hz, trans_NofD);  //metropolis on vertices = plaquette gs
                    else if (exc==4) {
                        update_object = std::make_shared<interaction_metropolis_aniso>(seed, n, h, spins, plaqs, NofD, hz, trans_NofD);  //metropolis on plaquettes  = vertex gs
                    }
                }
                break;
        case 5: {
                    assert (measure == 3);
                    for (auto s : spins)
                        s->copy_neighbors_internally(exc) ;
                    if (exc==3)
                        update_object = std::make_shared<interaction_wolff>(seed, n, h, spins, verts, NofD);  //metropolis on vertices = plaquette groundstate
                    else if (exc==4) {
                        update_object = std::make_shared<interaction_wolff>(seed, n, h, spins, plaqs, NofD);  //metropolis on plaquettes  = vertex groundstate
                    }
                }
                break;
        case 6: {
                    assert (measure == 5);
                    for (auto s : spins)
                        s->copy_neighbors_internally(exc) ;
                    if (exc==3)
                        update_object = std::make_shared<interaction_wolff_aniso>(seed, n, h, spins, verts, NofD, hz, trans_NofD);  //metropolis on vertices = plaquette gs
                    else if (exc==4) {
                        update_object = std::make_shared<interaction_wolff_aniso>(seed, n, h, spins, plaqs, NofD, hz, trans_NofD);  //metropolis on plaquettes  = vertex gs
                    }
                }
                break;
    }

    std::cout << "# L: " << L << " Steps: " << Nb_Steps  << " Spins: " <<spins.size()<<" Sites: "<<numsites<< std::endl;
    //print_information();

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

void toriccode::print_information() {
    for (int rep = 0; rep<n; ++rep) {
        for (sit = sites().first; sit!=sites().second; ++sit) {
            cout<<"On the level of ALPS lattice graph ------------------------------------------------------------- "<<endl;
            cout<<"Site "<<*sit<<" in replica "<<rep<<" has type "<<site_type(*sit)<<", and neighbors ";
            for (nit = neighbors(*sit).first; nit!= neighbors(*sit).second; ++nit) {
                cout<<*nit<<"  ";
            }
            cout<<endl;
            cout<<"On the level of the site object  ------------------------------------------------------------- "<<endl;
            if (site_type(*sit) == 6) {
                int siteindex = map_lat_to_vert[*sit + rep*numsites];
                cout<<"VERTEX, maps to No. "<<siteindex<<" which in turn has name "<<verts[siteindex]->get_name()<<endl;
                cout<<"Neighbors have names ";
                for (spit_t spit = verts[siteindex]->get_neighbors_begin(); spit != verts[siteindex]->get_neighbors_end(); ++spit) {
                    cout<<(*spit)->get_name()<<" ";
                }
                cout<<endl;
                cout<<"Counterpart in next replica has name "<<(verts[siteindex]->get_next())->get_name()<<endl;
                cout<<"Move in x direction yields vertex with name "<<(verts[siteindex]->move(0) )->get_name()<<endl;
                cout<<"Move in y direction yields vertex with name "<<(verts[siteindex]->move(1) )->get_name()<<endl;
                cout<<"Move in z direction yields vertex with name "<<(verts[siteindex]->move(2) )->get_name()<<endl;
            }
            else if (site_type(*sit) >= 3) {
                int siteindex = map_lat_to_plaq[*sit + rep*numsites];
                cout<<"PLAQUETTE with orientation "<<plaqs[siteindex]->get_orientation()<<", maps to No. "<<siteindex<<" which in turn has name "<<plaqs[siteindex]->get_name()<<endl;
                cout<<"Neighbors have names ";
                for (spit_t spit = plaqs[siteindex]->get_neighbors_begin(); spit != plaqs[siteindex]->get_neighbors_end(); ++spit) {
                    cout<<(*spit)->get_name()<<" ";
                }
                cout<<endl;
                cout<<"Counterpart in next replica has name "<<(plaqs[siteindex]->get_next())->get_name()<<endl;
            }
            else  {
                int siteindex = map_lat_to_spin[*sit + rep*numsites];
                cout<<"SPIN with orientation "<<spins[siteindex]->get_orientation()<<" with geometry "<<spins[siteindex]->get_geometry()<<", maps to No. "<<siteindex<<" which in turn has name "<<spins[siteindex]->get_name()<<endl;
                cout<<"Neighbors have names ";
                for (const_iit_t iit = spins[siteindex]->get_interaction_neighbors_begin(); iit != spins[siteindex]->get_interaction_neighbors_end(); ++iit) {
                    cout<<(*iit)->get_name()<<" ";
                }
                cout<<endl;
                cout<<"Initial weight is "<<update_object->get_weight_from_spin(spins[siteindex])<<endl;
                cout<<"Counterpart in next replica has name "<<(spins[siteindex]->get_next())->get_name()<<endl;
            }

            cout<<endl;
        }
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

