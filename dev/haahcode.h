#ifndef MYBASEMC_HPP
#define MYBASEMC_HPP
#include <alps/scheduler.h>
#include <string>
#include "updater.h" //includes site.h
#include "measurer.h" 


typedef alps::scheduler::LatticeMCRun<>::graph_type graph_type;

using namespace std;

class haahcode : public alps::scheduler::LatticeMCRun<graph_type>{

public :
    haahcode(const alps::ProcessList& where,const alps::Parameters& p,int node);
    static void print_copyright(std::ostream &);
    void save(alps::ODump& dump) const;
    void load(alps::IDump& dump);
    void dostep();
    bool is_thermalized() const;
    double work_done() const;
protected :
    alps::uint64_t Nb_Steps;
    alps::uint64_t Nb_Therm_Steps;
    alps::uint64_t Total_Steps;
    int L,N,numsites,NofD,exc,numspins,seed;
    double beta,ratio,h,hz;
    string IncStep;
    std::vector<int> geom;
    std::vector<int> map_lat_to_spin, map_lat_to_plaq, map_lat_to_vert;
    std::vector<spin_ptr> spins;
    std::vector<inter_ptr> plaqs;
    std::vector<inter_ptr> verts;
    spin_ptr candidate;
    std::shared_ptr<updater> update_object;
    std::shared_ptr<measurer> measurement_object;


    site_iterator sit;
    neighbor_iterator nit;


    void print_information() ;

};

typedef alps::scheduler::SimpleMCFactory<haahcode> HaahFactory;
    
#endif
