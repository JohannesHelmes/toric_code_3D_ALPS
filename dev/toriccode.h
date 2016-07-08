#ifndef MYBASEMC_HPP
#define MYBASEMC_HPP
#include <alps/scheduler.h>
#include <string>
#include "updater.h" //includes site.h
#include "measurer.h" 


typedef alps::scheduler::LatticeMCRun<>::graph_type graph_type;

using namespace std;

class toriccode : public alps::scheduler::LatticeMCRun<graph_type>{

public :
    toriccode(const alps::ProcessList& where,const alps::Parameters& p,int node);
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
    int L,W,N,numsites,start,n,NofD,exc,numspins,algo,seed,measure,trans_NofD;
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

    void heal_chain(int replica, std::deque<int> chain, std::deque<int> back_chain);
    void heal_chain_open(int replica, std::deque<int> chain, std::deque<int> back_chain);
    void heal_loop(int replica, std::deque<int> chain);
    int get_random_neighbor(int vertex);
    void flip_defect(int replica, int plaquette);

    site_iterator sit;
    neighbor_iterator nit;

    bool IsInA(site_descriptor);
    template<class T> typename deque<T>::reverse_iterator back_find(std::deque<T> *d,T item) {
        typename deque<T>::reverse_iterator rdit;
        rdit= std::find((*d).rbegin(),(*d).rend(),item);
        return rdit;
    };
    template<class T> typename deque<T>::iterator find(std::deque<T> *d,T item) {
        return std::find((*d).begin(),(*d).end(),item);
    };

    void print_information() ;

};

typedef alps::scheduler::SimpleMCFactory<toriccode> ToricFactory;
    
#endif
