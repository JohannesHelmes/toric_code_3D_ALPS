#ifndef MYBASEMC_HPP
#define MYBASEMC_HPP
#include <boost/random.hpp>
#include <alps/scheduler.h>
#include <alps/lattice/graph_helper.h>
#include <alps/alea.h>
#include <alps/alea/histogram.h>
#include <alps/scheduler/montecarlo.h>
#include <tuple>
#include <string>



// typedef boost::adjacency_list<boost::vecS,boost::vecS,boost::undirectedS, boost::property<alps::vertex_type_t,unsigned int>, boost::property<alps::edge_type_t,unsigned int, boost::property<alps::edge_index_t, unsigned int> > > graph_type; //doesn't work !!!
typedef alps::scheduler::LatticeMCRun<>::graph_type graph_type;

using namespace std;

struct defect {
    int vertex;
    int spin;

    defect() {vertex=-5; spin=-5;}
    defect(int a, int b) {vertex=a; spin=b;}
};

class toricFTSP : public alps::scheduler::LatticeMCRun<graph_type>{

public :
    toricFTSP(const alps::ProcessList& where,const alps::Parameters& p,int node);
    static void print_copyright(std::ostream &);
    void save(alps::ODump& dump) const;
    void load(alps::IDump& dump);
    void dostep();
    bool is_thermalized() const;
    double work_done() const;
protected :
    alps::uint64_t const Nb_Steps;
    alps::uint64_t const Nb_Therm_Steps;
    alps::uint64_t Total_Steps;
    int W,N,numspins,numsites,start,B,replica,NofD,ediff;
    int const L,n,d,exc;
    bool WasInA, InA, found;
    string IncStep;
    std::vector<bool> geom,edge;
    std::set<site_descriptor> edge_sites;
    std::vector<bool> spins_;
    std::vector<std::vector<bool> > vertex_defects;
    std::vector<int> vneighs;

    void heal_chain(int replica, std::deque<int> chain, std::deque<int> back_chain);
    void heal_chain_open(int replica, std::deque<int> chain, std::deque<int> back_chain);
    void heal_loop(int replica, std::deque<int> chain);
    int get_random_neighbor(int vertex);
    void flip_defect(int replica, int vertex);

    double beta,weight;
    std::vector<double> expmB;
    site_iterator sit;
    neighbor_iterator nit;
    std::deque<int>::iterator dit;
    std::vector<int>::iterator vit;
    std::vector<bool>::iterator vbit;

    void flip(int replica, int spin);

    virtual void do_measurements()=0;
    bool IsInA(site_descriptor) const;

    template<class T> typename deque<T>::reverse_iterator back_find(std::deque<T> *d,T item) {
        typename deque<T>::reverse_iterator rdit;
        rdit= std::find((*d).rbegin(),(*d).rend(),item);
        return rdit;
    };
    template<class T> typename deque<T>::iterator find(std::deque<T> *d,T item) {
        return std::find((*d).begin(),(*d).end(),item);
    };


};
#endif
