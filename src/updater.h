#include "site.h"
#include <boost/multi_array.hpp>
#include <boost/random/mersenne_twister.hpp> // for random_int, random_01
#include <boost/random/uniform_01.hpp>
#include <boost/random/uniform_int.hpp>
#include <boost/random/variate_generator.hpp>

typedef boost::random::mt19937 mt_rng;

class updater {
public:
    updater(int seed, int reps, double beta, std::vector<spin_ptr>& s);
    virtual void update()=0;
protected:
    
    //The order of these declaration is crucial because of dependencies in the initialization list!
    std::vector<spin_ptr> spins;
    const int N;
    mt_rng mtwister;
    boost::uniform_01<> real_dist;
    boost::variate_generator<mt_rng&, boost::uniform_01<> > random_01;
    std::vector<double> expmB;
    boost::multi_array<short, 4> option_dict{boost::extents[3][2][3][2]};

    boost::random::uniform_int_distribution<> int_dist;
    boost::variate_generator<mt_rng&, boost::random::uniform_int_distribution<> > random_int;
};

class single_spin_plaq : public updater {
public:
    single_spin_plaq(int seed, int reps, double beta, std::vector<spin_ptr>& s, std::vector<inter_ptr>& p, int& nofe);  //single spin update for plaquette Hamiltonian
    void update();
private:
    int cand_weight;
    std::vector<inter_ptr> plaqs;
    spin_ptr candidate;
    int &NofExc;
};
       
class single_spin_vert : public updater {
public:
    single_spin_vert(int seed, int reps, double beta, std::vector<spin_ptr>& s, std::vector<inter_ptr>& v, int& nofe);  //single spin update for vertex Hamiltonian
    void update();
private:
    int cand_weight;
    std::vector<inter_ptr> verts;
    spin_ptr candidate;
    int &NofExc;
};


class deconfined_vert : public updater {
public:
    deconfined_vert(int seed, int reps, double beta, std::vector<spin_ptr>& s, std::vector<inter_ptr>& v, int& nofe);  //single spin update for vertex Hamiltonian
    void update();
private:
    std::vector<inter_ptr> verts;
    const int N_verts_per_replica;
    int replica;
    boost::random::uniform_int_distribution<> int_dist_reps;
    boost::random::uniform_int_distribution<> int_dist_verts;
    boost::variate_generator<mt_rng&, boost::random::uniform_int_distribution<> > random_rep;
    boost::variate_generator<mt_rng&, boost::random::uniform_int_distribution<> > random_vert;
    int &NofExc;

    inter_ptr v_cand1, v_cand2, v_cand1_cpart, v_cand2_cpart, v_cand3;
    int r_vert;
    int label1,label2,label3;
    int cand_weight;

    void try_flip(inter_ptr& v1, inter_ptr& v2, int& NofD);
    void try_flip(inter_ptr& v1, inter_ptr& v2, inter_ptr& v3, inter_ptr& v4, int& NofD);
};


class interaction_metropolis : public updater {
public:
    interaction_metropolis(int seed, int reps, double h, std::vector<spin_ptr>& s, std::vector<inter_ptr>& v, int& total_magn);  //single vertex update for plaquette Hamiltonian
    void update();
private:
    std::vector<inter_ptr> verts;
    int const N_verts;
    boost::random::uniform_int_distribution<> int_dist_verts;
    boost::variate_generator<mt_rng&, boost::random::uniform_int_distribution<> > random_vert;
    int &TMagn;

    int cand_weight;
    inter_ptr cand, runner;
    spit_t nb_spin_it;
};
