#include "site.h"
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
    //std::unordered_map<vertex_pair_t, int> 

    boost::random::uniform_int_distribution<> int_dist;
    boost::variate_generator<mt_rng&, boost::random::uniform_int_distribution<> > random_int;
};

class single_spin_plaq : public updater {
public:
    single_spin_plaq(int seed, int reps, double beta, std::vector<spin_ptr>& s, std::vector<plaq_ptr>& p, int& nofe);  //single spin update for plaquette Hamiltonian
    void update();
private:
    int cand_weight;
    std::vector<plaq_ptr> plaqs;
    spin_ptr candidate;
    int &NofExc;
};
       
class single_spin_vert : public updater {
public:
    single_spin_vert(int seed, int reps, double beta, std::vector<spin_ptr>& s, std::vector<vert_ptr>& v, int& nofe);  //single spin update for vertex Hamiltonian
    void update();
private:
    int cand_weight;
    std::vector<vert_ptr> verts;
    spin_ptr candidate;
    int &NofExc;
};

class mix_spin_plaq_for_vert : public updater {
public:
    mix_spin_plaq_for_vert(int seed, int reps, double beta, std::vector<spin_ptr>& s, std::vector<plaq_ptr>& p, std::vector<vert_ptr>& v, int& nofe, double ratio);  //single spin update for vertex Hamiltonian
    void update();
private:
    int Nspinflips;
    int cand_weight;
    std::vector<vert_ptr> verts;
    std::vector<plaq_ptr> plaqs;
    spin_ptr candidate;
    int &NofExc;
};


class deconfined_vert : public updater {
public:
    deconfined_vert(int seed, int reps, double beta, std::vector<spin_ptr>& s, std::vector<vert_ptr>& v, int& nofe);  //single spin update for vertex Hamiltonian
    void update();
private:
    std::vector<vert_ptr> verts;
    const int N_verts_per_replica;
    int replica;
    boost::random::uniform_int_distribution<> int_dist_reps;
    boost::random::uniform_int_distribution<> int_dist_verts;
    boost::variate_generator<mt_rng&, boost::random::uniform_int_distribution<> > random_rep;
    boost::variate_generator<mt_rng&, boost::random::uniform_int_distribution<> > random_vert;
    int &NofExc;

    vert_ptr v_cand1, v_cand2, v_cand1_cpart, v_cand2_cpart, v_cand3,tmp;
    int r_vert, itmp;
    int label1,label2,label3;
    int cand_weight;

    void swap(vert_ptr& v1, vert_ptr& v2);
    void swap(int& i1, int& i2);
    void try_flip(vert_ptr& v1, vert_ptr& v2, int& NofD);
    void try_flip(vert_ptr& v1, vert_ptr& v2, vert_ptr& v3, vert_ptr& v4, int& NofD);
};
