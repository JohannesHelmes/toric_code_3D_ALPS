#include "site.h"
#include <boost/random/mersenne_twister.hpp> // for random_int, random_01
#include <boost/random/uniform_01.hpp>
#include <boost/random/uniform_int.hpp>
#include <boost/random/variate_generator.hpp>

typedef boost::random::mt19937 mt_rng;

class updater {
public:
    updater(int reps, double beta, std::vector<spin_ptr>& s);
    virtual void update()=0;
protected:
    //include boost rng here!

    mt_rng mtwister;
    boost::uniform_int<> int_dist;
    boost::uniform_01<> real_dist;
    boost::variate_generator<mt_rng, boost::uniform_int<> > random_int;
    boost::variate_generator<mt_rng&, boost::uniform_01<> > random_01;
    std::vector<spin_ptr>* spins;
    std::vector<double> expmB;
    int N;
};

class single_spin_plaq : updater {
public:
    single_spin_plaq(int reps, double beta, std::vector<spin_ptr>& s, std::vector<plaq_ptr>& p, int& nofe);  //single spin update for plaquette Hamiltonian
    void update();
private:
    int cand_weight;
    std::vector<plaq_ptr>* plaqs;
    spin_ptr candidate;
    int *NofExc;
};
       
class single_spin_vert : updater {
public:
    single_spin_vert(int reps, double beta, std::vector<spin_ptr>& s, std::vector<vert_ptr>& v, int& nofe);  //single spin update for vertex Hamiltonian
    void update();
private:
    int cand_weight;
    std::vector<vert_ptr>* verts;
    spin_ptr candidate;
    int *NofExc;
};
