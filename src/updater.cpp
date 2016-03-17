#include "updater.h"

/********** base class updater **************/
updater::updater(int reps, double beta, std::vector<spin_ptr>& s) : 
        spins(s),
        N(spins.size()),
        mtwister(42), //change this!!!
        int_dist(0,N-1),
        real_dist(),
        random_int(mtwister, int_dist),
        random_01(mtwister, real_dist)
{
    expmB.resize(8*reps+1);
    expmB[0]=1.0;
    for (int i=1; i<=8*reps; ++i) 
        expmB[i]=std::exp(-1.*i*beta);
}
    
/********** class single_spin **************/

single_spin_plaq::single_spin_plaq(int reps, double beta, std::vector<spin_ptr>& s, std::vector<plaq_ptr>& p, int& nofe) : 
        updater(reps, beta, s),
        plaqs(p),
        NofExc(nofe)
{
}

void single_spin_plaq::update() {

    for (int j=0; j<N; ++j) {
        candidate=spins[random_int()];
        cand_weight = candidate->get_weight_from_plaqs(); 

        if ((cand_weight>=0)||(random_01()<expmB[-2*cand_weight])) {
            NofExc -= 2*cand_weight; //is the old weight
            candidate->flip_and_flip_plaqs();
        }
    }
}

single_spin_vert::single_spin_vert(int reps, double beta, std::vector<spin_ptr>& s, std::vector<vert_ptr>& v, int& nofe) : 
        updater(reps, beta, s),
        verts(v),
        NofExc(nofe)
{
}

void single_spin_vert::update() {

    for (int j=0; j<N; ++j) { //generalize this !!!

        candidate=spins[random_int()];
        cand_weight = candidate->get_weight_from_verts();

        if ((cand_weight>=0)||(random_01()<expmB[-2*cand_weight])) {
            NofExc -= 2*cand_weight; //is the old weight
            candidate->flip_and_flip_verts();
        }
    }
}


mix_spin_plaq_for_vert::mix_spin_plaq_for_vert(int reps, double beta, std::vector<spin_ptr>& s, std::vector<plaq_ptr>& p, std::vector<vert_ptr>& v, int& nofe, double ratio) : 
        updater(reps, beta, s),
        Nspinflips((int)(ratio*N)),
        plaqs(p),
        verts(v),
        NofExc(nofe)
{
    std::cout<<"vertex MIX updater created, do "<<Nspinflips<<" single spin flips and "<<N-Nspinflips<<" plaquette flips"<<std::endl;
}

void mix_spin_plaq_for_vert::update() {

    for (int j=0; j<Nspinflips; ++j) { //generalize this !!!

        candidate=spins[random_int()];
        cand_weight = candidate->get_weight_from_verts();

        if ((cand_weight>=0)||(random_01()<expmB[-2*cand_weight])) {
            NofExc -= 2*cand_weight; //is the old weight
            candidate->flip_and_flip_verts();
        }
    }

    for (int j=Nspinflips; j<N; ++j) { //generalize this !!!
        plaqs[random_int()]->flip_neighbors();
    }
}

