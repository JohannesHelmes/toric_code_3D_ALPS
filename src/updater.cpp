#include "updater.h"
#include<cassert>

/********** base class updater **************/
updater::updater(int seed, int reps, double beta, std::vector<spin_ptr>& s) : 
        spins(s),
        N(spins.size()),
        mtwister(seed), //change this!!!
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

single_spin_plaq::single_spin_plaq(int seed, int reps, double beta, std::vector<spin_ptr>& s, std::vector<plaq_ptr>& p, int& nofe) : 
        updater(seed, reps, beta, s),
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

single_spin_vert::single_spin_vert(int seed, int reps, double beta, std::vector<spin_ptr>& s, std::vector<vert_ptr>& v, int& nofe) : 
        updater(seed, reps, beta, s),
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


mix_spin_plaq_for_vert::mix_spin_plaq_for_vert(int seed, int reps, double beta, std::vector<spin_ptr>& s, std::vector<plaq_ptr>& p, std::vector<vert_ptr>& v, int& nofe, double ratio) : 
        updater(seed, reps, beta, s),
        Nspinflips((int)(ratio*N)),
        plaqs(p),
        verts(v),
        NofExc(nofe)
{
    std::cout<<"vertex MIX updater created, do "<<Nspinflips<<" single spin flips and "<<N-Nspinflips<<" plaquette flips"<<std::endl;
    std::cout<<"plaqs size "<<p.size()<<" and spin size "<<s.size()<<std::endl;
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

deconfined_vert::deconfined_vert(int seed, int reps, double beta, std::vector<spin_ptr>& s, std::vector<vert_ptr>& v, int& nofe) :
        updater(seed, reps, beta, s),
        verts(v),
        N_verts_per_replica(v.size()/reps),
        NofExc(nofe),
        int_dist_reps(0,reps-1),
        int_dist_verts(0,N_verts_per_replica-1),
        random_rep(mtwister, int_dist_reps),
        random_vert(mtwister, int_dist_verts)
{
}

using namespace std;

void deconfined_vert::update() {
    replica=random_rep();
    cout<<"Picked replica "<<replica<<endl;
    r_vert = random_vert();
    v_cand1=verts[r_vert + replica * N_verts_per_replica];
    v_cand1_cpart=verts[r_vert + replica^1 * N_verts_per_replica];
    r_vert = random_vert();
    v_cand2=verts[r_vert + replica * N_verts_per_replica];
    v_cand2_cpart=verts[r_vert + replica^1 * N_verts_per_replica];

    if (v_cand1 == v_cand2)
        return;

    label1 = v_cand1->get_label();
    label2 = v_cand2->get_label();

    if (label2==1) {
        swap(v_cand1, v_cand2);
        swap(v_cand1_cpart, v_cand2_cpart);
    }
    else if (label1!=1) {
        if (label1 > label2) {
            swap(v_cand1, v_cand2);
            swap(v_cand1_cpart, v_cand2_cpart);
        }
        else if ((label1 == label2)&& (v_cand2->get_boundary())) {
            swap(v_cand1, v_cand2);
            swap(v_cand1_cpart, v_cand2_cpart);
        }
    }



    if (label1 == 1) {
        if ((label2 == 1)||(v_cand2->get_boundary())) {
            //flip both in A
            try_flip(v_cand1, v_cand2, v_cand1_cpart, v_cand2_cpart, NofExc);
        }
        else {
            v_cand3=verts[random_vert() + replica^1 * N_verts_per_replica]; 
            if (v_cand3->get_label() == label2) {
                try_flip(v_cand1, v_cand1_cpart, v_cand2, v_cand3, NofExc);
            }
        }
    }
    else if (label1 == label2) {
        if (v_cand2->get_boundary() ) { //implies v_cand1->get_boundary()==true
            assert(v_cand1->get_boundary()==true);
            try_flip(v_cand1, v_cand2, NofExc);
            try_flip(v_cand1_cpart, v_cand2_cpart, NofExc);
            //try both
        }
        else { //do update in B
            try_flip(v_cand1, v_cand2, NofExc);
        }
    }
    else {
        assert(label1==0);
        assert(label2==2);
        if ((v_cand1->get_boundary())&& (v_cand2->get_boundary())) {
            try_flip(v_cand1, v_cand2, v_cand1_cpart, v_cand2_cpart, NofExc);
        }
        else if ((v_cand1->get_boundary())^(v_cand2->get_boundary())) {
            v_cand3=verts[random_vert() + replica^1 * N_verts_per_replica]; 
            if ((v_cand1->get_boundary()) && (v_cand3->get_label() == label2)) {
                try_flip(v_cand1, v_cand1_cpart, v_cand2, v_cand3, NofExc);
            }
            if ((v_cand2->get_boundary()) && (v_cand3->get_label() == label1)) {
                try_flip(v_cand1, v_cand2, v_cand2_cpart, v_cand3, NofExc);
            }
        }
        else 
            std::cout<<"NO way "<<std::endl;
    }

    //otherwise, both are in B but disconnected -> abort
    
}

void deconfined_vert::swap(vert_ptr& v1, vert_ptr& v2) {
    tmp = v1;
    v1 = v2;
    v2 = tmp;
}

void deconfined_vert::try_flip(vert_ptr& v1, vert_ptr& v2, int& NofD) {
    cand_weight = v1->get_value() + v2->get_value() ;
    if ((cand_weight>=0)||(random_01()<expmB[-2*cand_weight])) {
        NofD -= 2*cand_weight; 
        v1->flip();
        v2->flip();
    }
}

void deconfined_vert::try_flip(vert_ptr& v1, vert_ptr& v2, vert_ptr& v3, vert_ptr& v4, int& NofD) {
    cand_weight = v1->get_value() + v2->get_value() + v3->get_value() + v4->get_value();
    if ((cand_weight>=0)||(random_01()<expmB[-2*cand_weight])) {
        NofD -= 2*cand_weight; 
        v1->flip();
        v2->flip();
        v3->flip();
        v4->flip();
    }
}
