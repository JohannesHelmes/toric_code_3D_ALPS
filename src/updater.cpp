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
    NofExc = -plaqs.size();
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
    NofExc = -verts.size();
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

using namespace std;

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
    NofExc = -verts.size();

    int counter;
    for (vit_t vit = verts.begin(); vit!=verts.end(); ++vit) {
        counter=0;
        for (const_spit_t spit = (*vit)->get_neighbors_begin(); spit != (*vit)->get_neighbors_end(); ++spit) {
            if ( (*spit)->get_geometry() != 1) {
                (*vit)->add_label( (*spit)->get_geometry() );
                (*vit)->set_boundary ( true );
                ++counter;
            }
            if (counter == 0) {
                (*vit)->add_label( 1 );
                (*vit)->set_boundary ( false );
            }
            else if (counter == 6)
                (*vit)->set_boundary ( false );
        }
    }


    /*
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
    */

    //pattern : option_dict[label1][isboundary1][label2][isboundary2]
    option_dict[1][0][1][0]=1;

    option_dict[0][0][0][0]=2;
    option_dict[2][0][2][0]=2;

    option_dict[0][1][0][0]=3;
    option_dict[0][0][0][1]=3;
    option_dict[2][1][2][0]=3;
    option_dict[2][0][2][1]=3;

    option_dict[0][1][0][1]=4;
    option_dict[2][1][2][1]=4;

    option_dict[1][0][0][0]=5;
    option_dict[1][0][2][0]=5;

    option_dict[0][0][1][0]=15;
    option_dict[2][0][1][0]=15;

    option_dict[1][0][0][1]=6;
    option_dict[1][0][2][1]=6;
    option_dict[0][1][1][0]=6;
    option_dict[2][1][1][0]=6;

    option_dict[0][0][2][0]=7;
    option_dict[2][0][0][0]=7;

    option_dict[0][1][2][0]=8;
    option_dict[2][0][0][1]=8;

    option_dict[2][1][0][0]=9;
    option_dict[0][0][2][1]=9;

    option_dict[0][1][2][1]=10;
    option_dict[2][1][0][1]=10;
    cout<<v.size()<<" vertices, "<<N_verts_per_replica<<" per replica"<<endl;
}


void deconfined_vert::update() {
    for (int j=0; j<=N_verts_per_replica; ++j) {
        replica=random_rep();
        r_vert = random_vert();
        v_cand1=verts[r_vert + replica * N_verts_per_replica];
        v_cand1_cpart=verts[r_vert + (replica^1) * N_verts_per_replica];
        //cout <<(r_vert+ replica * N_verts_per_replica)<<" and its counterpart "<<(r_vert + (replica^1) * N_verts_per_replica)<<endl;
        r_vert = random_vert();
        v_cand2=verts[r_vert + replica * N_verts_per_replica];
        v_cand2_cpart=verts[r_vert + (replica^1) * N_verts_per_replica];

        if (v_cand1 == v_cand2) 
            return;

        label1=v_cand1->get_label();
        label2=v_cand2->get_label();

        switch (option_dict[label1][v_cand1->get_boundary()][label2][v_cand2->get_boundary()]) {
            case  1:
            case  6: try_flip(v_cand1, v_cand2, v_cand1_cpart, v_cand2_cpart, NofExc); break;
            case  2: 
            case  3: try_flip(v_cand1, v_cand2, NofExc); break;
            case  4: try_flip(v_cand1, v_cand2, NofExc);
                     try_flip(v_cand1_cpart, v_cand2_cpart, NofExc); break;
            case  5:
            case  8: v_cand3=verts[random_vert() + (replica^1) * N_verts_per_replica]; 
                     if (v_cand3->get_label() == label2) {
                         try_flip(v_cand1, v_cand1_cpart, v_cand2, v_cand3, NofExc);
                     }
                     break;
            case  7: break;
            case  9:
            case 15: v_cand3=verts[random_vert() + (replica^1) * N_verts_per_replica]; 
                     if (v_cand3->get_label() == label1) {
                         try_flip(v_cand1, v_cand2_cpart, v_cand2, v_cand3, NofExc);
                     }
                     break;
        }
        
    }
    
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
