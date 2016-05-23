#include "updater.h"
#include<cassert>

using namespace std;



/********** base class updater **************/
updater::updater(int seed, int reps, double beta, std::vector<spin_ptr>& s) : 
        spins(s),
        beta(beta),
        N(spins.size()),
        mtwister(seed), //change this!!!
        int_dist(0,N-1),
        real_dist(),
        random_int(mtwister, int_dist),
        random_01(mtwister, real_dist)
{
    expmB.resize(24*reps+1);
    expmB[0]=1.0;
    for (int i=1; i<=24*reps; ++i) 
        expmB[i]=std::exp(-1.*i*beta);
}
    
/********** class single_spin_plaq **************/

single_spin_plaq::single_spin_plaq(int seed, int reps, double beta, std::vector<spin_ptr>& s, std::vector<inter_ptr>& p, int& nofe) : 
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

    //do_winding_update();

    loop_weight = 0;
    loop_set.clear();
    first_spin = spins[random_int()];
    orient = (random_int()%2 + first_spin->get_orientation() + 1) %3;

    ss_plaq_fill_loop(first_spin, loop_set, loop_weight);
    //cout<<"Loop has "<<loop_set.size()<<" elements and weight "<<loop_weight<<" and orientation "<<first_spin->get_orientation()<<", "<<orient<<" and geom "<<first_spin->get_geometry()<<endl;

    if ((loop_weight>=0)||(random_01()<exp(2*beta*loop_weight))) {
        for (auto l : loop_set)
            l->flip_and_flip_plaqs();
        NofExc -= 2*loop_weight;
    }

}

void single_spin_plaq::ss_plaq_fill_loop(spin_ptr spin, std::unordered_set<spin_ptr, std::my_hash>& l_set, int& weight) {
    if (l_set.find(spin) == l_set.end() ) {

        l_set.insert(spin);
        bool weightplaq;

        for (const_iit_t loop_iit = spin->get_interaction_neighbors_begin(); loop_iit != spin->get_interaction_neighbors_end(); ++loop_iit) {
            weightplaq = false;
            for (const_spit_t loop_sit = (*loop_iit)->get_neighbors_begin(); loop_sit != (*loop_iit)->get_neighbors_end(); ++loop_sit ) {
                if ( (*loop_sit)->get_orientation() ==orient) {
                    weight += (*loop_iit)->get_value(); //value determined from plaquettes
                    weightplaq = true;
                    //cout<<"Plaquette "<<*loop_iit<<" considered for weight "<<weight<<endl;
                    break;
                }
            }
            if (weightplaq)
                continue;

            //cout<<"Plaquette "<<*loop_iit<<" considered for propagation "<<endl;
            for (const_spit_t loop_sit = (*loop_iit)->get_neighbors_begin(); loop_sit != (*loop_iit)->get_neighbors_end(); ++loop_sit ) {
                if ( (*loop_sit != spin)&&( (*loop_sit)->get_orientation() == spin->get_orientation() ) )  {
                    ss_plaq_fill_loop(*loop_sit, l_set, weight);
                }
            }
        }
    }
}

/********** class single_spin_vert **************/

single_spin_vert::single_spin_vert(int seed, int reps, double beta, std::vector<spin_ptr>& s, std::vector<inter_ptr>& v, int& nofe) : 
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


/********** class deconfined_vert  **************/

deconfined_vert::deconfined_vert(int seed, int reps, double beta, std::vector<spin_ptr>& s, std::vector<inter_ptr>& v, int& nofe) :
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

    //label all connected regions and boundaries of vertices/plaquettes
    int counter;
    for (iit_t vit = verts.begin(); vit!=verts.end(); ++vit) {
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

void deconfined_vert::try_flip(inter_ptr& v1, inter_ptr& v2, int& NofD) {
    cand_weight = v1->get_value() + v2->get_value() ;
    if ((cand_weight>=0)||(random_01()<expmB[-2*cand_weight])) {
        NofD -= 2*cand_weight; 
        v1->flip();
        v2->flip();
    }
}

void deconfined_vert::try_flip(inter_ptr& v1, inter_ptr& v2, inter_ptr& v3, inter_ptr& v4, int& NofD) {
    cand_weight = v1->get_value() + v2->get_value() + v3->get_value() + v4->get_value();
    if ((cand_weight>=0)||(random_01()<expmB[-2*cand_weight])) {
        NofD -= 2*cand_weight; 
        v1->flip();
        v2->flip();
        v3->flip();
        v4->flip();
    }
}


/********** abstract class winding_updater  **************/

winding_updater::winding_updater(int seed, int reps, double beta, double h, std::vector<spin_ptr>& s, int& total_observable) :
        updater(seed, reps, beta, s),
        h(h),
        TObs(total_observable),
        dual(dual)
{
    TObs = -spins.size();
    cout<<"Winding updater initialized, h = "<<h<<endl;
}

void winding_updater::do_winding_update() {
    loop_weight = 0;
    loop_set.clear();
    first_spin = spins[random_int()];


    fill_loop(first_spin, loop_set, loop_weight);
    //cout<<"Loop has "<<loop_set.size()<<" elements and weight "<<loop_weight<<" and orientation "<<first_spin->get_orientation()<<endl;

    if ((loop_weight>=0)||(random_01()<exp(2*h*loop_weight))) {
        for (auto l : loop_set)
            l->flip();
        TObs -= 2*loop_weight;
    }
}

void winding_updater::fill_loop(spin_ptr spin, std::unordered_set<spin_ptr, std::my_hash>& l_set, int& weight) {
    if (l_set.find(spin) == l_set.end() ) {
        l_set.insert(spin);
        weight += spin->get_value();
        if (spin->get_geometry()==1)
            weight += spin->get_value();

        for (const_iit_t loop_iit = spin->get_interaction_neighbors_begin(); loop_iit != spin->get_interaction_neighbors_end(); ++loop_iit) {
            for (const_spit_t loop_sit = (*loop_iit)->get_neighbors_begin(); loop_sit != (*loop_iit)->get_neighbors_end(); ++loop_sit ) {
                if ( (*loop_sit != spin)&&( (*loop_sit)->get_orientation() == spin->get_orientation() ) )  {
                    fill_loop(*loop_sit, l_set, weight);
                }
            }
        }
    }
}


/*******************************************************/
/*                                                     */
/*          class interaction_metropolis               */
/*                                                     */
/*******************************************************/

interaction_metropolis::interaction_metropolis(int seed, int reps, double h, std::vector<spin_ptr>& s, std::vector<inter_ptr>& ia, int& total_magn) :
        winding_updater(seed, reps, h, h, s, total_magn),
        interactions(ia),
        N_interactions(ia.size()),
        TMagn(total_magn),
        int_dist_interactions(0,N_interactions-1),
        random_interaction(mtwister, int_dist_interactions)
{
    TMagn = -spins.size();

    //label all connected regions and boundaries of vertices/plaquettes
    int A_counter, B_counter;
    for (iit_t iit = interactions.begin(); iit!=interactions.end(); ++iit) {
        B_counter=0;
        A_counter=0;
        for (const_spit_t spit = (*iit)->get_neighbors_begin(); spit != (*iit)->get_neighbors_end(); ++spit) {
            if ( (*spit)->get_geometry() != 1) {
                (*iit)->add_label( (*spit)->get_geometry() );
                (*iit)->set_boundary ( true );
                ++B_counter;
            }
            else
                ++A_counter;
        }
        if (B_counter == 0) {
            (*iit)->add_label( 1 );
            (*iit)->set_boundary ( false );
        }
        else if (A_counter == 0) {
            (*iit)->set_boundary ( false );
        }
    }
    cout<<"Vertex metropolis initialized, vertices "<<N_interactions<<", magnetization "<<h<<endl;

}

void interaction_metropolis::update() {
    for (int j=0; j<N_interactions; ++j) {
        cand = interactions[random_interaction()];
        nb_spin_it= cand->get_neighbors_begin();
        cand_weight = 0;
        if ((cand->get_boundary()== true)||((*nb_spin_it)->get_geometry()==1)) { //the second option is for the case, that cand is completely in A
            //in other replicas, too
            runner = cand;
            do {
                for (nb_spin_it = runner->get_neighbors_begin(); nb_spin_it!=runner->get_neighbors_end(); ++nb_spin_it) {
                    cand_weight += (*nb_spin_it)->get_value();
                }
                runner = runner->get_next();
            } while (runner != cand);
            
            if ((cand_weight>=0)||(random_01()<expmB[-2*cand_weight])) {
                cand->flip_neighbors();
                runner = cand->get_next();
                while (runner != cand) { 
                    for (nb_spin_it = runner->get_neighbors_begin(); nb_spin_it!=runner->get_neighbors_end(); ++nb_spin_it) {
                        if ((*nb_spin_it)->get_geometry() != 1) {
                            assert(cand->get_boundary() == true);
                            (*nb_spin_it)->flip();
                        }
                    }
                    runner = runner->get_next();
                }
                TMagn -= 2*cand_weight; 
            }

        }
        else {
            for (nb_spin_it; nb_spin_it!=cand->get_neighbors_end(); ++nb_spin_it) {
                cand_weight += (*nb_spin_it)->get_value();
            }

            if ((cand_weight>=0)||(random_01()<expmB[-2*cand_weight])) {
                cand->flip_neighbors();
                TMagn -= 2*cand_weight; 
            }
        }

        /*
        int counted_magn=0;
        for (spit_t spit = spins.begin(); spit!=spins.end(); ++spit)
            counted_magn+= (*spit)->get_value();
        if (TMagn != counted_magn)
            cout<<"Error "<<TMagn<<" vs, "<<counted_magn<<endl;
            */
    }

    do_winding_update();


    // Try a single winding loop update
    //loop_weight = 0;
    //loop_set.clear();
    //first_spin = spins[random_int()];


    //fill_loop(first_spin, loop_set, loop_weight);
    ////cout<<"Loop has "<<loop_set.size()<<" elements and weight "<<loop_weight<<" and orientation "<<first_spin->get_orientation()<<endl;

    //if ((loop_weight>=0)||(random_01()<exp(2*beta*loop_weight))) {
    //    for (auto l : loop_set)
    //        l->flip();
    //    TMagn -= 2*loop_weight;
    //}
}


/*******************************************************/
/*                                                     */
/*          class interaction_wolff                    */
/*                                                     */
/*******************************************************/

interaction_wolff::interaction_wolff(int seed, int reps, double h, std::vector<spin_ptr>& s, std::vector<inter_ptr>& ia, int& total_magn) :
        winding_updater(seed, reps, h, h, s, total_magn),
        interactions(ia),
        N_interactions(ia.size()),
        TMagn(total_magn),
        int_dist_interactions(0,N_interactions-1),
        random_interaction(mtwister, int_dist_interactions)
{
    TMagn = -spins.size();

    //label all connected regions and boundaries of vertices/plaquettes
    int A_counter, B_counter;
    for (iit_t iit = interactions.begin(); iit!=interactions.end(); ++iit) {
        B_counter=0;
        A_counter=0;
        for (const_spit_t spit = (*iit)->get_neighbors_begin(); spit != (*iit)->get_neighbors_end(); ++spit) {
            if ( (*spit)->get_geometry() != 1) {
                (*iit)->add_label( (*spit)->get_geometry() );
                (*iit)->set_boundary ( true );
                ++B_counter;
            }
            else
                ++A_counter;
        }
        if (B_counter == 0) {
            (*iit)->add_label( 1 );
            (*iit)->set_boundary ( false );
        }
        else if (A_counter == 0) {
            (*iit)->set_boundary ( false );
        }
        //else {
         //   cout<<*iit<<" is boundary, Acounter= "<<A_counter<<endl;
        //}
    }
    cout<<"Wolff initialized, magnetization "<<h<<endl;
}


void interaction_wolff::flip_adjacents(inter_ptr the_inter) {
    the_inter->flip_neighbors();
    if (the_inter->get_boundary() ) {
        inter_ptr other = the_inter->get_next();
        //cout<<"  In Flipping boundary case "<<endl;
        for (const_spit_t spit = other->get_neighbors_begin(); spit != other->get_neighbors_end(); ++spit) {
            if ( (*spit)->get_geometry() != 1) {
                (*spit)->flip();
            }
        }
    }
}

void interaction_wolff::update() {

    for (int n=0; n<12; ++n) {
        cluster_members.clear();
        visited_vertices.clear();

        inter_ptr start_inter = interactions[random_interaction()];
        cluster_members.push_back(start_inter );
        flip_adjacents(start_inter);
        auto cluster_iterator = cluster_members.begin();

        //do the Wolff in while loop
        while ( cluster_iterator != cluster_members.end() ) {
            at_inter = *(cluster_iterator);

            if (visited_vertices.find(at_inter) != visited_vertices.end() ) {
                //cout<<" Failed !!!, break "<<endl;
                for (auto cm: cluster_members)
                    flip_adjacents(cm);
                break;
            }
            visited_vertices.insert(at_inter);


            //do the recursion
            for (spit = at_inter->get_neighbors_begin(); spit != at_inter->get_neighbors_end(); ++spit) {
                weight = (*spit)->get_geometry() == 1 ? 2*(*spit)->get_value() : (*spit)->get_value();
                if ( weight > 0 ) {
                    next_iter = (*spit)->get_dual_interaction_neighbors_begin();
                    while ( ( (*next_iter) == at_inter) || ( (*next_iter) == at_inter->get_next() ) ) {
                        ++next_iter;
                    }

                    //next_iter may or may not be the neighbor in the other replica 
                    if (random_01() < (1 - expmB[2*weight] ) ) {
                        cluster_members.push_back(*next_iter) ;
                        flip_adjacents(*next_iter);
                    }
                }
            }

            if (at_inter->get_boundary() ) {
                other = at_inter->get_next();
                for (spit = other->get_neighbors_begin(); spit != other->get_neighbors_end(); ++spit) {
                    if ( (*spit)->get_geometry() != 1) {
                        weight = (*spit)->get_value();
                        if ( weight > 0 ) {
                            next_iter = (*spit)->get_dual_interaction_neighbors_begin(); //plaquettes groundstate (=interaction) => vertices = dual interaction
                            if ( ( (*next_iter) == other)  )
                                ++next_iter;
                            if ( (random_01() < (1 - expmB[2*weight] ) ) )  {
                                cluster_members.push_back(*next_iter);
                                flip_adjacents(*next_iter);
                            }
                        }
                    }
                }
            }
            ++cluster_iterator;
        }
    }

    TMagn = 0;
    for (auto s : spins) {
        TMagn += s->get_value();
    }
    
    do_winding_update();
}
