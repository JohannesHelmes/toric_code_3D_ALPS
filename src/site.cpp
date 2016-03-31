#include "site.h"

/************* NON Members ********************************/
alps::ODump& operator<<(alps::ODump& dump, const std::vector<spin_ptr>& sp) {
    for (const_spit_t spit=sp.begin(); spit!=sp.end(); ++spit) {
        dump<<(*spit)->get_value();
    }
    return dump;
}

alps::IDump& operator>>(alps::IDump& dump, std::vector<spin_ptr>& sp) {
    int tmp;
    for (spit_t spit=sp.begin(); spit!=sp.end(); ++spit)  {
        dump>>tmp;
        (*spit)->set_value( (tmp==1)? true: false);
    }
    return dump;
}

alps::ODump& operator<<(alps::ODump& dump, const std::vector<plaq_ptr>& pl) {
    for (const_plit_t pit=pl.begin(); pit!=pl.end(); ++pit) 
        dump<<(*pit)->get_value();
    return dump;
}

alps::IDump& operator>>(alps::IDump& dump, std::vector<plaq_ptr>& pl) {
    bool tmp;
    for (plit_t pit=pl.begin(); pit!=pl.end(); ++pit)  {
        dump>>tmp;
        (*pit)->set_value(tmp);
    }
    return dump;
}

alps::ODump& operator<<(alps::ODump& dump, const std::vector<vert_ptr>& vt) {
    for (const_vit_t vit=vt.begin(); vit!=vt.end(); ++vit) 
        dump<<(*vit)->get_value();
    return dump;
}

alps::IDump& operator>>(alps::IDump& dump, std::vector<vert_ptr>& vt) {
    bool tmp;
    for (vit_t vit=vt.begin(); vit!=vt.end(); ++vit)  {
        dump>>tmp;
        (*vit)->set_value(tmp);
    }
    return dump;
}



/*******  base class site  ***********/
site::site() : value(false) { }

void site::flip() {
    value=!value;
}


/********  class spin  **************/
spin::spin(int geo): site(), geometry(geo), neighbor_in_next_replica(this)  { }

void spin::add_neighbor(plaq_ptr nb) {
    plaq_ptr new_nb(nb);
    p_neighbors.push_back(new_nb);
}

void spin::add_neighbor(vert_ptr nb) {
    vert_ptr new_nb(nb);
    v_neighbors.push_back(new_nb);
}

void spin::flip_and_flip_plaqs() {
    value=!value;
    for (plit=p_neighbors.begin(); plit!=p_neighbors.end(); ++plit) 
        (*plit)->flip();
}

void spin::flip_and_flip_verts() {
    value=!value;
    for (vit=v_neighbors.begin(); vit!=v_neighbors.end(); ++vit) 
        (*vit)->flip();
}

int spin::get_weight_from_plaqs() {
    energy = 0;
    for (plit=p_neighbors.begin(); plit!=p_neighbors.end(); ++plit) {
        energy+=(*plit)->get_value();
    }
    return energy;
}

int spin::get_weight_from_verts() {
    energy = 0;
    for (vit=v_neighbors.begin(); vit!=v_neighbors.end(); ++vit) {
        energy+=(*vit)->get_value();
    }
    return energy;
}

/********** class interaction *********/
interaction::interaction(): site() { }

void interaction::add_neighbor(spin_ptr nb) {
    spin_ptr new_nb(nb);
    neighbors.push_back(new_nb);
}

void interaction::flip_neighbors() {
    for (spit = neighbors.begin(); spit!=neighbors.end(); ++spit)
        (*spit)->flip();
}


