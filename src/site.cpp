#include "site.h"

/*******  base class site  ***********/
site::site() : value(false) { }



/********  class spin  **************/
spin::spin(): site()  { }

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

void interaction::flip() {
    value=!value;
}

