#include "site.h"

/*******  base class site  ***********/
site::site() : value(false) { }



/********  class spin  **************/
spin::spin(int inReps) : site(), weight(-4*inReps) { }

void spin::add_neighbor(plaq_ptr nb) {
    plaq_ptr new_nb(nb);
    neighbors.push_back(new_nb);
}

void spin::flip() {
    value=!value;
    for (plit=neighbors.begin(); plit!=neighbors.end(); ++plit) 
        (*plit)->flip();
}

void spin::mod_weight(bool plusminus) {
    weight += (plusminus)? 2: -2 ;
}

/********** class plaquette *********/
plaquette::plaquette(): site() { }

void plaquette::add_neighbor(spin_ptr nb) {
    spin_ptr new_nb(nb);
    neighbors.push_back(new_nb);
}

void plaquette::flip() {
    value=!value;
    for (spit=neighbors.begin(); spit!=neighbors.end(); ++spit)
        (*spit)->mod_weight(value);
}
