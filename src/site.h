#ifndef SITE_H
#define SITE_H

#include<vector>
#include<memory>
#include<alps/osiris/dump.h>

class interaction;
class spin;
class plaquette;
class vertexx; //vertex is somehow used by ALPS

//typedef std::shared_ptr<site> site_ptr;
//typedef std::vector<site_ptr>::iterator site_t;

typedef std::shared_ptr<interaction> inter_ptr;
typedef std::vector<inter_ptr>::iterator iit_t;
typedef std::vector<inter_ptr>::const_iterator const_iit_t;
typedef std::shared_ptr<spin> spin_ptr;
typedef std::vector<spin_ptr>::iterator spit_t;
typedef std::vector<spin_ptr>::const_iterator const_spit_t;
typedef std::shared_ptr<plaquette> plaq_ptr;
typedef std::vector<plaq_ptr>::iterator plit_t;
typedef std::vector<plaq_ptr>::const_iterator const_plit_t;
typedef std::shared_ptr<vertexx> vert_ptr;
typedef std::vector<vert_ptr>::iterator vit_t;
typedef std::vector<vert_ptr>::const_iterator const_vit_t;

/** must be declared as non-members **/
alps::ODump& operator<<(alps::ODump& dump, const std::vector<spin_ptr>& sp);
alps::ODump& operator<<(alps::ODump& dump, const std::vector<plaq_ptr>& sp);
alps::ODump& operator<<(alps::ODump& dump, const std::vector<vert_ptr>& sp);
alps::ODump& operator<<(alps::ODump& dump, const std::vector<inter_ptr>& sp);
alps::IDump& operator>>(alps::IDump& dump, std::vector<spin_ptr>& sp);
alps::IDump& operator>>(alps::IDump& dump, std::vector<plaq_ptr>& sp);
alps::IDump& operator>>(alps::IDump& dump, std::vector<vert_ptr>& sp);
alps::IDump& operator>>(alps::IDump& dump, std::vector<inter_ptr>& sp);



/**********************  BASE CLASS ****************************/
class site {
protected:
    bool value;
    const short orientation;
    const int name;
    void set_value(bool newvalue){ value = newvalue; }
    friend alps::IDump& operator>>(alps::IDump& dump, std::vector<spin_ptr>& sp);
    friend alps::IDump& operator>>(alps::IDump& dump, std::vector<plaq_ptr>& sp);
    friend alps::IDump& operator>>(alps::IDump& dump, std::vector<vert_ptr>& sp);
    friend alps::IDump& operator>>(alps::IDump& dump, std::vector<inter_ptr>& sp);

public:
    double get_value() const { return value? 1 : -1 ; }
    int get_orientation() const { return orientation; }
    int get_name() const { return name; }
    void flip();
    site(int name, short orientation);
};


/***********************  SPIN  ********************************/
class spin : public site {
    
private:
    int energy;
    const short geometry;
    std::vector<plaq_ptr> p_neighbors;  //ToDo : remove these two, work only with interaction_neighbors
    std::vector<vert_ptr> v_neighbors;
    std::vector<inter_ptr> interaction_neighbors;
    std::vector<inter_ptr> dual_interaction_neighbors;
    spin_ptr neighbor_in_next_replica;
    plit_t plit;
    vit_t vit;

public:
    spin(int name, short geo, short orientation); // orientation : x = 0, y = 1, z = 2
    double get_value() const { return value? 1 : -1 ; }

    void add_neighbor(plaq_ptr nb);
    void add_neighbor(vert_ptr nb);
    void copy_neighbors_internally(short excitation); //not very elegant: copies p_neighbors to interaction_neigbors and v_neighbors to dual_... or vv
    void set_ninr(spin_ptr ns) {neighbor_in_next_replica = ns; }
    void flip_and_flip_plaqs();
    void flip_and_flip_verts();
    int get_weight_from_plaqs(); 
    int get_weight_from_verts(); 
    int get_geometry() const { return geometry; }
    spin_ptr get_next() {return neighbor_in_next_replica; }
    spin_ptr move(short direction); //0 = x direction, 1 = y direction, 2 = z direction

    const_iit_t get_interaction_neighbors_begin() const {return interaction_neighbors.begin(); }
    const_iit_t get_interaction_neighbors_end() const {return interaction_neighbors.end(); }

    const_iit_t get_dual_interaction_neighbors_begin() const {return dual_interaction_neighbors.begin(); }
    const_iit_t get_dual_interaction_neighbors_end() const {return dual_interaction_neighbors.end(); }
};

/*************************  INTERACTIONS (many body terms)   *******************/
class interaction : public site { // plaquette or vertex
public:
    interaction(int name, short orientation=6);
    void add_neighbor(spin_ptr nb);
    void flip_neighbors();
    void add_label(int nlabel) {label = nlabel; }
    int get_label() { return label; }
    void set_boundary(bool nbound) {boundary = nbound; }
    bool get_boundary() { return boundary; }
    void set_ninr(inter_ptr nv) {neighbor_in_next_replica = nv; }
    void set_nindir(inter_ptr nv, short dir) {lattice_neighbor_in_dir[dir] = nv; } //sets neighbor in x,y, or z direction
    inter_ptr get_next() const {return neighbor_in_next_replica; }
    inter_ptr move(short direction) const {return lattice_neighbor_in_dir[direction]; } //0 = x direction, 1 = y direction, 2 = z direction

    const_spit_t get_neighbors_begin() const {return neighbors.begin(); }
    const_spit_t get_neighbors_end() const {return neighbors.end(); }
    spit_t get_neighbors_begin() {return neighbors.begin(); }
    spit_t get_neighbors_end() {return neighbors.end(); }
protected:
    std::vector<spin_ptr> neighbors;
    spit_t spit;
    int label;
    bool boundary;
    inter_ptr neighbor_in_next_replica;
    std::vector<inter_ptr> lattice_neighbor_in_dir;

    friend void spin::flip_and_flip_plaqs();
    friend void spin::flip_and_flip_verts();
};

/*******************  PLAQUETTE and VERTEX childs  **************************/
class plaquette : public interaction { using interaction::interaction; }; //inherit constructor
class vertexx : public interaction { using interaction::interaction; };

#endif  /* SITE_H */
