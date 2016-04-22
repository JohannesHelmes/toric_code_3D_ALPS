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
    void set_value(bool newvalue){ value = newvalue; }
    friend alps::IDump& operator>>(alps::IDump& dump, std::vector<spin_ptr>& sp);
    friend alps::IDump& operator>>(alps::IDump& dump, std::vector<plaq_ptr>& sp);
    friend alps::IDump& operator>>(alps::IDump& dump, std::vector<vert_ptr>& sp);
    friend alps::IDump& operator>>(alps::IDump& dump, std::vector<inter_ptr>& sp);

public:
    int get_value() const { return value? 1 : -1 ; }
    void flip();
    site();
};


/***********************  SITE  ********************************/
class spin : public site {
    
private:
    int energy;
    const int geometry;
    std::vector<plaq_ptr> p_neighbors;
    std::vector<vert_ptr> v_neighbors;
    spin_ptr neighbor_in_next_replica;
    plit_t plit;
    vit_t vit;

public:
    spin(int geo);
    void add_neighbor(plaq_ptr nb);
    void add_neighbor(vert_ptr nb);
    void set_ninr(spin_ptr ns) {neighbor_in_next_replica = ns; }
    void flip_and_flip_plaqs();
    void flip_and_flip_verts();
    int get_weight_from_plaqs(); 
    int get_weight_from_verts(); 
    int get_geometry() const { return geometry; }
    spin_ptr get_next() {return neighbor_in_next_replica; }

};

/*************************  INTERACTIONS (many body terms)   *******************/
class interaction : public site { // plaquette or vertex
public:
    interaction();
    void add_neighbor(spin_ptr nb);
    void flip_neighbors();
    void add_label(int nlabel) {label = nlabel; }
    int get_label() { return label; }
    void set_boundary(bool nbound) {boundary = nbound; }
    bool get_boundary() { return boundary; }
    void set_ninr(inter_ptr nv) {neighbor_in_next_replica = nv; }
    inter_ptr get_next() {return neighbor_in_next_replica; }

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

    friend void spin::flip_and_flip_plaqs();
    friend void spin::flip_and_flip_verts();
};

/*******************  PLAQUETTE and VERTEX childs  **************************/
class plaquette : public interaction { using interaction::interaction; }; //inherit constructor
class vertexx : public interaction { using interaction::interaction; };


/****************   NON ISOTROPIC VARIANTS  ************************/
class spin_z : public spin {
private:
    const double J;
public:
    spin_z(int geo, double nJ);
    int get_value() const  { return value? J : -J; }   //overwrites get_value() of site class
};
/*
class vertex_xxz;
typedef std::shared_ptr<vertex_xxz> vertxxz_ptr;
typedef std::vector<vertxxz_ptr>::iterator vxxzit_t;

class vertex_xxz : public site { 
public:
    vertex_xxz();
    void add_neighbor(spin_ptr nb, bool z_neighbor);
    void flip_neighbors();
    void add_label(int nlabel) {label = nlabel; }
    int get_label() { return label; }
    void set_boundary(bool nbound) {boundary = nbound; }
    bool get_boundary() { return boundary; }

    spit_t get_neighbors_xy_begin() {return neighbors_xy.begin(); }
    spit_t get_neighbors_xy_end() {return neighbors_xy.end(); }
protected:
    std::vector<spin_ptr> neighbors_xy, neighbors_z;
    spit_t spit;
    int label;
    bool boundary;

    friend void spin::flip_and_flip_plaqs();
    friend void spin::flip_and_flip_verts();
};
*/
#endif  /* SITE_H */
