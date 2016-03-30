#include<vector>
#include<memory>
#include<alps/osiris/dump.h>

class spin;
class plaquette;
class vertexx; //vertex is somehow used by ALPS

//typedef std::shared_ptr<site> site_ptr;
//typedef std::vector<site_ptr>::iterator site_t;

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
alps::IDump& operator>>(alps::IDump& dump, std::vector<spin_ptr>& sp);
alps::IDump& operator>>(alps::IDump& dump, std::vector<plaq_ptr>& sp);
alps::IDump& operator>>(alps::IDump& dump, std::vector<vert_ptr>& sp);



class site {
protected:
    bool value;
    void set_value(bool newvalue){ value = newvalue; }
    friend alps::IDump& operator>>(alps::IDump& dump, std::vector<spin_ptr>& sp);
    friend alps::IDump& operator>>(alps::IDump& dump, std::vector<plaq_ptr>& sp);
    friend alps::IDump& operator>>(alps::IDump& dump, std::vector<vert_ptr>& sp);

public:
    int get_value() const { return value? 1 : -1 ; }
    void flip();
    site();
    //IMPLEMENT Operator>> and Operator<< HERE!!!
};


class spin : public site {
    
private:
    int energy;
    const int geometry;
    std::vector<plaq_ptr> p_neighbors;
    std::vector<vert_ptr> v_neighbors;
    plit_t plit;
    vit_t vit;

public:
    spin(int geo);
    void add_neighbor(plaq_ptr nb);
    void add_neighbor(vert_ptr nb);
    void flip_and_flip_plaqs();
    void flip_and_flip_verts();
    int get_weight_from_plaqs(); 
    int get_weight_from_verts(); 
    int get_geometry() { return geometry; }

};

class interaction : public site { // plaquette or vertex
public:
    interaction();
    void add_neighbor(spin_ptr nb);
    void flip_neighbors();
    void add_label(int nlabel) {label = nlabel; }
    int get_label() { return label; }
    void set_boundary(bool nbound) {boundary = nbound; }
    bool get_boundary() { return boundary; }
    const_spit_t get_neighbors_begin() {return neighbors.begin(); }
    const_spit_t get_neighbors_end() {return neighbors.end(); }
protected:
    std::vector<spin_ptr> neighbors;
    spit_t spit;
    int label;
    bool boundary;

    friend void spin::flip_and_flip_plaqs();
    friend void spin::flip_and_flip_verts();
};

class plaquette : public interaction { using interaction::interaction; }; //inherit constructor
class vertexx : public interaction { using interaction::interaction; };


