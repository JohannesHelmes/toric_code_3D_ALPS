#include<vector>
#include<memory>

class spin;
class plaquette;
class vertexx; //vertex is somehow used by ALPS

//typedef std::shared_ptr<site> site_ptr;
//typedef std::vector<site_ptr>::iterator site_t;

typedef std::shared_ptr<spin> spin_ptr;
typedef std::vector<spin_ptr>::iterator spit_t;
typedef std::shared_ptr<plaquette> plaq_ptr;
typedef std::vector<plaq_ptr>::iterator plit_t;
typedef std::shared_ptr<vertexx> vert_ptr;
typedef std::vector<vert_ptr>::iterator vit_t;

class site {
protected:
    bool value;
public:
    site();
    //virtual void flip()=0;
};


class spin : site {
    
private:
    int energy;
    std::vector<plaq_ptr> p_neighbors;
    std::vector<vert_ptr> v_neighbors;
    plit_t plit;
    vit_t vit;
public:
    spin();
    void add_neighbor(plaq_ptr nb);
    void add_neighbor(vert_ptr nb);
    void flip_and_flip_plaqs();
    void flip_and_flip_verts();
    int get_weight_from_plaqs(); 
    int get_weight_from_verts(); 

};

class interaction : site { // plaquette or vertex
public:
    interaction();
    void add_neighbor(spin_ptr nb);
    int get_value(){ return value? 1 : -1 ; }
protected:
    std::vector<spin_ptr> neighbors;
    spit_t spit;
    void flip();

    friend void spin::flip_and_flip_plaqs();
    friend void spin::flip_and_flip_verts();
};

class plaquette : public interaction { using interaction::interaction; }; //inherit constructor
class vertexx : public interaction { using interaction::interaction; };


