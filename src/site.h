#include<vector>
#include<memory>

class site;
class spin;
class plaquette;

typedef std::shared_ptr<site> site_ptr;
typedef std::vector<site_ptr>::iterator site_t;

typedef std::shared_ptr<spin> spin_ptr;
typedef std::vector<spin_ptr>::iterator spit_t;
typedef std::shared_ptr<plaquette> plaq_ptr;
typedef std::vector<plaq_ptr>::iterator plit_t;

class site {
protected:
    bool value;
public:
    site();
    virtual void flip()=0;
};


class spin : site {
    
private:
    int weight, energy;
    std::vector<plaq_ptr> neighbors;
    plit_t plit;
public:
    spin(int inReps);
    void add_neighbor(plaq_ptr nb);
    void flip();
    void mod_weight(bool plusminus); //False = decrease, True = increase
    int const get_weight(); 
    int const num_neighbors(); 

};

class plaquette : site {
public:
    plaquette();
    void add_neighbor(spin_ptr nb);
    int get_value(){ return value? 1 : -1 ; }
private:
    std::vector<spin_ptr> neighbors;
    spit_t spit;
    void flip();

    friend void spin::flip();
};

class vertex : site {
    vertex();
private:
    void flip() { };
};

