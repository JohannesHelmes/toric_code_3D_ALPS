#include <alps/alea/observableset.h>

class measurer {
    
public:
    measurer(alps::ObservableSet& msmt);
    virtual void measure()=0;
protected:
    alps::ObservableSet& measurements;

};

class thermo_int : public measurer {
 
public:
    thermo_int(alps::ObservableSet& msmt,int& nofd,const int nspns);
    void measure();
private:
    int& NofD;
    int const numspins;
};

class switching : public measurer {

public:
    switching(alps::ObservableSet& msmt,bool& cnnctd);
    void measure();
private:
    bool& connected;
};
