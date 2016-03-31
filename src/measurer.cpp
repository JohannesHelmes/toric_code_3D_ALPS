#include "measurer.h"

using namespace std;

/* basic class measurer */
measurer::measurer(alps::ObservableSet& msmt) : measurements(msmt) {
}

/* class thermo_int */
thermo_int::thermo_int(alps::ObservableSet& msmt,int& nofd,const int nspns) : measurer(msmt), NofD(nofd), numspins(nspns) {
    measurements << alps::RealObservable("Energy");
    measurements << alps::RealObservable("Energy2");
}

void thermo_int::measure() {
    measurements["Energy"] << double(NofD)/(numspins);
    measurements["Energy2"] << double(NofD)/(numspins)*double(NofD)/(numspins);
}


/* class switching */
/* INCOMPLETE!!!! */
switching::switching(alps::ObservableSet& msmt, std::vector<spin_ptr>& spins, int N_spins_per_rep) : measurer(msmt)  {
    measurements << alps::RealObservable("EG"); //Ensemble glued at edge
    cout<<"Initializing switching measurement"<<endl;

    int ct =0;
    for (const_spit_t spit=spins.begin(); spit!=spins.end(); ++spit,++ct) {
        if ( (*spit)->get_geometry() == 2) {
            switcher.push_back(*spit);
        }
        if (ct >= N_spins_per_rep) 
            break;
    }
}


bool switching::isboth() {
    //int iterations ; int spin_no=0;
    for (const_spit_t spit=switcher.begin(); spit!=switcher.end(); ++spit) {
        reference = (*spit)->get_value();
        //cout<<"spin no "<<spin_no<<" has value "<<reference<<endl;
        runner = (*spit);
        runner = runner->get_next();
        //iterations=0;
        while (runner != *spit) {
            if ( runner->get_value () != reference) {
                return false;
            }
            runner = runner->get_next();
            //++iterations;
        }
        //cout<<"done "<<iterations<<" Its "<<endl;
    }
    return true;
}

void switching::measure() {
    if (isboth()) { //IMPLEMENT THIS!
        measurements["EG"] << 1.0;
        measurements["EG"] << 0.0;
    }
    else {
        measurements["EG"] << 0.0;
    }
}
