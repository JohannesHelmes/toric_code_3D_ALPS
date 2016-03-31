#include "measurer.h"

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

    int ct =0;
    for (const_spit_t spit=spins.begin(); spit!=spins.end(); ++spit) {
        if ( (*spit)->get_geometry() == 2)
            switcher.push_back(*spit);
        if (ct > N_spins_per_rep) 
            break;
    }
}


bool switching::isboth() {
    for (const_spit_t spit=switcher.begin(); spit!=switcher.end(); ++spit) {
        reference = (*spit)->get_value();
        runner = (*spit);
        runner->get_next();
        while (runner != *spit) {
            if ( runner->get_value () != reference)
                return false;
            runner->get_next();
        }
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
