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
switching::switching(alps::ObservableSet& msmt, std::vector<spin_ptr>& p_spins) : measurer(msmt), spins(p_spins) {
    measurements << alps::RealObservable("EG"); //Ensemble glued at edge
}

bool switching::isboth() {
    /*
    for (set<site_descriptor>::iterator it=edge_sites.begin();it!=edge_sites.end();++it) {
        for (int i=1; i<n; ++i) {
            if (spins_[*it]!=spins_[numsites*i+*it]) {
                return false;
            }
        }
    }
    */
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
