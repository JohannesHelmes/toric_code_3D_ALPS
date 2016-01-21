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
switching::switching(alps::ObservableSet& msmt, bool& cnncted) : measurer(msmt), connected(cnncted) {
}

void switching::measure() {
}
