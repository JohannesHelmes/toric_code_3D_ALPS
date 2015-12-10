/*********** Toric Code at finite temperature. Update: Single Spin flips only *****************/
/*********************         9.7.2015                         *******************************/

#include "FTSPRandomSwitch.h"

using namespace boost;


FTSPRandomSwitch::FTSPRandomSwitch(const alps::ProcessList& where,const alps::Parameters& p,int node) : toricFTSP(where,p,node) {

    sit=sites().first;
    for (int i=0; i<IncStep.length(); ++i,++sit) {
        if (sit==sites().second)
            break;
        while (site_type(*sit)!=0)
            ++sit;
        if (sit==sites().second)
            break;
        if (IncStep[i]=='1')
            geom[*sit]=true;
        else if (IncStep[i]=='2') {
            edge[*sit]=true;
            edge_sites.insert(*sit);
        }
    }

    cout<<"init ok"<<endl;

    measurements << alps::RealObservable("ED"); //Ensemble divided at edge
    measurements << alps::RealObservable("EG"); //Ensemble glued at edge
}


bool FTSPRandomSwitch::isboth() {
    for (set<site_descriptor>::iterator it=edge_sites.begin();it!=edge_sites.end();++it) {
        for (int i=1; i<n; ++i) {
            if (spins_[*it]!=spins_[numsites*i+*it]) {
                return false;
            }
        }
    }
    return true;
}


void FTSPRandomSwitch::do_measurements() {
    if (isboth()) {
        measurements["EG"] << 1.0;
        measurements["EG"] << 0.0;
    }
    else {
        measurements["EG"] << 0.0;
    }
}

