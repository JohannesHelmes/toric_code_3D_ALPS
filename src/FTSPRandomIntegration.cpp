#include "FTSPRandomIntegration.h"

using namespace boost;

FTSPRandomIntegration::FTSPRandomIntegration(const alps::ProcessList& where,const alps::Parameters& p,int node) : toricFTSP(where,p,node) {

    sit=sites().first;
    for (int i=0; i<IncStep.length(); ++i,++sit) {
        if (sit==sites().second)
            break;
        while (site_type(*sit)!=0)
            ++sit;
        if (sit==sites().second)
            break;
        if (IncStep[i]!='0')
            geom[*sit]=true;
    }

    /* print lattice connections
    for (sit=sites().first; sit!=sites().second; ++sit) {
        if (site_type(*sit)==2) {
            cout<<*sit<<": ";
            for (nit=neighbors(*sit).first; nit!=neighbors(*sit).second; ++nit) {
                cout<<*nit<<" ("<<geom[*nit]<<") ,";
            }
            cout<<endl;
        }
    }
    */

    
    measurements << alps::RealObservable("Energy");
    measurements << alps::RealObservable("Energy2");
    cout<<"init ok"<<endl;
}


void FTSPRandomIntegration::do_measurements() {
    /*
    int Nof_def=0;
    for (int i=0; i<n; ++i) {
        for (int j=0; j<numsites; ++j) {
            if (site_type(j)==2)
                Nof_def+=vertex_defects[i][j]? 1: -1;
        }
    }
    */
    measurements["Energy"] << double(NofD)/(numspins);
    measurements["Energy2"] << double(NofD*NofD)/(numspins*numspins);
}

