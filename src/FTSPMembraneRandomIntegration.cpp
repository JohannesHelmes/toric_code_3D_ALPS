#include "FTSPMembraneRandomIntegration.h"

using namespace boost;

FTSPMembraneRandomIntegration::FTSPMembraneRandomIntegration(const alps::ProcessList& where,const alps::Parameters& p,int node) : toricFTSPMembrane(where,p,node) {


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


void FTSPMembraneRandomIntegration::do_measurements() {
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
    measurements["Energy2"] << double(NofD)/(numspins)*double(NofD)/(numspins);

}

