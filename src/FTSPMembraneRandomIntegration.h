#ifndef MYMC_HPP
#define MYMC_HPP
#include "toric3DFTMembraneSingleSpin.h"

class FTSPMembraneRandomIntegration : public toricFTSPMembrane{

public :
    FTSPMembraneRandomIntegration(const alps::ProcessList& where,const alps::Parameters& p,int node);

private :
    void do_measurements();

};
typedef alps::scheduler::SimpleMCFactory<FTSPMembraneRandomIntegration> ToricFactory;
#endif
