#ifndef MYMC_HPP
#define MYMC_HPP
#include "toric3DFTSingleSpin.h"

class FTSPRandomIntegration : public toricFTSP{

public :
    FTSPRandomIntegration(const alps::ProcessList& where,const alps::Parameters& p,int node);

private :
    void do_measurements();

};
typedef alps::scheduler::SimpleMCFactory<FTSPRandomIntegration> ToricFactory;
#endif
