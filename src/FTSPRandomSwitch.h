#ifndef MYMC_HPP
#define MYMC_HPP
#include "toric3DFTSingleSpin.h"

class FTSPRandomSwitch : public toricFTSP{

public :
    FTSPRandomSwitch(const alps::ProcessList& where,const alps::Parameters& p,int node);
private :

    void do_measurements();
    bool isboth();

};
typedef alps::scheduler::SimpleMCFactory<FTSPRandomSwitch> ToricFactory;
#endif
