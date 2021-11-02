#ifndef mvrpmd_Esplit_ham_hpp
#define mvrpmd_Esplit_ham_hpp

#include "StateIndepPot.h"
#include "theta_Esplit.hpp"
#include "GTerm.h"
#include "mvrpmd_mixed.hpp"

class mvrpmd_Esplit_ham : public mvrpmd_mixed{

public:
    
    mvrpmd_Esplit_ham(double beta_num_beads,StateIndepPot &V0_In,
                      GTerm &GIn,theta_Esplit &thetaIn);

    double get_energy(const vector<double> &Q, const matrix<double> &x,
                      const matrix<double> &p);

    double get_energy();
    
    double get_energy_dyn(double mass,const vector<double> &Q, const vector<double> &P,
                          const matrix<double> &x,const matrix<double> &p);
    
    double get_sign();

private:

    /* Private data*/
    StateIndepPot * V0;
    GTerm * G;
    theta_Esplit * theta;
    double one_beta_num_beads; // 1/beta_num_beads
};

#endif
