#ifndef mvrpmd_Esplit_esti_hpp
#define mvrpmd_Esplit_esti_hpp

#include "StateIndepPot.h"
#include "theta_Esplit.hpp"
#include "theta_Esplit_dBeta.hpp"

class mvrpmd_Esplit_esti{
    
public:
    
    mvrpmd_Esplit_esti(int num_beads, double beta_num_beads,
                       StateIndepPot &V0In,
                       theta_Esplit &thetaIN, theta_Esplit_dBeta &dthetaIN);
    
    double get_estimator();
    
    double get_estimator(const vector<double> &Q,const matrix<double> &x,
                         const matrix<double> &p);
    
private:
    
/* Data */
    double ONE_HALF_beta_num_beads; // 0.5 * (1/beta_num_beads)
    double ONE_num_beads;
    
/* Objects*/
    StateIndepPot * V0;
    theta_Esplit * theta;
    theta_Esplit_dBeta * dtheta;
};

#endif
