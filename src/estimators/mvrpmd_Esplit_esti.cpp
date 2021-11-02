#include "mvrpmd_Esplit_esti.hpp"

mvrpmd_Esplit_esti::mvrpmd_Esplit_esti(int num_beads, double beta_num_beads,
                                       StateIndepPot &V0In,
                                       theta_Esplit &thetaIN,
                                       theta_Esplit_dBeta &dthetaIN)
    :ONE_HALF_beta_num_beads(0.5/beta_num_beads),
     ONE_num_beads(1.0/num_beads)
{
    V0 = &V0In;
    theta = &thetaIN;
    dtheta = &dthetaIN;
}
double mvrpmd_Esplit_esti::get_estimator(){
    
    double V0_temp = V0->get_V0();
    double theta_temp = theta->get_theta();
    double sgn_theta_temp = theta->get_signTheta();
    double dtheta_temp = dtheta->get_dtheta();
    

    
    return (ONE_HALF_beta_num_beads + ONE_num_beads*(V0_temp) -
            dtheta_temp/theta_temp)*sgn_theta_temp;
}
double mvrpmd_Esplit_esti::get_estimator(const vector<double> &Q,
                                        const matrix<double> &x,
                                        const matrix<double> &p){
    
    double V0_temp = V0->get_V0(Q);
    double theta_temp = theta->get_theta(Q,x,p);
    double sgn_theta_temp = theta->get_signTheta();
    double dtheta_temp = dtheta->get_dtheta();
   
    return (ONE_HALF_beta_num_beads + ONE_num_beads*(V0_temp)
            - dtheta_temp/theta_temp)*sgn_theta_temp;
}
