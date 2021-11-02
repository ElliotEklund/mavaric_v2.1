#include "rpmd_estimator.hpp"

rpmd_estimator::rpmd_estimator(int num_beads, double beta_num_beads,
                                    SpringEnergy &V_SpringIn, StateIndepPot &V0In)
    :ONE_HALF_beta_num_beads(0.5/beta_num_beads),
     ONE_num_beads(1.0/num_beads)
{
    V_spring = &V_SpringIn;
    V0 = &V0In;
}

double rpmd_estimator::get_estimator(){

    double V_spring_temp = V_spring->get_springEnergy();
    double V0_temp = V0->get_V0();

    return ONE_HALF_beta_num_beads + ONE_num_beads*(V0_temp - V_spring_temp);
}
double rpmd_estimator::get_estimator(const vector<double> &Q){

    double V_spring_temp = V_spring->get_springEnergy(Q);
    double V0_temp = V0->get_V0(Q);

    return ONE_HALF_beta_num_beads + ONE_num_beads*(V0_temp - V_spring_temp);
}
