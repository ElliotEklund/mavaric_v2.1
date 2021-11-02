#include "rpmd_ham.hpp"

rpmd_ham::rpmd_ham(int nuc_beads, double beta_num_beads,
                       SpringEnergy &V_springIn,StateIndepPot &V0_In)
    :nuc_beads(nuc_beads),
     one_beta_num_beads(1.0/beta_num_beads)
{
    V_spring = &V_springIn;
    V0 = &V0_In;
}
double rpmd_ham::get_energy(const vector<double> &Q){

    double V_spring_temp = V_spring->get_springEnergy(Q);
    double V0_temp = V0->get_V0(Q);
    return V_spring_temp + V0_temp;
}
double rpmd_ham::get_energy(){

    double V_spring_temp = V_spring->get_springEnergy();
    double V0_temp = V0->get_V0();
    return V_spring_temp + V0_temp;
}
double rpmd_ham::get_energy_dyn(double mass,const vector<double> &Q, const vector<double> &P){

    double kin_energy = 0.5 * inner_prod(P,P)/mass;
    double V_spring_temp = V_spring->get_springEnergy(Q);
    double V0_temp = V0->get_V0(Q);
    return kin_energy + V_spring_temp + V0_temp;
}
