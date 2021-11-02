#include "theta_Esplit_dQ.hpp"

theta_Esplit_dQ::theta_Esplit_dQ(int num_states, int elec_beads,
                                 C_Matrix &C_IN,M_Matrix &M_IN,
                                 dM_Matrix_dQ &M_dQ_IN)

    :num_states(num_states), elec_beads(elec_beads),

    f_chain(elec_beads,zero_matrix<std::complex<double> >(num_states,num_states)),
    b_chain(elec_beads,zero_matrix<std::complex<double> >(num_states,num_states)),

    theta_dQ_vec(1,0.0),

    temp1(num_states,num_states,0.0),temp2(num_states,num_states,0.0),
    fchain_temp(num_states,num_states,0.0),bchain_temp(num_states,num_states,0.0)
{
    C = &C_IN;
    M = &M_IN;
    M_dQ = & M_dQ_IN;
}
void theta_Esplit_dQ::update_theta_dQ(const vector<double> &Q){
    
    std::complex<double> tr (0,0);
     
    update_f_chain();
    update_b_chain();
    
    M_dQ->update_dM_dQ_vec(Q);
    
    vector<double> theta_dQ_temp(1,0);
    
    for (int bead=0; bead<elec_beads; bead++) {
        noalias(temp1) = prod(f_chain(bead), M_dQ->get_dMdQ_alpha(0));
        noalias(temp2) = prod(temp1, b_chain(bead));
        tr = trace<std::complex<double> >(temp2,num_states);
        theta_dQ_temp(0) += tr.real();
        tr = std::complex<double> (0,0);
    }
    
    theta_dQ_vec = theta_dQ_temp;
}

void theta_Esplit_dQ::update_f_chain(){

    f_chain[0] = C->get_C_alpha(0);

    for (int bead=1; bead<elec_beads; bead++) {
        noalias(fchain_temp) = prod(M->get_M_alpha(0), C->get_C_alpha(bead));
        f_chain[bead] = prod(f_chain[bead-1], fchain_temp);
    }
}
void theta_Esplit_dQ::update_b_chain(){

    b_chain[elec_beads-1] = identity_matrix<std::complex<double> > (num_states);

    for (int bead = elec_beads - 2; bead > -1; bead--) {
        noalias(bchain_temp) = prod(C->get_C_alpha(bead+1), M->get_M_alpha(0));
        b_chain[bead] = prod(bchain_temp, b_chain[bead+1]);
    }
}
vector<double> & theta_Esplit_dQ::get_theta_dQ_vec(){
    return theta_dQ_vec;
}

