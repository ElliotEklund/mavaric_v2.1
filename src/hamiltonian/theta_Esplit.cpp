#include "theta_Esplit.hpp"

theta_Esplit::theta_Esplit(int num_states, int elec_beads,C_Matrix &C_In, M_Matrix &M_In)
    :num_states(num_states),
     elec_beads(elec_beads),
     gamma_mat(identity_matrix<std::complex<double> > (num_states))
{
    C = &C_In;
    M = &M_In;
}
void theta_Esplit::update_gamma_mat(const vector<double> &Q,const matrix<double> &x,
                                   const matrix<double> &p){
    
    C->update_C_vec(x, p);
    M->update_M_vec(Q);
        
    /* Chain multiplicaton */
    gamma_mat = identity_matrix<std::complex<double> > (num_states);

    for(int bead=0; bead<elec_beads; bead++){
        gamma_mat = prod(gamma_mat,C->get_C_alpha(bead));
        gamma_mat = prod(gamma_mat,M->get_M_alpha(0));
    }
}
void theta_Esplit::update_theta(const vector<double> &Q,const matrix<double> &x,
                               const matrix<double> &p){
    
    update_gamma_mat(Q, x, p);
    std::complex<double> tr (0.0,0.0);
    tr = trace<std::complex<double> >(gamma_mat,num_states);
    theta = tr.real();
}
double theta_Esplit::get_theta(const vector<double> &Q,const matrix<double> &x,
                              const matrix<double> &p){
    update_theta(Q, x, p);
    return theta;
}
double theta_Esplit::get_theta(){return theta;}

double theta_Esplit::get_signTheta(){return sign(theta);}

double theta_Esplit::get_signTheta(const vector<double> &Q,const matrix<double> &x,
                                  const matrix<double> &p){
    
    update_theta(Q,x,p);
    return sign(theta);
}
matrix<std::complex<double> > theta_Esplit::get_gamm(){return gamma_mat;}
