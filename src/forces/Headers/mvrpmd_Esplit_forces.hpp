#ifndef mvrpmd_Esplit_forces_hpp
#define mvrpmd_Esplit_forces_hpp

#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>

#include "dStateIndep_dQ.hpp"
#include "theta_Esplit.hpp"
#include "theta_Esplit_dQ.hpp"
#include "theta_Esplit_dElec.hpp"
#include "mv_forces_temp.hpp"

using namespace boost::numeric::ublas;

class mvrpmd_Esplit_forces : public mv_forces_temp {
    
public:
    mvrpmd_Esplit_forces(int elec_beads, int num_states,
                        double mass, double betaIN, double alpha,
                        theta_Esplit &theta_IN, theta_Esplit_dQ &theta_dQ_IN,
                        theta_Esplit_dElec &theta_dElec_IN);
    
    /*
     Update dHdP to reflect the state of P
     P: vector of bead momentum
     */
    void update_dHdP(const vector<double> &P);
    
    /*
     Update dHdQ to reflect the state of Q,x,p
     Q: vector of bead positions
     x: matrix of x mapping variables
     p: matrix of p mapping variables
     */
    void update_dHdQ(const vector<double> &Q, const matrix<double> &x,
                        const matrix<double> &p);
    /*
     Update dHdx to reflect the state of Q,x,p
     Q: vector of bead positions
     x: matrix of x mapping variables
     p: matrix of p mapping variables
     */
    void update_dHdx(const vector<double> &Q, const matrix<double> &x,
                        const matrix<double> &p);
    /*
     Update dHdp to reflect the state of Q,x,p
     Q: vector of bead positions
     x: matrix of x mapping variables
     p: matrix of p mapping variables
     */
    void update_dHdp(const vector<double> &Q, const matrix<double> &x,
                        const matrix<double> &p);

    /* Update state of dHdP, dHdQ, dHdx, and dHdp.*/
    void update_Forces(const vector<double> &Q,const vector<double> &P,
                       const matrix<double> &x,const matrix<double> &p);
    
/* Mutators */
    /* Return dHdQ*/
    const vector<double> & get_dHdQ();
  
    /* Return dHdP*/
    const vector<double> & get_dHdP();
   
    /* Return dHdx*/
    const matrix<double> & get_dHdx();

    /* Return dHdp*/
    const matrix<double> & get_dHdp();
    
    double get_sign(const vector<double> &Q, const matrix<double> &x,
                        const matrix<double> &p);

/* Debugging*/
    void print_dHdQ();
    void print_dHdP();
    void print_dHdx();
    void print_dHdp();
    
private:
    
/* Data */
    double ONE_mass; //1.0/mass
    double beta; //beta/nuc_beads
    double ONE_beta; //1.0/beta_nuc_beads
    double TWO_beta; //2.0/beta_nuc_beads
    int elec_beads; //number of electronic beads
    int num_states; //number of electronic states
    double coeff_ONE_theta; //ONE_beta_nuc_beads /theta
    double alpha; //mapping variable prefactor
    double x_alpha, p_alpha; //x and p mapping variable prefactors

    vector<double> dHdQ, dHdP; //derivative of Hamiltonian wrt Q,P
    matrix<double> dHdx, dHdp; //derivative of Hamiltonian wrt x,p

    vector<double> dV0_dQ_vec;
    vector<double> dThetaMTS_dQ_vec;

    matrix<double> dThetaMTS_dx_vec;
    matrix<double> dThetaMTS_dp_vec;
    
/* Objects */
    dStateIndep_dQ dV0_dQ;
    
    theta_Esplit * theta;
    theta_Esplit_dQ * theta_dQ;
    theta_Esplit_dElec * theta_dElec;
};

#endif
