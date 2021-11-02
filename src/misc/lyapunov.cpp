#include "lyapunov.hpp"

lyapunov::lyapunov(int my_id,int num_procs, int root_proc)
    :my_id(my_id),
     num_procs(num_procs),
     root_proc(root_proc)
{}
void lyapunov::compute(){
    
    double dt = 0.001; //time step
    double time = 10; //total run time
    int num_steps = int(time/dt); //total number of steps
    int renorm = 100; //steps between renormalization
    int num_cycles = num_steps/renorm; //number of renormalizations
    
    vector<double> Q1(nuc_beads,0), Q2(nuc_beads,0);
    vector<double> P1(nuc_beads,0), P2(nuc_beads,0);
    matrix<double> x1(elec_beads,num_steps,0), x2(elec_beads,num_steps,0);
    matrix<double> p1(elec_beads,num_steps,0), p2(elec_beads,num_steps,0);
    
    vector<double> v1(2*nuc_beads + 2*elec_beads*num_states);
    vector<double> v2(2*nuc_beads + 2*elec_beads*num_states);
    vector<double> v_new(2*nuc_beads + 2*elec_beads*num_states);
    vector<double> dv(2*nuc_beads + 2*elec_beads*num_states);
    vector<double> d(num_cycles);
    double d0, dv_mag;

    
    /* Initialize Forces and Integrator*/
    C_Matrix C(elec_beads, num_states,alpha);
    M_Matrix M(num_states, elec_beads, beta_elec_beads);
    dM_Matrix_dQ dMdQ(elec_beads, num_states, beta_elec_beads, M);
    theta_mixed theta(num_states,nuc_beads,elec_beads,C,M);
    theta_mixed_dQ theta_dQ(num_states,nuc_beads,elec_beads,C,M,dMdQ);
    theta_mixed_dElec theta_dElec(num_states,elec_beads,alpha,C,M);

    mvrpmd_mixed_forces F(nuc_beads, elec_beads, num_states, mass,
                          beta_nuc_beads, alpha, theta, theta_dQ, theta_dElec);

    ABM_MVRPMD myABM1(F,dt,num_states,nuc_beads,elec_beads);
    ABM_MVRPMD myABM2(F,dt,num_states,nuc_beads,elec_beads);

    v1 = to_phase(Q1,P1,x1,p1);
    v2 = to_phase(Q2,P2,x2,p2);
    dv = v_diff(v2,v1);
    d0 = mag(dv);
    
    myABM1.initialize_rk4(Q1, P1, x1, p1);
    myABM2.initialize_rk4(Q2, P2, x2, p2);

    for (int cycle=0; cycle < num_cycles; cycle++) {
        for (int step=0; step < renorm; step++) {
            myABM1.take_step(Q1,P1,x1,p1);
            myABM2.take_step(Q2,P2,x2,p2);
        }
        v1 = to_phase(Q1,P1,x1,p1);
        v2 = to_phase(Q2,P2,x2,p2);
        dv = v_diff(v2,v1);
        dv_mag = mag(dv);
        dv = d0 * (dv/dv_mag);
        v_new = v1 + dv;
        from_phase(Q2,P2,x2,p2,v_new);
        d(cycle) = dv_mag;
        
        myABM2.initialize_rk4(Q2, P2, x2, p2);
    }
    
    vector<double> k(num_cycles,0);
    k = LCE(d,d0,renorm,num_cycles);
    
}

vector<double> lyapunov::v_diff(const vector<double> &v1,
                                const vector<double> &v2){
    return v1 - v2;
}

double lyapunov::mag(const vector<double> &v){

    return inner_prod(v,v);
}

vector<double> lyapunov::to_phase(const vector<double> &Q,const vector<double> &P,
                                 const matrix<double> &x,const matrix<double> &p){

    vector<double> v(2*nuc_beads + 2*elec_beads*num_states);

    for (int bead=0; bead<nuc_beads; bead++) {
        v(bead) = Q(bead);
        v(nuc_beads + bead) = P(bead);
    }

    for (int bead=0; bead<elec_beads; bead++) {
        for (int state=0; state<num_states; state++) {
            v(2*nuc_beads + bead*num_states + state) = x(bead,state);
            v(2*nuc_beads + elec_beads*num_states + bead*num_states + state) = p(bead,state);
        }
    }

    return v;
}

void lyapunov::from_phase(vector<double> &Q,vector<double> &P,
                                    matrix<double> &x,matrix<double> &p,
                                    const vector<double> &v){


    for (int bead=0; bead<nuc_beads; bead++) {
        Q(bead) = v(bead);
        P(bead) = v(nuc_beads + bead);
    }

    for (int bead=0; bead<elec_beads; bead++) {
        for (int state=0; state<num_states; state++) {
            x(bead,state) = v(2*nuc_beads + bead*num_states + state);
            p(bead,state) = v(2*nuc_beads + elec_beads*num_states + bead*num_states + state);
        }
    }
}

void lyapunov::set_system(int nuc_beadsIN, int elec_beadsIN, int num_statesIN,
                                 double massIN,double betaIN, double beta_nuc_beadsIN,
                                 double beta_elec_beadsIN, double alphaIN){
    nuc_beads = nuc_beadsIN;
    elec_beads = elec_beadsIN;
    num_states = num_statesIN;
    mass = massIN;
    beta = betaIN;
    beta_nuc_beads = beta_nuc_beadsIN;
    beta_elec_beads = beta_nuc_beadsIN;
    alpha = alphaIN;
}

vector<double> lyapunov::LCE(const vector<double> &d, double d0, double renorm, int num_cycles){
    
    vector<double> k(num_cycles,0);
    k(0) = log(fabs(d(0)/fabs(d0)))/renorm;
    
    for (int i=1; i<num_cycles; i++) {
        k(i) =  ((renorm * i * k(i-1)) + log(fabs(d(i)/fabs(d0))))/(renorm * (i+1));
    }
    
    return k;
}

void lyapunov::write_LCE(const vector <double> &k, int num_cycles, double stride_dt){
    
    std::ofstream myFile;
    myFile.open("lce");
    
    for (int i=0; i<num_cycles; i++) {
        myFile << i*stride_dt << " " << k(i) << std::endl;
     }
    
    myFile.close();
}
