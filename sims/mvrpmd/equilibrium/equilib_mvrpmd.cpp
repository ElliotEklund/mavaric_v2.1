#include "equilib_mvrpmd.hpp"

equilib_mvrpmd::equilib_mvrpmd(int my_id, int root_proc, int num_procs,
                               std::string root_path)
    :my_id(my_id),
     root_proc(root_proc),
     num_procs(num_procs),
     myRand(time(NULL) + my_id),
     helper(root_path,my_id,num_procs,root_proc),
     root_path(root_path)
{
    sys_set = false;
    files_set = false;
}
int equilib_mvrpmd::run(double nuc_ss, double x_ss, double p_ss,
                         unsigned long long num_steps,unsigned long long stride){
    
    if (!sys_set || !files_set) {
        if (my_id==root_proc) {
            std::cout << "ERROR: equilibrium mvrpmd variables not set!" << std::endl;
            std::cout << "Aborting calculation." << std::endl;
        }
        return -1;
    }

    //Declare vectors for monte carlo moves
    vector<double> Q(nuc_beads,0);
    matrix<double> x(elec_beads,num_states,0), p(elec_beads,num_states,0);

    //Initialize vectors
    gen_initQ(Q,nuc_beads,nuc_ss);
    gen_initElec(x,elec_beads,num_states,x_ss);
    gen_initElec(p,elec_beads,num_states,p_ss);
    
    //Over-write with saved vectors if requested
    if(readPSV){helper.read_PSV(nuc_beads, elec_beads, num_states, Q, x, p);}
    
    /* Assemble Hamiltonian and Estimator*/
    //SpringEnergy V_spring(nuc_beads,mass,beta/nuc_beads);
    StateIndepPot V0(nuc_beads,mass);
    GTerm G(elec_beads,num_states,alpha);
    C_Matrix C(elec_beads,num_states,alpha);
    M_Matrix M(num_states,1,beta/elec_beads);
    theta_Esplit theta(num_states,elec_beads,C,M);
    theta_Esplit_dBeta dtheta(elec_beads,num_states,beta/elec_beads,C,M);
    
    mvrpmd_Esplit_ham H(beta,V0,G,theta);
    mvrpmd_Esplit_esti Esti(1,beta,V0,theta,dtheta);
    
//    M_Matrix M2(num_states,nuc_beads,beta/elec_beads);
//    M_Matrix_MTS M_MTS(nuc_beads,elec_beads,num_states,M2);
//    Theta_MTS thetaMTS(num_states,elec_beads,C,M_MTS);
//    dTheta_MTS_dBeta thetaMTS_dBeta(nuc_beads,elec_beads,num_states,beta/elec_beads,
//                                    C,M2,M_MTS);
//
//    MVRPMD_MTS_Hamiltonian H2(beta/nuc_beads,V_spring,V0,G,thetaMTS);
//    MVRPMD_MTS_Estimator Esti2(nuc_beads,beta/nuc_beads,V_spring,V0,thetaMTS,
//                         thetaMTS_dBeta);
    
    /* Initialize energy and estimator variables*/
    int esti_samples = 0;
    double estimator_total(0), sgn_total(0);
    double energy(0), estimator(0), sgnTheta(0);
    vector<double> estimator_t(num_steps/stride,0);

    energy = H.get_energy(Q,x,p);
    estimator = Esti.get_estimator();

    /* Over-write estimator with saved value if requested*/
    if (readData) {helper.read_MC_data(sgn_total, estimator_total);}

     //Initialize nuclear stepping procedure
    system_step nuc_stepper(my_id,num_procs,root_proc,nuc_beads,beta);
    nuc_stepper.set_nuc_ss(nuc_ss);
    nuc_stepper.set_hamiltonian(H);

     //Initialize electronic stepping procedure
    elec_step elec_stepper(my_id,num_procs,root_proc,elec_beads,num_states,beta);
    elec_stepper.set_hamiltonian(H);
    elec_stepper.set_energy(energy);
    elec_stepper.set_ss(x_ss,p_ss);

    //r is a ratio that determines how often to sample each sub-system
    double r = double(nuc_beads)/ double(nuc_beads + num_states*elec_beads);
    

    std::ofstream progress;
    int ten_p = floor(num_steps/10.0); //ten percent of steps to take
    double get_prog = true;
    
    if (ten_p == 0) {
        get_prog = false;
        if (my_id==root_proc) {
            std::cout << "Warning: equilibrium progress will not be tracked"
            " because simulation uses small number of steps." <<std::endl;
        }
    }
    
    if (my_id == root_proc) {
        std::string file_name = root_path + "Output/equil_progress";
        progress.open(file_name.c_str());
    }

    for (int step=0; step<num_steps; step++) {
        if (myRand.doub() < r) {
            /* Sample Nuclear Coordinates*/
            nuc_stepper.step(energy,Q,x,p);
            energy = nuc_stepper.get_energy();
        }
        else{
            /* Sample Electronic Coordinates*/
            if (myRand.doub() >= 0.5) {
                /* Sample x*/
                elec_stepper.step_x(energy,Q,x,p);
                energy = elec_stepper.get_energy();
            }
            else{
                /* Sample p*/
                elec_stepper.step_p(energy,Q,x,p);
                energy = elec_stepper.get_energy();
            }
        }
        
        estimator = Esti.get_estimator(Q,x,p);
        sgnTheta = theta.get_signTheta();
        
        if (is_NaN(estimator)) {
            std::cout << "Bad value for estimator encountered on" << std::endl;
            std::cout << "process " << my_id << " at step " << step << std::endl;
        }
        else{
            estimator_total += estimator;
        }
        
        if (is_NaN(sgnTheta)) {
            std::cout << "Bad value for sgnTheta encountered on" << std::endl;
            std::cout << "process " << my_id << " at step " << step << std::endl;
        }
        else{
            sgn_total += sgnTheta;
        }
        
        if(step % stride == 0){
            estimator_t[esti_samples] = estimator_total/(sgn_total + 1);
            esti_samples += 1;
        }
        if (get_prog) {
            if (step % ten_p == 0) {
                if (my_id==root_proc) {
                    progress << 100 * step/double(num_steps) << "%" << std::endl;
                }
            }
        }
    }
    
    if (my_id == root_proc) {
        progress.close();
    }

    /* Retrive acceptance ratio information */
    unsigned long long nuc_steps_total = nuc_stepper.get_steps_total();
    unsigned long long nuc_steps_accpt = nuc_stepper.get_steps_accepted();

    unsigned long long x_steps_accpt = elec_stepper.get_x_steps_accpt();
    unsigned long long x_steps_total = elec_stepper.get_x_steps_total();

    unsigned long long p_steps_accpt = elec_stepper.get_p_steps_accpt();
    unsigned long long p_steps_total = elec_stepper.get_p_steps_total();

    /* Print Monte Carlo simulation Data */
    helper.set_sys_ratio(nuc_steps_total, nuc_steps_accpt);
    helper.set_x_ratio(x_steps_total, x_steps_accpt);
    helper.set_p_ratio(p_steps_total, p_steps_accpt);
    helper.set_average_energy(estimator_total, sgn_total);
    helper.write_estimator(estimator_t,stride);
    
    helper.final_report(nuc_beads,elec_beads,num_states,beta,num_steps,nuc_ss,
                        x_ss,p_ss);

    /* Write Monte Carlo information to file if requested*/
    if(writePSV){helper.write_PSV(nuc_beads, elec_beads, num_states, Q, x, p);}
    if (writeData) {helper.write_MC_data(sgn_total, estimator_total);}
    
    return 0;
}
void equilib_mvrpmd::gen_initQ(vector<double> &Q, int num_beads, double step_size){
    for (int bead=0; bead<num_beads; bead++) {
        Q(bead) = step_dist(myRand.doub(),step_size);
    }
}
void equilib_mvrpmd::gen_initElec(matrix<double> &v, int num_beads, int num_states,
                                   double step_size){
    for (int bead=0; bead<num_beads; bead++) {
        for (int state=0; state<num_states; state++) {
            v(bead,state) = step_dist(myRand.doub(),step_size);
        }
    }
}
void equilib_mvrpmd::initialize_system(int nuc_beads_IN,int elec_beadsIN,
                                       int num_statesIN,double massIN,double betaIN,
                                       double alphaIN){
    nuc_beads = nuc_beads_IN;
    elec_beads = elec_beadsIN;
    num_states = num_statesIN;
    mass = massIN;
    beta = betaIN;
    alpha = alphaIN;
    sys_set = true;
}
void equilib_mvrpmd::initialize_files(bool writePSV_IN, bool readPSV_IN,
                                      bool writeData_IN, bool readData_IN){
    writePSV = writePSV_IN;
    readPSV = readPSV_IN;
    writeData = writeData_IN;
    readData = readData_IN;
    files_set = true;
}
inline double equilib_mvrpmd::step_dist(const double rn, double step_size){
    return (rn * 2.0 * step_size) - step_size;
}
void equilib_mvrpmd::set_write_PSV(bool set_In){writePSV = set_In;}

void equilib_mvrpmd::set_read_PSV(bool set_In){readPSV = set_In;}

void equilib_mvrpmd::set_read_Data(bool set_In){readData = set_In;}

void equilib_mvrpmd::set_write_Data(bool set_In){writeData = set_In;}
