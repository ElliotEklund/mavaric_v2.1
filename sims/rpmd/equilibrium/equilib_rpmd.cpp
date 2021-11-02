#include "equilib_rpmd.hpp"

equilib_rpmd::equilib_rpmd(int my_id, int root_proc, int num_procs,
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
int equilib_rpmd::run(double nuc_ss, unsigned long long num_steps,
                      unsigned long long stride){

    if (!sys_set || !files_set) {
        if (my_id==root_proc) {
            std::cout << "ERROR: equilibrium rpmd variables not set!" << std::endl;
            std::cout << "Aborting calculation." << std::endl;
        }
        return -1;
    }

    //Declare and intialize vectors for monte carlo moves
    vector<double> Q(nuc_beads,0);
    gen_initQ(Q,nuc_beads,nuc_ss);

    //Over-write with saved vectors if requested
    if(readPSV){helper.read_PSV(nuc_beads, Q);}

    /* Assemble Hamiltonian and Estimator*/
    SpringEnergy V_spring(nuc_beads,mass,beta/nuc_beads);
    StateIndepPot V0(nuc_beads,mass);
    rpmd_ham H(nuc_beads,beta/nuc_beads,V_spring,V0);
    rpmd_estimator Esti(nuc_beads,beta/nuc_beads,V_spring,V0);

    /* Initialize energy and estimator variables*/
    int esti_samples = 0;
    double estimator_total(0);
    double energy(0), estimator(0);
    vector<double> estimator_t(num_steps/stride,0);
    energy = H.get_energy(Q);
    estimator = Esti.get_estimator();

    /* Over-write estimator with saved value if requested*/
    if (readData) {helper.read_MC_data(estimator_total);}

    //Initialize nuclear stepping procedure
    rpmd_system_step nuc_stepper(my_id,num_procs,root_proc,nuc_beads,beta);
    nuc_stepper.set_nuc_ss(nuc_ss);
    nuc_stepper.set_hamiltonian(H);

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
        /* Sample Nuclear Coordinates*/
        nuc_stepper.step(energy,Q);
        energy = nuc_stepper.get_energy();
        estimator = Esti.get_estimator(Q);
        estimator_total += estimator;

        if(step % stride == 0){
            estimator_t[esti_samples] = estimator_total;
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

    /* Print Monte Carlo simulation Data */
    helper.set_sys_ratio(nuc_steps_total, nuc_steps_accpt);
    helper.set_average_energy(estimator_total, num_steps);
    helper.final_report(nuc_beads,beta,num_steps,nuc_ss);

    /* Write Monte Carlo information to file if requested*/
    if(writePSV){helper.write_PSV(nuc_beads,Q);}
    if (writeData) {
      helper.write_MC_data(estimator_total);
      helper.write_estimator(estimator_t,stride);
    }

    return 0;
}
void equilib_rpmd::gen_initQ(vector<double> &Q, int num_beads, double step_size){
    for (int bead=0; bead<num_beads; bead++) {
        Q(bead) = step_dist(myRand.doub(),step_size);
    }
}
void equilib_rpmd::initialize_system(int nuc_beads_IN,double massIN,double betaIN){
    nuc_beads = nuc_beads_IN;
    mass = massIN;
    beta = betaIN;
    sys_set = true;
}
void equilib_rpmd::initialize_files(bool writePSV_IN, bool readPSV_IN,
                                      bool writeData_IN, bool readData_IN){
    writePSV = writePSV_IN;
    readPSV = readPSV_IN;
    writeData = writeData_IN;
    readData = readData_IN;
    files_set = true;
}
inline double equilib_rpmd::step_dist(const double rn, double step_size){
    return (rn * 2.0 * step_size) - step_size;
}
void equilib_rpmd::set_write_PSV(bool set_In){writePSV = set_In;}
void equilib_rpmd::set_read_PSV(bool set_In){readPSV = set_In;}
void equilib_rpmd::set_read_Data(bool set_In){readData = set_In;}
void equilib_rpmd::set_write_Data(bool set_In){writeData = set_In;}
