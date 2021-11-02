#include "sampling_rpmd.hpp"

sampling_rpmd::sampling_rpmd(int my_id, int root_proc, int num_procs)
    :my_id(my_id),
     root_proc(root_proc),
     num_procs(num_procs),
     myRand(time(NULL) + my_id),
     helper(my_id,num_procs,root_proc)
{
    sys_set = false;
    files_set = false;
}
int sampling_rpmd::run(double nuc_ss, unsigned long long num_trajs,
                        unsigned long long decorr){

    if (!sys_set || !files_set) {
        if (my_id==root_proc) {
            std::cout << "ERROR: sampling mvrpmd variables not set!" << std::endl;
            std::cout << "Aborting calculation." << std::endl;
        }
        return -1;
    }

    unsigned long long num_trajs_local = 0;
    num_trajs_local = get_trajs_local(num_trajs);

    /* Declare vectors for monte carlo moves*/
    vector<double> Q(nuc_beads,0);
    gen_initQ(Q,nuc_beads,nuc_ss);

    if (readPSV) {
        helper.read_PSV(nuc_beads,Q);
    }

    /* Declare vectors to store sampled trajectories*/
    vector<double> Q_trajs(num_trajs_local*nuc_beads,0);
    vector<double> P_trajs(num_trajs_local*nuc_beads,0);

    /* Assemble Hamiltonian*/
    SpringEnergy V_spring(nuc_beads,mass,beta/nuc_beads);
    StateIndepPot V0(nuc_beads,mass);
    rpmd_ham H(nuc_beads,beta/nuc_beads,V_spring,V0);
    double energy = H.get_energy(Q);
    double energy_prop = energy;

    /* Initialize nuclear stepping procedure */
    rpmd_system_step nuc_stepper(my_id,num_procs,root_proc,nuc_beads,beta);
    nuc_stepper.set_nuc_ss(nuc_ss);
    nuc_stepper.set_hamiltonian(H);
    nuc_stepper.set_energy(energy);

    int ten_p = floor(num_trajs_local/10.0); //ten percent of trajectories
    bool get_report = true;

    if (ten_p == 0) {
        if (my_id==root_proc) {
            std::cout << "WARNING: Sampling progress report will not be generated "
            "because num_trajs_local < 10" << std::endl;
            get_report = false;
        }
    }

    std::ofstream progress;
    if (my_id == root_proc) {
        std::string my_file = rootFolder + "Output/samp_progress";
        progress.open(my_file.c_str());
        if (!progress.is_open()) {
            std::cout << "ERROR: could not open file " << my_file << std::endl;
        }
    }

    /* Main Algorithm*/
    for (int traj=0; traj<num_trajs_local; traj++) {
        for (int step=0; step<decorr; step++) {
          /* Sample Nuclear Coordinates*/
          nuc_stepper.step(energy,Q);
          energy = nuc_stepper.get_energy();
        }
        for (int bead=0; bead<nuc_beads; bead++) {
            Q_trajs(traj*nuc_beads+bead) = Q(bead);
        }
        if (my_id == root_proc && get_report) {
            if (traj % ten_p == 0) {
                progress <<  100 * double(traj)/num_trajs_local << "%" << std::endl;
            }
        }
    }

    if (my_id==root_proc) {progress.close();}

    double stdev = sqrt(nuc_beads*mass/beta);
    //double stdev = sqrt(mass/(beta*nuc_beads));

    //Sample momentum distribution from Gaussian(mu, sigma, seed)
    Normaldev_BM momentum(0, stdev, myRand.int32());

    for(int i=0; i<num_trajs_local*nuc_beads; i++){
        P_trajs(i) = momentum.dev();
    }

    if (saveTrajs) {
      save_trajs(Q_trajs,"Output/Trajectories/Q",nuc_beads,beta,num_trajs,decorr);
      save_trajs(P_trajs,"Output/Trajectories/P",nuc_beads,beta,num_trajs,decorr);
    }

    unsigned long long nuc_steps_total = nuc_stepper.get_steps_total();
    unsigned long long nuc_steps_accpt = nuc_stepper.get_steps_accepted();
    helper.print_sys_accpt(nuc_steps_total,nuc_steps_accpt,"Nuclear");

    return 0;
}
void sampling_rpmd::gen_initQ(vector<double> &Q, int num_beads, double step_size){
    for (int bead=0; bead<num_beads; bead++) {
        Q(bead) = step_dist(myRand.doub(),step_size);
    }
}
inline double sampling_rpmd::step_dist(const double rn, double step_size){
    return (rn * 2.0 * step_size) - step_size;
}
void sampling_rpmd::initialize_system(int nuc_beads_IN,double massIN,double betaIN){
    nuc_beads = nuc_beads_IN;
    mass = massIN;
    beta = betaIN;
    sys_set = true;
}
void sampling_rpmd::initialize_files(bool readPSVIN,bool saveTrajsIN,
                                       std::string rootFolderIN){

    readPSV = readPSVIN;
    rootFolder = rootFolderIN;
    // helper.set_root(rootFolderIN);
    rootFolder = rootFolderIN;
    saveTrajs = saveTrajsIN;
    files_set = true;
}
void sampling_rpmd::save_trajs(vector<double> &v,std::string name, int nuc_beadsIN,
                                double betaIN,unsigned long long num_trajs_totalIN,
                                double decorrIN){

    std::string fileName = rootFolder + name;
    mpi_wrapper myWrap(num_procs,my_id,root_proc);
    myWrap.write_vector(v,fileName,nuc_beadsIN,betaIN,num_trajs_totalIN,decorrIN);
}
unsigned long long sampling_rpmd::get_trajs_local(unsigned long long num_trajs_totalIN){
    unsigned long long num_trajs_local = 0; //number of local trajectories
    /* Check for more trajs than procs*/
    if (num_trajs_totalIN < num_procs) {
        if (my_id == root_proc) {
            std::cout << "ERROR: num_trajs is small than num_procs.!" << std::endl;
        }
        return 0;
    }
    else{
        unsigned long long d = num_trajs_totalIN/num_procs;
        unsigned long long r = num_trajs_totalIN%num_procs;

        if (my_id < num_procs - r) {
            num_trajs_local = d;
        }
        else{
            num_trajs_local = d + 1;
        }
        return num_trajs_local;
    }
}
