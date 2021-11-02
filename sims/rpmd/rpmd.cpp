#include <iostream>
#include <fstream>
#include "mpi.h"
#include <vector>

#include "input_rpmd.hpp"
#include "equilib_rpmd.hpp"
#include "sampling_rpmd.hpp"
#include "dynamics_rpmd.hpp"

int main(int argc, char *argv[]) {

    int num_procs = 1; //number of processors program is distributed over
    int my_id = 0; //unique id of each processor
    int root_process = 0; //processor 0 is default root_pathprocess

    MPI_Status status;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD,&my_id);
    MPI_Comm_size(MPI_COMM_WORLD,&num_procs);
    MPI_Comm comm = MPI_COMM_WORLD;

    std::string root = argv[1];
    std::string sim = argv[2];
    std::string root_path = root + "sims/" + sim + "/";

                        /* BEGIN PROCESS 1 */
    /* This process reads in all parameters stored in
     InputFiles and distributes them to their appropriate
     variables.*/
    /* Vectors used to store parameters from InputFiles */
    std::vector<double> sys_parameters;
    std::vector<double> MC_parameters;
    std::vector<double> Samp_parameters;
    std::vector<double> Dyn_parameters;

    input_rpmd myInput;
    int abort = myInput.input_file_handler(root_path,sys_parameters,MC_parameters,
                                           Samp_parameters,Dyn_parameters);

    if (abort == -1){return -1;}

    /* Physical parameters.*/
    double temp, mass;
    int nuc_beads;
    double beta, beta_nuc;

    /* Monte Carlo parameters. */
    unsigned long long num_steps, esti_rate;
    double nuc_ss;
    bool writePSV, readPSV;
    bool readData, writeData;
    bool runMC;

    /* Sampling parameters.*/
    bool runSamp, saveTrajs;
    unsigned long long num_trajs;
    unsigned long long decor_len;

    /* Dynamics Variables */
    bool runDyn;
    double dt;
    double total_time;
    double tol;
    bool run_ac;
    int energy_stride;
    bool run_energ_conserv;
    double interval;

    /* From MonteCarlo */
    runMC = MC_parameters[0];
    num_steps = MC_parameters[1];
    esti_rate = MC_parameters[2];
    writePSV = MC_parameters[3];
    writeData = MC_parameters[4];
    readPSV = MC_parameters[5];
    readData = MC_parameters[6];

    /* From SystemParameters*/
    mass = sys_parameters[0];
    nuc_beads = sys_parameters[1];
    temp = sys_parameters[2];
    nuc_ss = sys_parameters[3];

    /* From Sampling */
    runSamp = Samp_parameters[0];
    num_trajs = Samp_parameters[1];
    decor_len = Samp_parameters[2];
    saveTrajs = Samp_parameters[3];

    /* From Dynamics */
    runDyn = Dyn_parameters[0];
    total_time = Dyn_parameters[1];
    dt = Dyn_parameters[2];
    double junk = Dyn_parameters[3];
    run_ac = Dyn_parameters[4];
    run_energ_conserv = Dyn_parameters[5];
    tol = Dyn_parameters[6];
    interval = Dyn_parameters[7];

    /* Set derived variables */
    beta = 1.0/temp;
    beta_nuc = beta/nuc_beads;

                        /* END PROCESS 1 */
    /* /////////////////////////////////////////////////////////*/

                          /* BEGIN PROCESS 2 */
      /* This process runs the Monte Carlo simulation if requested.*/
    if (runMC) {
        if (my_id == root_process) {
            std::cout << std::endl << std::endl;
            std::cout << "Begin Monte Carlo Simulation" << std::endl;
            std::cout << std::endl;
        }
        equilib_rpmd equilibrator(my_id,root_process,num_procs,root_path);
        equilibrator.initialize_system(nuc_beads,mass,beta);
        equilibrator.initialize_files(writePSV,readPSV,writeData,readData);

        clock_t start = clock();
        equilibrator.run(nuc_ss,num_steps,esti_rate);
        clock_t end = clock();
        double time_taken = double(end - start) / double(CLOCKS_PER_SEC);

        if (my_id == root_process) {
            std::cout << "\t Monte Carlo simulation time: " << time_taken <<
            std::endl << std::endl;
            std::cout << "End Monte Carlo Simulation" << std::endl;
            std::cout << std::endl;
        }
    }
                           /* END PROCESS 2 */
    /* /////////////////////////////////////////////////////////// */

    /* /////////////////////////////////////////////////////////// */
                          /* BEGIN PROCESS 3 */
      /* This process runs Sampling  if requested.*/
    if(runSamp){
        if (my_id == root_process) {
            std::cout << "Begin Sampling Simulation" << std::endl;
            std::cout << std::endl;
        }
        sampling_rpmd sampler(my_id,root_process,num_procs);
        sampler.initialize_system(nuc_beads,mass,beta);
        sampler.initialize_files(readPSV,saveTrajs,root_path);

        clock_t start = clock();
        sampler.run(nuc_ss,num_trajs,decor_len);
        clock_t end = clock();
        double time_taken = double(end - start) / double(CLOCKS_PER_SEC);

        if (my_id == root_process) {
            std::cout << "\t Sampling simulation time: " << time_taken <<
            std::endl << std::endl;

            std::cout << "End Sampling Simulation" << std::endl;
            std::cout << std::endl << std::endl;
        }
    }
                            /* END PROCESS 3 */
  /* /////////////////////////////////////////////////////////// */

    /* /////////////////////////////////////////////////////////// */
                          /* BEGIN PROCESS 4 */
    /* This process runs the Dynamics simulation if requested.*/
    if(runDyn){
        if (my_id == root_process) {
            std::cout << "Begin Dynamics Simulation" << std::endl;
            std::cout << std::endl << std::endl;
        }

        dynamics_rpmd my_dynamics(my_id,num_procs,root_process);
        my_dynamics.set_system(nuc_beads,mass,beta,beta_nuc);
        my_dynamics.set_time(dt,total_time);
        my_dynamics.set_trajs(num_trajs,root_path);

        std::string input_dir = root_path;
        std::string output_dir = root_path;

        clock_t start = clock();
        if(run_ac){
          std::cout << "Running auto-correlation function calculation." << std::endl;
          my_dynamics.compute_ac(interval,input_dir,output_dir);
        }
        if(run_energ_conserv){
          std::cout << "Running energy conservation calculation." << std::endl;
          my_dynamics.energy_conserve(tol,interval,input_dir,output_dir);
        }
        clock_t end = clock();
        double time_taken = double(end - start) / double(CLOCKS_PER_SEC);

        if (my_id == root_process) {
            std::cout << "Dynamics simulation time:" << time_taken << std::endl;
            std::cout << "End Dynamics Simulation" << std::endl;
        }
    }
                            /* END PROCESS 4 */
    /* /////////////////////////////////////////////////////////// */
    MPI_Finalize();
    return 0;
}
