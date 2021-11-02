#include "dynamics_rpmd.hpp"

dynamics_rpmd::dynamics_rpmd(int my_id, int num_procs, int root_proc)
    :my_id(my_id),
     num_procs(num_procs),
     root_proc(root_proc)
{
    is_sys_set = false;
    is_time_set = false;
    is_trajs_set = false;
}
int dynamics_rpmd::pre_comp(){
    int abort = 0;

    if (!is_sys_set) {
        abort = -1;
        if (my_id==root_proc) {
            std::cout << "ERROR: System variables are not set!" << std::endl;
            std::cout << "Aborting calculation." << std::endl;
        }
    }
    if (!is_time_set) {
        abort = -1;
        if (my_id==root_proc) {
            std::cout << "ERROR: Time variables are not set!" << std::endl;
            std::cout << "Aborting calculation." << std::endl;
        }
    }
    if (!is_trajs_set) {
        abort = -1;
        if (my_id==root_proc) {
            std::cout << "ERROR: Trajectory variables are not set!" << std::endl;
            std::cout << "Aborting calculation." << std::endl;
        }
    }
    return abort;
}
int dynamics_rpmd::compute_ac(int pac_stride,std::string input_dir,
                                 std::string output_dir){

    if(pre_comp()==-1){return -1;}
    rpmd_auto_corr my_ac(my_id,num_procs,root_proc,root_path);
    my_ac.compute_position_auto_corr(num_trajs_global,nuc_beads,dt,total_time,mass,
                                    beta,pac_stride);
    return 0;
}
int dynamics_rpmd::energy_conserve(double tol, int energy_stride,
                                      std::string input_dir,
                                      std::string output_dir){

    if(pre_comp()==-1){return -1;}
    rpmd_energy_conserv my_conserv(my_id,num_procs,root_proc,root_path);
    my_conserv.compute_energy(num_trajs_global,nuc_beads,dt,total_time,mass,beta,
                                                   energy_stride,tol);

    return 0;
}
void dynamics_rpmd::set_system(int nuc_beadsIN, double massIN,double betaIN,
                              double beta_nuc_beadsIN){
    nuc_beads = nuc_beadsIN;
    mass = massIN;
    beta = betaIN;
    beta_nuc_beads = beta_nuc_beadsIN;
    is_sys_set = true;
}
void dynamics_rpmd::set_time(double dtIN, double total_timeIN){
    dt = dtIN;
    total_time = total_timeIN;
    is_time_set = true;
}
void dynamics_rpmd::set_trajs(unsigned long long num_trajs_globalIN,
                                std::string root_pathIN){

    num_trajs_global = num_trajs_globalIN;
    num_trajs_local = 0; //number of local trajectories
    root_path = root_pathIN;

    /* Check for more trajs than procs*/
    if (num_trajs_global < num_procs) {
        if (my_id == root_proc) {
            std::cout << "ERROR: num_trajs is smaller than num_procs.!" << std::endl;
        }
    }
    else{
        unsigned long long d = num_trajs_global/num_procs;
        unsigned long long r = num_trajs_global%num_procs;

        if (my_id < num_procs - r) {
            num_trajs_local = d;
        }
        else{
            num_trajs_local = d + 1;
        }
        is_trajs_set = true;
    }
}
