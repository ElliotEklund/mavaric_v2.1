#include "rpmd_system_step.hpp"

rpmd_system_step::rpmd_system_step(int my_id,int num_procs,int root_proc,int num_beads,
                         double betaIN)
    :my_id(my_id),
     num_procs(num_procs),
     root_proc(root_proc),
     myRand(time(NULL) + my_id),
     num_beads(num_beads)
{
    Q_prop.resize(num_beads,0);
    set_beta(betaIN);
    steps_total = 0;
    steps_accepted = 0;
}
void rpmd_system_step::step(double energyIN, vector<double> &Q){

    Q_prop = Q;
    energy = energyIN;
    double mcMove=0;
    double energy_prop = 0;
    bool accept_move = false;

    for (int bead=0; bead<num_beads; bead++) {
        mcMove = rand_bead(myRand.int64(),num_beads);
        Q_prop(mcMove) = Q(mcMove) + step_dist(myRand.doub(),nuc_ss);
        energy_prop = H->get_energy(Q_prop);
        accept_move = check_move(energy,energy_prop);

        if (accept_move) {
            Q(mcMove) = Q_prop(mcMove);
            energy = energy_prop;
            steps_accepted += 1;
        }
        else{
            Q_prop(mcMove) = Q(mcMove);
        }
        steps_total += 1;
    }
}
inline double rpmd_system_step::step_dist(const double rn, double step_size){
    return (rn * 2.0 * step_size) - step_size;
}
inline int rpmd_system_step::rand_bead(const Ullong rn, int num_beads){
    return rn % num_beads;
}
bool rpmd_system_step::check_move(double energy, double energy_prop){

    double delta_energy = energy_prop - energy;
    double accept_move = false;
    /* Accept new system moves if energ_prop < energy*/
    if(delta_energy < 0){
        accept_move = true;
    }
    /* Accept new system moves if inequality is met*/
    else if (myRand.doub() <= exp(-beta * delta_energy/num_beads)){
        accept_move = true;
    }
    return accept_move;
}
void rpmd_system_step::set_hamiltonian(rpmd_ham &H_IN){H = &H_IN;}
void rpmd_system_step::set_nuc_ss(double nuc_ssIN){nuc_ss = nuc_ssIN;}
void rpmd_system_step::set_beta(double betaIN){beta = betaIN;}
void rpmd_system_step::set_energy(double energyIN){energy = energyIN;}
double rpmd_system_step::get_energy(){return energy;}
unsigned long long rpmd_system_step::get_steps_total(){return steps_total;}
unsigned long long rpmd_system_step::get_steps_accepted(){return steps_accepted;}
