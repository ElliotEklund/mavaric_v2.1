#ifndef rpmd_ham_hpp
#define rpmd_ham_hpp

#include "SpringEnergy.h"
#include "StateIndepPot.h"

class rpmd_ham{

public:
    rpmd_ham(int nuc_beads, double beta_num_beads,
               SpringEnergy &V_springIn,StateIndepPot &V0_In);

/* Mutators */
    double get_energy(const vector<double> &Q);
    double get_energy();
    double get_energy_dyn(double mass,const vector<double> &Q, const vector<double> &P);

private:
/* Data*/
    int nuc_beads; //nuclear  ring polymer beads
    double one_beta_num_beads; // 1/beta_num_beads

/* Objects */
    SpringEnergy * V_spring;
    StateIndepPot * V0;
};

#endif
