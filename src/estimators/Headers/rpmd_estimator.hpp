#ifndef rpmd_estimator_hpp
#define rpmd_estimator_hpp

#include "SpringEnergy.h"
#include "StateIndepPot.h"

class rpmd_estimator{

public:

    rpmd_estimator(int num_beads, double beta_num_beads,
                     SpringEnergy &V_SpringIn, StateIndepPot &V0In);

    /* Return energy estimator. This function assumes that all
     parts of the Hamiltonian have been updated. */
    double get_estimator();

    double get_estimator(const vector<double> &Q);


private:

    /* Private data.*/
    double ONE_HALF_beta_num_beads; // 0.5 * (1/beta_num_beads)
    double ONE_num_beads; //1.0/num_beads

    SpringEnergy * V_spring;
    StateIndepPot * V0;
};

#endif
