#include "StateIndepPot.h"

StateIndepPot::StateIndepPot(int num_beads, double mass)
    :num_beads(num_beads), mass(mass),one_vec(num_beads,1.0),
    Q_2(num_beads,0),Q_3(num_beads,0)
{}
void StateIndepPot::update_V0(const vector<double> &Q){

  Q_2 = element_prod(Q,Q);
  // Q_3 = element_prod(Q_2,Q);

  double QQ = inner_prod(Q_2,one_vec);
  // double QQQ = inner_prod(Q_3,one_vec);
  // double QQQQ = inner_prod(Q_2,Q_2);

    energy = 0.5*QQ;
    // energy = 0.5*QQ + 0.1*QQQ + 0.01*QQQQ;
}
double& StateIndepPot::get_V0(const vector<double> &Q){

    update_V0(Q);
    return energy;
}
double & StateIndepPot::get_V0(){return energy;}
