#include "dStateIndep_dQ.hpp"

dStateIndep_dQ::dStateIndep_dQ(int nuc_beads,double mass)
    :nuc_beads(nuc_beads), mass(mass),
     force(nuc_beads,0.0),QQ(nuc_beads,0.0),QQQ(nuc_beads,0.0)
{}

const vector<double> & dStateIndep_dQ::get_dStateIndep_dQ(const vector<double> &Q){

    // QQ = element_prod(Q,Q);
    // QQQ = element_prod(Q,QQ);

    // force = Q + 0.3*QQ + 0.04*QQQ;
    force = Q;

    return force;
}
