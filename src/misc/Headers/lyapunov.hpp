#ifndef lyapunov_hpp
#define lyapunov_hpp

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>

#include "SpringEnergy.h"
#include "StateIndepPot.h"
#include "GTerm.h"

#include "M_Matrix.h"
#include "C_Matrix.h"
#include "dM_Matrix_dQ.hpp"

#include "M_Matrix_MTS.hpp"
#include "dM_Matrix_dQ.hpp"

#include "theta_mixed.hpp"
#include "theta_mixed_dQ.hpp"
#include "theta_mixed_dElec.hpp"

#include "Forces_MTS.hpp"
#include "mvrpmd_mixed_forces.hpp"
#include "mvrpmd_mixed_ham.hpp"
#include "ABM_MVRPMD.hpp"

#include <fstream>
#include <string>
#include "mpi.h"
#include <math.h>
#include "trajs_io.hpp"
#include <list>
#include <algorithm>

using namespace boost::numeric::ublas;

class lyapunov{
    
public:
    
    lyapunov(int my_id, int num_procs, int root_proc);
    
    void compute();

    
private:
    
    int my_id, num_procs, root_proc; //mpi data
    int nuc_beads, elec_beads, num_states; //system data
    double mass, beta, beta_nuc_beads, beta_elec_beads; //system data
    double alpha; //mapping variable pre-factor
    
    void from_phase(vector<double> &Q,vector<double> &P,
                              matrix<double> &x,matrix<double> &p,
                              const vector<double> &v);
    
    vector<double> to_phase(const vector<double> &Q,const vector<double> &P,
                            const matrix<double> &x,const matrix<double> &p);
    
    vector<double> v_diff(const vector<double> &v1, const vector<double> &v2);

    double mag(const vector<double> &v);
    
    void set_system(int nuc_beadsIN, int elec_beadsIN, int num_statesIN,
                    double massIN,double betaIN, double beta_nuc_beadsIN,
                    double beta_elec_beadsIN, double alphaIN);
    
    vector<double> LCE(const vector<double> &d, double d0, double renorm, int num_cycles);
    
    void write_LCE(const vector <double> &k, int num_cycles, double stride_dt);

};

#endif
