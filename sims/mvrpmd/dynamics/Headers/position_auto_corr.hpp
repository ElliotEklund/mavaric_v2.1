#ifndef position_auto_corr_hpp
#define position_auto_corr_hpp

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>

#include <math.h>
#include <sstream>
#include <list>
#include <fstream>
#include <string>
#include "mpi.h"

#include "M_Matrix.h"
#include "C_Matrix.h"
#include "dM_Matrix_dQ.hpp"
#include "theta_mixed.hpp"
#include "theta_mixed_dQ.hpp"
#include "theta_mixed_dElec.hpp"

#include "Forces_MTS.hpp"
#include "mvrpmd_mixed_forces.hpp"
#include "ABM_MVRPMD.hpp"
#include "trajs_io.hpp"

class position_auto_corr{

public:
    position_auto_corr(int my_id, int root_proc, int num_procs);

    void set_system(int nuc_beadsIN,int elec_beadsIN,int num_statesIN,
        double massIN, double beta_nuc_beadsIN,double beta_elec_beadsIN,
        double alphaIN);

    void set_time(double dtIN, double total_timeIN);

    void compute(unsigned long long num_trajs,std::string input_dir,
        std::string output_dir,int interval);

    double centroid(const vector<double> & Q);

    void write_data(const matrix<double> &cQQ_final,const matrix<double> &sign_final,
      std::string output_dir,int num_trajs,int interval,int num_samples);


private:

  int nuc_beads, elec_beads, num_states; //system information
  double mass, beta_nuc_beads, beta_elec_beads; //system information
  double alpha; //mapping variable pre-factor
  bool is_sys_set;

  double dt;
  double total_time;
  bool is_time_set;

  int my_id, root_proc, num_procs;

};

#endif
