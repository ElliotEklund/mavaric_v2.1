#ifndef rpmd_energy_conserv_hpp
#define rpmd_energy_conserv_hpp

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>

#include <string>
#include "mpi.h"
#include <fstream>
#include <stdlib.h>
#include "trajs_io.hpp"
#include "rpmd_vv.hpp"
#include "rpmd_ham.hpp"

using namespace boost::numeric::ublas;

class rpmd_energy_conserv{

public:

    rpmd_energy_conserv(int my_id, int num_procs, int root_procs,std::string root_path);

    int compute_energy(int num_trajs,int nuc_beads,double dt,
                                    double total_time,double mass, double beta,
                                    int stride,double tol);

private:

  int read_in_trajs(int num_trajs,int nuc_beads);
  double centroid(const vector<double> &Q,int nuc_beads);
  int write_ac_data(int num_trajs,int num_samples,const matrix<double> &ac_data);

  /* Data */
  int my_id, num_procs, root_proc; //mpi data
  std::string root_path;

  vector<double> P, Q;

  };
  #endif
