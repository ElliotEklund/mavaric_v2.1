#ifndef equilib_rpmd_hpp
#define equilib_rpmd_hpp

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/storage.hpp>

#include <stdlib.h>
#include <time.h>
#include <string>

#include "nr3.h"
#include "ran.h"
#include "mpi.h"

#include "mpi_wrapper.hpp"

#include "MonteCarloHelper.h"
#include "rpmd_ham.hpp"
#include "rpmd_estimator.hpp"
#include "rpmd_system_step.hpp"
// #include "functions.hpp"

class equilib_rpmd{

public:
    equilib_rpmd(int my_id, int root_proc, int num_procs, std::string root_path);

    /* Runs Monte Carlo simulation.
     nuc_ss: nuclear step size used by nuclear part of system; used for proposing
              monte carlo steps
     x_ss: x mapping variable step size used by x part of system; used for proposing
              monte carlo steps
     p_ss: p mapping variable step size used by p part of system; used for proposing
              monte carlo steps
     num_steps: number of monte carlo steps simulation will take
     stride: length between successive collection of energy estimator
     */
    int run(double nuc_ss, unsigned long long num_steps, unsigned long long stride);

// /* Mutators */
    void initialize_system(int nuc_beads_IN,double massIN,double betaIN);

    void initialize_files(bool writePSV_IN, bool readPSV_IN,
                          bool writeData_IN, bool readData_IN);

    void set_write_PSV(bool set_In);
    void set_read_PSV(bool set_In);
    void set_read_Data(bool set_In);
    void set_write_Data(bool set_In);

private:

/* Data */
    int my_id, num_procs, root_proc; //mpi data
    std::string root_path;
    bool writePSV, readPSV; // determine how Phase Spave Variables are processed
    bool writeData, readData; // determine how MC data is processed

    /* Number of nuclear beads*/
    int nuc_beads;
    double mass, beta; //system mass; 1/kb T
    bool sys_set, files_set; //true only once initialize_system(initialize_files) called

/* Objects */
    Ran myRand; //NR3 random number generator
    MonteCarloHelper helper;//helper object

/* Functions */
    /* Initialize Q with random vaules.
     Q: vector of bead positions
     num_beads: number of beads in Q
     step_size: sets range of individual elements of Q; Q[i] ~ [-step_size,step_size]*/
    void gen_initQ(vector<double> &Q, int num_beads, double step_size);

    /* Returns a uniform random number between [-step_size,step_size]
     rn: random double
     step_size: specifices range of uniform random distribution*/
    inline double step_dist(const double rn, double step_size);
};

#endif
