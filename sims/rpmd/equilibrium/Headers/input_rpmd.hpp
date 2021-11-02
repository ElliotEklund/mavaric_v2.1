#ifndef INPUT_RPMD_H
#define INPUT_RPMD_H

#include <algorithm>
#include <fstream>
#include <iostream>
#include <math.h>
#include <string>
#include <stdlib.h>
#include <vector>

#include "system_in.h"
#include "monte_carlo_in.h"
#include "sampling_in.h"
#include "rpmd_dynamics_in.h"

class input_rpmd{

public:
    input_rpmd();

    /* Read in the parameters from the InputFile "fileName". Return
     * these parameters in a vector.*/
    std::vector<double> get_parameters(std::string fileName);

    /* Each of the five input vectors well be filled with their corresponding values after
     * input_file_handler is called. These vectors will only be returned if the function does
     * not find any corrupt files or errors.*/
    int input_file_handler(std::string root, std::vector<double> &sysParams,
                                  std::vector<double> &MCParams,
                                  std::vector<double> &sampParams,
                                  std::vector<double> &dynParams);

};

#endif
