cmake_minimum_required(VERSION 3.15.4)

set(MACHINE "ASTRA")

if(${MACHINE} STREQUAL "MAC")
  set(ROOT "/Users/ellioteklund/Desktop/MAVARIC_v2.0/MAVARIC")
  set(BOOST "/Users/ellioteklund/Downloads")
elseif(${MACHINE} STREQUAL "ASTRA")
  set(ROOT "/home/fs01/ece52/MAVARIC-MTS/MAVARIC")
  set(BOOST_DIR $ENV{BOOST_INC})
endif()

# set the project name
project(mavaric)
set(CMAKE_CXX_STANDARD 98)
set(CMAKE_CXX_COMPILER mpic++)

#if(${MACHINE} STREQUAL "MAC")
    set(EQUILIBRIUM_DIR "${ROOT}/sims/rpmd/equilibrium")
    set(SAMPLING_DIR "${ROOT}/sims/rpmd/sampling")
    set(DYNAMICS_DIR "${ROOT}/sims/rpmd/dynamics")
    set(HAMILTONIAN_DIR "${ROOT}/src/hamiltonian")
    set(ESTIMATORS_DIR "${ROOT}/src/estimators")
    set(FORCES_DIR "${ROOT}/src/forces")
    set(INTEGRATORS_DIR "${ROOT}/src/integrators")
    set(MISC_DIR "${ROOT}/src/misc")
    set(MONTECARLO_DIR "${ROOT}/src/monte_carlo")
    set(IO_DIR "${ROOT}/src/IO")
    set(EXECUTABLE_OUTPUT_PATH "${ROOT}/bin")
    set(CMAKE_CXX_FLAGS "-DMAC")

add_definitions( "-DNDEBUG -O3" )

# add the executable

add_executable(rpmd rpmd.cpp

    #equilibrium directory
        ${EQUILIBRIUM_DIR}/MonteCarloHelper.cpp
        ${EQUILIBRIUM_DIR}/equilib_rpmd.cpp
        ${EQUILIBRIUM_DIR}/input_rpmd.cpp

        #sampling directory
        ${SAMPLING_DIR}/sampling_rpmd.cpp
        ${SAMPLING_DIR}/SamplingHelper.cpp

        #dynamics directory
        ${DYNAMICS_DIR}/dynamics_rpmd.cpp
        ${DYNAMICS_DIR}/rpmd_auto_corr.cpp
        ${DYNAMICS_DIR}/rpmd_energy_conserv.cpp

        #hamiltonian directory
        ${HAMILTONIAN_DIR}/rpmd_ham.cpp
        ${HAMILTONIAN_DIR}/SpringEnergy.cpp
        ${HAMILTONIAN_DIR}/StateIndepPot.cpp

        #estimators directory
        ${ESTIMATORS_DIR}/rpmd_estimator.cpp

        #forces directory
        ${FORCES_DIR}/dSpring_dQ.cpp
        ${FORCES_DIR}/dStateIndep_dQ.cpp
        ${FORCES_DIR}/rpmd_force.cpp

        #integrators directory
        ${INTEGRATORS_DIR}/rpmd_vv.cpp

        #misc directory
        ${MISC_DIR}/mpi_wrapper.cpp

        #monte carlo directory
        ${MONTECARLO_DIR}/rpmd_system_step.cpp
    )

include_directories(${BOOST_DIR})
include_directories(${EQUILIBRIUM_DIR}/Headers)
include_directories(${SAMPLING_DIR}/Headers)
include_directories(${DYNAMICS_DIR}/Headers)
include_directories(${HAMILTONIAN_DIR}/Headers)
include_directories(${ESTIMATORS_DIR}/Headers)
include_directories(${FORCES_DIR}/Headers)
include_directories(${INTEGRATORS_DIR}/Headers)
include_directories(${MISC_DIR}/Headers)
include_directories(${MONTECARLO_DIR}/Headers)
include_directories(${IO_DIR}/Headers)

