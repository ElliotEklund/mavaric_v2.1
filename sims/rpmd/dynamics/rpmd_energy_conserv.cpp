#include "rpmd_energy_conserv.hpp"

rpmd_energy_conserv::rpmd_energy_conserv(int my_id, int num_procs, int root_proc,
                               std::string root_path)
    :my_id(my_id),
     num_procs(num_procs),
     root_proc(root_proc),
     root_path(root_path)
{}
int rpmd_energy_conserv::compute_energy(int num_trajs,int nuc_beads,double dt,
                                              double total_time,double mass, double beta,
                                              int stride, double tol){

  read_in_trajs(num_trajs,nuc_beads);
  int num_steps = total_time/dt;
  int num_samples = num_steps/stride;

  vector<double> Qtraj(nuc_beads,0);
  vector<double> Ptraj(nuc_beads,0);

  rpmd_vv vv(nuc_beads,mass,beta,dt);
  SpringEnergy V_spring (nuc_beads,mass,beta/nuc_beads);
  StateIndepPot V0 (nuc_beads,mass);
  rpmd_ham H(nuc_beads,beta/nuc_beads,V_spring,V0);
  int num_broken = 0;

  for (int traj=0; traj<num_trajs; traj++){
    /* Load new trajectory */
    for (int bead=0; bead<nuc_beads; bead++){
      Qtraj(bead) = Q(traj*nuc_beads + bead);
      Ptraj(bead) = P(traj*nuc_beads + bead);
    }

    bool broken = false;
    double energy_0 = H.get_energy_dyn(mass,Qtraj,Ptraj);
    int sample = 0;
    while(sample < num_samples && !broken){
    // for (int sample=0; sample<num_samples; sample++){
      double energy_t = H.get_energy_dyn(mass,Qtraj,Ptraj);
      double fidelity = 100*abs((energy_t - energy_0)/energy_0);
      sample += 1;

      if (fidelity > tol){
        num_broken += 1;
        broken = true;
      }

      for (int step=0; step<stride; step++){
        vv.step(Qtraj,Ptraj);
      }
    }
  }

  std::cout << "Percent Broken: " << 100 * num_broken/num_trajs << " %" << std::endl;

  return 0;
}
int rpmd_energy_conserv::read_in_trajs(int num_trajs,int nuc_beads){

  std::string P_traj_file = root_path + "Output/Trajectories/P";
  std::string Q_traj_file = root_path + "Output/Trajectories/Q";

  P = load_file(P_traj_file, num_trajs*nuc_beads);
  Q = load_file(Q_traj_file, num_trajs*nuc_beads);

  return 0;
}
double rpmd_energy_conserv::centroid(const vector<double> &Q,int nuc_beads){
  double Qcent = 0;
  for (int bead=0; bead<nuc_beads; bead++){
    Qcent += Q(bead);
  }
  return Qcent/nuc_beads;
}
int rpmd_energy_conserv::write_ac_data(int num_trajs,int num_samples,
                                  const matrix<double> &ac_data){

  std::string myFile = root_path + "Output/auto_corr_data.txt";
  std::ofstream myStream;
  myStream.open(myFile.c_str());

  for(int traj=0; traj<num_trajs; traj++){
    for(int sample=0; sample<num_samples; sample++){
      myStream << ac_data(sample,traj) << ",";
    }
    myStream << std::endl;
  }

  myStream.close();
  return 0;
}
