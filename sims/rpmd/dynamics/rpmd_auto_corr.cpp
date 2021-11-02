#include "rpmd_auto_corr.hpp"

rpmd_auto_corr::rpmd_auto_corr(int my_id, int num_procs, int root_proc,
                               std::string root_path)
    :my_id(my_id),
     num_procs(num_procs),
     root_proc(root_proc),
     root_path(root_path)
{}
int rpmd_auto_corr::compute_position_auto_corr(int num_trajs,int nuc_beads,double dt,
                                              double total_time,double mass, double beta,
                                              int stride){

  read_in_trajs(num_trajs,nuc_beads);
  int num_steps = total_time/dt;
  int num_samples = num_steps/stride;

  vector<double> Qtraj(nuc_beads,0);
  vector<double> Ptraj(nuc_beads,0);
  matrix<double> ac_data(num_samples,num_trajs,0);

  rpmd_vv vv(nuc_beads,mass,beta,dt);
  double Qcent_0 = 0;
  double Qcent_t = 0;

  for (int traj=0; traj<num_trajs; traj++){
    for (int bead=0; bead<nuc_beads; bead++){
      Qtraj(bead) = Q(traj*nuc_beads + bead);
      Ptraj(bead) = P(traj*nuc_beads + bead);
    }

    Qcent_0 = centroid(Qtraj,nuc_beads);

    for (int sample=0; sample<num_samples; sample++){
      Qcent_t = centroid(Qtraj,nuc_beads);
      ac_data(sample,traj) = Qcent_0*Qcent_t;
      for (int step=0; step<stride; step++){
        vv.step(Qtraj,Ptraj);
      }
    }
  }

  write_ac_data(num_trajs,num_samples,ac_data);
  return 0;
}
int rpmd_auto_corr::read_in_trajs(int num_trajs,int nuc_beads){

  std::string P_traj_file = root_path + "Output/Trajectories/P";
  std::string Q_traj_file = root_path + "Output/Trajectories/Q";

  P = load_file(P_traj_file, num_trajs*nuc_beads);
  Q = load_file(Q_traj_file, num_trajs*nuc_beads);

  return 0;
}
double rpmd_auto_corr::centroid(const vector<double> &Q,int nuc_beads){
  double Qcent = 0;
  for (int bead=0; bead<nuc_beads; bead++){
    Qcent += Q(bead);
  }
  return Qcent/nuc_beads;
}
int rpmd_auto_corr::write_ac_data(int num_trajs,int num_samples,
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
