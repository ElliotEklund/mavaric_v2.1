#include "position_auto_corr.hpp"


position_auto_corr::position_auto_corr(int my_id, int root_proc, int num_procs)
    :my_id(my_id),
     root_proc(root_proc),
     num_procs(num_procs)
{}
void position_auto_corr::compute(unsigned long long num_trajs,std::string input_dir,
    std::string output_dir,int interval){

  /* Initialize Forces and Integrator*/
  C_Matrix C(elec_beads, num_states,alpha);
  M_Matrix M(num_states, elec_beads, beta_elec_beads);
  dM_Matrix_dQ dMdQ(elec_beads, num_states, beta_elec_beads, M);
  theta_mixed theta(num_states,nuc_beads,elec_beads,C,M);
  theta_mixed_dQ theta_dQ(num_states,nuc_beads,elec_beads,C,M,dMdQ);
  theta_mixed_dElec theta_dElec(num_states,elec_beads,alpha,C,M);

  mvrpmd_mixed_forces F(nuc_beads, elec_beads, num_states, mass,
                        beta_nuc_beads, alpha, theta, theta_dQ, theta_dElec);

  ABM_MVRPMD myABM(F,dt,num_states,nuc_beads,elec_beads);


  vector<double> Q_traj (nuc_beads), P_traj (nuc_beads);
  matrix<double> x_traj (elec_beads,num_states), p_traj (elec_beads,num_states);

  vector<vector<double> > Q(num_trajs,zero_vector<double>(nuc_beads));
  vector<vector<double> > P(num_trajs,zero_vector<double>(nuc_beads));
  vector<matrix<double> > x(num_trajs,zero_matrix<double>
                                (elec_beads,num_states));
  vector<matrix<double> > p(num_trajs,zero_matrix<double>
                                (elec_beads,num_states));

  std::string Q_file = input_dir + "Q";
  std::string P_file = input_dir + "P";
  std::string x_file = input_dir + "xElec";
  std::string p_file = input_dir + "pElec";

  Q = get_trajs_reformat(Q_file,num_trajs*nuc_beads,
                         num_trajs*nuc_beads,my_id,num_procs,
                         root_proc,num_trajs,nuc_beads);

  P = get_trajs_reformat(P_file,num_trajs*nuc_beads,
                         num_trajs*nuc_beads,my_id,num_procs,
                         root_proc,num_trajs,nuc_beads);

  x = get_trajs_reformat(x_file,num_trajs*elec_beads*num_states,
                         num_trajs*elec_beads*num_states,my_id,num_procs,
                         root_proc,num_trajs,elec_beads,num_states);

  p = get_trajs_reformat(p_file,num_trajs*elec_beads*num_states,
                         num_trajs*elec_beads*num_states,my_id,num_procs,
                         root_proc,num_trajs,elec_beads,num_states);


  double sgnTheta = 0; //sign of Theta for a trajectory
  int num_steps = floor(total_time/dt);
  int num_samples = num_steps/interval;

  matrix<double> cQQ_final(num_trajs,num_samples);
  matrix<double> sign_final(num_trajs,num_samples);


  for (int traj=0; traj<num_trajs; traj++){

    /* Load new trajecty*/
    Q_traj = Q(traj);
    P_traj = P(traj);
    x_traj = x(traj);
    p_traj = p(traj);

    myABM.initialize_rk4(Q_traj, P_traj, x_traj, p_traj);
    sgnTheta = F.get_sign(Q_traj,x_traj,p_traj);

    double cent_0 = centroid(Q_traj);
    double cent_t = cent_0;
    double cQQ = cent_0 * cent_t * sgnTheta;

    for (int sample=0; sample<num_samples; sample++){
      cQQ_final(traj,sample) = cQQ;
      sign_final(traj,sample) = sgnTheta;

      for (int step=0; step<interval; step++){
        myABM.take_step(Q_traj, P_traj, x_traj, p_traj);
        cent_t = centroid(Q_traj);
        cQQ = cent_0 * cent_t * sgnTheta;
      }
    }
  }
  write_data(cQQ_final,sign_final,output_dir,num_trajs,interval,num_samples);
}
void position_auto_corr::set_system(int nuc_beadsIN, int elec_beadsIN,
      int num_statesIN, double massIN, double beta_nuc_beadsIN,
      double beta_elec_beadsIN, double alphaIN){

    nuc_beads = nuc_beadsIN;
    elec_beads = elec_beadsIN;
    num_states = num_statesIN;
    mass = massIN;
    beta_nuc_beads = beta_nuc_beadsIN;
    beta_elec_beads = beta_elec_beadsIN;
    alpha = alphaIN;
    is_sys_set = true;
}
void position_auto_corr::set_time(double dtIN, double total_timeIN){
    dt = dtIN;
    total_time = total_timeIN;
    is_time_set = true;
}
double position_auto_corr::centroid(const vector<double> & Q){
  double cent = 0;

  for (int i=0; i<nuc_beads; i++){
    cent += Q[i];
  }
  cent = cent/nuc_beads;
  return cent;
}
void position_auto_corr::write_data(const matrix<double> &cQQ_final,
    const matrix<double> &sign_final,std::string output_dir,int num_trajs,
    int interval, int num_samples){

  std::ofstream myFile;
  std::string fileName = output_dir + "position_auto_corr";
  myFile.open(fileName.c_str());

  if (!myFile.is_open()) {
      std::cout << "ERROR: Could not open " << fileName << std::endl;
  }

  myFile << "#dt:" << dt << std::endl;
  myFile << "#num_trajs:" << num_trajs << std::endl;

  for (int traj=0; traj<num_trajs; traj++){
    for(int sample=0; sample<num_samples; sample++){
      myFile << cQQ_final(traj,sample) << " ";
    }
    myFile << std::endl;
  }
  myFile.close();

  fileName = output_dir + "sign";
  myFile.open(fileName.c_str());

  //sign file
  if (!myFile.is_open()) {
      std::cout << "ERROR: Could not open " << fileName << std::endl;
  }

  myFile << "#dt:" << dt << std::endl;
  myFile << "#num_trajs:" << num_trajs << std::endl;

  for (int traj=0; traj<num_trajs; traj++){
    myFile << sign_final(traj,0) << " ";
    myFile << std::endl;
  }
  myFile.close();
}
