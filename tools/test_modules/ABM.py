import numpy as np
from RK4 import RK4
#from forces import Forces
from forces_mts import Forces


class ABM:
    def __init__(self,nuc_beads,elec_beads,num_states,beta,delta,mass,dt):
     
        self.nuc_beads = nuc_beads
        self.elec_beads = elec_beads
        self.num_states = num_states
        self.beta = beta
        self.delta = delta
        self.mass = mass
        self.dt = dt
    
        self.rk4 = RK4(nuc_beads,elec_beads,num_states,beta,delta,mass,-dt)

        self.F = Forces(self.nuc_beads,self.elec_beads,self.num_states,self.beta,
                        self.delta,self.mass)

        self.h1_p = self.dt * 55.0/24.0
        self.h2_p = self.dt*59.0/24.0
        self.h3_p = self.dt*37.0/24.0
        self.h4_p = self.dt*9.0/24.0
        self.h1_c = self.dt*9.0/24.0
        self.h2_c = self.dt*19.0/24.0
        self.h3_c = self.dt*5.0/24.0
        self.h4_c = self.dt*1.0/24.0

        self.f_Q_3 = np.zeros(self.nuc_beads)
        self.f_Q_2 = np.zeros(self.nuc_beads)
        self.f_Q_1 = np.zeros(self.nuc_beads)
        self.f_Q_0 = np.zeros(self.nuc_beads)
        self.f_Q_p1 = np.zeros(self.nuc_beads)
        
        self.f_P_3 = np.zeros(self.nuc_beads)
        self.f_P_2 = np.zeros(self.nuc_beads)
        self.f_P_1 = np.zeros(self.nuc_beads)
        self.f_P_0 = np.zeros(self.nuc_beads)
        self.f_P_p1 = np.zeros(self.nuc_beads)

        self.f_x_3 = np.zeros((self.elec_beads,self.num_states))
        self.f_x_2 = np.zeros((self.elec_beads,self.num_states))
        self.f_x_1 = np.zeros((self.elec_beads,self.num_states))
        self.f_x_0 = np.zeros((self.elec_beads,self.num_states))
        self.f_x_p1 = np.zeros((self.elec_beads,self.num_states))

        self.f_p_3 = np.zeros((self.elec_beads,self.num_states))
        self.f_p_2 = np.zeros((self.elec_beads,self.num_states))
        self.f_p_1 = np.zeros((self.elec_beads,self.num_states))
        self.f_p_0 = np.zeros((self.elec_beads,self.num_states))
        self.f_p_p1 = np.zeros((self.elec_beads,self.num_states))

        self.Q_pred = np.zeros(self.nuc_beads)
        self.P_pred = np.zeros(self.nuc_beads)
        self.x_pred = np.zeros((self.elec_beads,self.num_states))
        self.p_pred = np.zeros((self.elec_beads,self.num_states))
        
        
    def initialize_rk4(self,Q,P,x,p):
        Q_copy = np.zeros(self.nuc_beads)
        P_copy = np.zeros(self.nuc_beads)
        x_copy = np.zeros((self.elec_beads,self.num_states))
        p_copy = np.zeros((self.elec_beads,self.num_states))

        Q_copy[:] = Q
        P_copy[:] = P
        x_copy[:] = x
        p_copy[:] = p

        self.rk4.step(Q,P,x,p)
        self.update_f_1(Q,P,x,p)
        
        self.rk4.step(Q,P,x,p)
        self.update_f_2(Q,P,x,p)
        
        self.rk4.step(Q,P,x,p)
        self.update_f_3(Q,P,x,p)
        
        Q[:] = Q_copy
        P[:] = P_copy
        x[:] = x_copy
        p[:] = p_copy
            
    def step(self,Q,P,x,p):
        self.predict(Q,P,x,p)
        self.correct(Q,P,x,p)
        
    def predict(self,Q,P,x,p):
    
        self.update_f_0(Q,P,x,p)
        
        self.pred_Q = Q + (self.h1_p * self.f_Q_0) - (self.h2_p * self.f_Q_1) + (self.h3_p * self.f_Q_2) - (self.h4_p * self.f_Q_3);
        self.pred_P = P + (self.h1_p * self.f_P_0) - (self.h2_p * self.f_P_1) + (self.h3_p * self.f_P_2) - (self.h4_p * self.f_P_3);
        self.pred_x = x + (self.h1_p * self.f_x_0) - (self.h2_p * self.f_x_1) + (self.h3_p * self.f_x_2) - (self.h4_p * self.f_x_3);
        self.pred_p = p + (self.h1_p * self.f_p_0) - (self.h2_p * self.f_p_1) + (self.h3_p * self.f_p_2) - (self.h4_p * self.f_p_3);

    def correct(self,Q,P,x,p):
    
          self.update_f_p1(self.pred_Q,self.pred_P,self.pred_x,self.pred_p);
          
          Q  += (self.h1_c * self.f_Q_p1)  + (self.h2_c * self.f_Q_0)  - (self.h3_c * self.f_Q_1)  + (self.h4_c * self.f_Q_2);
          P  += (self.h1_c * self.f_P_p1)  + (self.h2_c * self.f_P_0)  - (self.h3_c * self.f_P_1)  + (self.h4_c * self.f_P_2);
          x  += (self.h1_c * self.f_x_p1)  + (self.h2_c * self.f_x_0)  - (self.h3_c * self.f_x_1)  + (self.h4_c * self.f_x_2);
          p  += (self.h1_c * self.f_p_p1)  + (self.h2_c * self.f_p_0)  - (self.h3_c * self.f_p_1)  + (self.h4_c * self.f_p_2);
          
          self.f_Q_3 = self.f_Q_2;
          self.f_Q_2 = self.f_Q_1;
          self.f_Q_1 = self.f_Q_0;

          self.f_P_3 = self.f_P_2;
          self.f_P_2 = self.f_P_1;
          self.f_P_1 = self.f_P_0;
        
          self.f_x_3 = self.f_x_2;
          self.f_x_2 = self.f_x_1;
          self.f_x_1 = self.f_x_0;

          self.f_p_3 = self.f_p_2;
          self.f_p_2 = self.f_p_1;
          self.f_p_1 = self.f_p_0;

    def update_f_3(self,Q,P,x,p):
        self.f_Q_3 =   self.F.dHdP(P)
        self.f_P_3 = - self.F.dHdQ(Q,x,p)
        self.f_x_3 =   self.F.dHdp(Q,x,p)
        self.f_p_3 = - self.F.dHdx(Q,x,p)

    def update_f_2(self,Q,P,x,p):
        self.f_Q_2 =   self.F.dHdP(P)
        self.f_P_2 = - self.F.dHdQ(Q,x,p)
        self.f_x_2 =   self.F.dHdp(Q,x,p)
        self.f_p_2 = - self.F.dHdx(Q,x,p)
        
    def update_f_1(self,Q,P,x,p):
        self.f_Q_1 =   self.F.dHdP(P)
        self.f_P_1 = - self.F.dHdQ(Q,x,p)
        self.f_x_1 =   self.F.dHdp(Q,x,p)
        self.f_p_1 = - self.F.dHdx(Q,x,p)

    def update_f_0(self,Q,P,x,p):
        self.f_Q_0 =   self.F.dHdP(P)
        self.f_P_0 = - self.F.dHdQ(Q,x,p)
        self.f_x_0 =   self.F.dHdp(Q,x,p)
        self.f_p_0 = - self.F.dHdx(Q,x,p)
        
        
    def update_f_p1(self,Q,P,x,p):
        self.f_Q_p1 =   self.F.dHdP(P)
        self.f_P_p1 = - self.F.dHdQ(Q,x,p)
        self.f_x_p1 =   self.F.dHdp(Q,x,p)
        self.f_p_p1 = - self.F.dHdx(Q,x,p)
