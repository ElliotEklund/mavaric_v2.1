import numpy as np
from Theta_Esplit import Theta
from potentials import V0

class Forces:
    def __init__(self,nuc_beads,elec_beads,num_states,beta,delta,mass):
        self.nuc_beads = nuc_beads
        self.elec_beads = elec_beads
        self.num_states = num_states
        self.beta = beta
        self.delta = delta
        self.mass = mass
        
        self.v = V0(self.nuc_beads,self.num_states,self.delta)
        self.theta = Theta(self.nuc_beads,self.elec_beads,self.num_states,self.beta,self.delta)

    def dHdP(self,P):

        return P/self.mass
    
    def dHdQ(self,Q,x,p):

        v_force = self.v.dV0_dQ(Q)
        th = self.theta.Theta(Q,x,p)
        theta_force = self.theta.grad_theta_dQ(Q,x,p)

        return v_force - theta_force/(self.beta*th)

    def dHdx(self,Q,x,p):

        th = self.theta.Theta(Q,x,p)
        theta_force = self.theta.grad_theta_dx(Q,x,p)

        return 2.0*x/self.beta - theta_force/(self.beta*th)

    def dHdp(self,Q,x,p):

        th = self.theta.Theta(Q,x,p)
        theta_force = self.theta.grad_theta_dp(Q,x,p)
    
        return 2.0*p/self.beta - theta_force/(self.beta*th)
