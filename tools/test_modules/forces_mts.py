import numpy as np
from Theta import Theta
from potentials import V0
from spring import spring

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
        self.Spring = spring(self.nuc_beads,self.beta,self.mass)
        
        #Derived constants
        self.beta_nuc = self.beta/self.nuc_beads
        self.one_beta_nuc = 1.0/(self.beta_nuc) #1.0/(beta/nuc_beads)
        self.one_nuc = 1.0/self.nuc_beads

    def dHdP(self,P):

#        print("dHdP",P/self.mass)
        return P/self.mass
    
    def dHdQ(self,Q,x,p):

        v_force = self.v.dV0_dQ(Q)
        th = self.theta.Theta(Q,x,p)
        theta_force = self.theta.grad_theta_dQ(Q,x,p)

#
#        print("dHdQ",self.Spring.spring_dQ(Q) + v_force - self.one_beta_nuc*theta_force/th)
#        print("spring:",self.Spring.spring_dQ(Q))
#        print("pot:",v_force)
#        print("theta:",th)
#        print("theta_force:",theta_force)
        
        return self.Spring.spring_dQ(Q) + v_force - self.one_beta_nuc*theta_force/th

    def dHdx(self,Q,x,p):

        th = self.theta.Theta(Q,x,p)
        theta_force = self.theta.grad_theta_dx(Q,x,p)

#        print("dHdx",(2.0*x - theta_force/th)*self.one_beta_nuc)
#        print("first term:",2.0*x*self.one_beta_nuc)
#        print("theta_dx:",theta_force)
        
        
        return (2.0*x - theta_force/th)*self.one_beta_nuc

    def dHdp(self,Q,x,p):

        th = self.theta.Theta(Q,x,p)
        theta_force = self.theta.grad_theta_dp(Q,x,p)
    
#        print("dHdp",(2.0*p - theta_force/th)*self.one_beta_nuc)
#        print("first term:",2.0*p*self.one_beta_nuc)
#        print("theta_dp:",theta_force)
        
        return (2.0*p - theta_force/th)*self.one_beta_nuc
