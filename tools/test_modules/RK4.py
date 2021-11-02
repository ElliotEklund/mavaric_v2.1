import numpy as np
#from forces import Forces
from forces_mts import Forces

class RK4:

    def __init__(self,nuc_beads,elec_beads,num_states,beta,delta,mass,dt):

        self.nuc_beads = nuc_beads
        self.elec_beads = elec_beads
        self.num_states = num_states
        self.beta = beta
        self.delta = delta
        self.mass = mass
        self.dt = dt
        
        self.F = Forces(self.nuc_beads,self.elec_beads,self.num_states,self.beta,
                        self.delta,self.mass)
    
        self.dt_half = 0.5*self.dt
        self.c1 = self.dt/6.0
        self.c2=2.0*self.dt/6.0
        self.c3=self.c2
        self.c4=self.c1
    
        self.k1Q = np.zeros(self.nuc_beads)
        self.k2Q = np.zeros(self.nuc_beads)
        self.k3Q = np.zeros(self.nuc_beads)
        self.k4Q = np.zeros(self.nuc_beads)
    
        self.k1P = np.zeros(self.nuc_beads)
        self.k2P = np.zeros(self.nuc_beads)
        self.k3P = np.zeros(self.nuc_beads)
        self.k4P = np.zeros(self.nuc_beads)

        self.k1x = np.zeros((self.elec_beads,self.num_states))
        self.k2x = np.zeros((self.elec_beads,self.num_states))
        self.k3x = np.zeros((self.elec_beads,self.num_states))
        self.k4x = np.zeros((self.elec_beads,self.num_states))
    
        self.k1p = np.zeros((self.elec_beads,self.num_states))
        self.k2p = np.zeros((self.elec_beads,self.num_states))
        self.k3p = np.zeros((self.elec_beads,self.num_states))
        self.k4p = np.zeros((self.elec_beads,self.num_states))

        
    def update_k1(self,Q,P,x,p):
        self.k1Q = self.F.dHdP(P)
        self.k1P = -self.F.dHdQ(Q,x,p)
        self.k1x = self.F.dHdp(Q,x,p)
        self.k1p = -self.F.dHdx(Q,x,p)
        
    def update_k2(self,Q,P,x,p):
        self.k2Q = self.F.dHdP(P)
        self.k2P = -self.F.dHdQ(Q,x,p)
        self.k2x = self.F.dHdp(Q,x,p)
        self.k2p = -self.F.dHdx(Q,x,p)
        
    def update_k3(self,Q,P,x,p):
        self.k3Q = self.F.dHdP(P)
        self.k3P = -self.F.dHdQ(Q,x,p)
        self.k3x = self.F.dHdp(Q,x,p)
        self.k3p = -self.F.dHdx(Q,x,p)
        
    def update_k4(self,Q,P,x,p):
        self.k4Q = self.F.dHdP(P)
        self.k4P = -self.F.dHdQ(Q,x,p)
        self.k4x = self.F.dHdp(Q,x,p)
        self.k4p = -self.F.dHdx(Q,x,p)

    def update_final(self,Q,P,x,p):
        
        Q += (self.c1*self.k1Q) + (self.c2*self.k2Q) + (self.c3*self.k3Q) + (self.c4*self.k4Q)
        P += (self.c1*self.k1P) + (self.c2*self.k2P) + (self.c3*self.k3P) + (self.c4*self.k4P)
        x += (self.c1*self.k1x) + (self.c2*self.k2x) + (self.c3*self.k3x) + (self.c4*self.k4x)
        p += (self.c1*self.k1p) + (self.c2*self.k2p) + (self.c3*self.k3p) + (self.c4*self.k4p)

    def step(self,Q,P,x,p):

        self.update_k1(Q,P,x,p)
        
        self.update_k2(Q+self.dt_half*self.k1Q, P+self.dt_half*self.k1P,
                  x+self.dt_half*self.k1x, p+self.dt_half*self.k1p)
                  
        self.update_k3(Q+self.dt_half*self.k2Q, P+self.dt_half*self.k2P,
                  x+self.dt_half*self.k2x, p+self.dt_half*self.k2p)
              
        self.update_final(Q,P,x,p)
