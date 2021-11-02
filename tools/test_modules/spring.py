import numpy as np

class spring:
    def __init__(self,nuc_beads,beta,mass):
        self.nuc_beads = nuc_beads
        self.beta = beta
        self.mass = mass
        
        #Derived variables
        self.beta_nuc = self.beta/self.nuc_beads
        self.energ_pf = self.mass/(2.0*self.beta_nuc*self.beta_nuc)
        self.deriv_pf = self.mass/(self.beta_nuc*self.beta_nuc)


    def energy(self,Q):
    
        Q_shift = np.roll(Q,1)
        diff = Q - Q_shift
        return self.energ_pf*np.dot(diff,diff)

    def spring_dQ(self,Q):
    
        Q_up = np.roll(Q,1)
        Q_down = np.roll(Q,-1)
        
        return self.deriv_pf*(2.0*Q - Q_up - Q_down)

