import numpy as np

class V0:
    def __init__(self,nuc_beads,num_states,delta):
        self.nuc_beads = nuc_beads
        self.num_states = num_states
        self.delta = delta

    def V(self,Q):

        """State Independent Potential
        Q: nuclear position vector"""

        v = 0.5*np.dot(Q,Q)
        return v

    def dV0_dQ(self,Q):
    
        """Derivative of State Independent Potential w.r.t Q
        Q: nuclear position vector"""
    
        return Q

    def Vmat(self,Q):

        """Diabatic potential matrix for a 2 state system with
        constant off-diagonal coupling.
        Q: nuclear position"""

        v = np.zeros((self.num_states,self.num_states))
        v[0,0] = Q + 3.0
        v[1,1] = -Q

        v[0,1] = self.delta
        v[1,0] = self.delta

        return v

    def Vmat_dQ(self,Q):

        """Derivative of Diabatic potential matrix w.r.t Q for a 2 state
        system with constant off-diagonal coupling.
        Q: nuclear position"""

        v = np.zeros((self.num_states,self.num_states))
        v[0,0] = 1.0
        v[1,1] = -1.0

        return v
