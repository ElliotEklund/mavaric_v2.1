import numpy as np
from potentials import V0

class M:
    def __init__(self,nuc_beads,elec_beads,num_states,beta,delta):
        self.elec_beads = elec_beads
        self.num_states = num_states
        self.beta = beta
        self.delta = delta
        self.beta_Ne = self.beta/self.elec_beads
        
        self.v = V0(nuc_beads,num_states,delta)
        
    def M_mat(self,Q):

        """M Matrix from MV-RPMD formulation.
            Q: nuclear position"""

        M = np.zeros((self.num_states,self.num_states))
        vmat = self.v.Vmat(Q)

        M[0,0] = np.exp(- self.beta_Ne * vmat[0,0])
        M[1,1] = np.exp(- self.beta_Ne * vmat[1,1])
        M[0,1] = -self.beta_Ne * vmat[0,1] * np.exp(-self.beta_Ne * vmat[0,0])
        M[1,0] = -self.beta_Ne * vmat[1,0] * np.exp(-self.beta_Ne * vmat[1,1])

        return M

    def M_mat_dQ(self,Q):

        """Derivative of M Matrix w.r.t Q.
            Q: nuclear position"""

        M_dQ = np.zeros((self.num_states,self.num_states))
        vmat = self.v.Vmat(Q)
        v_dQ = self.v.Vmat_dQ(Q)

        M_dQ[0,0] = - self.beta_Ne * v_dQ[0,0] * np.exp(-self.beta_Ne * vmat[0,0])
        M_dQ[1,1] = - self.beta_Ne * v_dQ[1,1] * np.exp(-self.beta_Ne * vmat[1,1])

        M_dQ[0,1] = - self.beta_Ne * v_dQ[0,1] * np.exp(-self.beta_Ne * vmat[0,0]) -\
                 (self.beta_Ne * vmat[0,1] ) * M_dQ[0,0]

        M_dQ[1,0] = - self.beta_Ne * v_dQ[1,0] * np.exp(-self.beta_Ne * vmat[1,1])  -\
                 (self.beta_Ne * vmat[1,0] ) * M_dQ[1,1]

        return M_dQ
