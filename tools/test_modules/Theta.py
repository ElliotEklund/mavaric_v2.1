import numpy as np
from MMatrix import M
from CMatrix import C

class Theta:
    def __init__(self,nuc_beads,elec_beads,num_states,beta,delta):
        self.nuc_beads = nuc_beads
        self.elec_beads = elec_beads
        self.num_states = num_states
        self.beta = beta
        self.delta = delta
        
        self.m = M(self.nuc_beads,self.elec_beads,self.num_states,self.beta,self.delta)
        self.c = C(self.num_states)
        

    def Gamma(self,Q,x,p):

        """Gamma Matrix from MV-RPMD formulation.
            x,p: mapping variables vectors of a given bead: both are length num_states
            Q: nuclear position"""

        gamma = np.identity(self.num_states)
        ratio = int(self.nuc_beads/self.elec_beads)

        for bead in range(self.elec_beads):
            M = self.m.M_mat(Q[bead*ratio])
            C = self.c.C_mat(x[bead,:],p[bead,:])
    
            gamma = np.matmul(gamma,C)
            gamma = np.matmul(gamma,M)

        return gamma


    def Gamma_dQ(self,Q,x,p,alpha):

        """Derivative of Gamma Matrix w.r.t Beta from MV-RPMD formulation using Electronic MTS
            x,p: mapping variables vectors of a given bead: both are length num_states
            Q: nuclear position"""

        ratio = int(self.nuc_beads/self.elec_beads)

        W = np.zeros((self.elec_beads,self.nuc_beads))
        gamma_dQ = np.identity(self.num_states)

        for i in range (self.elec_beads):
            W[i,i*ratio] = 1.0

        Q_trans = np.matmul(W,Q)

        for bead in range(self.elec_beads):
            C = self.c.C_mat(x[bead,:],p[bead,:])
            gamma_dQ = np.matmul(gamma_dQ,C)

            if alpha == bead:
                M_dQ = self.m.M_mat_dQ(Q_trans[bead])
                gamma_dQ = np.matmul(gamma_dQ,M_dQ)

            else:
                M = self.m.M_mat(Q_trans[bead])
                gamma_dQ = np.matmul(gamma_dQ,M)

        return gamma_dQ

    def Gamma_dx(self,Q,x,p,alpha,state):
        """Derivative wrt bead alpha and state"""

        gamma_dx = np.identity(self.num_states)
        ratio = int(self.nuc_beads/self.elec_beads)

        for bead in range(self.elec_beads):

            if alpha == bead:
                C_dx = self.c.C_mat_dx(x[bead,:],p[bead,:],state)
                gamma_dx = np.matmul(gamma_dx,C_dx)

            else:
                C = self.c.C_mat(x[bead,:],p[bead,:])
                gamma_dx = np.matmul(gamma_dx,C)

            M = self.m.M_mat(Q[ratio*bead])
            gamma_dx = np.matmul(gamma_dx,M)

        return gamma_dx


    def Gamma_dp(self,Q,x,p,alpha,state):
        """Derivative wrt bead alpha and state"""

        gamma_dp = np.identity(self.num_states)
        ratio = int(self.nuc_beads/self.elec_beads)

        for bead in range(self.elec_beads):

            if alpha == bead:
                C_dp = self.c.C_mat_dp(x[bead,:],p[bead,:],state)
                gamma_dp = np.matmul(gamma_dp,C_dp)

            else:
                C = self.c.C_mat(x[bead,:],p[bead,:])
                gamma_dp = np.matmul(gamma_dp,C)

            M = self.m.M_mat(Q[ratio*bead])
            gamma_dp = np.matmul(gamma_dp,M)

        return gamma_dp


    def Theta(self,Q,x,p):

        gamma = self.Gamma(Q,x,p)
        theta = np.trace(gamma)

        return np.real(theta)

    def grad_theta_dQ(self,Q,x,p):

        ratio = int(self.nuc_beads/self.elec_beads)
        v = np.zeros((self.nuc_beads,self.elec_beads))
        grad = np.zeros(self.elec_beads)

        for i in range (self.elec_beads):
            v[i*ratio,i] = 1.0

        for bead in range(self.elec_beads):
            grad[bead] = np.real(np.trace(self.Gamma_dQ(Q,x,p,bead)))

        return np.matmul(v,grad)

    def grad_theta_dx(self,Q,x,p):

        grad = np.zeros((self.elec_beads,self.num_states))

        for bead in range(self.elec_beads):
            for state in range (self.num_states):
                grad[bead,state] = np.real(np.trace(self.Gamma_dx(Q,x,p,bead,state)))
        return grad

    def grad_theta_dp(self,Q,x,p):

        grad = np.zeros((self.elec_beads,self.num_states))

        for bead in range(self.elec_beads):
            for state in range (self.num_states):
                grad[bead,state] = np.real(np.trace(self.Gamma_dp(Q,x,p,bead,state)))
        return grad
