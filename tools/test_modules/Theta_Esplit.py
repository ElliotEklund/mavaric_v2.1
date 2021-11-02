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

        """Gamma Matrix from MV-RPMD formulation using Electronic MTS
        x,p: mapping variables vectors of a given bead: both are length num_states
        Q: nuclear position"""

        gamma = np.identity(self.num_states)
        mmat = self.m.M_mat(Q[0])

        for bead in range(self.elec_beads):
            cmat = self.c.C_mat(x[bead,:],p[bead,:])
                
            gamma = np.matmul(gamma,cmat)
            gamma = np.matmul(gamma,mmat)
            
        return gamma
        
        
    def Gamma_dQ(self,Q,x,p):
        
        dMdQ = self.m.M_mat_dQ(Q)
        mmat = self.m.M_mat(Q)
        gamma_sum = np.zeros((self.num_states,self.num_states),dtype=complex)

        for i in range(self.elec_beads):
            gamma_dQ = np.identity(self.num_states)

            for j in range(self.elec_beads):
                cmat = self.c.C_mat(x[j,:],p[j,:])
                gamma_dQ = np.matmul(gamma_dQ,cmat)

                if i==j:
                    gamma_dQ = np.matmul(gamma_dQ,dMdQ)
                else:
                    gamma_dQ = np.matmul(gamma_dQ,mmat)
               
            gamma_sum += gamma_dQ

        return gamma_sum

    def Gamma_dx(self,Q,x,p,alpha,state):

        mmat = self.m.M_mat(Q)
        gamma_dx = np.identity(self.num_states)

        for i in range(self.elec_beads):
            cmat = self.c.C_mat(x[i,:],p[i,:])
    
            if alpha == i:
                dCdx = self.c.C_mat_dx(x[alpha,:],p[alpha,:],state)
                gamma_dx = np.matmul(gamma_dx,dCdx)
        
            else:
                gamma_dx = np.matmul(gamma_dx,cmat)
    
            gamma_dx = np.matmul(gamma_dx,mmat)

        return gamma_dx


    def Gamma_dp(self,Q,x,p,alpha,state):

        mmat = self.m.M_mat(Q)
        gamma_dp = np.identity(self.num_states)

        for i in range(self.elec_beads):
            cmat = self.c.C_mat(x[i,:],p[i,:])
    
            if alpha == i:
                dCdp = self.c.C_mat_dp(x[alpha,:],p[alpha,:],state)
                gamma_dp = np.matmul(gamma_dp,dCdp)
        
            else:
                gamma_dp = np.matmul(gamma_dp,cmat)
    
            gamma_dp = np.matmul(gamma_dp,mmat)

        return gamma_dp


    def Theta(self,Q,x,p):

        gamma = self.Gamma(Q,x,p)
        theta = np.trace(gamma)

        return np.real(theta)

    def grad_theta_dQ(self,Q,x,p):

        gamma = self.Gamma_dQ(Q,x,p)

        return np.real(np.trace(gamma))
        
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
