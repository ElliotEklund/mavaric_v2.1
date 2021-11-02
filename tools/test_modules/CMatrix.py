import numpy as np

class C:
    def __init__(self,num_states):
        self.num_states = num_states

    def k_delta(self,n,m):

        """Standard Kronecker delta function
            n,m: both integers"""
        if n == m:
            return 1.0
        else :
            return 0.0

    def C_mat(self,x,p):

        """C Matrix from MV-RPMD formulation.
            x,p: mapping variable vectors of a given bead; both are length num_states """

        x_p_p = x + 1j*p #x plus i p
        x_m_p = x - 1j*p #x minus i p

        return np.outer(x_p_p,x_m_p) - 0.5*np.identity(self.num_states)

    def C_mat_dx(self,x,p,alpha):

        """Derivative of C Matrix w.r.t x of a particular bead and state
            x,p: mapping variables vectors of a given bead: both are length num_states
            alpha: the state of x we are taking the derivative w.r.t
            num_states: number of electronic states"""

        C_dx = np.zeros((self.num_states,self.num_states),dtype='complex')

        for i in range(self.num_states):
            for j in range(self.num_states):
                C_dx[i,j] = self.k_delta(i,alpha)*x[j] + self.k_delta(j,alpha)*x[i] + \
                    1j*(p[i]*self.k_delta(j,alpha) - self.k_delta(i,alpha)*p[j])

        return C_dx


    def C_mat_dp(self,x,p,alpha):

        """Derivative of C Matrix w.r.t p of a particular bead and state
            x,p: mapping variables vectors of a given bead: both are length num_states
            alpha: the state of p we are taking the derivative w.r.t
            num_states: number of electronic states"""

        C_dp = np.zeros((self.num_states,self.num_states),dtype='complex')

        for n in range(self.num_states):
            for m in range(self.num_states):
                C_dp[n,m] = self.k_delta(n,alpha)*p[m] + self.k_delta(m,alpha)*p[n] + \
                    1j*(x[m]*self.k_delta(n,alpha) - x[n]*self.k_delta(m,alpha))

        return C_dp
