import numpy as np
import potentials as V
import Theta as theta

def energy(Q,P,x,p,nuc_beads,elec_beads,num_states,beta,delta):

    T = np.dot(P,P)/2.0
    v = V.V0(Q,nuc_beads)
    G = np.dot(x[0,:],x[0,:]) + np.dot(p[0,:],p[0,:])
    th = theta.theta(Q,x,p,nuc_beads,elec_beads,num_states,beta,delta)

    energ = T + v + G - np.log(np.absolute(th))
    return energ
