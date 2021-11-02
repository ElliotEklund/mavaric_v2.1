import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import time
from ABM import ABM
from Theta import Theta

def semi(x,p):

    xsq_psq = np.square(x) + np.square(p) - 1.0
    pops = np.sum(xsq_psq,axis=0)/(2.0*elec_beads)

    return pops

def boltz(Q,x,p):
    gam = theta.Gamma(Q,x,p)
    denom = np.trace(np.real(gam))
    s1 = np.real(gam[0,0])
    s2 = np.real(gam[1,1])
    pops = np.array([s1,s2])/denom

    return pops
    
    
#main ##########################################################################
dir = 'test_trajs/Model4/QQ/'

#Physical Parameters
nuc_beads = int(4)
elec_beads = int(4)
delta = 0.1

num_states = int(2)
temp = 1.0
beta = 1.0/temp
dt = 0.001
mass = 1.0

#Simulation Parameters
total_time = 5
num_cycles = int(6)
stride = int(10)
num_steps = int(total_time/dt)
num_samples = int(num_steps/stride)

#Read in and initialize trajectories
Qin = pd.read_csv(dir + 'Q',dtype=np.double).values
Pin = pd.read_csv(dir + 'P',dtype=np.double).values
xin = pd.read_csv(dir + 'xElec',dtype=np.double).values
pin = pd.read_csv(dir + 'pElec',dtype=np.double).values

num_trajs = int(Qin.size/nuc_beads)

Qrf = np.resize(Qin,(num_trajs,nuc_beads))
Prf = np.resize(Pin,(num_trajs,nuc_beads))
xrf = np.resize(xin,(num_trajs,elec_beads,num_states))
prf = np.resize(pin,(num_trajs,elec_beads,num_states))

#Define vectors for collecting data
cqq = np.zeros(num_samples)
css = np.zeros((num_samples,num_states))
cbb = np.zeros((num_samples,num_states))
csb = np.zeros((num_samples,num_states))
cbs = np.zeros((num_samples,num_states))

Qtraj = np.zeros((num_samples,nuc_beads))
xtraj = np.zeros((num_cycles,num_samples,elec_beads,num_states))
ptraj = np.zeros((num_cycles,num_samples,elec_beads,num_states))

timess = np.zeros(num_samples)

count = int(0)
theta_sum = 0

abm = ABM(nuc_beads,elec_beads,num_states,beta,delta,mass,dt)
theta = Theta(nuc_beads,elec_beads,num_states,beta,delta)

t0 = time.time()

for traj in range(num_cycles):

    Q = Qrf[traj,:]
    P = Prf[traj,:]
    x = xrf[traj,:]
    p = prf[traj,:]

    Q0 = np.sum(Q)/nuc_beads
    semi_0 = semi(x,p)
    boltz_0 = boltz(Q,x,p)
    sgn_theta = np.sign(theta.Theta(Q,x,p))
 
    theta_sum += sgn_theta

    abm.initialize_rk4(Q,P,x,p)

    count = int(0)

    for step in range(num_steps):
        abm.step(Q,P,x,p)

        if step % stride == 0:
            Qtraj[count,:] = Q
            xtraj[traj,count,:,:] = x
            ptraj[traj,count,:,:] = p
         
            Qcent = np.sum(Q)/nuc_beads
            semi_t = semi(x,p)
            boltz_t = boltz(Q,x,p)

            cqq[count]   += Q0 * Qcent * sgn_theta
            css[count,:] += semi_0 * semi_t * sgn_theta
            cbb[count,:] += boltz_0 * boltz_t * sgn_theta
            csb[count,:] += semi_0 * boltz_t * sgn_theta
            cbs[count,:] += boltz_0 * semi_t * sgn_theta

            timess[count] = step * dt
            count += 1

t1 = time.time()
print("Time:",t1-t0)

#fig, axs = plt.subplots(1)

cqq = cqq/theta_sum
css = css/theta_sum
cbb = cbb/theta_sum
csb = csb/theta_sum
cbs = cbs/theta_sum

np.savetxt(dir + 'cqq',cqq)
np.savetxt(dir + 'css',css)
np.savetxt(dir + 'cbb',cbb)
np.savetxt(dir + 'csb',csb)
np.savetxt(dir + 'cbs',cbs)
