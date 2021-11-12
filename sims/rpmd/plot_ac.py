import numpy as np
from numpy import genfromtxt
import matplotlib.pyplot as plt

file_name = "/Users/ellioteklund/Desktop/mavaric_v2.1/sims/rpmd/Output/auto_corr_data.txt"
tot_time = 5
num_trajs = 2000


my_data = genfromtxt(file_name, delimiter=',')[:][0:-1]
clean_data = np.delete(my_data,-1,1)
ac = np.sum(clean_data,axis=0)/num_trajs
num_samples = ac.size

tgrid = np.linspace(0,tot_time,num_samples)

plt.plot(tgrid,ac)
plt.show()
