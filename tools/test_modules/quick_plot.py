import numpy as np
import matplotlib.pyplot as plt


dir = "test_trajs/Model4/QQ/"

cqq = np.loadtxt(dir + "cqq")
css = np.loadtxt(dir + "css")
cbb = np.loadtxt(dir + "cbb")
csb = np.loadtxt(dir + "csb")
cbs = np.loadtxt(dir + "cbs")

fig, axs = plt.subplots(3)

axs[0].plot(cqq[:])

axs[1].plot(css[:,0])
axs[1].plot(css[:,1])
axs[1].plot(cbs[:,0])
axs[1].plot(cbs[:,1])

axs[2].plot(cbb[:,0])
axs[2].plot(cbb[:,1])
axs[2].plot(csb[:,0])
axs[2].plot(csb[:,1])

plt.show()

