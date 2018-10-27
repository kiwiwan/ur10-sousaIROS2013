# Initializations
datafolder = 'data/'
tmpfolder = 'tmp/'

from sympy import init_printing
init_printing()

import pickle
import sympybotics
import numpy


with open(tmpfolder +  'robotmodels/ur10_model.pkl', 'rb' ) as file:
          rbt = pickle.load( file )

with open(datafolder +  'trajectories/traj_shwfl_abq0.pkl', 'rb' ) as file:
          s_h_wf_l, a_b_q0 = pickle.load(file)   
          

h = 0.001
decimate = 10
h_plot = h*decimate

with open(datafolder +  'trajectories/traj.dat', 'r' ) as file:
          q_ref_orig = numpy.loadtxt(file)
s = q_ref_orig.shape[0] / decimate

q_ref = numpy.zeros((s, rbt.dof))
for i in range(s):
    q_ref[i, :] = q_ref_orig[i*decimate, :]


si = 1/h_plot
sf = si + 20/h_plot
t = numpy.arange(sf-si) * h_plot

from matplotlib import pyplot as plt
plt.close()

for i in range(rbt.dof):
    plt.plot(t,q_ref[si:sf,i], label="$q_%d$"%(i+1))
plt.legend()

plt.xlabel("Time (s)")
plt.ylabel("Joint positions (rad)")

plt.show()



