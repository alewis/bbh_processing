import horizons_tools
import matplotlib.pyplot as plt
theta = horizons_tools.h_theta()
time = horizons_tools.Times()

plt.plot(time,theta)
plt.show()

#trimmed_time, omega_theta = horizons_tools.omega_theta()

#k_phi_r = horizons_tools.k_from_w()
#omega_r = horizons_tools.omega_r_from_w()

#k_theta_r = 0

#print omega_theta.shape
#print omega_r[:,1].shape

