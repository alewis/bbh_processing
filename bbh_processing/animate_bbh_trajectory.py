"""Produce an animation of a kerr trajectory.
"""
import matplotlib.pyplot as plt
import bbh_processing.hfile_tools as hfile

from matplotlib import animation
import mpl_toolkits.mplot3d.axes3d as p3
import numpy as np
def lims(array):
    return [np.amin(array), np.amax(array)]

def update(i, rspace, time, rmag, theta, lines, factor):
    p = factor*i
    lines[0].set_data(rspace[:p, 0], rspace[:p, 1])
    lines[0].set_3d_properties(rspace[:p, 2])

    lines[1].set_data(rmag[:p], theta[:p])
    lines[2].set_data(time[:p], rmag[:p])
    lines[3].set_data(time[:p], theta[:p])
    return lines




def make_animation(fname):
    rvec = kerr.loadfile(fname=fname)[:25000, :]
    time = rvec[:, 0]
    rspace = rvec[:, 1:]
    # x = rvec[:, 1]
    # y = rvec[:, 2]
    # z = rvec[:, 3]
    fig = plt.figure()

    ax = fig.add_subplot(2, 2, 1, projection='3d') 
    ax.set_xlim3d(lims(rvec[:,1]))
    ax.set_xlabel('X')
    ax.set_ylim3d(lims(rvec[:,2]))
    ax.set_ylabel('Y')
    ax.set_zlim3d(lims(rvec[:,3]))
    ax.set_zlabel('Z')
    ax.set_title('Trajectory')

    rmag = kerr.rmag(rvec[:,1:])
    theta = kerr.theta(rvec[:,1:], [0, 0, 1])
    ax2 = fig.add_subplot(2, 2, 2) 
    ax2.set_xlabel("r")
    ax2.set_ylabel("theta")

    ax3 = fig.add_subplot(2, 2, 3)
    ax3.set_xlabel("time")
    ax3.set_ylabel("r")
    
    ax4 = fig.add_subplot(2, 2, 4)
    ax4.set_xlabel("time")
    ax4.set_ylabel("theta")
    
    lines = [ax.plot(rvec[0:1, 1], rvec[0:1, 2], rvec[0:1, 3])[0],
             ax2.plot(rmag, theta)[0],
             ax3.plot(time, rmag)[0],
             ax4.plot(time, theta)[0]]
    line_ani = animation.FuncAnimation(fig, update, 10000,
                                       fargs=(rspace, time, 
                                              rmag, theta, lines),
                                       interval=50,
                                       blit=False)
    line_ani.save('kerr_orbit_animation.mp4', fps=30, extra_args=['-vcodec', 'libx264'])
    #plt.show()

if __name__=="__main__":
    make_animation("/mnt/raid-project/nr/adlewis/pythontiming/Data/AaronsEccSeries2/Ecc_pt2_inc_40")
