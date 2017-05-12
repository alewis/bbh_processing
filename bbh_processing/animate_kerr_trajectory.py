"""Produce an animation of a kerr trajectory.
"""
import matplotlib.pyplot as plt
import analyze_kerr_data as kerr
import hfile_tools as hfile
from matplotlib import animation
import mpl_toolkits.mplot3d.axes3d as p3
from matplotlib.pyplot import cm
from mpl_toolkits.mplot3d.art3d import Line3DCollection
from matplotlib.collections import LineCollection
import numpy as np
def lims(array):
    return [np.amin(array), np.amax(array)]

def update(i, rspace, time, rmag, theta, lines, factor):
    p = factor*i
    #c = next(colour)
    lines[0].set_data(rspace[:p, 0], rspace[:p, 1])
    lines[0].set_3d_properties(rspace[:p, 2])
    

    lines[1].set_data(rmag[:p], theta[:p])
    lines[2].set_data(time[:p], rmag[:p])
    lines[3].set_data(time[:p], theta[:p])
    # for l in lines:
      # l.set_color(c)
    return lines


def extract_kerr_data(fname):
    rvec = kerr.loadfile(fname=fname)[:25000, :]
    time = rvec[:, 0]
    rspace = rvec[:, 1:]
    rmag = kerr.rmag(rvec[:,1:])
    theta = kerr.theta(rvec[:,1:], [0, 0, 1])
    return time, rspace, rmag, theta

def extract_bbh_data(fname):
    rspace = hfile.get_coord("rvec", fname)
    time = hfile.get_coord("t", fname)
    rmag = hfile.get_coord("r", fname)
    theta = hfile.get_coord("theta", fname)
    return time, rspace, rmag, theta

def make_animation(fname, ftype="bbh", factor=100):
    if ftype=="bbh":
        time, rspace, rmag, theta = extract_bbh_data(fname)
    elif ftype=="kerr":
        time, rspace, rmag, theta = extract_kerr_data(fname)
    else:
        raise ValueError("ftype must be 'bbh' or 'kerr")
    # x = rvec[:, 1]
    # y = rvec[:, 2]
    # z = rvec[:, 3]
    fig = plt.figure()

    ax = fig.add_subplot(2, 2, 1, projection='3d') 
    ax.set_xlim3d(lims(rspace[:,0]))
    ax.set_xlabel('X')
    ax.set_ylim3d(lims(rspace[:,1]))
    ax.set_ylabel('Y')
    ax.set_zlim3d(lims(rspace[:,2]))
    ax.set_zlabel('Z')
    ax.set_title('Trajectory')

    ax2 = fig.add_subplot(2, 2, 2) 
    ax2.set_xlabel("r")
    ax2.set_ylabel("theta")

    ax3 = fig.add_subplot(2, 2, 3)
    ax3.set_xlabel("time")
    ax3.set_ylabel("r")
    
    ax4 = fig.add_subplot(2, 2, 4)
    ax4.set_xlabel("time")
    ax4.set_ylabel("theta")

    # lines = []
    # mycmap=cm.Paired(np.arange(0, 1, float(factor)/len(time)))
    # rpoints = 
    # lines.append(Line3DCollection([], cmap=mycmap))
    # lines[0].set_array(rspace[0:1, 0], rspace[0:1, 1], rspace[0:1, 2])
    # ax.add_collection(lines[0])
    lines = [ax.plot(rspace[0:1, 0], rspace[0:1, 1], rspace[0:1, 2])[0],
             ax2.plot(rmag, theta)[0],
             ax3.plot(time, rmag)[0],
             ax4.plot(time, theta)[0]]

    line_ani = animation.FuncAnimation(fig, update, 
                                       fargs=(rspace, time, 
                                              rmag, theta, lines, factor),
                                       interval=50,
                                       save_count=len(time),
                                       blit=True)
    line_ani.save('kerr_orbit_animation.mp4', fps=30, extra_args=['-vcodec', 'libx264'])
    #plt.show()

if __name__=="__main__":
    make_animation("/home/adlewis/m/ResonanceAnalysis/Notes/GenericBBH/ExactOrbits/long_orbits/Kerr_p20_ecc9_zm6_a9.dat",
                  ftype="kerr", factor=100)
    #bbhpath = "/mnt/raid-project/nr/adlewis/pythontiming/Data/AaronsEccSeries2/Ecc_pt2_inc_40/Horizons.h5"
    #make_animation(bbhpath, ftype="bbh", factor=10)
