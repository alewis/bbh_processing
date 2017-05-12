import plot_omega_theta as omth
import plot_omega_phi as omphi
import plot_omega_r as omr

def make_datfiles(
    hname="Horizons.h5",
    omname="Omega.pkl"):
  omphi.make_plot(omname)
  omth.make_plot(hname)
  omr.make_plot(omname)

if __name__ == "__main__":
  make_datfiles()
