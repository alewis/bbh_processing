#This class represents a Cartesian basis which may rotate with respect to some
#other basis over time.

#INITIALIZATION
#Given a time series of 3D vectors V=(t,Vx,Vy,Vz) (stored as a numpy array
#such that V[:,0] is t, etc), the constructor should generate a series of axes
#
#The rotation is about the cross product between V and hat(z).

#METHODS
#t() -> returns a view into Vector[:,0]
#V() -> returns a view into Vector[:,1:]
#transform(R,T=None) -> input is a 3-column numpy array of vectors expressed in the
#same basis as Axis, plus, optionally, a 1-column numpy array of time points. If
#T is None we assume R to be on the same time series as 

#MEMBERS
#t    -> a 1-column numpy array storing the times
#vector -> a 3-column numpy array storing [Vx,Vy,Vz] at each t.
#mag    -> a 1-column numpy array storing |V| at each t.
#length -> length(t); just for convenience
#M -> the transformation matrices performing the rotation from the basis in
#          which Vector is expressed to that in which Vector(t) = [t,0,0,|V_z|]

import numpy as np
import numpy.linalg as la
from numpy.lib import stride_tricks
import utilities as bbh
class CartesianBasisSeries:
    #CONSTRUCTION METHODS
    #Two constructors to avoid confusion about how time is dealt with. 
    def __init__(self,V):
      if V.shape[1] != 4:
        raise ValueError("this constructor expects a time-series")
      self._set_values(t,V[:,1:]) 

    def __init__(self,t,V):
      if V.shape[1] != 3:
        raise ValueError("this constructor expects a time and a series of vectors")
      self._set_values(t,V)
      
    #Store data for construction.
    def _set_values(self,t,V):
      self.t=np.copy(t)
      self.length = len(t)
      zaxis = stride_tricks.as_strided(np.array([0.,0.,1.]),
          strides=(0,1*8),shape=(self.length,3))
      self.sinangles=bbh.sinang3D(V,zaxis) 
      self.cosangles=bbh.cosang(V,zaxis) 
      axes_raw = np.cross(V,zaxis) 
      self.axes = bbh.unit_vector(axes_raw) 
    

    #PUBLIC METHODS
    #Rotate r about the rotation axis by angle theta. 
    def rotate(self, r, t=None, debug=False):
      if t!=None:
        raise NotImplementedError("Interpolation onto new t-axis not implemented")
      else: 
        if len(r[:,0])!=self.length:
          raise ValueError("Vector must have same time sampling as basis.")
    
      x = r[:,0]
      y = r[:,1]
      z = r[:,2]

      u = self.axes[:,0]
      v = self.axes[:,1]
      w = self.axes[:,2]

      #Implicitly perform multiplication by matrix performing rotations
      #about unit vector from origin u,v,w by angle theta.
      #see inside.mines.edu/fs_home/gmurray/ArbitraryAxisRotation/ ss 5.2
      cos = self.cosangles
      sumterm = (u*x+v*y+w*z)*(1.-cos)
      sin = self.sinangles
      
      rnew = np.zeros(r.shape)
      rnew[:,0] = u*sumterm + x*cos + (-w*y + v*z)*sin
      rnew[:,1] = v*sumterm + y*cos + ( w*x - u*z)*sin
      rnew[:,2] = w*sumterm + z*cos + (-v*x + u*y)*sin

      if debug:
        rmag = la.norm(r,axis=1)
        rmag_new = la.norm(rnew,axis=1)
        err = np.sum((rmag-rmag_new)**2)
        print "(rotation debug) Total error: ", err 
        print "(rotation debug) Avg. error:", err/self.length
      return rnew

