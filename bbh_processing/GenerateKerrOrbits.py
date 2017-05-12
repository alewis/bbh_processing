#This translates Aaron's mathematica notebook into Python.


#Dimensionless functions
#Various thins for convenience.


#Converting Conserved Quantities
#Functions to transform (p,ecc,zm) into (E,Lz,Q).

def Lz(p,ecc,am,a):
  r1 = p /(1. + ecc)
  Enx = En(p,ecc,zm,a)


