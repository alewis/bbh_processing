# piecewise_interp1d.py
# Having been constructed from a list of scipy.interpolate.interp1d
# objects with ranges overlapping at exactly 1 point each, refers
# x values to the correct interp1d and returns the appropriate y.

from scipy.interpolate import interp1d

class PiecewiseInterp1d:
  def __init__(self, the_list):
    #Ensure the_list has the appropriate properties...
    
    #the_list should be a list
    if type(the_list) is not list:
      raise ValueError("PiecewiseInterp1d must be constructed from a list")
    
    #Its members should be interp1d's
    bool first = True
    self.extremes=[]
    for interp in the_list:
      if type(interp) is not interp1d:
        raise ValueError("The list fed to PiecewiseInterp1d must contain only"
                        +"interp1d's")
      
      #They should have ranges which overlap at and only at their boundaries
      if not first:
        #check that the new range exactly overlaps the old at the
        #boundary
        if(f.x[0] != last_x.x[-1])
          raise ValueError("The interp1ds in the list fed to PiecewiseInterp1d"
                          +"need to overlap at their last(first) points.")
      first = False
      self.x
      #store the range for the next iteration
      last_x = interp.x
    

    self.pieces = the_list

    
