#Write a header to the central file
def write_header_to_file(
    ofilename,
    headerstrings):
  the_header="" 
  for i in range(0,len(headerstrings)):
    the_header+="# ["+str(i+1)+"] = "+ headerstrings[i]+"\n" 

  print "\n"
  print "*****************************************************************"
  print "We will write the header \n %r." % the_header
  print "To file %r. " % ofilename
  print "This will overwrite the old header!"
  print "If you don't want that, hit CTRL-C (^C)."
  print "If you do want that, hit RETURN."
  print "*****************************************************************"
  raw_input("?")
  print "Very well then; overwriting the header..."
  
  #Read in the file, except for the old header.
  the_file="" 
  ofile = file(ofilename,'r')
  for line in ofile:
    if line[0] != "#":
      the_file+=line
  ofile.close()

  #Write in the new header, then the rest of the file.
  ofile=file(ofilename,'w') 
  ofile.write(the_header)
  ofile.write(the_file) 
  ofile.close()
#write_np_array_to_central_xmgr_file
#Adds the data in the numpy array 'array_to_write'
#to the file "outfilename", in an Xmgrblock-friendly
#format.

#Before doing this, makes a copy of the file
#"{filename}_COPY.dat" in the same directory,
#in case something goes wrong.

#Reads the name of the file

#This assumes the header has aleady been written,
#and does not check whether the data have a different
#number of dimensions than the array.

#central_file_location_path points to a file containing the path to the
#central file.


def write_np_array_to_central_xmgr_file(
  array_to_write,
  headerstrings,
  write_header=False,
  central_file_location_path="~/.centralfilepath"
  ):
  import numpy as np
  import argparse
  from shutil import copy2 

  
  parser=argparse.ArgumentParser(
      description=
   "Adds the data in the numpy array 'array_to_write'\n"+
   "to the file 'outfilename', in an Xmgrblock-friendly format.\n\n'"+
   "Before doing this, makes a copy of the file '{filename}_COPY.dat'\n"+
   "in the same directory, in case something goes wrong.\n\n"+
   "Does not check whether the data have a different\n"+
   " number of dimensions than the array.")


  
  
  #with open(expanded_path,'r') as f:
  #  outfilename = f.readline()

  #Write in the new data
  if (write_header): write_header_to_file(outfilename,headerstrings) 
  ofile = file(outfilename, 'a')
  np.savetxt(ofile,np.transpose(array_to_write))
  ofile.close()

