from module1    import *
from param_bins import *


"""
#===================================================================================================

Initializing the objects  par, containing all the parameters written in '~/py-modules/input.params' .
                          BB , containing the bins used in the different data processings           .

#===================================================================================================
"""


par,BB,attrList = read_params()



print ""

print "------parameters----------"
print "par.    "
print dir(par)

print ""

print "------bins----------------"
print "BB.    "
print dir(BB)


