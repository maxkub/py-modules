from module1    import *
from param_bins import *
from load_mod   import *


"""
#===================================================================================================
                                                                                                    
Initializing the objects  par, containing all the parameters written in '~/py-modules/input.params' .
                          BB , containing the bins used in the different data processings           .
                          var, containing the physical variables extracted from the GADGET data     . 
                                                                                                    
#===================================================================================================
"""



def initialize(t1,t2,Thres):
    par,BB,attrList = read_params()

    #var    = select_load(par)
    var    = select_load2(par,t1,t2,Thres)




    print ""

    print "------parameters----------"
    print "par.    "
    print dir(par)

    print ""

    print "------bins----------------"
    print "BB.    "
    print dir(BB)

    print ""

    print "------variables-----------"
    print "var.    "
    print dir(var)

    return par,BB,var




