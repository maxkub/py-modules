#!/usr/bin/env python

"""param_bins.py: - load parameters and their name with : par,attrL = read_params()
                  - build bins : BB = build_bins(par,attrL) 
                  - find_file  : findfile('name')  """


from pylab import *
from os import path



#==========================================================================
#==========================================================================

import os, sys
 
 
 
 
#============================================================================================================
 
def findfile(path):
    """Find the file named path in the sys.path.
       Returns the full path name if found, None if not found
    """ 
    for dirname in sys.path:
        possible = os.path.join(dirname, path)
        if os.path.isfile(possible):
            return possible
    return None


#==========================================================================
#==========================================================================


class Variables():
    """ contain all the variables describing particles moves """

    pass


#============================================================================================


def read_params(input_name,input_path):
    """ read parameters in input.params """

    home=path.expanduser('~/')
    
    param = findfile(input_name)
    ppath = findfile(input_path) 
    
    #print param
    #print ppath
    
    dat      = loadtxt( param )
    datpaths = loadtxt( ppath ,dtype='string')


    par = Variables()

    attrList=['nthing','thingmin','thingmax',
        'nrad' ,'rmin','rmax',
        'nDrad','Dradmin','Dradmax',
        'ntheta','theta_min','theta_max',
        'nz','zmin','zmax',
        'nt','tmin','tmax',
        'nL','Lmin','Lmax',
        'nDL','DLmin','DLmax',
        'nvelt','velmin_t','velmax_t',
        'nvelr','velmin_r','velmax_r',
        'nvelz','velmin_z','velmax_z',
        'nomega','omegmin','omegmax',
        'met','metmin','metmax',
        'ntime',
        'nmet',
        'nalpha',
        'nage',
        'power',
        'save',
        'fsize',
        'fsizel',
        'theFiles',
        ]




    k = 0
    for att in attrList :
        setattr( par , att , dat[k] )
        k = k + 1

    

    par.colorv    = ['b','g','r','c','m','y','k']
    par.loadpath  = datpaths[0]
    par.savepath  = home + datpaths[1]


    #bins = Bins(par,attrList)

    rcParams['font.size'] = par.fsize

    return par , attrList

#=====================================================================================

def build_bins(par,parList):
    """ create class containing the different bins """


    #print " create bins"

    bins = Variables()

    binsList=['thing',
        'rbins', 
        'Drbins',
        'thetabins',
        'zbins',
        'Tbins',
        'Lbins',
        'DLbins',
        'tvelbins',
        'rvelbins',
        'zvelbins',
        'omegbins',
        'metbins'
        ]


    i = 0
    for bi in binsList :

        N   = int( getattr( par , parList[i] ) )
        Min = getattr( par , parList[i+1] )
        Max = getattr( par , parList[i+2] )

        ex = zeros( N + 1)
        for j in range( N + 1):
            ex[j] = Min + (Max - Min)/float(N)*float(j) 

        i = i + 3

        setattr( bins , bi , ex ) 


    return bins


#=====================================================================================
#=====================================================================================






