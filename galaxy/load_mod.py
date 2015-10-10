from pylab import *
from galactic_mod import omeg_bar, omeg_mean




"""
============================================================================================
This module is used to load data from GADGET analysis scripts,
and treat the data to get information describing particles cinematics
and dynamics.

to be used in angmom* scripts

============================================================================================
"""






#============================================================================================
#============================================================================================

class Variables():
    """ contain all the variables describing particles moves """

    pass

#============================================================================================


def select_load(par):


    N = par.theFiles


    dat = Variables()


    if N == 1 :
        dat.dat1 = loadtxt(par.loadpath+'migrangm3_starsr0.0-30.5t1.0-2.0Th200.0_1.res')
        dat.dat2 = loadtxt(par.loadpath+'migrangm3_starsr0.0-30.5t1.0-2.0Th200.0_2.res')
        dat.dat3 = loadtxt(par.loadpath+'deltas_starsr0.0-30.5t1.0-2.0.res')

    if N == 2 :
        dat.dat1 = loadtxt(par.loadpath+'migrangm3_starsr0.0-30.5t2.0-3.0Th200.0_1.res')
        dat.dat2 = loadtxt(par.loadpath+'migrangm3_starsr0.0-30.5t2.0-3.0Th200.0_2.res')
        dat.dat3 = loadtxt(par.loadpath+'deltas_starsr0.0-30.5t2.0-3.0.res')

    if N == 3 :
        dat.dat1 = loadtxt(par.loadpath+'migrangm3_starsr0.0-30.5t4.0-5.0Th200.0_1.res')
        dat.dat2 = loadtxt(par.loadpath+'migrangm3_starsr0.0-30.5t4.0-5.0Th200.0_2.res')
        dat.dat3 = loadtxt(par.loadpath+'deltas_starsr0.0-30.5t4.0-5.0.res')

    if N == 4 :
        dat.dat1 = loadtxt(par.loadpath+'migrangm3_starsr0.0-30.5t9.0-10.0Th200.0_1.res')
        dat.dat2 = loadtxt(par.loadpath+'migrangm3_starsr0.0-30.5t9.0-10.0Th200.0_2.res')
        dat.dat3 = loadtxt(par.loadpath+'deltas_starsr0.0-30.5t9.0-10.0.res')

    if N == 5 :
        dat.dat1 = loadtxt(par.loadpath+'migrangm3_starsr0.0-30.5t4.0-5.0Th400.0_1.res')
        dat.dat2 = loadtxt(par.loadpath+'migrangm3_starsr0.0-30.5t4.0-5.0Th400.0_2.res')
        dat.dat3 = loadtxt(par.loadpath+'deltas_starsr0.0-30.5t4.0-5.0.res')

    if N == 6 :
        dat.dat1 = loadtxt(par.loadpath+'migrangm3_starsr0.0-30.5t4-5_1.res')
        dat.dat2 = loadtxt(par.loadpath+'migrangm3_starsr0.0-30.5t4-5_2.res')
        dat.dat3 = loadtxt(par.loadpath+'deltas_starsr0.0-30.5t4.0-5.0.res')




    print " load done "
    
    var = data_treatment(dat)

    return var

#=================================================================================
#=================================================================================



def select_load2(par,t1,t2,Thres):

    N = par.theFiles

    par.time1 = t1
    par.time2 = t2 


    dat = Variables()



    dat.dat1 = loadtxt(par.loadpath+'migrangm3_starsr0.0-30.5t'+str(t1)+'-'+str(t2)+'Th'+str(Thres)+'_1.res')
    dat.dat2 = loadtxt(par.loadpath+'migrangm3_starsr0.0-30.5t'+str(t1)+'-'+str(t2)+'Th'+str(Thres)+'_2.res')
    dat.dat3 = loadtxt(par.loadpath+'deltas_starsr0.0-30.5t'+str(t1)+'-'+str(t2)+'.res')



    print " load done "
    
    var = data_treatment(dat)

    var.omeg_bar1 = omeg_bar(par,t1)
    var.omeg_bar2 = omeg_bar(par,t2)
    var.omeg_mean = omeg_mean(par,t1,t2)
    


    return var

#=================================================================================
#=================================================================================



def data_treatment(dat) :


    var = Variables()


    var.idmig      = dat.dat1[:,1]
    var.age        = dat.dat1[:,2]
    var.btime      = dat.dat1[:,3]
    var.agemig     = dat.dat1[:,4]

    var.posmig     = dat.dat1[:,5:7]
    var.velmig     = dat.dat1[:,8:10]

    var.angm_migr  = var.posmig[:,0]*var.velmig[:,1] - var.posmig[:,1]*var.velmig[:,0]

    var.posinit    = dat.dat2[:,2:4]
    var.posfin     = dat.dat2[:,5:7]

    var.velinit    = dat.dat2[:,8:10]
    var.velfin     = dat.dat2[:,11:13]

    var.angminit   = var.posinit[:,0]*var.velinit[:,1] - var.posinit[:,1]*var.velinit[:,0] 
    var.angmfin    = var.posfin[:,0]*var.velfin[:,1] - var.posfin[:,1]*var.velfin[:,0] 

    var.delta_angm = var.angmfin - var.angminit
    var.delta_angm_migr = var.angm_migr - var.angminit


    var.rmig       = sqrt( var.posmig[:,0]**2 + var.posmig[:,1]**2 ) 

    var.ri         = sqrt( var.posinit[:,0]**2 + var.posinit[:,1]**2 ) 
    var.rf         = sqrt( var.posfin[:,0]**2 + var.posfin[:,1]**2 )

    var.delta_r    = var.rf - var.ri

    var.delta_gc   = var.ri/var.angminit*var.delta_angm

    var.mtime      = var.agemig + var.btime

    # tangential velocities
    #var.tveli     = (posinit[:,1]*velinit[:,0] + posinit[:,0]*velinit[:,1])/ri[:]
    #var.tvelf     = (posfin[:,1]*velfin[:,0] + posfin[:,0]*velfin[:,1])/rf[:]

    var.tveli      = (-var.posinit[:,1]*var.velinit[:,0] + var.posinit[:,0]*var.velinit[:,1])/var.ri[:]
    var.tvelf      = (-var.posfin[:,1]*var.velfin[:,0] + var.posfin[:,0]*var.velfin[:,1])/var.rf[:]
    var.tvelmig    = (-var.posmig[:,1]*var.velmig[:,0] + var.posmig[:,0]*var.velmig[:,1])/var.rmig[:]


    var.omegai     = var.tveli / var.ri
    var.omegaf     = var.tvelf / var.rf

    var.delta_vt   = var.tvelf - var.tveli

    #var.delta_gc2 = delta_angm/tveli - ri/(tveli)*delta_vt
    var.delta_gc2  = var.delta_angm/220.

    # radial velocities
    #var.rveli     = (posinit[:,1]*velinit[:,1] - posinit[:,0]*velinit[:,0])/ri[:]
    #var.rvelf     = (posfin[:,1]*velfin[:,1] - posfin[:,0]*velfin[:,0])/rf[:]

    var.rveli      = (var.posinit[:,1]*var.velinit[:,1] + var.posinit[:,0]*var.velinit[:,0])/var.ri[:]
    var.rvelf      = (var.posfin[:,1]*var.velfin[:,1] + var.posfin[:,0]*var.velfin[:,0])/var.rf[:]
    var.rvelmig    = (var.posmig[:,1]*var.velmig[:,1] + var.posmig[:,0]*var.velmig[:,0])/var.rmig[:]





    #### for all stars######################################################################
    var.idtot      = dat.dat3[:,1]

    var.positot    = dat.dat3[:,2:4]
    var.posftot    = dat.dat3[:,5:7]

    var.ritot      = sqrt( var.positot[:,0]**2 + var.positot[:,1]**2 ) 
    var.rftot      = sqrt( var.posftot[:,0]**2 + var.posftot[:,1]**2 )

    var.velitot    = dat.dat3[:,8:10]
    var.velftot    = dat.dat3[:,11:13]

    var.angmitot   = dat.dat3[:,14]
    var.angmftot   = dat.dat3[:,15]
    #var.angmftot  = posftot[:,0]*velftot[:,1]-posftot[:,1]*velftot[:,0]

    var.delta_angmtot = var.angmftot - var.angmitot

    var.delta_rtot    = var.rftot - var.ritot

    var.delta_gctot   = var.ritot#/angmitot#*delta_angmtot


    # tangential velocities
    var.tvelitot      = (-var.positot[:,1]*var.velitot[:,0] + var.positot[:,0]*var.velitot[:,1])/var.ritot[:]
    var.tvelftot      = (-var.posftot[:,1]*var.velftot[:,0] + var.posftot[:,0]*var.velftot[:,1])/var.rftot[:]


    var.omegaitot = var.tvelitot / var.ritot
    var.omegaftot = var.tvelftot / var.rftot

    var.delta_vttot   = var.tvelftot - var.tvelitot

    var.delta_gc2tot  = var.delta_angmtot/var.tvelitot - var.ritot/(var.tvelitot)*var.delta_vttot

    # radial velocities
    var.rvelitot      = (var.positot[:,1]*var.velitot[:,1] + var.positot[:,0]*var.velitot[:,0])/var.ritot[:]
    var.rvelftot      = (var.posftot[:,1]*var.velftot[:,1] + var.posftot[:,0]*var.velftot[:,0])/var.rftot[:]





    print " treatment done "

    return var


#==============================================================================================
#==============================================================================================
