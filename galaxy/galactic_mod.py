from pylab      import *
from module1    import *
from param_bins import *
from matplotlib.colors import LogNorm






#===============================================================================

#===============================================================================

def bar_radius(tt):
    """ determine the radius of the bar at time t"""

    par,attrL = read_params('input.galaxy_params','paths.galaxy_params')
   

    dat = loadtxt(par.loadpath+'bar_radius_0.0:10.0.res')

    n = dat[:,0]
    t = dat[:,1]
    r = dat[:,2]

    for i in range(len(t)):
        if t[i]>= tt :
            rad = r[i]
            break

    #rad  = r[ t == tt ]


    return rad

#====================================================================================
#====================================================================================

def omeg_bar(par,tt):
    """ determine the angular speed of the bar at time t"""

    par,attrL = read_params('input.galaxy_params','paths.galaxy_params')



    dat = loadtxt(par.loadpath+'bar_angle_0.:10._2.res')

    n = dat[:,0]
    t = dat[:,1]
    a = dat[:,2]


    dt = t[1]-t[0]

    omeg = []
    tv   = []
    nv   = []
    for i in range(len(t)-3):
        i = i + 1
        omeg.append( (a[i+1] - a[i-1]) / (2.*dt) )
        
        tv.append(t[i])
        nv.append(n[i-1])


    omeg = asarray(omeg)
    tv   = asarray(tv)
    nv   = asarray(nv)

    speed  = omeg[ tv == tt ]

    return speed

#==================================================================================
#==================================================================================

def omeg_mean(par,t1,t2):
    """ determine the mean angular speed of the bar between times t1 and t2"""


    par,attrL = read_params('input.galaxy_params','paths.galaxy_params')



    dat = loadtxt(par.loadpath+'bar_angle_0.:10._2.res')

    n = dat[:,0]
    t = dat[:,1]
    a = dat[:,2]


    dt = t[1]-t[0]

    omeg = []
    tv   = []
    nv   = []
    for i in range(len(t)-3):
        i = i + 1
        omeg.append( (a[i+1] - a[i-1]) / (2.*dt) )
        
        tv.append(t[i])
        nv.append(n[i-1])


    omeg = asarray(omeg)
    tv   = asarray(tv)
    nv   = asarray(nv)

    tot = omeg[ (tv>=t1) & (tv<=t2) ].sum()
    N   = nv[ (tv>=t1) & (tv<=t2) ]

    N1   = N[0]
    N2   = N[-1]

    print N1,N2,N2-N1,tot

    mean = tot/float(N2-N1)
    print 'mean=',mean

    return mean

#==================================================================================
#==================================================================================


def corot_rad(t) :
    """ compute corotation radius at time t """

    par,attrL = read_params('input.galaxy_params','paths.galaxy_params')
    path = par.loadpath

    dat     = loadtxt(path+'average_tang_speed.res')
    ta      = dat[:,0]
    tvel    = dat[:,1]

    angl   = loadtxt(par.loadpath+'bar_angle_0.:10._2.res')
    to     = angl[:,1]
    thetav = angl[:,2]


    speed = tvel[ ta >= t ]
 
    Rco = speed[0]/ang_speed(t)

 
    return Rco




#==================================================================================
#==================================================================================




def bar_angle(t) :
    """compute bar radius at time t """


    par,attrL = read_params('input.galaxy_params','paths.galaxy_params')
    path = par.loadpath


    dat     = loadtxt(path+'average_tang_speed.res')
    ta      = dat[:,0]
    tvel    = dat[:,1]

    angl   = loadtxt(par.loadpath+'bar_angle_0.:10._2.res')
    to     = angl[:,1]
    thetav = angl[:,2]


    for i,tt in enumerate(to):
        if ( tt >= t ) :
            thetab = thetav[i]
            break


    return thetab



#==================================================================================
#==================================================================================



def ang_speed(t):
    """ compute angular speed of the bar at time t
        expression from analyse with EUREKA software
    """
    return (154.8 - t)/(t + 2.508**(1.107 - (0.2419*t)**((0.5799*t)**(10.3/t)) - 0.513*t))

#==================================================================================
#==================================================================================

def coord_rot(x,y,theta):
    """ compute the rotated coordinates x and y by angle thetab """

    x_rot =  cos(theta) *x + sin( theta ) * y
    y_rot = -sin(theta) *x + cos( theta ) * y

    return x_rot , y_rot



#==================================================================================
#==================================================================================

def resonance_order_m(m , **kwargs ):
    """ compute the radius of the resonance of order m, as a function of time 
        Or at a specific time
    
    """
    
    opt = { 'comp'  : 'stars'    ,
            'ncell' : '200'      ,
            'time'  : None                
           }
           
    opt.update(kwargs)

    
    par,attrL = read_params('input.galaxy_params','paths.galaxy_params')
    path      = par.loadpath


    #inputfile = 'rotfreq_'+opt['comp']+'_'+opt['ncell']+'_kappa2'  
    #inputfile = 'rotfreq_1001_wins35_kappa2'  
    inputfile = 'kappa2_lia_wins11'
    dat1  = loadtxt( path + inputfile +'.pyres')
    rad  = dat1[0, 1:]    
    time = dat1[1:, 0]
    
    
    #inputfile = 'rotfreq_'+opt['comp']+'_'+opt['ncell']+'_smooth'  
    #inputfile = 'rotfreq_1001_wins35'  
    inputfile = 'rotfreq_lia'  
    dat2  = loadtxt( path + inputfile +'.pyres')

        
    
    rr = []
    tt = []

    for i in xrange(len(time)):
        
        #print dat2[i+1,2:]
        #print dat1[i+1,1:]
            
        diff = dat2[i+1,2:]- ang_speed(time[i]) + sqrt(dat1[i+1,1:])/m
        
        diff0 = diff[0]
        for j in xrange(len(diff) -1 ):
            diff1 = diff[j+1]

            if diff1/diff0 <= 0. :
                rr.append( (rad[j+1]+rad[j])/2. )
                tt.append( time[i] )
                break

            diff0 = diff1
        
        
    rr = asarray( rr )
    tt = asarray( tt )
    
    if opt['time'] != None :
        for i in range(len(tt)):
            if tt[i] >= opt['time'] :
                st = tt[i]
                sr = rr[i]
                break
                
    if opt['time'] == None :
        return tt, rr
    else:
        return st, sr

#==================================================================================
#==================================================================================

def average_vtang(t) :
    """ compute average tangential speed of stars in the disk at time t """

    par,attrL = read_params('input.galaxy_params','paths.galaxy_params')
    path = par.loadpath

    dat     = loadtxt(path+'average_tang_speed.res')
    ta      = dat[:,0]
    tvel    = dat[:,1]

 
    for i in range(len(ta)):
        if ta[i] >= t :
            st = ta[i]
            sv = tvel[i]
            break
   
 
    return st , sv

#==================================================================================
#==================================================================================

def vcirc(t , r) :
    """ gives the circular velocity of stars at given radius and time
        Uses the file from Lia
    """
    par,attrL = read_params('input.galaxy_params','paths.galaxy_params')
    path = par.loadpath
    
    
    dat = loadtxt(path+'vcirc_lia.pyres')
    
    rad  = dat[0, 2:]    
    time = dat[1:, 1]
    
    vcirc = dat[1:,2:]
    
    for i in range(len(time)):

        if time[i]>= t :
            
            for j in range(len(rad)):
                
                if rad[j] >= r :
                    
                    sr = rad[j]
                    st = time[i]
                    svcirc = vcirc[i,j]
                    
                    break
                    
                    
    return st, sr, svcirc
    
#==================================================================================
#==================================================================================

def vcirc2(t , r) :
    """ gives the circular velocity of stars at given radius and time
        Uses the file from Lia
        
        !!!!!!!   IS NOT FASTER THAN vcirc() !!!!!!!!!!!!!!!!!!
        
    """
    par,attrL = read_params('input.galaxy_params','paths.galaxy_params')
    path = par.loadpath
    
    
    dat = loadtxt(path+'vcirc_lia.pyres')
    
    rad  = dat[0, 2:]    
    time = dat[1:, 1]
    
    vcirc = dat[1:,2:]
    
    index_t = bin_index(time,6,1,t)
    index_r = bin_index(rad,6,1,r+rad[0])
    
    sr = rad[index_r]
    st = time[index_t]
    svcirc = vcirc[index_t,index_r]
                    
                    
    return st, sr, svcirc
    
#==================================================================================================
#==================================================================================================
    
def get_cod(ain,acomp,time):
      
    """
    extract the time and position and velocity of the mass center of the galaxy

    inputs
    character(len=80),intent(in) :: ain, acomp
    real(kind=wp),intent(in)     :: time

    !! outputs
    real(kind=wp),dimension(7)   :: cod
    """ 

    if ain == '/nethome/kubryk/works/uns_stuffs/list001.txt' :

        path = '/home/picardan1NS/kubryk/mwg/mwg001/analysis/'

        dat = loadtxt( path + 'mwg001.'+acomp+'.cod')

    if ain == '/nethome/kubryk/works/uns_stuffs/list003.txt' :

        path = '/home/picardan1NS/kubryk/mwg/mwg003/analysis/'

        dat = loadtxt( path + 'mwg003.'+acomp+'.cod')


    i=0
    while  1 :
        i=i+1

        if  i > 2001  : 
            print '>> ERROR in get_cod <<'

        if  dat.item(i,0) >= round(time,4)  : 
            print 'time cod', time, dat.item(i,0)
            break
                 
            
    return dat[i,:]


                
    


