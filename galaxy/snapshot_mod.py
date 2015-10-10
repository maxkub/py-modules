from pylab        import *
from module1      import *
from binning      import *
from param_bins   import *
from load_mod     import *
from galactic_mod import *
from matplotlib.colors import LogNorm




#==============================================================================
#==============================================================================

def solar_vel_distrib(comp,t , **kwargs ):

    opt = { 'power'     : 1.0        ,
            'vmin'      : None       ,  
            'vmax'      : None       ,
            'bar_label' : 'Number'   ,
            'fsize'     : 20         ,
            'cmap'      : 'jet'      ,
            'scale'     : 'lin'      ,
            'bar_orient': 'vertical' ,
           }

    opt.update(kwargs)
    print opt

    par , attrL = read_params('input.galaxy_params','paths.galaxy_params')
    BB          = build_bins(par,attrL)
    path        = par.loadpath

    figure()
    #title(str(t)+' Gyr- '+comp,fontsize=par.fsize)
    title('time '+str(t)+'Gyr - '+comp,fontsize=par.fsize)

    ylabel(r'$\bf V\ (km/s) $',fontsize=par.fsizel)
    xlabel(r'$\bf U\ (km/s) $',fontsize=par.fsizel)



    thetab = bar_angle(t)
    
    dat=loadtxt(par.loadpath+'snapshot_'+comp+'_'+str(t)+'.res')

    x = dat[:,2]
    y = dat[:,3]

    vx = dat[:,5]
    vy = dat[:,6]
    
    rad = sqrt( x**2 + y**2 )
    
    
    rvel = - ( y*vy + x*vx )/rad  ## 'minus' to orient the positive Radial speeds to galactic center
    tvel = ( - y*vx + x*vy )/rad
    
    for i,r in enumerate(BB.rbins) :
        if r >= 8. :
            nrsun = i
            break
    
   
   
    ## Local Standard of Rest :
    rvel_sun = rvel[ (rad>=7.8) & (rad<=8.2) ]
    rvel_mean = rvel_sun.sum()/float(len(rvel_sun))
    
    tvel_sun = tvel[ (rad>=7.8) & (rad<=8.2) ]
    tvel_mean = tvel_sun.sum()/float(len(tvel_sun))
    
    rvel_sun = rvel_sun - rvel_mean
    tvel_sun = tvel_sun - tvel_mean




    bins = binning_2D( BB.rvelbins, rvel_sun, 5,1, BB.tvelbins, tvel_sun, 5,2, 1. )
   
    
       
    plot_pcolor(BB.rvelbins,BB.tvelbins,bins , invalid = 0. )




#==========================================================================



def snapshot(comp,t , **kwargs ):

    opt = { 'power'     : 1.0        ,
            'vmin'      : None       ,  
            'vmax'      : None       ,
            'bar_label' : 'Number'   ,
            'fsize'     : 20         ,
            'cmap'      : 'jet'      ,
            'scale'     : 'lin'      ,
            'bar_orient': 'vertical' ,
            'ref'       : 'bar'
           }

    opt.update(kwargs)
    

    par,attrL = read_params('input.galaxy_params','paths.galaxy_params')
    BB        = build_bins(par,attrL)
    path = par.loadpath


    
    #figure()
    #title(str(t)+' Gyr- '+comp,fontsize=par.fsize)
    #title('time '+str(t)+'Gyr - '+comp,fontsize=par.fsize)

    #ylabel(r'$\bf y\ (kpc) $',fontsize=par.fsizel)
    #xlabel(r'$\bf x\ (kpc) $',fontsize=par.fsizel)



    thetab = bar_angle(t)

    dat1=loadtxt(par.loadpath+'snapshot_'+comp+'_%3.2f.res'%(t) )

    x=dat1[:,2]
    y=dat1[:,3]


    if opt['ref'] == 'bar' :
        xr , yr = coord_rot( x , y , thetab )
        
    elif opt['ref'] == 'inertial' :
        xr = x
        yr = y
   
    bins = binning_2D(BB.rbins,xr,5,1,BB.rbins,yr,5,1,1.)

 
    plot_pcolor(BB.rbins,BB.rbins,bins , **opt )


    theta = arange(0.,2*pi,0.01)
    Rc = corot_rad(t)
    print 'Rcorot = ',Rc,' kpc'
    plot(Rc*cos(theta),Rc*sin(theta),'k')

    axes().set_aspect('equal')
    xlim(par.rmin,par.rmax)
    ylim(par.rmin,par.rmax)


#=======================================================================================================

def snapshot_mean(comp,t1,t2,dt, **kwargs):
    """  """

    opt = { 
            'rad_min' : 0.
           }

    opt.update(kwargs)
    


    par,attrL = read_params('input.galaxy_params','paths.galaxy_params')
    BB        = build_bins(par,attrL)
    path = par.loadpath

    figure()
    #title(str(t)+' Gyr- '+comp,fontsize=par.fsize)
    title('time'+str(t1)+'-'+str(t2)+'Gyr '+comp,fontsize=par.fsize)

    ylabel(r'$\bf y\ (kpc) $',fontsize=par.fsizel)
    xlabel(r'$\bf x\ (kpc) $',fontsize=par.fsizel)



    Nt = int( (t2-t1)/dt ) 
    print Nt

    bins = zeros( (len(BB.rbins),len(BB.rbins)) )

    for i in range(Nt) :

        t = t1 + i*dt
        print 't',t

        thetab = bar_angle(t)

        dat1=loadtxt(par.loadpath+'snapshot_'+comp+'_%3.2f.res' %(t))

        x=dat1[:,2]
        y=dat1[:,3]

        xr , yr = coord_rot( x , y , thetab )
  
        binst = binning_2D(BB.rbins,xr,5,1,BB.rbins,yr,5,1,1.)

        bins = bins + binst

        

    bins = bins/float(Nt)

    if opt['rad_min'] != 0. :

        k = where( (BB.rbins>-opt['rad_min']) & (BB.rbins<opt['rad_min']) )

        bins[ k[0][0]:k[0][-1]+1 , k[0][0]:k[0][-1]+1 ] = 0.


    plot_pcolor(BB.rbins,BB.rbins,bins)



    axes().set_aspect('equal')
    xlim(par.rmin,par.rmax)
    ylim(par.rmin,par.rmax)

    return BB.rbins , bins


#===========================================================================================================

def snapshot_levels(t,**kwargs):
   

    opt = { 'lv1'     : 0.03       ,
            'lv2'     : 30.        ,
            'nlevels' : 200        ,
            'comp'    : 'stars'    ,
            'ref'     : 'inertial' ,
            'cmap'    : 'jet'
           }

    opt.update(kwargs)


    lv1     = opt['lv1']
    lv2     = opt['lv2']
    comp    = opt['comp']
    nlevels = opt['nlevels']

    par,attrL = read_params('input.galaxy_params','paths.galaxy_params')
    BB        = build_bins(par,attrL)

    path=par.loadpath
    fsizel = par.fsizel



    thetab = bar_angle(t)

    #-------------------------------------------------------------------------

    dat1=loadtxt(path+'snapshot_'+comp+'_%3.2f.res'%(t) )
    

    pos = dat1[:,2:4]

    if opt['ref'] == 'bar' :
        x_rot,y_rot = coord_rot( pos[:,0] , pos[:,1] , thetab )
        
    elif opt['ref'] == 'inertial' :
        x_rot = pos[:,0]
        y_rot = pos[:,1]
        
    
    N_bins = binning_2D(BB.rbins,x_rot,7,1,BB.rbins,y_rot,7,1,1.)


    N_bins[ N_bins >= N_bins.max()*lv1 ] = nan  # at 2Gyr # lv1 = 0.03  # at 2.5Gyr  lv1 = 0.015
    N_bins[ N_bins <= lv2 ] = nan                         # lv2 = 30.                lv2 = 20.

    N_bins = ma.masked_invalid(N_bins)
    #cmap = get_cmap('YlGn')
    #cmap = get_cmap('Reds')
    cmap = get_cmap(opt['cmap'])
    cmap.set_bad(color = 'w', alpha = 1.)



    def func(x,y):
      return N_bins[x,y]

    I=range(len(BB.rbins))
    J=range(len(BB.rbins))
    X,Y = meshgrid(I, J)



    Z = func(X, Y)

    #PP = pcolor(r_bins[X], r_bins[Y], Z, norm=LogNorm(vmin=Z.min(), vmax=Z.max()) , cmap = cmap )#, vmin=0.0 ,v)
    PP = pcolor(BB.rbins[X], BB.rbins[Y], Z, norm=LogNorm(vmin=1., vmax=Z.max()) , cmap = cmap )#, vmin=0.0 ,v)
    cb = colorbar(PP, orientation='vertical')
    cb.set_label(r'$\bf Number $',fontsize=fsizel)


#================================================================================================

def snapshot_met(par,BB,comp,t,omega):
    """ see shape, comp = 'stars' or 'gas' """

    figure()
    #title(str(t)+' Gyr - '+comp,fontsize=18)
    title( 'time step '+str(t)+' - '+comp, fontsize = par.fsize )

    ylabel( r'$\bf y\ (kpc) $', fontsize = par.fsizel )
    xlabel( r'$\bf x\ (kpc) $', fontsize = par.fsizel )



    dat=loadtxt( par.loadpath+'snapshot_'+comp+'_'+str(t)+'.res' )

    pos  = dat[:,2:4]
    met  = dat[:,8]
    mass = dat[:,9]
   
    uu = met * mass


    bins1 = binning_2D( BB.rbins, pos[:,0], 5, 1, BB.rbins, pos[:,1], 5, 1, uu)
    bins2 = binning_2D( BB.rbins, pos[:,0], 5, 1, BB.rbins, pos[:,1], 5, 1, mass)


    bins = 0.*bins1

    for i in xrange(len(bins[:,0])):
        for j in xrange(len(bins[0,:])):

            if (bins2[i,j] != 0.):
                bins[i,j] = bins1[i,j] / bins2[i,j]
            else:
                bins[i,j] = 0.

            if (bins[i,j] != 0. ):
                bins[i,j] = log10(bins[i,j]/0.0142)
            else:
                bins[i,j] = -1000. 



    def func(x,y):
      return bins[x,y]

    I=range(len(BB.rbins))
    J=range(len(BB.rbins))
    X,Y = meshgrid(I, J)


    #Z = (func(X, Y))**(power)
    Z = func(X, Y)
    pcolor(BB.rbins[X], BB.rbins[Y], Z, vmin=-1.5, vmax=1. )
    cb = colorbar(format='%.1f')
    cb.set_label(r'$\bf [O/H] $',fontsize=16)



    theta = arange(0.,2*pi,0.01)
    R = 220./omega
    print 'Rcorot = ',R,' kpc'
    plot(R*cos(theta),R*sin(theta),'k')



    axes().set_aspect('equal')
    xlim(par.rmin,par.rmax)
    ylim(par.rmin,par.rmax)



#================================================================================================


