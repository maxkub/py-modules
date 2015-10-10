from pylab import * 
from module1 import *
from scipy.interpolate import griddata







#===========================================================================================================

def effective_pot(par,comp,t,lv1,lv2):
    """ show effective potential in frame rotating with the bar """

    fig = figure()

    title(str(t)+'Gyr')


    dat = loadtxt(par.loadpath+'potential_halo,disk,stars_'+str(t)+'Gyr200.res')
    #dat = loadtxt(path+'potential_halo,disk,stars_'+str(t)+'Gyr400.res')

    # cartesian coordinate
    xdat = dat[0,3:]
    ydat = xdat.copy()

    pot = - dat[ 1 : ,3:]



    N = len(xdat)



    def ang_speed(t):
        # from analyse with EUREKA software
        return (154.8 - t)/(t + 2.508**(1.107 - (0.2419*t)**((0.5799*t)**(10.3/t)) - 0.513*t))


    eff_pot = 0. * pot

    # effective potential Binney & tremaine : pot(x,y) - 0.5*omega_bar**2 * (x**2 + y**2)
    for i in range(N-1):
        for j in range(N-1):
            # using the unit mass in the simulation K = 10e10 Msun --> K*G = 4.498e4 kpc3.Gyr-2
            eff_pot[i,j] = pot[i,j]*4.498e4 - 0.5*ang_speed(t)**2*(xdat[i]**2 + ydat[j]**2) 



    """
    dat = loadtxt(par.loadpath+'bar_angle_0.:10._2.res')

    tb     = dat[:,1]
    angleb = dat[:,2]

    k=where( (tb > t-0.01) & (tb<t+0.01) )
    angle = angleb[k].sum()/float(len(angleb[k]))

  
    print 'angle',angle , angle%2.*3.14

    x =  cos( angle ) *xdat + sin( angle ) * ydat
    y = -sin( angle ) *xdat + cos( angle ) * ydat



    points = []
    values = []
    grid = []
    grid_y = []
    for i in range(len(x)):

        gx = ones( len(x) )
        gx = gx * xdat[i]
        points = append( points , gx , 0 )
        points = append( points , ydat , 1 )

        gx = ones( len(x) )
        gx = gx * x[i]
        grid = append( grid , gx , 0 )
        grid = append( grid , ydat , 1 )

        for j in range(len(x)):

            values.append( eff_pot[i,j] )




    grid_pot_rot = griddata(points, values, grid , method='cubic')

    print grid_pot_rot
    """



    #value = -290000#-290000.
    #eff_pot[ eff_pot<value ] = nan

    #eff_pot = ma.masked_invalid(eff_pot)
    cmap = get_cmap('jet')
    cmap.set_bad(color = 'w', alpha = 1.)

    X, Y = meshgrid(xdat, ydat)

    I=range(N-1)
    J=range(N-1)
    X,Y = meshgrid(I, J)


    def func(x,y):
        return eff_pot[x,y]

    Z = (func(X, Y))**(1.)
    #pcolor(x[X], x[Y], Z)#, vmin=0.0 ,v)
    contour(xdat[X], ydat[Y], Z, 400, cmap=cmap)#, vmin=-290000.0 )#,v)
    cb = colorbar(format='%.0f')
    cb.set_label(r'$\bf potential\ (kpc^{2}.Gyr^{-2}) $',fontsize=16)



    #----------------------------------------------------------------------
    """

    dat1=loadtxt(path+'snapshot_'+comp+'_'+str(t)+'.res')#_codallcomptest.res')
    

    pos = dat1[:,2:4]



    r_bins=zeros( nrad +1)

    for i in range(nrad+1):
        r_bins[i] = rmin + (rmax - rmin)/float(nrad)*float(i) 



    N_bins = binning_2D(r_bins,pos[:,1],10,1,r_bins,pos[:,0],10,1,1.)


    def func(x,y):
      return N_bins[x,y]

    I=range(len(r_bins))
    J=range(len(r_bins))
    X,Y = meshgrid(I, J)



    Z = func(Y, X)
    #PP = pcolor(r_bins[X], r_bins[Y], Z, norm=LogNorm(vmin=1., vmax=Z.max()) )#, vmin=0.0 ,v)
    #cb = colorbar(PP, orientation='vertical')
    #cb.set_label(r'$\bf Number $',fontsize=fsizel)
    
    PP = contour(r_bins[X], r_bins[Y], Z, levels=[lv1,lv2],colors = 'k' )  # levels = [50,200]



    theta = arange(0.,2*pi,0.01)
    R = 220./ang_speed(t)
    print 'Rcorot = ',R,' kpc'
    plot(R*cos(theta),R*sin(theta),'k--')
   

    grid('on')

    xlim(rmin,rmax)
    ylim(rmin,rmax)
    axes().set_aspect('equal')

    """


