
from scipy.optimize    import leastsq
from pylab             import *
from matplotlib.colors import LogNorm
from binning           import *





#==============================================================================================================

def color_cycle(n , **kwargs) :
    
    opt = { 'colors'    : ['b' , 'g' , 'r' , 'c' , 'm' , 'y' , 'k' ]  }

    opt.update(kwargs)
    
    colors = opt['colors']
    
    return colors[int(n%len(colors))]



#==============================================================================================================
# from scipy cookbook 


def savitzky_golay(y, window_size, order, deriv=0, rate=1):
    r"""Smooth (and optionally differentiate) data with a Savitzky-Golay filter.
    The Savitzky-Golay filter removes high frequency noise from data.
    It has the advantage of preserving the original shape and
    features of the signal better than other types of filtering
    approaches, such as moving averages techniques.
    Parameters
    ----------
    y : array_like, shape (N,)
        the values of the time history of the signal.
    window_size : int
        the length of the window. Must be an odd integer number.
    order : int
        the order of the polynomial used in the filtering.
        Must be less then `window_size` - 1.
    deriv: int
        the order of the derivative to compute (default = 0 means only smoothing)
    Returns
    -------
    ys : ndarray, shape (N)
        the smoothed signal (or it's n-th derivative).
    Notes
    -----
    The Savitzky-Golay is a type of low-pass filter, particularly
    suited for smoothing noisy data. The main idea behind this
    approach is to make for each point a least-square fit with a
    polynomial of high order over a odd-sized window centered at
    the point.
    Examples
    --------
    t = np.linspace(-4, 4, 500)
    y = np.exp( -t**2 ) + np.random.normal(0, 0.05, t.shape)
    ysg = savitzky_golay(y, window_size=31, order=4)
    import matplotlib.pyplot as plt
    plt.plot(t, y, label='Noisy signal')
    plt.plot(t, np.exp(-t**2), 'k', lw=1.5, label='Original signal')
    plt.plot(t, ysg, 'r', label='Filtered signal')
    plt.legend()
    plt.show()
    References
    ----------
    .. [1] A. Savitzky, M. J. E. Golay, Smoothing and Differentiation of
       Data by Simplified Least Squares Procedures. Analytical
       Chemistry, 1964, 36 (8), pp 1627-1639.
    .. [2] Numerical Recipes 3rd Edition: The Art of Scientific Computing
       W.H. Press, S.A. Teukolsky, W.T. Vetterling, B.P. Flannery
       Cambridge University Press ISBN-13: 9780521880688
    """
    import numpy as np
    from math import factorial

    try:
        window_size = np.abs(np.int(window_size))
        order = np.abs(np.int(order))
    except ValueError, msg:
        raise ValueError("window_size and order have to be of type int")
    if window_size % 2 != 1 or window_size < 1:
        raise TypeError("window_size size must be a positive odd number")
    if window_size < order + 2:
        raise TypeError("window_size is too small for the polynomials order")
    order_range = range(order+1)
    half_window = (window_size -1) // 2
    # precompute coefficients
    b = np.mat([[k**i for i in order_range] for k in range(-half_window, half_window+1)])
    m = np.linalg.pinv(b).A[deriv] * rate**deriv * factorial(deriv)
    # pad the signal at the extremes with
    # values taken from the signal itself
    firstvals = y[0] - np.abs( y[1:half_window+1][::-1] - y[0] )
    lastvals = y[-1] + np.abs(y[-half_window-1:-1][::-1] - y[-1])
    y = np.concatenate((firstvals, y, lastvals))
    return np.convolve( m[::-1], y, mode='valid')




#==============================================================================================================
# from scipy cookbook 


def sgolay2d ( z, window_size, order, derivative=None):
    """
    hx, hy added by max, is the interval in x and y
    """
    # number of terms in the polynomial expression
    n_terms = ( order + 1 ) * ( order + 2)  / 2.0

    if  window_size % 2 == 0:
        raise ValueError('window_size must be odd')

    if window_size**2 < n_terms:
        raise ValueError('order is too high for the window size')

    half_size = window_size // 2

    # exponents of the polynomial. 
    # p(x,y) = a0 + a1*x + a2*y + a3*x^2 + a4*y^2 + a5*x*y + ... 
    # this line gives a list of two item tuple. Each tuple contains 
    # the exponents of the k-th term. First element of tuple is for x
    # second element for y.
    # Ex. exps = [(0,0), (1,0), (0,1), (2,0), (1,1), (0,2), ...]
    exps = [ (k-n, n) for k in range(order+1) for n in range(k+1) ]

    # coordinates of points
    ind = arange(-half_size, half_size+1, dtype=float64)
    dx = repeat( ind, window_size )
    dy = tile( ind, [window_size, 1]).reshape(window_size**2, )

    # build matrix of system of equation
    A = empty( (window_size**2, len(exps)) )
    for i, exp in enumerate( exps ):
        A[:,i] = (dx**exp[0]) * (dy**exp[1])

    # pad input array with appropriate values at the four borders
    new_shape = z.shape[0] + 2*half_size, z.shape[1] + 2*half_size
    Z = zeros( (new_shape) )

    """
    # top band
    band = z[0, :]
    Z[:half_size, half_size:-half_size] =  band -  abs( flipud( z[1:half_size+1, :] ) - band )
    # bottom band
    band = z[-1, :]
    Z[-half_size:, half_size:-half_size] = band  + abs( flipud( z[-half_size-1:-1, :] )  -band )
    # left band
    band = tile( z[:,0].reshape(-1,1), [1,half_size])
    Z[half_size:-half_size, :half_size] = band - abs( fliplr( z[:, 1:half_size+1] ) - band )
    # right band
    band = tile( z[:,-1].reshape(-1,1), [1,half_size] )
    Z[half_size:-half_size, -half_size:] =  band + abs( fliplr( z[:, -half_size-1:-1] ) - band )
    # central band
    Z[half_size:-half_size, half_size:-half_size] = z
    """

    # top band
    band = z[0, :]
    Z[:half_size, half_size:-half_size] =  band + band -  flipud( z[1:half_size+1, :] )
    # bottom band
    band = z[-1, :]
    Z[-half_size:, half_size:-half_size] = band + band  - flipud( z[-half_size-1:-1, :] ) 
    # left band
    band = tile( z[:,0].reshape(-1,1), [1,half_size])
    Z[half_size:-half_size, :half_size] = band + band -  fliplr( z[:, 1:half_size+1] )
    # right band
    band = tile( z[:,-1].reshape(-1,1), [1,half_size] )
    Z[half_size:-half_size, -half_size:] =  band + band -  fliplr( z[:, -half_size-1:-1] ) 
    # central band
    Z[half_size:-half_size, half_size:-half_size] = z

    """
    # top left corner
    band = z[0,0]
    Z[:half_size,:half_size] = band - abs( flipud(fliplr(z[1:half_size+1,1:half_size+1]) ) - band )
    # bottom right corner
    band = z[-1,-1]
    Z[-half_size:,-half_size:] = band + abs( flipud(fliplr(z[-half_size-1:-1,-half_size-1:-1]) ) - band )
    """

    # top left corner
    band = z[0,0]
    Z[:half_size,:half_size] = band + band - flipud(fliplr(z[1:half_size+1,1:half_size+1]) ) 
    # bottom right corner
    band = z[-1,-1]
    Z[-half_size:,-half_size:] = band + band - flipud(fliplr(z[-half_size-1:-1,-half_size-1:-1]) )

    """    
    # top right corner
    band = Z[half_size,-half_size:]
    Z[:half_size,-half_size:] = band - abs( flipud(Z[half_size+1:2*half_size+1,-half_size:]) - band )
    # bottom left corner
    band = Z[-half_size:,half_size].reshape(-1,1)
    Z[-half_size:,:half_size] = band - abs( fliplr(Z[-half_size:, half_size+1:2*half_size+1]) - band )
    """

    # top right corner
    band = z[0,-1]
    Z[:half_size,-half_size:] = band + band -  flipud(fliplr(z[1:half_size+1,-half_size-1:-1]) )
    # bottom left corner
    band = z[-1,0]
    Z[-half_size:,:half_size] = band + band -  flipud(fliplr(z[-half_size-1:-1,1:half_size+1]) ) 
    


    #matshow(Z)

    #print 'test' , linalg.pinv(A)[0]

    # solve system and convolve
    if derivative == None:
        m = linalg.pinv(A)[0].reshape((window_size, -1))
        return signal.fftconvolve(Z, m, mode='valid')
    elif derivative == 'col':
        c = linalg.pinv(A)[1].reshape((window_size, -1))
        return signal.fftconvolve(Z, -c, mode='valid')
        #return signal.fftconvolve(Z, c/hx, mode='valid')
    elif derivative == 'row':
        r = linalg.pinv(A)[2].reshape((window_size, -1))
        return signal.fftconvolve(Z, -r, mode='valid')
        #return signal.fftconvolve(Z, r/hy, mode='valid')
    elif derivative == 'both':
        c = linalg.pinv(A)[1].reshape((window_size, -1))
        r = linalg.pinv(A)[2].reshape((window_size, -1))
        return signal.fftconvolve(Z, -r, mode='valid'), signal.fftconvolve(Z, -c, mode='valid')
        #return signal.fftconvolve(Z, r/hy, mode='valid'), signal.fftconvolve(Z, c/hx, mode='valid')
    elif derivative == 'all':
        m = linalg.pinv(A)[0].reshape((window_size, -1))
        c = linalg.pinv(A)[1].reshape((window_size, -1))
        r = linalg.pinv(A)[2].reshape((window_size, -1))
        return signal.fftconvolve(Z, m, mode='valid'), signal.fftconvolve(Z, -r, mode='valid'), signal.fftconvolve(Z, -c, mode='valid')



#====================================================================================
# use to smooth data - works with 2D data 
# >>>>>>>>>>>>> see test_smooth.py <<<<<<<<<<<<<<<<<<<
# smooth-2d.py
#
# purpose:  Smooth 2D data (horizontal/vertical sections)
# author:   Filipe P. A. Fernandes
# e-mail:   ocefpaf@gmail
# web:      http://ocefpaf.tiddlyspot.com/
# created:  27-Oct-2010
# modified: Wed 27 Oct 2010 11:43:39 AM EDT
#
# obs: http://www.scipy.org/Cookbook/SignalSmooth
#  check also cookb_signalsmooth.py


def gauss_kern(size, sizey=None):
    """ Returns a normalized 2D gauss kernel array for convolutions """
    size = int(size)
    if not sizey:
        sizey = size
    else:
        sizey = int(sizey)
    x, y = mgrid[-size:size+1, -sizey:sizey+1]
    g = exp(-(x**2/float(size)+y**2/float(sizey)))
    return g / g.sum()

def blur_image(im, n, ny=None) :
    """ blurs the image by convolving with a gaussian kernel of typical
        size n. The optional keyword argument ny allows for a different
        size in the y direction.
    """
    g = gauss_kern(n, sizey=ny)
    improc = signal.convolve(im,g, mode='valid')
    return(improc)


#==================================================================================

def coords(x,y):
   """ find the 3 first significant number of the coords """
   x = x/10.     # 10. because radius is never greater than 99 kpc
   px1 = int(x)
   px2 = int( (x-px1)*10.)  
   px3 = int( ((x-px1)*10.-px2)*10. + 0.5)  # +0.5 to get the right troncature
   px = px1*100 + px2*10 + px3

   y = y/10.
   py1 = int(y)
   py2 = int( (y-py1)*10.)
   py3 = int( ((y-py1)*10.-py2)*10. + 0.5)
   py = py1*100 + py2*10 + py3

   return px,py
  





#======================================================================


def gauss_fit(xi,yi,params):
    """ gaussian fit :  fit(x) = p[0]*exp( -(x-p[1])**2/(2*p[2]**2) ) 

        inputs :  x      : antecedants 
                  y      : image y(x)
                  params : array with initial parameters [ p[0],p[1],p[2] ] 
        output :  resfit : array with computed [ p[0],p[1],p[2] ]  """

    def residuals(p, y, x):
        err = y-peval(x,p)
        return err
    
    def peval(x, p):
        return p[0]*exp( -(x-p[1])**2/(2*p[2]**2) )    # gaussian


     
 
    def fits(x,y,A1,x01,sig1):

        pname = (['A','x0','sig'])

        p0 = array([A1,x01,sig1])

        p , plsq = leastsq(residuals, p0, args=(y, x), maxfev=600000)

        print "Final parameters for gaussian fit"
        for l in range(len(pname)):
             print '%s = %.5f' % (pname[l],p[l])
         
        return p


    resfit = fits(xi,yi,params[0],params[1],params[2])

    return resfit


#======================================================================


def gaussian(x,params):
    """ gaussian function :  y(x) = p[0]*exp( -(x-p[1])**2/(2*p[2]**2) ) 

        inputs :  x      : antecedants 
                  params : array with parameters [ p[0],p[1],p[2] ] 
        output :  y      : array with computed y(x) """



    y = p[0]*exp( -(x-p[1])**2/(2*p[2]**2) ) 

    return y







#====================================================================================
# doesn't work yet 


def residuals(p, y, x , cod):
    err = y-peval(x,p,cod)
    return err
    
def peval(x, p , cod):

    if (cod == 'g+g') :  # gaussian + gaussian
       return abs(p[0])*exp(-(x-p[1])**2/(2.*p[2]**2))+  abs(p[3])*exp(-(x-p[4])**2/(2.*p[5]**2)) 

    if (cod == 'exp') : 
       return p[0]*exp(-(x-p[1])**2/p[2])                  

    if (cod == 's') :   # sech
       return p[0]*2./(exp(-x/p[1])+exp(x/p[1]))          

    if (cod == 's2'):   # sech2 
       return p[0]*4./(exp(-x/p[1])+exp(x/p[1]))**2    

    if (cod == 'log(exp)') : # log10(exp) 
       return p[0]-x/p[1]                 

    if (cod == 'log(s2)') : # log10(sech2) 
       return p[0]-2.*log10(exp(-x/p[1])+exp(x/p[1]))                 

  
     
 
def fits(x,y,pname,p0,cod):
   """ inputs : x , y , pname , p0 , code

       fits with : gaussian+gaussian, code = 'g+g'
                   exponential      , code = 'exp'
                   sech             , code = 's'
                   sech**2          , code = 's2' 
                   log10(exp)       , code = 'log(exp)' 
                   log10(sech**2)   , code = 'log(s2)'  """
   #A1=0.1
   #A2=0.1
   #n1=1.
   #n2=1.
  
   #pname = (['A1','x01','sig1','A2','x02','sig2'])

   #p0 = array([A1,x01,sig1, A2,x02,sig2 ])


   p , plsq = leastsq(residuals, p0, args=(y, x), maxfev=2000000)

   print
   print "Final parameters for fit"
   for i in range(len(pname)):
       print '%s = %.5f' % (pname[i],p[i])

   return p



#============================================================================================

def plot_pcolor(xbins,ybins,bins,**kwargs):
    """ pcolor plot z particules associated to values x and y, are put in xbins and ybins"""



    opt = { 'power'     : 1.0        ,
            'vmin'      : None       ,  
            'vmax'      : None       ,
            'bar_label' : 'Number'   ,
            'fsize'     : 20         ,
            'cmap'      : 'jet'      ,
            'scale'     : 'lin'      ,
            'bar_orient': 'vertical' ,
            'position'  : 'center'   ,
            'invalid'   : nan 
           }

    opt.update(kwargs)
    print opt

    if opt['invalid'] != nan :
        bins[ bins == opt['invalid'] ] = nan


    bins = ma.masked_invalid(bins)
    cmap = get_cmap(opt['cmap'])
    cmap.set_bad(color = 'w', alpha = 1.)


    def func(x,y):
      return bins[x,y]

    I=range(len(xbins))
    J=range(len(ybins))
    X,Y = meshgrid(I, J)



    Z = func(X, Y)**(opt['power'])


    if opt['vmin'] == None :
        vmin = Z.min()
    else :
        vmin = opt['vmin']
        
        
    if opt['vmax'] == None :
        vmax = Z.max()
    else :
        vmax = opt['vmax']
        
        
    if opt['position'] == 'center' :
        x_offset = 0.5 * (xbins[1]-xbins[0])
        y_offset = 0.5 * (ybins[1]-ybins[0])
    else :
        x_offset = 0.
        y_offset = 0.



    if opt['scale'] == 'log' :

        PP = pcolor(xbins[X] - x_offset , ybins[Y] - x_offset , Z, norm=LogNorm(vmin=1.0, vmax=vmax) , cmap = cmap )#, vmin=0.0 ,v)

    else :

        #PP = pcolor(xbins[X] - x_offset , ybins[Y] - x_offset , Z, vmin=vmin, vmax=vmax , cmap = cmap )#, vmin=0.0 ,v)
        PP = pcolor(xbins[X] , ybins[Y] , Z, vmin=vmin, vmax=vmax , cmap = cmap )#, vmin=0.0 ,v)

   

    cb = colorbar(PP, orientation=opt['bar_orient'])
    cb.set_label(r'$\bf '+opt['bar_label']+' $',fontsize=opt['fsize'])
  

    xlim(xbins[0],xbins[-1])
    ylim(ybins[0],ybins[-1])
   
   
   
   
#============================================================================================

def plot_contour(xbins,ybins,bins,**kwargs):
    """ pcolor plot z particules associated to values x and y, are put in xbins and ybins"""



    opt = { 'power'     : 1.0        ,
            'vmin'      : None       ,  
            'vmax'      : None       ,
            'bar_label' : 'Number'   ,
            'fsize'     : 20         ,
            'cmap'      : 'jet'      ,
            'scale'     : 'lin'      ,
            'bar_orient': 'vertical' ,
            'position'  : 'center'   ,
            'invalid'   : nan        ,
            'levels'    : None
           }

    opt.update(kwargs)
    print opt

    if opt['invalid'] != nan :
        bins[ bins == opt['invalid'] ] = nan


    bins = ma.masked_invalid(bins)
    cmap = get_cmap(opt['cmap'])
    cmap.set_bad(color = 'w', alpha = 1.)


    def func(x,y):
      return bins[x,y]

    I=range(len(xbins))
    J=range(len(ybins))
    X,Y = meshgrid(I, J)



    Z = func(X, Y)**(opt['power'])


    if opt['vmin'] == None :
        vmin = Z.min()
    else :
        vmin = opt['vmin']
        
        
    if opt['vmax'] == None :
        vmax = Z.max()
    else :
        vmax = opt['vmax']
        
        
    if opt['position'] == 'center' :
        x_offset = 0.5 * (xbins[1]-xbins[0])
        y_offset = 0.5 * (ybins[1]-ybins[0])
    else :
        x_offset = 0.
        y_offset = 0.



    if opt['scale'] == 'log' :

        PP = contour(xbins[X] - x_offset , ybins[Y] - x_offset , Z, norm=LogNorm(vmin=1.0, vmax=vmax) , cmap = cmap ,levels=opt['levels'] )#, vmin=0.0 ,v)

    else :

        PP = contour(xbins[X] - x_offset , ybins[Y] - x_offset , Z, vmin=vmin, vmax=vmax , cmap = cmap , levels=opt['levels'] )#, vmin=0.0 ,v)

   
    if opt['bar_label'] == None :
        pass
    else :
        cb = colorbar(PP, orientation=opt['bar_orient'])
        cb.set_label(r'$\bf '+opt['bar_label']+' $',fontsize=opt['fsize'])
  

    xlim(xbins[0],xbins[-1])
    ylim(ybins[0],ybins[-1])
   

#================================================================================  

def maxval(l):

   mx = l[0]
   for i in range(len(l)):
       if l[i] > mx :
          mx = l[i]
          

   return mx


#======================================================================================

#========================================================================================


def fourier_m2(rbins,x,y,nsr,pr):
    """ ==========================================================================

        compute m=2 fourier mode in cylindrical coordinates, 

        2D optimized binning
        inputs :  rbins    : the radial bins 
                  x        : x coordinate of particles to radialy bin
                  y        : y coordinate of particles to radialy bin
                  nsr      : number of signif nums to extract from each radius[i]
                  pr       : max( p | 10**p < max(radius[:]) )
   
        output :  nbins    : array of integer or reals , number of particules in each bin

        x and y have to be of same dimension

        work with the function 'signif_nums' 
        ===========================================================================
     """



    rad = sqrt( x**2 + y**2 )


    """

    in each radial bin :

      N = number of particles

      coeff a = sum( cos(2*theta[i]), i=1..N )
      coeff b = sum( sin(2*theta[i]), i=1..N )

    coeff c = sqrt( a**2 + b**2 )  -->  stength of mode m=2






    cos(2*a) = 1 - 2*sin**2(a)
             = 1 - 2* y**2 / ( x**2 + y**2 )

    sin(2*a) = 2*sin(a)*cos(a)
             = 2* y * x / ( x**2 + y**2 )

    """


    rminbins = min(rbins[:])  ; rmaxbins =max(rbins[:])
    rwidth = rmaxbins - rminbins
    rintwidth = signif_nums ( rwidth , nsr , pr)


    dnr = rintwidth/(len(rbins)-1)



    nbins      = zeros ( (len(rbins)) )
    c          = zeros ( (len(rbins)) )

    nbins_temp = zeros ( (len(rbins),2) )

    for i in range(len(rad)):
        rb = signif_nums(rad[i] - rminbins, nsr , pr)

        rb=rb/dnr

        if (rb>=0) and (rb<=len(rbins)-1):

           a = 1. - 2. * y[i]**2 / ( x[i]**2 + y[i]**2 )
           b = 2. * y[i] * x[i]  / ( x[i]**2 + y[i]**2 )


           nbins_temp[rb,0] = nbins_temp[rb,0] + a
           nbins_temp[rb,1] = nbins_temp[rb,1] + b

           nbins[rb] = nbins[rb] + 1.



    c[:] =  sqrt( nbins_temp[:,0]**2 + nbins_temp[:,1]**2 )



    return c,nbins





#========================================================================================


def fourier_m(m,x,y):
    """ ==========================================================================

        compute fourier mode m in cylindrical coordinates, 

        
        inputs :  m        : the desired fourier mode
                  x        : x coordinate of data points
                  y        : y coordinate of data points

   
        output :  c        : strength of mode m

        x and y have to be of same dimension

 
        ===========================================================================
     
        in each radial bin :

        N = number of particles

        coeff a = 2*sum( y[i]*cos(m*x[i]), i=1..N )
        coeff b = 2*sum( y[i]*sin(m*x[i]), i=1..N )

        coeff c = sqrt( a**2 + b**2 )  -->  stength of mode m
    """


    am = 0.
    bm = 0.

    for i in xrange(len(x)):
    
        if ( type(y).__name__ == 'float' ):

            am = am + cos( m * x[i] )
            bm = bm + sin( m * x[i] )

        else :

            am = am + y[i]*cos( m * x[i] )
            bm = bm + y[i]*sin( m * x[i] )

    am = 2.*am/float(len(x))
    bm = 2.*bm/float(len(x))

    c =  sqrt( am**2 + bm**2 )



    return c,am,bm





#=====================================================================================

def lin_interp_2D(xbins,ybins,values,x,y):
    """ do a 2D linear interpolation 
        We know the values on a grid (xbins,ybins);
        we get the interpolated value res at point (x,y) 
     """


    bins = binning_2D(xbins,x,7,1,ybins,y,7,1,1.)

    
    res = values[ where( bins != 0. ) ]

    return res

























