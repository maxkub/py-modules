from pylab   import *

#================================================================================
#================================================================================

def signif_nums(x,n,to, **kwargs):

   """ find the n first significant number of x 
       inputs : x : number to work with
                n : number of signif nums to extract from x
                to: max(10**(to)) such as 10**(to) < x
       outputs: out: integer made of the n signif nums """

   opt = { 'print'   : False  # outputs to test the results   
           }
           

   opt.update(kwargs)

   pwo = []

   xtemp = x/10.**(to)    
   pwo.append( int(xtemp) )

   for i in xrange(n-1):
       xtemp = (xtemp - pwo[i]) * 10.
       #pwo[i+1] = int(xtemp)
       pwo.append( int(xtemp) )

  
   xtemp = (xtemp - pwo[n-2]) * 10. + 0.5
   #pwo[n-1] = int(xtemp)
   pwo.append( int(xtemp) )


   out = 0
   for i in xrange(n):
       out = out + pwo[i]*10**(n-1-i)
       
   if opt['print'] :
      print 'signif nums' , out

   return out
   
#======================================================================


def binning_1D(xbins,x,ns,p,z ,  **kwargs):
    """ 1D optimized binning
        inputs :  xbins  : the bins 
                  x      : values to put in bins
                  ns     : number of signif nums to extract from each x[i]
                  p      : max( p | 10**p < max(x[:]) )
                  z      : value associated with each x[i], to bin (can be 1 to count) 
        output :  nbins  : array of integer , number of particules in each bin

        work with the function 'signif_nums' """
        
        
        
    opt = { 'print'   : False  # outputs to test the results   
           }
           

    opt.update(kwargs)

    minbins = min(xbins[:])  ; maxbins =max(xbins[:])
    #xwidth = abs(minbins) + abs(maxbins)
    xwidth = maxbins - minbins
    intwidth = signif_nums ( xwidth , ns , p)

    dn = intwidth/(len(xbins)-1)

    nbins = zeros (len(xbins))


    if opt['print']:
       print 'xwidth' , xwidth
       print 'intwidth' , intwidth
       print 'dn' , dn


    if ( type(z).__name__ == 'int' ) or ( type(z).__name__ == 'float' ):
      
       for i in range(len(x)):
           a = signif_nums(x[i] - minbins, ns , p)

           a=a/dn

           if (a>=0) and (a<=len(xbins)-1):
              nbins[a] = nbins[a]+z

    else:

       for i in xrange(len(x)):
           a = signif_nums(x[i] - minbins, ns , p)

           a=a/dn

           if (a>=0) and (a<=len(xbins)-1):
              nbins[a] = nbins[a]+z[i]


    return nbins
#========================================================================================


def binning_2D(xbins,x,nsx,px,ybins,y,nsy,py,z):
    """ 2D optimized binning
        inputs :  xbins    : the x bins 
                  x        : values to put in xbins
                  nsx      : number of signif nums to extract from each x[i]
                  px       : max( p | 10**p < max(x[:]) )
                  ybins    : the y bins 
                  y        : values to put in ybins
                  nsy      : number of signif nums to extract from each y[i]
                  py       : max( p | 10**p < max(y[:]) )
                  z        : value associated with each particules, to bin (can be 1 to count) 
        output :  nbins    : array of integer or reals , number of particules in each bin

        x and y have to be of same dimension

        work with the function 'signif_nums' """



    xminbins = min(xbins[:])  ; xmaxbins =max(xbins[:])
    #xwidth = abs(xminbins) + abs(xmaxbins)
    xwidth = xmaxbins - xminbins
    xintwidth = signif_nums ( xwidth , nsx , px)

    yminbins = min(ybins[:])  ; ymaxbins =max(ybins[:])
    #ywidth = abs(yminbins) + abs(ymaxbins)
    ywidth = ymaxbins - yminbins
    yintwidth = signif_nums ( ywidth , nsy , py)


    dnx = xintwidth/(len(xbins)-1)
    dny = yintwidth/(len(ybins)-1)

    nbins = zeros ( (len(xbins),len(ybins)) )


    if ( type(z).__name__ == 'int' ) or ( type(z).__name__ == 'float' ):

        for i in xrange(len(x)):
            a = signif_nums(x[i] - xminbins, nsx , px)
            b = signif_nums(y[i] - yminbins, nsy , py)

            #print 'a-b',a,b

            a=a/dnx
            b=b/dny
            
            #print 'dnx - dny', dnx,dny, 'a,b',a,b

            if (a>=0) and (a<=len(xbins)-1) and (b>=0) and (b<=len(ybins)-1):
               nbins[a,b] = nbins[a,b]+z

    else:

        for i in xrange(len(x)):
            a = signif_nums(x[i] - xminbins, nsx , px)
            b = signif_nums(y[i] - yminbins, nsy , py)

            a=a/dnx
            b=b/dny

            if (a>=0) and (a<=len(xbins)-1) and (b>=0) and (b<=len(ybins)-1):
               nbins[a,b] = nbins[a,b]+z[i]


    return nbins



#========================================================================================


def binning_3D(xbins,x,nsx,px,ybins,y,nsy,py,zbins,z,nsz,pz,cpt):
    """ 3D optimized binning
        inputs :  xbins    : the x bins 
                  x        : values to put in xbins
                  nsx      : number of signif nums to extract from each x[i]
                  px       : max( p | 10**p < max(x[:]) )
                  ybins    : the y bins 
                  y        : values to put in ybins
                  nsy      : number of signif nums to extract from each y[i]
                  py       : max( p | 10**p < max(y[:]) )
                  zbins    : the z bins 
                  z        : values to put in zbins
                  nsz      : number of signif nums to extract from each z[i]
                  pz       : max( p | 10**p < max(z[:]) )
                  cpt        : value associated with each particules, to bin (can be 1 to count) 
        output :  nbins    : array of integer or reals , number of particules in each bin

        x and y have to be of same dimension

        work with the function 'signif_nums' """



    xminbins = min(xbins[:])  ; xmaxbins =max(xbins[:])
    xwidth = xmaxbins - xminbins
    xintwidth = signif_nums ( xwidth , nsx , px)

    yminbins = min(ybins[:])  ; ymaxbins =max(ybins[:])
    ywidth = ymaxbins - yminbins
    yintwidth = signif_nums ( ywidth , nsy , py)

    zminbins = min(zbins[:])  ; zmaxbins =max(zbins[:])
    zwidth = zmaxbins - zminbins
    zintwidth = signif_nums ( zwidth , nsz , pz)


    dnx = xintwidth/(len(xbins)-1)
    dny = yintwidth/(len(ybins)-1)
    dnz = zintwidth/(len(zbins)-1)

    nbins = zeros ( (len(xbins),len(ybins),len(zbins)) )


    if ( type(cpt).__name__ == 'int' ) or ( type(cpt).__name__ == 'float' ):

        for i in xrange(len(x)):
            a = signif_nums(x[i] - xminbins, nsx , px)
            b = signif_nums(y[i] - yminbins, nsy , py)
            c = signif_nums(z[i] - zminbins, nsz , pz)

            a=a/dnx
            b=b/dny
            c=c/dnz

            if (a>=0) and (a<=len(xbins)-1) and (b>=0) and (b<=len(ybins)-1)  and (c>=0) and (c<=len(zbins)-1):
               nbins[a,b,c] = nbins[a,b,c]+cpt

    else:

        for i in xrange(len(x)):
            a = signif_nums(x[i] - xminbins, nsx , px)
            b = signif_nums(y[i] - yminbins, nsy , py)
            c = signif_nums(z[i] - zminbins, nsz , pz)

            a=a/dnx
            b=b/dny
            c=c/dnz

            if (a>=0) and (a<=len(xbins)-1) and (b>=0) and (b<=len(ybins)-1) and (c>=0) and (c<=len(zbins)-1):
               nbins[a,b,c] = nbins[a,b,c]+cpt[i]


    return nbins




#========================================================================================


def binning_4D(xbins,x,nsx,px,ybins,y,nsy,py,zbins,z,nsz,pz,wbins,w,nsw,pw,cpt):
    """ 4D optimized binning
        inputs :  xbins    : the x bins 
                  x        : values to put in xbins
                  nsx      : number of signif nums to extract from each x[i]
                  px       : max( p | 10**p < max(x[:]) )
                  ybins    : the y bins 
                  y        : values to put in ybins
                  nsy      : number of signif nums to extract from each y[i]
                  py       : max( p | 10**p < max(y[:]) )
                  zbins    : the z bins 
                  z        : values to put in zbins
                  nsz      : number of signif nums to extract from each z[i]
                  pz       : max( p | 10**p < max(z[:]) )
                  wbins    : the z bins 
                  w        : values to put in zbins
                  nsw      : number of signif nums to extract from each z[i]
                  pw       : max( p | 10**p < max(z[:]) )
                  cpt        : value associated with each particules, to bin (can be 1 to count) 
        output :  nbins    : array of integer or reals , number of particules in each bin

        x and y have to be of same dimension

        work with the function 'signif_nums' """



    xminbins = min(xbins[:])  ; xmaxbins =max(xbins[:])
    xwidth = xmaxbins - xminbins
    xintwidth = signif_nums ( xwidth , nsx , px)

    yminbins = min(ybins[:])  ; ymaxbins =max(ybins[:])
    ywidth = ymaxbins - yminbins
    yintwidth = signif_nums ( ywidth , nsy , py)

    zminbins = min(zbins[:])  ; zmaxbins =max(zbins[:])
    zwidth = zmaxbins - zminbins
    zintwidth = signif_nums ( zwidth , nsz , pz)

    wminbins = min(wbins[:])  ; wmaxbins =max(wbins[:])
    wwidth = wmaxbins - wminbins
    wintwidth = signif_nums ( wwidth , nsw , pw)


    dnx = xintwidth/(len(xbins)-1)
    dny = yintwidth/(len(ybins)-1)
    dnz = zintwidth/(len(zbins)-1)
    dnw = wintwidth/(len(wbins)-1)
    
    print dnx,dny,dnz,dnw

    nbins = zeros ( (len(xbins),len(ybins),len(zbins),len(wbins)) )


    if ( type(cpt).__name__ == 'int' ) or ( type(cpt).__name__ == 'float' ):

        for i in xrange(len(x)):
            a = signif_nums(x[i] - xminbins, nsx , px)
            b = signif_nums(y[i] - yminbins, nsy , py)
            c = signif_nums(z[i] - zminbins, nsz , pz)
            d = signif_nums(w[i] - wminbins, nsw , pw)

            a=a/dnx
            b=b/dny
            c=c/dnz
            d=d/dnw

            if (a>=0) and (a<=len(xbins)-1) and (b>=0) and (b<=len(ybins)-1)  and (c>=0) and (c<=len(zbins)-1) and (d>=0) and (d<=len(wbins)-1):
               nbins[a,b,c,d] = nbins[a,b,c,d]+cpt

    else:

        for i in xrange(len(x)):
            a = signif_nums(x[i] - xminbins, nsx , px)
            b = signif_nums(y[i] - yminbins, nsy , py)
            c = signif_nums(z[i] - zminbins, nsz , pz)
            d = signif_nums(w[i] - wminbins, nsw , pw)
            
            print 'a',a,b,c,d

            a=a/dnx
            b=b/dny
            c=c/dnz
            d=d/dnw

            if (a>=0) and (a<=len(xbins)-1) and (b>=0) and (b<=len(ybins)-1) and (c>=0) and (c<=len(zbins)-1) and (d>=0) and (d<=len(wbins)-1):
               nbins[a,b,c,d] = nbins[a,b,c,d]+cpt[i]


    return nbins


#========================================================================================


def binning_ND( N , **kwargs ):
    """ ND optimized binning
        inputs :  N        : number of parameter to bin
        
        in opt['ins'] xbins: the x bins 
                  x        : values to put in xbins
                  nsx      : number of signif nums to extract from each x[i]
                  px       : max( p | 10**p < max(x[:]) )
                  ybins    : the y bins 
                  y        : values to put in ybins
                  nsy      : number of signif nums to extract from each y[i]
                  py       : max( p | 10**p < max(y[:]) )
                  zbins    : the z bins 
                  z        : values to put in zbins
                  nsz      : number of signif nums to extract from each z[i]
                  pz       : max( p | 10**p < max(z[:]) )
                  wbins    : the z bins 
                  w        : values to put in zbins
                  nsw      : number of signif nums to extract from each z[i]
                  pw       : max( p | 10**p < max(z[:]) )
                  cpt      : value associated with each particules, to bin (can be 1 to count) 
                  
        output :  nbins    : array of integer or reals , number of particules in each bin

        x and y have to be of same dimension

        work with the function 'signif_nums' """

    opt={ 'ins'  :  []  }

    opt.update( kwargs)
    
    
    minbins  = zeros( N )
    maxbins  = zeros( N )
    width    = zeros( N )
    intwidth = zeros( N )
    dn       = zeros( N )
    num      = zeros( N )
    
    for i in range(N) :
        minbins[i]  = min(opt['ins'][4*i])  ; maxbins[i] =max(opt['ins'][4*i])
        width[i]    = maxbins[i] - minbins[i]
        intwidth[i] = signif_nums( width[i] , opt['ins'][2+4*i] , opt['ins'][3+4*i] )

        dn[i] = int( intwidth[i]/(len(opt['ins'][4*i])-1) )
        print i, dn[i]
 

    

    nbins = zeros ( list( len(opt['ins'][4*i]) for i in range(N) ) )#, dtype = 'i4')



    if ( type( opt['ins'][-1] ).__name__ == 'int' ) or ( type( opt['ins'][-1] ).__name__ == 'float' ):

        for i in xrange(len( opt['ins'][1] )):
           
            for j in range(N) :
               
                try :
                    num[j] = signif_nums( opt['ins'][1+4*j][i] - minbins[j] , opt['ins'][2+4*j] , opt['ins'][3+4*j] )
                except IndexError:
                    print 'IndexError',i,j


            num = [ [int(num[j]/dn[j])] for j in xrange(N) ]
 
            #print num
            #print nbins[num]
            
            try :
                nbins[num] = nbins[num]+ opt['ins'][-1]
            except IndexError :
                print 'IndexError', i, num
                
            #print nbins[num]


    else:

        for i in xrange(len(opt['ins'][1])):
           
           
            for j in range(N) :
                num[j] = signif_nums( opt['ins'][1+4*j][i] - minbins[j] , opt['ins'][2+4*j] , opt['ins'][3+4*j] )
                

            num = [ [num[j]/dn[j]] for j in range(N) ]
 
            try :
                nbins[num] = nbins[num]+ opt['ins'][-1][i]
            except IndexError :
                pass


    return nbins



#====================================================================================


def bin_index(xbins,ns,p,value):
    """ give the index associated to a value, in the bins
        inputs :  xbins  : the bins 
                  ns     : number of signif nums to extract from each x[i]
                  p      : max( p | 10**p < max(x[:]) )
                  value  : value whom we are looking for the corresponding index in the bins 
        output :  index  : index corresponding to the value (integer)

        work with the function 'signif_nums' """



    minbins = min(xbins[:])  ; maxbins =max(xbins[:])
    #xwidth = abs(minbins) + abs(maxbins)
    xwidth = maxbins - minbins
    intwidth = signif_nums ( xwidth , ns , p)

    dn = intwidth/(len(xbins)-1)

    nbins = zeros (len(xbins))

      
   
    a = signif_nums(value - minbins, ns , p)
    a=a/dn

    if (a>=0) and (a<=len(xbins)-1):
        index = a
        return index





    return nbins



