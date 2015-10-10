



"""
count_angm()            # sum of the DeltaL for migrating stars and all stars
plot_tot_angm()         # histogram with final and initial angular momentum for all stars
plot_migr_angm()        # histogram with final and initial angular momentum for selected stars
plot_delta_angm()       # histogram with delta angular momentum 
plot_delta_r()          # histogram with delta r 
plot_delta_gc()         # test of guiding center determination
plot_agemig()           # histogram with the age at wich the stars started migrating

plot_rveli()             # histogram with initial radial velocities
plot_tveli()             # histogram with initial tangential velocities

plot_dvelt_migr()        #histogram with delta tangential velocities of migrating particules
plot_dvelt_tot()        #histogram with delta tangential velocities of all particules


plot_hist_ritot()       # histogram with R initial for all stars

plot_dl_linit_migr()    # plot of delta_L vs Linit  for migrating stars


plot_dl_linit_tot()     # plot of delta_L vs Linit  for all stars

plot_dl_linit_tot_tvel()     # plot of delta_L vs Linit vs tangential velocity for all stars
plot_dl_linit_tot_rvel()     # plot of delta_L vs Linit vs radial velocity for all stars

plot_dl_linit_tot_disptvel() # plot of delta_L vs Linit vs tangential velocity dispersion for all stars
plot_dl_linit_tot_disprvel() # plot of delta_L vs Linit vs radial velocity dispersion for all stars

plot_dl_linit_tot_rf()        # plot of delta_L vs Linit vs R final for all stars
plot_dl_linit_tot_ri()        # plot of delta_L vs Linit vs R initial for all stars
plot_dl_linit_tot_dr()       # plot of delta_L vs Linit vs delta R for all stars
plot_dl_linit_tot_dispdr()   # plot of delta_L vs Linit vs delta R dispersion for all stars


plot_dgc_rinit_migr()   # plot of delta_guiding center vs Rinit  for migrating stars


plot_dgc_rinit_tot()   # plot of delta_guiding center vs Rinit  for all stars


plot_dl_lfin_migr()     # plot of delta_L vs Lfin  for migrating stars

plot_dl_lfin_tot()      # plot of delta_L vs Lfin  for all stars
plot_dl_Rfin_tot()      # plot of delta_L vs Rfin  for all stars

plot_dl_mt()            # plot of delta_L vs time of migration beginning for selected stars

plot_dr_mt()            # plot of delta_R vs time of migration beginning for selected stars
plot_dr_rinit_migr()    # plot of delta_R vs Rinit of migration beginning for selected stars


plot_dr_rinit_tot()     # plot of delta_R vs Rinit for all stars

plot_dr_rinit_tot_dLtot()


plot_posmig(t1,t2)      # point plot x,y of when the stars have started migrating


plot_tang_veli_r()       # plot of initial tangential vel vs radius  for migrating stars

plot_rad_veli_r()        # plot of initial radial vel vs radius  for migrating stars


plot_tang_velitot_r()    # plot of initial tangential vel vs radius  for all stars

plot_rad_velitot_r()     # plot of initial radial vel vs radius  for all stars


plot_tang_velf_r()      # plot of final tangential vel vs radius  for migrating stars
plot_rad_velf_r()       # plot of final radial vel vs radius  for migrating stars

plot_tang_velftot_r()   # plot of final tangential vel vs radius  for all stars
plot_rad_velftot_r()    # plot of final radial vel vs radius  for all stars


plot_num_rad_tot()      # plot of init and final number vs radius for all stars


plot_i_f_pos_tot()      # plot initial and final position (2D) for all stars
plot_i_f_pos_migr()     # plot initial and final position (2D) for migrating stars




plot_ri_rf_tot()       # plot initial vs final radius of particles for all particles
plot_ri_rf_migr()      # plot initial vs final radius of particles for migrating particles
plot_ri_rf_nomigr()    # plot initial vs final radius of particles for non migrating particles



plot_angmi_angmf_tot()       # plot initial vs final angular momentum of particles for all particles
plot_angmi_angmf_migr()      # plot initial vs final angular momentum of particles for migrating particles
plot_angmi_angmf_nomigr()      # plot initial vs final angular momentum of particles for non migrating particles

plot_angmi_angmf_nomigr_fold()      # plot initial vs final angular momentum of particles for non migrating particles and fold it with gaussian
"""




#==============================================================================

#### find the line of idd 

def find_l(idd):
   """ find the line of the star number idd"""
   #idd

   cpt=132109
   find=0
   while (find ==0):
      if (idn[cpt] == idd):
         print str(idd), 'l=', cpt
         find=1
      cpt = cpt +1
   return
       

#==============================================================================


def plot_hist_ritot():


   figure()

   xlabel('radius (kpc)')

   title('Rinit distrib')

   hist(ritot,50)


   print 'name', plot_hist_ritot.__name__+'.png'
  
   print inspect.stack()[2]

   if save==1:
      #savefig(savepath+'hist_Ritot.png')
      savefig(savepath+'hist_Ritot.png')




#=============================================================================

def count_dir():
   """ count the number of stars migrating inward and outward"""
   plus = 0
   minus = 0
   for i in range(len(rf)):
      if (rf[i] - ri[i] > 0.):
         plus = plus + 1
      else:
         minus = minus + 1

   print 'minus=',minus
   print 'plus=',plus
 

#==============================================================================

def count_angm():
   """ compute the total z-angular momentum variation """

   angvar = 0.
   for i in range(len(delta_angm)):
      angvar = angvar + delta_angm[i]

   angvart = 0.
   for i in range(len(positot)):
      angvart = angvart + delta_angmtot[i]


   print 'delta angular momentum all migr=', angvar
   print 'delta angular momentum all =', angvart
   print
   


#===========================================================================================
def plot_delta_gc():
   """ histogram with delta guiding center """
   figure()

   xlabel(r'$\Delta r$',fontsize=16)
   ylabel(r'$normed\ distribution$',fontsize=16)

   title('delta guiding center and delta radius')

   #hist(delta_gc,200,range=(-8,8),label='guid. cent.')
   #hist(delta_r,200,histtype='step',color='r',linewidth=2,label='radius')


   dr_bins=[]
   for i in range(nDrad+1):
       dr_bins.append( Dradmin + (Dradmax - Dradmin)*float(i)/float(nDrad) ) 


   bins1 = binning_1D(dr_bins,delta_gc,5,1,1.)

   bins2 = binning_1D(dr_bins,delta_r,5,1,1.)


   max1 =max(bins1[:])
   bins1 =bins1/max1
   bins2 =bins2/max1


   plot(dr_bins,bins1,linewidth=2,label='guid. cent.')
   plot(dr_bins,bins2,linewidth=2,label='radius')

   

   ####folds with gaussian

   bins3 =zeros(len(bins2))

   dr =(Dradmax - Dradmin)/float(nDrad)
   sfold = 1.

   for j in range(len(bins3)):
       for i in range(len(bins3)):
           bins3[j]=bins3[j]+bins1[i]/(sfold*sqrt(2.*pi))*exp(-(dr_bins[i]-dr_bins[j])**2/(2.*sfold**2))*dr

   #bins3=bins3/dr  
   #max_val=max(bins3)
   #bins3=bins3/max_val

   #bins3 = bins3/max1

   plot(dr_bins,bins3,'k--',linewidth=1,label='folded')



  


   def residuals(p, y, x):
       err = y-peval(x,p)
       return err
    
   def peval(x, p):
       return abs(p[0])*exp(-(x-p[1])**2/(2.*p[2]**2))+  abs(p[3])*exp(-(x-p[4])**2/(2.*p[5]**2))  
       #return p[0]*exp(-(x-p[1])**2/p[2])+  0.000009*exp(-(x-12.2)**2/3.)                     # gaussian + gaussian


     
 
   def fits(x,y,A1,x01,sig1,A2,x02,sig2):
      #A1=0.1
      #A2=0.1
      #n1=1.
      #n2=1.
  
      pname = (['A1','x01','sig1','A2','x02','sig2'])

      p0 = array([A1,x01,sig1, A2,x02,sig2 ])

      #pname = (['A1','x01','n1'])

      #p0 = array([A1 ,x01,n1])


      p , plsq = leastsq(residuals, p0, args=(y, x), maxfev=2000000)

      print
      print "Final parameters for fit"
      for i in range(len(pname)):
          print '%s = %.5f' % (pname[i],p[i])

      return p



   res2 = fits(dr_bins,bins2,0.1,-1.,1.,0.1,1.,1.)
   plot(dr_bins,peval(dr_bins,res2),'r--',linewidth=1,label='fit rad')


   res3 = fits(dr_bins,bins3,0.1,-1.,1.,0.1,1.,1.)
   plot(dr_bins,peval(dr_bins,res3),'r--',linewidth=1,label='fit fold')




   grid('on')
   legend(loc='best')


   if save==1:
      savefig(savepath+'Dgc.png')



#===========================================================================================
def plot_tveli():
   """ histogram with initial tangential velocities """
   figure()

   xlabel('velocities')

   title('tangential velocities distrib')

   hist(tveli,50,label='migr')
  
   legend(loc='best')

   if save==1:
      savefig(savepath+'Tvel.png')


#===========================================================================================
def plot_dvelt_migr():
   """ histogram with delta tangential velocities of migrating particules"""
   figure()

   xlabel('velocities')

   title('delta tangential velocities / migr stars')

   hist(delta_vt,50,label='migr')
  
   legend(loc='best')

   if save==1:
      savefig(savepath+'DTvel-migr.png')

#===========================================================================================
def plot_dvelt_tot():
   """ histogram with delta tangential velocities of all particules"""
   figure()

   xlabel('velocities')

   title('delta tangential velocities / all stars')

   hist(delta_vttot,50,label='all')
  
   legend(loc='best')

   if save==1:
      savefig(savepath+'DTvel-all.png')

#===========================================================================================
def plot_rveli():
   """ histogram with initial radial velocities """
   figure()

   xlabel('velocities')

   title('radial velocities distrib')

   hist(rveli,50,label='migr')
  
   legend(loc='best')

   if save==1:
      savefig(savepath+'init-Rvel.png')

#===========================================================================================
def plot_delta_angm():
   """ histogram with delta angular momentum """
   figure()
   title(r'$\bf \Delta angular\ momentum$',fontsize=18)

   xlabel(r'$\Delta ang.\ momentum$',fontsize=16)
   ylabel(r'$normed \ distribution$',fontsize=16)


   hist(delta_angmtot,50,range=(DLmin,DLmax),normed='True',histtype='step',color='r',linewidth=2,label='tot')
   hist(delta_angm,50,range=(DLmin,DLmax),normed='true',alpha=0.5,label='migr')

 
   legend(loc='best')

   if save==1:
      savefig(savepath+'Dangm.png')


#===========================================================================================
def plot_tot_angm():
   """ histogram with final and initial angular momentum for all stars"""
   figure()
   title(r'$\bf tot\ angular\ momentum\ distrib$',fontsize=18)

   xlabel(r'$ ang.\ momentum$',fontsize=16)
   ylabel(r'$normed \ distribution$',fontsize=16)


   #hist(delta_angmtot,200,range=(-1000,1000),normed='True',histtype='step',color='r',linewidth=2,label='tot')
   #hist(delta_angm,200,range=(-1000,1000),normed='true',alpha=0.5,label='migr')

   
   #nb = 50
   angm_bins=zeros(nL+1)
   for i in range(nL+1):
       angm_bins[i] = Lmin + (Lmax - Lmin)*float(i)/float(nL)


   bins1 = binning_1D(angm_bins,angmitot,7,3,1)

   bins2 = binning_1D(angm_bins,angmftot,7,3,1)


   max1 =max(bins1[:])
   bins1 =bins1/max1
   bins2 =bins2/max1

   plot(angm_bins,bins1,linewidth=2,label='init tot')
   plot(angm_bins,bins2,linewidth=2,label='fin tot')
   #plot(dangm_bins,bins,label='tot-migr')
   

   legend(loc='best')
   grid('on')

   if save==1:
      savefig(savepath+'i-f-angm-all.png')
#===========================================================================================
def plot_migr_angm():
   """ histogram with final and initial angular momentum for selected stars"""
   figure()
   title(r'$\bf migr\ angular\ momentum\ distrib$',fontsize=18)

   xlabel(r'$ ang.\ momentum$',fontsize=16)
   ylabel(r'$normed \ distribution$',fontsize=16)


   #hist(delta_angmtot,200,range=(-1000,1000),normed='True',histtype='step',color='r',linewidth=2,label='tot')
   #hist(delta_angm,200,range=(-1000,1000),normed='true',alpha=0.5,label='migr')

 
   #nb = 50
   angm_bins=zeros(nL+1)
   for i in range(nL+1):
       angm_bins[i] = Lmin + (Lmax - Lmin)*float(i)/float(nL)


   bins1 = binning_1D(angm_bins,angminit,5,1,1)

   bins2 = binning_1D(angm_bins,angmfin,5,1,1)


   max1 =max(bins1[:])
   bins1 =bins1/max1
   bins2 =bins2/max1

   plot(angm_bins,bins1,linewidth=2,label='init migr')
   plot(angm_bins,bins2,linewidth=2,label='fin migr')
   #plot(dangm_bins,bins,label='tot-migr')


   legend(loc='best')
   grid('on')

   if save==1:
      savefig(savepath+'i-f-angm-migr.png')


#===========================================================================================

def plot_delta_r():
   """ histogram with delta r """
   figure()
   title(r'$\bf \Delta r$',fontsize=18)

   xlabel(r'$\Delta r$',fontsize=16)
   ylabel(r'$normed\ distribution$',fontsize=16)


   #h1 = hist(delta_rtot,200,range=(-10,10),normed='True',histtype='step',color='r',linewidth=2,label='tot')
   #h2 = hist(delta_r,200,range=(-10,10),normed='True',alpha=0.5,label='migr')


   #nb = 50
   dr_bins=[]
   for i in range(nDrad+1):
       dr_bins.append( Dradmin + (Dradmax - Dradmin)*float(i)/float(nDrad) ) 


   bins1 = binning_1D(dr_bins,delta_rtot,5,1,1)

   bins2 = binning_1D(dr_bins,delta_r,5,1,1)

   max1 =max(bins1[:])
   bins1 =bins1/max1
   bins2 =bins2/max1


   plot(dr_bins,bins1,linewidth=2,label='tot')
   plot(dr_bins,bins2,linewidth=2,label='migr')
   #plot(dr_bins,bins,label='tot-migr')





   def residuals(p, y, x):
       err = y-peval(x,p)
       return err
    
   def peval(x, p):
       return abs(p[0])*exp(-(x-p[1])**2/(2.*p[2]**2))+  abs(p[3])*exp(-(x-p[4])**2/(2.*p[5]**2))  
       #return p[0]*exp(-(x-p[1])**2/p[2])+  0.000009*exp(-(x-12.2)**2/3.)                     # gaussian + gaussian


     
 
   def fits(x,y,A1,x01,sig1,A2,x02,sig2):
      #A1=0.1
      #A2=0.1
      #n1=1.
      #n2=1.
  
      pname = (['A1','x01','sig1','A2','x02','sig2'])

      p0 = array([A1,x01,sig1, A2,x02,sig2 ])

      #pname = (['A1','x01','n1'])

      #p0 = array([A1 ,x01,n1])


      p , plsq = leastsq(residuals, p0, args=(y, x), maxfev=2000000)

      print
      print "Final parameters for fit"
      for i in range(len(pname)):
          print '%s = %.5f' % (pname[i],p[i])

      return p



   res = fits(dr_bins,bins2,0.1,-1.,1.,0.1,1.,1.)


   plot(dr_bins,peval(dr_bins,res),'r--',linewidth=1,label='fit')


   #res2 = fits(dr_bins,bins1,0.1,-0.,1.,0.,0.,0.1)

   #plot(dr_bins,peval(dr_bins,res2),'g--',linewidth=1,label='fit2')


   xlim(-10,10)
   legend(loc='best')
   grid('on')

   if save==1:
      savefig(savepath+'Dr.png')


#===========================================================================================

def plot_agemig():
   """ histogram with the age at wich the stars started migrating"""
   figure()
   title('time of migration')
   xlabel('time (Gyr)')
   ylabel('number')

   
   hist(agemig,200,range=(0,1))
#   hist(mtime,200,range=(4,5))
#   hist(mtime,50)
   
   #nage = 50
   time_bins = zeros(nage+1)
   for i in range(nage+1):
       time_bins[i] = 1. + 1./float(nage)*i


#   bins = binning_1D(time_bins, mtime, 5 , 0 , 1)

   #plot(age_bins,bins)
#   semilogy(time_bins,bins)

   grid('on')
   #xlim(4,5)

   if save==1:
      savefig(savepath+'agemig.png')


#===========================================================================================

def plot_posmig(t1,t2):
   """ plot x,y of when the stars have started migrating"""
   figure()
  
   title( str(t1)+'-'+str(t2) )

   w = 3.  # rad/Gyr 
   for i in range(len(mtime)):
     if (mtime[i]>t1) and (mtime[i]<t2):
       a = posmig[i,0]
       b = posmig[i,1]
       #posmigi[i,0] = a *cos(w*mtime[i]) - b *sin(w*mtime[i])
       #posmigi[i,1] = b *cos(w*mtime[i]) + a *sin(w*mtime[i])

       plot(posmig[i,0],posmig[i,1],'+b')
   grid('on')

   if save==1:
      savefig(savepath+'posmig.png')
  

#==============================================================================================



   
#==========================================================================







# plot of delta_L vs Linit  for stars detected as migrating ones

def plot_dl_linit_migr():

    figure()
    title(r'$\bf \Delta L\ vs\ L_{init.}\ -\ migr\ stars$',fontsize=18)

    ylabel(r'$\Delta L$',fontsize=16)
    xlabel(r'$L_{init.}$',fontsize=16)

    
    Lbins=zeros( nL +1)
    DLbins=zeros( nDL +1)

    for i in range(nL+1):
        Lbins[i] = Lmin + (Lmax - Lmin)/float(nL)*float(i) 

    for i in range(nDL+1):
        DLbins[i] = DLmin + (DLmax - DLmin)/float(nDL)*float(i) 


    bins = binning_2D(Lbins,angminit,5,3,DLbins,delta_angm,5,3,1)

    def bbb(x,y):
        return bins[x,y]


    I=range(len(Lbins))
    J=range(len(DLbins)) 
    X,Y = meshgrid(I, J)


    Z  = (bbb(X, Y))**(power)
    CS =contourf(Lbins[X], DLbins[Y], Z , 10)#, vmin=0.0 ,vmax=0.15)
    cb = colorbar(format='%.3f')
    cb.set_label(r'$\bf number^{%.2f}$'%(power),fontsize=16)

    axhline(y= 400.,color='r',linestyle='--')
    axhline(y=-400.,color='r',linestyle='--')
    axhline(y=0.,color='w',linestyle='--')


    line = 3000-2.*Lbins

    plot(Lbins,line,'k--')


    """
    open(savepath+'corot_stars.pyout','w').close()
    pyout = zeros(10)
    find=0
    i = 0
    while find<10 :    #and i<len(angminit):
        if (angminit[i] > 1100.) and (angminit[i] < 1400.) and (delta_angm[i] > 400.) and (delta_angm[i] < 500.) :
           pyout[find] = idmig[i]
  
           plot(angminit[i],delta_angm[i],'b+',markersize=10,linewidth=2)
           find = find +1
           print find

        i = i + 1 

    savetxt(path+'corot_stars.pyout',pyout,fmt='%i')



    open(savepath+'chaot_stars.pyout','w').close()
    pyout = zeros(10)
    find=0
    i = 0
    while find<10 and i<len(angminit):
        if (angminit[i] > 2071.) and (angminit[i] < 2280.) and (delta_angm[i] < 0.) and (delta_angm[i] > -140.) :
           pyout[find] = idmig[i]
           plot(angminit[i],delta_angm[i],'b+',markersize=10,linewidth=2)
           find = find +1
           print find


        i = i +1

    savetxt(path+'chaot_stars.pyout',pyout,fmt='%i')
    """


    ylim(DLmin,DLmax)

    if save==1:
      savefig(savepath+'DL-Linit-migr.png')


#==========================================================================
# plot of delta_guiding center vs Rinit  for stars detected as migrating ones

def plot_dgc_rinit_migr():

    figure()
    title(r'$\bf \Delta gc\ vs\ R_{init.}\ -\ migr\ stars$',fontsize=18)

    ylabel(r'$\Delta gc\ (kpc)$',fontsize=16)
    xlabel(r'$R_{init.}\ (kpc)$',fontsize=16)

    
    #nb = 50
    dr_bins=zeros(nDrad+1)
    for i in range(nDrad+1):
        dr_bins[i] = Dradmin + (Dradmax -Dradmin)*float(i)/float(nDrad)  


    #nrad = 50
    r_bins=zeros( nrad +1)
    for i in range(nrad+1):
        r_bins[i] = rmin + (rmax -rmin)/float(nrad)*float(i) 


    bins = binning_2D(r_bins,ri,5,1,dr_bins,delta_gc2,5,1,1.)

    for i in range(len(bins[:,0])):
        bins[i,:] = bins[i,:]/sum(bins[i,:])

#    for i in range(len(bins[:,0])):
#        bins[i,:] = bins[i,:]/max(bins[i,:])



    def bbb(x,y):
        return bins[x,y]


    I=range(len(r_bins))
    J=range(len(dr_bins)) 
    X,Y = meshgrid(I, J)


    Z  = (bbb(X, Y))**(power)
    CS =contourf(r_bins[X], dr_bins[Y], Z , 10)#, vmin=0.0 ,vmax=0.15)
    cb = colorbar(format='%.3f')
    cb.set_label(r'$\bf fraction^{%.2f}$'%(power),fontsize=16)

    #axhline(y= 400.,color='r',linestyle='--')
    #axhline(y=-400.,color='r',linestyle='--')
    axhline(y=0.,color='w',linestyle='--')

    #xlim(4.5,15)
    #ylim(-5,5)


    if save==1:
      savefig(savepath+'Dgc-Ri-migr.png')

#==========================================================================
# plot of delta_L vs Linit  for all stars

def plot_dl_linit_tot():
    

    figure()
    title(r'$\bf \Delta L\ vs\ L_{init.}\ -\ total\ stars$',fontsize=18)

    ylabel(r'$ \Delta L$',fontsize=16)
    xlabel(r'$L_{init.}$',fontsize=16)


    Lbins=zeros( nL +1)
    DLbins=zeros( nDL +1)

    for i in range(nL+1):
        Lbins[i] = Lmin + (Lmax - Lmin)/float(nL)*float(i) 

    for i in range(nDL+1):
        DLbins[i] = DLmin + (DLmax - DLmin)/float(nDL)*float(i) 


    bins = binning_2D(Lbins,angmitot,5,3,DLbins,delta_angmtot,5,3,1)

    def bbb(x,y):
        return bins[x,y]

    I=range(len(Lbins))
    J=range(len(DLbins)) 
    X,Y = meshgrid(I, J)


    Z  = (bbb(X, Y))**(power)
    CS =contourf(Lbins[X], DLbins[Y], Z , 10, vmin=0.0 ,vmax=3.6)
    cb = colorbar(format='%.3f')
    cb.set_label(r'$\bf number^{%.2f}$'%(power),fontsize=16)

    axhline(y=400.,color='r',linestyle='--')
    axhline(y=-400.,color='r',linestyle='--')
    axhline(y=0.,color='w',linestyle='--')


    line = 3000-2.*Lbins

    plot(Lbins,line,'k--')



    dl_mean=zeros(len(Lbins))

    for i in range(len(Lbins)):
        mm=0.
        for j in range(len(DLbins)-1):
            dl_mean[i] = dl_mean[i] + bins[i,j]*(DLbins[j]+DLbins[j+1])/2.
            mm=mm+bins[i,j]
        dl_mean[i]=dl_mean[i]/mm


    plot(Lbins,dl_mean,'w--',linewidth=2,label='average DL')




    ylim(DLmin,DLmax)

    if save==1:
      savefig(savepath+'DL-Li-all.png')





#==========================================================================
# plot of delta_L vs Linit  vs tangential velocity for all stars

def plot_dl_linit_tot_tvel():
    

    figure()
    title(r'$\bf \Delta L\ vs\ L_{init.}\ tvel\ -\ total\ stars$',fontsize=18)

    ylabel(r'$ \Delta L$',fontsize=16)
    xlabel(r'$L_{init.}$',fontsize=16)


    Lbins    = zeros( nL +1)
    DLbins   = zeros( nDL +1)
    vel_bins = zeros( nvel +1)


    for i in range(nL+1):
        Lbins[i] = Lmin + (Lmax - Lmin)/float(nL)*float(i) 

    for i in range(nDL+1):
        DLbins[i] = DLmin + (DLmax - DLmin)/float(nDL)*float(i) 

    for i in range(nvel+1):
        vel_bins[i] = velmin_t+ (velmax_t - velmin_t)/float(nvel)*float(i) 


    #bins = binning_2D(Lbins,angmitot,5,3,DLbins,delta_angmtot,5,3,1)
    bins = binning_3D(Lbins,angmitot,5,3,DLbins,delta_angmtot,5,3,vel_bins,tvelitot,5,2,1.)
    #bins2 = binning_3D(Lbins,angmitot,5,3,DLbins,delta_angmtot,5,3,vel_bins,tvelitot,5,2,tvelitot)


    def bbb(x,y):
        a = 0.
        tot = bins[x,y,:].sum(axis=2)
      
        for i in range( len ( bins[x,y,:] ) ):
            a = a + vel_bins[i] * bins[x,y,i] / tot

        return a
        #return bins[x,y]

    I=range(len(Lbins))
    J=range(len(DLbins)) 
    X,Y = meshgrid(I, J)


    Z  = (bbb(X, Y))**(1)
    CS =contourf(Lbins[X], DLbins[Y], Z , 10, vmin=0.0 ,vmax=3.6)
    cb = colorbar(format='%.3f')
    cb.set_label(r'$\bf tang.\ vel\ (km.s^{-1})$',fontsize=16)

    axhline(y=400.,color='r',linestyle='--')
    axhline(y=-400.,color='r',linestyle='--')
    axhline(y=0.,color='w',linestyle='--')


    line = 3000-2.*Lbins

    plot(Lbins,line,'k--')


    """
    dl_mean=zeros(len(Lbins))

    for i in range(len(Lbins)):
        mm=0.
        for j in range(len(DLbins)-1):
            dl_mean[i] = dl_mean[i] + bins[i,j]*(DLbins[j]+DLbins[j+1])/2.
            mm=mm+bins[i,j]
        dl_mean[i]=dl_mean[i]/mm


    plot(Lbins,dl_mean,'w--',linewidth=2,label='average DL')
    """



    ylim(DLmin,DLmax)

    if save==1:
      savefig(savepath+'dl_linit_tot_tvel.png')



#==========================================================================
# plot of delta_L vs Linit  vs radial velocity for all stars

def plot_dl_linit_tot_rvel():
    

    figure()
    title(r'$\bf \Delta L\ vs\ L_{init.}\ radveli\ -\ total\ stars$',fontsize=18)

    ylabel(r'$ \Delta L$',fontsize=16)
    xlabel(r'$L_{init.}$',fontsize=16)


    Lbins    = zeros( nL +1)
    DLbins   = zeros( nDL +1)
    vel_bins = zeros( nvel +1)


    for i in range(nL+1):
        Lbins[i] = Lmin + (Lmax - Lmin)/float(nL)*float(i) 

    for i in range(nDL+1):
        DLbins[i] = DLmin + (DLmax - DLmin)/float(nDL)*float(i) 

    for i in range(nvel+1):
        vel_bins[i] = velmin_r+ (velmax_r - velmin_r)/float(nvel)*float(i) 


    #bins = binning_2D(Lbins,angmitot,5,3,DLbins,delta_angmtot,5,3,1)
    bins = binning_3D(Lbins,angmitot,5,3,DLbins,delta_angmtot,5,3,vel_bins,rvelitot,5,1,1.)
    #bins2 = binning_3D(Lbins,angmitot,5,3,DLbins,delta_angmtot,5,3,vel_bins,tvelitot,5,2,tvelitot)


    def bbb(x,y):

        a = 0.
        tot = bins[x,y,:].sum(axis=2)

        for i in range( len ( bins[x,y,:] ) ):
            a = a + vel_bins[i] * bins[x,y,i] / tot

        return a
        #return bins[x,y]

    I=range(len(Lbins))
    J=range(len(DLbins)) 
    X,Y = meshgrid(I, J)


    Z  = (bbb(X, Y))**(1)
    CS =contourf(Lbins[X], DLbins[Y], Z , 10, vmin=0.0 ,vmax=3.6)
    cb = colorbar(format='%.3f')
    cb.set_label(r'$\bf rad.\ vel\ (km.s^{-1})$',fontsize=16)

    axhline(y=400.,color='r',linestyle='--')
    axhline(y=-400.,color='r',linestyle='--')
    axhline(y=0.,color='w',linestyle='--')


    line = 3000-2.*Lbins

    plot(Lbins,line,'k--')


    """
    dl_mean=zeros(len(Lbins))

    for i in range(len(Lbins)):
        mm=0.
        for j in range(len(DLbins)-1):
            dl_mean[i] = dl_mean[i] + bins[i,j]*(DLbins[j]+DLbins[j+1])/2.
            mm=mm+bins[i,j]
        dl_mean[i]=dl_mean[i]/mm


    plot(Lbins,dl_mean,'w--',linewidth=2,label='average DL')
    """



    ylim(DLmin,DLmax)

    if save==1:
      savefig(savepath+'dl_linit_tot_rvel.png')






#==========================================================================
# plot of delta_L vs Linit  vs tangential velocity dispersion for all stars

def plot_dl_linit_tot_disptvel():
    

    figure()
    title(r'$\bf \Delta L\ vs\ L_{init.}\ tvel\ disp\ -\ total\ stars$',fontsize=18)

    ylabel(r'$ \Delta L$',fontsize=16)
    xlabel(r'$L_{init.}$',fontsize=16)


    Lbins    = zeros( nL +1)
    DLbins   = zeros( nDL +1)
    vel_bins = zeros( nvel +1)


    for i in range(nL+1):
        Lbins[i] = Lmin + (Lmax - Lmin)/float(nL)*float(i) 

    for i in range(nDL+1):
        DLbins[i] = DLmin + (DLmax - DLmin)/float(nDL)*float(i) 

    for i in range(nvel+1):
        vel_bins[i] = velmin_t+ (velmax_t - velmin_t)/float(nvel)*float(i) 


    #bins = binning_2D(Lbins,angmitot,5,3,DLbins,delta_angmtot,5,3,1)
    bins = binning_3D(Lbins,angmitot,5,3,DLbins,delta_angmtot,5,3,vel_bins,tvelitot,5,2,1.)
    #bins2 = binning_3D(Lbins,angmitot,5,3,DLbins,delta_angmtot,5,3,vel_bins,tvelitot,5,2,tvelitot)


    def bbb(x,y):
        a = 0.
        b = 0.
        tot = bins[x,y,:].sum(axis=2)
      
        for i in range( len ( bins[x,y,:] ) ):
            a = a + vel_bins[i] * bins[x,y,i] / tot

        for i in range( len ( bins[x,y,:] ) ):
            b = b + vel_bins[i]**2 * bins[x,y,i] / tot


        c = sqrt ( b - a**2 )

        return c
        #return bins[x,y]

    I=range(len(Lbins))
    J=range(len(DLbins)) 
    X,Y = meshgrid(I, J)


    Z  = (bbb(X, Y))**(1)
    CS =contourf(Lbins[X], DLbins[Y], Z , 10, vmin=0.0 ,vmax=3.6)
    cb = colorbar(format='%.3f')
    cb.set_label(r'$\bf tang.\ vel\ disp\ (km.s^{-1})$',fontsize=16)

    axhline(y=400.,color='r',linestyle='--')
    axhline(y=-400.,color='r',linestyle='--')
    axhline(y=0.,color='w',linestyle='--')


    line = 3000-2.*Lbins

    plot(Lbins,line,'k--')


    """
    dl_mean=zeros(len(Lbins))

    for i in range(len(Lbins)):
        mm=0.
        for j in range(len(DLbins)-1):
            dl_mean[i] = dl_mean[i] + bins[i,j]*(DLbins[j]+DLbins[j+1])/2.
            mm=mm+bins[i,j]
        dl_mean[i]=dl_mean[i]/mm


    plot(Lbins,dl_mean,'w--',linewidth=2,label='average DL')
    """



    ylim(DLmin,DLmax)

    if save==1:
      savefig(savepath+'dl_linit_tot_disptvel.png')









#==========================================================================
# plot of delta_L vs Linit  vs radial velocity dispersion for all stars

def plot_dl_linit_tot_disprvel():
    

    figure()
    title(r'$\bf \Delta L\ vs\ L_{init.}\ rvel\ disp\ -\ total\ stars$',fontsize=18)

    ylabel(r'$ \Delta L$',fontsize=16)
    xlabel(r'$L_{init.}$',fontsize=16)


    Lbins    = zeros( nL +1)
    DLbins   = zeros( nDL +1)
    vel_bins = zeros( nvel +1)


    for i in range(nL+1):
        Lbins[i] = Lmin + (Lmax - Lmin)/float(nL)*float(i) 

    for i in range(nDL+1):
        DLbins[i] = DLmin + (DLmax - DLmin)/float(nDL)*float(i) 

    for i in range(nvel+1):
        vel_bins[i] = velmin_r+ (velmax_r - velmin_r)/float(nvel)*float(i) 


    #bins = binning_2D(Lbins,angmitot,5,3,DLbins,delta_angmtot,5,3,1)
    bins = binning_3D(Lbins,angmitot,5,3,DLbins,delta_angmtot,5,3,vel_bins,rvelitot,5,2,1.)
    #bins2 = binning_3D(Lbins,angmitot,5,3,DLbins,delta_angmtot,5,3,vel_bins,tvelitot,5,2,tvelitot)


    def bbb(x,y):
        a = 0.
        b = 0.
        tot = bins[x,y,:].sum(axis=2)
      
        for i in range( len ( bins[x,y,:] ) ):
            a = a + vel_bins[i] * bins[x,y,i] / tot

        for i in range( len ( bins[x,y,:] ) ):
            b = b + vel_bins[i]**2 * bins[x,y,i] / tot


        c = sqrt ( b - a**2 )

        return c
        #return bins[x,y]

    I=range(len(Lbins))
    J=range(len(DLbins)) 
    X,Y = meshgrid(I, J)


    Z  = (bbb(X, Y))**(1)
    CS =contourf(Lbins[X], DLbins[Y], Z , 10, vmin=0.0 ,vmax=3.6)
    cb = colorbar(format='%.3f')
    cb.set_label(r'$\bf rad.\ vel\ disp\ (km.s^{-1})$',fontsize=16)

    axhline(y=400.,color='r',linestyle='--')
    axhline(y=-400.,color='r',linestyle='--')
    axhline(y=0.,color='w',linestyle='--')


    line = 3000-2.*Lbins

    plot(Lbins,line,'k--')


    """
    dl_mean=zeros(len(Lbins))

    for i in range(len(Lbins)):
        mm=0.
        for j in range(len(DLbins)-1):
            dl_mean[i] = dl_mean[i] + bins[i,j]*(DLbins[j]+DLbins[j+1])/2.
            mm=mm+bins[i,j]
        dl_mean[i]=dl_mean[i]/mm


    plot(Lbins,dl_mean,'w--',linewidth=2,label='average DL')
    """



    ylim(DLmin,DLmax)

    if save==1:
      savefig(savepath+'dl_linit_tot_disprvel.png')

#==========================================================================
# plot of delta_L vs Linit  vs R final for all stars

def plot_dl_linit_tot_rf():
    

    figure()
    title(r'$\bf \Delta L\ vs\ L_{init.}\ R\ fin\ -\ total\ stars$',fontsize=18)

    ylabel(r'$ \Delta L$',fontsize=16)
    xlabel(r'$L_{init.}$',fontsize=16)


    Lbins    = zeros( nL +1)
    DLbins   = zeros( nDL +1)
    r_bins=zeros(nrad+1)


    for i in range(nL+1):
        Lbins[i] = Lmin + (Lmax - Lmin)/float(nL)*float(i) 

    for i in range(nDL+1):
        DLbins[i] = DLmin + (DLmax - DLmin)/float(nDL)*float(i) 

    for i in range(nrad+1):
        r_bins[i] = rmin + (rmax -rmin)*float(i)/float(nrad)  



    #bins = binning_2D(Lbins,angmitot,5,3,DLbins,delta_angmtot,5,3,1)
    bins = binning_3D(Lbins,angmitot,5,3,DLbins,delta_angmtot,5,3,r_bins,rftot,5,1,1.)
    #bins2 = binning_3D(Lbins,angmitot,5,3,DLbins,delta_angmtot,5,3,vel_bins,tvelitot,5,2,tvelitot)


    def bbb(x,y):
        a = 0.
        tot = bins[x,y,:].sum(axis=2)
      
        for i in range( len ( bins[x,y,:] ) ):
            a = a + r_bins[i] * bins[x,y,i] / tot

        return a
        #return bins[x,y]

    I=range(len(Lbins))
    J=range(len(DLbins)) 
    X,Y = meshgrid(I, J)


    Z  = (bbb(X, Y))**(1)
    CS =contourf(Lbins[X], DLbins[Y], Z , 10, vmin=0.0 ,vmax=3.6)
    cb = colorbar(format='%.3f')
    cb.set_label(r'$\bf R\ fin\ (kpc)$',fontsize=16)

    axhline(y=400.,color='r',linestyle='--')
    axhline(y=-400.,color='r',linestyle='--')
    axhline(y=0.,color='w',linestyle='--')


    line = 3000-2.*Lbins

    plot(Lbins,line,'k--')


    """
    dl_mean=zeros(len(Lbins))

    for i in range(len(Lbins)):
        mm=0.
        for j in range(len(DLbins)-1):
            dl_mean[i] = dl_mean[i] + bins[i,j]*(DLbins[j]+DLbins[j+1])/2.
            mm=mm+bins[i,j]
        dl_mean[i]=dl_mean[i]/mm


    plot(Lbins,dl_mean,'w--',linewidth=2,label='average DL')
    """



    ylim(DLmin,DLmax)

    if save==1:
      savefig(savepath+'dl_linit_tot_rf.png')







#==========================================================================


#==========================================================================
# plot of delta_L vs Linit  vs R initial for all stars

def plot_dl_linit_tot_ri():
    

    figure()
    title(r'$\bf \Delta L\ vs\ L_{init.}\ R\ init\ -\ total\ stars$',fontsize=18)

    ylabel(r'$ \Delta L$',fontsize=16)
    xlabel(r'$L_{init.}$',fontsize=16)


    Lbins    = zeros( nL +1)
    DLbins   = zeros( nDL +1)
    r_bins=zeros(nrad+1)


    for i in range(nL+1):
        Lbins[i] = Lmin + (Lmax - Lmin)/float(nL)*float(i) 

    for i in range(nDL+1):
        DLbins[i] = DLmin + (DLmax - DLmin)/float(nDL)*float(i) 

    for i in range(nrad+1):
        r_bins[i] = rmin + (rmax -rmin)*float(i)/float(nrad)  



    #bins = binning_2D(Lbins,angmitot,5,3,DLbins,delta_angmtot,5,3,1)
    bins = binning_3D(Lbins,angmitot,5,3,DLbins,delta_angmtot,5,3,r_bins,ritot,5,1,1.)
    #bins2 = binning_3D(Lbins,angmitot,5,3,DLbins,delta_angmtot,5,3,vel_bins,tvelitot,5,2,tvelitot)


    def bbb(x,y):
        a = 0.
        tot = bins[x,y,:].sum(axis=2)
      
        for i in range( len ( bins[x,y,:] ) ):
            a = a + r_bins[i] * bins[x,y,i] / tot

        return a
        #return bins[x,y]

    I=range(len(Lbins))
    J=range(len(DLbins)) 
    X,Y = meshgrid(I, J)


    Z  = (bbb(X, Y))**(1)
    CS =contourf(Lbins[X], DLbins[Y], Z , 10, vmin=0.0 ,vmax=3.6)
    cb = colorbar(format='%.3f')
    cb.set_label(r'$\bf R\ init\ (kpc)$',fontsize=16)

    axhline(y=400.,color='r',linestyle='--')
    axhline(y=-400.,color='r',linestyle='--')
    axhline(y=0.,color='w',linestyle='--')


    line = 3000-2.*Lbins

    plot(Lbins,line,'k--')


    """
    dl_mean=zeros(len(Lbins))

    for i in range(len(Lbins)):
        mm=0.
        for j in range(len(DLbins)-1):
            dl_mean[i] = dl_mean[i] + bins[i,j]*(DLbins[j]+DLbins[j+1])/2.
            mm=mm+bins[i,j]
        dl_mean[i]=dl_mean[i]/mm


    plot(Lbins,dl_mean,'w--',linewidth=2,label='average DL')
    """



    ylim(DLmin,DLmax)

    if save==1:
      savefig(savepath+'dl_linit_tot_ri.png')







#==========================================================================
# plot of delta_L vs Linit  vs delta R for all stars

def plot_dl_linit_tot_dr():
    

    figure()
    title(r'$\bf \Delta L\ vs\ L_{init.}\ dr\ -\ total\ stars$',fontsize=18)

    ylabel(r'$ \Delta L$',fontsize=16)
    xlabel(r'$L_{init.}$',fontsize=16)


    Lbins    = zeros( nL +1)
    DLbins   = zeros( nDL +1)
    dr_bins=zeros(nDrad+1)


    for i in range(nL+1):
        Lbins[i] = Lmin + (Lmax - Lmin)/float(nL)*float(i) 

    for i in range(nDL+1):
        DLbins[i] = DLmin + (DLmax - DLmin)/float(nDL)*float(i) 

    for i in range(nDrad+1):
        dr_bins[i] = Dradmin + (Dradmax -Dradmin)*float(i)/float(nDrad)  



    #bins = binning_2D(Lbins,angmitot,5,3,DLbins,delta_angmtot,5,3,1)
    bins = binning_3D(Lbins,angmitot,5,3,DLbins,delta_angmtot,5,3,dr_bins,delta_rtot,5,1,1.)
    #bins2 = binning_3D(Lbins,angmitot,5,3,DLbins,delta_angmtot,5,3,vel_bins,tvelitot,5,2,tvelitot)


    def bbb(x,y):
        a = 0.
        tot = bins[x,y,:].sum(axis=2)
      
        for i in range( len ( bins[x,y,:] ) ):
            a = a + dr_bins[i] * bins[x,y,i] / tot

        return a
        #return bins[x,y]

    I=range(len(Lbins))
    J=range(len(DLbins)) 
    X,Y = meshgrid(I, J)


    Z  = (bbb(X, Y))**(1)
    CS =contourf(Lbins[X], DLbins[Y], Z , 10, vmin=0.0 ,vmax=3.6)
    cb = colorbar(format='%.3f')
    cb.set_label(r'$\bf \Delta R\ (kpc)$',fontsize=16)

    axhline(y=400.,color='r',linestyle='--')
    axhline(y=-400.,color='r',linestyle='--')
    axhline(y=0.,color='w',linestyle='--')


    line = 3000-2.*Lbins

    plot(Lbins,line,'k--')


    """
    dl_mean=zeros(len(Lbins))

    for i in range(len(Lbins)):
        mm=0.
        for j in range(len(DLbins)-1):
            dl_mean[i] = dl_mean[i] + bins[i,j]*(DLbins[j]+DLbins[j+1])/2.
            mm=mm+bins[i,j]
        dl_mean[i]=dl_mean[i]/mm


    plot(Lbins,dl_mean,'w--',linewidth=2,label='average DL')
    """



    ylim(DLmin,DLmax)

    if save==1:
      savefig(savepath+'dl_linit_tot_dr.png')







#==========================================================================
# plot of delta_L vs Linit  vs radial velocity dispersion for all stars

def plot_dl_linit_tot_dispdr():
    

    figure()
    title(r'$\bf \Delta L\ vs\ L_{init.}\ dr\ disp\ -\ total\ stars$',fontsize=18)

    ylabel(r'$ \Delta L$',fontsize=16)
    xlabel(r'$L_{init.}$',fontsize=16)


    Lbins    = zeros( nL +1)
    DLbins   = zeros( nDL +1)
    dr_bins=zeros(nDrad+1)


    for i in range(nL+1):
        Lbins[i] = Lmin + (Lmax - Lmin)/float(nL)*float(i) 

    for i in range(nDL+1):
        DLbins[i] = DLmin + (DLmax - DLmin)/float(nDL)*float(i) 

    for i in range(nDrad+1):
        dr_bins[i] = Dradmin + (Dradmax -Dradmin)*float(i)/float(nDrad)  



    #bins = binning_2D(Lbins,angmitot,5,3,DLbins,delta_angmtot,5,3,1)
    bins = binning_3D(Lbins,angmitot,5,3,DLbins,delta_angmtot,5,3,dr_bins,delta_rtot,5,1,1.)
    #bins2 = binning_3D(Lbins,angmitot,5,3,DLbins,delta_angmtot,5,3,vel_bins,tvelitot,5,2,tvelitot)


    def bbb(x,y):
        a = 0.
        b = 0.
        tot = bins[x,y,:].sum(axis=2)
      
        for i in range( len ( bins[x,y,:] ) ):
            a = a + dr_bins[i] * bins[x,y,i] / tot

        for i in range( len ( bins[x,y,:] ) ):
            b = b + dr_bins[i]**2 * bins[x,y,i] / tot

        c = sqrt( b - a**2 )

        return c
        #return bins[x,y]

    I=range(len(Lbins))
    J=range(len(DLbins)) 
    X,Y = meshgrid(I, J)


    Z  = (bbb(X, Y))**(1)
    CS =contourf(Lbins[X], DLbins[Y], Z , 10, vmin=0.0 ,vmax=3.6)
    cb = colorbar(format='%.3f')
    cb.set_label(r'$\bf \Delta R\ disp\ (kpc)$',fontsize=16)

    axhline(y=400.,color='r',linestyle='--')
    axhline(y=-400.,color='r',linestyle='--')
    axhline(y=0.,color='w',linestyle='--')


    line = 3000-2.*Lbins

    plot(Lbins,line,'k--')


    """
    dl_mean=zeros(len(Lbins))

    for i in range(len(Lbins)):
        mm=0.
        for j in range(len(DLbins)-1):
            dl_mean[i] = dl_mean[i] + bins[i,j]*(DLbins[j]+DLbins[j+1])/2.
            mm=mm+bins[i,j]
        dl_mean[i]=dl_mean[i]/mm


    plot(Lbins,dl_mean,'w--',linewidth=2,label='average DL')
    """



    ylim(DLmin,DLmax)

    if save==1:
      savefig(savepath+'dl_linit_tot_dispdr.png')




#==========================================================================
# plot of delta_guiding center vs Rinit  for all stars 

def plot_dgc_rinit_tot():

    figure()
    title(r'$\bf \Delta gc\ vs\ R_{init.}\ -\ total\ stars$',fontsize=18)

    ylabel(r'$\Delta gc\ (kpc)$',fontsize=16)
    xlabel(r'$R_{init.}\ (kpc)$',fontsize=16)

    
    #nb = 50
    dr_bins=zeros(nDrad+1)
    for i in range(nDrad+1):
        dr_bins[i] = Dradmin + (Dradmax - Dradmin)*float(i)/float(nDrad)  


    #nrad = 50
    r_bins=zeros( nrad +1)
    for i in range(nrad+1):
        r_bins[i] = rmin + (rmax -rmin)/float(nrad)*float(i) 


    bins = binning_2D(r_bins,ritot,5,1,dr_bins,delta_gc2tot,5,1,1.)

    for i in range(len(bins[:,0])):
        bins[i,:] = bins[i,:]/sum(bins[i,:])

#    for i in range(len(bins[:,0])):
#        bins[i,:] = bins[i,:]/max(bins[i,:])


    def bbb(x,y):
        return bins[x,y]


    I=range(len(r_bins))
    J=range(len(dr_bins)) 
    X,Y = meshgrid(I, J)


    Z  =(bbb(X, Y))**(power)
    CS =contourf(r_bins[X], dr_bins[Y], Z , 10, vmin=0.0 ,vmax=0.1)
    cb = colorbar(format='%.3f')
    cb.set_label(r'$\bf fraction^{%.2f}$'%(power),fontsize=16)

    #axhline(y= 400.,color='r',linestyle='--')
    #axhline(y=-400.,color='r',linestyle='--')
    axhline(y=0.,color='w',linestyle='--')

    #xlim(4.5,15)
    #ylim(-10,10)


    if save==1:
      savefig(savepath+'Dgc-Ri-all.png')


#==========================================================================
# plot of delta_L vs Lfin  for stars detected as migrating ones

def plot_dl_lfin_migr():

    figure()
    title(r'$\bf \Delta L\ vs\ L_{fin.}\ -\ migr\ stars$',fontsize=18)

    ylabel(r'$\Delta L$',fontsize=16)
    xlabel(r'$L_{fin.}$',fontsize=16)

    
    Lbins=zeros( nL +1)
    DLbins=zeros( nDL +1)

    for i in range(nL+1):
        Lbins[i] = Lmin + (Lmax - Lmin)/float(nL)*float(i) 

    for i in range(nDL+1):
        DLbins[i] = DLmin + (DLmax - DLmin)/float(nDL)*float(i) 


    bins = binning_2D(Lbins,angmfin,5,3,DLbins,delta_angm,5,3,1)

    def bbb(x,y):
        return bins[x,y]


    I=range(len(Lbins))
    J=range(len(DLbins)) 
    X,Y = meshgrid(I, J)


    Z  = (bbb(X, Y))**(power)
    CS =contourf(Lbins[X], DLbins[Y], Z , 10)#, vmin=0.0 ,vmax=0.15)
    cb = colorbar(format='%.3f')
    cb.set_label(r'$\bf number^{%.2f}$'%(power),fontsize=16)

    axhline(y= 400.,color='r',linestyle='--')
    axhline(y=-400.,color='r',linestyle='--')
    axhline(y=0.,color='w',linestyle='--')

    if save==1:
      savefig(savepath+'DL-Lf-migr.png')


#==========================================================================
# plot of delta_L vs Lfin  for all stars

def plot_dl_lfin_tot():
    

    figure()
    title(r'$\bf \Delta L\ vs\ L_{fin.}\ -\ total\ stars$',fontsize=18)

    ylabel(r'$ \Delta L$',fontsize=16)
    xlabel(r'$L_{fin.}$',fontsize=16)


    Lbins=zeros( nL +1)
    DLbins=zeros( nDL +1)

    for i in range(nL+1):
        Lbins[i] = Lmin + (Lmax - Lmin)/float(nL)*float(i) 

    for i in range(nDL+1):
        DLbins[i] = DLmin + (DLmax - DLmin)/float(nDL)*float(i) 


    bins = binning_2D(Lbins,angmftot,5,3,DLbins,delta_angmtot,5,3,1)

    def bbb(x,y):
        return bins[x,y]

    I=range(len(Lbins))
    J=range(len(DLbins)) 
    X,Y = meshgrid(I, J)


    Z  = (bbb(X, Y))**(power)
    CS =contourf(Lbins[X], DLbins[Y], Z , 10, vmin=0.0 ,vmax=3.6)
    cb = colorbar(format='%.3f')
    cb.set_label(r'$\bf number^{%.2f}$'%(power),fontsize=16)

    axhline(y=400.,color='r',linestyle='--')
    axhline(y=-400.,color='r',linestyle='--')
    axhline(y=0.,color='w',linestyle='--')


    if save==1:
      savefig(savepath+'DL-Lf-all.png')




#==========================================================================
# plot of delta_L vs Rfin  for all stars
def plot_dl_Rfin_tot():
    

    figure()
    title(r'$\bf \Delta L\ vs\ R_{fin.}\ -\ total\ stars$',fontsize=18)

    ylabel(r'$ \Delta L$',fontsize=16)
    xlabel(r'$R_{fin.}$',fontsize=16)


    Rbins=zeros( nrad +1)
    DLbins=zeros( nDL +1)

    for i in range(nrad+1):
        Rbins[i] = rmin + (rmax - rmin)/float(nrad)*float(i) 

    for i in range(nDL+1):
        DLbins[i] = DLmin + (DLmax - DLmin)/float(nDL)*float(i) 


    bins = binning_2D(Rbins,rftot,5,1,DLbins,delta_angmtot,5,3,1)

    def bbb(x,y):
        return bins[x,y]

    I=range(len(Rbins))
    J=range(len(DLbins)) 
    X,Y = meshgrid(I, J)


    Z  = (bbb(X, Y))**(power)
    #CS =contourf(Rbins[X], DLbins[Y], Z , 10, vmin=0.0 ,vmax=3.6)
    CS =pcolor(Rbins[X], DLbins[Y], Z ,norm=LogNorm(vmin=1., vmax=Z.max() ) )
    cb = colorbar(format='%.3f')
    cb.set_label(r'$\bf number^{%.2f}$'%(power),fontsize=16)

    axhline(y=400.,color='r',linestyle='--')
    axhline(y=-400.,color='r',linestyle='--')
    axhline(y=0.,color='w',linestyle='--')


    if save==1:
      savefig(savepath+'DL-Rf-all.png')






#==========================================================================

#==========================================================================
# plot of delta_L vs time of migration beginning for selected stars

def plot_dl_mt():
    

    figure()
    title(r'$\bf \Delta L\ vs\ migr\ time\ -\ migr stars$',fontsize=18)

    ylabel(r'$ \Delta L$',fontsize=16)
    xlabel(r'$migr\ time (Gyr)$',fontsize=16)


    Lbins=zeros( nL +1)
    DLbins=zeros( nDL +1)

    #nage = 50
    age_bins = zeros(nage+1)
    for i in range(nage+1):
        age_bins[i] = 0. + 1./float(nage)*i


    for i in range(nL+1):
        Lbins[i] = Lmin + (Lmax - Lmin)/float(nL)*float(i) 

    for i in range(nDL+1):
        DLbins[i] = DLmin + (DLmax - DLmin)/float(nDL)*float(i) 


    bins = binning_2D(age_bins,agemig,5,0,DLbins,delta_angm,5,3,1)

   
    #bins = bins/float(len(ri))
   
   
    def bbb(x,y):
        return bins[x,y]

    I=range(len(age_bins))
    J=range(len(DLbins)) 
    X,Y = meshgrid(I, J)


    Z  = (bbb(X, Y))**(power)
    CS =contourf(age_bins[X], DLbins[Y], Z )#, vmin=0.0 ,vmax=3.6)
    cb = colorbar(format='%.3f')
    cb.set_label(r'$\bf number^{%.2f}$'%(power),fontsize=16)

    axhline(y=200.,color='r',linestyle='--')
    axhline(y=-200.,color='r',linestyle='--')
    axhline(y=0.,color='w',linestyle='--')

    line1 = 2000-1600.*age_bins
    line2 = 3000-3400.*age_bins
    line3 = -2000+1600.*age_bins
    line4 = -3000+3400.*age_bins

    plot(age_bins,line1,'w--')
    plot(age_bins,line2,'w--')
    plot(age_bins,line3,'w--')
    plot(age_bins,line4,'w--')

    ylim(DLmin,DLmax)

    if save==1:
      savefig(savepath+'DL-mtime-migr.png')

#==========================================================================
# plot of delta_R vs time of migration beginning for selected stars

def plot_dr_mt():
    

    figure()
    title(r'$\bf \Delta R\ vs\ migr\ time\ -\ migr stars$',fontsize=18)

    ylabel(r'$ \Delta R$',fontsize=16)
    xlabel(r'$migr\ time (Gyr)$',fontsize=16)



    #nage = 50
    age_bins = zeros(nage+1)
    for i in range(nage+1):
        age_bins[i] = 0. + 1./float(nage)*i


    #nb = 50
    dr_bins=zeros(nDrad+1)
    for i in range(nDrad+1):
        dr_bins[i] = Dradmin + (Dradmax - Dradmin)*float(i)/float(nDrad)  



    bins = binning_2D(age_bins,agemig,7,0,dr_bins,delta_r,5,1,1)
   
   
    def bbb(x,y):
        return bins[x,y]

    I=range(len(age_bins))
    J=range(len(dr_bins)) 
    X,Y = meshgrid(I, J)


    Z  = (bbb(X, Y))**(power)
    CS =contourf(age_bins[X], dr_bins[Y], Z )#, vmin=0.0 ,vmax=3.6)
    cb = colorbar(format='%.3f')
    cb.set_label(r'$\bf number^{%.2f}$'%(power),fontsize=16)

    axhline(y=0.,color='w',linestyle='--')

    if save==1:
      savefig(savepath+'Dr-mtime-migr.png')

#==========================================================================
# plot of delta_R vs Rinit  for stars detected as migrating ones

def plot_dr_rinit_migr():

    figure()
    title(r'$\bf \Delta R\ vs\ R_{init.}\ -\ migr\ stars$',fontsize=18)

    ylabel(r'$\Delta R\ (kpc)$',fontsize=16)
    xlabel(r'$R_{init.}\ (kpc)$',fontsize=16)

   

    #nb = 50
    dr_bins=zeros( nDrad + 1)
    for i in range(nDrad+1):
        dr_bins[i] =  Dradmin + (Dradmax - Dradmin)*float(i)/float(nDrad) 

    #nR = 50
    r_bins = zeros(nrad+1)
    for i in range(nrad+1):
        r_bins[i] =  rmin + (rmax - rmin)*float(i)/float(nrad) 



    bins = binning_2D(r_bins,ri,5,1,dr_bins,delta_r,5,1,1)


    for i in range(len(bins[:,0])):
        bins[i,:] = bins[i,:]/sum(bins[i,:])

    def bbb(x,y):
        return bins[x,y]


    I=range(len(r_bins))
    J=range(len(dr_bins)) 
    X,Y = meshgrid(I, J)


    Z  = (bbb(X, Y))**(power)
    CS =contourf(r_bins[X], dr_bins[Y], Z , 10)#, vmin=0.0 ,vmax=0.15)
    cb = colorbar(format='%.3f')
    cb.set_label(r'$\bf fraction^{%.2f} $'%(power),fontsize=16)

    #axhline(y= 400.,color='r',linestyle='--')
    #axhline(y=-400.,color='r',linestyle='--')
    axhline(y=0.,color='w',linestyle='--')

    ylim(Dradmin,Dradmax)
    xlim(5,rmax)


    if save==1:
      savefig(savepath+'Dr-Ri-migr.png')

#==========================================================================
# plot of delta_R vs Rinit  for all stars

def plot_dr_rinit_tot():
    

    figure()
    title(r'$\bf \Delta R\ vs\ R_{init.}\ -\ total\ stars$',fontsize=18)

    ylabel(r'$ \Delta R\ (kpc)$',fontsize=16)
    xlabel(r'$R_{init.}\ (kpc)$',fontsize=16)


    #nb = 50
    dr_bins=zeros( nDrad + 1)
    for i in range(nDrad+1):
        dr_bins[i] =  Dradmin + (Dradmax - Dradmin)*float(i)/float(nDrad) 

    #nR = 50
    r_bins = zeros(nrad+1)
    for i in range(nrad+1):
        r_bins[i] =  rmin + (rmax - rmin)*float(i)/float(nrad) 


    bins = binning_2D(r_bins,ritot,5,1,dr_bins,delta_rtot,5,0,1.)

    for i in range(len(bins[:,0])):
        bins[i,:]=bins[i,:]/sum(bins[i,:])

    def bbb(x,y):
        return bins[x,y]

    I=range(len(r_bins))
    J=range(len(dr_bins)) 
    X,Y = meshgrid(I, J)


    Z  = (bbb(X, Y))**(power)
    CS =contourf(r_bins[X], dr_bins[Y], Z , 10)#, vmin=0.0 ,vmax=3.6)
    cb = colorbar(format='%.3f')
    cb.set_label(r'$\bf fraction^{%.2f}$'%(power),fontsize=16)

    #axhline(y=400.,color='r',linestyle='--')
    #axhline(y=-400.,color='r',linestyle='--')
    axhline(y=0.,color='w',linestyle='--')


    ylim(Dradmin,Dradmax)
    xlim(5,rmax)


    if save==1:
      savefig(savepath+'dr_rinit_tot.png')




#====================================================================================



#==========================================================================
# plot of delta_R vs Rinit  for all stars

def plot_dr_rinit_tot_dLtot():
    

    figure()
    title(r'$\bf \Delta R\ vs\ R_{init.}\ -\ total\ stars$',fontsize=18)

    ylabel(r'$ \Delta R\ - \Delta L/220.\ (kpc)$',fontsize=16)
    xlabel(r'$R_{init.}\ (kpc)$',fontsize=16)


    #nb = 50
    dr_bins=zeros( nDrad + 1)
    for i in range(nDrad+1):
        dr_bins[i] =  Dradmin + (Dradmax - Dradmin)*float(i)/float(nDrad) 

    #nR = 50
    r_bins = zeros(nrad+1)
    for i in range(nrad+1):
        r_bins[i] =  rmin + (rmax - rmin)*float(i)/float(nrad) 


    bins1 = binning_2D(r_bins,ritot,5,1,dr_bins,delta_rtot,5,0,1.)
    #bins2 = binning_2D(r_bins,ritot,5,1,dr_bins,delta_gc2tot,5,0,1.)
    bins2 = binning_2D(r_bins,ritot,5,1,dr_bins,delta_angmtot/220.,5,0,1.)

    for i in range(len(bins1[:,0])):
        bins1[i,:]=bins1[i,:]/sum(bins1[i,:])
        bins2[i,:]=bins2[i,:]/sum(bins2[i,:])



    bins = bins1 - bins2


    #for i in range(len(bins[:,0])):
    #    bins[i,:]=bins[i,:]/sum(bins[i,:])

    def bbb(x,y):
        return bins[x,y]

    I=range(len(r_bins))
    J=range(len(dr_bins)) 
    X,Y = meshgrid(I, J)


    Z  = (bbb(X, Y))**(power)
    CS =contourf(r_bins[X], dr_bins[Y] , Z , 10 )#, vmin=0.0 , vmax=0.1 )
    cb = colorbar(format='%.3f')
    cb.set_label(r'$\bf fraction^{%.2f}$'%(power),fontsize=16)

    #axhline(y=400.,color='r',linestyle='--')
    #axhline(y=-400.,color='r',linestyle='--')
    axhline(y=0.,color='w',linestyle='--')


    ylim(Dradmin,Dradmax)
    xlim(5,rmax)


    if save==1:
      savefig(savepath+'dr_rinit_tot_dLtot.png')




#====================================================================================



# plot of initial tangential vel vs radius  for stars detected as migrating ones

def plot_tang_veli_r():

    figure()
    title(r'$\bf init.\ tangential\ vel. -\ migr\ stars$',fontsize=18)

    ylabel(r'$velocity\ (km/s)$',fontsize=16)
    xlabel(r'$radius (kpc)$',fontsize=16)

    
    r_bins=zeros( nrad +1)
    vel_bins=zeros( nvel +1)

    for i in range(nrad+1):
        r_bins[i] = rmin + (rmax - rmin)/float(nrad)*float(i) 

    for i in range(nvel+1):
        vel_bins[i] = velmin_t+ (velmax_t - velmin_t)/float(nvel)*float(i) 


    bins = binning_2D(r_bins,ri,5,0,vel_bins,tveli,5,2,1)

    def bbb(x,y):
        return bins[x,y]


    I=range(len(r_bins))
    J=range(len(vel_bins)) 
    X,Y = meshgrid(I, J)


    Z  = (bbb(X, Y))**(power)
    CS =contourf(r_bins[X], vel_bins[Y], Z , 10)#, vmin=0.0 ,vmax=0.15)
    cb = colorbar(format='%.3f')
    cb.set_label(r'$\bf number^{%.2f}$'%(power),fontsize=16)




    v_mean=zeros(len(r_bins))

    for i in range(len(r_bins)):
        mm=0.
        for j in range(len(vel_bins)-1):
            v_mean[i] = v_mean[i] + bins[i,j]*(vel_bins[j]+vel_bins[j+1])/2.
            mm=mm+bins[i,j]
        v_mean[i]=v_mean[i]/mm


    plot(r_bins,v_mean,'w--',linewidth=2,label='average velocity')

    legend(loc='best')


    #xlim(4.5,15)



    if save==1:
      savefig(savepath+'i-Tvel-migr.png')

#====================================================================================

# plot of final tangential vel vs radius  for stars detected as migrating ones

def plot_tang_velf_r():

    figure()
    title(r'$\bf final.\ tangential\ vel. -\ migr\ stars$',fontsize=18)

    ylabel(r'$velocity\ (km/s)$',fontsize=16)
    xlabel(r'$radius (kpc)$',fontsize=16)

    
    r_bins=zeros( nrad +1)
    vel_bins=zeros( nvel +1)

    for i in range(nrad+1):
        r_bins[i] = rmin + (rmax - rmin)/float(nrad)*float(i) 

    for i in range(nvel+1):
        vel_bins[i] = velmin_t+ (velmax_t - velmin_t)/float(nvel)*float(i) 


    bins = binning_2D(r_bins,rf,5,1,vel_bins,tvelf,5,2,1)

    def bbb(x,y):
        return bins[x,y]


    I=range(len(r_bins))
    J=range(len(vel_bins)) 
    X,Y = meshgrid(I, J)


    Z  = (bbb(X, Y))**(power)
    CS =contourf(r_bins[X], vel_bins[Y], Z , 10)#, vmin=0.0 ,vmax=0.15)
    cb = colorbar(format='%.3f')
    cb.set_label(r'$\bf number^{%.2f}$'%(power),fontsize=16)




    v_mean=zeros(len(r_bins))

    for i in range(len(r_bins)):
        mm=0.
        for j in range(len(vel_bins)-1):
            v_mean[i] = v_mean[i] + bins[i,j]*(vel_bins[j]+vel_bins[j+1])/2.
            mm=mm+bins[i,j]
        v_mean[i]=v_mean[i]/mm


    plot(r_bins,v_mean,'w--',linewidth=2,label='average velocity')

    legend(loc='best')


    if save==1:
      savefig(savepath+'f-Tvel-migr.png')
 

#====================================================================================

# plot of initial radial vel vs radius  for stars detected as migrating ones

def plot_rad_veli_r():

    figure()
    title(r'$\bf init.\ radial\ vel. -\ migr\ stars$',fontsize=18)

    ylabel(r'$velocity\ (km/s)$',fontsize=16)
    xlabel(r'$radius (kpc)$',fontsize=16)

    
    r_bins=zeros( nrad +1)
    vel_bins=zeros( nvel +1)

    for i in range(nrad+1):
        r_bins[i] = rmin + (rmax - rmin)/float(nrad)*float(i) 

    for i in range(nvel+1):
        vel_bins[i] = velmin_r+ (velmax_r - velmin_r)/float(nvel)*float(i) 



    bins = binning_2D(r_bins,ri,5,1,vel_bins,rveli,5,1,1)

    def bbb(x,y):
        return bins[x,y]


    I=range(len(r_bins))
    J=range(len(vel_bins)) 
    X,Y = meshgrid(I, J)


    Z  = (bbb(X, Y))**(power)
    CS =contourf(r_bins[X], vel_bins[Y], Z , 10)#, vmin=0.0 ,vmax=0.15)
    cb = colorbar(format='%.3f')
    cb.set_label(r'$\bf number^{%.2f}$'%(power),fontsize=16)


    v_mean=zeros(len(r_bins))

    for i in range(len(r_bins)):
        mm=0.
        for j in range(len(vel_bins)-1):
            v_mean[i] = v_mean[i] + bins[i,j]*(vel_bins[j]+vel_bins[j+1])/2.
            mm=mm+bins[i,j]
        v_mean[i]=v_mean[i]/mm


    plot(r_bins,v_mean,'w--',linewidth=2,label='average velocity')

    legend(loc='best')



    if save==1:
      savefig(savepath+'i-Rvel-migr.png')

#====================================================================================

# plot of final radial vel vs radius  for stars detected as migrating ones

def plot_rad_velf_r():

    figure()
    title(r'$\bf final.\ radial\ vel. -\ migr\ stars$',fontsize=18)

    ylabel(r'$velocity\ (km/s)$',fontsize=16)
    xlabel(r'$radius (kpc)$',fontsize=16)

    
    r_bins=zeros( nrad +1)
    vel_bins=zeros( nvel +1)

    for i in range(nrad+1):
        r_bins[i] = rmin + (rmax - rmin)/float(nrad)*float(i) 

    for i in range(nvel+1):
        vel_bins[i] = velmin_r+ (velmax_r - velmin_r)/float(nvel)*float(i) 



    bins = binning_2D(r_bins,rf,5,1,vel_bins,rvelf,5,1,1)

    def bbb(x,y):
        return bins[x,y]


    I=range(len(r_bins))
    J=range(len(vel_bins)) 
    X,Y = meshgrid(I, J)


    Z  = (bbb(X, Y))**(power)
    CS =contourf(r_bins[X], vel_bins[Y], Z , 10)#, vmin=0.0 ,vmax=0.15)
    cb = colorbar(format='%.3f')
    cb.set_label(r'$\bf number^{%.2f}$'%(power),fontsize=16)




    v_mean=zeros(len(r_bins))

    for i in range(len(r_bins)):
        mm=0.
        for j in range(len(vel_bins)-1):
            v_mean[i] = v_mean[i] + bins[i,j]*(vel_bins[j]+vel_bins[j+1])/2.
            mm=mm+bins[i,j]
        v_mean[i]=v_mean[i]/mm


    plot(r_bins,v_mean,'w--',linewidth=2,label='average velocity')

    legend(loc='best')


    if save==1:
      savefig(savepath+'f-Rvel-migr.png')

#====================================================================================

# plot of initial tangential vel vs radius  for all stars
def plot_tang_velitot_r():

    figure()
    title(r'$\bf init.\ tangential\ vel. -\ all\ stars$',fontsize=18)

    ylabel(r'$velocity\ (km/s)$',fontsize=16)
    xlabel(r'$radius (kpc)$',fontsize=16)

    
    r_bins=zeros( nrad +1)
    vel_bins=zeros( nvel +1)

    for i in range(nrad+1):
        r_bins[i] = rmin + (rmax - rmin)/float(nrad)*float(i) 

    for i in range(nvel+1):
        vel_bins[i] = velmin_t+ (velmax_t - velmin_t)/float(nvel)*float(i) 



    bins = binning_2D(r_bins,ritot,5,0,vel_bins,tvelitot,5,2,1)

    def bbb(x,y):
        return bins[x,y]


    I=range(len(r_bins))
    J=range(len(vel_bins)) 
    X,Y = meshgrid(I, J)


    Z  = (bbb(X, Y))**(power)
    CS =contourf(r_bins[X], vel_bins[Y], Z , 10)#, vmin=0.0 ,vmax=0.15)
    cb = colorbar(format='%.3f')
    cb.set_label(r'$\bf number^{%.2f}$'%(power),fontsize=16)



    v_mean=zeros(len(r_bins))

    for i in range(len(r_bins)):
        mm=0.
        for j in range(len(vel_bins)-1):
            v_mean[i] = v_mean[i] + bins[i,j]*(vel_bins[j]+vel_bins[j+1])/2.
            mm=mm+bins[i,j]
        v_mean[i]=v_mean[i]/mm


    plot(r_bins,v_mean,'w--',linewidth=2,label='average velocity')

    legend(loc='best')



    if save==1:
      savefig(savepath+'i-Tvel-all.png')

#====================================================================================

# plot of final tangential vel vs radius  for all stars

def plot_tang_velftot_r():

    figure()
    title(r'$\bf final.\ tangential\ vel. -\ all\ stars$',fontsize=18)

    ylabel(r'$velocity\ (km/s)$',fontsize=16)
    xlabel(r'$radius (kpc)$',fontsize=16)

    
    r_bins=zeros( nrad +1)
    vel_bins=zeros( nvel +1)

    for i in range(nrad+1):
        r_bins[i] = rmin + (rmax - rmin)/float(nrad)*float(i) 

    for i in range(nvel+1):
        vel_bins[i] = velmin_t+ (velmax_t - velmin_t)/float(nvel)*float(i) 



    bins = binning_2D(r_bins,rftot,5,1,vel_bins,tvelftot,5,2,1)

    def bbb(x,y):
        return bins[x,y]


    I=range(len(r_bins))
    J=range(len(vel_bins)) 
    X,Y = meshgrid(I, J)


    Z  = (bbb(X, Y))**(power)
    CS =contourf(r_bins[X], vel_bins[Y], Z , 10)#, vmin=0.0 ,vmax=0.15)
    cb = colorbar(format='%.3f')
    cb.set_label(r'$\bf number^{%.2f}$'%(power),fontsize=16)




    v_mean=zeros(len(r_bins))

    for i in range(len(r_bins)):
        mm=0.
        for j in range(len(vel_bins)-1):
            v_mean[i] = v_mean[i] + bins[i,j]*(vel_bins[j]+vel_bins[j+1])/2.
            mm=mm+bins[i,j]
        v_mean[i]=v_mean[i]/mm


    plot(r_bins,v_mean,'w--',linewidth=2,label='average velocity')

    legend(loc='best')


    if save==1:
      savefig(savepath+'f-Tvel-all.png')
 

#====================================================================================

# plot of initial radial vel vs radius  for all stars

def plot_rad_velitot_r():

    figure()
    title(r'$\bf init.\ radial\ vel. -\ all\ stars$',fontsize=18)

    ylabel(r'$velocity\ (km/s)$',fontsize=16)
    xlabel(r'$radius (kpc)$',fontsize=16)

    
    r_bins=zeros( nrad +1)
    vel_bins=zeros( nvel +1)

    for i in range(nrad+1):
        r_bins[i] = rmin + (rmax - rmin)/float(nrad)*float(i) 

    for i in range(nvel+1):
        vel_bins[i] = velmin_r+ (velmax_r - velmin_r)/float(nvel)*float(i) 


    bins = binning_2D(r_bins,ritot,5,1,vel_bins,rvelitot,5,1,1)

    def bbb(x,y):
        return bins[x,y]


    I=range(len(r_bins))
    J=range(len(vel_bins)) 
    X,Y = meshgrid(I, J)


    Z  = (bbb(X, Y))**(power)
    CS =contourf(r_bins[X], vel_bins[Y], Z , 10)#, vmin=0.0 ,vmax=0.15)
    cb = colorbar(format='%.3f')
    cb.set_label(r'$\bf number^{%.2f}$'%(power),fontsize=16)


    v_mean=zeros(len(r_bins))

    for i in range(len(r_bins)):
        mm=0.
        for j in range(len(vel_bins)-1):
            v_mean[i] = v_mean[i] + bins[i,j]*(vel_bins[j]+vel_bins[j+1])/2.
            mm=mm+bins[i,j]
        v_mean[i]=v_mean[i]/mm


    plot(r_bins,v_mean,'w--',linewidth=2,label='average velocity')

    legend(loc='best')


    if save==1:
      savefig(savepath+'i-Rvel-all.png')

#====================================================================================

# plot of final radial vel vs radius  for all stars

def plot_rad_velftot_r():

    figure()
    title(r'$\bf final.\ radial\ vel. -\ all\ stars$',fontsize=18)

    ylabel(r'$velocity\ (km/s)$',fontsize=16)
    xlabel(r'$radius (kpc)$',fontsize=16)

    
    r_bins=zeros( nrad +1)
    vel_bins=zeros( nvel +1)

    for i in range(nrad+1):
        r_bins[i] = rmin + (rmax - rmin)/float(nrad)*float(i) 

    for i in range(nvel+1):
        vel_bins[i] = velmin_r+ (velmax_r - velmin_r)/float(nvel)*float(i) 


    bins = binning_2D(r_bins,rftot,5,1,vel_bins,rvelftot,5,1,1)

    def bbb(x,y):
        return bins[x,y]


    I=range(len(r_bins))
    J=range(len(vel_bins)) 
    X,Y = meshgrid(I, J)





    Z  = (bbb(X, Y))**(power)
    CS =contourf(r_bins[X], vel_bins[Y], Z , 10)#, vmin=0.0 ,vmax=0.15)
    cb = colorbar(format='%.3f')
    cb.set_label(r'$\bf number^{%.2f}$'%(power),fontsize=16)




    v_mean=zeros(len(r_bins))

    for i in range(len(r_bins)):
        mm=0.
        for j in range(len(vel_bins)-1):
            v_mean[i] = v_mean[i] + bins[i,j]*(vel_bins[j]+vel_bins[j+1])/2.
            mm=mm+bins[i,j]
        v_mean[i]=v_mean[i]/mm


    plot(r_bins,v_mean,'w--',linewidth=2,label='average velocity')

    legend(loc='best')


    if save==1:
      savefig(savepath+'f-Rvel-all.png')



#==========================================================================
# plot number vs Radius  for all stars

def plot_num_rad_tot():
    

    figure()
    title(r'$\bf Number\ vs\ Radius\ -\ total\ stars$',fontsize=18)

    ylabel(r'$\bf Number $',fontsize=16)
    xlabel(r'$\bf Radius\ (kpc) $',fontsize=16)



    r_bins=zeros( nrad +1 )

    for i in range(nrad+1):
        r_bins[i] = rmin + (rmax - rmin)/float(nrad)*float(i) 


    bins1 = binning_1D(r_bins,ritot,5,1,1.)

    bins2 = binning_1D(r_bins,rftot,5,1,1.)



    plot( r_bins,bins1,linewidth=2,label='init')
    plot( r_bins,bins2,linewidth=2,label='fin' )


    grid('on')
    legend(loc='best')


    if save==1:
      savefig(savepath+'num_rad_tot.png')

#==========================================================================
# plot number vs Radius  for all stars

def plot_i_f_pos_tot():
    

    figure(figsize=(16,6))
    subplot(121)
    title(r'$\bf initial\ -\ total\ stars$',fontsize=18)

    ylabel(r'$\bf y\ (kpc) $',fontsize=16)
    xlabel(r'$\bf x\ (kpc) $',fontsize=16)



    nbins = 100
    rr = zeros( nbins +1)
    for i in range(nbins+1):
        rr[i] = -15. + 30./float(nbins)*float(i) 


    binsi = binning_2D(rr,positot[:,0],5,1,rr,positot[:,1],5,1,1.)



    def func(x,y):
      return binsi[x,y]

    I=range(len(rr))
    J=range(len(rr))
    X,Y = meshgrid(I, J)


    Z = (func(X, Y))**(1)
    pcolor(rr[X], rr[Y], Z)#, vmin=0.0 ,v)
    cb = colorbar(format='%.3f')
    #cb.set_label(r'$\bf number $',fontsize=16)


    subplot(122)
    title(r'$\bf final\ -\ total\ stars$',fontsize=18)

    ylabel(r'$\bf y\ (kpc) $',fontsize=16)
    xlabel(r'$\bf x\ (kpc) $',fontsize=16)

    binsf = binning_2D(rr,posftot[:,0],5,1,rr,posftot[:,1],5,1,1.)



    def func(x,y):
      return binsf[x,y]

    I=range(len(rr))
    J=range(len(rr))
    X,Y = meshgrid(I, J)


    Z = (func(X, Y))**(1)
    pcolor(rr[X], rr[Y], Z)#, vmin=0.0 ,v)
    cb = colorbar(format='%.3f')
    cb.set_label(r'$\bf number $',fontsize=16)



    if save==1:
      savefig(savepath+'i_f_pos_tot.png')



#==========================================================================
# plot number vs Radius  for migrating stars

def plot_i_f_pos_migr():
    

    figure(figsize=(16,6))
    subplot(121)
    title(r'$\bf initial\ -\ migr\ stars$',fontsize=18)

    ylabel(r'$\bf y\ (kpc) $',fontsize=16)
    xlabel(r'$\bf x\ (kpc) $',fontsize=16)



    nbins = 100
    rr = zeros( nbins +1)
    for i in range(nbins+1):
        rr[i] = -15. + 30./float(nbins)*float(i) 


    binsi = binning_2D(rr,posinit[:,0],5,1,rr,posinit[:,1],5,1,1.)



    def func(x,y):
      return binsi[x,y]

    I=range(len(rr))
    J=range(len(rr))
    X,Y = meshgrid(I, J)


    Z = (func(X, Y))**(1)
    pcolor(rr[X], rr[Y], Z)#, vmin=0.0 ,v)
    cb = colorbar(format='%.3f')
    #cb.set_label(r'$\bf number $',fontsize=16)


    subplot(122)
    title(r'$\bf final\ -\ migr\ stars$',fontsize=18)

    ylabel(r'$\bf y\ (kpc) $',fontsize=16)
    xlabel(r'$\bf x\ (kpc) $',fontsize=16)

    binsf = binning_2D(rr,posfin[:,0],5,1,rr,posfin[:,1],5,1,1.)



    def func(x,y):
      return binsf[x,y]

    I=range(len(rr))
    J=range(len(rr))
    X,Y = meshgrid(I, J)


    Z = (func(X, Y))**(1)
    pcolor(rr[X], rr[Y], Z,edgecolors = 'None')#, vmin=0.0 ,v)
    cb = colorbar(format='%.3f')
    cb.set_label(r'$\bf number $',fontsize=16)



    if save==1:
      savefig(savepath+'i_f_pos_migr.png')


#=====================================================================================================


def plot_ri_rf_tot():
    figure()  # Rinit vs Rfin for total stars 


    title('R(birth-time) vs. R(10 Gyr) - tot stars')#,fontsize=12)

    ylabel('R initial (kpc)')
    xlabel('R final (kpc)')



    r_bins=zeros( nrad +1)

    for i in range(nrad+1):
        r_bins[i] = rmin + (rmax - rmin)/float(nrad)*float(i) 



    n_bins = binning_2D(r_bins,rftot,5,1,r_bins,ritot,5,1,1.)

    def bbb(x,y):
        return n_bins[x,y]/n_bins[:,y].sum()*100.

    I=range(len(r_bins))
    J=range(len(r_bins))
    X,Y = meshgrid(I, J)


    Z  = (bbb(X, Y))#**(power)
    CS =pcolor(r_bins[X], r_bins[Y], Z,   norm=LogNorm(vmin=0.00001, vmax=Z.max()) )#, vmin=0.0 ,vmax=0.15)
    cb = colorbar()
    #cb.set_label(r'$\bf number \ of \  stars/pc^{2}$',fontsize=fsizel)
    cb.set_label(r'$\bf \% \ of \ tot \ mass \ at \ R_{fin}$',fontsize=fsizel)


    plot(r_bins,r_bins,color='k')


    r_mean=zeros(len(r_bins))

    for i in range(len(r_bins)):
        mm=0.
        for j in range(len(r_bins)-1):
            r_mean[i] = r_mean[i] + n_bins[i,j]*(r_bins[j]+r_bins[j+1])/2.
            mm=mm+n_bins[i,j]
        r_mean[i]=r_mean[i]/mm


    plot(r_bins,r_mean,'k',linewidth=2,label='average birth radius')


    legend(loc='best')


    #ylim(0.,18.)
    #xlim(0.,18.)


    """
    open(path+'mean_rad.res','w').close()

    f = open(path+'mean_rad.res','w')

    print >>f, '# Date (m/y) 04/2012'
    print >>f, '# mean R initial (kpc), for each R final(kpc) at final time (10 Gyr)'
    print >>f, '# results from GADGET simulation mwg001'
    print >>f, ' '
    print >>f, '# Rf(kpc)  ,  Ri_mean(kpc)'
 


    for i in range(len(rad_bins)):
        print >> f, '%10.3f %10.3f' % (rad_bins[i],r_mean[i])


    f.close()
    """



    if save==1:
        savefig(home+'works/res-works/ri_rf_tot.png')
#============================================================================================



def plot_ri_rf_migr():
    figure()  # Rinit vs Rfin for total stars 


    title('R(birth-time) vs. R(10 Gyr) - migr stars')#,fontsize=12)

    ylabel('R initial (kpc)')
    xlabel('R final (kpc)')



    r_bins=zeros( nrad +1)

    for i in range(nrad+1):
        r_bins[i] = rmin + (rmax - rmin)/float(nrad)*float(i) 



    n_bins = binning_2D(r_bins,rf,5,1,r_bins,ri,5,1,1.)

    def bbb(x,y):
        return n_bins[x,y]/n_bins[:,y].sum()*100.

    I=range(len(r_bins))
    J=range(len(r_bins))
    X,Y = meshgrid(I, J)


    Z  = (bbb(X, Y))#**(power)
    CS =pcolor(r_bins[X], r_bins[Y], Z,   norm=LogNorm(vmin=0.00001, vmax=Z.max()) )#, vmin=0.0 ,vmax=0.15)
    cb = colorbar()
    #cb.set_label(r'$\bf number \ of \  stars/pc^{2}$',fontsize=fsizel)
    cb.set_label(r'$\bf \% \ of \ tot \ mass \ at \ R_{fin}$',fontsize=fsizel)


    plot(r_bins,r_bins,color='k')


    r_mean=zeros(len(r_bins))

    for i in range(len(r_bins)):
        mm=0.
        for j in range(len(r_bins)-1):
            r_mean[i] = r_mean[i] + n_bins[i,j]*(r_bins[j]+r_bins[j+1])/2.
            mm=mm+n_bins[i,j]
        r_mean[i]=r_mean[i]/mm


    plot(r_bins,r_mean,'k',linewidth=2,label='average birth radius')


    legend(loc='best')


    #ylim(0.,18.)
    #xlim(0.,18.)





    if save==1:
        savefig(home+'works/res-works/ri_rf_migr.png')
#============================================================================================


def plot_ri_rf_nomigr():
    figure()  # Rinit vs Rfin for non migrating stars 


    title('R(birth-time) vs. R(10 Gyr) - no migr stars')#,fontsize=12)

    ylabel('R initial (kpc)')
    xlabel('R final (kpc)')



    r_bins=zeros( nrad +1)

    for i in range(nrad+1):
        r_bins[i] = rmin + (rmax - rmin)/float(nrad)*float(i) 



    n_binstot = binning_2D(r_bins,rftot,5,1,r_bins,ritot,5,1,1.)

    n_binsmigr = binning_2D(r_bins,rf,5,1,r_bins,ri,5,1,1.)

    n_bins = n_binstot - n_binsmigr

    def bbb(x,y):
        return n_bins[x,y]/n_bins[:,y].sum()*100.

    I=range(len(r_bins))
    J=range(len(r_bins))
    X,Y = meshgrid(I, J)


    Z  = (bbb(X, Y))#**(power)
    CS =pcolor(r_bins[X], r_bins[Y], Z,   norm=LogNorm(vmin=0.00001, vmax=Z.max()) )#, vmin=0.0 ,vmax=0.15)
    cb = colorbar()
    #cb.set_label(r'$\bf number \ of \  stars/pc^{2}$',fontsize=fsizel)
    cb.set_label(r'$\bf \% \ of \ tot \ mass \ at \ R_{fin}$',fontsize=fsizel)


    plot(r_bins,r_bins,color='k')


    r_mean=zeros(len(r_bins))

    for i in range(len(r_bins)):
        mm=0.
        for j in range(len(r_bins)-1):
            r_mean[i] = r_mean[i] + n_bins[i,j]*(r_bins[j]+r_bins[j+1])/2.
            mm=mm+n_bins[i,j]
        r_mean[i]=r_mean[i]/mm


    plot(r_bins,r_mean,'k',linewidth=2,label='average birth radius')


    legend(loc='best')


    #ylim(0.,18.)
    #xlim(0.,18.)





    if save==1:
        savefig(home+'works/res-works/ri_rf_nomigr.png')
#============================================================================================


def plot_angmi_angmf_migr():
    figure()  # Linit vs Lfin for migr stars 


    title('Linit vs. Lfin - migr stars')#,fontsize=12)

    ylabel('L initial (kpc2/Gyr)')
    xlabel('L final (kpc2/Gyr)')



    Lbins=zeros( nL +1)

    for i in range(nL+1):
        Lbins[i] = Lmin + (Lmax - Lmin)/float(nL)*float(i) 




    n_bins = binning_2D(Lbins,angmfin,5,3,Lbins,angminit,5,3,1.)

    def bbb(x,y):
        return n_bins[x,y]/n_bins[:,y].sum()*100.

    I=range(len(Lbins))
    J=range(len(Lbins))
    X,Y = meshgrid(I, J)


    Z  = (bbb(X, Y))#**(power)
    CS =pcolor(Lbins[X], Lbins[Y], Z,   norm=LogNorm(vmin=0.00001, vmax=Z.max()) )#, vmin=0.0 ,vmax=0.15)
    cb = colorbar()
    #cb.set_label(r'$\bf number \ of \  stars/pc^{2}$',fontsize=fsizel)
    cb.set_label(r'$\bf \% \ of \ tot \ mass \ at \ R_{fin}$',fontsize=fsizel)


    plot(Lbins,Lbins,color='k')


    r_mean=zeros(len(Lbins))

    for i in range(len(Lbins)):
        mm=0.
        for j in range(len(Lbins)-1):
            r_mean[i] = r_mean[i] + n_bins[i,j]*(Lbins[j]+Lbins[j+1])/2.
            mm=mm+n_bins[i,j]
        r_mean[i]=r_mean[i]/mm


    plot(Lbins,r_mean,'k',linewidth=2,label='average birth radius')


    legend(loc='best')


    #ylim(0.,18.)
    #xlim(0.,18.)




    if save==1:
        savefig(home+'works/res-works/angmi_angmf_migr.png')
#============================================================================================





def plot_angmi_angmf_tot():
    figure()  # Linit vs Lfin for tot stars 


    title('Linit vs. Lfin - tot stars')#,fontsize=12)

    ylabel('L initial (kpc2/Gyr)')
    xlabel('L final (kpc2/Gyr)')



    Lbins=zeros( nL +1)

    for i in range(nL+1):
        Lbins[i] = Lmin + (Lmax - Lmin)/float(nL)*float(i) 




    n_bins = binning_2D(Lbins,angmitot,5,3,Lbins,angmftot,5,3,1.)

    def bbb(x,y):
        return n_bins[x,y]/n_bins[:,y].sum()*100.

    I=range(len(Lbins))
    J=range(len(Lbins))
    X,Y = meshgrid(I, J)


    Z  = (bbb(X, Y))#**(power)
    CS =pcolor(Lbins[X], Lbins[Y], Z,  norm=LogNorm(vmin=0.00001, vmax=Z.max()) )#, vmin=0.0 ,vmax=0.15)
    cb = colorbar()
    #cb.set_label(r'$\bf number \ of \  stars/pc^{2}$',fontsize=fsizel)
    cb.set_label(r'$\bf \% \ of \ tot \ mass \ at \ R_{fin}$',fontsize=fsizel)


    plot(Lbins,Lbins,color='k')


    r_mean=zeros(len(Lbins))

    for i in range(len(Lbins)):
        mm=0.
        for j in range(len(Lbins)-1):
            r_mean[i] = r_mean[i] + n_bins[i,j]*(Lbins[j]+Lbins[j+1])/2.
            mm=mm+n_bins[i,j]
        r_mean[i]=r_mean[i]/mm


    plot(Lbins,r_mean,'k',linewidth=2,label='average birth radius')


    legend(loc='best')


    #ylim(0.,18.)
    #xlim(0.,18.)




    if save==1:
        savefig(home+'works/res-works/angmi_angmf_tot.png')
#===========================================================================================








def plot_angmi_angmf_nomigr():
    figure()  # Linit vs Lfin for non migrating stars 


    title('Linit vs. Lfin - non migr stars')#,fontsize=12)

    ylabel('L initial (kpc2/Gyr)')
    xlabel('L final (kpc2/Gyr)')



    Lbins=zeros( nL +1)

    for i in range(nL+1):
        Lbins[i] = Lmin + (Lmax - Lmin)/float(nL)*float(i) 




    n_binstot  = binning_2D(Lbins,angmitot,5,3,Lbins,angmftot,5,3,1.)

    n_binsmigr = binning_2D(Lbins,angminit,5,3,Lbins,angmfin,5,3,1.)

    n_bins     = n_binstot - n_binsmigr


    def bbb(x,y):
        return n_bins[x,y]/n_bins[:,y].sum()*100.

    I=range(len(Lbins))
    J=range(len(Lbins))
    X,Y = meshgrid(I, J)


    Z  = (bbb(X, Y))**(1.)
    CS =pcolor(Lbins[X], Lbins[Y], Z, norm=LogNorm(vmin=0.00001, vmax=5e-3))#Z.max()) )#, vmin=0.0 ,vmax=0.15)
    cb = colorbar()
    #cb.set_label(r'$\bf number \ of \  stars/pc^{2}$',fontsize=fsizel)
    cb.set_label(r'$\bf \% \ of \ tot \ mass \ at \ R_{fin}$',fontsize=fsizel)


    plot(Lbins,Lbins,color='k')


    r_mean=zeros(len(Lbins))

    for i in range(len(Lbins)):
        mm=0.
        for j in range(len(Lbins)-1):
            r_mean[i] = r_mean[i] + n_bins[i,j]*(Lbins[j]+Lbins[j+1])/2.
            mm=mm+n_bins[i,j]
        r_mean[i]=r_mean[i]/mm


    plot(Lbins,r_mean,'k',linewidth=2,label='average birth radius')


    legend(loc='best')


    #ylim(0.,18.)
    #xlim(0.,18.)




    if save==1:
        savefig(home+'works/res-works/angmi_angmf_nomigr.png')
#===========================================================================================


def plot_angmi_angmf_nomigr_fold():
    figure()  # Linit vs Lfin for non migrating stars 


    title('Linit vs. Lfin - non migr stars')#,fontsize=12)

    ylabel('L initial (kpc2/Gyr)')
    xlabel('L final (kpc2/Gyr)')



    Lbins=zeros( nL +1)

    for i in range(nL+1):
        Lbins[i] = Lmin + (Lmax - Lmin)/float(nL)*float(i) 




    n_binstot  = binning_2D(Lbins,angmitot,5,3,Lbins,angmftot,5,3,1.)

    n_binsmigr = binning_2D(Lbins,angminit,5,3,Lbins,angmfin,5,3,1.)

    n_bins     = n_binstot - n_binsmigr


    def bbb(x,y):
        return n_bins[x,y]/n_bins[:,y].sum()*100.


    fold_bins = 0.*n_bins


    for i in range( len( Lbins ) ):
        for j in range( len( Lbins ) ):
            fold_bins[i,j] = fold_bins[i,j] + n_bins 


    I=range(len(Lbins))
    J=range(len(Lbins))
    X,Y = meshgrid(I, J)


    Z  = (bbb(X, Y))**(1.)
    CS =pcolor(Lbins[X], Lbins[Y], Z, norm=LogNorm(vmin=0.00001, vmax=Z.max()) )#, vmin=0.0 ,vmax=0.15)
    cb = colorbar()
    #cb.set_label(r'$\bf number \ of \  stars/pc^{2}$',fontsize=fsizel)
    cb.set_label(r'$\bf \% \ of \ tot \ mass \ at \ R_{fin}$',fontsize=fsizel)


    plot(Lbins,Lbins,color='k')


    r_mean=zeros(len(Lbins))

    for i in range(len(Lbins)):
        mm=0.
        for j in range(len(Lbins)-1):
            r_mean[i] = r_mean[i] + n_bins[i,j]*(Lbins[j]+Lbins[j+1])/2.
            mm=mm+n_bins[i,j]
        r_mean[i]=r_mean[i]/mm


    plot(Lbins,r_mean,'k',linewidth=2,label='average birth radius')


    legend(loc='best')


    #ylim(0.,18.)
    #xlim(0.,18.)




    if save==1:
        savefig(home+'works/res-works/angmi_angmf_nomigr_fold.png')













