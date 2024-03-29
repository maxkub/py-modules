########################################################################
# program: fft.py
# author: Tom Irvine
# Email: tom@vibrationdata.com
# version: 1.9
# date: September 12, 2013
# description:  
#    
#  Determine the Fast Fourier transform of a signal.
#  The file must have two columns: time(sec) & amplitude.
#              
########################################################################

from __future__ import print_function
    
import sys

if sys.version_info[0] == 2:
    print ("Python 2.x")
    import Tkinter as tk
    from tkFileDialog import asksaveasfilename

           
if sys.version_info[0] == 3:
    print ("Python 3.x")    
    import tkinter as tk 
    from tkinter.filedialog import asksaveasfilename    
    



from math import pi,atan2,log
from numpy import zeros,argmax,mean

from scipy.fftpack import fft

import matplotlib.pyplot as plt


########################################################################

class FFT:

    def __init__(self,a,b,imr):
        self.a=a
        self.b=b
        self.imr=imr
        pass
        
    def fft_data(self):
        
#   Truncate to 2**n

        num=len(self.b)

        noct=int(log(num)/log(2))

        num_fft=2**noct

        bb=self.b[0:num_fft]
        
        if(self.imr==1):
            bb=bb-mean(bb)

        dur_fft=self.a[num_fft-1]-self.a[0]

        df=1/dur_fft
      
        
        z =fft(bb)

        nhalf=num_fft/2

        print (" ")
        print (" %d samples used for FFT " %num_fft)
        print ("df = %8.4g Hz" %df)

        zz=zeros(nhalf,'f')
        ff=zeros(nhalf,'f')
        ph=zeros(nhalf,'f')

        freq=zeros(num_fft,'f')

        z/=float(num_fft)

        for k in range(0,int(num_fft)):
            freq[k]=k*df
    
        ff=freq[0:nhalf]
        

    
        for k in range(0,int(nhalf)):    

            if(k > 0):			 
                zz[k]=2.*abs(z[k])
            else:    
                zz[k]= abs(z[k])

            ph[k]=atan2(z.real[k],z.imag[k])
  

        idx = argmax(abs(zz))        
 
        return idx,freq,ff,z,zz,ph,nhalf,df,num_fft    
    

########################################################################


#print (" ")
#print (" Remove mean:  1=yes  2=no ")

#imr = GetInteger2()
          
########################################################################
 
#idx,freq,ff,z,zz,ph,nhalf,df,num_fft=FFT(a,b,imr).fft_data() 

########################################################################


