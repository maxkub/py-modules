import os


#====================================================================
#====================================================================


def video(case='', identification='',vcodec='',fps=''):
    command = ('mencoder',
               'mf://'+case+'/'+identification+'*.png',
               '-mf',
               'type=png:w=800:h=600:fps='+fps,   #25
               '-ovc',
               'lavc',
               '-lavcopts',
               'vcodec='+vcodec,   #mpeg4,mjpeg
               '-oac',
               'copy',
               '-o',
               case + '/'+identification+'_video.avi')

    
   
    os.spawnvp(os.P_WAIT, 'mencoder', command)
   
   


#===================================================================
