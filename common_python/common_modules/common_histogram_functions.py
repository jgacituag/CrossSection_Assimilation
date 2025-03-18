import numpy as np
import datetime as dt
import os.path

import sys
sys.path.append('../../common_python/common_functions/')
sys.path.append('../../common_python/common_modules/')
import ctl_reader as ctlr
import binary_io as bio



def read_histogram(filenamehist,filenamemax,filenamemin,nx,ny,nbins,ctl_dict,dtypein) :

   my_hist=dict() #Initialize my_hist dictionary

   #Max_val and Min_val are already dictionaries.
   max_val=kld=ctlr.read_data_grads(filenamemax,ctl_dict,masked=False)
   min_val=kld=ctlr.read_data_grads(filenamemin,ctl_dict,masked=False)

   nz=(np.max( ctl_dict['end_record'] )+1) * nbins #Total number of layers times the number of bins

   undef=np.float32( ctl_dict['undef'] )
   #The total number of records is the total number of variables in the ctl times the number of bins in the histogram.
   nrecords = ( ctl_dict['end_record'] + 1 ) * nbins

   tmp_data=bio.read_data_direct(filenamehist,nx,ny,nz,dtypein=dtypein,undef_in=undef,undef_out=undef,seq_acces=False)  #Read the data.

   print(np.max(tmp_data),np.min(tmp_data))
   #Loop over variables to create the dictionary. 
   ini_record = ctl_dict['ini_record']

   ivar=0

   for my_var in ctl_dict['var_list']   :

     nlevs = int( ctl_dict['var_size'][ivar] )
     if nlevs == 0   :
        nlevs =  int(1)

     my_hist[my_var]=dict()

     my_hist[my_var]['hist']=np.ones((nx,ny,nlevs,nbins)) * undef
 
     my_hist[my_var]['minval']=min_val[my_var]

     my_hist[my_var]['maxval']=max_val[my_var]

     for iz in range(0,nlevs) :

        my_hist[my_var]['hist'][:,:,iz,:]=tmp_data[:,:,(ini_record[ivar]+iz)*nbins:(ini_record[ivar]+iz+1)*nbins]

     ivar=ivar+1

   return my_hist 

def analyze_histogram_fun( my_hist , thresholdmin_input )  :

   output=dict()

   for key in my_hist :

     output[key]=dict()
 
     tmp_hist=my_hist[key]['hist']

     [nx , ny , nz , nbins] = np.shape( tmp_hist )

     enssize=np.max(np.max(np.max(np.sum( tmp_hist , 3 ))))

     tmp_hist = tmp_hist / enssize

     #Smooth histogram

     filter_length=2
 
     tmp_hist_s=tmp_hist

     for ii in range(0,nbins) :

       maxi=ii+filter_length

       mini=ii-filter_length

       if(maxi > nbins) :
            maxi=nbins
       if(mini < 1    ) : 
             mini=1

       tmp_hist_s[:,:,:,ii]=np.nanmean( tmp_hist[:,:,:,mini:maxi] , 3)

     max1val=np.zeros((nx,ny,nz),dtype=int16)
     max1loc=np.zeros((nx,ny,nz),dtype=int16)
     max1mas=np.ones((nx,ny,nz),dtype=bool)

     max2val=np.zeros((nx,ny,nz),dtype=int16)
     max2loc=np.zeros((nx,ny,nz),dtype=int16)
     max2mas=np.zeros((nx,ny,nz),dtype=bool)

     min1mas=np.zeros((nx,ny,nz),dtype=bool)
     min1val=99999*np.ones((nx,ny,nz),dtype=int16)
     min1loc=np.zeros((nx,ny,nz),dtype=int16)

     for i in range(2,nbins) :

        tmpbin=tmp_hist_s[:,:,:,i]

        #Searching for first pdf maximum (usually the only mode of the pdf)

        tmpmask=np.logical_and( tmpbin > max1val , max1mas )

        max1val[ tmpmask ]=tmpbin( tmpmask )

        max1loc[ tmpmask ]=i

        #If pdf values drops more than thresholdmin_input % 

        tmpmask= np.logical_and( ( max1val - tmpbin ) > thresholdmin_input , max1mas )

        min1mas[ tmpmask ]=True    #Start searching for local minima.
        max1mas[ tmpmask ]=False   #Stop searching for first maximum.

        tmpmask= np.logical_and( min1val > tmpbin , min1mas )
        min1val[ tmpmask ]=tmpbin( tmpmask )
        min1loc[ tmpmask ]=i

        #If pdf raises again after first minimum more than thresholdmin_input % then start searching for second maximum.

        tmpmask=np.logical_and( ( tmpbin - min1val ) > thresholdmin_input , min1mas )

        max2mas[ tmpmask ]=True;   #Start searching for second maximum.
        min1mas[ tmpmask ]=False;  #Stop searching for local minima.

        tmpmask=np.logical_and( tmpbin > max2val , max2mas )

        max2val[ tmpmask ]=tmpbin( tmpmask )
        max2loc[ tmpmask ]=i


     min1loc[ max2loc == 0 ]=0   #If we fail to find a second maximum then ignore the local minimum. 

     min1val[ max2loc == 0 ]=0

     output[key]['max1val']=max1val
     output[key]['max1loc']=max1loc
     output[key]['max1mas']=max1mas

     output[key]['max2val']=max2val
     output[key]['max2loc']=max2loc
     output[key]['max2mas']=max2mas

     output[key]['min1mas']=min1mas
     output[key]['min1val']=min1val
     output[key]['min1loc']=min1loc

   return output

       




