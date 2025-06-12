# -*- coding: utf-8 -*-
"""
This is a python driver for the simple_letkf_wloc fortran routine.
This code provides a way to run simple assimilation experiments with
realistic priors. 

This version uses tempering (iterative) LETKF to compute the update
 
@author:
"""
import sys
sys.path.append('common_python/common_functions/')
sys.path.append('common_python/common_modules/')
sys.path.append('common_python/common_letkf/')

import numpy as np
import os

from cletkf_wloc      import common_da        as cda


root_data_path='/home/jorge.gacitua/experimentos/CrossSection_Assimilation/Data'

qg_index = 0         #Variable index corresponding to QG
qr_index = 1         #Variable index corresponding to QR
qs_index = 2         #Variable index corresponding to QS
tt_index = 3         #Variable index corresponding to TEMP
pp_index = 4         #Variable index corresponding to PRES

case = 'mean' #Case to be considered. Options are: 'mean', 'above' and 'below'
if case == 'mean':
#### closest to mean
    truth_ens_member = 22 #Which ensemble member will be considered as the true state.
    obs_loc_x = [30]     #List of x coordinate of observations
    obs_loc_y = [0]      #List of y coordinate of observations
    obs_loc_z = [18]     #List of z coordinate of observations
elif case == 'above':
#### farthest above mean
    truth_ens_member = 9 #Which ensemble member will be considered as the true state.
    obs_loc_x = [30]     #List of x coordinate of observations
    obs_loc_y = [0]      #List of y coordinate of observations
    obs_loc_z = [18]     #List of z coordinate of observations
elif case == 'below':
#### farthest below mean
    truth_ens_member = 19 #Which ensemble member will be considered as the true state.
    obs_loc_x = [30]     #List of x coordinate of observations
    obs_loc_y = [0]      #List of y coordinate of observations
    obs_loc_z = [18]     #List of z coordinate of observations

obs_error = 5.0      #Standard deviation of the observation error.

loc_scales = np.array([5,5,5])  #Localization scales in x,y and z.

range_ntemp  = np.arange(1,4)   #Number of tempering iterations 
range_alpha  = [1,2]            #Parameter controling the size of tempering steps.

#=========================================================
#  READ DATA
#=========================================================

input_ens = np.load("/home/jorge.gacitua/experimentos/CrossSection_Assimilation/Data/ensemble_cross_sections_39-5S.npz")["cross_sections"]

for NTemp in range_ntemp:
    for Alpha in range_alpha:

        #=========================================================
        #  SEPARATE THE TRUE STATE AND THE FORECAST ENSEMBLE
        #=========================================================

        tmp_mask = np.zeros( input_ens.shape[3] ).astype(bool) 
        tmp_mask[truth_ens_member] = True

        true_state = input_ens[:,:,:,tmp_mask,:][:,:,:,0,:]
        xf = input_ens[:,:,:,~tmp_mask,:]

        #Get the size of the forecast ensemble (as will be used in the DA)
        [nx , ny , nz , nbv , nvar] = xf.shape

        #=========================================================
        #  GET THE OBSERVATIONS
        #=========================================================

        nobs = len(obs_loc_x)
        yo = np.zeros( nobs )
        hxf = np.zeros( ( nobs , nbv ) )
        obs_error = obs_error * np.ones( nobs )
        
        for ii in range( nobs )  : 
        
            ox = obs_loc_x[ii]
            oy = obs_loc_y[ii]
            oz = obs_loc_z[ii]
            qr = true_state[ox,oy,oz,qr_index]
            qg = true_state[ox,oy,oz,qg_index]
            qs = true_state[ox,oy,oz,qs_index]
            tt = true_state[ox,oy,oz,tt_index]
            pp = true_state[ox,oy,oz,pp_index]

            yo[ii] = cda.calc_ref( qr , qs , qg , tt , pp )

        #=========================================================
        #  COMPUTE THE TEMPERING STEPS
        #=========================================================

        #NTemp is the number of tempering steps to be performed.
        #Alpha is a slope coefficient. Larger alpha means only a small part of the information
        #will be assimilated in the first step (and the largest part will be assimilated in the last step).

        dt=1.0/float(NTemp+1)
        steps = np.exp( 1.0 * Alpha / np.arange( dt , 1.0-dt/100.0 , dt ) )
        steps = steps / np.sum(steps)
        #Compute normalized pseudo_time tempering steps:
        steps = ( 1.0 / steps ) /  np.sum( 1.0 / steps )

        #=========================================================
        #  ALLOCATE THE SPACE FOR THE ANALYSIS STEPS
        #=========================================================

        xatemp = np.zeros( (nx,ny,nz,nbv,nvar,NTemp+1) )
        xatemp[:,:,:,:,:,0] = xf[:]

        #=========================================================
        #  APPLY OBSERVATION OPERATOR
        #=========================================================

        for it in range(NTemp) :
            for ii in range( nobs ) :
                #Loop over the ensemble members
                for jj in range( nbv ) :

                    qr = xatemp[ox,oy,oz,jj,qr_index,it]
                    qg = xatemp[ox,oy,oz,jj,qg_index,it]
                    qs = xatemp[ox,oy,oz,jj,qs_index,it]
                    tt = xatemp[ox,oy,oz,jj,tt_index,it]
                    pp = xatemp[ox,oy,oz,jj,pp_index,it]

                    hxf[ii,jj] = cda.calc_ref( qr , qs , qg , tt , pp )

            #Get the ensemble mean observation departure.
            dep = yo - np.mean( hxf , 1 ) 

            #=========================================================
            #  COMPUTE THE ANALYSIS UPDATE
            #=========================================================

            #Compute the simple analysis update
            print('Computing the letkf update')


            xatemp = np.asfortranarray( xatemp.astype('float32') )   #Transform data type 
            hxf = np.asfortranarray( hxf.astype('float32') )
            dep = dep.astype('float32')
            obs_error = obs_error.astype('float32')

            obs_error_temp = obs_error * ( 1.0 / steps[it] )

            xatemp[:,:,:,:,:,it+1]=cda.simple_letkf_wloc(nx=nx,ny=ny,nz=nz,
                                    nbv=nbv,nvar=nvar,nobs=nobs,
                                    hxf=hxf,xf=xatemp[:,:,:,:,:,it],
                                    dep=dep,ox=obs_loc_x,
                                    oy=obs_loc_y ,oz=obs_loc_z,
                                    locs=loc_scales,
                                    oerr=obs_error_temp
                                    ).astype('float32')


        #Write the analysis for the update variables
        print('Writing data')
        np.savez_compressed(f'{root_data_path}/output_Ntemp{NTemp}_alpha{Alpha}_8var_Rloc_5_{case}.npz',xatemp=xatemp,xf=xf,yo=yo,
                            hxf=hxf,obs_error=obs_error,obs_error_temp=steps,
                            obs_loc_x=obs_loc_x,
                            obs_loc_y=obs_loc_y,
                            obs_loc_z=obs_loc_z,
                            true_state=true_state)

        print ( "We are done" )