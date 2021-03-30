#------------------------------------------------------------------------- 
#  
# Shallow Water Model
#
# Danilo Couto de Souza
# Universidade de São Paulo (USP)
# Instituto de Astronomia, Geociências e Ciências Atmosféricas
#-------------------------------------------------------------------------

"""
Created on Tue Dec 22 12:17:12 2020

@author: Danilo
"""

# Import modules
import numpy as np
import time
import matplotlib.pyplot as plt 
import pylab as pl

# Import subroutines
from forcing import (gauss_space,decay)
from grid import (create_grid, coriolis)
from diagnostic_eqs import (fluxU, fluxV, Bournelli, csi)
from filters import (RAW_filter, difusion_x, difusion_y,
                     RA_filter_x, RA_filter_y, RA_filter_t, time_filter)
import radiation_boundary_constant as RB_cte
import radiation_boundary_estimatev2 as RB_est
import radiation_boundary_Pedro as RB_p

def integrate(xmax, xmin, ymax, ymin, delta, H, tmax,
              rotation, forcing, boundary, radiation_type, D, gamma,
              frecplot):
    
    
    print('----------------------------------------')
    print("Running Shallow Water Model using C-grid")
    print(" ")
    print("Parameters: ")
    print('tmax = '+str(tmax))
    print('rotation = '+str(rotation))
    print('forcing type = '+str(forcing))
    print('----------------------------------------')
    
    
    # ---------------------------------------------
    # Data to be expoerted at the end of the model run  
    
    start_time = time.time() # store initial time   
    result_u = []
    result_v = []
    result_h = []
    result_csi = []
    # for checking conservation of mass and energy
    Vol = [] # volume
    Ep = []  # potential energy
    Ek = []  # kinetic energy
    ape = [] # absolute potential enstrophy
    # data from one grid point to be stored at all time steps
    histgrid1 = [] 
    histgrid2 = []  
    histgrid3 = []   
    histgrid4 = [] 
    
    # ---------------------------------------------
    # Create grids
    
    dx, dy = delta, delta
    
    tmp = create_grid(xmax, xmin, ymax, ymin, delta)
    
    lonsu,latsu = tmp[0]['lonsu'],tmp[0]['latsu']
    lonsv,latsv = tmp[0]['lonsv'],tmp[0]['latsv']
    lonsz,latsz = tmp[0]['lonsz'],tmp[0]['latsz']    
    nx,ny = tmp[1],tmp[2]
    dt = tmp[3]
    
    fu = coriolis(latsu,rotation)
    fv = coriolis(latsv,rotation)
    fz = coriolis(latsz,rotation)

    g = 9.8
    
    # ---------------------------------------------
    # Define forcing    
    
    if forcing == 2:
        a,cx,cy,nrx,nry,dx = 1,0,0,10,5,delta
    if forcing == 3:
        a,cx,cy,nrx,nry,dx = .1,0,0,15,15,delta            
        
    gauss = gauss_space(xmin, xmax, ymin, ymax, nx, ny, a, cx, cy, nrx, nry, dx, dy)
    
    if forcing == 1:
        alpha = 0.8
    if forcing == 2 or forcing ==3:
        alpha = 0.02
        
    # ---------------------------------------------        
    # Initial conditions  

    u_next =  np.zeros((ny,nx+1))
    v_next =  np.zeros((ny+1,nx))
    h_next =  np.zeros((ny,nx)) + H
    U_next = fluxU(u_next,h_next,nx,ny)
    V_next = fluxV(v_next,h_next,nx,ny)    

    # Update prognostic matrices       
    
    u_prev,v_prev,h_prev = u_next*np.nan,v_next*np.nan,h_next*np.nan              
    u_curr,v_curr,h_curr = u_next,v_next,h_next
    u_next,v_next,h_next = u_curr*np.nan,v_curr*np.nan,h_curr*np.nan 
    
    U_prev,V_prev = u_next*np.nan,v_next*np.nan             
    U_curr,V_curr = U_next,V_next
    U_next,V_next = u_curr*np.nan,v_curr*np.nan
    
    # ---------------------------------------------
    # Time-loop:
    for t in range(1,tmax):
        
        frc = gauss*decay(alpha,t,H)
        
        # Calculate diagnostic matrices
        U_next = fluxU(u_curr,h_curr,nx,ny)
        V_next = fluxV(v_curr,h_curr,nx,ny)
        B_next = Bournelli(u_curr,v_curr,h_curr,nx,ny) 
        csi_next = csi(u_curr,v_curr,h_curr,fz,nx,ny,dx,dy)

        if t > 1:            
            if boundary == 'radiation':
                if radiation_type == 'constant':
                    U_next = RB_cte.west(U_next,U_curr,V_curr,fv,H,dx,dt)
                    V_next = RB_cte.north(V_next,V_curr,U_curr,fu,H,dy,dt)         
                    V_next = RB_cte.south(V_next,V_curr,U_curr,fu,H,dy,dt)         
                elif radiation_type == 'estimate':
                    U_next = RB_est.west(U_next,U_curr,U_prev,V_curr,fv,H,dx,dt)
                    V_next = RB_est.north(V_next,V_curr,V_prev,U_curr,fu,H,dy,dt)         
                    V_next = RB_est.south(V_next,V_curr,V_prev,U_curr,fu,H,dy,dt)         
                elif radiation_type == 'Pedro': 
                    U_next[:,-1] = 0
                    U_next = RB_p.west(U_next,U_curr,U_prev,V_curr,fv,H,dx,dt)
                    V_next = RB_p.north(V_next,V_curr,V_prev,U_curr,fu,H,dy,dt)           
                    V_next = RB_p.south(V_next,V_curr,V_prev,U_curr,fu,H,dy,dt)               
    
        # First time-step:
        # discretize using Euler forward in time and centered in space: 
        if t == 1:
                                        
            # calculate u field
            V_mean = .5 * ((csi_next[:-1,1:-1] * \
                            .5 * (V_next[:-1,:-1] + V_next[:-1,1:])) +
                           (csi_next[1:,1:-1] * \
                            .5 * (V_next[1:,:-1] + V_next[1:,1:])))
            delta_Bu = (B_next[:,1:] - B_next[:,:-1])/delta
            u_next[:,1:-1] = u_curr[:,1:-1] + dt * (V_mean-delta_Bu)           
            
            # calculate v field
            U_mean = .5 * ((csi_next[1:-1,1:] * \
                            .5 * (U_next[1:,1:] + U_next[:-1,1:])) +
                           (csi_next[1:-1,:-1] * \
                            .5 * (U_next[1:,:-1] + U_next[:-1,:-1])))
            delta_Bv = (B_next[:-1] - B_next[1:])/delta            
            v_next[1:-1] = v_curr[1:-1] + dt * (-U_mean-delta_Bv)
            
            # calculate h field
            delta_U =  (U_next[:,1:] - U_next[:,:-1])/delta
            delta_V = (V_next[:-1,:] - V_next[1:,:])/delta            
            h_next = h_curr + dt *  (-delta_U-delta_V)
        
        # Others time-steps
        # discretize using Leapfrog scheme:  
          
        else:           
        
            # calculate u field
            V_mean = .5 * ((csi_next[:-1,1:-1] * \
                            .5 * (V_next[:-1,:-1] + V_next[:-1,1:])) +
                           (csi_next[1:,1:-1] * \
                            .5 * (V_next[1:,:-1] + V_next[1:,1:])))
            delta_Bu = (B_next[:,1:] - B_next[:,:-1])/delta
            u_next[:,1:-1] =  u_prev[:,1:-1] + 2*dt * (V_mean-delta_Bu)
                 
            # calculate v field
            U_mean = .5 * ((csi_next[1:-1,1:] * \
                            .5 * (U_next[1:,1:] + U_next[:-1,1:])) +
                           (csi_next[1:-1,:-1] * \
                            .5 * (U_next[1:,:-1] + U_next[:-1,:-1])))
            delta_Bv = (B_next[:-1] - B_next[1:])/delta     
            v_next[1:-1] = v_prev[1:-1] + 2*dt * (-U_mean-delta_Bv)
            
            # calculate h field
            delta_U =  (U_next[:,1:] - U_next[:,:-1])/delta
            delta_V = (V_next[:-1,:] - V_next[1:,:])/delta 
            h_next = h_prev + 2*dt *  (-delta_U-delta_V)        

        # ---------------------------------------------    
        # Add forcing 
        if forcing == 2: # add forcing
            u_next[1:-1,1:-1] =  u_next[1:-1,1:-1] + frc[1:-1,1:-1]        
  
        if forcing == 3: # add forcing
            h_next[1:-1,1:-1] =  h_next[1:-1,1:-1] + .5*(frc[1:-1,2:-1] + frc[1:-1,1:-2] )     
                        
        # ---------------------------------------------    
        # Boundaries
        
        # west is always fixed        
        if boundary == 'fixed':
            u_next[:,-1], u_next[:,0],v_next[0],v_next[-1] = 0,0,0,0
            
        if boundary == 'radiation':

            if radiation_type == 'constant':
                u_next[:,-1] = 0
                v_next = RB_cte.north(v_next,v_curr,u_curr,fu,H,dy,dt)         
                v_next = RB_cte.south(v_next,v_curr,u_curr,fu,H,dy,dt)
                u_next = RB_cte.west(u_next,u_curr,v_curr,fv,H,dx,dt)
                
            elif radiation_type == 'estimate':
                    u_next[:,-1] = 0
                    u_next = RB_est.west(u_next,u_curr,u_prev,v_curr,fv,H,dx,dt)
                    v_next = RB_est.north(v_next,v_curr,v_prev,u_curr,fu,H,dy,dt)         
                    v_next = RB_est.south(v_next,v_curr,v_prev,u_curr,fu,H,dy,dt)                
            elif radiation_type == 'Pedro':
                    u_next[:,-1] = 0
                    u_next = RB_p.west(u_next,u_curr,u_prev,v_curr,fv,H,dx,dt)
                    v_next = RB_p.north(v_next,v_curr,v_prev,u_curr,fu,H,dy,dt)         
                    v_next = RB_p.south(v_next,v_curr,v_prev,u_curr,fu,H,dy,dt)                             

        # ---------------------------------------------                
        # Apply filters             
        if t > 3:
            
            if D > 0:
            
                # Difusion in x
                u_curr  = difusion_x(u_curr,u_prev,D,dx,dt)
                v_curr  = difusion_x(v_curr,v_prev,D,dx,dt)
                h_curr  = difusion_x(h_curr,h_prev,D,dx,dt)
    
#                # Difusion in y            
#                u_curr  = difusion_y(u_curr,u_prev,D,dx,dt)
#                v_curr  = difusion_y(v_curr,v_prev,D,dx,dt)
#                h_curr  = difusion_y(h_curr,h_prev,D,dx,dt)
            
            if gamma > 0:
            
                # Robert-Asselin in time  
                u_curr  = RA_filter_t(u_next,u_curr,u_prev,gamma)
                v_curr  = RA_filter_t(v_next,v_curr,v_prev,gamma)
                h_curr  = RA_filter_t(h_next,h_curr,h_prev,gamma)
                
                # Robert-asselin-Willians
                u_next  = RAW_filter(u_next,u_curr,u_prev,gamma,alpha)[0]
                v_next  = RAW_filter(v_next,v_curr,v_prev,gamma,alpha)[0]
                h_next  = RAW_filter(h_next,h_curr,h_prev,gamma,alpha)[0]
            
                u_curr  = RAW_filter(u_next,u_curr,u_prev,gamma,alpha)[1]
                v_curr  = RAW_filter(v_next,v_curr,v_prev,gamma,alpha)[1]
                h_curr  = RAW_filter(h_next,h_curr,h_prev,gamma,alpha)[1]
    
                # Time-filter          
                u_curr  = time_filter(u_next,u_curr,u_prev)
                v_curr  = time_filter(v_next,v_curr,v_prev)
                h_curr  = time_filter(h_next,h_curr,h_prev)
            
                # Robert-Asselin in x-direction          
                u_curr  = RA_filter_x(u_curr,gamma)
                v_curr  = RA_filter_x(v_curr,gamma)
                h_curr  = RA_filter_x(h_curr,gamma)
                
#                # Robert-Asselin in y-direction                  
#                u_curr  = RA_filter_y(u_curr,gamma)
#                v_curr  = RA_filter_y(v_curr,gamma)
#                h_curr  = RA_filter_y(h_curr,gamma)

        # ---------------------------------------------                 
        # Update prognostic matrices        
        u_prev,v_prev,h_prev = u_curr,v_curr,h_curr           
        u_curr,v_curr,h_curr = u_next,v_next,h_next
        u_next,v_next,h_next = u_curr*np.nan,v_curr*np.nan,h_curr*np.nan
        
        # Update diagnostic matrices
        U_prev,V_prev = U_curr,V_curr                   
        U_curr,V_curr,csi_curr = U_next,V_next,csi_next        
        U_next,V_next,csi_next = U_curr*np.nan,V_curr*np.nan,csi_curr*np.nan
        
        # ---------------------------------------------    
        ## Computate conservative properties
        
        # total volume        
        Vol.append(((h_curr)*delta**2).sum())
        
        # potential energy
        Ep.append(((h_curr)**2*delta**2).sum()*g/2)
        
        # kinectic energy
        EKu = (u_curr[:,:-1]+u_curr[:,1:])/2
        EKv = (v_curr[:-1]+v_curr[1:])/2
        h_meanxy = h_curr.mean()
        Ek.append((((EKu**2 + EKv**2)*h_meanxy/2*(delta**2)).sum()))
        
        # enstrophy
        h_meanx = (h_curr[:,:-1] + h_curr[:,1:])/2
        h_meanxy = (h_meanx[:-1] + h_meanx[1:])/2 
        tmp = csi_curr*np.nan
        tmp[1:-1,1:-1] = (csi_curr[1:-1,1:-1]**2)*h_meanxy
        tmp[0,:-1] = (csi_curr[0,:-1]**2)*h_curr[0]
        tmp[-1,1:] = (csi_curr[-1,1:]**2)*h_curr[-1]
        tmp[1:,0] = (csi_curr[1:,0]**2)*h_curr[:,0]
        tmp[:-1,-1] = (csi_curr[:-1,-1]**2)*h_curr[:,-1]
        ape.append((tmp*(delta**2)).sum()/2)

        # ---------------------------------------------            
        # Sotre middle points
        histgrid1.append(h_curr[40,40])
        histgrid2.append(u_curr[40,41])
        histgrid3.append(v_curr[41,40])
        histgrid4.append(csi_curr[41,41])
        
        # ---------------------------------------------            
        # Track mean CFL
        cmax = np.sqrt(np.amax(h_curr)*g)
        cflmax = cmax * dt/dx
        print('timestep = '+str(t)+', cfl max: '+str(cflmax))
 
        # ---------------------------------------------            
        # Store results at some time steps
        if t % frecplot == 0:    
            print('')            
            print('Storing results, timestep = '+str(t)+\
                  ' (time = '+str(round(t*dt/60/60/24,2))+'d, '+ \
                  str(round(t*dt/60/60,2))+'h)...')
            result_u.append(u_curr)
            result_v.append(v_curr)
            result_h.append(h_curr)
            result_csi.append(csi_curr)
            
            mx1,mx2,mx3,mx4 = round(np.amax(h_curr)),round(np.amax(u_curr)),round(np.amax(v_curr)),round(np.amax(csi_curr[1:-1,1:-1] )),
            m1,m2,m3,m4 = round(np.mean(h_curr)),round(np.mean(u_curr)),round(np.mean(v_curr)),round(np.mean(csi_curr[1:-1,1:-1] ))  
            mn1,mn2,mn3,mn4 = round(np.amin(h_curr)),round(np.amin(u_curr)),round(np.amin(v_curr)),round(np.amin(csi_curr[1:-1,1:-1] ))             
            
            print('max values of h, u, v and csi: '+str(mx1),str(mx2),str(mx3),str(mx4))
            print('mean values of h, u, v and csi: '+str(m1),str(m2),str(m3),str(m4))
            print('min values of h, u, v and csi: '+str(mn1),str(mn2),str(mn3),str(mn4))
            print('')

    # ---------------------------------------------                                
    # After finishing model integration:
    #   plot h, u and v at the center of each grid
    tx = np.arange(0,len(histgrid1)*dt/60/60/24,dt/60/60/24)
    fig,axs = plt.subplots(4,figsize=(10, 10), constrained_layout=True)
    axs[0].plot(tx,histgrid1,linewidth=6,color='b')
    axs[0].set_title( 'h at point 40,40',color='b', fontsize = 22)
    axs[0].tick_params(axis='both', which='major', labelsize=16)
    axs[1].plot(tx,histgrid2,linewidth=6,color='g')
    axs[1].set_title('u at point 40,41',color='g', fontsize = 22)
    axs[1].tick_params(axis='both', which='major', labelsize=16)    
    axs[2].plot(tx,histgrid3,linewidth=6,color='r')
    axs[2].set_title('v at point 41,40',color='r', fontsize = 22)
    axs[2].tick_params(axis='both', which='major', labelsize=16)
    axs[2].set_xlabel('time (days)', fontsize=18)   
    axs[3].plot(tx,histgrid4,linewidth=6,color='y')
    axs[3].set_title('csi at point 41,41',color='y', fontsize = 22)
    axs[3].tick_params(axis='both', which='major', labelsize=16)
    axs[3].set_xlabel('time (days)', fontsize=18)
    pl.savefig('histgrid.png')     
     
    
    # plot energy and mass
    fig,axs = plt.subplots(4,figsize=(10, 10), constrained_layout=True)
    axs[0].plot(tx,Vol,linewidth=6,color='b')
    axs[0].set_title('Vol',color='b', fontsize = 22)
    axs[0].tick_params(axis='both', which='major', labelsize=16)    
    axs[1].plot(tx,Ep,linewidth=6,color='g')
    axs[1].set_title('Ep',color='g', fontsize = 22)
    axs[1].tick_params(axis='both', which='major', labelsize=16)
    axs[2].plot(tx,Ek,linewidth=6,color='r')
    axs[2].set_title('Ek',color='r', fontsize = 22)      
    axs[2].tick_params(axis='both', which='major', labelsize=16)
    axs[3].plot(tx,ape,linewidth=6,color='y')
    axs[3].set_title('Abs. P. Entrophy',color='y', fontsize = 22)      
    axs[3].tick_params(axis='both', which='major', labelsize=16)
    pl.savefig('properties.png')     
            
    endtime = (time.time() - start_time)
    print(" ")
    print('Finished!')
    print('Total time elapsed = '+str(endtime)+' seconds')
    print(" ")
            
    return result_u, result_v, result_h, result_csi, [lonsu,lonsz,lonsv], [latsu,latsz,latsv],\
            endtime, Vol, Ep, Ek, ape, histgrid1, histgrid2, histgrid3, dt     
    
         