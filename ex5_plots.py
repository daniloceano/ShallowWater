#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov  5 13:41:26 2020

@author: Danilo
"""
import os
import glob

import numpy as np
import math
import pylab as pl

import matplotlib.pyplot as plt 
import matplotlib.colors as colors
from mpl_toolkits.axes_grid1 import make_axes_locatable
from mpl_toolkits.mplot3d import Axes3D 
import cmocean.cm as cmo
from matplotlib.colors import LightSource

import metpy.calc
from metpy.units import units


from filme import anim
from core import integrate


#-------------------------------------------------------------------------
# Animation
#-------------------------------------------------------------------------

# animate through time loop
def animate(u,v,h,csi,H,dt):

    # specify max and min contours
    minh, maxh, meanh = np.amin(h),np.amax(h),np.mean(h)
    minu, maxu, meanu = np.amin(u),np.amax(u),np.mean(u)
    minv, maxv, meanv = np.amin(v),np.amax(v),np.mean(v)
    mincsi, maxcsi, meancsi = np.amin(csi),np.amax(csi),np.mean(csi)
    normh = colors.TwoSlopeNorm(vmin=minh,
           vcenter=meanh,
            vmax=maxh)
    normu = colors.TwoSlopeNorm(vmin=minu,
           vcenter=meanu,
            vmax=maxu)  
    normv = colors.TwoSlopeNorm(vmin=minv,
           vcenter=meanv,
            vmax=maxv)
    normcsi = colors.TwoSlopeNorm(vmin=mincsi,
           vcenter=meancsi,
            vmax=maxcsi)    
    
    clevh = np.linspace(minh,maxh,10)
    clevu = np.linspace(minu,maxu,10)
    clevv = np.linspace(minv,maxv,10)
    

    # skip every 2 vectors
    skip=(slice(None,None,2),slice(None,None,2))
    
    # loops    
    for t in range(len(u)):
        
        # specify fig
        fig, axs = plt.subplots(2,2,figsize=(22, 10))
        plt.rcParams["axes.axisbelow"] = False
        
        u_mean = (u[t][:,:-1]+u[t][:,1:])/2
        v_mean = (v[t][:-1]+v[t][1:])/2
        csi_mean = (csi[t][:-1,:-1]+csi[t][1:,1:])/2
        
        # vorticitiy
        # assign units to fields for computations
        ums = u * units.meter / units.second
        vms = v * units.meter / units.second
        ums = (ums[t][:,1:] + ums[t][:,:-1])/2
        vms = (vms[t][1:] + vms[t][:-1])/2    
        delta = 100e3 * units.meter
        
        # take the last time step to stablish the vort and div levs
        vort = metpy.calc.vorticity(ums,vms,delta,delta)
        maxvort = np.amax(vort)
        clevort = np.linspace(-maxvort,maxvort,10)
               
#        # specify max and min contours
#        minh, maxh, meanh = np.amin(h[t]),np.amax(h[t]),np.mean(h[t])
#        minu, maxu, meanu = np.amin(u[t]),np.amax(u),np.mean(u[t])
#        minv, maxv, meanv = np.amin(v[t]),np.amax(v[t]),np.mean(v[t])
#        mincsi, maxcsi, meancsi = np.amin(csi[t]),np.amax(csi[t]),np.mean(csi[t])
#        normh = colors.TwoSlopeNorm(vmin=minh,
#               vcenter=meanh,
#                vmax=maxh)
#        normu = colors.TwoSlopeNorm(vmin=minu,
#               vcenter=meanu,
#                vmax=maxu)  
#        normv = colors.TwoSlopeNorm(vmin=minv,
#               vcenter=meanv,
#                vmax=maxv)
#        normcsi = colors.TwoSlopeNorm(vmin=mincsi,
#               vcenter=meancsi,
#                vmax=maxcsi)    
#        
#        clevh = np.linspace(minh,maxh,10)
#        clevu = np.linspace(minu,maxu,10)
#        clevv = np.linspace(minv,maxv,10)        
        
        
        fig.suptitle('t = '+str(round(t*frecplot*dt/60/60,2))+'h,' \
           ' ('+str(round(t*frecplot*dt/60/60/24,2))+'d)',fontsize=22)
        
        ploth = axs[0,0].pcolormesh(lonsh,latsh,h[t], 
                    cmap= 'coolwarm',norm=normh)
        axs[0,0].contour(lonsh, latsh, h[t], clevh, colors='k')        
        axs[0,0].quiver(lonsh[::5],latsh[::5],
                      u_mean[::5],v_mean[::5],scale = 0.2)
        axs[0,0].text(-3700000,4100000,
           'Pertubation height h(x,y,t) (m) and velocity field ' + r'$\vec{V}$'+' (m/s)' ,
           fontsize=18) 
        axs[0,0].grid(True, color="grey", linestyle='--',lw=2)        
        
        plotcsi = axs[0,1].pcolormesh(lonsh,latsh,csi_mean,
                   cmap= 'RdBu',norm=normcsi)
        axs[0,1].contour(lonsh, latsh, vort, clevort, colors='k')
        if t > 0:
            plotcsi = axs[0,1].contourf(lonsh, latsh, vort,
               levels=clevort, cmap='seismic',alpha=0.5)
        axs[0,1].text(-3700000,4100000,'Abs. potential (shaded) and vertical (contour) vorticity '  + r'$\xi$'+' (u,v,t) (m)',fontsize=18)
        axs[0,1].grid(True, color="grey", linestyle='--',lw=2)

        
        plotu = axs[1,0].pcolormesh(lonsu,latsu,u[t],
                   cmap= 'Spectral',norm=normu)
        axs[1,0].contour(lonsu, latsu, u[t], clevu, colors='k')
        axs[1,0].text(-3700000,4100000,'Zonal wind ' + r'$\vec{u}$'+' (m/s)',
           fontsize=18)
        axs[1,0].grid(True, color="grey", linestyle='--',lw=2)
        
        
        plotv = axs[1,1].pcolormesh(lonsv,latsv,v[t], 
                   cmap= 'BrBG',norm=normv)
        axs[1,1].contour(lonsv, latsv, v[t], clevv, colors='k')
        axs[1,1].text(-3700000,4100000,'Meridional wind '+ r'$\vec{v}$'+' (m/s)',
           fontsize=18)
        axs[1,1].grid(True, color="grey", linestyle='--',lw=2)
               
        plots = [ploth,plotcsi,plotu,plotv]
        
        i = 0
        for row in range(2):
            for col in range(2):
                divider = make_axes_locatable(axs[row,col])
                cax = divider.append_axes('right', size='5%', pad=0.05)      
                fig.colorbar(plots[i],cax = cax) 
                i = i+1                
        
        pl.savefig('./h_wnd/'+str(t*frecplot)+'.png')  
        print('saved figure: ./h_wnd/'+str(t*frecplot)+'.png')   


#-------------------------------------------------------------------------

# animate through time loop
def animate_vort_div(u,v,h,H,dt):
    
    print('Saving vorticity and divergence figures...')
    
    # assign units to fields for computations
    ums = u * units.meter / units.second
    vms = v * units.meter / units.second
    delta = 100e3 * units.meter
    
#    # specify max and min for pcolormesh
#    minh, maxh, meanh = np.amin(h),np.amax(h), np.mean(h)
#    minvort, maxvort, meanvort = np.amin(vort),np.amax(vort), np.mean(vort)
#    mindiv, maxdiv, meandiv = np.amin(div),np.amax(div), np.mean(div)  
    minh, maxh, meanh = [], [], []
    minvort, maxvort, meanvort = [], [], []
    mindiv, maxdiv, meandiv = [], [], []
    for t in range(len(u)):
        vort = metpy.calc.vorticity(ums[t][:,:-1],vms[t][:-1],delta,delta)
        div = metpy.calc.divergence(ums[t][:,:-1],vms[t][:-1],delta,delta)
        uni = vort.units
        vort,div = vort/uni,div/uni
        
        minh_, maxh_, meanh_ = np.amin(h),np.amax(h), np.mean(h)
        minvort_, maxvort_, meanvort_ = np.amin(vort),np.amax(vort), np.mean(vort)
        mindiv_, maxdiv_, meandiv_ = np.amin(div),np.amax(div), np.mean(div)

        minh.append(minh_), maxh.append(maxh_), meanh.append(meanh_)
        minvort.append(minvort_), maxvort.append(maxvort_), meanvort.append(meanvort_)
        mindiv.append(mindiv_), maxdiv.append(maxdiv_), meandiv.append(meandiv_)

    minh, maxh, meanh = np.amin(minh),np.amax(maxh), np.mean(meanh)
    minvort, maxvort, meanvort = np.amin(minvort),np.amax(maxvort), np.mean(meanvort)
    mindiv, maxdiv, meandiv = np.amin(mindiv),np.amax(maxdiv), np.mean(meandiv)
        
        
     # skip every 2 vectors
    skip=(slice(None,None,2),slice(None,None,2))
    
    
    normvort = colors.TwoSlopeNorm(vmin=minvort,
           vcenter=meanvort,
            vmax=maxvort)
    normdiv = colors.TwoSlopeNorm(vmin=mindiv,
           vcenter=meandiv,
            vmax=maxdiv)
    normh = colors.TwoSlopeNorm(vmin=minh,
           vcenter=meanh,
            vmax=maxh) 
    
    # specify contours
    clevh = np.linspace(minh,maxh,10)
    clevvort = np.linspace(minvort,maxvort,10)
    clevdiv = np.linspace(mindiv,maxdiv,10)
    
    # loops    
    for t in range(len(u)):
    
        # create fig
        fig, axs = plt.subplots(3,1,figsize=(22, 16))
#        plt.tight_layout()
        
        # global title
        fig.suptitle('t = '+str(round(t*frecplot*dt/60/60,2))+'h,' \
           ' ('+str(round(t*frecplot*dt/60/60/24,2))+'d)',fontsize=22)
        
        # assign units to fields for computations
        u_ = u[t]
        v_ = v[t]
        
        # interpolate u and v into h grid
        u_mean = (u_[:,:-1]+u_[:,1:])/2
        v_mean = (v_[:-1]+v_[1:])/2
        hgt = h[t]
        
        # compute vort and div for the time step
        vort = metpy.calc.vorticity(ums[t][:,:-1],vms[t][:-1],delta,delta)
        div = metpy.calc.divergence(ums[t][:,:-1],vms[t][:-1],delta,delta)
        uni = vort.units
        vort,div = vort/uni,div/uni
#        
#        minh, maxh, meanh = np.amin(h[t]),np.amax(h[t]), np.mean(h[t])
#        minvort, maxvort, meanvort = np.amin(vort[t]),np.amax(vort[t]), np.mean(vort[t])
#        mindiv, maxdiv, meandiv = np.amin(div[t]),np.amax(div[t]), np.mean(div[t])
#            
            
#         # skip every 2 vectors
#        skip=(slice(None,None,2),slice(None,None,2))
#        
#        
#        normvort = colors.TwoSlopeNorm(vmin=minvort,
#               vcenter=meanvort,
#                vmax=maxvort)
#        normdiv = colors.TwoSlopeNorm(vmin=mindiv,
#               vcenter=meandiv,
#                vmax=maxdiv)
#        normh = colors.TwoSlopeNorm(vmin=minh,
#               vcenter=meanh,
#                vmax=maxh) 
        
        # specify contours
        clevh = np.linspace(minh,maxh,10)
        clevvort = np.linspace(minvort,maxvort,10)
        clevdiv = np.linspace(mindiv,maxdiv,10)                

        #plot height
        ploth = axs[0].pcolormesh(lonsh,latsh,hgt, 
                    cmap='coolwarm',norm=normh)
        axs[0].quiver(lonsh[::5],latsh[::5],
                      u_mean[::5],v_mean[::5], scale = 0.4)
        axs[0].contour(lonsh, latsh, hgt, clevh, colors='k')
        axs[0].text(-4000000,4100000,'Surface elevation and Velocity field' \
           +' u(x,y,t)',
           fontsize=18) 
        
        # plot vorticity
        plotvort = axs[1].pcolormesh(lonsh,latsh,vort,
                   cmap= 'coolwarm',norm=normvort)
        axs[1].contour(lonsh, latsh, vort, clevvort, colors='k')
        axs[1].text(-4000000,4100000,'Vertical vorticity field $\zeta$ (1/s)',
           fontsize=18)
        axs[1].grid(True, color="grey", linestyle='--',lw=1)

        
        # plot divergence
        plotdiv = axs[2].pcolormesh(lonsh,latsh,div,
                   cmap= 'Spectral',norm=normdiv)
        axs[2].contour(lonsh, latsh, div, clevdiv, colors='k')
        axs[2].text(-4000000,4100000,'Wind divergence field ∇⋅'+ r'$\vec{v}$'+ ' (1/s)',
           fontsize=18)
        axs[2].grid(True, color="grey", linestyle='--',lw=1)

        
        plots = [ploth,plotvort,plotdiv]
        for row in range(3):
            divider = make_axes_locatable(axs[row])
            cax = divider.append_axes('right', size='5%', pad=0.05)      
            fig.colorbar(plots[row],cax = cax) 
    
        pl.savefig('./div_vort/'+str(t*frecplot)+'.png')  
        print('saved figure: ./div_vort/'+str(t*frecplot)+'.png')    
#-------------------------------------------------------------------------
# 3D Animation
#-------------------------------------------------------------------------
def animate3D(u,v,h,H,dt):
    
    print('Saving 3D figures...')

    # specify max and min contours
    minh, maxh = np.amin(h),np.amax(h)      
    norm_h = colors.TwoSlopeNorm(vmin=minh,
           vcenter=H,
            vmax=maxh)
    
    # colormap
    cmap = 'cmo.deep'
    
    # add light for shadings
    ls = LightSource(azdeg=315, altdeg=45)
    
    # loop to animate
    for t in range(round(tmax/frecplot)+1):
        
        # fig params
        fig = plt.figure(figsize=(20, 10))
        
        # 3D plot
        ax = fig.add_subplot(1, 2, 1, projection='3d')
        ax.set_zlim(minh, maxh)
        # view zimuth and angle
        ax.view_init(50, 260)
        # labels in km
        ax.set_xlabel('X (km)', fontweight='bold')
        ax.set_ylabel('Y (km)', fontweight='bold')
        ax.set_zlabel('Height (m)', fontweight='bold')
        
        # plot surface
        z = h[t]
        rgb = ls.shade(z, cmap=cmo.deep, vert_exag=0.1, blend_mode='soft')
        ax.plot_surface(lonsh/1e3,latsh/1e3, z,
                               rstride=1, cstride=1, facecolors=rgb,
                               cmap = cmap,
                               linewidth=0,
                               antialiased=False,
                               shade=False,
                               norm=norm_h) 
        ax.set_title('Surface elevation h(x,y,t) after t ='\
           +str(round(t*dt*frecplot/60/60/24,2))+'d, '+ \
                  str(round(t*dt/60/60,2))+'h')
        # plot vectors
        ax = fig.add_subplot(1, 2, 2)        
        u_mean = (u[t][:,:-1]+u[t][:,1:])/2
        v_mean = (v[t][:-1]+v[t][1:])/2
        ax.quiver(lonsh[::5]/1e3,latsh[::5]/1e3, u_mean[::5],v_mean[::5],scale=0.5)
        ax.set_title('Velocity field u(x,y,t) after t ='\
           +str(round(t*dt*frecplot/60/60/24,2))+'d, '+ \
                  str(round(t*dt/60/60,2))+'h')
        
        print('animating timestep = '+str(t))
        pl.savefig('./3Dfigs/'+str(t)+'.png') 

#-------------------------------------------------------------------------
# Plot forcing
#-------------------------------------------------------------------------
def anim_forcing(H,lons,lats,dt):
    # Pertubation (2D Gaussian in space):
    # a: amplitude of the pertubation
    # cx: displacement of the peak on x dims
    # cy: displacement of the peak on y dims
    # nrx: width of the pertubation on x dims
    # nry: width of the pertubation on y dis
    # dx: dims spacing
    nx = 81 # (number of grid points on x)
    ny = 81 # (number of grid points on y)
    delta = 100e3 #  grid spacing (m) 
    xmax,xmin = 4e6,-4e6
    ymax,ymin = 4e6,-4e6
    a,cx,cy,nrx,nry,dx = 1,40,41,11,5,delta
    xgauss = np.linspace(xmin, xmax,nx+1)
    ygauss = np.linspace(ymax, ymin,ny)    
    longa,latga = np.meshgrid(xgauss,ygauss)     
    gauss_space = a * np.exp(- (((longa-cx)**2)/(nrx*dx)**2) -\
                          (((latga-cy)**2)/(nry*dx)**2))
    
    normf = colors.TwoSlopeNorm(vmin=0,
           vcenter=0.05,
            vmax=0.1)
    
    alpha = 0.01
    for t in range(120):
        fig, axs = plt.subplots(figsize=(15, 10))
        #    # Time component of the pertubation
        gauss_time = (H/2)*(alpha**3)*(t**2)*math.exp(-alpha*t)
        forcing = (gauss_space * gauss_time)
        # plot
        p1 = axs.pcolormesh(lons,lats,forcing[:,:-1],cmap = 'Reds')
        axs.text(-4000000,4000000,'t = '\
           +str(round(t*dt/60/60,2))+'h',fontsize=18)
        axs.grid(True, color="grey", linestyle='--',lw=2)        
        divider = make_axes_locatable(axs)
        cax = divider.append_axes('right', size='5%', pad=0.05)
        plt.colorbar(p1,cax = cax)
        
        
        print('animating timestep = '+str(t))
        pl.savefig('./forcing/'+str(t)+'.png')     
    
    


## --------------------------
## get model data for given params
xmax, xmin, ymax, ymin, delta, H, = 40e5,-40e5, 40e5,-40e5, 100e3, 1
#D, gamma = 0, 0
D, gamma = 2e3, 0.01      
tmax = 1501
H = 1
rotation = 'on'
forcing = 2
boundary = 'radiation'
radiation_type = 'estimate'
#frecplot = round(tmax/25)
frecplot = 10
# data
#data = integrate(xmax, xmin, ymax, ymin, delta, H, tmax,
#              rotation, forcing, boundary, radiation_type, D, gamma,
#              frecplot)
u, v, h, csi = data[0], data[1], data[2], data[3]
histh,histu,histv = data[11],data[12],data[13]
lonsu, latsu = data[4][0],data[5][0]
lonsz, latsz = data[4][1],data[5][1]
lonsv, latsv = data[4][2],data[5][2]
lonsh = (lonsz[1:,1:] + lonsz[:-1,:-1])/2
latsh = (latsz[1:,1:] + latsz[:-1,:-1])/2
dt = data[-1]


## # # # # # # # #
## animate h, u, v  
## # # # # # # # #
## delete pre existing files before making new ones
#fileList = glob.glob('./h_wnd/*', recursive=True)
#for filePath in fileList:
#    try:
#        os.remove(filePath)
#    except OSError:
#        print("Error while deleting file")
## make figures
#animate(u,v,h,csi,H,dt)
## make video
#anim('./h_wnd/')
#    
#
# # # # # # # # # # # # #
# animate vort and div 
# # # # # # # # # # # # #
# delete pre existing files before making new ones
fileList = glob.glob('./div_vort/*', recursive=True)
for filePath in fileList:
    try:
        os.remove(filePath)
    except OSError:
        print("Error while deleting file")
## make figures
animate_vort_div(u,v,h,H,dt)
# make video
anim('./div_vort/')

## # # # # # # 
## animate 3D 
## # # # # # # 
## delete pre existing files before making new ones
#fileList = glob.glob('./3Dfigs/*', recursive=True)
#for filePath in fileList:
#    try:
#        os.remove(filePath)
#    except OSError:
#        print("Error while deleting file")
### make figures
#animate3D(u,v,h,H,dt)
## make video
#anim('./3Dfigs/')
#        

        