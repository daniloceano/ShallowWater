#------------------------------------------------------------------------- 
#  
# Shallow Water Model
#
# Danilo Couto de Souza
# Universidade de São Paulo (USP)
# Instituto de Astronomia, Geociências e Ciências Atmosféricas
#-------------------------------------------------------------------------

import numpy as np
import math


def gauss_space(xmin, xmax, ymin, ymax, nx, ny, a, cx, cy, nrx, nry, dx, dy):
 # --------------------------
    # Pertubation (2D Gaussian in space):
        # a: amplitude of the pertubation
        # cx: displacement of the peak on x dims (in meters)
        # cy: displacement of the peak on y dims(in meters)
        # nrx: width of the pertubation on x dims (in grid points)
        # nry: width of the pertubation on y dis (in grid points)
        # xmax, xmin,ymax,ymin: endpoints of the grids
        # dx: dims spacing
        # dy: dims spacing
    xgauss = np.linspace(xmin, xmax,nx+1)
    ygauss = np.linspace(ymax, ymin,ny)  
    
    longa,latga = np.meshgrid(xgauss,ygauss)     
    gauss_space = a * np.exp(- (((longa-cx)**2)/(nrx*dx)**2) -\
                          (((latga-cy)**2)/(nry*dx)**2))
    
    return gauss_space


def decay(alpha,t,H):   
#    # Time component of the pertubation
#    #   this will make the gaussian to decay
    decay = (H/2)*(alpha**3)*(t**2)*math.exp(-alpha*t)
    
    return decay


## See forcing
#H = 1   
#frc = []
#gs = gauss_space(xmin, xmax, ymin, ymax, nx, ny, a, cx, cy, nrx, nry, dx, dy)
#for t in range (1000):   
#    gt = decay(0.02,t,H)
#    tmp = gs*gt
#    frc.append(tmp[40,40])
#    
#plt.plot(frc)
