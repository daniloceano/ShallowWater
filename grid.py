#------------------------------------------------------------------------- 
#  
# Shallow Water Model
#
# Danilo Couto de Souza
# Universidade de São Paulo (USP)
# Instituto de Astronomia, Geociências e Ciências Atmosféricas
#-------------------------------------------------------------------------


import numpy as np

def create_grid(xmax, xmin, ymax, ymin, delta):
    """
    Create an Arakawa-C grid for the model
    
    The zonal wind (u) is located in the right edges,
    the meridional wind (v) is located at the top edge and
    the absolute potential vorticity (csi or z) is located at the top right vertex
    """
    
#    # ----------
#    # Grid params
#    delta = 100e3 #  grid spacing (m) 
#    xmax,xmin = 40e5,-40e5
#    ymax,ymin = 40e5,-40e5
    nx = int((xmax+delta-xmin)/delta) # (number of grid points on x)
    ny = int((ymax+delta-ymin)/delta) # (number of grid points on y)
    
    # grid for u
    xu = np.arange(-40, 42,1)*delta
    yu = np.arange(40, -41,-1)*delta
    lonsu,latsu = np.meshgrid(xu,yu) 
    # grid for csi   
    xz = np.arange(-40, 42.,1.)*delta
    yz = (np.arange(41, -41,-1)-.5)*delta
    lonsz,latsz = np.meshgrid(xz,yz) 
    # grid for v    
    xv = (np.arange(-40, 41.,1.)+.5)*delta
    yv = (np.arange(41, -41,-1)-.5)*delta
    lonsv,latsv = np.meshgrid(xv,yv)
    
    # ----------   
    # Constants
    g = 9.8 # gravity (m/s-2)
    mu = 0.2
    H = 1
    dt = (mu*delta/(np.sqrt(g*H)))
    
    grids = {'lonsu':lonsu,
             'latsu':latsu,
             'lonsz':lonsz,
             'latsz':latsz,
             'lonsv':lonsv,
             'latsv':latsv,}
    
    return grids,nx,ny,dt

def coriolis(lats,rotation):
    """
    Create matrices for the Coriolis force for each prognostic variable
    """
    # ----------
    # Rotation effects
    if rotation == 'on':
        # calculate Coriollis parameter using Beta plane approximation
        beta = (2*7.29*1e-5*np.cos(0*np.pi/180)/6400000)
    if rotation == 'off':
        # turn off Coriollis
        beta = 0
    # Coriolis force exerced by each wind component
    f = lats*beta 
    
    return f
    
    
    
    
