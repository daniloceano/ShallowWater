#------------------------------------------------------------------------- 
#  
# Shallow Water Model
#
# Danilo Couto de Souza
# Universidade de São Paulo (USP)
# Instituto de Astronomia, Geociências e Ciências Atmosféricas
#-------------------------------------------------------------------------

import numpy as np

g = 9.8

def west(u_next,u_curr,u_prev,v_curr,fv,H,dx,dt):
    
    fv_mean  = 0     
    cpo = -( -fv_mean + (u_curr[:,1]-u_prev[:,1])/dt)/((u_prev[:,1]-u_prev[:,2])/dx)
    cpo[cpo>dx/dt]=dx/dt 
    cpo[np.isnan(cpo)==True]=0 
    cpo[cpo==-np.inf]=0 
    cpo[cpo==np.inf]=0
    cpo[cpo<0]=0 
    u_next[:,0] = u_curr[:,0]-cpo*(dt/dx)*(u_curr[:,0]-u_curr[:,1])
                
    return u_next


def south(v_next,v_curr,v_prev,u_curr,fu,H,dy,dt):
    
    fu_mean = ((fu*u_curr)[-2,1:] +(fu*u_curr)[-3,1:] +(fu*u_curr)[-2,0:-1] +(fu*u_curr)[-3,0:-1,])/4
    cpn = -(fu_mean + (v_curr[-2,:]-v_prev[-2,:])/dt)/((v_prev[-2,:]-v_prev[-3,:])/dy)
    cpn[cpn>dy/dt]=dy/dt 
    cpn[np.isnan(cpn)==True]=0 
    cpn[cpn==-np.inf]=0
    cpn[cpn==np.inf]=0
    cpn[cpn<0]=0 
    fu_mean = ((fu*u_curr)[-2,1:] +(fu*u_curr)[-1,1:] +(fu*u_curr)[-2,0:-1] +(fu*u_curr)[-1,0:-1,])/4
    v_next[-1,:] = v_curr[-1,:]-dt*fu_mean-cpn*(dt/dy)*(v_curr[-1,:]-v_curr[-2,:]) 
             
    return v_next


def north(v_next,v_curr,v_prev,u_curr,fu,H,dy,dt):
    
    fu_mean = ((fu*u_curr)[2,1:] +(fu*u_curr)[1,1:] +(fu*u_curr)[2,0:-1] +(fu*u_curr)[1,0:-1,])/4
    cps = -(fu_mean + (v_curr[1,:]-v_prev[1,:])/dt)/((v_prev[1,:]-v_prev[2,:])/dy)
    cps[cps>dy/dt]=dy/dt 
    cps[np.isnan(cps)==True]=0 
    cps[cps==-np.inf]=0
    cps[cps==np.inf]=0
    cps[cps<0]=0
    fu_mean = (((fu*u_curr)[0,1:] +(fu*u_curr)[1,1:] +(fu*u_curr)[0,0:-1] +(fu*u_curr)[1,0:-1,])/4)
    v_next[0,:] = v_curr[0,:] -fu_mean*dt -cps*(dt/dy)*(v_curr[0,:]-v_curr[1,:])
              
    return v_next
 
