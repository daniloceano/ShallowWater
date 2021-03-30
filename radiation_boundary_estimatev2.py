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
                     

#    fv_mean = .25*((fv[i,0] * (v_curr[i,0] + v_curr[i,1])) +
#                     (fv[i+1,0] *(v_curr[i+1,0] + v_curr[i+1,1])))
    fv_mean = 0
    c = -((((-u_curr[:,1] - u_prev[:,1])/dt) - fv_mean)/
    ((u_prev[:,1] - u_prev[:,2])/dx))
    c[np.isinf(c)]=np.sqrt(g*H)
    c[np.isnan(c)]=np.sqrt(g*H)
    c[np.absolute(c)>(dx/dt)]=np.sqrt(g*H)
    c[c<0]=0
    u_next[:,0] = u_curr[:,0] - dt*fv_mean -c*(dt/dx) * (u_curr[:,0] - u_curr[:,1])       


            
    return u_next

def south(v_next,v_curr,v_prev,u_curr,fu,H,dy,dt):

    fu_mean = ((fu*u_curr)[-2,1:] +(fu*u_curr)[-1,1:] +(fu*u_curr)[-2,:-1] +(fu*u_curr)[-1,:-1])/4   
    c =  -((((v_next[-2] - v_curr[-2])/dt) + fu_mean)/ 
    ((v_curr[-2] - v_curr[-3])/dy))
    c[np.isinf(c)]=np.sqrt(g*H)
    c[np.isnan(c)]=np.sqrt(g*H)
    c[np.absolute(c)>(dy/dt)]=np.sqrt(g*H)
    c[c<0]=0 
    v_next[-1,:] = v_curr[-1,:]-dt*fu_mean-c*(dt/dy)*(v_curr[-1,:]-v_curr[-2,:]) 
   
    return v_next

def north(v_next,v_curr,v_prev,u_curr,fu,H,dy,dt):
            
    fu_mean = (((fu*u_curr)[0,1:] +(fu*u_curr)[1,1:] +(fu*u_curr)[0,0:-1] +(fu*u_curr)[1,0:-1,])/4)
    c = -((((v_next[1] - v_curr[1])/dt) + fu_mean)/
     ((v_curr[1] - v_curr[2])/dy)) 
    c[np.isinf(c)]=np.sqrt(g*H)
    c[np.isnan(c)]=np.sqrt(g*H)
    c[c>(dy/dt)]=np.sqrt(g*H)  
    c[c<0]=0            
    v_next[0] = v_curr[0] - dt * (
            (c * ((v_curr[0] - v_curr[1])/dy)) + fu_mean)

                    
    return v_next

 
