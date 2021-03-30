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
         
    c = np.zeros(u_next[:,0].shape)
    
    for i in range(len(u_next[:,0])):
                
        if u_next[i,1] < 0:
            
            fv_mean = .25*((fv[i,0] * (v_curr[i,0] + v_curr[i,1])) +
                             (fv[i+1,0] *(v_curr[i+1,0] + v_curr[i+1,1])))
            fv_mean = 0
            c[i] = -((((u_curr[i,1] - u_prev[i,1])/dt) - fv_mean)/
            ((u_prev[i,1] - u_prev[i,0])/dx))            
            c[np.isinf(c)]=np.sqrt(g*H)
            c[np.isnan(c)]=np.sqrt(g*H)
            c[np.absolute(c)>(dx/dt)]=np.sqrt(g*H)
            c[c<0]=0
#            fv_mean = .25 * ((fv[i,0] * v_curr[i,0]) +
#                             (fv[i+1,0] * v_curr[i+1,0]))
            u_next[i,0] = u_curr[i,0] - dt * (
                    (c[i]  *  ((u_curr[i,1] - u_curr[i,0])/dx)) - fv_mean)
                    
        else: 
            u_next[i,0] = u_next[i,1]
    
    return u_next


def north(v_next,v_curr,v_prev,u_curr,fu,H,dy,dt):
       
    c = np.zeros(v_next[0].shape)
    
    for i in range(len(v_next[0])):
        
        if v_next[1,i] > 0:
            
            fu_mean = .25*((fu[0,i] * (u_curr[0,i] + u_curr[0,i+1])) +
                             (fu[1,i] *(u_curr[1,i] + u_curr[1,i+1])))
            c[i] = -((((v_curr[1,i] - v_prev[1,i])/dt) + fu_mean)/
             ((v_prev[1,i] - v_prev[2,i])/dy)) 
            c[np.isinf(c)]=np.sqrt(g*H)
            c[np.isnan(c)]=np.sqrt(g*H)
            c[c>(dy/dt)]=np.sqrt(g*H)  
            c[c<0]=0            
#            fu_mean = .25 * fu[0,i] * (u_curr[0,i] +  u_curr[0,i+1])
            v_next[0,i] = v_curr[0,i] - dt * (
                    (c[i] * ((v_curr[0,i] - v_curr[1,i])/dy)) + fu_mean)

        else: 
             v_next[0,i] = v_next[1,i]
                    
    return v_next


def south(v_next,v_curr,v_prev,u_curr,fu,H,dy,dt):

    
    c = np.zeros(v_next[0].shape)
    
    for i in range(len(v_next[-1])):  
        
        if v_next[-2,i] < 0:  
            
#            fu_mean = .25*((fu[-2,i] * (u_curr[-2,i] + u_curr[-2,i+1])) +
#                             (fu[-3,i] *(u_curr[-3,i] + u_curr[-3,i+1])))
            fu_mean = .25*((fu[-1,i] * (u_curr[-1,i] + u_curr[-1,i+1])) +
                             (fu[-2,i] *(u_curr[-2,i] + u_curr[-2,i+1])))             
            c[i] =  -((((v_curr[-2,i] - v_prev[-2,i])/dt) + fu_mean)/ 
            ((v_prev[-3,i] - v_prev[-2,i])/dy))
            c[np.isinf(c)]=-np.sqrt(g*H)
            c[np.isnan(c)]=-np.sqrt(g*H)
            c[np.absolute(c)>(dy/dt)]=-np.sqrt(g*H)
            c[c>0]=0
#            fu_mean = .25 * fu[-1,i] * (u_curr[-1,i] + u_curr[-1,i+1])
#            fu_mean = .25*((fu[-1,i] * (u_curr[-1,i] + u_curr[-1,i+1])) +
#                             (fu[-2,i] *(u_curr[-2,i] + u_curr[-2,i+1])))            
            v_next[-1,i] = v_curr[-1,i] - dt * (
                    (c[i] * ((v_curr[-2,i] - v_curr[-1,i])/dy)) + fu_mean) 
             
        else: 
            v_next[-1,i] = v_next[-2,i]
            
            
    return v_next
 
