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

def west(u_next,u_curr,v_curr,fv,H,dx,dt):
      
    
    c = np.zeros(u_next[:,0].shape)
    
    for i in range(len(u_next[:,0])):
                
        if u_next[i,1] < 0 :
            
            c[i] = c[i] + np.sqrt(g*H)
               
            fv_mean = .5*fv[i,0] * (v_curr[i,0] + v_curr[i+1,0])
            u_next[i,0] = u_curr[i,0] + dt * (
                    (c[i]  *  ((u_curr[i,1] - u_curr[i,0])/dx)) + fv_mean)
        
        else: 
#            u_next[i,0] = u_curr[i,0]
            u_next[i,0] = u_next[i,1]

    u_next[np.isnan(u_next)] = 0    
    
    return u_next


def north(v_next,v_curr,u_curr,fu,H,dy,dt):
    
    
    c = np.zeros(v_next[0].shape)
    
    for i in range(len(v_next[0])):
        
        if v_next[1,i] > 0:

            c[i] = c[i] + np.sqrt(g*H) 
           
            fu_mean = .5 * fu[0,i] * (u_curr[0,i] +  u_curr[0,i+1])
#            fu_mean = vn*0                
            v_next[0,i] = v_curr[0,i] + dt * (
                    (-c[i] * ((v_curr[0,i] - v_curr[1,i])/dy)) - fu_mean)
        
        else: 
#            v_next[0,i] = v_curr[0,i]
             v_next[0,i] = v_next[1,i]
            
#    v_next[np.isnan(v_next)] = 0    
        
    return v_next


def south(v_next,v_curr,u_curr,fu,H,dy,dt):

    
    c = np.zeros(v_next[0].shape)
    
    for i in range(len(v_next[-1])):  
        
        if v_next[-2,i] < 0:  
          
            c[i] = c[i] + np.sqrt(g*H)  
          
            fu_mean = .5 * fu[-1,i] * (u_curr[-1,i] + u_curr[-1,i+1])
#            fu_mean = vs*0                            
            v_next[-1,i] = v_curr[-1,i] + dt * (
                    (c[i] * ((v_curr[-2,i] - v_curr[-1,i])/dy)) - fu_mean) 
        
        else: 
#            v_next[-1,i] = v_curr[-1,i]
            v_next[-1,i] = v_next[-2,i]
            
            
#    v_next[np.isnan(v_next)] = 0                

    return v_next
 

    