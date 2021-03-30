#------------------------------------------------------------------------- 
#  
# Shallow Water Model
#
# Danilo Couto de Souza
# Universidade de São Paulo (USP)
# Instituto de Astronomia, Geociências e Ciências Atmosféricas
#-------------------------------------------------------------------------


def difusion_x(var_curr,var_prev,D,dx,dt):
    
    F = D*dt/(dx**2)
    var_curr[:,1:-1] = var_curr[:,1:-1] + F * (
    var_prev[:,2:] - 2*var_prev[:,1:-1] + var_prev[:,:-2])
    
    return var_curr 
    
def difusion_y(var_curr,var_prev,D,dy,dt):
    
    F = D*dt/(dy**2)
    var_curr[1:-1] = var_curr[1:-1] + F * (
    var_prev[:-2] - 2*var_prev[1:-1] + var_prev[2:])
    
    return var_curr     

def RA_filter_x(var,gamma):              
                
    var[:,1:-1] = var[:,1:-1] + gamma * (var[:,:-2] - 2*var[:,1:-1] + var[:,2:])

    return var

def RA_filter_y(var,gamma):              

    var[1:-1] = var[1:-1] + gamma * (var[2:] - 2*var[1:-1] + var[:-2])
    
    return var   

def RA_filter_t(var_next,var_curr,var_prev,gamma):  
    
    var_curr = var_curr + gamma * (var_prev - 2*var_curr + var_next)
    
    return var_curr
    
def RAW_filter(var_next,var_curr,var_prev,gamma,alpha):
                
    f1 = gamma*(1-alpha)/2
    f2 = gamma*alpha/2
    #north
    var_next = var_next - f1 * (var_next - 2*var_curr + var_prev)
    var_curr = var_curr + f2 * (var_next - 2*var_curr + var_prev)
    
    return var_next,var_curr

def time_filter(var_next,var_curr,var_prev):
    
    var_curr = (var_prev + 2*var_curr + var_next)/4
    
    return var_curr
