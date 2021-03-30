#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec 22 14:57:15 2020

@author: Danilo
"""

g = 9.8

import numpy as np

def fluxU(u_curr,h_curr,nx,ny):

    # flux U:
    U_next = np.zeros((ny,nx+1))
    U_next[:,1:-1] = u_curr[:,1:-1] * .5 * (h_curr[:,:-1] + h_curr[:,1:])
    U_next[:,0] = u_curr[:,0] * h_curr[:,0]
    U_next[:,-1] = u_curr[:,-1] * h_curr[:,-1]
   
    return U_next

def fluxV(v_curr,h_curr,nx,ny):

    # flux V:
    V_next = np.zeros((ny+1,nx))  
    V_next[1:-1] = v_curr[1:-1] * .5 * (h_curr[1:] + h_curr[:-1])
    V_next[0] = v_curr[0] * h_curr[0]
    V_next[-1] = v_curr[-1] * h_curr[-1]
    
   
    return V_next

def Bournelli(u_curr,v_curr,h_curr,nx,ny):
        
    # calculte Bournelli:  
    B_next = np.zeros((ny,nx))         
    u_mean = .5 * (u_curr[:,1:]**2 + u_curr[:,:-1]**2)
    v_mean = .5 * (v_curr[:-1]**2 + v_curr[1:]**2)
    B_next = g*(h_curr) + .5*(u_mean+v_mean)
    
    return B_next

def csi(u_curr,v_curr,h_curr,fz,nx,ny,dx,dy):
    
    csi_next = np.zeros((ny+1,nx+1))*np.nan
    delta_v = (v_curr[1:-1,1:]-v_curr[1:-1,:-1])/dx
    delta_u = (u_curr[:-1,1:-1] - u_curr[1:,1:-1])/dy
    h_mean = .25*(h_curr[1:,:-1] + h_curr[1:,1:] + h_curr[:-1,:-1] + h_curr[:-1,1:])
    csi_next[1:-1,1:-1] = (fz[1:-1,1:-1] + delta_v - delta_u)/h_mean
    
	#NORTE E SUL
    csi_next[0,1:-1] = (fz[0,1:-1] + (v_curr[0,1:]-v_curr[0,0:-1])/dx)/((h_curr[0,:-1]+h_curr[0,1:])/2) 
    csi_next[-1,1:-1] = (fz[-1,1:-1] + (v_curr[-1,1:]-v_curr[-1,0:-1])/dx)/((h_curr[-1,0:-1]+h_curr[-1,1:])/2) 
	#LESTE E OESTE
    csi_next[1:-1,0] = (fz[1:-1,0] - (u_curr[1:,0]-u_curr[0:-1,0])/dy)/((h_curr[0:-1,0]+h_curr[1:,0])/2) 
    csi_next[1:-1,-1] = (fz[1:-1,-1] - (u_curr[1:,-1]-u_curr[0:-1,-1])/dy)/((h_curr[0:-1,-1]+h_curr[1:,-1])/2)    
    
    csi_next[np.isnan(csi_next)]=0
    
    return csi_next

