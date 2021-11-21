# -*- coding: utf-8 -*-
"""
Created on Fri May 17 07:53:39 2019

@author: Emeric Villette
"""

import numpy as np
import matplotlib.pyplot as plt
import copy
import time

start = time.time()

c=0.5
c2=c**2

L=50
T=1000



M=400
N=500
dt=T/(N+1)
dx=L/(M+1)


g=c2*dt**2/dx**2



A=(1+g)*np.eye(M)-(g/2)*np.diag(np.ones(M-1),-1)-(g/2)*np.diag(np.ones(M-1),1)
invA=np.linalg.inv(A)

def alpha(t) :
    #The alpha function represent the movement at the beginning (left) of the rope
    #In our case it is a sin function that has an amplitude that is decreasing over time
    global L
    return  (L/(t+1))*np.sin(3*np.pi*t/(L))

def beta(t) :#fin de la corde
    #The alpha function represent the movement at the end (right) of the rope
    #In our case it is a sin function that has an amplitude that is decreasing over time
    global L
    return  (L/(t+1))*np.sin(3*np.pi*t/(L))


x=np.linspace(dx,L-dx,M)
t=np.linspace(0,T,N+1)

a=alpha(t)
b=beta(t)

#a=np.zeros(N+1)
#b=np.zeros(N+1)

u0=np.zeros((M+2,1))
u0[0,0]=alpha(0)
u0[M+1,0]=beta(0)

"""
Here we have the equation that our rope will have overall
If we put that alpha and beta return 0 (no movement at the left or right)
the rope is going to have a sinusoidal movement (up and down movement)
"""
for i in range(1,M+1) :
        u0[i,0]=np.sin(dx*i*np.pi/L)
        
        
u1=np.zeros((M+2,1))        
for i in range(1,M+1) :
        u1[i,0]=u0[i,0]+(g/2)*(u0[i+1,0]-2*u0[i,0]+u0[i-1,0])
        
        
u0c=u0[1:M+1,0][np.newaxis].T        


u1c=u1[1:M+1,0][np.newaxis].T         

q=np.zeros((M,1))
q[0,0]=1

r=np.zeros((M,1))
r[M-1,0]=1

u=2*(invA  @ u1c) - u0c+g*((0.5*a[0]+0.5*a[2])*(invA @ q)+(0.5*b[0]+0.5*b[2])*(invA @ r))

v=u1c
d=u


import matplotlib.animation as manimation

FFMpegWriter = manimation.writers['ffmpeg']
metadata = dict(title='Corde', artist='Guillaume Couffignal',
                comment='TD4_5')
writer = FFMpegWriter(fps=24, metadata=metadata)

fig = plt.figure()
l, = plt.plot([], [], 'b')

plt.xlim(0, L)
plt.ylim(-L/2, L/2)

with writer.saving(fig, "Corde_simple.mp4", N):
    for n in range(N):
        plt.axis([0, L, -L/2, L/2])  
        l.set_data(x, u)
        plt.title("Pour : c=%s, T=%s, M=%s, N=%s, n = %s" % (c, T, M, N, n))
        u=2*(invA  @ d) - v +g*(0.5*a[n-1]+0.5*a[n+1])*(invA @ q)+g*(0.5*b[n-1]+0.5*b[n+1])*(invA @ r)
        v=copy.deepcopy(d)
        d=copy.deepcopy(u)
        writer.grab_frame()
        
end = time.time()

print("Time to run (seconds) : ", end-start)