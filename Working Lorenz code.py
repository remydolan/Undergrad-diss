#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 24 16:56:44 2020

@author: remydolan
"""

from numpy import *
from matplotlib.pyplot import *   
i=100001

#stepsize + upper t limit 
tmax=40
h=tmax/float(i-1)  

h1=5.962
v1=1
v2=4
def f1(t,y):
    
    f1 = np.array([h1*y[1]-v1*y[0]-y[1]*y[2],y[0]*y[2]+h1*y[0]-v2*y[1],y[0]*y[1]-y[2]])
    return f1

y0=0.1
y1=0
y2=0







#equilibria are (0,0,0),(1.392955961,4.032096673,5.616533094),(-1.392955961,-4.032096673,5.616533094)

#0,0,0
eq1=np.zeros([i,3]) 
#1.3,4,5.6
eq2=np.zeros([i,3]) 
#-1.3,-4,5.6
eq3=np.zeros([i,3]) 

ip=np.zeros([i,3]) 

fp=np.zeros([i,3])
y4=np.zeros([i,3]) 
y5=np.zeros([i,3])
y4[0]=np.array([y0,y1,y2])
y5[0]=np.array([y0,y1,y2])
tc=zeros([i])

Errest=zeros([i,3])
eq2[0:i]=np.array([1.392955961,4.032096673,5.616533094])
eq3[0:i]=np.array([-1.392955961,-4.032096673,5.616533094])
ip[0:i]=np.array([y0,y1,y2])

def phi4(t,y,f1,h):
    k1=h*f1(t,y)
    k2=h*f1(t+(h/4),y+(k1/4))
    k3=h*f1(t+(3*h/8),y+(3*k1/32)+(9*k2/32))
    k4=h*f1(t+(12*h/13),y+(1932*k1/2197)-(7200*k2/2197)+(7296*k3/2197))
    k5=h*f1(t+h,y+(439*k1/216)-(8*k2)+(3680*k3/513)-(845*k4/4104))
    k6=h*f1(t+(h/2),y-(8*k1/27)+(2*k2)-(3544*k3/2565)+(1859*k4/4104)-(11*k5/40))
    phi4=(25*k1/216)+(1408*k3/2565)+(2197*k4/4104)-(k5/5)
    return phi4

def phi5(t,y,f1,h):
    k1=h*f1(t,y)
    k2=h*f1(t+(h/4),y+(k1/4))
    k3=h*f1(t+(3*h/8),y+(3*k1/32)+(9*k2/32))
    k4=h*f1(t+(12*h/13),y+(1932*k1/2197)-(7200*k2/2197)+(7296*k3/2197))
    k5=h*f1(t+h,y+(439*k1/216)-(8*k2)+(3680*k3/513)-(845*k4/4104))
    k6=h*f1(t+(h/2),y-(8*k1/27)+(2*k2)-(3544*k3/2565)+(1859*k4/4104)-(11*k5/40))
    phi5=(16*k1/135)+(6656*k3/12825)+(28561*k4/56430)-(9*k5/50)+(2*k6/55)
    return phi5

t=0
n=1
while t<tmax:
    y4[n]=y4[n-1]+phi4(t,y4[n-1],f1,h)
    y5[n]=y5[n-1]+phi5(t,y5[n-1],f1,h)
    Errest[n]=abs(y4[n]-y5[n])

    while max(Errest[n])>0.00001:
        h=h/2
        y4[n]=y4[n-1]+phi4(t,y4[n-1],f1,h)
        y5[n]=y5[n-1]+phi5(t,y5[n-1],f1,h)
        Errest[n]=abs(y4[n]-y5[n])

        print(Errest[n],n,h,t)
    if max(Errest[n])<1e-17:
        h=2*h
    if h>=1:
        h=1
    if h<=1e-10:
        h=1e-10
    
    y4[n]=y4[n-1]+phi4(t,y4[n-1],f1,h)
    y5[n]=y5[n-1]+phi5(t,y5[n-1],f1,h)
    tc[n-1]=t
    t=t+h

    if n==i-1:
        break
    n=n+1

print(t,h,Errest[n-1],n-1,y5[n-1])
#print(z[n-2])
#print(t,h,Errest[n],n,y5[n])
#1.392955961,4.032096673,5.616533094
fp[0:i]=np.array([y5[n-1,0],y5[n-1,1],y5[n-1,2]]) 
maxus=n-1
#empty list for x coordinate of positive steady state
ssxp=zeros([maxus])
#reassigning all the zeros
ssxp[0:maxus]=1.392955961
#same for negative steady state
ssxn=zeros([maxus])
ssxn[0:maxus]=-1.392955961
ssyp=zeros([maxus])
ssyp[0:maxus]=4.032096673
ssyn=zeros([maxus])
ssyn[0:maxus]=-4.032096673
ssz=zeros([maxus])
ssz[0:maxus]=5.616533094
plot(tc[0:maxus],y5[0:maxus,0],label="How x changes with time")
plot(tc[0:maxus],ssxp)
plot(tc[0:maxus],ssxn)
annotate('jump from one orbit to another',xy=(23,3.5),xytext=(1.4,4),arrowprops=dict(facecolor='black', shrink=0.05))
legend()
show()
plot(tc[0:maxus],y5[0:maxus,1],label="How y changes with time")
plot(tc[0:maxus],ssyp)
plot(tc[0:maxus],ssyn)
legend()
show()
plot(tc[0:maxus],y5[0:maxus,2],label="How z changes with time")
plot(tc[0:maxus],ssz)
legend()
show()
from mpl_toolkits import mplot3d
import numpy as np
import matplotlib.pyplot as plt
ax=plt.axes(projection='3d')
#maxus=n-1
ax.plot3D(y5[0:maxus,0],y5[0:maxus,1],y5[0:maxus,2])
ax.plot3D(eq1[0:maxus,0],eq1[0:maxus,1],eq1[0:maxus,2],'-x')
ax.plot3D(eq2[0:maxus,0],eq2[0:maxus,1],eq2[0:maxus,2],'-x')
ax.plot3D(eq3[0:maxus,0],eq3[0:maxus,1],eq3[0:maxus,2],'-x')
ax.plot3D(ip[0:maxus,0],ip[0:maxus,1],ip[0:maxus,2],'-o')
ax.plot3D(fp[0:maxus,0],fp[0:maxus,1],fp[0:maxus,2],'-o')
#annotate('loss of accuracy, less data points',xyz=(0,0,0),xyztext=(0,0,0),arrowprops=dict(facecolor='black', shrink=0.05))
xlabel('x')
ylabel('y')
#zlabel('z')
show()
#annotate('loss of accuracy, less data points',xy=(4,-2.19),xytext=(5,-2.3),arrowprops=dict(facecolor='black', shrink=0.05))

#ax.view_init(0,180)


