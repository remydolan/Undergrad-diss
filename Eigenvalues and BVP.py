#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sun Feb  2 14:56:32 2020

@author: remydolan
"""

from numpy import *
from matplotlib.pyplot import *

mrange=200
i=101
tmax=1.0
h=tmax/float(i-1) 
beta=0
#range of guesses for lamda^2
lamdag=np.arange(0,mrange,1)



def f1(t,y):
    
    f1 = np.array([y[1],-y[2]*y[0],0])
    return f1


t=linspace(0,tmax,i)


 
y1i=1
y0i=0
node=4
#variable to be changed is y2

y=np.zeros([i,3])



yrb=zeros([mrange])


#code to look for sign difference, calculates y using different guesses for lamda,
#so can get 2 boundaries where solution lies
for k in range (0,mrange):
  
    y[0]=np.array([y0i,y1i,lamdag[k]])
    
    for n in range (1,i): 
        y[n]=y[n-1]+h*f1(t[n-1]+(h/2),y[n-1]+(h/2)*f1(t[n-1],y[n-1]))
        
        
#creates list with all the different values of y for corresponding lamdag
# from 0 to 200
    yrb[k]=y[i-1,0]-beta


#searches list for sign difference, stops when sign difference detected
k=0
for j in range (0,mrange):
    
    if yrb[j]<0 and yrb[j+1]>=0:
        k=k+1 
#        print pr[j],a[j],pr[j+1],a[j+1],y1i,"-ve first"
#        break
        
    if yrb[j]>0 and yrb[j+1]<=0:
        k=k+1
    if k==node:
        break
#        print pr[j],a[j],pr[j+1],a[j+1],y1i,"positive first"
#        break
        

#assigns initial guesses for where solution lies
l=zeros([6])
l[1]=lamdag[j]  
l[0]=lamdag[j+1] 


#E=yestimate-boundary condition, so want E=0
#empty lists to be filled
E=zeros([6])
delkx=zeros([6])
E[0]=yrb[j+1]




for k in range(1,6):
    #array for y[0] with all the different values of l(lamda^2), want y approx equal to 0
    y[0]=np.array([y0i,y1i,l[k]])
#   for each l, recalculate y, get closer and closer to 0 with each iteration
    if k>4:
        break
    for n in range (1,i):
        y[n]=y[n-1]+h*f1(t[n-1]+(h/2),y[n-1]+(h/2)*f1(t[n-1],y[n-1]))
#    #sec method code   
    E[k]=y[i-1,0]-beta
    delkx[k-1]=-E[k]*(l[k]-l[k-1])/(E[k]-E[k-1])
    l[k+1]=l[k]+delkx[k-1]
   
    print(l[k],E[k])
print(pi*pi*node*node)




ki=zeros([6])
plot(E,'-x') 
ylabel('E')
xlabel('k iteration')
title('iterations homing in on solution')
plot(ki)
show()
#plotting function with approximate for lamda^2
y[0]=np.array([y0i,y1i,l[k]])
for n in range (1,i): 
    y[n]=y[n-1]+h*f1(t[n-1]+(h/2),y[n-1]+(h/2)*f1(t[n-1],y[n-1]))
tmin=0
plot(t[tmin:i],y[tmin:i,0],label="Modified Euler approximate")
title("BVP plot")
ylabel("y")
xlabel("t")
legend()
show()
