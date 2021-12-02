#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 20 20:52:06 2020

@author: remydolan
"""

from numpy import *
from matplotlib.pyplot import *
#y(0)=1,y(1)=0 y''+y=0, solution is y(x) = cos(x) - cot(1) sin(x)
#y[0]=x(x is gradient of y),y[1]=y
i=101
tmax=1.0
h=tmax/float(i-1) 
#beta is desired BC on RHS
beta=0
#gg is guess for gradient
gg=np.arange(-10,11,1)
print(gg)

#function to be solved
def f1(t,y):
    
    f1 = np.array([-y[1],y[0]])
    return f1
#empty real solution list of zeros
rs=np.zeros([i])
#range in time
t=linspace(0,tmax,i)


#initial value for y 
y1=1

#empty approximate list of zeros
y=np.zeros([i,2])
#IC for real solution
rs[0]=np.array([1])
#yrb is a list of the y approximates at the Right Boundary
yrb=zeros([21])


#code to look for sign difference, calculates y using different gradients,
#so can get 2 boundaries where solution lies
for k in range (0,21):
  
    y[0]=np.array([gg[k],y1])
    
    for n in range (1,i): 
        y[n]=y[n-1]+h*f1(t[n-1]+(h/2),y[n-1]+(h/2)*f1(t[n-1],y[n-1]))
        
        
#creates list with all the different values of y for corresponding gradients
# from -10 to 10
    yrb[k]=y[100,0]-beta
print(yrb)

#searches list for sign difference, stops when sign difference detected
for j in range (0,21):
    if yrb[j]<0 and yrb[j+1]>=0:
        

        break
        
    if yrb[j]>0 and yrb[j+1]<=0:
        

        break
        
print(yrb[j],yrb[j+1])
print(gg[j],gg[j+1])
#assigns initial guesses for where solution lies
ydash=zeros([6])
ydash[1]=gg[j]  
ydash[0]=gg[j+1] 

#E=yestimate-boundary condition, so want E=0
#empty lists to be filled
E=zeros([6])
delkx=zeros([6])
E[0]=yrb[j+1]
print(E[0])



for k in range(1,6):
    #array for y[0] with all the different gradients ydash, want y approx equal to 0
    y[0]=np.array([ydash[k],y1])
#   for each ydash, recalculate y, get closer and closer to 0 with each iteration
    if k>4:
        break
    for n in range (1,i):
        y[n]=y[n-1]+h*f1(t[n-1]+(h/2),y[n-1]+(h/2)*f1(t[n-1],y[n-1]))
    #sec method code, from book   
    E[k]=y[i-1,1]-beta
    delkx[k-1]=-E[k]*(ydash[k]-ydash[k-1])/(E[k]-E[k-1])
    ydash[k+1]=ydash[k]+delkx[k-1]
   
    print(ydash[k],E[k],-cos(1)/sin(1))
print(-cos(1)/sin(1))
print(ydash[k])    
        
#plot of iterative scheme homing in on solution
#ki is zero line
ki=zeros([6])
plot(E,'-x') 
ylabel('E')
xlabel('k iteration')
title('iterations homing in on solution')
plot(ki)
show()

#plotting function with approximate for gradient
y[0]=np.array([ydash[k],y1])
for n in range (1,i): 
    y[n]=y[n-1]+h*f1(t[n-1]+(h/2),y[n-1]+(h/2)*f1(t[n-1],y[n-1]))
    rs[n]=cos(t[n])-cos(1)*sin(t[n])/sin(1)
tmin=0
plot(t[tmin:i],y[tmin:i,1],label="Modified Euler approximate")
plot(t[tmin:i],rs[tmin:i],label="Real solution")
title("Approximate compared to real solution(analytic gradient)")
ylabel("y")
xlabel("t")
legend()
show()

#plot(log(E),'-x') 
#ylabel('E')
#xlabel('k iteration')
#title('iterations homing in on solution')
#plot(ki)
#show()
#print(E)
