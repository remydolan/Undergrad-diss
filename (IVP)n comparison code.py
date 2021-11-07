#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sun Jan 19 19:18:42 2020

@author: remydolan
"""

### pair od diff equations derived from x''=x+t
#x(0)=1,x'(0)=-1
from numpy import *
from matplotlib.pyplot import *
#number of steps

i=1001

#stepsize + upper t limit
tmax=10
h=tmax/float(i-1)

#2nd order ode to be solved, x''=x+t, (y,x+t)
def f1(t,y):
    
    f1 = np.array([y[1],y[0]+t])
    return f1
#empty list for real solution
rs=np.zeros([i])

#trange
t=linspace(0,tmax,i)
#initial conditions
tnot=0 
y0=1
y1=-1 
#empty approximate list
#y=np.zeros([i,2])
ye=zeros([i,2])
yme=zeros([i,2])
yrk=zeros([i,2])
#initial conditions
#y[0]=np.array([y0,y1])
rs[0]=np.array([tnot])
ye[0]=np.array([y0,y1])
yme[0]=np.array([y0,y1])
yrk[0]=np.array([y0,y1])

def phi(x,y,f1,h):
    k1=f1(x,y)
    k2=f1(x+h/2,y+k1*h/2)
    k3=f1(x+h/2,y+h*k2/2)
    k4=f1(x+h,y+h*k3) 
    phi=(k1+(2*k2)+(2*k3)+k4)/6 
    return phi
#ivp^n, using modified euler
for n in range (1,i): 
    ye[n]=ye[n-1]+h*f1(t[n-1],ye[n-1])
    yme[n]=yme[n-1]+h*f1(t[n-1]+(h/2),yme[n-1]+(h/2)*f1(t[n-1],yme[n-1]))
    yrk[n]=yrk[n-1]+h*phi(t[n-1],yrk[n-1],f1,h)
    rs[n]=-t[n]+0.5*(exp(t[n])+exp(-t[n]))
   
   
#plots
tmin=1
#n=np.arange(tmin,i,1)
rs[n]=-t[n]+0.5*(exp(t[n])+exp(-t[n]))
plot(t[tmin:i],log(abs((ye[tmin:i,0]-rs[tmin:i])/rs[tmin:i])),label="Euler")
plot(t[tmin:i],log(abs((yme[tmin:i,0]-rs[tmin:i])/rs[tmin:i])),label="Modified Euler")
plot(t[tmin:i],log(abs((yrk[tmin:i,0]-rs[tmin:i])/rs[tmin:i])),label="Runge-Kutta")
#title("plot of the absolute difference between y approximate and real solution for 2nd order ode")
xlabel("t")
ylabel("log(abs(y-rs)/rs)")
legend()
show()
#plot(t[tmin:i],ye[tmin:i,0],label="Euler approximate")
plot(t[tmin:i],yme[tmin:i,0],label="Modified Euler approximate")
#plot(t[tmin:i],yrk[tmin:i,0],label="Runge-Kutta approximate")
plot(t[tmin:i],rs[tmin:i],label="Real solution")
title("Approximate compared to real solution")
ylabel("y")
xlabel("t")
legend()
show()

