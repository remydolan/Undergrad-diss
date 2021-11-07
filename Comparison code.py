#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sun Nov 24 19:21:44 2019

@author: remydolan
"""

from numpy import *
from matplotlib.pyplot import *
#1st order ode
#exact solution=x+2exp(-x)-1, y=2/e when x=1
#dy/dx+y=x,y(0)=1
#dy/dx=x+y,y(0)=1, es=-x+2exp(x)-1
def f1(x,y):
    f1 = x+y
    return f1

#initial conditions
x0=0 
y0=1
#upper x limit
xmax=10
#number of steps
i=10001 
#step size
h=xmax/float(i-1)

#empty lists to be filled
ye=zeros([i])
yme=zeros([i])
yrk=zeros([i])
rs=zeros([i])
rs[0]=0
ye[0]=y0
yme[0]=y0
yrk[0]=y0
#phi for RK
def phi(x,y,f1,h):
    k1=f1(x,y)
    k2=f1(x+h/2,y+k1*h/2)
    k3=f1(x+h/2,y+h*k2/2)
    k4=f1(x+h,y+h*k3) 
    phi=(k1+(2*k2)+(2*k3)+k4)/6 
    return phi
#xrange 
x=linspace(0,xmax,i)

for n in range (1,i):
    #euler method
    ye[n]=ye[n-1]+h*f1(x[n-1],ye[n-1])
    #exact solution
    rs[n]=-x[n]+2*exp(x[n])-1
    #mod euler method
    yme[n]=yme[n-1]+h*f1(x[n-1]+(h/2),yme[n-1]+(h/2)*f1(x[n-1],yme[n-1]))
    #RK method
    yrk[n]=yrk[n-1]+h*phi(x[n-1],yrk[n-1],f1,h)
    

xmin=8

#plots
title("log10 of absolute difference")
plot(x[xmin:i],log10(abs((ye[xmin:i]-rs[xmin:i])/rs[xmin:i])), label="Euler")
plot(x[xmin:i],log10(abs((yme[xmin:i]-rs[xmin:i])/rs[xmin:i])),label="Modified Euler")
plot(x[xmin:i],log10(abs((yrk[xmin:i]-rs[xmin:i])/rs[xmin:i])), label="Runge-Kutta")

legend()
#title("log of absolute difference between approximate and real solution")
xlabel("x")
ylabel("log10(abs(y-rs)/rs)")
#xlim(0,100)
#ylim(-60,0) 
show()

#plot(x[xmin:i],ye[xmin:i],'b', label="Euler")
#plot(x[xmin:i],yme[xmin:i],'go',label="Modified Euler")
#plot(x[xmin:i],yrk[xmin:i],'rx', label="Runge-Kutta")
#plot(x[xmin:i],rs[xmin:i],'k', label="Analytic solution")
plot(x[xmin:i],ye[xmin:i],'b', label="Euler")
plot(x[xmin:i],yme[xmin:i],'y--',label="Modified Euler")
plot(x[xmin:i],yrk[xmin:i],'r-.', label="Runge-Kutta")
plot(x[xmin:i],rs[xmin:i],'', label="Analytic solution")
title("Plot of approximate solutions of f1")
ylabel("y")
xlabel("x")
#xlim(0,100)
#ylim(-60,0)
legend()
show()