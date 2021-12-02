#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sat Mar 21 16:26:13 2020

@author: remydolan
"""

#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 16 22:35:03 2020

@author: remydolan
"""
from numpy import *
from matplotlib.pyplot import * 
Rmax=1e-2
x0=0
y0=1
i=101
y4=zeros([i]) 
y5=zeros([i])
y4[0]=y0
y5[0]=y0
tmax=10
x=linspace(2,12,i)
h=tmax/float(i-1)
Errest=zeros([i])
rs=zeros([i])
rs[0]=1
xc=zeros([i])
xc[0]=2
def f1(x,y):
    f1 =-x*y**2
    return f1

#def phi(x,y,f1,h):
#    k1=f1(x,y)
#    k2=f1(x+h/2,y+k1*h/2)
#    k3=f1(x+h/2,y+h*k2/2)
#    k4=f1(x+h,y+h*k3) 
#    phi=(k1+(2*k2)+(2*k3)+k4)/6 
#    return phi
def phi4(x,y,f1,h):
    k1=h*f1(x,y)
    k2=h*f1(x+(h/4),y+(k1/4))
    k3=h*f1(x+(3*h/8),y+(3*k1/32)+(9*k2/32))
    k4=h*f1(x+(12*h/13),y+(1932*k1/2197)-(7200*k2/2197)+(7296*k3/2197))
    k5=h*f1(x+h,y+(439*k1/216)-(8*k2)+(3680*k3/513)-(845*k4/4104))
    k6=h*f1(x+(h/2),y-(8*k1/27)+(2*k2)-(3544*k3/2565)+(1859*k4/4104)-(11*k5/40))
    phi4=(25*k1/216)+(1408*k3/2565)+(2197*k4/4104)-(k5/5)
    return phi4

def phi5(x,y,f1,h):
    k1=h*f1(x,y)
    k2=h*f1(x+(h/4),y+(k1/4))
    k3=h*f1(x+(3*h/8),y+(3*k1/32)+(9*k2/32))
    k4=h*f1(x+(12*h/13),y+(1932*k1/2197)-(7200*k2/2197)+(7296*k3/2197))
    k5=h*f1(x+h,y+(439*k1/216)-(8*k2)+(3680*k3/513)-(845*k4/4104))
    k6=h*f1(x+(h/2),y-(8*k1/27)+(2*k2)-(3544*k3/2565)+(1859*k4/4104)-(11*k5/40))
    phi5=(16*k1/135)+(6656*k3/12825)+(28561*k4/56430)-(9*k5/50)+(2*k6/55)
    return phi5

#def Errest(x,y,f1,h):
#    k1=h*f1(x,y)
#    k2=h*f1(x+(h/4),y+(k1/4))
#    k3=h*f1(x+(3*h/8),y+(3*k1/32)+(9*k2/32))
#    k4=h*f1(x+(12*h/13),y+(1932*k1/2197)-(7200*k2/2197)+(7296*k3/2197))
#    k5=h*f1(x+h,y+(439*k1/216)-(8*k2)+(3680*k3/513)-(845*k4/4104))
#    k6=h*f1(x+(h/2),y-(8*k1/27)+(2*k2)-(3544*k3/2565)+(1859*k4/4104)-(11*k5/40))
#    Errest=(k1/360)-(128*k3/4275)-(2097*k4/75240)+(k5/50)+(2*k6/55)
#    return Errest
x=2
n=1
while x<10:
    y4[n]=y4[n-1]+phi4(x,y4[n-1],f1,h)
    y5[n]=y5[n-1]+phi5(x,y5[n-1],f1,h)
    Errest[n]=abs(y4[n]-y5[n])

    while Errest[n]>1e-5:
        h=h/2
        y4[n]=y4[n-1]+phi4(x,y4[n-1],f1,h)
        y5[n]=y5[n-1]+phi5(x,y5[n-1],f1,h)
        Errest[n]=abs(y4[n]-y5[n])
        print(Errest[n],n,h,x)
    if Errest[n]<1e-10:
        h=2*h
    if h>=1:
        h=1
    if h<=1e-10:
        h=1e-10
    
    y4[n]=y4[n-1]+phi4(x,y4[n-1],f1,h)
    y5[n]=y5[n-1]+phi5(x,y5[n-1],f1,h)
    
    rs[n-1]=2/(x**2-2)
  
    xc[n-1]=x
    x=x+h
    if n==i-1:
        break
    n=n+1
    
#print(y5[1],rs[1],Errest[1],1)
xmin=0
xmax=n-1
plot(xc[xmin:xmax],rs[xmin:xmax],linestyle='solid',label="Analytic Solution")
plot(xc[xmin:xmax],y5[xmin:xmax],'x',label="Approximate")
title("Plot of approximate solution of f1 using RKF")
ylabel("y")
xlabel("x")
legend()
show()
#plot(z[xmin:n-1],log10(abs((y4[xmin:n-1]-rs[xmin:n-1])/rs[xmin:n-1])), label="RKF")
#title("log of absolute difference between approximate and real solution")
#ylabel("log(abs(y-rs)/rs)")
#xlabel("x")
#annotate('loss of accuracy, less data points',xy=(4,-2.19),xytext=(5,-2.3),arrowprops=dict(facecolor='black', shrink=0.05))
#legend()
#show()

#plot(Ratio[0:100])
#print(rs,z,y)
print(xc,rs,y5,n)       
#print(z[n],rs[n],y5[n])
