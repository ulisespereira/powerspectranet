import numpy as np
import math as mt
from scipy.integrate import odeint as myode
import matplotlib.pyplot as plt
#parameters

phi=0.067
V1=-1.2
V2=18.
V3=12.
V4=17.4
eca=120.
ek=-84.
el=-60.
gca=4.
gk=8.
gl=2.

cm=20.

eadim=-el
gadim=gl


g1=gk/gadim
g2=gca/gadim
g3=gl/gadim
u1=ek/eadim
u2=eca/eadim
u3=el/eadim

lam=cm/gadim

a1=eadim/V4
b1=-V3/V4
a2=eadim/V2
b2=-V1/V2
c=(1/(phi*lam))
cfast=0.1*c

Izap=0.02
#I0=0.28
fm=1.
period=1600.
tstart=200.



##################################
#######functions
##################################
#dynamics

def minf(u):
	return 0.5*(1.+np.tanh(a2*u+b2))

def ninf(u):
	return 0.5*(1.+np.tanh(a1*u+b1))

def tau(u):
	return c/np.cosh((a1*u+b1)/2)

def Dminf(u):
	return 0.5*a2*np.cosh(b2+a2*u)**(-2)

def Dninf(u):
	return 0.5*a1*np.cosh(b1+a1*u)**(-2)

def Df(u):
	return g1*ninf(u)+g2*minf(u)+g3+g1*Dninf(u)*(u-u1)+g2*Dminf(u)*(u-u2)

def bet1(u):
	return -Dninf(u)

def m1(u):
	return g1*(u-u1)

def a(u,w):
	return (m1(u)*bet1(u)*tau(u)*tau(u))/(w*w*tau(u)*tau(u)+1)

def b(u,w):
	return 1-(m1(u)*bet1(u)*tau(u))/(w*w*tau(u)*tau(u)+1)

def z(u,w):
	return 1/np.sqrt(w*w*b(u,w)*b(u,w)+(w*w*a(u,w)-Df(u))**2)



def zap(t):
	if t>tstart and t<tstart+period:
		return Izap*np.sin(2.*mt.pi*(fm*(t-tstart)**2)/(2*period))
	else:
		return 0.

def field(x,t):
	field1=I0+zap(t)-g1*x[1]*(x[0]-u1)-g2*minf(x[0])*(x[0]-u2)-g3*(x[0]-u3)
	field2=(1/tau(x[0]))*(ninf(x[0])-x[1])
	return np.array([field1,field2])

####################################################
##########Plotting
###################################################

#integrating

myfreq=np.fft.fftfreq(16000,d=0.1)
t=np.linspace(0,2000,20000)

figure=plt.figure()
spectrum=figure.add_subplot(121)
dynamics=figure.add_subplot(122)
for i in range(1):
	I0=0.3*(i/10.)
	x0=np.array([-1.,ninf(-1.)])
	sol=myode(field,x0,t)
	myfft=np.fft.fft(sol[tstart*10:(tstart+period)*10,0])
	myfftNorm=np.sqrt(np.power(np.real(myfft),2)+np.power(np.imag(myfft),2))
#plotting
	spectrum.plot(myfreq,10*z(0.7,myfreq),'o')
	spectrum.plot(myfreq,myfftNorm,'o')
	dynamics.plot(t[tstart*10:(tstart+period)*10],sol[tstart*10:(tstart+period)*10,0])
	print(i)
#spectrum.set_xlim([0,0.5])
#spectrum.set_ylim([1,7])
plt.show()











#field

