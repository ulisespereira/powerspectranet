import numpy as np
import math as mt
from scipy.integrate import odeint as myode
import matplotlib.pyplot as plt
import scipy.fftpack

#parameters




tstart=0.
freq=50.
period=400.
Izap=0.1
def zap(t):
	if t>tstart and t<tstart+period:
		return Izap*np.sin(2.*mt.pi*freq*t)
	else:
		return 0.

def field(x,t):
	field1=zap(t)-x
	return field1


#number of sample points
N=800
#sample spacing
dt=1/1000.
#
t=np.linspace(0,N*dt,N)
#y=10.+np.sin(50.0*2.*mt.pi*t)+0.5*np.sin(80.0*2.0*np.pi*t)
x0=0
sol=myode(field,x0,t)
yf = scipy.fftpack.fft(sol)
xf = np.linspace(0.0, 1.0/(2.0*dt), N/2)




figure=plt.figure()
spectrum=figure.add_subplot(121)
dynamics=figure.add_subplot(122)
spectrum.plot(xf[0:N/2],(2.0/N)*(np.abs(yf[0:N/2])))
dynamics.plot(t,sol)

#x0=0.
#freq=freq+i*1.
#myfft=np.fft.fft(sol[tstart:(tstart+period)])
#myfftNorm=np.sqrt(np.power(np.real(myfft),2)+np.power(np.imag(myfft),2))
#spectrum.plot(myfftNorm,'o')
#dynamics.plot(t[tstart*10:(tstart+period)*10],sol[tstart*10:(tstart+period)*10,0])
#print(sol)
#spectrum.set_xlim([0.1,1000])
	#spectrum.set_ylim([1,7])
plt.show()











#field

