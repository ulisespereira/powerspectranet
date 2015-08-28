from mynet import *
from scipy import signal
from integrator import *
from scipy.integrate import odeint
import numpy as np
import matplotlib.pyplot as plt

K=10
X0=-1*np.ones((K,2))
dt=0.1
fs=1/dt
Nsteps=10000
myintegrator=integrator(dt,Nsteps)
mynet=net(K)

#stimulation function
def mysin(x):
	if x>tstart and x<=period+tstart:
		return amp*np.sin(2*np.pi*freq*(x-tstart)**2/(2*period))+I0
	else:
		return I0


amp=0.1 #ampplitude of stimulation
tstart=100.  #begining of stimulation
period=400.
freq=0.4
I0=0. #dc
time=np.arange(Nsteps)*dt
nstart=int(tstart/dt)
nstop=int((tstart+period)/dt)


mynet.stim(mysin)
sol=myintegrator.rk4(mynet.field,X0)
current_stim=np.array([mysin(t) for t in time[nstart:nstop]])

mean_voltages=[]
power_spectras=[]
for i in range(0,K):	
	f_v,P_v=signal.periodogram(sol[nstart:nstop,i],fs,scaling='spectrum')
	f_I,P_I=signal.periodogram(current_stim,fs,scaling='spectrum')
	mean_v=(np.mean(sol[nstart:-1,i]))
	mean_voltages.append(mean_v)
	power_spectras.append(P_v/P_I)


figure=plt.figure()
zap=figure.add_subplot(121)
spectra=figure.add_subplot(122)

for i in range(0,K):
	zap.plot(sol[:,i])
for i in range(0,K):
	spectra.plot(f_v,power_spectras[i])

spectra.set_xlabel('Normalized Freq')
spectra.set_ylabel('Normalized Impedance')
spectra.set_xlim([0.01,0.4])	
spectra.set_ylim([0.1,1.8])	
plt.show()
