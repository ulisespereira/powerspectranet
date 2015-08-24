from scipy import signal
import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import odeint
from morrislecar import *
from testmodel import *



''' This code is to calculate the power spectra of a given model 
    There are  an stimulation function which is a sin wave of a given frequencty
    after stimulate for a certain number of cycles of a single mode the 
    power spectra is calculated and extracted the maximum. This is 
    repeated in the loop to create a complete power spectra'''



#stimulation function
def mysin(x):
	if x>tstart:
		return amp*np.sin(2*np.pi*freq*(x-tstart)**2/(2*period))+I0
	else:
		return I0

#intializing the model
#morris lecar
myneuron=morrislecar()
myneuron.thestim(mysin)
x0=np.array([-1.,myneuron.ninf(-1.)])

#initializing the simulation parameters
dt=1e-2   
fs=1/dt
amp=0.005 #ampplitude of stimulation
tstart=100.  #begining of stimulation
nstart=tstart/dt 
I0=0. #dc
period=3000
freq=0.4
delta_freq=0.01
n_sim=1




N=nstart+period/dt
time=np.arange(N)/fs
figure=plt.figure()
zap=figure.add_subplot(121)
spectra=figure.add_subplot(122)
thefreq=0.001*np.arange(400)
for I in 0.05*np.arange(5):
	I0=I
	current_stim=np.array([mysin(t) for t in time[nstart:-1]])
	
	sol=odeint(myneuron.field,x0,time)
	#sol_test=odeint(mytest.field,x0,time)
	zap.plot(myneuron.lam*time[nstart:-1],myneuron.eadim*sol[nstart:-1,0])

	#mean_v.append(np.mean(sol[nstart:-1,0]))
	mean_v=(np.mean(sol[nstart:-1,0]))

	f_v,P_v=signal.periodogram(sol[nstart:-1,0],fs,scaling='spectrum')
	f_I,P_I=signal.periodogram(current_stim,fs,scaling='spectrum')
	
	theo_curve=np.array([myneuron.z(mean_v,2*np.pi*f) for f in thefreq])

	spectra.plot(f_v,P_v/P_I,'o')
	spectra.plot(thefreq,theo_curve,'b')
	print(I)
zap.set_xlabel('ms')
zap.set_ylabel('mV')
spectra.set_xlabel('Normalized Freq')
spectra.set_ylabel('Normalized Impedance')
spectra.set_xlim([0.01,0.4])	
spectra.set_ylim([0.1,2])	
plt.show()











#field

