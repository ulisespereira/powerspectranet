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
		return amp*np.sin(2*np.pi*freq*x)+I0
	else:
		return 0

#intializing the model
#morris lecar
myneuron=morrislecar()
myneuron.thestim(mysin)
x0=np.array([-1.,myneuron.ninf(-1.)])
#test for the impedance spectra analysis
mytest=testmodel()
mytest.thestim(mysin)
x0_test=0.

#initializing the simulation parameters
power_spectra=[] #here comes the power spectra
mean_v=[] #the mean voltage
dt=1e-2   
fs=1/dt
amp=0.1 #ampplitude of stimulation
tstart=100.  #begining of stimulation
nstart=tstart/dt 
I0=0. #dc
freq0=0.01
freq=freq0
delta_freq=0.01
n_sim=50
for i in delta_freq*np.arange(n_sim):
	N=nstart+10.*(1/freq)*fs
	time=np.arange(N)/fs
	current_stim=np.array([mysin(t) for t in time[nstart:-1]])
	#sol=odeint(myneuron.field,x0,time)
	sol_test=odeint(mytest.field,x0,time)
	plt.plot(time[nstart:-1],sol_test[nstart:-1])
	plt.show()
	#mean_v.append(np.mean(sol[nstart:-1,0]))
	mean_v.append(np.mean(sol_test[nstart:-1]))

	#f_v,P_v=signal.periodogram(sol[nstart:-1,0],fs,'flattop',scaling='spectrum')
	f_v,P_v=signal.periodogram(sol_test[nstart:-1],fs,'flattop',scaling='spectrum')
	f_I,P_I=signal.periodogram(current_stim,fs,'flattop',scaling='spectrum')
	
	P_max_v=np.amax(P_v)
	P_max_I=np.amax(P_I)
#	print(f_v[P_v==P_max_v])
#	print(f_I[P_I==P_max_I])
#	print(freq)
	print(P_max_v,P_max_I)
	power_spectra.append(P_max_v/P_max_I)

#	plt.plot(P_v)
#	plt.xlim([0,10])	
	plt.show()
	freq=freq+i
#	myneuron.thestim(mysin)
	print(i)

the_mean_v=np.mean(mean_v)
thefreq=(delta_freq*np.arange(n_sim)+freq0)
plt.plot(thefreq,power_spectra,'or')
#plt.plot(thefreq,np.array([myneuron.z(the_mean_v,2*np.pi*f) for f in thefreq]),'b')
plt.plot(thefreq,np.array([mytest.z(2*np.pi*f) for f in thefreq]),'b')
plt.show()
##integrating


#	x0=np.array([-1.,ninf(-1.)])
#	myfft=np.fft.fft(sol[tstart*10:(tstart+period)*10,0])
#	myfftNorm=np.sqrt(np.power(np.real(myfft),2)+np.power(np.imag(myfft),2))
##plotting
#	spectrum.plot(myfreq,10*z(0.7,myfreq),'o')
#	spectrum.plot(myfreq,myfftNorm,'o')
#	dynamics.plot(t[tstart*10:(tstart+period)*10],sol[tstart*10:(tstart+period)*10,0])
#	print(i)
##spectrum.set_xlim([0,0.5])
##spectrum.set_ylim([1,7])
#plt.show()
#










#field

