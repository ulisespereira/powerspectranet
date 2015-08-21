from scipy import signal
import matplotlib.pyplot as plt
import numpy as np
from morrislecar import *
from scipy.integrate import odeint
import timeit




#parameters simulations
dt=1e-2
fs=1/dt
N=2e4
time=np.arange(N)/fs

#parameters stimuli
amp=0.1
freq=0.
tstart=100.
nstart=tstart/dt
I0=0.
def mysin(x):
	if x>tstart:
		return amp*np.sin(2*np.pi*freq*x)+I0
	else:
		return 0


#intializing the model
myneuron=morrislecar()
myneuron.thestim(mysin)
x0=np.array([-1.,myneuron.ninf(-1.)])
power_spectra=[]
mean_v=[]

freq0=0.01
n_sim=50
freq=freq0
delta_freq=0.05
for i in delta_freq*np.arange(n_sim):
	current_stim=np.array([mysin(t) for t in time[nstart:-1]])
	sol=odeint(myneuron.field,x0,time)
	plt.plot(time[nstart:-1],sol[nstart:-1,0])
	plt.show()


	f_v,P_v=signal.periodogram(sol[nstart:-1,0],fs,'flattop',scaling='spectrum')
	f_I,P_I=signal.periodogram(current_stim,fs,'flattop',scaling='spectrum')
	
	P_max_v=np.amax(P_v)
	P_max_I=np.amax(P_I)
#	print(f_v[P_v==P_max_v])
#	print(f_I[P_I==P_max_I])
#	print(freq)
	power_spectra.append(P_max_v/P_max_I)

	#plt.plot(f,Pxx)
	#plt.xlim([0,10])	
	#plt.show()
	freq=freq+i
#	myneuron.thestim(mysin)
#	print(i)

thefreq=(delta_freq*np.arange(n_sim)+freq0)
plt.plot(thefreq,power_spectra,'or')
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

