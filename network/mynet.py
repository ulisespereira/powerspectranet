from morrislecar import *
import numpy as np

''' This is a class  network for a network of Morris-Lecar neurons
This class is very limmited since there is no heterogeneity in the neurons
and there is stimulation to just one neuron.
'''

class net:
	def __init__(self,K):
		self.K=K
		self.f=(1.5/K)*np.ones(self.K)
		self.g=(1.5/K)*np.ones(self.K)
		self.connectivity=np.outer(self.f,self.g)
		self.myneuron=morrislecar()
		self.thestim=lambda x: 0.

	def stim(self,mystim):
		self.thestim=lambda x: mystim(x)

	def field(self,X,t):
		# X is an Kx2 array where the first two are u and x
		field=np.array([self.myneuron.field(x,t) for x in X])
	#	print(field)
		field[0,0]=field[0,0]+self.thestim(t) #stimulation
		field[:,0]=field[:,0]+np.dot(self.connectivity,X[:,0])
	#	print(field)
		return field

	def z2(self,u,w):
		impedance_spectra=[]
		SNz2=[]
		SNzre=[]
		SNzim=[]
		ReZ=0.
		ImZ=0.
		for i in range(self.K):
			self.myneuron.synapses=self.f[i]*(np.sum(self.g))
			SNz2.append(self.myneuron.z2(u[i],w))
			SNzre.append(self.myneuron.zre(u[i],w))
			SNzim.append(self.myneuron.zim(u[i],w))
			ReZ=ReZ+SNzre[i]*self.f[i]*self.g[i]
			ImZ=ImZ+SNzim[i]*self.f[i]*self.g[i] 
		ReZ0=ReZ-SNzre[0]*self.f[0]*self.g[0]
		ImZ0=ImZ-SNzim[0]*self.f[0]*self.g[0]
		denominator=1+2*ReZ+np.multiply(ReZ,ReZ)+np.multiply(ImZ,ImZ)
		numerator=1+2*ReZ0+np.multiply(ReZ0,ReZ0)+np.multiply(ImZ0,ImZ0)
		
		#impedance_spectra.append(np.multiply(SNz2[0],np.multiply(numerator1,1./denominator)))
		impedance_spectra.append(np.multiply(SNz2[0],np.multiply(numerator,1/denominator)))
		for i in range(1,self.K):
			impedance_spectra.append(np.multiply(SNz2[0]*self.g[0]**2,np.multiply(SNz2[0]*self.f[i]**2,1./denominator)))
		return impedance_spectra
		
	
