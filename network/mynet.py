from morrislecar import *
import numpy as np

''' This is a class  network for a network of Morris-Lecar neurons
This class is very limmited since there is no heterogeneity in the neurons
and there is stimulation to just one neuron.
'''

class net:
	def __init__(self,K):
		self.K=K
		self.f=0.1*np.ones(self.K)
		self.g=0.1*np.ones(self.K)
		self.connectivity=np.outer(self.f,self.g)
		self.myneuron=morrislecar()
		self.f=lambda x: 0.

	def stim(self,mystim):
		self.f=lambda x: mystim(x)

	def field(self,X,t):
		# X is an Kx2 array where the first two are u and x
		field=np.array([self.myneuron.field(x,t) for x in X])
		field[0,0]=field[0,0]+self.f(t) #stimulation
		field[:,0]=field[:,0]+np.dot(self.connectivity,X[:,0])
		return field

		
	