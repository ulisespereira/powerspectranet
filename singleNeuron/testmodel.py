import numpy as np
import math as mt

""" 
This module gives you  the Morris Lecar neuron model 
and its power spectra

"""


class testmodel:
	def __init__(self):
		#parameters
		self.tau=1.0

		self.f=lambda x: 0.                # stimuli funct

		#functions

	def z(self,w):
		return (self.tau*self.tau)/(w*w*self.tau*self.tau+1)
		
		# setting up the stimulation function
	def thestim(self,elstim):
		self.f=lambda x: elstim(x)

	def field(self,x,t):
		return self.f(t)-(1/self.tau)*x

