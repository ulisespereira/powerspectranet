import numpy as np
import math as mt

""" 
This module gives you  the Morris Lecar neuron model 
and its power spectra

"""


class morrislecar:
	def __init__(self):
		#parameters
		self.phi=0.067
		self.V1=-1.2
		self.V2=18.
		self.V3=12.
		self.V4=17.4
		self.eca=120.
		self.ek=-84.
		self.el=-60.
		self.gca=4.
		self.gk=8.
		self.gl=1.
		self.cm=20.
		#nondimensionalization
		self.eadim=-self.el
		self.gadim=self.gl

		self.g1=self.gk/self.gadim
		self.g2=self.gca/self.gadim
		self.g3=self.gl/self.gadim
		self.u1=self.ek/self.eadim
		self.u2=self.eca/self.eadim
		self.u3=self.el/self.eadim
		self.lam=self.cm/self.gadim #real time step
		self.a1=self.eadim/self.V4
		self.b1=-self.V3/self.V4
		self.a2=self.eadim/self.V2
		self.b2=-self.V1/self.V2
		self.c=(1/(self.phi*self.lam))
		self.connections=0
		self.f=lambda x: 0.                # stimuli funct

		#functions

	def minf(self,u):
		return 0.5*(1.+np.tanh(self.a2*u+self.b2))

	def ninf(self,u):
		return 0.5*(1.+np.tanh(self.a1*u+self.b1))

	def tau(self,u):
		return self.c/np.cosh((self.a1*u+self.b1)/2)
	
	def Dminf(self,u):
		return 0.5*self.a2*np.cosh(self.b2+self.a2*u)**(-2)

	def Dninf(self,u):
		return 0.5*self.a1*np.cosh(self.b1+self.a1*u)**(-2)

	def Df(self,u):
		return self.g1*self.ninf(u)+self.g2*self.minf(u)+self.g3+self.g1*self.Dninf(u)*(u-self.u1)+self.g2*self.Dminf(u)*(u-self.u2)+self.connections

	def bet1(self,u):
		return -self.Dninf(u)

	def m1(self,u):
		return self.g1*(u-self.u1)

	def a(self,u,w):
		return (self.m1(u)*self.bet1(u)*self.tau(u)*self.tau(u))/(w*w*self.tau(u)*self.tau(u)+1)

	def b(self,u,w):
		return 1+(self.m1(u)*self.bet1(u)*self.tau(u))/(w*w*self.tau(u)*self.tau(u)+1)

	def z(self,u,w):
		return 1/(w*w*self.b(u,w)*self.b(u,w)+(w*w*self.a(u,w)-self.Df(u))**2)
		
		# setting up the stimulation function
	def thestim(self,elstim):
		self.f=lambda x: elstim(x)

	def field(self,x,t):
		field1=self.f(t)-self.g1*x[1]*(x[0]-self.u1)-self.g2*self.minf(x[0])*(x[0]-self.u2)-self.g3*(x[0]-self.u3)
		field2=(1/self.tau(x[0]))*(self.ninf(x[0])-x[1])
		return np.array([field1,field2])

