import numpy as np

""" this is a integrator """


class integrator:
	def __init__(self,dt,nsteps):
		self.dt=dt
		self.nsteps=nsteps
	def rk4(self,field,x0):
		myx=x0
		mysol=[myx.flatten('F')]
		for t in np.arange(self.nsteps)*self.dt:
			k1=field(myx,t)
			k2=field(myx+(self.dt/2)*k1,t+self.dt/2)
			k3=field(myx+(self.dt/2)*k2,t+self.dt/2)
			k4=field(myx+self.dt*k3,t+self.dt)
			myx=myx+(self.dt/6)*(k1+2*k2+2*k3+k4)
			mysol.append(myx.flatten('F'))
			print t
		return np.array(mysol)

