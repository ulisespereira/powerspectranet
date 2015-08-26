from mynet import *
from integrator import *
from scipy.integrate import odeint
import numpy as np
import matplotlib.pyplot as plt

K=4
X0=-0.1*np.ones((K,2))
dt=0.1
Nsteps=1000
myintegrator=integrator(dt,Nsteps)

mynet=net(K)
sol=myintegrator.rk4(mynet.field,X0)

print(sol[0])
plt.plot(sol[0,:])
plt.show()
