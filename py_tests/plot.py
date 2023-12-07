import numpy as np
import matplotlib.pyplot as plt

data = np.loadtxt(f'C:\ESC\ASC-ODE\ASC-ODE\py_tests\output.txt')


plt.plot(data[:,0], data[:,1],'r--',label="x")
plt.plot(data[:,0], data[:,2],'b--',label="y")
plt.legend()
plt.grid()
plt.show()

plt.plot(data[:,1], data[:,2])
plt.show()

