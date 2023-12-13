import numpy as np
import matplotlib.pyplot as plt

data_ie = np.loadtxt(f'C:\ESC\ASC-ODE\ASC-ODE\py_tests\output_ie.txt')
data_ee = np.loadtxt(f'C:\ESC\ASC-ODE\ASC-ODE\py_tests\output_ee.txt')


plt.plot(data_ie[:,0], data_ie[:,1],'r-',label="Implicit Euler")
plt.plot(data_ee[:,0], data_ee[:,1],'b-',label="Explicit Euler")
plt.legend()
plt.grid()
plt.show()


