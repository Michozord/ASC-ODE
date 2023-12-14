import numpy as np
import matplotlib.pyplot as plt

data = np.loadtxt(f'C:\ESC\ASC-ODE\ASC-ODE\py_tests\output_alpha_2.txt')

plt.plot(data[:,1], data[:,2],'b-',label="Masse 1")
plt.plot(data[:,3], data[:,4],'r-',label="Masse 2")
plt.legend()
plt.grid()
plt.show()


