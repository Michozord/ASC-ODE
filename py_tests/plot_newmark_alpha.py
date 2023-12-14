import numpy as np
import matplotlib.pyplot as plt

data_n = np.loadtxt(f'C:\ESC\ASC-ODE\ASC-ODE\py_tests\output_newmark.txt')
data_a = np.loadtxt(f'C:\ESC\ASC-ODE\ASC-ODE\py_tests\output_alpha.txt')


plt.plot(data_n[:,0], data_n[:,1],'r-',label="Newmark Method")
plt.legend()
plt.grid()
plt.show()
plt.plot(data_a[:,0], data_a[:,1],'b-',label="$\\alpha$ Method")
plt.legend()
plt.grid()
plt.show()


