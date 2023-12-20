import numpy as np
import matplotlib.pyplot as plt

data_ie = np.loadtxt(f'C:\ESC\ASC-ODE\ASC-ODE\py_tests\output_ie.txt')
data_ee = np.loadtxt(f'C:\ESC\ASC-ODE\ASC-ODE\py_tests\output_ee.txt')
# data_cn = np.loadtxt(f'C:\ESC\ASC-ODE\ASC-ODE\py_tests\output_cn.txt')
data_rk = np.loadtxt(f'C:\ESC\ASC-ODE\ASC-ODE\py_tests\output_rk.txt')
data_rad = np.loadtxt(f'C:\ESC\ASC-ODE\ASC-ODE\py_tests\output_rad.txt')

plt.plot(data_ie[:,0], data_ie[:,1],'r-',label="Implicit Euler")
plt.plot(data_ee[:,0], data_ee[:,1],'b-',label="Explicit Euler")
# plt.plot(data_cn[:,0], data_cn[:,1],'g-',label="Crank-Nicolson")
plt.plot(data_rk[:,0], data_rk[:,1],'c-',label="RK2/impl. midpoint")
plt.plot(data_rad[:,0], data_rad[:,1],'m-',label="Radau IIA")
plt.legend()
plt.grid()
plt.show()


