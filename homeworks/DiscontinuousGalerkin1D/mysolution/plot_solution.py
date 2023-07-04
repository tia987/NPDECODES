from matplotlib.pyplot import figure, plot, savefig, xlabel, ylabel
from numpy import genfromtxt
from sys import argv
from mpl_toolkits import mplot3d
import matplotlib.pyplot as plt
import numpy as np

input_file = str(argv[1])
output_file = str(argv[2])

data = genfromtxt(input_file, delimiter=',')
x_t = data[0]
N = 10
x_t = np.tile(x_t, (N,1))
x = x_t.reshape(-1)
z_t = data[1:]
z = z_t.reshape(-1)
 
y = np.linspace(0, 5, 10)
y = np.repeat(y,81)

fig = plt.figure()
ax = plt.axes(projection='3d')

ax.scatter3D(x, y, z, c=z, cmap='Blues')
ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_zlabel('z')
savefig(output_file)
print('Generated ' + output_file)
