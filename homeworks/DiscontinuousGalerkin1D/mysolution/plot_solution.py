from matplotlib.pyplot import figure, plot, savefig, xlabel, ylabel
from numpy import genfromtxt
from sys import argv
from mpl_toolkits import mplot3d

input_file = str(argv[1])
output_file = str(argv[2])

data = genfromtxt(input_file, delimiter=',')
x = data[1]
y = data[0]

#figure()
#plot(x, y)
#xlabel('x')
#ylabel('u(x,1)')
#savefig(output_file)
#
#print('Generated ' + output_file)

fig = plt.figure()
ax = plt.axes(projection='3d')

ax.plot3D(xline, yline, zline, 'gray')
