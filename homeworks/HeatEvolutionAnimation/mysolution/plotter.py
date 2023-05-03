import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

M = 100 # number of frames
dpi = 100 # dots per inch

data = np.loadtxt("solution.txt")
axis_length = int(np.sqrt(data.shape[1]))
data = data.reshape((data.shape[0], axis_length, axis_length))

step = int(data.shape[0] / M)

fig, ax = plt.subplots()
im = ax.imshow(data[0], vmin=0.)
plt.colorbar(im)
plt.axis('off')
fig.tight_layout()

def update(i):
    im.set_data(data[i*step])
    return im,

ani = FuncAnimation(fig, update, frames=M, interval=200, blit=True, repeat=False)

ani.save("animation.gif", dpi=dpi, writer='pillow')
plt.show()
