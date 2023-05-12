import matplotlib.pyplot as plt
import numpy as np

stable = np.loadtxt("upwind_quadrature_solution_stable.txt")
unstable = np.loadtxt("upwind_quadrature_solution_unstable.txt")

shape = (int(np.sqrt(stable.size)),)*2

fig, (ax1, ax2) = plt.subplots(1, 2)
im1 = ax1.imshow(np.reshape(stable, shape), origin="lower")
ax1.set_title("Upwind Quadrature")
fig.colorbar(im1)

im2 = ax2.imshow(np.reshape(unstable, shape), origin="lower")
ax2.set_title("Standard Quadrature")
fig.colorbar(im2)

plt.savefig("solution.png")
plt.show()
