import numpy as np
import matplotlib.pyplot as plt


plt.loglog([250, 500, 1000, 2000, 4000], [1321, 3539, 16405, 91830, 628768], label="runtime")
plt.loglog([250, 500, 1000, 2000, 4000], [1321, 4*1321, 4*4*1321, 4*4*4*1321, 4*4*4*4*1321], label="quadratically")
plt.loglog([250, 500, 1000, 2000, 4000], [1321, 8*1321, 8*8*1321, 8*8*8*1321, 8*8*8*8*1321], label="cubically")

plt.legend(loc="upper left")

plt.savefig("runtime_invert.png")