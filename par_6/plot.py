import numpy as np
import matplotlib.pyplot as plt


plt.loglog([32, 128, 512, 2048, 8192, 32768], [0.001, 0.02, 0.25, 4.3, 64, 343])
plt.loglog([32, 128, 512, 2048, 8192, 32768], [0.03, 0.17, 2.8, 90, 355, 467])

plt.show()