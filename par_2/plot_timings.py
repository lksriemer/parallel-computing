import numpy as np
import matplotlib.pyplot as plt

timings1 = np.array([9760,9426,10331,9441,9491,9388,9948,9472,9365,9387])
timings2 = np.array([94993,94322,94730,94426,94982,94403,97287,94805,94688,94654])
timings3 = np.array([943699,951331,946759,950015,945436,944827,947136,994202,952314,953473])

mean1 = np.mean(timings1)
mean2 = np.mean(timings2)
mean3 = np.mean(timings3)
std1 = np.std(timings1)
std2 = np.std(timings2)
std3 = np.std(timings3)


plt.yscale('log')
plt.xscale('log')


plt.fill_between([1000, 10000, 100000], np.array([mean1, mean2, mean3]) + np.array([std1, std2, std3]), np.array([mean1, mean2, mean3]) - np.array([std1, std2, std3]), alpha=0.5)


plt.show()