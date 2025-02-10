import numpy as np
import matplotlib.pyplot as plt

timings1 = np.array([9824,10897,9552,9710,9585,9575,9632,9713,9777,9560])
timings2 = np.array([100529,100740,102453,101109,101115,100976,101176,100968,101120,101097])
timings3 = np.array([1091842,1090154,1094093,1093855,1097918,1103877,1116322,1093820,1086232,1087975])

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