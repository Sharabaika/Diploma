import matplotlib.pyplot as plt
import numpy as np

x = np.arange(0.01, 10, 0.1)
f = lambda t : 1.0/np.tanh(t)-1.0/t
y = [f(t) for t in x]

plt.xlabel("ξ")
plt.ylabel("M")
plt.plot(x,y)
plt.savefig("SavedPlots/Ланжевен.png")