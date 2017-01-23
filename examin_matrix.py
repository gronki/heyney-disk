import numpy as np
import matplotlib.pyplot as plt


M = np.load('M.npy')
A = np.load('A.npy')

fig,ax = plt.subplots(figsize=(7,7),dpi=141)
ax.pcolor(M != 0, cmap='gray')
[ (plt.axvline(a,color='white'),plt.axhline(a,color='white')) for a in [0,16,32,48] ]
ax.set_aspect(1)
plt.show()

fig,ax = plt.subplots(figsize=(7,7),dpi=141)
ax.pcolor(np.arcsinh(M), cmap='coolwarm')
[ (plt.axvline(a,color='white'),plt.axhline(a,color='white')) for a in [0,16,32,48] ]
ax.set_aspect(1)
plt.show()
