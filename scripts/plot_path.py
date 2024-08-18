import sys

import numpy as np
from matplotlib import pyplot as plt

pts = []
with open(sys.argv[1], 'r') as f:
    for line in f.readlines():
        pts.append([float(i) for i in line.split(',')])

pts = np.array(pts)
fig = plt.figure()
ax = fig.add_subplot(projection='3d')
ax.plot(pts[:, 0], pts[:, 1], pts[:, 2])
plt.savefig('test.png', dpi=400)