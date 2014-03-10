import numpy as np
import matplotlib
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt

from test_main import *

def plot3d(x,y,z, xl='X', yl='Y', zl='Z', c='#cc0000'):
    ax.scatter(x, y, z, c=c)
    ax.set_xlabel(xl)
    ax.set_ylabel(yl)
    ax.set_zlabel(zl)

l = LedaFits('vis_00.uvfits')
xyz = l.d_array_geometry["STABXYZ"]
x,y,z = np.split(xyz, 3, axis=1)

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
plot3d(x, y, z, 'X', 'Y', 'Z', c='#00cc00')
plt.show()

bls = coords.computeBaselineVectors(xyz)

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
uvw = coords.computeUVW(bls, H=0, d=np.rad2deg(0)) * 1e6
u,v,w = np.split(uvw, 3, axis=1)

plot3d(u, v, w, 'U', 'V', 'W')

uvw = coords.computeUVW(bls, H=0, d=np.deg2rad(34.07)) * 1e6
u,v,w = np.split(uvw, 3, axis=1)
plot3d(u, v, w, 'U', 'V', 'W', c='#00cc00')

plt.show()