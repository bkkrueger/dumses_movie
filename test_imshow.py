import copy
import matplotlib as mpl
from matplotlib.colors import colorConverter
import matplotlib.pyplot as pl
from matplotlib.ticker import MaxNLocator
from mpl_toolkits.axes_grid1 import make_axes_locatable
import numpy as np

Nx = 125*2+1
Ny = 55*2+1

x = (np.arange(Nx, dtype=float) - 0.5*(Nx-1)) / 10.0
y = (np.arange(Ny, dtype=float) - 0.5*(Ny-1)) / 10.0

xx = x.reshape((len(x),1))
yy = y.reshape((1,len(y)))

data = -1.0 \
       + 3.0 * np.exp(-((xx-10.3)**2 + (yy+1.2)**2)/7.7**2) \
       + 1.2 * np.exp(-((xx+9.1)**2 + (yy-3.7)**2)/4.8**2)

mask = np.zeros(data.shape)
mask[np.abs(data) < 0.15] = 1
masked_data = np.ma.masked_where(mask, data)

# transform to imshow orientation
plot_data = np.flipud(masked_data)
xx = y
yy = x

my_cmap = copy.copy(mpl.cm.seismic)
my_cmap.set_bad(color='k', alpha=1)

fig = pl.figure()
ax = fig.add_axes([0,0.05,1,0.9])

divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="20%", pad=0.05)

img = ax.imshow(plot_data, extent=[xx.min(), xx.max(), yy.min(), yy.max()],
      aspect=1.0, cmap=my_cmap)
cbar = pl.colorbar(img, cax=cax)
ax.contour(xx, yy, np.abs(x.reshape((len(x),1)).repeat(len(y),1)),
      levels=[4], colors="#00FF00")
ax.xaxis.set_major_locator(MaxNLocator(6))
ax.yaxis.set_major_locator(MaxNLocator(10))
cbar.locator = MaxNLocator(nbins=10)

