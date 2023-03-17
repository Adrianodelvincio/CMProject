import numpy as np
from matplotlib import pyplot as plt

x = np.random.uniform(0,1,10000)
y = np.random.uniform(0,1,10000)
x = 1 - x*x
y = 2*(1 - y*y)

edge_x = np.linspace(0,1,20)
edge_y = np.linspace(0,2,40)

plt.figure(1)
nn , xedges, yedges , image = plt.hist2d(x,y, bins = [edge_x,edge_y], alpha = 0.5)

plt.figure(2)
plt.hist(x, bins = 20)
plt.hist(y,bins = 20)

plt.figure(3)
nn = nn /360
plt.title("Non-transposed x,y matrix")
plt.pcolor(yedges,xedges,nn,cmap = 'seismic')
plt.colorbar()

plt.figure(4)
plt.errorbar(x,y,marker = '.', linestyle = '', markersize = 1)

plt.figure(5)
plt.title("transposed x,y matrix")
nn = np.transpose(nn)
plt.pcolormesh(xedges,yedges,nn,cmap = 'seismic')

plt.show()
