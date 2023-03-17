import ROOT
import numpy as np
from matplotlib import pyplot as plt

fileUp = "Evt_Selected/ordered_up.root"
fileDown = "Evt_Selected/ordered_down.root"
rdf_up = ROOT.RDataFrame ("DecayTree", fileUp)
rdf_dw = ROOT.RDataFrame ("DecayTree", fileDown)

#Save values inside arrays
lowup = np.asarray(rdf_up.Filter("Bcharge == +1").Take['double']("roLow"))
lowdw = np.asarray(rdf_dw.Filter("Bcharge == +1").Take['double']("roLow"))
nlowup = np.asarray(rdf_up.Filter("Bcharge == -1").Take['double']("roLow"))
nlowdw = np.asarray(rdf_dw.Filter("Bcharge == -1").Take['double']("roLow"))

highup = np.asarray(rdf_up.Filter("Bcharge == +1").Take['double']("roHigh"))
highdw = np.asarray(rdf_dw.Filter("Bcharge == +1").Take['double']("roHigh"))
nhighup = np.asarray(rdf_up.Filter("Bcharge == -1").Take['double']("roHigh"))
nhighdw = np.asarray(rdf_dw.Filter("Bcharge == -1").Take['double']("roHigh"))

#print("lowup : ", lowup.shape, " lowdw " , lowdw.shape, " nlowup ", nlowup.shape, " nlowdw ", nlowdw.shape)

#Now we create a color 2d histogram for the asymmetry, to study the local asymmetries
rolow = np.concatenate((lowup,lowdw),axis=None)
nrolow = np.concatenate((nlowup,nlowdw),axis=None)
rohigh = np.concatenate((highup,highdw),axis=None)
nrohigh = np.concatenate((nhighup,nhighdw),axis=None)
#print(rolow.shape,rohigh.shape,nrolow.shape,nrohigh.shape)


print("total number of event B+ : " , len(rolow), " total number of event B- : " , len(nrolow))
#now we create a 2d histogram
xbinning = 10
ybinning = 14

#define the minimum statistics required for each bin from the formula of the error of the asymmetry
sigma_expected = 10/100 # max error percentage required for each bin
Nmin = int(1/(2*sigma_expected**2))
print("Minimum statistics for each bin : %d" % Nmin )
Base = 6
Cmin = Nmin
rolow = rolow/1e6 ; nrolow = nrolow/1e6 ; rohigh = rohigh/1e6 ; nrohigh = nrohigh/1e6

xmin = min(rolow)
xmax = max(rolow)
ymin = min(rohigh)
ymax = max(rohigh)

xmin1 = min(nrolow)
xmax1 = max(nrolow)
ymin1 = min(nrohigh)
ymax1 = max(nrohigh)

xmin = np.minimum(xmin,xmin1)
xmax = np.maximum(xmax,xmax1)
ymin = np.minimum(ymin,ymin1)
ymax = np.maximum(ymax,ymax1)

edge_x = np.logspace(np.emath.logn(Base,xmin),np.emath.logn(Base,xmax), num = xbinning, base = Base, endpoint = True)
edge_y = np.linspace(ymin,ymax, num = ybinning, endpoint = True)

plt.figure(1)
plt.title(r"Matter $\rho_{low}$ vs $\rho_{high}$ ")
matter, xedges, yedges, Image = plt.hist2d(rolow,rohigh, bins = [edge_x,edge_y], alpha = 0.5, cmin = Cmin)
plt.xlabel(r'$\rho_{low}$')
plt.ylabel(r'$\rho_{high}$')


plt.figure(2)
plt.title(r"antiMatter $\rho_{low}$ vs $\rho_{high}$ ")
antimatter, xedges, yedges, Image = plt.hist2d(nrolow,nrohigh, bins = [edge_x,edge_y], alpha = 0.5, cmin = Cmin)
plt.xlabel(r'$\rho_{low}$')
plt.ylabel(r'$\rho_{high}$')
plt.colorbar()

plt.figure(3)
plt.errorbar(rolow,rohigh,linestyle = '', marker = '.', markersize = 0.5, color = "black")
plt.errorbar(nrolow,nrohigh, linestyle = '', marker = '.', markersize = 0.5, color = "yellow")

plt.figure(4)
plt.title(r"$\rho_{low}$")
plt.hist(rolow, bins = 20, histtype = 'step')
#now we can compute the asymmetry in the following way:
asymmetry = np.zeros(matter.shape)

#Mask for the bins with zero asymmetry
sum = matter + antimatter
difference = antimatter - matter
mask = (sum == 0)
mask2 = (sum != 0)
#print("sum equal to 0 : " , sum[mask])
#print("mask : ", mask)

#mask 2d arrays
asymmetry[mask] = 0
np.putmask(asymmetry, mask2 , difference / sum)

#print("asymmetry : ", asymmetry)
print(asymmetry.shape, xedges.shape, yedges.shape)
asymmetry = np.transpose(asymmetry)

plt.figure(5)
plt.errorbar(rolow,rohigh,linestyle = '', marker = '.', alpha = 0.25,markersize = 0.5, color = 'black')
plt.errorbar(nrolow,nrohigh, linestyle = '', marker = '.', alpha = 0.25, markersize = 0.5, color = 'yellow')
plt.xlabel(r"$\rho_{low}$ [$GeV^{2}]$")
plt.ylabel(r"$\rho_{high}$ [$GeV^{2}$]")
plt.pcolormesh(xedges,yedges,asymmetry,cmap = 'seismic', edgecolors = 'k', linewidths = '0.5')
plt.colorbar()
plt.show()
