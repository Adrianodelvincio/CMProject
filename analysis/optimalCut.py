import ROOT
import numpy as np
from matplotlib import pyplot as plt

ROOT.gSystem.Load("AutoDict_Functions_cxx.so") #load the library

#Mass of B particles
Bmass = 5279

x = np.linspace(1,300.,50)
y = np.linspace(1,300.,50)
x = Bmass - x
y = Bmass + y


cpp_code = """
// Function definition
double fexp(double *x, double *par ) { return par[0]*TMath::Exp(-par[1]*(x[0]-5000)); }
"""
# Inject the code in the ROOT interpreter
ROOT.gInterpreter.ProcessLine(cpp_code)

parameterFourBody = np.array([1.356e04,15.26,2.31e04,5024.,7.962e-07])
parameterCruijff = np.array([16.81,15.74, 5285.,0.09649,0.1015,1959.])
parameterFexp = np.array([32.24,0.0006275])

def Integral(x,y, n = 1000):
	xx = np.linspace(x,y,n)
	interval_width = (y - x)/n
	xx = interval_width/2 + xx
	fourbody = [ROOT.fourBodybackground( c,parameterFourBody ) for c in xx]
	signal = [ROOT.Cruijff(c,parameterCruijff ) for c in xx ]
	combinatorial = [ROOT.fexp(c, parameterFexp ) for c in xx]
	signal = np.sum(signal)*interval_width
	combinatorial = np.sum(combinatorial)*interval_width
	fourbody = np.sum(fourbody)*interval_width
	return signal, fourbody, combinatorial


Significance = np.zeros((len(x),len(y)))

def Value(signal,Combinatorial,fourbody):

	return signal/(np.sqrt(signal + Combinatorial + fourbody))


for i in range(0,len(x)):
	for j in range(0,len(y)):
		Events, background4body, Combinatorial = Integral(x[i],y[j])
		Significance[i][j] = Value(Events,Combinatorial,background4body)

indx , indy = np.unravel_index(np.argmax(Significance, axis = None), Significance.shape)

print("index : ", indx, " ", indy, " optimal x : ", x[indx], " optimal y ", y[indy])

plt.figure(1)
plt.title("Significance")
plt.plot(x, Significance[:][indy], color = "black", label = "fixed y = ymax")
plt.plot(y, Significance[indx][:], color = "green", label = r"fixed x = xmax")
plt.legend();
plt.figure(2)
plt.plot()
xx = np.linspace(5200,5400,100)
signal = [ROOT.Cruijff(c,parameterCruijff ) for c in xx ]
plt.plot(xx,signal)
plt.show()
