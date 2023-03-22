Welcome!
This project is an analysis to study the CP asymmetry in the decay of B+/- meson in 3 kaons.
The data were acquired at LHCb detector, one of the 4 largest experiments of Cern, in 2011.
The analysis is developed using ROOT data Frame, so ROOT is required to run the code.

Briefly an explanation of the project:

 - In display.cc some Histogram are defined, for a first visualization of the Data. The data files contain 25 variables for the analysis divided in column, only 20 are analyzed. The Data can be loaded [here][https://opendata.cern.ch/record/4900]

 - In invmass.cc the first analysis is performed: the invariant mass of the B+/- candidates is computed with the three kaons momenta and some histograms are filled. Then two new datafiles are created (processed_Up.root and processed_Down.root).

 - In globalAsymmetry.cpp some cuts are defined and applied to discard events with Muons and Pions. The recostructed invariant mass is fitted, using the Cruijff function (a gaussian function with asymmetric tails) for the signal and an exponential + 4body function for the background. In the end a global Asymmetry between B+/B- decays is computed.  

 - In optimalCut.py an optimisation is performed to define the best cut for the B+/- candidates. Only the particle that have a invariant mass close to the known value (5279 MeV/c^2) will be analyzed later. The accepted interval (5279 - x, 5279 + y) is selected looking at the maximum of the function (s are the expeted signal events and b is the expected background):
  $$\frac{s(x,y)}{\sqrt{s(x,y) + b(x,y)}}$$ 

 - In Dalitz.cpp Dalitz plots are defined, with the same cuts used before in globalAsymmetry. The invariant masses for each pair of kaons is computed, to study local asymmetries, then a new dataset, with the only events selected for final analysis, is produced. 
 
 - in 2dPlot.py we study the local CP asymmetry and a color map is generated. The axis are the M_low and M_high invariant masses and for each bin the asymmetry between B+/B- deacy is computed. Also the statistical significance.


To Run the code, you can use the following command in the directory of the project : root analysis/scriptname -l
The analysis follow this order:
 - display.cpp (To visualize the data)
 - invmass.cpp
 - globalAsymmetry.cpp
 - simuldata.cpp
 - dalitz.cpp
 - optimalCut.py
 - 2dPlot.py

debugpcolor.py is a scripts used to understand how the pcolormash function in matplotlib library produce the colormap, cruijff.cpp to visualize the shape of expected signal distribution. They are not part of the main analysis
