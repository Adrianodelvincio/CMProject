Welcome!
This project is an analysis to study the CP asymmetry in the decay of B+/- meson in 3 kaons.
The data were acquired at LHCb detector, one of the 4 largest experiments of Cern, in 2011.
The analysis is developed using ROOT data Frame, so ROOT is required to run the code.

Briefly an explanation of the project:

 - In display.cc some Histogram are defined, for a first visualization of the Data. The data files contain 25 variables for the analysis divided in column, only 20 are analyzed. The Data can be loaded [here][https://opendata.cern.ch/record/4900]
 - In invmass.cc the first analysis is performed: the invariant mass of the B+/- candidates is computed with the three kaons momenta and then some histograms are filled. Then two new datafiles are created (processed_Up.root and processed_Down.root), with only the relevant quantities for the analysis.
 - In globalAsymmetry.cpp some cuts are defined and applied to discard Muons and Pions, the recostructed invariant mass is fitted, using the Cruijff function (a gaussian function with asymmetric tails) for the signal and an exponential + 4body function for the background. Then the global Asymmetry is computed.  
 - In search.py an optimisation is performed to define the best cut for the B+/- candidates. Only the particle that have a invariant mass close to the known value (5279 MeV/c^2) will be analyzed later, to study the CP asymmetry
 - In Dalitz.cpp Dalitz plots are defined, with the same cuts used before in globalAsymmetry. The invariant masses for each pair of kaons is computed, to study loal asymmetries. 
 - in 2dPlot.py we study the local CP asymmetry and pcolor map is generated. The axis are the M_low and M_high invariant masses and for each bin the asymmetry between B+/B- deacy is computed. 


	[x] finire alla massa dei candidati B+/- (Crystal Ball function)
	[x] capire se le masse invarianti calcolate sono giuste o sbagliate
	[x] capire come funziona la color map di matplotlib
	[ ] organizzare la search per ottimizzare i tagli (probK and probPi)

