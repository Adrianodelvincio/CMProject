#include <iostream>
#include <Math/Vector4D.h>

int Dalitz(){
        //load the two datafile
        auto fileDown = "B2HHH_MagnetDown.root";
        auto fileUp = "B2HHH_MagnetUp.root";

        ROOT::RDataFrame rdf_up("DecayTree", fileUp);
        ROOT::RDataFrame rdf_down("DecayTree", fileDown);

	//for the momenta of a particle
	auto Momentum = [] (double px, double py, double pz){ return sqrt(px*px + py*py + pz*pz);};

	//for the invariant mass
	auto invMass = [] (double p1x, double p1y, double p1z, double p2x, double p2y, double p2z, double p3x, double p3y, double p3z){
	double KaonMass = 493.677; // MeV/c**2
	ROOT::Math::PxPyPzMVector k1(p1x,p1y,p1z, KaonMass);
	ROOT::Math::PxPyPzMVector k2(p2x,p2y,p2z, KaonMass);
	ROOT::Math::PxPyPzMVector k3(p3x,p3y,p3z, KaonMass);
	double invariant = k1.M2() + k2.M2() + k3.M2() + 2*k1.Dot(k2) + 2*k1.Dot(k3) + 2*k2.Dot(k3);
	return TMath::Sqrt(invariant);};

	//for the energy of the particle
	auto energyKaon = [] (double p1x, double p1y, double p1z){
	double KaonMass = 493.677; // MeV/c**2
	ROOT::Math::PxPyPzMVector k(p1x,p1y,p1z, KaonMass);
	return k.energy();};

	auto invB = [](double px,double py,double pz,double Energy){
	ROOT::Math::PxPyPzEVector k(px,py,pz,Energy);
	return k.M();};

	//Define Energy of Kaons, momenta and Energy of meson B+/-
	auto rdf_up1 = rdf_up.Define("invMass", invMass, {"H1_PX","H1_PY","H1_PZ","H2_PX","H2_PY","H2_PZ","H3_PX","H3_PY","H3_PZ"})
			.Define("H1_E",energyKaon, {"H1_PX","H1_PY","H1_PZ"})
			.Define("H2_E",energyKaon, {"H2_PX","H2_PY","H2_PZ"})
			.Define("H3_E",energyKaon, {"H3_PX","H3_PY","H3_PZ"})
			.Define("B_PX","H1_PX + H2_PX + H3_PX")
			.Define("B_PY","H1_PY + H2_PY + H3_PY")
			.Define("B_PZ","H1_PZ + H2_PZ + H3_PZ")
			.Define("B_E","H1_E + H2_E + H3_E")
			.Filter("H1_E <= 1e6 && H2_E <= 1e6 && H3_E <= 1e6")
			.Filter("!H1_isMuon && !H2_isMuon && !H3_isMuon")
			.Filter("B_E  <= 1e6"); // ???

	auto rdf_down1 = rdf_down.Define("invMass", invMass, {"H1_PX","H1_PY","H1_PZ","H2_PX","H2_PY","H2_PZ","H3_PX","H3_PY","H3_PZ"})
			.Define("H1_E",energyKaon, {"H1_PX","H1_PY","H1_PZ"})
			.Define("H2_E",energyKaon, {"H2_PX","H2_PY","H2_PZ"})
			.Define("H3_E",energyKaon, {"H3_PX","H3_PY","H3_PZ"})
			.Define("B_PX","H1_PX + H2_PX + H3_PX")
			.Define("B_PY","H1_PY + H2_PY + H3_PY")
			.Define("B_PZ","H1_PZ + H2_PZ + H3_PZ")
			.Define("B_E","H1_E + H2_E + H3_E")
			.Filter("!H1_isMuon && !H2_isMuon && !H3_isMuon")
			.Filter("H1_E <= 1e6 && H2_E <= 1e6 && H2_E <= 1e6")
			.Filter("B_E  <= 1e6"); // ???

	auto Bselection_up = rdf_up1.Define("MassB",invB, {"B_PX","B_PY","B_PZ","B_E"});
	auto Bselection_down = rdf_down1.Define("MassB",invB, {"B_PX","B_PY","B_PZ","B_E"});

	//Now produce the Dalitz plot to study the decay of the B+/- particles

	return 0;
}
