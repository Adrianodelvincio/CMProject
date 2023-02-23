#include <iostream>
#include <Math/Vector4D.h>

int simuldata(){
        //load the two datafile
        auto fileDown = "PhaseSpaceSimulation.root";
        //auto fileUp = "B2HHH_MagnetUp.root";

        ROOT::RDataFrame rdf_up("PhaseSpaceTree", fileDown);
       // ROOT::RDataFrame rdf_down("DecayTree", fileDown);

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
			.Filter("B_E  <= 1e6");


	auto Bselection_up = rdf_up1.Define("MassB",invB, {"B_PX","B_PY","B_PZ","B_E"});


	rdf_up1.Display({"H1_PX", "H1_PY", "H1_PZ","H1_E","invMass"})->Print();
	double mass = 493.677; double max = 400e3;

	auto histMassup = rdf_up1.Histo1D({"mass up","invariant Mass",128u,5270,5300},"invMass");
	auto histE1up = rdf_up1.Histo1D({"energy 1","Energy",128u,mass,max},"H1_E");
	auto histE2up = rdf_up1.Histo1D({"energy 2","Energy",128u,mass,max},"H2_E");
	auto histE3up = rdf_up1.Histo1D({"energy 3","Energy",128u,mass,max},"H3_E");
	auto BE_up = rdf_up1.Histo1D({"energy","Energy B+/-",128u,-100,1000e3},"B_E");
	auto Bmass_up = Bselection_up.Histo1D({"mass B","B mass",128u,4900,6900},"MassB");

	//Define the Dalitz Plot for the simulated Data, to study the decay of B+/- particles

	auto mxy = [](double bx, double by, double bz, double bE, double kx, double ky, double kz, double kE){
	ROOT::Math::PxPyPzEVector k(kx,ky,kz,kE);
	ROOT::Math::PxPyPzEVector b(bx,by,bz,bE);
	return (b - k).M();
};

	auto dalitz = rdf_up1.Define("Bcharge","(H1_Charge + H2_Charge + H3_Charge)")
				.Filter("Bcharge == +1")
				.Define("m12", mxy ,{"B_PX","B_PY", "B_PZ", "B_E", "H3_PX","H3_PY","H3_PZ", "H3_E"})
				.Define("m23", mxy ,{"B_PX","B_PY", "B_PZ", "B_E", "H1_PX","H1_PY","H1_PZ", "H1_E"})
				.Define("m13", mxy, {"B_PX","B_PY", "B_PZ", "B_E", "H2_PX","H2_PY","H2_PZ", "H2_E"});

	auto hist_ro3 = dalitz.Filter("H3_Charge == -1").Histo1D({"mm","K1 K2 invariant mass",128u,880., 9000.},"m12");
	auto hist_ro1 = dalitz.Filter("H1_Charge == -1").Histo1D({"", "", 128u, 880., 9000.},"m23");
	auto hist_ro2 = dalitz.Filter("H2_Charge == -1").Histo1D({"", "", 128u, 880., 9000.},"m13");

	int total = 3000;
	auto graph_ro3 = dalitz.Filter("H3_Charge == -1").Range(total).Graph("m12","m23");
	auto graph_ro1 = dalitz.Filter("H1_Charge == -1").Range(total).Graph("m23","m12");
	auto graph_ro2 = dalitz.Filter("H2_Charge == -1").Range(total).Graph("m13","m23");

	hist_ro3->Add(hist_ro1.GetPtr(),1); hist_ro3->Add(hist_ro2.GetPtr(),1);

	dalitz.Display({"m12", "m23", "m13","invMass"},10,10)->Print();

	//Plot the invariant mass
	auto c = new TCanvas("c","c",900,900);
	auto p = new TPad("p","p",0,0,1,1);
	p->Divide(2,2,0.01,0.01);
	p->Draw();
	p->cd(1);
	histMassup->DrawClone();
	p->cd(2);
	histE1up->DrawClone();
	histE2up->DrawClone("Same");
	histE3up->DrawClone("Same");
	p->cd(3);
	hist_ro3->SetTitle("K1 K2 invariant mass");
	hist_ro3->SetXTitle("invariant mass k+ k+ [MeV**2]");
	hist_ro3->DrawClone();
	p->cd(4);
	graph_ro3->SetTitle("Dalitz plot k+ k+ k- ; invariant mass m12 ++ ; invariant mass m23 +-");
	float marksize = 0.2;
	graph_ro3->SetMarkerSize(marksize);
	graph_ro2->SetMarkerSize(marksize);
	graph_ro1->SetMarkerSize(marksize);

	graph_ro3->SetMarkerStyle(20);
	graph_ro2->SetMarkerStyle(20);
	graph_ro1->SetMarkerStyle(20);
	graph_ro3->SetMarkerColor(2);
	graph_ro2->SetMarkerColor(2);
	graph_ro1->SetMarkerColor(2);
	graph_ro3->DrawClone("APL");
	graph_ro2->DrawClone("same AP");
	graph_ro1->DrawClone("same AP");
	// B related quantities
	auto bcanvas = new TCanvas("bcavas","bcanvas",900,900);
	auto bpad = new TPad("bpad","bpad",0,0,1,1);
	bpad->Divide(2,2,0.01,0.01);
	bpad->Draw();
	bpad->cd(1);
	BE_up->SetYTitle("Counts");
	BE_up->SetXTitle("Energy [MeV]");
	BE_up->DrawClone();
	bpad->cd(2);
	Bmass_up->SetYTitle("Counts");
	Bmass_up->SetXTitle("Mass [MeV]");
	Bmass_up->DrawClone();
	bpad->cd(4);
	return 0;
}
