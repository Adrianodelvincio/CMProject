#include <iostream>
#include <Math/Vector4D.h>


double Cruijff(double *x, double *par){
	double arg = 0;
	double arg2 = (x[0] - par[3]);
	if (par[1] != 0 && par != 0){
		if( arg2 <= 0){
			arg = TMath::Exp(TMath::Power(arg2,2) / (2*par[1]*par[1] + par[4] *TMath::Power(arg2,2)));
		}
		else{
			arg = TMath::Exp(TMath::Power(arg2,2) / (2*par[2]*par[2] + par[5] *TMath::Power(arg2,2)));
		}
	}
	return arg;
}

int invmass(){
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

	rdf_up1.Display({"H1_PX", "H1_PY", "H1_PZ","H1_E","invMass"})->Print();
	double mass = 493.677; double max = 400e3;
	auto histMassdown = rdf_down1.Histo1D({"mass down","invariant Mass",128u,4500,6500},"invMass");
	auto histMassup = rdf_up1.Histo1D({"mass up","invariant Mass",128u,4500,6500},"invMass");

	auto histE1up = rdf_up1.Histo1D({"energy 1","Energy",128u,mass,max},"H1_E");
	auto histE2up = rdf_up1.Histo1D({"energy 2","Energy",128u,mass,max},"H2_E");
	auto histE3up = rdf_up1.Histo1D({"energy 3","Energy",128u,mass,max},"H3_E");
	auto histE1down = rdf_down1.Histo1D({"energy 1","Energy",128u,mass,max},"H1_E");
	auto histE2down = rdf_down1.Histo1D({"energy 2","Energy",128u,mass,max},"H2_E");
	auto histE3down = rdf_down1.Histo1D({"energy 3","Energy",128u,mass,max},"H3_E");

	auto BE_up = rdf_up1.Histo1D({"energy","Energy B+/-",128u,0,800e3},"B_E");
	auto BE_down = rdf_down1.Histo1D({"energy","Energy B+/-",128u,0,800e3},"B_E");

	auto Bmass_up = Bselection_up.Histo1D({"mass B","B mass",128u,4900,6900},"MassB");
	auto Bmass_down = Bselection_down.Histo1D({"mass B","B mass",128u,4900,6900},"MassB");
	//fit to the invariant mass of the B meson



	//Plot the invariant mass
	auto c = new TCanvas("c","c",900,900);
	auto p = new TPad("p","p",0,0,1,1);
	p->Divide(2,2,0.01,0.01);
	p->Draw();
	p->cd(1); histMassup->DrawClone("e");
	p->cd(2); histMassdown->DrawClone("e");
	p->cd(3);
	histE1up->DrawClone();
	histE2up->DrawClone("Same");
	histE3up->DrawClone("Same");
	p->cd(4);
	histE1down->DrawClone();
	histE2down->DrawClone("Same");
	histE3down->DrawClone("Same");
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
	BE_down->SetYTitle("Counts");
	BE_down->SetXTitle("Energy [MeV]");
	BE_down->DrawClone();
	bpad->cd(3);
	Bmass_up->SetYTitle("Counts");
	Bmass_up->SetXTitle("Mass B magnet UP [MeV]");
	Bmass_up->DrawClone();
	bpad->cd(4);
	Bmass_down->SetYTitle("Counts");
	Bmass_down->SetXTitle("Mass B magnet DOWN [MeV]");
	Bmass_down->DrawClone();

	//Save the files with the new data
	Bselection_up.Snapshot("DecayTree", "processed_Up.root",{"invMass","H1_ProbK","H2_ProbK","H3_ProbK",
	"H1_ProbPi","H2_ProbPi","H3_ProbPi", "H1_Charge", "H2_Charge","H3_Charge"});
	Bselection_down.Snapshot("DecayTree", "processed_Down.root",{"invMass","H1_ProbK","H2_ProbK","H3_ProbK",
	"H1_ProbPi","H2_ProbPi","H3_ProbPi", "H1_Charge", "H2_Charge","H3_Charge"});
	//print number of events after the filters
	auto rdf_up2 = rdf_up.Filter("!H1_isMuon && !H2_isMuon && !H3_isMuon").Count();
	auto rdf_down2 = rdf_down.Filter("!H1_isMuon && !H2_isMuon && !H3_isMuon").Count();
	std::cout << "Magnet UP events: " << rdf_up2.GetValue()  << std::endl;
	std::cout << "Magnet DOWN events: " << rdf_down2.GetValue() << std::endl;
	return 0;
}
