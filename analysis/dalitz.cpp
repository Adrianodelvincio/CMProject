#include <iostream>
#include <string>
#include <Math/Vector4D.h>

int dalitz(){
   	ROOT::EnableImplicitMT(2);
        //load the two datafile

        auto fileDown = "B2HHH_MagnetDown.root";
        auto fileUp = "B2HHH_MagnetUp.root";
	ROOT::RDataFrame rdf_up("DecayTree", fileUp);
        ROOT::RDataFrame rdf_down("DecayTree", fileDown);

	//auto fileUp = "PhaseSpaceSimulation.root";
        //ROOT::RDataFrame rdf_up("PhaseSpaceTree", fileUp);

	//The invariant mass
	auto invMass = [] (double p1x, double p1y, double p1z, double p2x, double p2y, double p2z, double p3x, double p3y, double p3z){
		double KaonMass = 493.677; // MeV/c**2
		ROOT::Math::PxPyPzMVector k1(p1x,p1y,p1z, KaonMass);
		ROOT::Math::PxPyPzMVector k2(p2x,p2y,p2z, KaonMass);
		ROOT::Math::PxPyPzMVector k3(p3x,p3y,p3z, KaonMass);
		double invariant = k1.M2() + k2.M2() + k3.M2() + 2*k1.Dot(k2) + 2*k1.Dot(k3) + 2*k2.Dot(k3);
	return TMath::Sqrt(invariant);};

	//Energy of the particle
	auto energyKaon = [] (double p1x, double p1y, double p1z){
		double KaonMass = 493.677; // MeV/c**2
		ROOT::Math::PxPyPzMVector k(p1x,p1y,p1z, KaonMass);
	return k.energy();};

	//Invariant mass of two particles
	auto mxy = [](double bx, double by, double bz, double bE, double kx, double ky, double kz, double kE){
		ROOT::Math::PxPyPzEVector k(kx,ky,kz,kE);
		ROOT::Math::PxPyPzEVector b(bx,by,bz,bE);
		auto c = b - k;
	return (b-k).M2();};

	//Define Energy of Kaons, momenta and Energy of meson B+/-
	auto rdfKinematicUp = rdf_up.Define("invMass", invMass, {"H1_PX","H1_PY","H1_PZ","H2_PX","H2_PY","H2_PZ","H3_PX","H3_PY","H3_PZ"})
			.Define("H1_E",energyKaon, {"H1_PX","H1_PY","H1_PZ"})
			.Define("H2_E",energyKaon, {"H2_PX","H2_PY","H2_PZ"})
			.Define("H3_E",energyKaon, {"H3_PX","H3_PY","H3_PZ"})
			.Define("B_PX","H1_PX + H2_PX + H3_PX")
			.Define("B_PY","H1_PY + H2_PY + H3_PY")
			.Define("B_PZ","H1_PZ + H2_PZ + H3_PZ")
			.Define("B_E","H1_E + H2_E + H3_E")
			.Filter("H1_E <= 1e6 && H2_E <= 1e6 && H3_E <= 1e6")
			.Filter("!H1_isMuon && !H2_isMuon && !H3_isMuon")
			.Filter("B_E  <= 1e6");
	
	auto rdfKinematicDown = rdf_down.Define("invMass", invMass, {"H1_PX","H1_PY","H1_PZ","H2_PX","H2_PY","H2_PZ","H3_PX","H3_PY","H3_PZ"})
			.Define("H1_E",energyKaon, {"H1_PX","H1_PY","H1_PZ"})
			.Define("H2_E",energyKaon, {"H2_PX","H2_PY","H2_PZ"})
			.Define("H3_E",energyKaon, {"H3_PX","H3_PY","H3_PZ"})
			.Define("B_PX","H1_PX + H2_PX + H3_PX")
			.Define("B_PY","H1_PY + H2_PY + H3_PY")
			.Define("B_PZ","H1_PZ + H2_PZ + H3_PZ")
			.Define("B_E","H1_E + H2_E + H3_E")
			.Filter("!H1_isMuon && !H2_isMuon && !H3_isMuon")
			.Filter("H1_E <= 1e6 && H2_E <= 1e6 && H2_E <= 1e6")
			.Filter("B_E  <= 1e6");
	
	//Cuts and event selection
	auto Pselect_up = rdfKinematicUp.Filter("H1_ProbPi <= 0.6")
				.Filter("H2_ProbPi <= 0.6")
				.Filter("H3_ProbPi <= 0.6")
				.Filter("H1_ProbK >= 0.5")
				.Filter("H2_ProbK >= 0.5")
				.Filter("H3_ProbK >= 0.5");

	auto Pselect_down =  rdfKinematicDown.Filter("H1_ProbPi <= 0.6")
				.Filter("H2_ProbPi <= 0.6")
				.Filter("H3_ProbPi <= 0.6")
				.Filter("H1_ProbK >= 0.5")
				.Filter("H2_ProbK >= 0.5")
				.Filter("H3_ProbK >= 0.5");

	auto mass_up_selection = Pselect_up.Filter("invMass <= 5347")
				.Filter("invMass >= 5223")
				.Define("Bcharge","H1_Charge + H2_Charge + H3_Charge");
	auto mass_down_selection = Pselect_down.Filter("invMass <= 5347")
				.Filter("invMass >= 5223")
				.Define("Bcharge","H1_Charge + H2_Charge + H3_Charge");

	//Define invariant masses
	//m12 -> m13 ; m23 -> m12 ; m13 -> m23
	auto dalitz_up = mass_up_selection
				.Define("m13", mxy ,{"B_PX","B_PY", "B_PZ", "B_E", "H3_PX","H3_PY","H3_PZ", "H3_E"})
				.Define("m12", mxy ,{"B_PX","B_PY", "B_PZ", "B_E", "H1_PX","H1_PY","H1_PZ", "H1_E"})
				.Define("m23", mxy, {"B_PX","B_PY", "B_PZ", "B_E", "H2_PX","H2_PY","H2_PZ", "H2_E"})
				.Define("m12r","TMath::Sqrt(m12)")
				.Define("m13r","TMath::Sqrt(m13)")
				.Define("m23r","TMath::Sqrt(m23)");

	auto dalitz_down = mass_down_selection
                                .Define("m13", mxy ,{"B_PX","B_PY", "B_PZ", "B_E", "H3_PX","H3_PY","H3_PZ", "H3_E"})
                                .Define("m12", mxy ,{"B_PX","B_PY", "B_PZ", "B_E", "H1_PX","H1_PY","H1_PZ", "H1_E"})
                                .Define("m23", mxy, {"B_PX","B_PY", "B_PZ", "B_E", "H2_PX","H2_PY","H2_PZ", "H2_E"})
				.Define("m12r","TMath::Sqrt(m12)")
				.Define("m13r","TMath::Sqrt(m13)")
				.Define("m23r","TMath::Sqrt(m23)");

	//invariant mass Histograms B+ -> K+K+K-
	auto hist_m12_up   = dalitz_up.Filter("Bcharge == 1").Histo1D({"mm","K^{+}K^{+}", 128u, 50e3, 30e6},"m12");
	auto hist_m12_down = dalitz_down.Filter("Bcharge == 1").Histo1D({"mm","K^{+}K^{+}", 128u, 50e3, 30e6},"m12");
	auto hist_m13_up   = dalitz_up.Filter("Bcharge == 1").Histo1D({"","",128u,0,23e6},"m13");
	auto hist_m13_down = dalitz_down.Filter("Bcharge == 1").Histo1D({"","",128u,0,23e6},"m13");
	auto hist_m23_up   = dalitz_up.Filter("Bcharge == 1").Histo1D({"","",128u,0,30e6},"m23");
	auto hist_m23_down = dalitz_down.Filter("Bcharge == 1").Histo1D({"","",128u,0,30e6},"m23");

	//Dalitz plots
	auto M12M13_up = dalitz_up.Filter("Bcharge == 1").Filter("H1_Charge == -1").Graph("m12","m13");
	auto M12M13_down = dalitz_down.Filter("Bcharge == 1").Filter("H1_Charge == -1").Graph("m12","m13");
	auto M12M23_up = dalitz_up.Filter("Bcharge == 1").Graph("m12","m23");
	auto M12M23_down = dalitz_down.Filter("Bcharge == 1").Graph("m12","m23");
	auto M12M13_up_neg = dalitz_up.Filter("Bcharge == -1").Filter("H1_Charge == +1").Graph("m12","m13");
	auto M12M13_down_neg = dalitz_down.Filter("Bcharge == -1").Filter("H1_Charge == +1").Graph("m12","m13");
	auto M12M23_up_neg = dalitz_up.Filter("Bcharge == -1").Graph("m12","m23");
	auto M12M23_down_neg = dalitz_down.Filter("Bcharge == -1").Graph("m12","m23");
	auto M13M23_up = dalitz_up.Filter("Bcharge == 1").Filter("H1_Charge == -1").Graph("m13","m23");
        auto M13M23_down = dalitz_down.Filter("Bcharge == 1").Filter("H1_Charge == -1").Graph("m13","m23");
	auto M13M23_up_neg = dalitz_up.Filter("Bcharge == -1").Graph("m13","m23");
        auto M13M23_down_neg = dalitz_down.Filter("Bcharge == -1").Graph("m13","m23");

	//B- -> K+K-K-
	auto hist_m12_up_neg = dalitz_up.Filter("Bcharge == -1").Histo1D({"hist_m12_up_neg","K^{-}K^{-}", 128u, 50e3, 30e6},"m12");
	auto hist_m12_down_neg = dalitz_down.Filter("Bcharge == -1").Histo1D({"hist_m12_down_neg","K^{-}K^{-}", 128u, 50e3, 30e6},"m12");
	auto hist_m13_up_neg = dalitz_up.Filter("Bcharge == -1").Histo1D({"hist_m13_up_neg","",128u,0,23e6},"m13");
	auto hist_m13_down_neg = dalitz_down.Filter("Bcharge == -1").Histo1D({"hist_m13_down_neg","",128u,0,23e6},"m13");
	auto hist_m23_up_neg = dalitz_up.Filter("Bcharge == -1").Histo1D({"hist_m23_up_neg","",128u,0,30e6},"m23");
	auto hist_m23_down_neg = dalitz_down.Filter("Bcharge == -1").Histo1D({"hist_m23_down_neg","",128u,0,30e6},"m23");

	//Plots and histograms
	gStyle->SetOptStat(0);
	auto Dalitz = new TCanvas("Dalitz","Dalitz M12M23 M12M23",900,900);
	auto pDalitz = new TPad("pDalitz","p",0,0,1,1);
	pDalitz->Divide(2,1,0.02,0.02);
	pDalitz->Draw();
	float marksize = 1;
	pDalitz->cd(1);
	M12M23_up->SetTitle("m^{2}_{12} vs m^{2}_{13}; invariant mass m12 ; invariant mass m23 ");
	M12M23_up->SetMarkerColorAlpha(kBlack, 0.35);
	M12M23_up->DrawClone("AP");
	M12M23_up_neg->SetMarkerColorAlpha(kBlack, 0.35);
	M12M23_up_neg->DrawClone("P");
	M12M23_down->SetMarkerColorAlpha(kBlack, 0.35);
	M12M23_down->DrawClone("P");
	M12M23_down_neg->SetMarkerColorAlpha(kBlack, 0.35);
	M12M23_down_neg->DrawClone("P");
	pDalitz->cd(2);
	M12M13_up->SetTitle("m^{2}_{12} vs m^{2}_{13} ; invariant mass m12 ; invariant mass m13");
	M12M13_up->SetMarkerColorAlpha(kBlack, 0.35);
	M12M13_up->DrawClone("AP");
	M12M13_up_neg->SetMarkerColorAlpha(kBlack, 0.35);
	M12M13_up_neg->DrawClone("P");
        M12M13_down->SetMarkerColorAlpha(kBlack, 0.35);
        M12M13_down->DrawClone("P");
        M12M13_down_neg->SetMarkerColorAlpha(kBlack, 0.35);
        M12M13_down_neg->DrawClone("P");

	auto c = new TCanvas("c","Dalitz M13M23",900,900);
        auto p = new TPad("p","p",0,0,1,1);
        p->Divide(1,1,0.02,0.02);
        p->Draw();
	p->cd(1);
	M13M23_up->SetTitle("m^{2}_{13} VS m^{2}_{23}; invariant mass m13 ; invariant mass m23 ");
        M13M23_up->SetMarkerColorAlpha(kBlack, 0.35);
        M13M23_up->DrawClone("AP");
        M13M23_up_neg->SetMarkerColorAlpha(kBlack, 0.35);
        M13M23_up_neg->DrawClone("P");
        M13M23_down->SetMarkerColorAlpha(kBlack, 0.35);
        M13M23_down->DrawClone("P");
        M13M23_down_neg->SetMarkerColorAlpha(kBlack, 0.35);
        M13M23_down_neg->DrawClone("P");

	auto d = new TCanvas("d","B+ -> k+ k+ k-",900,900);
	auto q = new TPad("q","q",0,0,1,1);

	q->Divide(2,2,0.02,0.02);
	q->Draw();

	q->cd(1);
	hist_m12_down->SetLineColor(2);
	hist_m12_down->DrawClone();
	hist_m12_down->GetXaxis()->SetTitle("m12 [MeV^{2}] ");
	hist_m12_up->SetLineColor(1);
	hist_m12_up->DrawClone("same");

	q->cd(2);

	hist_m13_down->SetLineColor(2);
	hist_m13_down->DrawClone();
	hist_m13_down->GetXaxis()->SetTitle("m13 [MeV^{2}] ");
	hist_m13_up->SetLineColor(1);
	hist_m13_up->DrawClone("SAME");

	q->cd(3);
	hist_m23_down->GetXaxis()->SetTitle("m23 [MeV^{2}] ");
	hist_m23_down->SetLineColor(2);
	hist_m23_down->DrawClone();
	hist_m23_up->SetLineColor(1);
	hist_m23_up->DrawClone("same");


	auto e = new TCanvas("e","B- -> k- k- k+",900,900);
	auto f = new TPad("f","f",0,0,1,1);

	f->Divide(2,2,0.02,0.02);
	f->Draw();

	f->cd(1);
	hist_m12_down_neg->SetLineColor(2);
	hist_m12_down_neg->DrawClone();
	hist_m12_up_neg->GetXaxis()->SetTitle(" m_{12} [MeV^{2}] ");
	hist_m12_up_neg->SetLineColor(1);
	hist_m12_up_neg->DrawClone("same");


	auto legend3 = new TLegend(0.1,0.7,0.48,0.9);
	legend3->SetHeader("m_{12} histogram");
	legend3->AddEntry("hist_m12_up_neg","m_{12} MAGNET UP", "l");
	legend3->AddEntry("hist_m12_down_neg","m_{12} MAGNET DOWN", "l");
	legend3->SetTextSize(0.03);
	legend3->Draw();

	f->cd(2);
	hist_m13_up_neg->GetXaxis()->SetTitle(" m_{13} [MeV^{2}] ");
	hist_m13_up_neg->SetLineColor(1);
	hist_m13_up_neg->DrawClone();
	hist_m13_down_neg->SetLineColor(2);
	hist_m13_down_neg->DrawClone("same");

	auto legend4 = new TLegend(0.1,0.7,0.48,0.9);
	legend4->SetHeader("m_{13} histogram");
	legend4->AddEntry("hist_m13_up_neg","m_{13} MAGNET UP", "l");
	legend4->AddEntry("hist_m13_down_neg","m_{13} MAGNET DOWN", "l");
	legend4->SetTextSize(0.03);
	legend4->Draw();

	f->cd(3);
	hist_m23_down_neg->GetXaxis()->SetTitle(" m23 [MeV^{2}] ");
	hist_m23_down_neg->SetLineColor(2);
	hist_m23_down_neg->DrawClone();
	hist_m23_up_neg->SetLineColor(1);
	hist_m23_up_neg->DrawClone("same");
	auto legend5 = new TLegend(0.1,0.7,0.48,0.9);
	legend5->SetHeader("m_{13} histogram");
	legend5->AddEntry("hist_m23_up_neg","m_{23} MAGNET UP", "l");
	legend5->AddEntry("hist_m23_down_neg","m_{23} MAGNET DOWN", "l");
	legend5->SetTextSize(0.03);
	legend5->Draw();


	//ordering the masses of the two resonances

	auto ordered_up = dalitz_up.Filter("(m13 > 3474496 + 100000) || (m13 < 3474496 - 100000)")
				.Filter("(m23 > 3474496 + 100000) || (m23 < 3474496 - 100000)")
				.Filter("(m13 > 11600000 + 100000) || (m13 < 11600000 - 100000)")
				.Filter("(m23 > 11600000 + 100000) || (m23 < 11600000 - 100000)")
				.Define("roLow" , "std::min(m23,m13)")
				.Define("roHigh", "std::max(m23,m13)");
	auto ordered_down = dalitz_down.Filter("(m13 > 3474496 + 1000) || (m13 < 3474496 - 1000)")
				.Filter("(m23 > 3474496 + 100000) || (m23 < 3474496 - 100000)")
				.Filter("(m13 > 11600000 + 100000) || (m13 < 11600000 - 100000)")
				.Filter("(m23 > 11600000 + 100000) || (m23 < 11600000 - 100000)")
				.Define("roLow" , "std::min(m23,m13)")
				.Define("roHigh", "std::max(m23,m13)");

	auto filterGraph = ordered_up.Graph("m13","m23");
	auto s = new TCanvas("s","Dalitz plot without resonance");
        auto j = new TPad("j","j",0,0,1,1);
        j->Divide(1,1,0.01,0.01);
        j->Draw();
        j->cd(1);
	filterGraph->DrawClone("AP");

	auto Rolow = ordered_up.Filter("Bcharge == 1").Histo1D({"#rho_{low+}","",20u,0.8e6,28e5},"roLow");
	auto Rolow_neg = ordered_up.Filter("Bcharge == -1").Histo1D({"#rho_{low-}","",20u,0.8e6,28e5},"roLow");
	auto aa = ordered_down.Filter("Bcharge == 1").Histo1D({"","",20u,0.8e6,28e5},"roLow");
	auto bb = ordered_down.Filter("Bcharge == -1").Histo1D({"","",20u,0.8e6,28e5},"roLow");

	Rolow->Add(aa.GetPtr(),1); Rolow_neg->Add(bb.GetPtr(),1);

	auto RoHigh = ordered_up.Filter("Bcharge == 1").Histo1D({"#rho_{High+}","",26u, 2e6, 15e6},"roHigh");
	auto RoHigh_neg = ordered_up.Filter("Bcharge == -1").Histo1D({"#rho_{High-}","",26u, 2e6, 15e6},"roHigh");
	auto cc = ordered_down.Filter("Bcharge == 1").Histo1D({"","",26u,2e6,15e6},"roHigh");
	auto dd = ordered_down.Filter("Bcharge == -1").Histo1D({"","",26u,2e6,15e6},"roHigh");

	RoHigh->Add(cc.GetPtr(),1); RoHigh_neg->Add(dd.GetPtr(),1);
	//ordered_up.Display({"m13","m12","roHigh", "roLow"})->Print();

	auto ordGraph_up = ordered_up.Filter("Bcharge == 1").Graph("roLow","roHigh");
	auto ordGraph_down = ordered_down.Filter("Bcharge == 1").Graph("roLow","roHigh");
	auto ordGraph_negup = ordered_up.Filter("Bcharge == -1").Graph("roLow","roHigh");
	auto ordGraph_negDown = ordered_down.Filter("Bcharge == -1").Graph("roLow","roHigh");


	auto g = new TCanvas("g", "Mass Order", 1000, 1000);
	auto h = new TPad("h","h",0,0,1,1);
	h->Divide(2,2,0.01,0.01);
	h->Draw();
	h->cd(1);
	//histogram K+k- (Low)
	Rolow_neg->SetLineColor(kRed);
	Rolow_neg->GetXaxis()->SetTitle(" #rho_{low} [MeV**2] ");
	Rolow_neg->DrawClone();
	Rolow->SetLineColor(kBlack);
	Rolow->DrawClone("same");

	auto legend1 = new TLegend(0.1,0.7,0.48,0.9);
	legend1->SetHeader("#rho_{low} histogram");
	legend1->AddEntry("#rho_{low+}","#rho_{low} from B +", "l");
	legend1->AddEntry("#rho_{low-}", "#rho_{low} from B -", "l");
	legend1->SetTextSize(0.03);
	legend1->Draw();

	h->cd(2);
	//histogram K+K- (High)
	RoHigh->SetLineColor(kBlack);
	RoHigh->GetXaxis()->SetTitle(" #rho_{high} [MeV**2] ");
	RoHigh->DrawClone();
	RoHigh_neg->SetLineColor(kRed);
	RoHigh_neg->DrawClone("same");

	auto legend2 = new TLegend(0.1,0.7,0.48,0.9);
	legend2->SetHeader("#rho_{High} histogram");
	legend2->AddEntry("#rho_{High+}","#rho_{high} from B+", "l");
	legend2->AddEntry("#rho_{High-}","#rho_{high} from B-", "l");
	legend2->SetTextSize(0.03);
	legend2->Draw();

	h->cd(3);
	ordGraph_up->SetTitle(" Orderer dalitz magnet UP; invariant K+K- (Low); K+K- (High) ");
	ordGraph_up->DrawClone("AP");
	ordGraph_down->DrawClone("P");
	h->cd(4);
	ordGraph_negup->SetTitle(" Orderer dalitz magnet Down; invariant K+K-(Low); K+K- (High) ");
	ordGraph_negup->DrawClone("AP");
	ordGraph_negDown->DrawClone("P");

	//Save the files with the new data
	ordered_up.Snapshot("DecayTree", "Evt_Selected/ordered_up.root",{"invMass","H1_Charge", "H2_Charge","H3_Charge","Bcharge","roLow",
		"roHigh","m12" , "m13" , "m23"});
	ordered_down.Snapshot("DecayTree", "Evt_Selected/ordered_down.root",{"invMass","H1_Charge", "H2_Charge","H3_Charge","Bcharge","roLow",
		"roHigh","m12" , "m13" , "m23"});


        auto histBpos1 = ordered_up.Filter("Bcharge == 1").Histo1D({"histB1","recostructed Mass B+ meson",128u,4900,6900},"invMass");
        auto histBneg1 = ordered_up.Filter("Bcharge == -1").Histo1D({"histB2","recostructed Mass B- meson",128u,4900,6900},"invMass");
        auto histBpos2 = ordered_down.Filter("Bcharge == 1").Histo1D({"histB3","recostructed Mass B+ meson",128u,4900,6900},"invMass");
        auto histBneg2 = ordered_down.Filter("Bcharge == -1").Histo1D({"histB4","recostructed Mass B- meson",128u,4900,6900},"invMass");


        auto positiveBup = histBpos1->GetEntries();
        auto negBup = histBneg1->GetEntries();
        auto positiveBdown = histBpos2->GetEntries();
        auto negBdown = histBneg2->GetEntries();

        //Compute the asymmetry for the events selected
        float asymmetry1 = (positiveBup - negBup)/(positiveBup + negBup);
        float asymmetry2 = (positiveBdown - negBdown)/(positiveBdown + negBdown);

        float sigmaUp = TMath::Sqrt((1- TMath::Power(asymmetry1,2))/(positiveBup + negBup));
        float sigmaDown = TMath::Sqrt((1- TMath::Power(asymmetry2,2))/(positiveBdown + negBdown));

        std::cout << "- - -" << std::endl;
	std::cout << "Global asymmetry rejecting charm resonance " << std::endl;
        std::cout << "Asymmetry for MAGNET UP: " << asymmetry1 << " +/- " << sigmaUp << std::endl;
        std::cout << "Asymmetry for MAGNET DOWN: " << asymmetry2 << " +/- " << sigmaDown <<std::endl;

	return 0;
}

