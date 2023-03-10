#include <iostream>
#include <Math/Vector4D.h>
#include "fit.h"


int globalAsymmetry(){
        //load the two datafile
        auto fileDown = "../processed_Down.root";
        auto fileUp = "../processed_Up.root";

        ROOT::RDataFrame rdf_up("DecayTree", fileUp);
        ROOT::RDataFrame rdf_down("DecayTree", fileDown);

	auto colnames = rdf_up.GetColumnNames();
	for(auto &&colname : colnames) std::cout<< colname << std::endl;

	//After the Muons
	int entriesUp = 2505220 ; int entriesDown = 3806297;
	std::cout << "Number of Entries magnet UP: " << entriesUp << std::endl;
	std::cout << "Number of Entries magnet DOWN: " << entriesDown << std::endl;

	//Now start with the selection
	auto SelectionUp = rdf_up.Filter("H1_ProbPi <= 0.5")
	.Filter("H2_ProbPi <= 0.5")
	.Filter("H3_ProbPi <= 0.5")
	.Filter("H1_ProbK >= 0.5")
	.Filter("H2_ProbK >= 0.5")
	.Filter("H3_ProbK >= 0.5").Count();
	auto SelectionDown = rdf_down.Filter("H1_ProbPi <= 0.5")
	.Filter("H2_ProbPi <= 0.5")
	.Filter("H3_ProbPi <= 0.5")
	.Filter("H1_ProbK >= 0.5")
	.Filter("H2_ProbK >= 0.5")
	.Filter("H3_ProbK >= 0.5").Count();

	//after filtering the pions
	auto pifilterUp = SelectionUp.GetValue();
	auto pifilterDown = SelectionDown.GetValue();
	std::cout << "After Pi filter magnet UP: " << pifilterUp << std::endl;
	std::cout << "After Pi filter magnet DOWN: " << pifilterDown << std::endl;
	auto lost = [](int total, int final ){ return (float) final/total;};
	std::cout << "Events lost Magnet UP: " << 1 - lost(entriesUp,pifilterUp) << std::endl;
	std::cout << "Events lost Magnet DOWN: " << 1 - lost(entriesDown,pifilterDown) << std::endl;

	//plot the probabilities
	auto upKaon1 = rdf_up.Histo1D({"Pion","Prob. is pion",128u,0.,1.},"H1_ProbPi");
	auto upKaon2 = rdf_up.Histo1D({"Pion","Prob. is pion",128u,0.,1.},"H2_ProbPi");
	auto upKaon3 = rdf_up.Histo1D({"Pion","Prob. is pion",128u,0.,1.},"H3_ProbPi");

	auto downKaon1 = rdf_up.Histo1D({"Pion","Prob. is pion",128u,0.,1.},"H1_ProbPi");
	auto downKaon2 = rdf_up.Histo1D({"Pion","Prob. is pion",128u,0.,1.},"H2_ProbPi");
	auto downKaon3 = rdf_up.Histo1D({"Pion","Prob. is pion",128u,0.,1.},"H3_ProbPi");
	/*
	auto b = new TCanvas("b", "b", 900, 900);
	auto d = new TPad("d","d",0,0,1,1);
	d->Draw();

	d->Divide(2,2,0.01,0.01);
	d->cd(1);
	upKaon1->DrawClone("PLC");
	upKaon2->DrawClone("same PLC");
	upKaon3->DrawClone("same PLC");
	d->cd(2);
	downKaon1->DrawClone("PLC");
	downKaon2->DrawClone("same PLC");
	downKaon3->DrawClone("same PLC");
	*/

	//Now select the kaon event and plot invariant mass after the selection
	auto selectUp = rdf_up.Filter("H1_ProbPi <= 0.5")
	.Filter("H2_ProbPi <= 0.5")
	.Filter("H3_ProbPi <= 0.5")
	.Filter("H1_ProbK >= 0.5")
	.Filter("H2_ProbK >= 0.5")
	.Filter("H3_ProbK >= 0.5");

	auto selectDown =  rdf_down.Filter("H1_ProbPi <= 0.5")
	.Filter("H2_ProbPi <= 0.5")
	.Filter("H3_ProbPi <= 0.5")
	.Filter("H1_ProbK >= 0.5")
	.Filter("H2_ProbK >= 0.5")
	.Filter("H3_ProbK >= 0.5");

	auto selectUp1 = rdf_up.Filter("H1_ProbPi <= 0.2")
	.Filter("H2_ProbPi <= 0.2")
	.Filter("H3_ProbPi <= 0.2")
	.Filter("H1_ProbK >= 0.5")
	.Filter("H2_ProbK >= 0.5")
	.Filter("H3_ProbK >= 0.5");

	auto selectDown1 =  rdf_down.Filter("H1_ProbPi <= 0.2")
	.Filter("H2_ProbPi <= 0.2")
	.Filter("H3_ProbPi <= 0.2")
	.Filter("H1_ProbK >= 0.5")
	.Filter("H2_ProbK >= 0.5")
	.Filter("H3_ProbK >= 0.5");

	auto histUp = selectUp.Histo1D({"Bup","B mass",128u,4700,6400},"invMass");
	auto histDown = selectDown.Histo1D({"Bdown","B mass",128u,4700,6400},"invMass");

	//Now fit the invariant mass of the B meson
	//Use the cristal ball function to fit the mass!!!
	gInterpreter->GenerateDictionary("Functions", "fit.h");
	TF1 *func1 = new TF1("fit_Bmass", Cruijff,4700,6400,6);
	TF1 *func2 = new TF1("fit_Bmass2",CrystallBall,4700,6400,5);
	TF1 *func3 = new TF1("fit_Bmass3",fourBodybackground,4700,6400,5);
	//set the initial values
	func1->SetParameters(histUp->GetRMS(),histUp->GetRMS(),histUp->GetMean(),0,0,100);
	func1->SetParNames("sigma_L", "sigma_R", "mean", "alpha_L", "alpha_R", "norm");
	func2->SetParameters(1,2,histUp->GetRMS(),histUp->GetMean(),500);
	func2->SetParNames("alpha","n","sigma","mean","norm");
	func3->SetParameters(900,10,1,5000,1e-5);
	func3->SetParNames("m0","c","p");

	//Now remove the event above and below the invariant mass recostructed, and select B+ B-
	// it's important to define a way to select a good range
	double mesonBmass = 5279.15; // MeV/c^2

	auto mass_up_selection = selectUp.Filter("invMass <= 5279.15 + 100")
				.Filter("invMass >= 5279.15 - 100")
				.Define("Bcharge","(H1_Charge + H2_Charge + H3_Charge)");
	auto mass_down_selection = selectDown.Filter("invMass <= 5279.15 + 100")
				.Filter("invMass >= 5279.15 - 100")
				.Define("Bcharge","(H1_Charge + H2_Charge + H3_Charge)");

	//PLOTS

	auto c = new TCanvas("c","Invariant mass B meson",900,900);
	auto p = new TPad("p","p",0,0,1,1);

	p->Divide(2,1,0.02,0.02);
	p->Draw();

	p->cd(1);
	gStyle->SetOptStat(0);
	histUp->Fit("fit_Bmass3","L","same",5000,5200);
	gStyle->SetOptFit(1);
	histUp->DrawClone();
	//func3->DrawClone();
	p->cd(2);
	gStyle->SetOptStat(0);
	histDown->Fit("fit_Bmass3","L","same",5000,5200);
	gStyle->SetOptFit(1);
	histDown->DrawClone();

	//auto e = new TCanvas("e","Invariant mass B+/B-",900,900);
        //auto f = new TPad("f","f",0,0,1,1);
        //f->Divide(2,2,0.01,0.01);
        //f->Draw();

	//Identify B+ and B- by the kaon charge

	auto histBpos1 = mass_up_selection.Filter("Bcharge == 1").Histo1D({"histB1","recostructed Mass B+ meson",128u,4900,6900},"invMass");
	auto histBneg1 = mass_up_selection.Filter("Bcharge == -1").Histo1D({"histB2","recostructed Mass B- meson",128u,4900,6900},"invMass");
	auto histBpos2 = mass_down_selection.Filter("Bcharge == 1").Histo1D({"histB3","recostructed Mass B+ meson",128u,4900,6900},"invMass");
	auto histBneg2 = mass_down_selection.Filter("Bcharge == -1").Histo1D({"histB4","recostructed Mass B- meson",128u,4900,6900},"invMass");

	//f->cd(1); histBpos1->DrawClone();
	//f->cd(2); histBneg1->DrawClone();

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
	std::cout << "Asymmetry for MAGNET UP: " << asymmetry1 << " +/- " << sigmaUp << std::endl;
	std::cout << "Asymmetry for MAGNET DOWN: " << asymmetry2 << " +/- " << sigmaDown <<std::endl;

	//Now compute the uncertainty of the asymmetry
	return 0;
}
