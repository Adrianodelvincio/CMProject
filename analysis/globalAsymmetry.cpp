#include <iostream>
#include <Math/Vector4D.h>
#include "fit.h"


int globalAsymmetry(){
        //load the two datafile
        auto fileDown = "processed_Down.root";
        auto fileUp = "processed_Up.root";

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

	//Use the cristal ball/cruijff function as a model
	gInterpreter->GenerateDictionary("Functions", "analysis/fit.h");
	TF1 *func1 = new TF1("fit_Bmass", Cruijff,5000,5400,6);
	TF1 *func2 = new TF1("4body",fourBodybackground,5039,5300,5);
	TF1 *f = new TF1("Model",model,5000,6400,13);
	TF1 *fexp = new TF1("fexp", "[0]*TMath::Exp(-[1]*(x-5000))",5038,6400);


	//set the initial values
	func1->SetParNames("sigma_L", "sigma_R", "mean", "alpha_L", "alpha_R", "Norm");
	func2->SetParNames("m0","c","p","leftBoundary","Norm");
	f->SetParNames("#sigma_{L}", "#sigma_{R}", "mean", "#alpha_{L}", "#alpha_{R}", "N_{Cruijff}","m_{0}","c","p","leftBoundary","N_{4b}");
	f->SetParName(11,"#lambda");
	f->SetParName(12,"N_{exp}");
	f->SetLineStyle(1);
	f->SetLineColor(kRed);
	fexp->SetParNames("norm","lambda");
	Double_t params[13] = 	{18.61,  //sigma_L
				12.32,	 //sigma_R
				5288,	 //mean
				0.0988,	 //alpha_L
				0.2111,	 //alpha_R
				1900,	 //Normalization Cruijff
				1.524e04,//m0
				15.08,	 //c
				1.861e04,//p
				5008,	 //leftBoundary
				0.6e-06, //Normalization 4b
				0.0006,  //lambda
				50};	 //Normalization exponential
	f->SetParameters(params);


	//Now remove the event above and below the invariant mass recostructed, and select B+ B-
	//f->FixParameter();

	//PLOTS
	auto c = new TCanvas("c","Invariant mass B meson",1100,900);
	auto p = new TPad("p","p",0,0,1,1);
	p->Divide(2,1,0.02,0.02);
	p->Draw();

	p->cd(1);
	gStyle->SetOptStat(0);
	gStyle->SetOptFit(1);

	histUp->Fit("Model","","same",5035,6300);
	//histUp->DrawClone("same");
	histUp->DrawClone("E1 X0");

	TF1 *fitUp = histUp->GetFunction("Model");

	double p0 = fitUp->GetParameter(0);
	double p1 = fitUp->GetParameter(1);
	double p2 = fitUp->GetParameter(2);
	double p3 = fitUp->GetParameter(3);
	double p4 = fitUp->GetParameter(4);
	double p5 = fitUp->GetParameter(5);
	func1->SetParameters(p0,p1,p2,p3,p4,p5);
	func1->SetLineColor(kBlue);
	func1->SetLineStyle(2);
	func1->DrawClone("same");

	double p6 = fitUp->GetParameter(6);
	double p7 = fitUp->GetParameter(7);
	double p8 = fitUp->GetParameter(8);
	double p9 = fitUp->GetParameter(9);
	double p10 = fitUp->GetParameter(10);
	func2->SetParameters(p6,p7,p8,p9,p10);
	func2->SetLineColor(kGreen);
	func2->SetLineStyle(2);
	func2->DrawClone("same");

	double p11 = fitUp->GetParameter(11);
	double p12 = fitUp->GetParameter(12);
	fexp->SetParameters(p12,p11);
	fexp->SetLineColorAlpha(kOrange,1);
	fexp->SetLineStyle(2);
	fexp->DrawClone("same");

	p->cd(2);
	gStyle->SetOptStat(0);
	gStyle->SetOptFit(1);
	histDown->Fit("Model","","same",5035,6300);
	//histDown->DrawClone();
	histDown->DrawClone("E1 X0");
	TF1 *fitDw = histDown->GetFunction("Model");
        double pp11 = fitDw->GetParameter(11);
        double pp12 = fitDw->GetParameter(12);
        fexp->SetParameters(pp12,pp11);
        fexp->SetLineColorAlpha(kOrange,1);
        fexp->SetLineStyle(2);
        fexp->DrawClone("same");

	double pp0 = fitDw->GetParameter(0);
        double pp1 = fitDw->GetParameter(1);
        double pp2 = fitDw->GetParameter(2);
        double pp3 = fitDw->GetParameter(3);
        double pp4 = fitDw->GetParameter(4);
        double pp5 = fitDw->GetParameter(5);
        func1->SetParameters(pp0,pp1,pp2,pp3,pp4,pp5);
        func1->SetLineColor(kBlue);
        func1->SetLineStyle(2);
        func1->DrawClone("same");

        double pp6 = fitDw->GetParameter(6);
        double pp7 = fitDw->GetParameter(7);
        double pp8 = fitDw->GetParameter(8);
        double pp9 = fitDw->GetParameter(9);
        double pp10 = fitDw->GetParameter(10);
        func2->SetParameters(pp6,pp7,pp8,pp9,pp10);
        func2->SetLineColor(kGreen);
        func2->SetLineStyle(2);
        func2->DrawClone("same");


	double mesonBmass = 5279.15; // MeV/c^2

	auto mass_up_selection = selectUp.Filter("invMass <= 5279.15 + 100")
				.Filter("invMass >= 5279.15 - 100")
				.Define("Bcharge","(H1_Charge + H2_Charge + H3_Charge)");
	auto mass_down_selection = selectDown.Filter("invMass <= 5279.15 + 100")
				.Filter("invMass >= 5279.15 - 100")
				.Define("Bcharge","(H1_Charge + H2_Charge + H3_Charge)");

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
