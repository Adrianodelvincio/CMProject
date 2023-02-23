#include <iostream>
#include <Math/Vector4D.h>

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

	auto histUp = selectUp.Histo1D({"mass B","B mass",128u,4900,6900},"invMass");
	auto histDown = selectDown.Histo1D({"mass B","B mass",128u,4900,6900},"invMass");
	// try stronger selection
	auto histUp1 = selectUp1.Histo1D({"mass B","B mass",128u,4900,6900},"invMass");
	auto histDown1 = selectDown1.Histo1D({"mass B","B mass",128u,4900,6900},"invMass");

	//Now remove the event above and below the invariant mass recostructed, and select B+ B-
	// it's important to define a way to select a good range
	double mesonBmass = 5279.15; // MeV/c^2

	auto mass_up_selection = selectUp.Filter("invMass <= 5279.15 + 100")
				.Filter("invMass >= 5279.15 - 100")
				.Define("Bcharge","(H1_Charge + H2_Charge + H3_Charge)");
	auto mass_down_selection = selectDown.Filter("invMass <= 5279.15 + 100")
				.Filter("invMass >= 5279.15 - 100")
				.Define("Bcharge","(H1_Charge + H2_Charge + H3_Charge)");

	auto histUp2 = mass_up_selection.Histo1D({"mass B","B mass",128u,4900,6900},"invMass");
	auto histDown2 = mass_down_selection.Histo1D({"mass B","B mass",128u,4900,6900},"invMass");

	auto c = new TCanvas("c","c",900,900);
	auto p = new TPad("p","p",0,0,1,1);
	p->Divide(2,2,0.01,0.01);
	p->Draw();
	p->cd(1); histUp->DrawClone("PLC");
	histUp1->DrawClone("same PLC");
	p->cd(2); histDown->DrawClone("PLC");
	histDown1->DrawClone("same PLC");
	p->cd(3);
	histUp2->DrawClone();
	p->cd(4);
	histDown2->DrawClone();

	//Identify B+ and B- by the kaon charge

	auto histBpos1 = mass_up_selection.Filter("Bcharge == 1").Histo1D({"histB","recostructed Mass B+ meson",128u,4900,6900},"invMass");
	auto histBneg1 = mass_up_selection.Filter("Bcharge == -1").Histo1D({"histB","recostructed Mass B- meson",128u,4900,6900},"invMass");
	auto histBpos2 = mass_down_selection.Filter("Bcharge == 1").Histo1D({"histB","recostructed Mass B+ meson",128u,4900,6900},"invMass");
	auto histBneg2 = mass_down_selection.Filter("Bcharge == -1").Histo1D({"histB","recostructed Mass B- meson",128u,4900,6900},"invMass");


	auto e = new TCanvas("e","e",900,900);
	auto f = new TPad("f","f",0,0,1,1);
	f->Divide(2,2,0.01,0.01);
	f->Draw();
	f->cd(1); histBpos1->DrawClone();
	f->cd(2); histBneg1->DrawClone();

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
