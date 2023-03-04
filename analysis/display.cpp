#include <iostream>


int display(){
	//load the two datafile
	auto fileDown = "B2HHH_MagnetDown.root";
	auto fileUp = "B2HHH_MagnetUp.root";

	ROOT::RDataFrame rdf_up("DecayTree", fileUp);
	ROOT::RDataFrame rdf_down("DecayTree", fileDown);
	float range = 4000;

	//Print column name of the data tree
	auto colNames =  rdf_up.GetColumnNames();
	for (auto &&colname : colNames) std::cout << colname << std::endl;

	//look at the histogram
	auto h_upx = rdf_up.Histo1D({"h_upx","First particle px",64u,-range,range},"H1_PX");
	auto h_upy = rdf_up.Histo1D({"h_upy","First particle py",64u,-100.,100.},"H1_PY");
	auto h_downy = rdf_down.Histo1D({"h_downy","First particle py",64u,-100.,100.},"H1_PY");
	auto h_downx = rdf_down.Histo1D({"h_downx","First particle px",64u,-range,range},"H1_PX");
	h_upx->SetMarkerSize(0.5);

	h_upx->SetYTitle("Counts");
	h_upy->SetYTitle("Counts");
	h_downx->SetYTitle("Counts");
	h_downy->SetYTitle("Counts");

	auto c = new TCanvas("c","c", 600, 600);
	int nx = 2, ny = 2;
	auto p = new TPad("p","",0,0,1,1);
	p->Divide(nx,ny,0.01,0.01);
	p->Draw();
	p->cd(1);
	h_upx->DrawClone();
	p->cd(2);
	h_upy->DrawClone();
	p->cd(3);
	h_downx->DrawClone();
	p->cd(4);
	h_downy->DrawClone();

	auto BselectionUp = rdf_up.Define("Bcharge","(H1_Charge + H2_Charge + H3_Charge)");
	auto BselectionDown = rdf_down.Define("Bcharge","(H1_Charge + H2_Charge + H3_Charge)");

	BselectionUp.Display({"H1_Charge", "H2_Charge","H3_Charge", "Bcharge"},25,10)->Print();
        BselectionDown.Display({"H1_Charge", "H2_Charge","H3_Charge","Bcharge"},25,10)->Print();


	//look at the histogram
	auto h1_upx = rdf_up.Histo1D({"h_upx","First particle px",64u,-range,range},"H2_PX");
	auto h1_upy = rdf_up.Histo1D({"h_upy","First particle py",64u,-100.,100.},"H2_PY");
	auto h1_downy = rdf_down.Histo1D({"h_downy","First particle py",64u,-100.,100.},"H2_PY");
	auto h1_downx = rdf_down.Histo1D({"h_downx","First particle px",64u,-range,range},"H2_PX");
	h1_upx->SetMarkerSize(0.5);

	h1_upx->SetYTitle("Counts");
	h1_upy->SetYTitle("Counts");
	h1_downx->SetYTitle("Counts");
	h1_downy->SetYTitle("Counts");

	auto c1 = new TCanvas("c1","c1", 600, 600);
	auto p1 = new TPad("p1","",0,0,1,1);
	p1->Divide(nx,ny,0.01,0.01);
	p1->Draw();
	p1->cd(1);
	h1_upx->DrawClone();
	p1->cd(2);
	h1_upy->DrawClone();
	p1->cd(3);
	h1_downx->DrawClone();
	p1->cd(4);
	h1_downy->DrawClone();
	return 0;
}
