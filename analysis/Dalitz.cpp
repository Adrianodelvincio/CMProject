#include <iostream>
#include <Math/Vector4D.h>

int Dalitz(){
        //load the two datafile
        auto fileDown = "Processed_Down";
        auto fileUp = "Processed_Up";

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
	.Filter("H3_ProbPi <= 0.5").Count();
	auto SelectionDown = rdf_down.Filter("H1_ProbPi <= 0.5")
	.Filter("H2_ProbPi <= 0.5")
	.Filter("H3_ProbPi <= 0.5").Count();

	//after filtering the pions
	auto pifilterUp = SelectionUp.GetValue();
	auto pifilterDown = SelectionDown.GetValue();
	std::cout << "After Pi filter magnet UP: " << pifilterUp << std::endl;
	std::cout << "After Pi filter magnet DOWN: " << pifilterDown << std::endl;
	auto lost = [](int total, int final ){ return (float) final/total;};
	std::cout << "Events lost Magnet UP: " << 1 - lost(entriesUp,pifilterUp) << std::endl;
	std::cout << "Events lost Magnet DOWN: " << 1 - lost(entriesDown,pifilterDown) << std::endl;
	return 0;
	}
