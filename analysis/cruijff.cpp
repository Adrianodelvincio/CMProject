#include<iostream>
#include <ROOT>

double Cruijff(double *x, double *par){
        double arg = 0;
        double arg2 = (x[0] - par[2]);

        if (par[0] != 0 && par[1] != 0){
                if( arg2 <= 0){
                        arg = par[5]*TMath::Exp(-TMath::Power(arg2,2) / (2*par[0]*par[0] + par[3] *TMath::Power(arg2,2)));
                }
                else{
                        arg = par[5]*TMath::Exp(-TMath::Power(arg2,2) / (2*par[1]*par[1] + par[4] *TMath::Power(arg2,2)));
                }
        }
        return arg;
}


void cruijff(){
	TF1 *func = new TF1("fit_Bmass",Cruijff,5200,5400,6);
	func->SetParameters(200,200,5300,1,10,10);
	func->DrawClone();
}
