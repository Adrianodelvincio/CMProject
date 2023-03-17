#ifndef FIT_H
#define FIT_H
#include <TMath.h>
#include <Math/PdfFuncMathCore.h>
#include <Math/ProbFuncMathCore.h>
#include <iostream>

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

//par[0] -> alpha, par[1]-> n, par[2]-> sigma, par[3] -> mu, par[4]-> norm
double CrystallBall(double *x, double *par){
	return par[4]*ROOT::Math::crystalball_function(x[0],par[0] , par[1], par[2], par[3]);
}


// par[0]-> m0, par[1] = c, par[2] = p, par[3] initial, par[4]-> normalization
double fourBodybackground(double *x, double *par){
	double xx = x[0] - par[3];
	double arg1 = 0;
	if(x[0] > par[3]){
		arg1 = xx*TMath::Power(1 - TMath::Power(xx/par[0],2), par[2]) * TMath::Exp(par[1] * (1 - TMath::Power(xx/par[0],2)));
	}
	return par[4]*arg1;
}

// par[0] -> sigma_L
// par[1] -> sigma_R
// par[2] -> mean
// par[3] -> alpha_L
// par[4] -> alpha_R
// par[5] -> Norm_cruijff
// par[6] -> m0
// par[7] -> c
// par[8] -> p
// par[9] -> leftBoundary
// par[10] -> Norm 4b
// par[11] -> lambda
// par[12] -> N exp

double model(double *x,double *par){
	//Cruijff
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
	//4body
	double xx = x[0] - par[9];
        double arg1 = 0;
        arg1 = xx*TMath::Power(1 - TMath::Power(xx/par[6],2), par[8]) * TMath::Exp(par[7] * (1 - TMath::Power(xx/par[6],2)));
	//exponential
	double arg3 = 0;
	arg3 = par[12]*TMath::Exp(-par[11]*xx);
	std::cout << "arg3 : " << arg3 << " par[12] : " << par[12] << " par[11] : " << par[11] << " xx " << xx << " TMath::Exp(-par[11]*xx): "<< TMath::Exp(-par[11]*xx)  <<std::endl;
        return (arg + par[10]*arg1 +arg3);
}

#endif


