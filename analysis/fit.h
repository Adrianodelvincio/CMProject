#ifndef FIT_H
#define FIT_H
#include <TMath.h>
#include <Math/PdfFuncMathCore.h>
#include <Math/ProbFuncMathCore.h>

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
	double arg2 = 0;
	double arg3 = 0;
	if(x[0] > par[3] && x[0] < par[3] + par[0]){
		arg1 = xx*TMath::Power(1 - TMath::Power(xx/par[0],2), par[2]) * TMath::Exp(par[1] * (1 - TMath::Power(xx/par[0],2)));
		arg2 = 1;
		arg3 = 1;
	}
	return par[4]*arg3*arg1*arg2;
}


/*
//par[0]-> alpha, par[1] -> sigma , par[2]-> mean, par[3]-> n, par[4]-> normalization
double CrystallBall(double *x, double *par){
	double arg = 0;
	double arg2 = (x[0] - par[2]);
	if(par[1] != 0){
		if(arg2/par[1] > -par[0]){
			arg = TMath::Exp(-TMath::Power(arg2,2)/(2*par[1]*par[1]));
		}
		else{
			double a = TMath::Power((par[3]/TMath::Abs(par[0])),par[3])*TMath::Exp(-TMath::Abs(par[0]*par[0])/2);
			double b = par[3]/TMath::Abs(par[0]) - TMath::Abs(par[0]);
			double c = (par[3]/TMath::Abs(par[0]))*(1/(par[3]-1))*TMath::Exp(-TMath::Abs(par[0]*par[0]/2));
			double d = TMath::Sqrt(TMath::Pi()/2)*(1 + TMath::Erf(TMath::Abs(par[0]/TMath::Sqrt(2))));
			double N = 1/(par[2]*(c + d));
			arg = a*TMath::Power((b - arg2/par[1]),-par[3]);
		}
	}
	return par[4]*arg;
}
*/

#endif
