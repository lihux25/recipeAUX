// Source file for the Oxbridge Stransverse Mass Library -- oxbridgekinetics.
// See http://www.hep.phy.cam.ac.uk/~lester/mt2/index.html
// Authors: Christopher Lester and Alan Barr

#include "recipeAUX/OxbridgeMT2/interface/Mt2Calculators.h"
#include "recipeAUX/OxbridgeMT2/interface/Mt2Units.h"
#include "Minuit2/MnUserParameterState.h"
#include "Minuit2/MnSimplex.h"
#include "Minuit2/MnMinimize.h"
#include "Minuit2/MnStrategy.h"
#include "Minuit2/FunctionMinimum.h"
#include "recipeAUX/OxbridgeMT2/interface/mT2Fcn_4441_a.h"

using ROOT::Minuit2::MnUserParameters;
using ROOT::Minuit2::MnSimplex;
using ROOT::Minuit2::MnMinimize;
using ROOT::Minuit2::MnStrategy;
using ROOT::Minuit2::FunctionMinimum;

using std::cout;
using std::endl;

namespace Mt2 {

  template <class T> T smallestPositive(const T a, const T b) {
    if (a<=b) {
      if (a>=0) {
	return a;
      }
      return b; // b may not be positive, but it's the most positive thing remaining!
    } else {
      // b<=a
      if (b>=0) {
	return b;
      }
      return a; // a may not be positive, but it's the most positive thing remaining!
    }
  }

  

  double Basic_Mt2_4441_Calculator::
  mt2_4441(const Mt2::LorentzVector& pxpypzeAlpha,  // 4 d.o.f. 
	   const Mt2::LorentzVector& pxpypzeBeta,   // 4 d.o.f. 
	   const Mt2::LorentzVector& pxpypzeG,      // 4 d.o.f.
	   const double rootS,    /* eg 14 TeV */ // 1 d.o.f
	   const double mChi,
	   const double sign) {
    MnUserParameters upar;
    upar.Add("mpqlam", 0.5, 0.1, // next come the limits
	     0.0, 1.0);
    upar.Add("mhlam",  0.5, 0.1, // next come the limits
	     0.0, 1.0);
    const int highQuality=2;
    MnStrategy mnStrategy(highQuality);

    mT2Fcn_4441_a theFCN_4441_a(pxpypzeAlpha,     
				pxpypzeBeta,     
				pxpypzeG,         
				rootS,
				mChi,
				sign,
				this->debug());

    MnMinimize minimize(theFCN_4441_a, upar, mnStrategy);
    FunctionMinimum min = minimize(0, Mt2::GeV);

    const double mt2Candidate = min.Fval();
    return mt2Candidate;
  }

double Basic_Mt2_4441_Calculator::
  mt2_4441(const Mt2::LorentzVector& pxpypzeAlpha,      // 4 d.o.f. 
	   const Mt2::LorentzVector& pxpypzeBeta,       // 4 d.o.f. 
	   const Mt2::LorentzVector& pxpypzeG,          // 4 d.o.f.
	   const double rootS,    /* eg 14 TeV */     // 1 d.o.f
	   const double mChi) {
    
    const Mt2::LorentzVector pxpypzeAllVisible = pxpypzeG + pxpypzeAlpha + pxpypzeBeta;
    static const Mt2::LorentzVector pxpypzeLambda(Mt2::LorentzVector::InitEPxPyPz(1,0,0,0));
    const Mt2::LorentzVector pxpypzeF = 
      pxpypzeLambda*rootS - pxpypzeAlpha - pxpypzeBeta -  pxpypzeG;


    if (pxpypzeAlpha.m2()<0 ||
	pxpypzeBeta.m2()<0 ||
	pxpypzeG.m2()<0 ||
	pxpypzeAllVisible.m2() >= rootS*rootS ||
	pxpypzeF.m2()<0) {
      // there is definitely something bad about the input values
      return -1;
    }

//    double mt2Candidate[1];
    double mt2Candidate[2];
    for (int i=0; i<=1; ++i) {
      const double sign = (i==0 ? +1 : -1);
      mt2Candidate[i] = mt2_4441(pxpypzeAlpha,     
				 pxpypzeBeta,     
				 pxpypzeG,         
				 rootS,
				 mChi,
				 sign);
    }
    
    if (this->debug()) {
      cout << "mt2+ mt2- " << mt2Candidate[0] << " " << mt2Candidate[1] << endl;
    }

    const double mt2 = smallestPositive(mt2Candidate[0],
					mt2Candidate[1]);
    
    return mt2;
  }
  
} // end of Mt2 Namespace
