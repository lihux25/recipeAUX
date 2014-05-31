// Source file for the Oxbridge Stransverse Mass Library -- oxbridgekinetics.
// See http://www.hep.phy.cam.ac.uk/~lester/mt2/index.html
// Authors: Christopher Lester and Alan Barr


#include "recipeAUX/OxbridgeMT2/interface/SUSYPhys_Mt2_222_Calculator.h"
#include "recipeAUX/OxbridgeMT2/interface/Mt2Units.h"
#include "Minuit2/MnUserParameterState.h"
#include "Minuit2/MnMigrad.h"
#include "Minuit2/FunctionMinimum.h"

#include <math.h>
#include <string>

using ROOT::Minuit2::MnUserParameters;
using ROOT::Minuit2::MnMigrad;
using ROOT::Minuit2::FunctionMinimum;

namespace Mt2 {

  /*
    
  stealing parts of files in
  
  http://isscvs.cern.ch/cgi-bin/cvsweb.cgi/offline/PhysicsAnalysis/SUSYPhys/SUSYPhysUser/src/?cvsroot=atlas
  
  like
  
  http://isscvs.cern.ch/cgi-bin/cvsweb.cgi/offline/PhysicsAnalysis/SUSYPhys/SUSYPhysUser/src/SusyUserExampleHistTool.cxx?rev=1.4;content-type=text%2Fplain;cvsroot=atlas
  
  and
  
  http://isscvs.cern.ch/cgi-bin/cvsweb.cgi/~checkout~/offline/PhysicsAnalysis/SUSYPhys/SUSYPhysUser/src/mT2Fcn.cxx?rev=1.1;content-type=text%2Fplain;cvsroot=atlas
  
  */
  double SUSYPhys_Mt2_222_Calculator::mt2_222(const TwoVector& visA, 
					      const TwoVector& visB,
					      const TwoVector& ptmiss, 
					      const double mEachInvisible){
 
    mT2Fcn theFCN(ptmiss.px(),
		  ptmiss.py(),
		  mEachInvisible,
		  visA.px(),
		  visA.py(),
		  visB.px(),
		  visB.py());
    
    double guessx = 0.5*(ptmiss.px());
    double guessy = 0.5*(ptmiss.py());

    MnUserParameters upar;
    upar.Add("etx", guessx, 100.0*Mt2::GeV); 
    upar.Add("ety", guessy, 100.0*Mt2::GeV);
    
    MnMigrad migrad(theFCN, upar, 2);
    
    FunctionMinimum min = migrad(0, Mt2::GeV);
    
    double best = min.Fval();

    if (best>=0) {
      return sqrt(best);
    } else {
      // something went wrong -- maybe just rounding errors;
      return sqrt(fabs(best));
    }
    
  }
  

  double  SUSYPhys_Mt2_222_Calculator::mT2Fcn::operator()(const std::vector<double>& par) const {
    
    double qT1x = par[0];
    double qT1y = par[1];
    
    double qT2x = theExmiss - qT1x;
    double qT2y = theEymiss - qT1y;
    
    double ETj1 = hypot(thePT1x,thePT1y);    // This is where this method assumes that the visible particles are massless!
    double ETj2 = hypot(thePT2x,thePT2y);    // This is where this method assumes that the visible particles are massless!    
    
    double ETchi1 = sqrt(qT1x*qT1x + qT1y*qT1y + theMchi*theMchi);
    double ETchi2 = sqrt(qT2x*qT2x + qT2y*qT2y + theMchi*theMchi);
    
    double mTsq1 = theMchi*theMchi + 2.0*(ETj1*ETchi1 - (thePT1x*qT1x + thePT1y*qT1y));  // This is where this method assumes that the visible particles are massless! Compare with Basic_Mt2_332_Calculator.cpp
    double mTsq2 = theMchi*theMchi + 2.0*(ETj2*ETchi2 - (thePT2x*qT2x + thePT2y*qT2y));  // This is where this method assumes that the visible particles are massless! Compare with Basic_Mt2_332_Calculator.cpp
    
    return fmax(mTsq1,mTsq2);
    
  }


} // end of Mt2 Namespace
