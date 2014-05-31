// Source file for the Oxbridge Stransverse Mass Library -- oxbridgekinetics.
// See http://www.hep.phy.cam.ac.uk/~lester/mt2/index.html
// Authors: Christopher Lester and Alan Barr

#include "recipeAUX/OxbridgeMT2/interface/Mt2Vectors.h"
#include <iomanip>
#include <sstream>

namespace Mt2{

Mt2Exception::Mt2Exception(const std::string & reason) throw() : m_reason(reason){
}

Mt2Exception::~Mt2Exception() throw() {}

const char* Mt2Exception::what() const throw() {
  return m_reason.c_str();
}

void LorentzTransverseVector::print(std::ostream& os) const{
  using std::setw;
  os << "Lor2Vec: x = " << setw(9) << px() 
     << ", y = " << setw(9) << py()
     << ", Et = " << setw(9) << Et()
     << ", mass = " << mass();
}

void TwoVector::print(std::ostream& os) const {
  using std::setw;
  os << "TwoVector: x = " << setw(9) << px()
     << ", y = " << setw(9) << py();
}

}

std::ostream& operator << (std::ostream& os , const Mt2::LorentzTransverseVector& v){
  v.print(os);
  return os;
}
                                                                                                    
std::ostream& operator << (std::ostream& os , const Mt2::TwoVector& v){
  v.print(os);
  return os;
}


