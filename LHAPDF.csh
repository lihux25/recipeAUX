#!/bin/csh

if ($?CMSSW_BASE) then 
  echo $CMSSW_BASE
  setenv LHAPDF_INCLUDE `scram tool info lhapdf | grep ROOT_INCLUDE_PATH | cut -d= -f2`
  setenv LHAPDF_LIBDIR `scram tool info lhapdf | grep LIBDIR | cut -d= -f2`
  setenv LHAPDF_DATA_PATH /cvmfs/sft.cern.ch/lcg/external/lhapdfsets/current/
  echo $LHAPDF_INCLUDE 
  echo $LHAPDF_LIBDIR
  echo $LHAPDF_DATA_PATH
  return 1
endif

set PDFList=(NNPDF30_nlo_as_0118 CT10nlo MMHT2014nlo68cl)
if ( ! -d LHAPDF) then
  wget  http://www.hepforge.org/archive/lhapdf/LHAPDF-6.1.5.tar.gz -O- | tar xz -C ./
  mv LHAPDF-6.1.5 LHAPDF
  cd LHAPDF
  mkdir build
  ./configure --disable-python --prefix=$PWD/build
  make -j 4
  make install
  cd -
endif
setenv LHAPDF_INCLUDE $PWD/LHAPDF/build/include
setenv LHAPDF_LIBDIR $PWD/LHAPDF/build/lib
setenv LHAPDF_DATA_PATH $PWD/LHAPDF/build/share
setenv PATH $PWD/build/bin:$PATH
foreach PDF ($PDFList)
  echo $PDF
  if (! -d $LHAPDF_DATA_PATH/$PDF) then
    wget http://www.hepforge.org/archive/lhapdf/pdfsets/6.1/${PDF}.tar.gz -O- | tar xz -C $LHAPDF_DATA_PATH
  endif
end
