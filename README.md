recipeAUX
===============
Works in CMSSW_7_1_0_pre8

*Setup LHAPDF*

- cd recipeAUX
- source LHAPDF.csh

It will setup enviromental variable:
- $LHAPDF\_INCLUDE : For include path of LHAPDF header
- $LHAPDF_LIBDIR : For shared library location of LHAPDF
- $LHAPDF_DATA_PATH : The dataset location of PDFsets. On LPC/Lxplus, it
included all the current PDFSet. On local machines, it downloads NNPDF30_nlo_as_0118,
CT10nlo and MMHT2014nlo68cl so far.
