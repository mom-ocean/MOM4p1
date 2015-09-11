#These are the notes on how to make the MOM4p1 public release.
#The notes are intended only for the GFDL people responsible for making the MOM4p1 public release.  
#
#check outs
mkdir src
cd src
cvs co -r riga_prerelease2 atmos_bgrid atmos_coupled atmos_ebm atmos_fv_dynamics
cvs co -r riga_prerelease2 atmos_null atmos_param atmos_shared atmos_spectral
cvs co -r riga_prerelease2 coupler 
cvs co -r riga_prerelease2 ice_sis ice_param
cvs co -r riga_prerelease2 land_lad land_null land_param
cvs co -r riga_prerelease2 shared postprocessing preprocessing  tools
cvs co -r riga_prerelease2 ocean_shared
cvs co -r riga_prerelease2 mom4p1

cvs up -r mom4p1_riga_14dec2009_smg mom4p1/
cvs up -r mom4p1_riga_16dec2009_smg mom4p1/doc

#Test, test, test and more test. 

#tag up everythying

cvs tag -b mom4p1_pubrel_dec2009_nnz atmos_ebm/ atmos_fv_dynamics/ atmos_null/ atmos_param/ atmos_shared/ atmos_spectral/ coupler/ ice_param/ ice_sis/ land_lad/ land_null/ land_param/ ocean_shared/ postprocessing/ preprocessing/ shared/ tools/  atmos_bgrid/ atmos_coupled/ mom4p1/

cvs up -r  mom4p1_pubrel_dec2009_nnz atmos_ebm/ atmos_fv_dynamics/ atmos_null/ atmos_param/ atmos_shared/ atmos_spectral/ coupler/ ice_param/ ice_sis/ land_lad/ land_null/ land_param/ ocean_shared/ postprocessing/ preprocessing/ shared/ tools/  atmos_bgrid/ atmos_coupled/ mom4p1/


#test the new tag
cd ..
mv src src.old
mkdir src
cd src

cvs co -r mom4p1_pubrel_dec2009_nnz atmos_ebm atmos_fv_dynamics atmos_null atmos_param atmos_shared atmos_spectral coupler ice_param ice_sis land_lad land_null land_param ocean_shared postprocessing preprocessing shared tools  atmos_bgrid atmos_coupled mom4p1

#test, test, test

#cleanup
find . -name '*.pdf'  -exec rm -rf {} \;
find . -name '*.ps' -exec rm -rf {} \;
find . -name '*.jpg'   -exec rm -rf {} \;
\rm -rf postprocessing/analysis

#GNUize all F90 and c files

/home/nnz/bin/GNULicense.pl -f --dir=. --recursive

#This will touch ALL code, so you have to test, test, test!

# cvs ci or not ci ,  I wouldn't do it GNU!!

#Documentation
#Note that .pdf and .ps files are too big to be part of the release bundle (GFORGE has a limit of 8M for each file).
So these big files must go on GForge separately.
But you should leave the .html files in the original mom4p1/doc directory because they have relative links to the src code.
#Work on the doc/README check it and move it up to the ROOT dir.
cd doc
#Work on the doc/quickstart_guide.xml
#Generate the html and pdf from xml
/home/arl/bin/mkdocbk mom4_manual.xml

cvs ci README quickstart_guide.xml quickstart_guide.html quickstart_guide.pdf
mv README ../../


#Date tag (sticky) the whole thing
cd src
cvs tag mom4p1_pubrel_18dec2009_nnz *

#Also you may want to remove the CVS directories at the very last step when sure there are no more changes.
#find . -name 'CVS'  -exec rm -rf {} \;



#
#tar up
#
cd ../../
tar cvf mom4p1_pubrel_dec2009.tar mom4p1_pubrel_dec2009/src
tar rvf mom4p1_pubrel_dec2009.tar mom4p1_pubrel_dec2009/bin
tar rvf mom4p1_pubrel_dec2009.tar mom4p1_pubrel_dec2009/exp
gzip mom4p1_pubrel_dec2009.tar


