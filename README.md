# GMOSS
This repository contains C codes that generate GMOSS: "Global MOdel for the radio Sky Spectrum". see "https://export.arxiv.org/abs/1607.07453"
Some notes: 
There are three sets of codes
1. patch_integrate_besselk_10oct16_1_500 -- This is the actual code that does the fitting of the GMOSS model to input sky maps using physics described in the paper.
   Note this code as of the first version on Github is hardcoded to use fixed set of sky maps, and loop over HEALPIX pixels 1 to 500 -- please change to 1 to 3072 if using the original input data sets.

2. GMOSS_8aug17.c which uses already optimized GMOSS physical parameters (in the file GMOSS_params_allpix.txt) to generate spectra over a range of frequencies as specified in an input file.

3. spcgen_eor_GMOSS_saras2_Timbaktu_for_Saurabh_noEoR_nonoise_v2.c which generates expectation of sky spectra with GMOSS as the sky model, when observed by a mock instrument looking up at the sky with a certain antenna beam, having some noise properties, at a certain location (lat,lon). The specifics of the mock observation need to be modified as per the original version pushed on github. Cleaner revisions will be updated in due course which are less hardcoded!

There are also some supporting files, including numerical recipes (author acknowledges the original Numerical Recipes in C and F), sky maps (sources cited in the paper). 
