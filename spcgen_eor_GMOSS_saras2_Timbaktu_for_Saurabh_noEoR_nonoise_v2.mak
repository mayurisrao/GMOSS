spcgen_eor_GMOSS_saras2_Timbaktu_for_Saurabh_noEoR_nonoise_v2: spcgen_eor_GMOSS_saras2_Timbaktu_for_Saurabh_noEoR_nonoise_v2.o \
	nrutil.o  adpint.o amoeba.o amotry.o besselk.o qromo.o midpnt.o polint.o fpoly.o lfit.o \
	qromo_alternate.o midpnt_alternate.o polint_alternate.o spline.o splint.o gaussj.o covsrt.o \
	refraction.o ran1.o cal_lst.o beam_definition_cos2_pattern.o precess.o gasdev.o
	gfortran -o spcgen_eor_GMOSS_saras2_Timbaktu_for_Saurabh_noEoR_nonoise_v2 \
	spcgen_eor_GMOSS_saras2_Timbaktu_for_Saurabh_noEoR_nonoise_v2.o nrutil.o besselk.o adpint.o fpoly.o lfit.o \
	amoeba.o amotry.o qromo.o midpnt.o polint.o qromo_alternate.o midpnt_alternate.o gaussj.o covsrt.o \
	polint_alternate.o spline.o splint.o refraction.o ran1.o cal_lst.o beam_definition_cos2_pattern.o precess.o gasdev.o\
	 -L/usr/local/lib \
	 -L/usr/lib \
	 -lm 
spcgen_eor_GMOSS_saras2_Timbaktu_for_Saurabh_noEoR_nonoise_v2.o : spcgen_eor_GMOSS_saras2_Timbaktu_for_Saurabh_noEoR_nonoise_v2.c
	gfortran -c spcgen_eor_GMOSS_saras2_Timbaktu_for_Saurabh_noEoR_nonoise_v2.c
nrutil.o : nrutil.c
	gfortran -c nrutil.c
cal_lst.o : cal_lst.c
	  gfortran -c cal_lst.c
amoeba.o : amoeba.c
	gfortran -c amoeba.c
amotry.o : amotry.c
	gfortran -c amotry.c
besselk.o : besselk.f
	gfortran -c besselk.f
adpint.o : adpint.f
	gfortran -c adpint.f
beam_definition_cos2_pattern.o : beam_definition_cos2_pattern.c
	gfortran -c beam_definition_cos2_pattern.c
precess.o : precess.f
	gfortran -c precess.f
qromo.o : qromo.c
	gfortran -c qromo.c
ran1.o : ran1.c
	gfortran -c ran1.c
refraction.o : refraction.c
	gfortran -c refraction.c
midpnt.o : midpnt.c
	gfortran -c midpnt.c
qromo_alternate.o : qromo_alternate.c
	gfortran -c qromo_alternate.c
midpnt_alternate.o : midpnt_alternate.c
	gfortran -c midpnt_alternate.c
polint.o : polint.c
	gfortran -c polint.c
polint_alternate.o : polint_alternate.c
	gfortran -c polint_alternate.c
spline.o : spline.c
	gfortran -c spline.c
splint.o : splint.c
	gfortran -c splint.c
gasdev.o : gasdev.c
	gfortran -c gasdev.c
fpoly.o : fpoly.c
	gfortran -c fpoly.c
lfit.o : lfit.c
	gfortran -c lfit.c
covsrt.o : covsrt.c
	gfortran -c covsrt.c
gaussj.o : gaussj.c
	gfortran -c gaussj.c

