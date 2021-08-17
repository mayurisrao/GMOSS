patch_integrate_besselk_10oct16_1_500: patch_integrate_besselk_10oct16_1_500.o \
	nrutil.o  adpint.o amoeba.o amotry.o besselk.o qromo.o midpnt.o polint.o \
	qromo_alternate.o midpnt_alternate.o polint_alternate.o spline.o splint.o
	gfortran -o patch_integrate_besselk_10oct16_1_500 \
	patch_integrate_besselk_10oct16_1_500.o nrutil.o besselk.o adpint.o \
	amoeba.o amotry.o qromo.o midpnt.o polint.o qromo_alternate.o midpnt_alternate.o \
	polint_alternate.o spline.o splint.o \
	 -L/usr/local/lib \
	 -L/usr/lib \
	 -lm 
patch_integrate_besselk_10oct16_1_500.o : patch_integrate_besselk_10oct16_1_500.c
	gfortran -c patch_integrate_besselk_10oct16_1_500.c
nrutil.o : nrutil.c
	gfortran -c nrutil.c
amoeba.o : amoeba.c
	gfortran -c amoeba.c
amotry.o : amotry.c
	gfortran -c amotry.c
besselk.o : besselk.f
	gfortran -c besselk.f
adpint.o : adpint.f
	gfortran -c adpint.f
qromo.o : qromo.c
	gfortran -c qromo.c
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

