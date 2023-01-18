import numpy as np
import scipy as sp
import matplotlib.pyplot as plt

def modbessik2(u):
	xnu = 5.0/3.0
	N = 1
	ARG = (1.0/u)
	ORDER = xnu
	mmbskr_(ARG,ORDER,N,BK,IER)
	xrk = BK[0]/np.exp(ARG)
    return xrk/(u*u)

def fofx1(gama):
	nu_c = (gama*gama*scale_gam_nu)/1.0e9
	x1 = 0.0
	xu = 1/x
	reps = float(inttolf)
	aeps = float(inttolf)
	nmax = nrecurs
    #adpint_(&rint,&xl,&xu,&reps,&aeps,&dif,modbessik2,&ier,&npt,&nmax)
    #rint = (double)qromo_alternate(modbessik2float,(float)xl,(float)xu,midpnt_alternate)
    rint = qromo_alternate(modbessik2,xl,xu,midpnt_alternate)
    p1 = (2*C1) - 3.0
    integ = rint*np.power(gama,-1.0*p1)*(x)
    return integ
