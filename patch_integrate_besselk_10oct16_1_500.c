/* This is a c code that ought to fit a spectrum to a set of five measurements 
   at 22, 45, 150, 408, 1420 and 23G
   The spectrum is a model containing a broken power law synthrotron 
   This version attempts to improve accuracy by scaling values
   This version changes from adaptive integration to qromo for writing out the spectrum 
   in the 2000-6000 MHz range only in steps of 10 MHz.

   This version fits using integration in gamma space rather than x=nu/nu_c

   The output file generated is in 2 columns with frequency and temperature.

   Output is read and processed for smoothness using fit_complex_spectrum_b.py

   Original from Mayuri 03 Dec 2015
   Modified by Ravi thereafter
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
//#include </usr/local/miriad/inc/maxdimc.h>
#include </home/mayuris/software/jnk/miriad/inc/linux64/maxdimc.h>
//#include </Users/mayuri/work/C/ANSI-C/nrutil.h>
#include </home/mayuris/work/HEALPy/codes/broken_pl_28dec15/nrutil.h>
//#include </Users/rsubrahm/CODESPACE/INCLUDE/nrutil.h>
//#include </usr/local/miriad/inc/miriad.h>
#include </home/mayuris/software/jnk/miriad/inc/miriad.h>
#include <string.h>

typedef double (*func1arg) (double);

void amoeba(float **p, float y[], int ndim, float ftol,
	    float (*funk)(float []), int *nfunk);
float amotry(float **p, float y[], float psum[], int ndim,
	     float (*funk)(float []), int ihi, float fac);

void mmbskr_(double *ARG,double *ORDER,int *N, double *BK, int *IER);

void adpint_(double *rint, double *xl, double *xu, double *reps, double *aeps,
	     double *dif, double (*func)(double), int *ier, int *npt, int *nmax);
double qromo(double (*func)(double), double a, double b,
	     double (*choose)(double(*)(double), double, double, int));
double midpnt(double (*func)(double), double a, double b, int n);
double qromo_alternate(double (*func)(double), double a, double b,
		       double (*choose)(double(*)(double), double, double, int));
double midpnt_alternate(double (*func)(double), double a, double b, int n);
void polint(double xa[], double ya[], int n, double x, double *y, double *dy);
void polint_alternate(double xa[], double ya[], int n, double x, double *y, double *dy);

/* functions for interpolation of recombination data 
   and conversion between brightness intensity and temperature */
void spline(float x[], float y[], int n, float yp1, float ypn, float y2[]);
void splint(float xa[], float ya[], float y2a[], int n, float x, float *y);

#define kk 1.3806488e-23 // m^2 kg s^-1
#define hh 6.62606957e-34 // m^2 kg s^-2 K^-1
#define cvel 2.99792458e+08 // m s^-1
#define m_e 9.1e-31
#define q_e 1.6e-19
#define sin_alph 1.0
#define Bmag 1e-9 // Tesla == 10 micro-Gauss
#define PI 3.14159265358979
#define FTOL 1.0e-6
#define T_e 8000.0 // Thermal electron temperature
#define TCMB 2.72548
#define NHPIX 3072
#define GSPAN 100.0

#define NN 5000 // length of array containing EoR spectrum

/* Global variables */
float *xx, *yy;
double xnu,xri,xrk,xrip,xrkp;
double nu,alpha1,alpha2,temp1,temp2,output,fnorm,fnorm1,fnorm2;
double AA,BB,C1,C2,D,ft,therm,fbreak;
int flag;
double nu_break;
float x,y;
float f1,f2;
float inttolf = 1.0e-8;
float difft[7];
double gama_min,gama_max,gama_break,scale_gam_nu;
int nrecurs = 10000000;

/* Evaulate bessel function of order 5/3, integrand  */
double modbessik2(double u)
{
  double ARG, ORDER, BK[1];
  int N, IER;
  xnu = 5.0/3.0;
  ORDER = (double) xnu;
  N = 1;
  ARG = (double)(1.0/u);
  mmbskr_(&ARG,&ORDER,&N,BK,&IER);
  xrk = BK[0]/exp(ARG);
  return (double)(xrk/(u*u));
}

/* Evaulate bessel function of order 5/3, integrand
   float modbessik2float(float u)
   {
   double ARG, ORDER, BK[1];
   int N, IER;
   xnu = 5.0/3.0;
   ORDER = (double) xnu;
   N = 1;
   ARG = (double)(1.0/u);
   mmbskr_(&ARG,&ORDER,&N,BK,&IER);
   xrk = BK[0]/exp(ARG);
   return (float)(xrk/(u*u));
   }
*/

double fofx1(double gama)
{
  double p1,nu_c;
  double integ;
    
  double x;
  double rint, xl, xu, reps, aeps, dif;
  int ier, npt, nmax;

  nu_c = (gama*gama*scale_gam_nu)/1.0e9; //GHz
  x = ((double)nu/nu_c);
    
  xl = (double)0.0;
  xu = (double)(1/x);
  reps = (double)inttolf;
  aeps = (double)inttolf;
  nmax = nrecurs;
    
  adpint_(&rint,&xl,&xu,&reps,&aeps,&dif,modbessik2,&ier,&npt,&nmax);
  //adapt_gq3(modbessik2,xl,xu,reps,&rint,&aeps,&npt);
  //  rint = (double)qromo(modbessik2float,(float)xl,(float)xu,midpnt);
    
  p1 = (double)((2*C1) - 3.0); 
  integ = rint*pow(gama,-1.0*p1)*(x);
  return integ;
}

double fofx1float(double gamafloat)
{
  double p1,nu_c,gama;
  double integ;

  double x;
  double rint, xl, xu, reps, aeps, dif;
  int ier, npt, nmax;
    
  gama = (double)gamafloat;
  nu_c = (gama*gama*scale_gam_nu)/1.0e9; //GHz
  x = ((double)nu/nu_c);
    
  xl = (double)0.0;
  xu = (double)(1/x);
  reps = (double)inttolf;
  aeps = (double)inttolf;
  nmax = nrecurs;
    
  //  adpint_(&rint,&xl,&xu,&reps,&aeps,&dif,modbessik2,&ier,&npt,&nmax);
  //  rint = (double)qromo_alternate(modbessik2float,(float)xl,(float)xu,midpnt_alternate);
  rint = qromo_alternate(modbessik2,xl,xu,midpnt_alternate);
    
  p1 = (double)((2*C1) - 3.0); 
  integ = rint*pow(gama,-1.0*p1)*(x);
  return integ;
}


float func(float pp[])
{
  int i;
  double FFIT,nu_min,nu_max,nu_break,Tx,Te,nu_t,extn,d_alph;
  float chisq,DDIF;

  double rint, rint1, rint2, xl, xu, xb, reps, aeps, dif;
  int ier, npt, nmax;

  /* Use global variables C1 for temperature spectral index to pass to fofx1 */
  fnorm = pow(10.0,(double)pp[1]);
  alpha1 = pow(10.0,(double)pp[2]);
  d_alph = pow(10.0,(double)pp[3]);
  alpha2 = alpha1+d_alph;
  nu_break = pow(10.0,(double)pp[4]);
  Tx = pow(10.0,(double)pp[5]);
  Te = pow(10.0,(double)pp[6]);
  nu_t = pow(10.0,(double)pp[7]);

  reps = (double)inttolf;
  aeps = (double)inttolf;
  nmax = nrecurs;
 
  /* Impose conditions on parameters, If violated chisquare explodes */
  if ( alpha1 < 2.0 || alpha1 > 3.0)
    { return 100000.0; }
  if ( alpha2 < 2.0 || alpha2 > 3.0)
    { return 100000.0; }
  if(Te < 0.0 || Te > 10000.0)
    { return 100000.0; }

  /* Compute Chisquare */
  chisq=0.0;
  for (i=1;i<=6;i++)
    {
      nu = (double)xx[i];		
      nu_min = nu*1e9/(double)GSPAN;
      nu_max = nu*1e9*(double)GSPAN;
      gama_min = sqrt((double)(nu_min)/scale_gam_nu);
      gama_max = sqrt((double)(nu_max)/scale_gam_nu);
      gama_break = sqrt((double)(nu_break)/scale_gam_nu);
      xl = (double)gama_min;
      xu = (double)gama_max;
      xb = (double)gama_break;
      
      if(xl > xb)
	{
	  C1 = alpha2;
	  //adapt_gq3(fofx1,xl,xu,reps,&rint,&aeps,&npt);
	  adpint_(&rint,&xl,&xu,&reps,&aeps,&dif,fofx1,&ier,&npt,&nmax);
	  rint *= pow((double)(gama_break),(double)(2*C1-3));
	}
      else if(xu < xb)
	{
	  C1 = alpha1;
	  //adapt_gq3(fofx1,xl,xu,reps,&rint,&aeps,&npt);
	  adpint_(&rint,&xl,&xu,&reps,&aeps,&dif,fofx1,&ier,&npt,&nmax);
	  rint *= pow((double)(gama_break),(double)(2*C1-3));
	}
      else
	{
	  xu = xb;
	  C1 = alpha1;
	  //adapt_gq3(fofx1,xl,xu,reps,&rint1,&aeps,&npt);
	  adpint_(&rint1,&xl,&xu,&reps,&aeps,&dif,fofx1,&ier,&npt,&nmax);
	  rint1 *= pow((double)(gama_break),(double)(2*C1-3));
	  xl = xb;
	  xu = (double)gama_max;
	  C1 = alpha2;
	  //adapt_gq3(fofx1,xl,xu,reps,&rint2,&aeps,&npt);
	  adpint_(&rint2,&xl,&xu,&reps,&aeps,&dif,fofx1,&ier,&npt,&nmax);
	  rint2 *= pow((double)(gama_break),(double)(2*C1-3));
	  rint = rint1 + rint2;
	}
      extn = exp((double)-1.0*pow((nu_t/nu),(double)2.1));
      FFIT = fnorm*((pow(nu,-2.0)*rint) + Tx*pow(nu,-2.1))*extn + Te*((double)1.0-extn);
      DDIF = (yy[i] - (float)FFIT)/yy[i];
      difft[i]=(DDIF);
      if (i <= 6) {chisq += (DDIF*DDIF);}
    }
  chisq /= 6.0;
  return chisq;
}

float func1(float pp[])
{

  int i;
  double FFIT,d_alph,temp,Tx,Te,extn,nu_t;
  float chisq,DDIF;

  /* Use global variables C1 for temperature spectral index to pass to fofx1 */  
  
  fnorm1 = pow(10.0,(double)pp[1]);
  alpha1 = pow(10.0,(double)pp[2]);
  d_alph = pow(10.0,(double)pp[3]);
  fnorm2 = pow(10.0,(double)pp[4]);
  alpha2 = alpha1-d_alph;
  Tx = pow(10.0,(double)pp[5]);
  Te = pow(10.0,(double)pp[6]);
  nu_t = pow(10.0,(double)pp[7]);

  /* Impose conditions on parameters, If violated chisquare explodes */
  if ( alpha1 < 2.0 || alpha1 > 3.0)
    { return 100000.0; }
  if ( alpha2 < 2.0 || alpha2 > 3.0)
    { return 100000.0; }  
  if(Te < 0.0 || Te > 10000.0)
    { return 100000.0; }

  /* Compute Chisquare */
  chisq=0.0;
  for (i=1;i<=6;i++)
    {
      extn = exp((double)-1.0*pow((nu_t/xx[i]),2.1));
      temp = pow((double)xx[i],-1.0*alpha1) + fnorm2*pow((double)xx[i],-1.0*alpha2);
      FFIT = fnorm1*(temp + Tx*pow(xx[i],-2.1))*extn + Te*((double)1.0-extn);
      DDIF = (yy[i] - (float)FFIT)/yy[i];
      difft[i]=DDIF;
      if (i <= 6) {chisq += (DDIF*DDIF);}
    }
  chisq /= 6.0;
  return chisq;
}

main()
{
  static const char filename_22[] = "../../data/all_sky_maps_5deg/map_22_r4_5deg_nested_galactic_Kelvin.txt";
  static const char filename_45[] = "../../data/all_sky_maps_5deg/map_45_r4_5deg_nested_galactic_Kelvin.txt";
  static const char filename_150[] = "../../data/all_sky_maps_5deg/map_150_r4_5deg_nested_galactic_Kelvin.txt";
  static const char filename_408[] = "../../data/all_sky_maps_5deg/map_408_r4_5deg_nested_galactic_Kelvin.txt";
  static const char filename_1420[] = "../../data/all_sky_maps_5deg/map_1420_r4_5deg_nested_galactic_Kelvin.txt";
  static const char filename_23G[] = "../../data/all_sky_maps_5deg/map_23000_r4_5deg_nested_galactic_Kelvin.txt";  
  static const char filename_fit_params[]="../../data/eor_sim/total_model_params_10oct16_1_500.txt";
  static const char filename_pixel_spec[]="../../data/eor_sim/total_model_spec_10oct16_1_500.txt";
  static const char filename_eor[]= "../../data/eor_sim/eor_Tb_values.dat";

  double rint, rint1, rint2, xl, xu, xb, reps, aeps, dif,Tx,Te,extn, nu_t;
  int ier, npt, nmax;

  int npoints1;
  float *xarray,*yarray,*y2;
  float read1,read2;
  float ffmin,ffmax;
  float yp1,ypn,yf;
    
  /* Variable declarations */
  float *yamoeba,*afit,*aafit,*cof,**pfun,*yfit,*sky_22,
    *sky_45,*sky_150,*sky_408,*sky_1420,*sky_23000;
  FILE *fptr, *fptr1,*fptr2;
  float read_value;
  float break_freq;
  float correction_150offset = 21.4; // Kelvin, Patra et al. 2015
  float correction_150scaling = 1.05; // Patra et al. 2015
  char inbuf[200];
  float T_150,T_408,T_1420,T_23G,gama_1,gama_2,gama_3,nu_min,nu_max;
  float chisq;
  int i,ii,jj,j,nu_min_MHz,nu_max_MHz,nfunc;
  xx = vector(1,6);
  yy = vector(1,6);
  yamoeba = vector(1,8);
  afit = vector(1,7);
  aafit = vector(1,7);
  cof = vector(1,7);
  pfun = matrix(1,8,1,7);
  sky_22 = vector(1,NHPIX);
  sky_45 = vector(1,NHPIX);
  sky_150 = vector(1,NHPIX);
  sky_408 = vector(1,NHPIX);
  sky_1420 = vector(1,NHPIX);
  sky_23000 = vector(1,NHPIX);

  // read in the sky maps and coordinates; subtract CMB where necessary;

  fptr=fopen(filename_22,"r");
  i=1;
  while (fgets (inbuf, sizeof inbuf, fptr) != NULL)
    {
      sscanf(inbuf,"%g",&read_value); 
      sky_22[i]=read_value;
      ++i;
    }
  fclose(fptr);
  //printf(" Read %d values from 22 MHz image \n",i);

  fptr=fopen(filename_45,"r");
  i=1;
  while (fgets (inbuf, sizeof inbuf, fptr) != NULL)
    {
      sscanf(inbuf,"%g",&read_value); 
      sky_45[i]=read_value;
      ++i;
    }
  fclose(fptr);
  //printf(" Read %d values from 45 MHz image \n",i);

  fptr=fopen(filename_150,"r");
  i=1;
  while (fgets (inbuf, sizeof inbuf, fptr) != NULL)
    {
      sscanf(inbuf,"%g",&read_value); 
      sky_150[i]=(read_value - correction_150offset)*correction_150scaling; 
      sky_150[i]= sky_150[i]-TCMB;
      ++i;
    }
  fclose(fptr);
  //printf(" Read %d values from 150 MHz image \n",i);
		
  fptr=fopen(filename_408,"r");
  i=1;
  while (fgets (inbuf, sizeof inbuf, fptr) != NULL)
    {
      sscanf(inbuf,"%g",&read_value);
      sky_408[i]=read_value-TCMB;
      ++i;
    }
  fclose(fptr);
  //printf(" Read %d values from 408 MHz image \n",i);
		
  fptr=fopen(filename_1420,"r");
  i=1;
  while (fgets (inbuf, sizeof inbuf, fptr) != NULL)
    {
      sscanf(inbuf,"%g",&read_value);
      sky_1420[i]=read_value-TCMB;
      ++i;   
    }
  fclose(fptr);
  //printf(" Read %d values from 1420 MHz image \n",i);

  fptr=fopen(filename_23G,"r");
  i=1;
  while (fgets (inbuf, sizeof inbuf, fptr) != NULL)
    {
      sscanf(inbuf,"%g",&read_value);
      sky_23000[i]=read_value;
      ++i;
    }
  fclose(fptr);
  //printf(" Read %d values from 23000 MHz image \n",i);

  int fflush(FILE *fptr); 

  fptr  = fopen(filename_fit_params,"a+");  
  fptr2  = fopen(filename_pixel_spec,"a+");  
  scale_gam_nu = (double)(3.0*q_e*Bmag*sin_alph)/(4.0*PI*m_e*cvel);

  for(ii=60; ii<76; ii++)
    {
      nu = powf(1.05,(float)ii)/1e3;
      fprintf(fptr2,"%f ",nu);
    }
  for (ii=40; ii<201; ii++)
    {
      nu = (float)(ii)/1e3;
      fprintf(fptr2,"%f ",nu);
    }
  for(ii=109; ii<210; ii++)
    {
      nu = powf(1.05,(float)ii)/1e3;
      fprintf(fptr2,"%f ",nu);
    }
  fprintf(fptr2," \n");


  /* Frequency in GHz */
  xx[1]=0.022;
  xx[2]=0.045;
  xx[3]=0.150;
  xx[4]=0.408;
  xx[5]=1.420;
  xx[6]=22.690;

  reps = (double)inttolf;
  aeps = (double)inttolf;
  nmax = nrecurs;
  //  printf("Pixel    alpha1     alpha2\n");
  for (j=1;j<=500;++j)
    {
      printf("\n\n****************\n\nPixel number: %d\n",j);             
	{
	  yy[1] = sky_22[j];
	  yy[2] = sky_45[j];
	  yy[3] = sky_150[j];
	  yy[4] = sky_408[j];
	  yy[5] = sky_1420[j];
	  yy[6] = sky_23000[j];

	  printf("temperatures : %f %f %f %f %f %f \n",yy[1],yy[2],yy[3],yy[4],yy[5],yy[6]);

	  /* Spectral index between 45 MHz and 150 MHz  */
	  alpha1 = (double)(log10f(yy[2]) - log10f(yy[3]))/(log10f(xx[3]) - log10f(xx[2]));
	  if(alpha1 < 2.0) alpha1 = 2.0001;
	  if(alpha1 > 3.0) alpha1 = 2.9999;
	  /* Spectral index between 408 MHz and 1420 MHz  */
	  alpha2 = (double)(log10f(yy[4]) - log10f(yy[5]))/(log10f(xx[5]) - log10f(xx[4]));
	  if(alpha2 < 2.0) alpha2 = 2.0001;
	  if(alpha2 > 3.0) alpha2 = 2.9999;
      
	  //	  printf("%d     %f      %f\n",j,alpha1,alpha2);
	  nu_break = sqrt(0.150*0.408)*1e9;
	  Te = 8000.0;
	  nu_t = 0.001;
	  extn = (double)exp((double)-1.0*pow((nu_t/xx[5]),(double)2.1));

	  /*  If data requires alpha2 steeper than alpha1, then model as synchrotron with break */
	  if(alpha1 < alpha2)
	    {
	      printf("Alpha2 steeper than alpha1! \n");
	      /* Initial computation of normalization */
	      nu = (double)1.420;
	      flag = 0;
	      nu_min = nu*1e9/GSPAN;
	      nu_max = nu*1e9*GSPAN;
	      gama_min = sqrt((double)(nu_min)/scale_gam_nu);
	      gama_max = sqrt((double)(nu_max)/scale_gam_nu);
	      gama_break = sqrt((double)(nu_break)/scale_gam_nu);
      
	      xb = gama_break;
	      xl = gama_min;
	      xu = gama_max;
      
	      if(xl > xb)
		{
		  C1 = alpha2;
		  //	      adapt_gq3(fofx1,xl,xu,reps,&rint,&aeps,&npt);
		  adpint_(&rint,&xl,&xu,&reps,&aeps,&dif,fofx1,&ier,&npt,&nmax);
		  rint *= pow((double)(gama_break),(double)(2*C1-3));
		}
	      else if(xu < xb)
		{
		  C1 = alpha1;
		  //adapt_gq3(fofx1,xl,xu,reps,&rint,&aeps,&npt);
		  adpint_(&rint,&xl,&xu,&reps,&aeps,&dif,fofx1,&ier,&npt,&nmax);
		  rint *= pow((double)(gama_break),(double)(2*C1-3));
		}
	      else
		{
		  xu = xb;
		  C1 = alpha1;
		  //adapt_gq3(fofx1,xl,xu,reps,&rint1,&aeps,&npt);
		  adpint_(&rint1,&xl,&xu,&reps,&aeps,&dif,fofx1,&ier,&npt,&nmax);
		  rint1 *= pow((double)(gama_break),(double)(2*C1-3));
		  xl = xb;
		  xu = (double)gama_max;
		  C1 = alpha2;
		  //adapt_gq3(fofx1,xl,xu,reps,&rint2,&aeps,&npt);
		  adpint_(&rint2,&xl,&xu,&reps,&aeps,&dif,fofx1,&ier,&npt,&nmax);
		  rint2 *= pow((double)(gama_break),(double)(2*C1-3));
		  rint = rint1 + rint2;
		}
	      fnorm = (double)(yy[5] - (Te*((double)1.0-(double)extn)))/((powf(nu,-2.0)*(float)rint)*extn);
	      printf("fnorm is %f\n",fnorm);
	      printf("rint is %f\n",rint);
	      printf("extn is %f\n",extn);

	      nu = (double)22.690;
	      nu_min = nu*1e9/GSPAN;
	      nu_max = nu*1e9*GSPAN;
	      gama_min = sqrt((double)(nu_min)/scale_gam_nu);
	      gama_max = sqrt((double)(nu_max)/scale_gam_nu);
	      gama_break = sqrt((double)(nu_break)/scale_gam_nu);
      
	      xb = gama_break;
	      xl = gama_min;
	      xu = gama_max;
      
	      if(xl > xb)
		{
		  C1 = alpha2;
		  adpint_(&rint,&xl,&xu,&reps,&aeps,&dif,fofx1,&ier,&npt,&nmax);
		  rint *= pow((double)(gama_break),(double)(2*C1-3));
		}
	      else if(xu < xb)
		{
		  C1 = alpha1;
		  adpint_(&rint,&xl,&xu,&reps,&aeps,&dif,fofx1,&ier,&npt,&nmax);
		  rint *= pow((double)(gama_break),(double)(2*C1-3));
		}
	      else
		{
		  xu = xb;
		  C1 = alpha1;
		  adpint_(&rint1,&xl,&xu,&reps,&aeps,&dif,fofx1,&ier,&npt,&nmax);
		  rint1 *= pow((double)(gama_break),(double)(2*C1-3));
		  xl = xb;
		  xu = (double)gama_max;
		  C1 = alpha2;
		  adpint_(&rint2,&xl,&xu,&reps,&aeps,&dif,fofx1,&ier,&npt,&nmax);
		  rint2 *= pow((double)(gama_break),(double)(2*C1-3));
		  rint = rint1 + rint2;
		}
	      extn = exp((double)-1.0*pow((nu_t/nu),2.1));
	      temp1 = fnorm*extn;
	      Tx = (((yy[6] - Te*(1.0-extn))/temp1)-(powf(nu,-2.0)*(float)rint))/pow(xx[5],-2.1);
	      if (Tx <= 0.0)
		{
		  Tx = (double)1.0e-10;
		  printf("Tx negative!\n");
		}
	      printf("Te is %f\n",Te);
	      printf("Initial values of fnorm, alpha1, alpha2 fbreak (MHz), Tx, Te, nu_t: %e %f %f %f %e %f %f\n",fnorm,alpha1,alpha2,nu_break/1e6,Tx,Te,nu_t);
      
	      afit[1] = (float)log10(fnorm);
	      afit[2] = (float)log10(alpha1);
	      afit[3] = (float)log10(alpha2 - alpha1);
	      afit[4] = (float)log10(nu_break);
	      afit[5] = (float)log10(Tx);
	      afit[6] = (float)log10(Te);
	      afit[7] = (float)log10(nu_t);
	      chisq = func(afit);
	      printf("Chisq and sqrt(chisq) before optimization = %e %e \n",chisq,sqrtf(chisq));
	      printf(" differences: %f %f %f %f %f %f\n",
		     difft[1],difft[2],difft[3],difft[4],difft[5],difft[6]);

	      for(ii=1;ii<=8;ii++)
		{
		  for(jj=1;jj<=7;jj++)
		    {
		      aafit[jj]=pfun[ii][jj]=(ii==(jj+1) ? afit[jj]*1.5 : afit[jj]);
		    }
		  yamoeba[ii]=func(aafit);
		}
	      printf("Going to optimize\n");
	      amoeba(pfun,yamoeba,7,FTOL,func,&nfunc);

	      printf("\n OPTIMIZATION DONE\n\n");

	      cof[1]=pfun[1][1];
	      cof[2]=pfun[1][2];
	      cof[3]=pfun[1][3];
	      cof[4]=pfun[1][4];
	      cof[5]=pfun[1][5];
	      cof[6]=pfun[1][6];
	      cof[7]=pfun[1][7];

	      afit[1] = cof[1];
	      afit[2] = cof[2];
	      afit[3] = cof[3];
	      afit[4] = cof[4];
	      afit[5] = cof[5];
	      afit[6] = cof[6];
	      afit[7] = cof[7];

	      fnorm = powf(10.0,cof[1]);
	      printf("Fit variables fnorm, alpha1, alpha2, nu_break (MHz), Tx, Te, nu_t: %e %f %f %f %f %f %f\n",
		     fnorm,powf(10.0,cof[2]),powf(10.0,cof[3])+powf(10.0,cof[2]),powf(10.0,cof[4])/1e6,powf(10.0,cof[5]),powf(10.0,cof[6]),powf(10.0,cof[7]));
	
	      chisq = func(afit);
	      printf("Chisq and sqrt(chisq) after optimization = %e %e \n",chisq,sqrtf(chisq));
	      printf(" differences: %f %f %f %f %f %f \n",difft[1],difft[2],difft[3],difft[4],difft[5],difft[6]);
	      fprintf(fptr,"%d %d %f %f %f %f %f %f %f %f\n",j,flag,cof[1],cof[2],cof[3],cof[4],cof[5],cof[6],cof[7],chisq);
	 
	      /* for(ii=60; ii<76; ii++) */
	      /* 	{ */
	      /* 	  nu = pow((double)1.05,(double)ii)/(double)1e3; */
	      /* 	  nu_min = nu*1e9/GSPAN; */
	      /* 	  nu_max = nu*1e9*GSPAN; */
	      /* 	  gama_min = sqrt((nu_min)/scale_gam_nu); */
	      /* 	  gama_max = sqrt((nu_max)/scale_gam_nu); */
	      /* 	  xl = (double)gama_min; */
	      /* 	  xu = (double)gama_max; */
	    
	      /* 	  fnorm = pow((double)10.0,(double)afit[1]); */
	      /* 	  alpha1 = pow(10.0,(double)afit[2]); */
	      /* 	  alpha2 = pow(10.0,(double)afit[3])+pow(10.0,(double)afit[2]); */
	      /* 	  nu_break = pow(10.0,(double)afit[4]); */
	      /* 	  gama_break = sqrt((nu_break)/scale_gam_nu); */
	      /* 	  xb = (double)gama_break; */
	      /* 	  Tx = pow((double)10.0,(double)afit[5]); */
	      /* 	  Te = pow((double)10.0,(double)afit[6]); */
	      /* 	  nu_t = pow((double)10.0,(double)afit[7]); */

	      /* 	  if(xl > xb) */
	      /* 	    { */
	      /* 	      C1 = alpha2; */
	      /* 	      adpint_(&rint,&xl,&xu,&reps,&aeps,&dif,fofx1,&ier,&npt,&nmax); */
	      /* 	      rint *= pow((double)(gama_break),(double)(2*C1-3)); */
	      /* 	    } */
	      /* 	  else if(xu < xb) */
	      /* 	    { */
	      /* 	      C1 = alpha1; */
	      /* 	      adpint_(&rint,&xl,&xu,&reps,&aeps,&dif,fofx1,&ier,&npt,&nmax); */
	      /* 	      rint *= pow((double)(gama_break),(double)(2*C1-3)); */
	      /* 	    } */
	      /* 	  else */
	      /* 	    { */
	      /* 	      xu = xb; */
	      /* 	      C1 = alpha1; */
	      /* 	      adpint_(&rint1,&xl,&xu,&reps,&aeps,&dif,fofx1,&ier,&npt,&nmax); */
	      /* 	      rint1 *= pow((double)(gama_break),(double)(2*C1-3)); */
	      /* 	      xl = xb; */
	      /* 	      xu = (double)gama_max; */
	      /* 	      C1 = alpha2; */
	      /* 	      adpint_(&rint2,&xl,&xu,&reps,&aeps,&dif,fofx1,&ier,&npt,&nmax); */
	      /* 	      rint2 *= pow((double)(gama_break),(double)(2*C1-3)); */
	      /* 	      rint = rint1 + rint2; */
	      /* 	    } */
	      /* 	  extn = exp((double)-1.0*pow((nu_t/nu),2.1)); */
	      /* 	  output = fnorm*(powf(nu,-2.0)*(float)rint + Tx*pow(nu,-2.1))*extn + Te*(1.0-extn); */
	      /* 	  fprintf(fptr2,"%20.15lf ",output); */
	      /* 	} */
	      for(ii=2000; ii<6001; ii=ii+10)
		{
		  //printf(" Computing output spectrum at freq %d MHz \n",ii);
		  nu = (double)(ii)/(double)1e3;
		  nu_min = nu*(double)1e9/(double)GSPAN;
		  nu_max = nu*(double)1e9*(double)GSPAN;
		  gama_min = sqrt((nu_min)/scale_gam_nu);
		  gama_max = sqrt((nu_max)/scale_gam_nu);
		  xl = (double)gama_min;
		  xu = (double)gama_max;
	  
		  fnorm = pow((double)10.0,(double)afit[1]);
		  alpha1 = pow(10.0,(double)afit[2]);
		  alpha2 = pow(10.0,(double)afit[3])+pow(10.0,(double)afit[2]);
		  nu_break = pow(10.0,(double)afit[4]);
		  gama_break = sqrt((nu_break)/scale_gam_nu);
		  xb = (double)gama_break;
		  Tx = pow((double)10.0,(double)afit[5]);
		  Te = pow((double)10.0,(double)afit[6]);
		  nu_t = pow((double)10.0,(double)afit[7]);

		  if(xl > xb)
		    {
		      C1 = alpha2;
		      //				adpint_(&rint,&xl,&xu,&reps,&aeps,&dif,fofx1,&ier,&npt,&nmax);
		      rint = (double)qromo(fofx1float,(float)xl,(float)xu,midpnt);
		      rint *= pow((double)(gama_break),(double)(2*C1-3));
		    }
		  else if(xu < xb)
		    {
		      C1 = alpha1;
		      //				adpint_(&rint,&xl,&xu,&reps,&aeps,&dif,fofx1,&ier,&npt,&nmax);
		      rint = (double)qromo(fofx1float,(float)xl,(float)xu,midpnt);
		      rint *= pow((double)(gama_break),(double)(2*C1-3));
		    }
		  else
		    {
		      xu = xb;
		      C1 = alpha1;
		      //				adpint_(&rint1,&xl,&xu,&reps,&aeps,&dif,fofx1,&ier,&npt,&nmax);
		      rint1 = (double)qromo(fofx1float,(float)xl,(float)xu,midpnt);
		      rint1 *= pow((double)(gama_break),(double)(2*C1-3));
		      xl = xb;
		      xu = (double)gama_max;
		      C1 = alpha2;
		      //				adpint_(&rint2,&xl,&xu,&reps,&aeps,&dif,fofx1,&ier,&npt,&nmax);
		      rint2 = (double)qromo(fofx1float,(float)xl,(float)xu,midpnt);
		      rint2 *= pow((double)(gama_break),(double)(2*C1-3));
		      rint = rint1 + rint2;
		    }
		  extn = exp((double)-1.0*pow((nu_t/nu),2.1));
		  output = fnorm*(powf(nu,-2.0)*(float)rint + Tx*pow(nu,-2.1))*extn + Te*(1.0-extn);/* + (double)yf */
		  fprintf(fptr2,"%20.15lf ",output);
		}
	      /* for(ii=109; ii<210; ii++) */
	      /* 	{ */
	      /* 	  nu = pow((double)1.05,(double)ii)/(double)1e3; */
	      /* 	  nu_min = nu*1e9/GSPAN; */
	      /* 	  nu_max = nu*1e9*GSPAN; */
	      /* 	  gama_min = sqrt((nu_min)/scale_gam_nu); */
	      /* 	  gama_max = sqrt((nu_max)/scale_gam_nu); */
	      /* 	  xl = (double)gama_min; */
	      /* 	  xu = (double)gama_max; */

	      /* 	  fnorm = pow((double)10.0,(double)afit[1]); */
	      /* 	  alpha1 = pow(10.0,(double)afit[2]); */
	      /* 	  alpha2 = pow(10.0,(double)afit[3])+pow(10.0,(double)afit[2]); */
	      /* 	  nu_break = pow(10.0,(double)afit[4]); */
	      /* 	  gama_break = sqrt((nu_break)/scale_gam_nu); */
	      /* 	  xb = (double)gama_break; */
	      /* 	  Tx = pow((double)10.0,(double)afit[5]); */
	      /* 	  Te = pow((double)10.0,(double)afit[6]); */
	      /* 	  nu_t = pow((double)10.0,(double)afit[7]); */


	      /* 	  if(xl > xb) */
	      /* 	    { */
	      /* 	      C1 = alpha2; */
	      /* 	      adpint_(&rint,&xl,&xu,&reps,&aeps,&dif,fofx1,&ier,&npt,&nmax); */
	      /* 	      rint *= pow((double)(gama_break),(double)(2*C1-3)); */
	      /* 	    } */
	      /* 	  else if(xu < xb) */
	      /* 	    { */
	      /* 	      C1 = alpha1; */
	      /* 	      adpint_(&rint,&xl,&xu,&reps,&aeps,&dif,fofx1,&ier,&npt,&nmax); */
	      /* 	      rint *= pow((double)(gama_break),(double)(2*C1-3)); */
	      /* 	    } */
	      /* 	  else */
	      /* 	    { */
	      /* 	      xu = xb; */
	      /* 	      C1 = alpha1; */
	      /* 	      adpint_(&rint1,&xl,&xu,&reps,&aeps,&dif,fofx1,&ier,&npt,&nmax); */
	      /* 	      rint1 *= pow((double)(gama_break),(double)(2*C1-3)); */
	      /* 	      xl = xb; */
	      /* 	      xu = (double)gama_max; */
	      /* 	      C1 = alpha2; */
	      /* 	      adpint_(&rint2,&xl,&xu,&reps,&aeps,&dif,fofx1,&ier,&npt,&nmax); */
	      /* 	      rint2 *= pow((double)(gama_break),(double)(2*C1-3)); */
	      /* 	      rint = rint1 + rint2; */
	      /* 	    } */
	      /* 	  extn = exp((double)-1.0*pow((nu_t/nu),2.1)); */
	      /* 	  output = fnorm*(powf(nu,-2.0)*(float)rint + Tx*pow(nu,-2.1))*extn + Te*(1.0-extn); */
	      /* 	  fprintf(fptr2,"%20.15lf ",output); */
	      /* 	} */
	      fprintf(fptr2,"\n");
	      /********************************************/
	    }
	  /*  If data requires alpha2 flatter than alpha1, then model as sum of power laws, i.e. steep and flat spectrum sources */
	  else
	    {
	      printf("alpha2 flatter than alpha1!\n");
	      flag = 1;
	      /* Initial computation of normalization of the form: T(nu) = fnorm1*(nu^-alpha1) + fnorm2*(nu^-alpha2) */
	      fnorm1 = yy[3]/(pow(xx[3],-1.0*alpha1));
	      fnorm2 = (yy[4]/(pow(xx[4],-1.0*alpha2)))/fnorm1;	     
	  
	      extn = exp((double)-1.0*pow((nu_t/xx[6]),2.1));
	      Tx = ((double)1.0/pow(xx[6],-2.1))*(((yy[6]-Te*((double)1.0-extn))/(fnorm1*extn))-((pow(xx[6],-1.0*alpha1))+(fnorm2*pow(xx[6],-1.0*alpha2))));
	  
	      if (Tx <=0.0)
		{
		  printf("Tx negative!\n");
		  Tx = (double)1.0e-10;
		}
	  
	      afit[1] = (float)log10(fnorm1);
	      afit[2] = (float)log10(alpha1);
	      afit[3] = (float)log10(alpha1 - alpha2);
	      afit[4] = (float)log10(fnorm2);
	      afit[5] = (float)log10(Tx);
	      afit[6] = (float)log10(Te);
	      afit[7] = (float)log10(nu_t);

	      printf("Initial values of fnorm1, alpha1, alpha1-alpha2 fnorm2: %e %f %f %f %f %f %f\n",fnorm1,alpha1,alpha2,fnorm2,Tx,Te,nu_t);
	      chisq = func1(afit);
	      printf("Chisq and sqrt(chisq) before optimization = %e %e \n",chisq,sqrtf(chisq));
	      printf(" differences: %f %f %f %f %f %f\n",
		     difft[1],difft[2],difft[3],difft[4],difft[5],difft[6]);

	      for(ii=1;ii<=8;ii++)
		{
		  for(jj=1;jj<=7;jj++)
		    {
		      aafit[jj]=pfun[ii][jj]=(ii==(jj+1) ? afit[jj]*1.5 : afit[jj]);
		    }
		  yamoeba[ii]=func1(aafit);
		}

	      amoeba(pfun,yamoeba,7,FTOL,func1,&nfunc);

	      printf("\n OPTIMIZATION DONE\n\n");
	  
	      cof[1]=pfun[1][1];
	      cof[2]=pfun[1][2];
	      cof[3]=pfun[1][3];
	      cof[4]=pfun[1][4];
	      cof[5]=pfun[1][5];
	      cof[6]=pfun[1][6];
	      cof[7]=pfun[1][7];

	      afit[1] = cof[1];
	      afit[2] = cof[2];
	      afit[3] = cof[3];
	      afit[4] = cof[4];
	      afit[5] = cof[5];
	      afit[6] = cof[6];
	      afit[7] = cof[7];

	      fnorm1 = pow(10.0,(double)cof[1]);
	      fnorm2 = pow(10.0,(double)cof[4]);
	      alpha1 = pow(10.0,(double)cof[2]);
	      alpha2 = alpha1 - pow(10.0,(double)cof[3]);
	      printf("Fit variables fnorm1, alpha1, alpha2, fnorm2: %e %f %f %f\n",
		     fnorm1,powf(10.0,cof[2]),powf(10.0,cof[2])-powf(10.0,cof[3]),fnorm2);
	      Tx = pow((double)10.0,(double)afit[5]);
	      Te = pow((double)10.0,(double)afit[6]);
	      nu_t = pow((double)10.0,(double)afit[7]);
	
	      chisq = func1(afit);
	      printf("Chisq and sqrt(chisq) after optimization = %e %e \n",chisq,sqrtf(chisq));
	      printf(" differences: %f %f %f %f %f %f \n",difft[0],difft[1],difft[2],difft[3],difft[4],difft[5]);
	      fprintf(fptr,"%d %d %f %f %f %f %f %f %f %f\n",j,flag,fnorm1,cof[2],cof[3],fnorm2,cof[5],cof[6],cof[7],chisq);
	  	  
	      /* for(ii=60; ii<76; ii++) */
	      /* 	{ */
	      /* 	  nu = powf(1.05,(float)ii)/1e3; */
	      /* 	  extn = exp((double)-1.0*pow((nu_t/nu),(double)2.1)); */
	      /* 	  output = fnorm1*(pow(nu,-1.0*alpha1) + fnorm2*(pow(nu,-1.0*alpha2)) + Tx*pow(nu,-2.1))*extn + Te*((double)1.0 - extn); */
	      /* 	  fprintf(fptr2,"%20.15lf ",output); */
	      /* 	} */
	      for(ii=2000; ii<6001; ii=ii+10)
		{
		  nu = (float)(ii)/1e3;
		  extn = exp((double)-1.0*pow((nu_t/nu),(double)2.1));
		  output = fnorm1*(pow(nu,-1.0*alpha1) + fnorm2*(pow(nu,-1.0*alpha2)) + Tx*pow(nu,-2.1))*extn + Te*((double)1.0 - extn);
		  fprintf(fptr2,"%20.15lf ",output);
		}
	      /* for(ii=109; ii<210; ii++) */
	      /* 	{ */
	      /* 	  nu = powf(1.05,(float)ii)/1e3; */
	      /* 	  extn = exp((double)-1.0*pow((nu_t/nu),(double)2.1)); */
	      /* 	  output = fnorm1*(pow(nu,-1.0*alpha1) + fnorm2*(pow(nu,-1.0*alpha2)) + Tx*pow(nu,-2.1))*extn + Te*((double)1.0 - extn); */
	      /* 	  fprintf(fptr2,"%20.15lf ",output); */
	      /* 	} */
	  
	      fprintf(fptr2," \n");

	    }
	}
    }
  fclose(fptr);
  fclose(fptr2);
  free_vector(sky_22,1,NHPIX);
  free_vector(sky_45,1,NHPIX);
  free_vector(sky_150,1,NHPIX);
  free_vector(sky_408,1,NHPIX);
  free_vector(sky_1420,1,NHPIX);
  free_vector(sky_23000,1,NHPIX);
  free_vector(xx,1,6);
  free_vector(yy,1,6);
  free_vector(yamoeba,1,7);
  free_vector(afit,1,7);
  free_vector(aafit,1,7);
  free_matrix(pfun,1,8,1,7);
  free_vector(cof,1,7);
}

