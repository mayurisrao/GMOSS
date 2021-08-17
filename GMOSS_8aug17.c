/* This file generates GMOSS spectra for every pixel out of 1-3072 in a HEALPix scheme
   CHANGE: PATH TO CURRENT DIRECTORY FROM HOME FOR nrutil.h

   Input file: 'frequecies_list.txt' should contain frequenceis in MHz in single column

   CHANGE: NFREQ to the number of frequencies in the list

   The output file generated is in 2 columns with frequency and temperature.
   NOTE: Before reading the output file with another code, remove any double/multiple
   spaces between the temperatures/frequencies for convenient columnwise reading.

   Original from Mayuri 03 Dec 2015
   Modified by Mayuri on 8 August 2017  
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include </home/mayuris/work/HEALPy/codes/broken_pl_28dec15/to_magendran/nrutil.h>
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
#define NFREQ 72 //NOTE CHANGE ACCORDING TO REQUIREMENT!!
#define GSPAN 100.0

#define NN 5000 // length of array containing EoR spectrum

/* Global variables */
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
    
  rint = qromo_alternate(modbessik2,xl,xu,midpnt_alternate);
    
  p1 = (double)((2*C1) - 3.0); 
  integ = rint*pow(gama,-1.0*p1)*(x);
  return integ;
}

main()
{  
  static const char filename_freq[] = "frequencies_list.txt";
  static const char filename_fit_params[]="GMOSS_params_allpix.txt";
  static const char filename_pixel_spec[]="output_spectra.txt";

  double rint, rint1, rint2, xl, xu, xb, reps, aeps, dif,Tx,Te,extn, nu_t;
  int ier, npt, nmax;
  double read_value;
  float read_value1,read_value2,read_value3,read_value4,read_value5,read_value6,read_value7,read_value8,read_value9,read_value10,read_value11,read_value12;
  int npoints1;
  float read1,read2;
  float ffmin,ffmax;
  float yp1,ypn,yf;
    
  /* Variable declarations */
  float *norm_param,*alpha1_param,*dalph_param,*nubrk_param,*Tx_param,*Te_param,*nut_param;
  double *freq_GHz;
  int *flag;
  FILE *fptr, *fptr1,*fptr2;
  float break_freq;
  char inbuf[200];
  float gama_1,gama_2,gama_3,nu_min,nu_max;
  float chisq;
  int i,ii,jj,j,nu_min_MHz,nu_max_MHz,nfunc;
  flag = ivector(1,NHPIX);
  norm_param = vector(1,NHPIX);
  alpha1_param = vector(1,NHPIX);
  dalph_param = vector(1,NHPIX);
  nubrk_param = vector(1,NHPIX);
  Tx_param = vector(1,NHPIX);
  Te_param = vector(1,NHPIX);
  nut_param = vector(1,NHPIX); 
  freq_GHz = dvector(1,NFREQ);
  
  // read in the parameters and frequencies;
  fptr=fopen(filename_freq,"r");
  printf("Opened file to read frequencies \n");
  i=1;
  while (fgets (inbuf, sizeof inbuf, fptr) != NULL)
    {
      sscanf(inbuf,"%lf",&read_value);
      freq_GHz[i]=(double)read_value/1e3;     
      ++i;         
    }
  printf("\n");
  fclose(fptr);

  fptr=fopen(filename_fit_params,"r");
  printf("Opened file to read parameters \n");
  i=1;
  while (fgets (inbuf, sizeof inbuf, fptr) != NULL)
    {
      sscanf(inbuf,"%f%f%f%f%f%f%f%f%f%f",&read_value1,&read_value2,&read_value3,&read_value4,&read_value5,&read_value6,&read_value7,&read_value8,&read_value9,&read_value10);
      flag[i]=(int)read_value2;
      norm_param[i] = read_value3;
      alpha1_param[i] = read_value4;
      dalph_param[i] = read_value5;
      nubrk_param[i] = read_value6;
      Tx_param[i] = read_value7;
      Te_param[i] = read_value8;
      nut_param[i] = read_value9;
      if(i == 10) printf("Printing parameters of pixel 10 as a consistency check! \n %f  %f  %f  %f  %f  %f  %f  %f  %f  %f\n\n",read_value1,read_value2,read_value3,read_value4,read_value5,read_value6,read_value7,read_value8,read_value9,read_value10);
      ++i;
    }
  fclose(fptr);
  int fflush(FILE *fptr); 
  int fflush(FILE *fptr);

  fptr2  = fopen(filename_pixel_spec,"a+");  
  printf("Opened file to write out spectra \n");
  scale_gam_nu = (double)(3.0*q_e*Bmag*sin_alph)/(4.0*PI*m_e*cvel);

  for(ii=1; ii<=NFREQ; ii++)
    {
      
      nu = (double)freq_GHz[ii];

      fprintf(fptr2,"%20.15lf ",nu);
    }
  fprintf(fptr2," \n");

  for (j=1;j<=3072;++j)
    {
      printf("\n\n****************\n\nPixel number: %d\n",j);             
	  if(flag[j]==0)
	    {
	      printf("alpha2 steeper than alpha1!\n");
	      for(ii=1; ii<=NFREQ; ii=ii+1)
		{
		  //printf(" Computing output spectrum at freq %d MHz \n",ii);
		  
		  nu = (double)freq_GHz[ii];
		  nu_min = nu*(double)1e9/(double)GSPAN;
		  nu_max = nu*(double)1e9*(double)GSPAN;
		  gama_min = sqrt((nu_min)/scale_gam_nu);
		  gama_max = sqrt((nu_max)/scale_gam_nu);
		  xl = (double)gama_min;
		  xu = (double)gama_max;
	  
		  fnorm = pow((double)10.0,(double)norm_param[j]);
		  alpha1 = pow(10.0,(double)alpha1_param[j]);
		  alpha2 = pow(10.0,(double)dalph_param[j])+pow(10.0,(double)alpha1_param[j]);
		  nu_break = pow(10.0,(double)nubrk_param[j]);
		  gama_break = sqrt((nu_break)/scale_gam_nu);
		  xb = (double)gama_break;
		  Tx = pow((double)10.0,(double)Tx_param[j]);
		  Te = pow((double)10.0,(double)Te_param[j]);
		  nu_t = pow((double)10.0,(double)nut_param[j]);

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

	      fprintf(fptr2,"\n");
	      /********************************************/
	    }
	  /*  If data requires alpha2 flatter than alpha1, then model as sum of power laws, i.e. steep and flat spectrum sources */
	  else
	    {
	      printf("alpha2 flatter than alpha1!\n");
	      
	      fnorm1 = norm_param[j];
	      fnorm2 = nubrk_param[j];
	      alpha1 = pow(10.0,(double)alpha1_param[j]);
	      alpha2 = alpha1 - pow(10.0,(double)dalph_param[j]);
	      Tx = pow((double)10.0,(double)Tx_param[j]);
	      Te = pow((double)10.0,(double)Te_param[j]);
	      nu_t = pow((double)10.0,(double)nut_param[j]);
	
	      for(ii=1; ii<=NFREQ; ii=ii+1)
		{
		  
		  nu = (double)freq_GHz[ii];
		  extn = exp((double)-1.0*pow((nu_t/nu),(double)2.1));
		  output = fnorm1*(pow(nu,-1.0*alpha1) + fnorm2*(pow(nu,-1.0*alpha2)) + Tx*pow(nu,-2.1))*extn + Te*((double)1.0 - extn);
		  fprintf(fptr2,"%20.15lf ",output);
		}
	      fprintf(fptr2," \n");

	    }
	}    
  fclose(fptr);
  fclose(fptr2);
  free_vector(norm_param,1,NHPIX);
  free_vector(alpha1_param,1,NHPIX);
  free_vector(dalph_param,1,NHPIX);
  free_vector(nubrk_param,1,NHPIX);
  free_vector(Tx_param,1,NHPIX);
  free_vector(Te_param,1,NHPIX);
  free_vector(nut_param,1,NHPIX);
  free_dvector(freq_GHz,1,NFREQ);
}

