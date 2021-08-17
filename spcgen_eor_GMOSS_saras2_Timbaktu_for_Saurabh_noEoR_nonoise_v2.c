/* spcgen.c

 *
 * model the sky, model the instrument, generate the measurement set
 *
 * spectra are written at each time assuming sky is stationary in each integration.
 * spectra are written at intervals of INTEGRATION_TIME
 * 
 * Ravi: original - 07 Dec 2012
 * Modification: 19may16 mayuri - changed number of channels to 167 to extend up to ~260 MHz with 5 extra log spaced 
 * freq points beyond 200 MHz
 */

// Comment/uncomment path to nrutil.h

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include </home/mayuris/software/jnk/miriad/inc/linux64/maxdimc.h>
//#include </Users/mayuri/work/C/ANSI-C/nrutil.h>
#include </home/mayuris/work/HEALPy/codes/broken_pl_28dec15/nrutil.h>
//#include </Users/rsubrahm/CODESPACE/INCLUDE/nrutil.h>
//#include </usr/local/miriad/inc/miriad.h>
#include </home/mayuris/software/jnk/miriad/inc/miriad.h>
#include <string.h>
//#include <cpgplot.h>

#define CHANNEL_WIDTH 1.0 /* MHz */
#define NUMBER_OF_CHANNELS 103 /* 103 depending on contents of frequency file */
#define START_FREQUENCY 0.04 /* 0.01 GHz = 10 MHz */
#define INTEGRATION_TIME 5.0 /*3600.0*//* sec 	This is the UT time interval in which
				   spectra are written out */ 

/* #define NOISE_INT_TIME 5000.0 *//* seconds   4 hour noise integration
				This is the time used for computation
				of noise added to the spectrum ... using time for just over 3 sigma*/ 
#define	SITE_LATITUDE +32.7908 // GBD: +13.01333 /* deg *// Chile: -23.0228// Hanle: +32.7908 // Timbaktu +14.2331114
#define	SITE_LONGITUDE 79.0002 //* deg */// GBD: 77.580833 Chile: -67.7550//  Hanle: 79.0002 // Timbaktu 77.6108055
#define site_altitude 4500 /* 200.0 */ /* meters *///4800//
#define NHPIX 3072
#define PI 3.14159265358979
#define NN 5000 // length of array containing recombination line spectra

#define Trx 50.0 // K Not really used anywhere further in the code..
#define T_n1 50.0 // K Receiver noise temperature in path 1 in cross-corr receiver
#define T_n2 50.0 // K Receiver noise temperature in path 2 in cross-corr receiver
#define T_ref 300.0 //K reference nosie source
#define T_atm 0.0113044 //K ( atmospheric correction - increase on including absorption 
//							and translating T_sys (Trx + atmosphere) to above atmosphere)
#define TCMB 2.72548
#define T_cold 2730.0
#define T_hot 3730.0 
#define d2r PI/180.0
#define kk 1.3806488e-23 // m^2 kg s^-1
#define hh 6.62606957e-34 // m^2 kg s^-2 K^-1
#define cvel 2.99792458e+08 // m s^-1
#define FTOL 1.0e-4
#define T_e 8000.0 // Thermal electron temperature
#define GSPAN 100.0
#define m_e 9.1e-31
#define q_e 1.6e-19
#define sin_alph 1.0
#define Bmag 1e-9 // Tesla == 10 micro-Gauss
#define NFREQ 103

//static const char filename_chisq[]= "eor_chisq_values_map_interp_23apr15.txt";

double cal_lst(double utc_julian_date, double longitude_radians);
void precess_(float *ra, float *dec, float *epoch1, float *epoch2);
double beam_definition_cos2_pattern(float frequency, double azimuth, double elevation);

//double beam_definition_saras2(float frequency, double azimuth, double elevation);

float gasdev(long *idum); 
float gammq(float a, float x);

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

/* Function to account for effects of atmospheric refraction */
double refraction (double rad_elevation, double m_altitude);

/* Functions used in generating SARAS beam*/
void fpoly(float x, float p[], int np);
void lfit(float x[], float y[], float sig[], int ndat, float a[], int ia[],
	  int ma, float **covar, float *chisq, void (*funcs)(float, float [], int));
void covsrt(float **covar, int ma, int ia[], int mfit);
void gaussj(float **a, int n, float **b, int m);

/* Global Variables */
float *xx, *yy;
double xnu,xri,xrk,xrip,xrkp;
double nu,alpha1,alpha2,output,fnorm,fnorm1,fnorm2;
double AA,BB,C1,C2,D,ft,therm,fbreak;
double nu_break;
float x,y;
float f1,f2;
double temp1,temp2;
float inttolf = 1.0e-8;
float difft[5];
double gama_min,gama_max,gama_break,scale_gam_nu;
int nrecurs = 10000000;

/* Evaulate bessel function of order 5/3, integrand  */

main()
{
  const char version[30]="spcgen: 21apr15";
  static const char filename_freq[] = "saras2/freq_revised_16.txt";
  static const char filename_coord[] = "../../data/all_sky_maps_5deg/PIXEL_LISTING_NESTED_R4_GALACTIC.txt";
  static const char filename_eor[]= "../../data/eor_sim/eor_Tb_values.dat";
  static const char filename_beam[]= "monopole_beam.txt";
  static const char filename_fit_params[]="../../data/eor_sim/total_model_params_11apr16_allpix.txt";
  static const char filename_pix_spec[]="../../data/eor_sim/total_model_spec_13jan17_1_3072_used.txt";
  static const char filename6[]="../../data/eor_sim/dummy_spectrum_14jan17.txt";  // corresponds to LST~0.25 hrs

  double rint, rint1, rint2, xl, xu, xb, reps, aeps, dif;
  double *freq_GHz;
  int ier, npt, nmax, ctt;
  float *yamoeba,*afit,*aafit,*cof,**pfun,*yfit,*beam_elevation,*beam_gain,*beam_sigma,*a_lfit,*p_fpoly;
  float T_150,T_408,T_1420,T_23G,gama_1,gama_2,gama_3,nu_min,nu_max;
  time_t utime;
  long seed;  
  float *xarray,*yarray,*y2;
  float read1,read2;
  float ffmin,ffmax;
  float yp1,ypn;
  double new_yf[NUMBER_OF_CHANNELS];
  double cmb_intensity[NUMBER_OF_CHANNELS],P_cold[NUMBER_OF_CHANNELS],P_hot[NUMBER_OF_CHANNELS],fg_int;
  double P_diff[NUMBER_OF_CHANNELS];
  double P_diff_temp[NUMBER_OF_CHANNELS];
  float yf;
  int npoints1,nfunc;
  char outfile[20];
  char inbuf[200];
  int i,j,k,ii,jj,jjj,itst;
  int nhpix,nbadpixels,beam_points,ma_lfit,ma_fpoly;
  int npolyfit;
  float read_value1,read_value2,read_value3,read_value4,read_value5,read_value6,read_value7,read_value8,read_value9,read_value10,read_value11,read_value12;
  double read_value;
  double cwt;
  double cfreq,cfreq1,cfreq0;
  double ctemp;
  float ra_precess,dec_precess,epoch1,epoch2;
  float *ll_coordinate, *bb_coordinate;
  int *flag,*ia_lfit;
  float **cc,**covar_lfit,chisq_lfit;
  double **ctemp_val;
  double intensity, cfreq1_hz;
  double final_temp; 
  double final_temp_1, final_temp_2;
  double cof0,cof1,cof2,cof3;
  FILE *fptr,*fptr1,*fptr2,*fptr3,*fptr4,*fptr5,*fptr6,*fptr7;
  double ft;

  /* variables for output vis file */

  char source_name[5];
  char velocity_type[9] = {'V','E','L','O','-','O','B','S','\0'};
  int tno;
  int var_ivalue[1];
  int npoints,nn;
  int antenna1,antenna2;
  int flags[NUMBER_OF_CHANNELS];
  float baseline_value[1];
  double preamble[4];
  double p1,p2,p3;
  double ant_az[2],ant_el[2];
  double freq_channel1[1];
  double freq_inc[1];
  double time_var[1];
  double coord_var[2];
  double site_latitude[1],site_longitude[1];
  float data[2*NUMBER_OF_CHANNELS];
  double sumwt[NUMBER_OF_CHANNELS],sumwt0[NUMBER_OF_CHANNELS],sumwt1[NUMBER_OF_CHANNELS];
  float var_veldop[1];
  float var_vsource[1];
  double var_restfreq[1];
  double sra[1],sdec[1];
  double lo1[1],lo2[1],freq[1],freqif[1];
  int mount[1];
  float evector[1];
  int nnpol;
  float jyperk[1],inttime[1],epoch[1];
  double antpos[6];
  float tpower[1];
   
  long int jjdd,yyyy,mmmm,dddd,hrhr,minmin;
  long int nspect;
  long int int_time;
  double julian_date, longitude_radians, lst;
  double cross_imag[NUMBER_OF_CHANNELS],cross_real[NUMBER_OF_CHANNELS],cross_real0[NUMBER_OF_CHANNELS],cross_real1[NUMBER_OF_CHANNELS],cross_real2[NUMBER_OF_CHANNELS];
  double sigma, variance; 
  float secsec;
  double c_ll, c_bb;
  double ra_b1950, dec_b1950;
  double ra_date, dec_date;
  double azaz, elel, new_elel,x_fpoly;
  double haha;

  double time_utc,time_utc_res;
  int time_utc_hh,time_utc_mm,julian_day;
  float time_utc_sec;	

  float AA,BB,D,FFIT,DDIF,chisq;
  float correction_150offset = 21.4; // Kelvin, Patra et al. 2015
  float correction_150scaling = 1.05; // Patra et al. 2015
  // Variable for pgplot
  int npointsp;
  float xmin,xmax,ymin,ymax,xxc[1],yyc[1];
  

  ll_coordinate = vector(1,NHPIX);
  bb_coordinate = vector(1,NHPIX);
  freq_GHz = dvector(1,NFREQ);
  flag = ivector(1,NHPIX);

  ctemp_val = dmatrix(1,NHPIX+1,1,NUMBER_OF_CHANNELS);
  nhpix = NHPIX;

  fptr=fopen(filename_freq,"r");
  printf("Opened file to read frequencies \n");
  i=1;
  while (fgets (inbuf, sizeof inbuf, fptr) != NULL)
    {
      sscanf(inbuf,"%lf",&read_value);
      freq_GHz[i]=(double)read_value/1e3;
      if(i == 10) printf("\n at %d %20.15lf\n",i,(double)read_value/1e3);
      if(i == NFREQ) printf("at %d %20.15lf\n",i,(double)read_value/1e3);
      ++i;         
    }
  printf("\n");
  fclose(fptr);
  
  fptr=fopen(filename_coord,"r");
  //  fgets (inbuf, sizeof inbuf, fptr);  // dummy read
  i=1;
  while (fgets (inbuf, sizeof inbuf, fptr) != NULL)
    {
      sscanf(inbuf,"%g%g",&read_value1,&read_value2);
      ll_coordinate[i]=read_value1*d2r;
      bb_coordinate[i]=read_value2*d2r;
      ++i;
    }
  fclose(fptr);
  printf(" Read %d coordinate values \n",i);

  fptr=fopen(filename_fit_params,"r");
  //fgets (inbuf, sizeof inbuf, fptr);  // dummy read
  i=1;
  while (fgets (inbuf, sizeof inbuf, fptr) != NULL)
    {
      sscanf(inbuf,"%f%f%f%f%f%f%f%f%f%f",&read_value1,&read_value2,&read_value3,&read_value4,&read_value5,&read_value6,&read_value7,&read_value8,&read_value9,&read_value10);
      flag[i]=(int)read_value2;
      if(i == 1000) printf("%f  %f  %f  %f  %f  %f  %f  %f  %f  %f\n\n",read_value1,read_value2,read_value3,read_value4,read_value5,read_value6,read_value7,read_value8,read_value9,read_value10);
      ++i;
    }
  fclose(fptr);
  printf(" Read %d flag values \n",i);

  int fflush(FILE *fptr);
  fptr=fopen(filename_pix_spec,"r");
  //fgets (inbuf, sizeof inbuf, fptr);  // dummy read
  j=1;
  while (fgets (inbuf, sizeof inbuf, fptr) != NULL)
    {
      //      printf("%d\n",j);
      i=1;
      while(i <= NUMBER_OF_CHANNELS)
	{
	  fscanf(fptr,"%lf",&ctemp_val[j][i]);
	  if(j == 1) printf("%20.15f ",ctemp_val[j][i]);
	  if(j == 3072) printf("%20.15f ",ctemp_val[j][i]);
	  //	  ctemp_val[j][i] = temp1;
	    ++i;
	} 
      if(j == 1) printf("\n ");
      if(j == 3072) printf("\n ");
      ++j;            
    }
  printf("\n");
  fclose(fptr);
  printf(" Read %d spectra values \n",j);


  // Open and read input file containing beam of SARAS2 antenna (as on 7 May 2016)
  fptr=fopen(filename_beam,"r");
  beam_points=0;
  beam_elevation=vector(1,20);
  beam_gain=vector(1,20);
  beam_sigma=vector(1,20);
  while (fgets (inbuf, sizeof inbuf, fptr) != NULL)
    {
      ++beam_points;
      sscanf(inbuf,"%g%g",&read_value1,&read_value2);
      beam_elevation[beam_points]=read_value1;
      beam_gain[beam_points]=read_value2;
  	  beam_sigma[beam_points]=0.1;
    }
  fclose(fptr);
  printf(" Read %d values from beam measurements file \n",beam_points);


  // Generate a functional form for the beam pattern to be used to get beam weights
  ia_lfit=ivector(1,7);
  a_lfit=vector(1,7);
  covar_lfit = matrix(1,7,1,7);
  p_fpoly = vector(1,7);
  
  ma_lfit=7;
  ma_fpoly=7;
  for (i=1;i<8;++i) ia_lfit[i]=1;
  lfit(beam_elevation, beam_gain, beam_sigma, beam_points, a_lfit, ia_lfit,
       ma_lfit, covar_lfit, &chisq_lfit, fpoly);
  
  free_vector(beam_elevation,1,20);
  free_vector(beam_gain,1,20);
  free_vector(beam_sigma,1,20);
  free_ivector(ia_lfit,1,7);
  free_matrix(covar_lfit,1,7,1,7);

// SET SITE LATITUDE AND LONGITUDE
  site_latitude[0] = (double)(SITE_LATITUDE * PI/180.0);
  site_longitude[0] = (double)(SITE_LONGITUDE * PI/180.0);

  // Open and initiatize the measurement set output file 

  /* set all flags to good */
  for(i=0;i<NUMBER_OF_CHANNELS;++i) flags[i]=1;

  npoints=NUMBER_OF_CHANNELS; /* number of spectral points in the records */

  printf("Type output vis file name : ");
  gets(inbuf);
  sscanf(inbuf,"%s",outfile);
	

  int_time=INTEGRATION_TIME;
  printf(" Type number of %ld-sec spectral records to be written: ",int_time);
  gets(inbuf);
  sscanf(inbuf,"%ld",&nspect);

  printf("Type start UTC yyyy, mm (unit offset), dd, hh, mm, ss.ss :");
  gets(inbuf);
  sscanf(inbuf,"%ld%ld%ld%ld%ld%f",&yyyy,&mmmm,&dddd,&hrhr,&minmin,&secsec);
  jjdd = ( 1461 * ( yyyy + 4800 + ( mmmm - 14 ) / 12 ) ) / 4 +
    ( 367 * ( mmmm - 2 - 12 * ( ( mmmm - 14 ) / 12 ) ) ) / 12 -
    ( 3 * ( ( yyyy + 4900 + ( mmmm - 14 ) / 12 ) / 100 ) ) / 4 +
    dddd - 32075;
  julian_date = (double)jjdd -0.5;
  julian_date += (double)secsec/(3600.0*24.0) + 
    (double)minmin/(60.0*24.0) + (double)hrhr/24.0;

  julian_date -= (double)INTEGRATION_TIME/(3600.0*24.0);  
  // backoff timestamp by one integration
  // time because time is incrememted first
  // off in the loop below
  	
  /* Read the file containing EoR spectrum information */
  
  ffmin=0.01;  // 0.01 GHz.. 10 MHz
  ffmax=0.20; // 0.2 GHz.. 200 MHz
  xarray=vector(1,NN);
  yarray=vector(1,NN);
  y2=vector(1,NN);
  npoints1=0;
  fptr1  = fopen(filename_eor,"r");
  while(fscanf(fptr1,"%g%g",&read1,&read2) != EOF)
    {
      if(read1 >= ffmin && read1 <= ffmax)
  	{
  	  ++npoints1;
  	  xarray[npoints1]=read1;
  	  yarray[npoints1]=read2;
  	}
    }
  fclose(fptr1);
  printf ("Got the EoR /del T_b values in range %f to %f GHz\n",ffmin,ffmax);

  /* /\* calculation of second order differential for spline interpolation *\/ */

  yp1=(yarray[2]-yarray[1])/(xarray[2]-xarray[1]);
  ypn=(yarray[npoints1]-yarray[npoints1-1])/(xarray[npoints1]-xarray[npoints1-1]);
  spline(xarray, yarray, npoints1, yp1, ypn, y2);
 
  //compute the coefficient arrays here for spectral interpolations fptr2 
  // to store spectral indices that are outliers, if any
  scale_gam_nu = (double)(3.0*q_e*Bmag*sin_alph)/(4.0*PI*m_e*cvel);

  // Open ascii file, 
  // ftpr1 to store all the simulated spectra for python processing

  int fflush(FILE *fptr1);

  fptr1  = fopen(filename6,"a+");
  for(i=0;i<NUMBER_OF_CHANNELS;++i)
    {
      cfreq = freq_GHz[i+1];      
      //      if (i==10) printf("Frequency check.. at index 10, frequency in GHz is %lf\n",cfreq);


      fprintf(fptr1,"%20.15f ",cfreq);
      /* fprintf(fptr3,"%f ",cfreq); */
      /* fprintf(fptr4,"%f ",cfreq); */
      //  fprintf(fptr2,"%f ",cfreq);
      /* compute the recombination spectrum intensity at the current frequency*/
      splint(xarray,yarray,y2,npoints1,cfreq,&yf);
      new_yf[i]=(double)yf;

      cfreq1_hz = cfreq*1e9;
      //uncomment below to include EoR
      /* cmb_intensity[i] = (2.0*(double)hh*cfreq1_hz*cfreq1_hz*cfreq1_hz)/ */
      /* 	     (((double)cvel*(double)cvel)* */
      /* 	(exp(((double)hh*(double)cfreq1_hz)/((double)kk*(double)(TCMB+new_yf[i])))-1.0)); */
      // uncomment below to not inclue EoR
      cmb_intensity[i] = (2.0*(double)hh*cfreq1_hz*cfreq1_hz*cfreq1_hz)/
      	(((double)cvel*(double)cvel)*
      	 (exp(((double)hh*(double)cfreq1_hz)/((double)kk*(double)(TCMB)))-1.0));
      P_cold[i]=2*kk*T_cold*(((hh*cfreq1_hz)/(kk*T_cold))/(exp((hh*cfreq1_hz)/(kk*T_cold))-1));
      P_hot[i]=2*kk*T_hot*(((hh*cfreq1_hz)/(kk*T_hot))/(exp((hh*cfreq1_hz)/(kk*T_hot))-1));
      P_diff[i]=P_hot[i]-P_cold[i];
    }
  
  fprintf(fptr1," \n");
  /* fprintf(fptr3," \n"); */
  /* fprintf(fptr4," \n"); */
  //  fprintf(fptr2," \n");

  /* initialize seed for random number generator for adding gaussian noise */
  time(&utime);
  seed = -(long)utime;
  /* seed_count += 1; */
  printf(" Start loop for writing nspect spectral records \n");

  for(ii=0;ii<nspect;++ii)
    {
      //fprintf(fptr2,"Spec %d\n",ii);
      // write the header record here for the current time
      
      source_name[0]='S';
      source_name[1]='P';
      source_name[2]='C';
      source_name[3]='0';
      source_name[4]='\0';
      //      uvputvr_c(tno,1,"source",source_name,4);

      preamble[0]=0.0;  /* u coordinate */
      preamble[1]=0.0;  /* v coordinate */

      // assemble the current JD in preamble[2] from UTC

      julian_date += (double)INTEGRATION_TIME/(3600.0*24.0);
      preamble[2] = julian_date;

      // Type UTC on terminal

      julian_day = (int)floor(julian_date + 0.5);
      time_utc = (julian_date+0.5) - (double)julian_day;
      if(time_utc >= 1.0) time_utc -= 1.0;
      if(time_utc < 0.0) time_utc += 1.0;
      time_utc *= 24.0;
      time_utc_hh = (int)floor(time_utc);
      time_utc_res = time_utc - (double)time_utc_hh;
      time_utc_mm = (int)floor(60.0*time_utc_res);    
      time_utc_sec = (double)(60.0*(60.0*time_utc_res-time_utc_mm));
      // printf("%d  ",ii+1);
      printf(" UTC: %2d %2d %5.2f\n",time_utc_hh,time_utc_mm,time_utc_sec);
      //      fprintf(fptr2," UTC: %2d %2d %5.2f\n",time_utc_hh,time_utc_mm,time_utc_sec);
      antenna1=1;
      antenna2=2;
      preamble[3]=(double)(256*antenna1+antenna2);

      // compute LST and record this variable

      longitude_radians = (double)site_longitude[0];
      lst=cal_lst(julian_date, longitude_radians);
      time_var[0]=lst;
      //      uvputvrd_c(tno,"lst",time_var,1);

      // prepare and write the spectral data

      for(i=0;i<NUMBER_OF_CHANNELS;++i)
  	{
  	  cross_real[i]=0.0;
	  cross_real2[i] = 0.0;
  	  cross_imag[i]=0.0;
  	  sumwt[i]=0.0;
  	}
      //      ctt = 0;
      for(j=1;j<=NHPIX;++j)
  	{
	  //	  printf("Pixel number : %d\n",j);
  	  c_ll=(double)ll_coordinate[j];
  	  c_bb=(double)bb_coordinate[j];
  	  cwt=0.0;
  	  ctt = 0;
	  // Convert ll,bb to ra,dec B1950.0 epoch
  	  // Wikipedia - Celestial coordinate systems
	  
  	  ra_b1950 = atan2( sin(c_ll - 123.0*PI/180.0),
  			    (cos(c_ll - 123.0*PI/180.0)*sin(27.4*PI/180.0)-tan(c_bb)*cos(27.4*PI/180.0)) )
  	    + 12.25*PI/180.0;
  	  dec_b1950 = asin( sin(c_bb)*sin(27.4*PI/180.0) +
  			    cos(c_bb)*cos(27.4*PI/180.0)*cos(c_ll-123.0*PI/180.0) );

  	  // Precess from B1950.0 to date

  	  ra_precess = (180.0/PI)*(float)ra_b1950;
  	  dec_precess = (180.0/PI)*(float)dec_b1950;
  	  epoch1=1950.0;
  	  epoch2=(float)yyyy + (float)(mmmm-1)/12.0 ;

  	  precess_(&ra_precess,&dec_precess,&epoch1,&epoch2);
  	  ra_date = (double)(ra_precess*PI/180.0);
  	  dec_date = (double)(dec_precess*PI/180.0);

  	  // Convert ra,dec date to az,el using the current LST
  	  // AZ is defined as the angle from N towards E

  	  haha = lst - ra_date; // Convert LST and RA to HA

  	  azaz = (double)(PI) + atan2( sin(haha) ,
  				       (cos(haha)*sin(site_latitude[0])-tan(dec_date)*cos(site_latitude[0])) );
  	  elel = asin( sin(site_latitude[0])*sin(dec_date)
  		       +cos(site_latitude[0])*cos(dec_date)*cos(haha) );

  	  /* Included calculation to correct for atmospheric refraction
  	     elel is in radians, new_elel is also in radians, site_altitude is in meters */
  	  new_elel = refraction(elel, site_altitude);

  	  // If az,el is within the beam weight by the beam pattern and accumulate the spectrum
  	  //	  printf("Pixel number: %d ",j);
  	  for(i=0;i<NUMBER_OF_CHANNELS;++i)
  	    {
  	      cfreq = freq_GHz[i+1];      
	      //	      if (i==10) printf("Frequency check.. at index 10, frequency in GHz is %lf\n",cfreq);
	      //cwt = beam_definition_cos2_pattern(cfreq,azaz,new_elel);
	     
	      x_fpoly = new_elel*180.0/3.14159265;
  	      cwt = 0.0;
	      fpoly(x_fpoly,p_fpoly,ma_fpoly);
	      for(k=1;k<=ma_fpoly;++k) {cwt += a_lfit[k]*p_fpoly[k];}
	      if (x_fpoly > 90.0 || x_fpoly < 5.0) cwt=0.0;
	     
  	      if(cwt > 0.0)
  		{

  		   if (ctt==0) 
  		    { 
		     /* printf("Pixel number: %d\n",j);  */
		      //		     fprintf(fptr2,"%d %d\n",flag[j],j);
  		     } 
		   ctt+=1; 
		  /* continue; */
  		  cfreq1=cfreq;
  		  cfreq1_hz = cfreq1*1e9;
  		  nu = cfreq1;
		  ctemp = (double)ctemp_val[j][i+1];
		   /* if (j == 257) */
		   /*   {printf("%f\n",ctemp);} */
  		  /* Convert sky temperature to intensity,
  		     add the recombination spectrum intensity to the computed sky brightness
  		     and record the final antenna temperature*/

  		  //uncomment below to convert ctemp in Planckian
  		  /* fg_int = (2.0*(double)hh*cfreq1_hz*cfreq1_hz*cfreq1_hz)/ */
  		  /*   (((double)cvel*(double)cvel)* */
  		  /*    (exp(((double)hh*(double)cfreq1_hz)/((double)kk*(ctemp)))-1.0)); */
  		  /* intensity = (cmb_intensity[i]+fg_int)*(cvel*cvel)/(cfreq1_hz*cfreq1_hz); */

  		  // uncomment below to convert ctemp in RJ
  		  intensity = (cmb_intensity[i]
  		  	       + (2.0*(cfreq1_hz)*(cfreq1_hz)*kk*ctemp/(cvel*cvel))
  		  	       )*(cvel*cvel)/(cfreq1_hz*cfreq1_hz);
		  
  		  // uncomment below for no CMB in spectra
  		  /* intensity =  ((2.0*(cfreq1_hz)*(cfreq1_hz)*kk*ctemp/(cvel*cvel)) */
  		  /* 		)*(cvel*cvel)/(cfreq1_hz*cfreq1_hz); */

  		  //ucnomment below for no foreground
  		  //intensity = (cmb_intensity[i])*(cvel*cvel)/(cfreq1_hz*cfreq1_hz);
		  
  		  
		  cross_real[i] += cwt*intensity;
		  sumwt[i] += cwt;		      		  
  		}
	           
  	    } // end loop for freq
	  //	  printf("%d  %f \n",j,cwt);
  	} // end loop for pixel

      for(i=0;i<NUMBER_OF_CHANNELS;++i)
      	{
      	  if(sumwt[i]>0.0)
      	    {
      	      cross_real[i] /= sumwt[i];
	      cross_real2[i]= (cross_real[i]/P_diff[i])*(T_hot-T_cold);
      	    }
	  /* fprintf(fptr3,"%20.15lf ",cross_real[i]); */
	  /* fprintf(fptr4,"%20.15lf ",sumwt[i]); */
      	}
 
      for(i=0;i<NUMBER_OF_CHANNELS;++i)
      	{
      	  data[2*i]=(float)cross_real2[i];
      	  data[2*i+1]=(float)cross_imag[i];
      	  fprintf(fptr1,"%20.15lf ",cross_real2[i]);
      	}

      fprintf(fptr1,"\n");

    } /* end of loop writing ispect records */
 
  fclose(fptr1);
  /* fclose(fptr3); */
  /* fclose(fptr4); */
  //  fclose(fptr2);   
  free_ivector(flag,1,NHPIX);
  free_vector(ll_coordinate,1,NHPIX);
  free_vector(bb_coordinate,1,NHPIX);
  free_dmatrix(ctemp_val,1,NHPIX+1,1,NUMBER_OF_CHANNELS);
  free_vector(xarray,1,NN);
  free_vector(yarray,1,NN);
  free_vector(y2,1,NN);
  free_dvector(freq_GHz,1,NFREQ);
  /* free_vector(a_lfit,1,7); */
  /* free_vector(p_fpoly,1,7); */
  return 0;
} /* end of main */



