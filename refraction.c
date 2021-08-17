/* Program to calculate the apparent altitude given the observed altitude and vice versa
   Author - Mayuri SRao
   Date - 15 Feb 2013
*/

/* #include<stdio.h> */
/* #include<math.h> */
/* # define alpha 0.0065  //temp lapse rate [deg C per meter] */
/* # define PI 3.14159265 */
/* double refraction (double rad_elevation, double m_altitude); */
/* int main() */
/* { */
/*   double elevation, new_elevation, altitude; */
/*   printf("Enter the true elevation in radians to be corrected for refraction\n"); */
/*   scanf("%lf",&elevation); */
/*   printf("Enter the altitude of the observing location in meters.. 0 being sea level\n"); */
/*   scanf("%lf",&altitude); */
/*   new_elevation = refraction(elevation, altitude); */
/*   printf("true_elevation = %g, \tobserved_elevation = %g\n",elevation,new_elevation); */
/*   return 0; */
/* } */


#include<stdio.h>
#include<math.h>
# define alpha 0.0065  //temp lapse rate [deg C per meter]
# define PI 3.14159265 
double refraction (double rad_elevation, double m_altitude)
{ 
  double d2r = PI/180.0;
  double deg_elevation;
  double tpcor,R,new_ele,pres,temp;
  pres = 1010.*pow((1-6.5/288000*m_altitude),5.255);
  if (m_altitude > 11000) 
    {temp = 211.5;}
  else 
    {temp = 283.0 - alpha*m_altitude;}
  deg_elevation = rad_elevation/d2r;
  //  printf("elevation in deg = %lf \n",deg_elevation);
  R = 1.02/tan((deg_elevation + (10.3/(deg_elevation + 5.11)))*d2r);//refraction correction in arcminutes
  if (deg_elevation == 90.0)
    { R += 0.0019279;}
  tpcor = pres/1010. * 283/temp;
  R = tpcor*(R/60);
  new_ele = deg_elevation + R;
  return new_ele*d2r;
}

