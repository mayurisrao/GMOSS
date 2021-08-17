/* cal_lst_2011.c
*	
* compute LST - the local apparent sidereal time - from UTC (as Julian date) and longitude
* assumes DUT1 = UT1-UTC = 0.3 sec which is correct late 2012
*
* using equations from aa.usno.navy.mil/faq/docs/GAST.php
*
*/

#include <stdio.h>
#include <math.h>
#define PI 3.14159265358979

double cal_lst(double utc_julian_date, double longitude_radians)
{
	int julian_day,day_of_year;
	long int yyyy,mmmm,dddd,jjdd;
	double ut1_julian_date,time_ut1,gmst,lst,gast;
	double dd,TT,HH,JD0,D0,DD;
	double ee, LL, OO, DP, eqeq;

	ut1_julian_date = utc_julian_date + (double)0.3/(3600.0*24.0); 
	/* assume DUT1 = UT1-UTC = 0.3 sec which is correct late 2012 */

	julian_day = (int)floor(ut1_julian_date + 0.5);
	time_ut1 = (ut1_julian_date+0.5) - (double)julian_day;
	if(time_ut1 >= 1.0) time_ut1 -= 1.0;
	if(time_ut1 < 0.0) time_ut1 += 1.0;
	JD0 = (double)julian_day - 0.5;
	HH = (ut1_julian_date - (double)JD0)*24.0;
	D0 = JD0 - 2451545.0;
	DD = ut1_julian_date - 2451545.0;
	TT = DD/36525.0;

	gmst = 6.697374558 + 0.06570982441908 * D0 + 1.00273790935 * HH + 0.000026 * TT*TT;
	while(gmst>=24.0) gmst -= 24.0;
	while(gmst<0.0) gmst += 24.0;

//	printf("gmst2 = %g\n",gmst);

	ee = 23.4393 - 0.0000004 * DD;
	LL = 280.47 + 0.98565 * DD;
	OO = 125.04 - 0.052954 * DD;
	
	ee *= PI/180.0;
	LL *= PI/180.0;
	OO *= PI/180.0;

	DP = -0.000319*sin(OO) - 0.000024*sin(2*LL);
	eqeq = DP * cos(ee);

	gast = gmst + eqeq;

	lst = gast + (longitude_radians*180.0/PI)*(24.0/360.0);
	if(lst>24.0) lst -= 24.0;
	if(lst<0.0) lst += 24.0;

	lst *= (360.0*PI)/(24.0*180.0);
	
	return lst;
}

