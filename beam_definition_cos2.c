/* beam_definition.c
 *
 * defines the beam pattern
 *
 * function gives values of the beam power pattern  
 *
 * inputs: 
 *         frequency in MHz
 *         azimuth in radians
 *         elevation in radians
 */

#include <stdio.h>
#include <math.h>

#define PI 3.14159265

double beam_definition_cos2(float frequency, double azimuth, double elevation)
{
	double beam;
	/* if(elevation<=(PI/4.0)) beam=0.0; */
	/* else */
	beam = cos(2.0*(elevation - PI/2))*cos(2.0*(elevation-PI/2));

    return(beam);
}
