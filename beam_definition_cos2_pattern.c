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

double beam_definition_cos2_pattern(float frequency, double azimuth, double elevation)
{
	double beam;
	if(elevation<=0.0) beam=0.0;
	else
	beam = sin(elevation)*sin(elevation);

    return(beam);
}
