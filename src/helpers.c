/* Helper functions */

#include <stdlib.h>    /* atoi()       */
#include <stdio.h>     /* printf() etc.*/
#include "type.h"
#include "def.h"
#include "global.h"
#include "helpers.h"

/*
* Random - Generator of the pseudo-random "doubles"
*/
double Random( double x )
{
	x = 11. * x + M_PI;
	return x - (double)((int)x);
} /* end Random */


/*
* Check Courant number - Perfmorms check of a maximal Courant number for time-stepping stability reasons.
*/
void checkCoNum( void ){

  unsigned int i, j, k;
  real maxCo,CoNum;

   maxCo = 0.0;

  /* Go trough inner field cells and perform check */
  for (k = 1; k < DEPP; k++) {
    for (j = 1; j < HIGG; j++) {
      for (i = 1; i < LENN; i++) {

 
        //  We will just check velocity in X-axis (streamwise) direction
        R = U1[i][j][k];
        U = U2[i][j][k]/R;
        maxCo = ( (CoNum = U*deltaT_X) > maxCo ) ? CoNum : maxCo;
        
      }
    }
  }
  printf( "Maximum Courant no. = %f\n", maxCo );

  /* Perform Courant number asjustment to timestep dt*/
  deltaT *= ( maxCoNum / maxCo); /* adjust time step so max Courant number doesn't exceed maxCoNum */
  printf( "New time-step size:  %f\n", deltaT );
  // other adjustments as well:
  deltaT_X = deltaT / deltaX;     deltaT_Y = deltaT / deltaY;     deltaT_Z = deltaT / deltaZ;

} /* end checkCoNum() */


/************
*  ARRAY3D  *   Memory allocation procedure (for flat memory isn'it inefficient?)
************/
real ***Array3D( int columns, int rows, int floors )
{
real ***x;
int i, j;

	if( (x=(real ***)malloc(columns*sizeof(real**))) == NULL ) {
	   puts( "Cannot allocate memory" );
	   exit( -1 );
	}
	for( i = 0; i < columns; i++ ) {
		if( (x[i]=(real **)malloc(rows*sizeof(real*))) == NULL ) {
		   puts( "Cannot allocate memory" );
		   exit(1);
		}
		for( j = 0; j < rows; j++ )
			if( (x[i][j]=(real *)malloc(floors*sizeof(real))) == NULL ) {
			   puts( "Cannot allocate memory" );
			   exit( -1 );
			}
	}
	return x;

} /* end Array3D() */