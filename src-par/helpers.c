#include "mpi.h"
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
void checkCoNum( int myid ){

  unsigned int i, j, k;
  real maxCo,CoNum,gloMaxCo;

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

  // Allreduce maximum for maxCo
  MPI_Allreduce( &maxCo, &gloMaxCo, 1, MPI_FLOAT, MPI_MAX, MPI_COMM_WORLD );

  if (myid == 0) fprintf(stdout, "Maximum Courant no. = %f\n", maxCo );

  maxCo = gloMaxCo;

  /* Perform Courant number asjustment to timestep dt*/
  deltaT *= ( maxCoNum / maxCo); /* adjust time step so max Courant number doesn't exceed maxCoNum */

  if (myid == 0) fprintf(stdout, "New time-step size:  %f\n", deltaT );

  // other adjustments as well:
  deltaT_X = deltaT / deltaX;     deltaT_Y = deltaT / deltaY;     deltaT_Z = deltaT / deltaZ;

} /* end checkCoNum() */



/************
*  ARRAY2D  *   Memory allocation procedure for 2D arrays
************/
real **Array2D(unsigned columns, unsigned rows)
{
real **x;
unsigned i;

	if((x=(real **)malloc(columns*sizeof(real*))) == NULL) {
		fprintf(stderr, "mpi_duct: can't allocate memory");
		MPI_Abort(MPI_COMM_WORLD, 1);
	}
	for (i = 0; i < columns; i++)
		if ((x[i]=(real *)malloc(rows*sizeof(real))) == NULL) {
			fprintf(stderr, "mpi_duct: can't allocate memory");
			MPI_Abort(MPI_COMM_WORLD, 1);
		}
	return x;

} // end Array2D()

/************
*  ARRAY3D  *   Memory allocation procedure for 3D arrays
************/
real ***Array3D(unsigned columns, unsigned rows, unsigned floors)
{
real ***x;
unsigned i, j;

	if( (x=(real ***)malloc(columns*sizeof(real**))) == NULL ) {
		fprintf(stderr, "mpi_duct: can't allocate memory");
		MPI_Abort(MPI_COMM_WORLD, 1);
	}
	for (i = 0; i < columns; i++) {
		if ((x[i]=(real **)malloc(rows*sizeof(real*))) == NULL) {
			fprintf(stderr, "mpi_duct: can't allocate memory");
			MPI_Abort(MPI_COMM_WORLD, 1);
		}
		for (j = 0; j < rows; j++)
			if ((x[i][j]=(real *)malloc(floors*sizeof(real))) == NULL) {
				fprintf(stderr, "mpi_duct: can't allocate memory");
				MPI_Abort(MPI_COMM_WORLD, 1);
			}
	}
	return x;

} // end Array3D()

  /*
 * Free Array 3D
 */
void free3D( real *** arr, int columns, int rows ){
  
  int i,j;

  for ( i = 0; i < columns; i++ ){
    for ( j = 0; j < rows; j++ ){
      free( arr[i][j] );
    }
    free( arr[i] );
    }
  free( arr );

} 


