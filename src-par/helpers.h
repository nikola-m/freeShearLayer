#ifndef HELPERS_H
#define HELPERS_H

/*
* Random - Generator of the pseudo-random "doubles"
*/
double Random( double x );

/*
* Check Courant number - Perfmorms check of a maximal Courant number for time-stepping stability reasons.
*/
void checkCoNum( int myid );


/************
*  ARRAY2D  *   Memory allocation procedure for 2D arrays
************/
real **Array2D(unsigned columns, unsigned rows);

/*
*  ARRAY3D - Memory allocation procedure (for flat memory isn'it inefficient?)
*/
real ***Array3D( unsigned columns, unsigned rows, unsigned floors );

/*
 * Free Array 3D
 */
void free3D( real *** arr, int columns, int rows );


#endif