#ifndef HELPERS_H
#define HELPERS_H

/*
* Random - Generator of the pseudo-random "doubles"
*/
double Random( double x );

/*
* Check Courant number - Perfmorms check of a maximal Courant number for time-stepping stability reasons.
*/
void checkCoNum( void );


/*
*  ARRAY3D - Memory allocation procedure (for flat memory isn'it inefficient?)
*/
real ***Array3D( int columns, int rows, int floors );

#endif