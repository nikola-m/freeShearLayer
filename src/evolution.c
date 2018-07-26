/*
*  EVOLUTION  
* Computes new (or intermediate) values of conservative variables
* using two-stage TVD Runge-Kutta algorithm
*/
#include "type.h"
#include "def.h"      /* Definitions, parameters */
#include "global.h"   /* global variables */

#include "evolution.h"

void Evolution( int Stage )
{
/*register*/ unsigned i, j, k, _i, _j, _k;

	if( Stage == 1 ) {

		for( i = 1, _i = 0; i < LENN; i++, _i++ ) {
			for( j = 1, _j = 0; j < HIGG; j++, _j++ ) { /* indices i-1, k-1 etc */
				for( k = 1, _k = 0; k < DEPP; k++, _k++ ) {

					U1p[i][j][k] = U1[i][j][k] 
					               + deltaT_X * ( Fx1[_i][_j][_k] - Fx1[i][_j][_k] )
								   + deltaT_Y * ( Fy1[_i][_j][_k] - Fy1[_i][j][_k] )
								   + deltaT_Z * ( Fz1[_i][_j][_k] - Fz1[_i][_j][k] );
					U2p[i][j][k] = U2[i][j][k] 
					               + deltaT_X * ( Fx2[_i][_j][_k] - Fx2[i][_j][_k] )
								   + deltaT_Y * ( Fy2[_i][_j][_k] - Fy2[_i][j][_k] )
								   + deltaT_Z * ( Fz2[_i][_j][_k] - Fz2[_i][_j][k] );
					U3p[i][j][k] = U3[i][j][k] 
					               + deltaT_X * ( Fx3[_i][_j][_k] - Fx3[i][_j][_k] )
								   + deltaT_Y * ( Fy3[_i][_j][_k] - Fy3[_i][j][_k] )
								   + deltaT_Z * ( Fz3[_i][_j][_k] - Fz3[_i][_j][k] );
					U4p[i][j][k] = U4[i][j][k] 
					               + deltaT_X * ( Fx4[_i][_j][_k] - Fx4[i][_j][_k] )
								   + deltaT_Y * ( Fy4[_i][_j][_k] - Fy4[_i][j][_k] )
								   + deltaT_Z * ( Fz4[_i][_j][_k] - Fz4[_i][_j][k] );
					U5p[i][j][k] = U5[i][j][k] 
					               + deltaT_X * ( Fx5[_i][_j][_k] - Fx5[i][_j][_k] )
								   + deltaT_Y * ( Fy5[_i][_j][_k] - Fy5[_i][j][_k] )
								   + deltaT_Z * ( Fz5[_i][_j][_k] - Fz5[_i][_j][k] );

				 } /* end for */
			} /* end for */
		} /* end for */

		/* next stage */
		Stage = 2;
		U1_ = U1p, U2_ = U2p, U3_ = U3p, U4_ = U4p, U5_ = U5p;

	} /* end if() First stage */
	else {

		for( i = 1, _i = 0; i < LENN; i++, _i++ ) {
			for( j = 1, _j = 0; j < HIGG; j++, _j++ ) {
				for( k = 1, _k = 0; k < DEPP; k++, _k++ ) {

					U1[i][j][k] = 0.5 * ( U1[i][j][k] + U1p[i][j][k]
								 + deltaT_X * ( Fx1[_i][_j][_k] - Fx1[i][_j][_k] )
								 + deltaT_Y * ( Fy1[_i][_j][_k] - Fy1[_i][j][_k] )
								 + deltaT_Z * ( Fz1[_i][_j][_k] - Fz1[_i][_j][k] )
								);
					U2[i][j][k] = 0.5 * ( U2[i][j][k] + U2p[i][j][k]
								 + deltaT_X * ( Fx2[_i][_j][_k] - Fx2[i][_j][_k] )
								 + deltaT_Y * ( Fy2[_i][_j][_k] - Fy2[_i][j][_k] )
								 + deltaT_Z * ( Fz2[_i][_j][_k] - Fz2[_i][_j][k] )
								);
					U3[i][j][k] = 0.5 * ( U3[i][j][k] + U3p[i][j][k]
								 + deltaT_X * ( Fx3[_i][_j][_k] - Fx3[i][_j][_k] )
								 + deltaT_Y * ( Fy3[_i][_j][_k] - Fy3[_i][j][_k] )
								 + deltaT_Z * ( Fz3[_i][_j][_k] - Fz3[_i][_j][k] )
								);
					U4[i][j][k] = 0.5 * ( U4[i][j][k] + U4p[i][j][k]
								 + deltaT_X * ( Fx4[_i][_j][_k] - Fx4[i][_j][_k] )
								 + deltaT_Y * ( Fy4[_i][_j][_k] - Fy4[_i][j][_k] )
								 + deltaT_Z * ( Fz4[_i][_j][_k] - Fz4[_i][_j][k] )
								);
					U5[i][j][k] = 0.5 * ( U5[i][j][k] + U5p[i][j][k]
								 + deltaT_X * ( Fx5[_i][_j][_k] - Fx5[i][_j][_k] )
								 + deltaT_Y * ( Fy5[_i][_j][_k] - Fy5[_i][j][_k] )
								 + deltaT_Z * ( Fz5[_i][_j][_k] - Fz5[_i][_j][k] )
								);

				 } /* end for */
			} /* end for */
		} /* end for */

		/* next stage */
		Stage = 1;
		U1_ = U1, U2_ = U2, U3_ = U3, U4_ = U4, U5_ = U5;

	} /* end else Second stage */

} /* end Evolution() */
