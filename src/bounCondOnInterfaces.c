#include "type.h"
#include "def.h"      /* Definitions, parameters */
#include "global.h"   /* global variables */

#include "bounCondOnInterfaces.h"

/***********************
* BOUNCONDONINTERFACES 
*  Completing the "outer" parameters' values at the domain boundary
***********************/
void BounCondOnInterfaces( void )
{

	/*--- x ---*/

	/*- left side - plain INFLOW */
	P = 100000.;
	U = 76.4; /* V = 0.; W = 0.; */
	R = 1.0;

	for( j = 0; j < HIG/2; j++ ) {
		for( k = 0; k < DEP; k++ ) {
			xU1[0][j][k] = R;
			xU2[0][j][k] = R * U;
			xU3[0][j][k] = 0./*R * V*/;
			xU4[0][j][k] = 0./*R * W*/;
			xU5[0][j][k] = P / K_1 + 0.5 * R * ( U * U /*+ V * V + W * W*/ );
		} /* end for() */
	} /* end for() */
	/*- left side - block of disturbed velocities on INFLOW */
	for( j = HIG/2; j < HIG/2 + BL_HIG; j++ ) {
		for( k = 0; k < DEP; k++ ) {
			if( k_max > k_min )
				if( k >= k_min && k < k_max ) {
					   U = 200. + Ud, V = Vd, W = Wd;
				}
				else { U = 200., V = W = 0.; }
			else
				if( k >= k_max && k < k_min ) {
					   U = 200., V = W = 0.;
				}
				else { U = 200. + Ud, V = Vd, W = Wd; }
			xU1[0][j][k] = R;
			xU2[0][j][k] = R * U;
			xU3[0][j][k] = R * V;
			xU4[0][j][k] = R * W;
			xU5[0][j][k] = P / K_1 + 0.5 * R * ( U * U + V * V + W * W );
		} /* end for() */
	} /* end for() */
	/*- left side - plain INFLOW */
	U = 200.; /* V = 0.; W = 0.; */
	for( j = HIG/2 + BL_HIG; j < HIG; j++ ) {
		for( k = 0; k < DEP; k++ ) {
			xU1[0][j][k] = R;
			xU2[0][j][k] = R * U;
			xU3[0][j][k] = 0./*R * V*/;
			xU4[0][j][k] = 0./*R * W*/;
			xU5[0][j][k] = P / K_1 + 0.5 * R * ( U * U /*+ V * V + W * W*/ );
		} /* end for() */
	} /* end for() */
	/* right side - OUTFLOW */
	P = 100000.;
	for( j = 0; j < HIG; j++ ) {
		for( k = 0; k < DEP; k++ ) {
			U1x[LEN][j][k] = R = xU1[LEN][j][k];
			U2x[LEN][j][k] = U = xU2[LEN][j][k];
			U3x[LEN][j][k] = V = xU3[LEN][j][k];
			U4x[LEN][j][k] = W = xU4[LEN][j][k];
			U5x[LEN][j][k] = P / K_1 + 0.5 * ( U * U + V * V + W * W ) / R;
		} /* end for() */
	} /* end for() */

	/*--- y ---*/

	for( i = 0; i < LEN; i++ ) {
		for( k = 0; k < DEP; k++ ) {
			/* bottom side - SLIP */
			yU1[i][0][k] =   U1y[i][0][k];
			yU2[i][0][k] =   U2y[i][0][k];
			yU3[i][0][k] = - U3y[i][0][k];
			yU4[i][0][k] =   U4y[i][0][k];
			yU5[i][0][k] =   U5y[i][0][k];
			/* top side - SLIP */
			U1y[i][HIG][k] =   yU1[i][HIG][k];
			U2y[i][HIG][k] =   yU2[i][HIG][k];
			U3y[i][HIG][k] = - yU3[i][HIG][k];
			U4y[i][HIG][k] =   yU4[i][HIG][k];
			U5y[i][HIG][k] =   yU5[i][HIG][k];
		} /* end for */
	}  /* end for */

	/*--- z ---*/

	for( i = 0; i < LEN; i++ ) {
		for( j = 0; j < HIG; j++ ) {
			/* back side - PERIODIC */
			zU1[i][j][0] = zU1[i][j][DEP];
			zU2[i][j][0] = zU2[i][j][DEP];
			zU3[i][j][0] = zU3[i][j][DEP];
			zU4[i][j][0] = zU4[i][j][DEP];
			zU5[i][j][0] = zU5[i][j][DEP];
			/* front side - PERIODIC */
			U1z[i][j][DEP] = U1z[i][j][0];
			U2z[i][j][DEP] = U2z[i][j][0];
			U3z[i][j][DEP] = U3z[i][j][0];
			U4z[i][j][DEP] = U4z[i][j][0];
			U5z[i][j][DEP] = U5z[i][j][0];
		} /* end for */
	} /* end for */
} /* end BounCondOnInterfaces( void ) */
