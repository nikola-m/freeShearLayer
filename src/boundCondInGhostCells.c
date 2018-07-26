#include "type.h"
#include "def.h"      /* Definitions, parameters */
#include "global.h"   /* global variables */

#include "bounCondInGhostCells.h"

/***********************
* BOUNCONDINGHOSTCELLS *   BC in the cells out of the domain
***********************/
void BounCondInGhostCells( void )
{
/*register*/ unsigned i, j, k;

		/*-- x  */

		/* left side - INFLOW */
		P = 100000.;
		U = 76.4; /* V = 0.; W = 0.; */
		R = 1.0;
		for( j = 1; j <= HIG/2; j++ ) {
			for( k = 1; k < DEPP; k++ ) {
				U1_[0][j][k] = R;
				U2_[0][j][k] = R * U;
				U3_[0][j][k] = 0.;/*R * V; */
				U4_[0][j][k] = 0.;/*R * W; */
				U5_[0][j][k] = P / K_1 + 0.5 * R * ( U * U /*+ V * V + W * W*/ );
			} /* end for() */
		} /* end for() */


		U = 200; /* V = 0.; W = 0.; */
		for( j = HIG/2+1; j < HIGG; j++ ) {
			for( k = 1; k < DEPP; k++ ) {
				U1_[0][j][k] = R;
				U2_[0][j][k] = R * U;
				U3_[0][j][k] = 0.;/*R * V; */
				U4_[0][j][k] = 0.;/*R * W; */
				U5_[0][j][k] = P / K_1 + 0.5 * R * ( U * U /*+ V * V + W * W*/ );
			} /* end for() */
		} /* end for() */


		/* right side - OUTFLOW */
		for( j = 1; j < HIGG; j++ ) {
			for( k = 1; k < DEPP; k++ ) {
				U1_[LENN][j][k] = R = U1_[LEN][j][k];
				U2_[LENN][j][k] = U = U2_[LEN][j][k];
				U3_[LENN][j][k] = V = U3_[LEN][j][k];
				U4_[LENN][j][k] = W = U4_[LEN][j][k];
				U5_[LENN][j][k] = P / K_1 + 0.5 * ( U * U + V * V + W * W ) / R;
			} /* end for() */
		} /* end for() */

		/*
		 * z - PERIODIC 
         * These is lateral direction.
		 */
		for( i = 0; i <= LENN; i++ ) {
			for( j = 1; j < HIGG; j++ ) {
				/* back side */
				U1_[i][j][0] = U1_[i][j][DEP];
				U2_[i][j][0] = U2_[i][j][DEP];
				U3_[i][j][0] = U3_[i][j][DEP];
				U4_[i][j][0] = U4_[i][j][DEP];
				U5_[i][j][0] = U5_[i][j][DEP];
				/* front side */
				U1_[i][j][DEPP] = U1_[i][j][1];
				U2_[i][j][DEPP] = U2_[i][j][1];
				U3_[i][j][DEPP] = U3_[i][j][1];
				U4_[i][j][DEPP] = U4_[i][j][1];
				U5_[i][j][DEPP] = U5_[i][j][1];
			} /* end for */
		} /* end for */

		/* 
		 * y - SLIP 
         * This is vertical direction
		 */
		for( i = 0; i <= LENN; i++ ) {
			for( k = 0; k <= DEPP; k++ ) {
				/* bottom side */
				U1_[i][0][k] =   U1_[i][1][k];
				U2_[i][0][k] =   U2_[i][1][k];
				U3_[i][0][k] = - U3_[i][1][k];
				U4_[i][0][k] =   U4_[i][1][k];
				U5_[i][0][k] =   U5_[i][1][k];
				/* top side */
				U1_[i][HIGG][k] =   U1_[i][HIG][k];
				U2_[i][HIGG][k] =   U2_[i][HIG][k];
				U3_[i][HIGG][k] = - U3_[i][HIG][k];
				U4_[i][HIGG][k] =   U4_[i][HIG][k];
				U5_[i][HIGG][k] =   U5_[i][HIG][k];
			} /* end for */
		}  /* end for */

} /* end BounCondInGhostCells() */
