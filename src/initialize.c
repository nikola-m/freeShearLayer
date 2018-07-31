#include <stdio.h>     /* printf() etc.*/
#include <stdlib.h>    /* atoi()       */
#include <math.h>      /* pow(), sqrt()   */

#include "type.h"
#include "def.h"      /* Definitions, parameters */
#include "global.h"   /* global variables */
#include "helpers.h"  /* helper functions */
#include "initialize.h"

/***************
*  INITIALIZE  *    Performs necessary initializstions
***************/
void Initialize( void )
{

	/*--- Loading data from "layer2.ini" file ---*/
	FILE *pF;
	if( ( pF = fopen( "layer2.ini", "r" ) ) == NULL ) {
	   printf( "Cannot open \"layer2.ini\" file");
	   exit( -1 );
	}

	int i = 0;
	do{
		if ( fgets( str, 79, pF ) != NULL ) {
			switch( i ) {
				case  0:  printf( "LEN    = %d\n",  LEN    = atoi( str ) ); break;
				case  1:  printf( "HIG    = %d\n",  HIG    = atoi( str ) ); break;
				case  2:  printf( "DEP    = %d\n",  DEP    = atoi( str ) ); break;
				case  3:  printf( "deltaX = %f\n",  deltaX = atof( str ) ); break;
				case  4:  printf( "deltaY = %f\n",  deltaY = atof( str ) ); break;
				case  5:  printf( "deltaZ = %f\n",  deltaZ = atof( str ) ); break;
				case  6:  printf( "deltaT = %f\n",  deltaT = atof( str ) ); break;
				case  7:  printf( "Answer = %d\n",  Answer = atoi( str ) ); break;
				case  8:  printf( "numstep= %d\n",  numstep= atoi( str ) ); break;
				case  9:  printf( "mu_L   = %f\n",  mu_L   = atof( str ) ); break;
				case 10:  printf( "Pr_L   = %f\n",  Pr_L   = atof( str ) ); break;
				case 11:  printf( "Pr_T   = %f\n",  Pr_T   = atof( str ) ); break;
				case 12:  printf( "Cs     = %f\n",  Cs     = atof( str ) ); break;
				case 13:  printf( "BL_HIG = %d\n",  BL_HIG = atoi( str ) ); break;
				case 14:  printf( "Ua     = %f\n",  Ua     = atof( str ) ); break;
				case 15:  printf( "Va     = %f\n",  Va     = atof( str ) ); break;
				case 16:  printf( "Wa     = %f\n",  Wa     = atof( str ) ); break;
				case 17:  printf( "Ns_min = %d\n",  Ns_min = atoi( str ) ); break;
				case 18:  printf( "Nst    = %d\n",  Nst    = atoi( str ) ); break;
				case 19:  printf( "f_step = %d\n",  f_step = atoi( str ) ); break;
				case 20:  printf( "nStages = %d\n", nStages  = atoi( str ) ); break;
				case 21:  printf( "maxCoNum = %f\n",maxCoNum = atof( str ) ); break;
				default: break;
			} /* end switch */
			i++;
		} /* end if() */
		else {
		   puts( "Error while reading \"layer2.ini\" file" );
		   exit( -1 );
		}
	} while( i < 22 );

	fclose( pF );

	printf("Total number of cells in computational domain: %d\n", LEN*HIG*DEP );

	/*--- other initializatons ---*/
		/* complexes with deltas */
	deltaT_X = deltaT / deltaX;     deltaT_Y = deltaT / deltaY;     deltaT_Z = deltaT / deltaZ;
	 _deltaX = 1. / deltaX;         _deltaY  = 1. / deltaY;          _deltaZ = 1. / deltaZ;
	_4deltaX = 1. / ( 4 * deltaX ); _4deltaY = 1. / ( 4 * deltaY ); _4deltaZ = 1. / ( 4 * deltaZ );
	_2deltaX = _4deltaX + _4deltaX; _2deltaY = _4deltaY + _4deltaY; _2deltaZ = _4deltaZ + _4deltaZ;
		/* sizes + 1 */
	LENN = LEN + 1; HIGG = HIG + 1; DEPP = DEP + 1;
		/* constant-value lambda_L */
	lambda_L = mu_L * Cv * K / Pr_L;
		/* complex cp_Pr_T */
	cp_Pr_T = Cv * K / Pr_T;
		/* complex Cd * delta * delta in Smagorinsky model */
	CsDD = Cs * pow( deltaX * deltaY * deltaZ, 0.66666667 );
	/*--- Dynamic memory allocation ---*/
	{
	   /* temporal array of pointers */
	real ***fp[10];
		/* total amount of memory required, for gas dynamics 5 */
	mem = 10 * ( (unsigned long)LEN + 2 ) * ( (unsigned long)HIG + 2 ) * ( (unsigned long)DEP + 2 )
		+ 10 * ( (unsigned long)LEN + 1 ) *   (unsigned long)HIG	   *   (unsigned long)DEP
		+ 10 *   (unsigned long)LEN       * ( (unsigned long)HIG + 1 ) *   (unsigned long)DEP
		+ 10 *   (unsigned long)LEN	      *   (unsigned long)HIG	   * ( (unsigned long)DEP + 1 );
	printf( "%lu bytes of dynamic memory required...", mem * sizeof( real ) );
		/* for quantities in cells */
	for( i = 0; i < 10; i++ ) fp[i] = Array3D( LEN+2, HIG+2, DEP+2 );
	U1 = fp[0]; U2 = fp[1]; U3 = fp[2]; U4 = fp[3]; U5 = fp[4];
	U1p= fp[5]; U2p= fp[6]; U3p= fp[7]; U4p= fp[8]; U5p= fp[9];
		/* for x fluxes */
	for( i = 0; i < 10; i++ ) fp[i] = Array3D( LENN, HIG, DEP );
	xU1 = fp[0]; xU2 = fp[1]; xU3 = fp[2]; xU4 = fp[3]; xU5 = fp[4];
	U1x = fp[5]; U2x = fp[6]; U3x = fp[7]; U4x = fp[8]; U5x = fp[9];
		/* for y fluxes */
	for( i = 0; i < 10; i++ ) fp[i] = Array3D( LEN, HIGG, DEP );
	yU1 = fp[0]; yU2 = fp[1]; yU3 = fp[2]; yU4 = fp[3]; yU5 = fp[4];
	U1y = fp[5]; U2y = fp[6]; U3y = fp[7]; U4y = fp[8]; U5y = fp[9];
		/* for x fluxes */
	for( i = 0; i < 10; i++ ) fp[i] = Array3D( LEN, HIG, DEPP );
	zU1 = fp[0]; zU2 = fp[1]; zU3 = fp[2]; zU4 = fp[3]; zU5 = fp[4];
	U1z = fp[5]; U2z = fp[6]; U3z = fp[7]; U4z = fp[8]; U5z = fp[9];
	printf(" allocated!\n" );
	} /* end block */
		/* array of probes */
	if( (probes=(real*)malloc( sizeof(real)*HIG)) == NULL ) {
	   puts( "Cannot allocate memory for probes" );
	   exit( -1 );
	}
	/*--- Initial conditions ---*/
	printf( "Answer = %d\n", Answer );
	if( Answer ) { /* restoring from backup */
		if ( ( pF = fopen( "backup.dat", "rb") ) == NULL ) {
		   puts( "Cannot open \"backup.dat\" file to read" );
		   exit( -1 );
		}
		for ( i = 1; i < LENN; i++ ) {
			for ( j = 1; j < HIGG; j++ ) {
				if ( fread( &U1[i][j][1], DEP * sizeof( U1[i][j][1] ), 1, pF ) == 1 &&
					 fread( &U2[i][j][1], DEP * sizeof( U2[i][j][1] ), 1, pF ) == 1 &&
					 fread( &U3[i][j][1], DEP * sizeof( U3[i][j][1] ), 1, pF ) == 1 &&
				     fread( &U4[i][j][1], DEP * sizeof( U4[i][j][1] ), 1, pF ) == 1 &&
					 fread( &U5[i][j][1], DEP * sizeof( U5[i][j][1] ), 1, pF ) == 1 ) {

			      // puts( "Reading fields succesfull." );

				} else {

				  puts( "Could not read the fields" );

				}
			}
		}
		if ( fread( &step, sizeof( step ), 1, pF ) == 1 ){
			puts( "Reading step value successful." );
		}
		fclose( pF );
	} /* end if restoring from backup */
	else { /* start new simulation */
		P = 100000.;
		R = 1.0;
		U = 76.4; /*V = W = 0.;*/
		for ( i = 1; i < LENN; i++ ) {
			for ( j = 1; j <= HIG/2; j++ ) {
				for ( k = 1; k < DEPP; k++ ) {
					U1[i][j][k] = R;
					U2[i][j][k] = R * U;
					U3[i][j][k] = 0.;/*R * V;*/
					U4[i][j][k] = 0.;/*R * W;*/
					U5[i][j][k] = P / K_1 + 0.5 * R * ( U * U /*+ V * V + W * W*/ );
				} /* end for() */
			} /* end for() */
		} /* end for() */
		U = 200;
		for ( i = 1; i < LENN; i++ ) {
			for ( j = HIG/2+1; j < HIGG; j++ ) {
				for ( k = 1; k < DEPP; k++ ) {
					U1[i][j][k] = R;
					U2[i][j][k] = R * U;
					U3[i][j][k] = 0.;/*R * V*/;
					U4[i][j][k] = 0.;/*R * W*/;
					U5[i][j][k] = P / K_1 + 0.5 * R * ( U * U /*+ V * V + W * W*/ );
				} /* end for() */
			} /* end for() */
		} /* end for() */
		step = 0;  /* null step number at the beginning */
	} /* end else start new simulation */

	Stage = 1;
	U1_ = U1, U2_ = U2, U3_ = U3, U4_ = U4, U5_ = U5;

} /* end Initialize( ) */