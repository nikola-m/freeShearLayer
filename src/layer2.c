/* layer2.c ***********************************\
*  -- Flow in Mixing Layer modeled by LES --   *
*    explicit two-stage time differencing      *
*  written by Andrei Chernousov  Dec 23, 2000  *
*        E-mail: <andrei99@iname.com>          *
*  http://www.geocities.com/andrei_chernousov  *
*  Modified by Nikola Mirkov 2018              *
*        E-mail: largeddysimulation@gmail.com  *
\**********************************************/

#include <stdio.h>     /* printf() etc.*/
#include <math.h>      /* sqrt()       */
#include <stdlib.h>    /* atoi()       */
#include <signal.h>    /* sigint function to stop the program execution */

#include "type.h"
#include "def.h"      /* Definitions, parameters */
#include "global.h"   /* Global variables */
#include "helpers.h"  /* Helper functions */ 

#include "initialize.h"
#include "bounCondOnInterfaces.h"
#include "bounCondInGhostCells.h"
#include "reconstruction.h"
#include "fluxes.h"
#include "evolution.h"
#include "output.h"
#include "probes.h"
#include "finalize.h"
#include "turbulence.h"


/* Definition of global variables */

real
	***U1_, ***U2_, ***U3_, ***U4_, ***U5_, /* substitutional pointers */
	  /* Cons. variables at the cells' centroids */
	***U1,  ***U2,  ***U3,  ***U4,  ***U5,
	***U1p, ***U2p, ***U3p, ***U4p, ***U5p,
	  /* Cons. variables at cell boundaries in x-direction */
	***xU1, ***xU2, ***xU3, ***xU4, ***xU5,
	***U1x, ***U2x, ***U3x, ***U4x, ***U5x,
	  /* Cons. variables at cell boundaries in y-direction */
	***yU1, ***yU2, ***yU3, ***yU4, ***yU5,
	***U1y, ***U2y, ***U3y, ***U4y, ***U5y,
	  /* Cons. variables at cell boundaries in y-direction */
	***zU1, ***zU2, ***zU3, ***zU4, ***zU5,
	***U1z, ***U2z, ***U3z, ***U4z, ***U5z,

	R, U, V, W, P, C, /* vector of primitive flow parameters */
	u1, u2, u3, u4, u5; /* conservative flow variables */

real 
	deltaX, deltaY, deltaZ, /* cell spacings */
	deltaT, /* time step */
	deltaT_X, deltaT_Y, deltaT_Z, /* convenient ratios */
	_deltaX, _deltaY, _deltaZ,
	_2deltaX, _2deltaY, _2deltaZ,
	_4deltaX, _4deltaY, _4deltaZ;


real/* transport coefficients */
	/* molecular */
	mu_L, lambda_L, Pr_L,
	/* subgrid-scale */
	mu_T, Pr_T, cp_Pr_T,
	/* effective */
	mu_E, lambda_E,
	/* Smagorinsky constant, complex CsDD = Cs * delta * delta, where delta = (dX*dY*dZ)~0.33333 */
	Cs, CsDD;

real maxCoNum;  /* maximum Courant number */
int nStages;    /* number of stages of Runge-Kutta algorithm */

unsigned
   step,    /* current time step */
   numstep, /* overall prescribed number of time steps */
   mem,     /* total amount of dynamic memory required (in bytes) */
   f_step,  /* frame taking step */
   
   Answer, /* solution continuation flag      */
   LEN, HIG, DEP,  /* cell numbers in x-, y-, z-directions, respectively for 3d box */
   LENN, HIGG, DEPP,
   i,  j,  k, i__,  j__,  k__, /* indices */
   l;
   
char
   Stage,    /* indicator of the current stage */
   str[80];  /* buffer string */
  
real *probes; /*--- Array of probes ---*/

/*--- for the specification of disturbed velocity block ---*/
unsigned BL_HIG,  /* height of the disturbing block            */
	  k_min, k_max;      /* its "minimal"-"maximal" indices      */
real Ud, Vd, Wd, /* curent disturbed velocities           */
	  Ua, Va, Wa;     /* their amplitude values                */
int Ns_min, Nst;  /* minimal period of existence of this disturbing block and its max-min */
int continFlag;
double X;
	 

/* SIGINT handler */
// #pragma argsused

void termHandler( int sigType )
{
   if( continFlag-- ); /* 1 -> 0 */
   puts( "User break, terminating..." );
}


/*********
*  MAIN  *   Main function
*********/
int main( void )
{

	int counter = 0;
    float totalTime = 0.;

	continFlag = 1; /* continuation flag (to be changed by SIGINT handler) */


	/*--- A handler for SIGINT for this program ---*/
	signal( SIGINT, termHandler );
	
	/*--- Get the seed of random number ---*/
    if (scanf("%lf", &X) == 1) {
        printf( "Random number [0...1]: " );
        printf( "X = %f\n", X );
    } else {
        printf("Failed to read random number.\n");
    }

	/*--- Initialization ---*/
	Initialize();

	/*--- Time stepping ---*/
	while( 1 ) {

		/*--- Check exit conditions ---*/
		if( step == numstep || !continFlag  ) break;

        /* Check Courant number at every timestep */
		checkCoNum();

		printf( "Timestep no.: %d, total time: %f sec.\n", ++step, (totalTime+=deltaT) );

		/*--- Define disturbing block */
		if( counter > 0 ) counter--;
		else {
			counter = (int)( Ns_min + Nst * ( X = Random( X ) ) );
			k_min   = (unsigned)( DEP * ( X = Random( X ) ) );
			k_max   = (unsigned)( k_min + DEP * ( X = Random( X ) ) );
			if( k_max > DEP ) k_max -= DEP;
			X = Random( X ); Ud = Ua * ( 1 - X - X );
			X = Random( X ); Vd = Va * ( 1 - X - X );
			X = Random( X ); Wd = Wa * ( 1 - X - X );
			printf( "k_min = %d, k_max = %d, Ud = %f, Vd = %f, Wd = %f\n", k_min, k_max, Ud, Vd, Wd );
		}

		/*=== Go trough stages ===*/
		for( Stage = 1; Stage <= nStages; Stage++){

			/*--- BC in ghost cells ---*/
			BounCondInGhostCells( );
			/*--- Parameters at the cell boundaries ---*/
			Reconstruction( );
			/*--- BC at cell boundaries at the domain boundary ---*/
			BounCondOnInterfaces( );
			/*--- Fluxes ---*/
			Fluxes( );
			/*--- Evolution ---*/
			Evolution( nStages, Stage );

		}

		/*--- Probes output ---*/
		Probes();
		/*--- Taking frame ---*/
		if( step % f_step == 0 ) Output();

	} /* end while */

	/*--- Output the flowfield ---*/
	Output();
	/*--- Backup solution---*/
	Finalize();

	puts( "Simulation finished!" );
	return 0;

} /* end main() */

/* end of layer2.c */