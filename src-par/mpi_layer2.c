/* layer2.c ***********************************\
*  -- Flow in Mixing Layer modelled by LES --  *
*    explicit two/three-stage time differencing*
*  written by Andrei Chernousov  Dec 23, 2000  *
*        E-mail: <andrei99@iname.com>          *
*  http://www.geocities.com/andrei_chernousov  *
*  Modified by Nikola Mirkov 2018-present      *
*        E-mail: largeddysimulation@gmail.com  *
\**********************************************/

#include "mpi.h"

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
	  /* SGS viscosity */
	***mu_SGS,

	R, U, V, W, P, C, /* vector of primitive flow parameters */
	u1, u2, u3, u4, u5; /* conservative flow variables */

real 
	deltaX, deltaY, deltaZ, /* cell spacings */
	deltaT, /* time step */
	deltaT_X, deltaT_Y, deltaT_Z, /* convenient ratios */
	_deltaX, _deltaY, _deltaZ,
	_2deltaX, _2deltaY, _2deltaZ,
	_4deltaX, _4deltaY, _4deltaZ,
	// array of probes
	*probes;

real/* transport coefficients */
	/* molecular */
	mu_L, lambda_L, Pr_L,
	/* subgrid-scale */
	mu_T, Pr_T, cp_Pr_T,
	/* effective */
	mu_E, lambda_E,
	/* Smagorinsky constant, complex CsDD = Cs * delta * delta, where delta = (dX*dY*dZ)~0.33333 */
	Cs, CsDD, DD;

real maxCoNum;  /* maximum Courant number */
int nStages;    /* number of stages of Runge-Kutta algorithm */

/* MPI Buffer */
// real *Buf;
// unsigned BufCountU, BufCountF; // BufCountU > BufCountF

float *Bufp, *Bufn, *pBuf, *nBuf;
int BufCountU, BufCountF;


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

//--- handles of "probes.myid" files
FILE *pFprobes;

//-- for randomized disturbing blocks	
int 
    j_base, j_base_max = 6,
    j_period, j_period_min = 5, j_period_maxmin = 8;

/*--- for the specification of disturbed velocity block ---*/
int BL_HIG,           /* height of the disturbing block            */
	    k_min, k_max; /* its "minimal"-"maximal" indices      */

real Ud, Vd, Wd,     /* curent disturbed velocities           */
	 Ua, Va, Wa;     /* their amplitude values                */
int Ns_min, Nst;     /* minimal period of existence of this disturbing block and its max-min */

int continFlag;
double X;
	 


void termHandler( int sigType )
{
   if( continFlag-- ); /* 1 -> 0 */
   puts( "User break, terminating..." );
}


/*********
*  MAIN  *   Main function
*********/
int main(int argc, char *argv[])
{
int numprocs; // number of processes
int myid;     // identifier of _this_ process

// char processor_name[MPI_MAX_PROCESSOR_NAME];
char a;

int counter = 0;
float totalTime = 0.;

continFlag = 1; /* continuation flag (to be changed by SIGINT handler) */

MPI_File fh;
MPI_Status status;


    //-- MPI Initialization block
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);

	/*--- Get the seed of random number ---*/
    X = 0.317391239817; // hardcoded instead of reading it from file.

    //-- Initializations for numerical scheme
    Initialize(myid);

    /*--- time stepping ---*/
    while (1) {

    	//-- Broadcast step from the process with the rank "root" to others
    	MPI_Bcast(&step, 1, MPI_INT, 0, MPI_COMM_WORLD);

	    // Only root process checks single-byte record in "stopfile": 'g' or 's'
		if(0 == myid) {
		    MPI_File_open(MPI_COMM_SELF, "stopfile",
				    MPI_MODE_CREATE | MPI_MODE_RDWR, MPI_INFO_NULL, &fh);
		    MPI_File_set_view(fh, 0, MPI_CHAR, MPI_CHAR, "native", MPI_INFO_NULL);
		    MPI_File_read(fh, &a, 1, MPI_CHAR, &status);
		    MPI_File_close(&fh);
		    //
		    if (a == 's') {
			a = 'g';
			continFlag = 0;
			MPI_File_open(MPI_COMM_SELF, "stopfile",
					MPI_MODE_CREATE | MPI_MODE_RDWR,
					MPI_INFO_NULL, &fh);
					MPI_File_set_view(fh, 0, MPI_CHAR, MPI_CHAR,
					"native", MPI_INFO_NULL);
			MPI_File_write(fh, &a, 1, MPI_CHAR, &status);
			MPI_File_close(&fh);
		    }
		}

		//-- Broadcast continFlag from the process with the rank "root" to others
		MPI_Bcast(&continFlag, 1, MPI_INT, 0, MPI_COMM_WORLD);

		//-- Every process may check exit conditions
		if(step == numstep || !continFlag) break;

	    /* Check Courant number at every timestep */
		checkCoNum( myid );

		//-- root process prints current step and defines disturbances
		if(0 == myid) {
	        fprintf(stdout, "Timestep no.: %d, total time: %f sec.\n", ++step, (totalTime+=deltaT));
		    if (counter > 0) counter--;
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
		}

		/*=== Go trough stages ===*/
		for( Stage = 1; Stage <= nStages; Stage++){

			/*--- BC in ghost cells ---*/
			BounCondInGhostCells( myid, numprocs );
			/*--- Parameters at the cell boundaries ---*/
			Reconstruction( );
			/*--- BC at cell boundaries at the domain boundary ---*/
			BounCondOnInterfaces( myid, numprocs );
			/*--- Fluxes ---*/
			Fluxes( myid, numprocs );
			/*--- Evolution ---*/
			Evolution( nStages, Stage );

		}

		//-- Probes output
		Probes(myid);

		//-- Take frame
		if (step%f_step == 0 && step!=0) Output(myid);

    } // end while(1)

	/*--- Output the flowfield ---*/
	Output();

    //-- Finalize: backup the solution and close probes.* files
    Finalize(myid);

    //---
    fprintf(stdout, "%d/%d process stopped!\n", myid+1, numprocs);

    //---
    MPI_Finalize();
    
    //---
    return 0;

} // end main()

//-- end of mpi_layer2.c ---