#include "mpi.h"
#include <stdio.h>     /* printf() etc.*/
#include <stdlib.h>    /* atoi()       */
#include <math.h>      /* pow(), sqrt()   */

#include "type.h"
#include "def.h"      /* Definitions, parameters */
#include "global.h"   /* global variables */
#include "helpers.h"  /* helper functions */
#include "initialize.h"

/***************
*  INITIALIZE  *    all necessary initializstions
***************/
void Initialize(int myid) // identifier of _this_ process
{
  // temporary arrays
float ***fp[10];
char filename[30];

unsigned i;
MPI_File fh;
MPI_Status status;
float buf;

	//--- root process (node w/rank 0) reads data from file "mpi_layer2.bin"
	//    and broadcasts it to others
	if(0 == myid) {
		// open the file...
		MPI_File_open(MPI_COMM_SELF, "mpi_layer2.bin", MPI_MODE_CREATE | MPI_MODE_RDWR, MPI_INFO_NULL, &fh);
		MPI_File_set_view(fh, 0, MPI_FLOAT, MPI_FLOAT, "native", MPI_INFO_NULL);
	}

	//--- ... and reads/broadcast data (if myid == 0), otherwise receives it
	i = 0;
	do{

		// syncronization
		MPI_Barrier(MPI_COMM_WORLD);

		// root reads the next data item
		if (0 == myid)
		    MPI_File_read(fh, &buf, 1, MPI_FLOAT, &status);

		// root process broadcasts the next data item
		MPI_Bcast(&buf, 1, MPI_FLOAT, 0, MPI_COMM_WORLD);

		// For debugging. May cause a lot of lines in stdout for larger no. of processes...
		// fprintf(stdout, "%d process: after broadcast %dth data item = %f\n",
		// 			myid+1, i, buf);

		// dispatch to specific parameters
		switch (i) {
			case  0: fprintf(stdout, "LEN    = %d\n", LEN    = buf); break;
			case  1: fprintf(stdout, "HIG    = %d\n", HIG    = buf); break;
			case  2: fprintf(stdout, "DEP    = %d\n", DEP    = buf); break;
			case  3: fprintf(stdout, "deltaX = %f\n", deltaX = buf); break;
			case  4: fprintf(stdout, "deltaY = %f\n", deltaY = buf); break;
			case  5: fprintf(stdout, "deltaZ = %f\n", deltaZ = buf); break;
			case  6: fprintf(stdout, "deltaT = %f\n", deltaT = buf); break;
			case  7: fprintf(stdout, "Answer = %d\n", Answer = buf); break;
			case  8: fprintf(stdout, "numstep= %d\n", numstep= buf); break;
			case  9: fprintf(stdout, "mu_L   = %f\n", mu_L   = buf); break;
			case 10: fprintf(stdout, "Pr_L   = %f\n", Pr_L   = buf); break;
			case 11: fprintf(stdout, "Pr_T   = %f\n", Pr_T   = buf); break;
			case 12: fprintf(stdout, "Cs     = %f\n", Cs     = buf); break;
			case 13: fprintf(stdout, "bl_hig = %d\n", BL_HIG = buf); break;
			case 14: fprintf(stdout, "Ua     = %f\n", Ua     = buf); break;						
			case 15: fprintf(stdout, "Va     = %f\n", Va     = buf); break;
			case 16: fprintf(stdout, "Wa     = %f\n", Wa     = buf); break;
			case 17: fprintf(stdout, "Ns_min = %d\n", Ns_min = buf); break;
			case 18: fprintf(stdout, "Nst    = %d\n", Nst    = buf); break;
			case 19: fprintf(stdout, "f_step = %d\n", f_step = buf); break;
			case 20: fprintf(stdout, "nStages = %d\n", nStages = buf); break;
			case 21: fprintf(stdout, "maxCoNum = %f\n", maxCoNum = buf); break;	
			default: break;
		} // end switch

		i++;

	} while (i < 22);

	//--- close file
	if (0 == myid)
	    MPI_File_close(&fh);

	//--- other initializatons
	// complexes with deltas
	deltaT_X = deltaT / deltaX;
	deltaT_Y = deltaT / deltaY;
	deltaT_Z = deltaT / deltaZ;
	 _deltaX = 1. / deltaX;
	 _deltaY = 1. / deltaY;
	 _deltaZ = 1. / deltaZ;
	_4deltaX = 1. / ( 4 * deltaX );
	_4deltaY = 1. / ( 4 * deltaY );
	_4deltaZ = 1. / ( 4 * deltaZ );
	_2deltaX = _4deltaX + _4deltaX;
	_2deltaY = _4deltaY + _4deltaY;
	_2deltaZ = _4deltaZ + _4deltaZ;
	// box sizes + 1
	LENN = LEN + 1;
	HIGG = HIG + 1;
	DEPP = DEP + 1;
	// constant-value lambda_L
	lambda_L = mu_L * Cv * K / Pr_L;
	// complex cp_Pr_T
	cp_Pr_T = Cv * K / Pr_T;
	/* complex Cd * delta * delta in Smagorinsky model */
	DD = pow( deltaX * deltaY * deltaZ, twoThirds ); // Filter width squared.
	CsDD = Cs * DD;

	//--- dynamic allocation of memory
	// total amount of memory required, for gas dynamics 5
	mem = 11 * (LEN + 2) * (HIG + 2) * (DEP + 2)
	    + 10 * (LEN + 1) *  HIG	 *  DEP
	    + 10 *  LEN      * (HIG + 1) *  DEP
	    + 10 *  LEN      *  HIG	 * (DEP + 1);
	fprintf(stdout, "%d process: %lu bytes of memory required\n",
				myid+1, mem*sizeof(float) );
	// for quantities in cells
	for (i = 0; i < 10; i++) 
	    fp[i] = Array3D(LEN+2, HIG+2, DEP+2);
	U1 = fp[0]; U1p= fp[5];
	U2 = fp[1]; U2p= fp[6];
	U3 = fp[2]; U3p= fp[7];
	U4 = fp[3]; U4p= fp[8];
	U5 = fp[4]; U5p= fp[9];
	// for x fluxes
	for (i = 0; i < 10; i++) 
	    fp[i] = Array3D(LENN, HIG, DEP);
	xU1 = fp[0]; U1x = fp[5];
	xU2 = fp[1]; U2x = fp[6];
	xU3 = fp[2]; U3x = fp[7];
	xU4 = fp[3]; U4x = fp[8];
	xU5 = fp[4]; U5x = fp[9];
	// for y fluxes
	for (i = 0; i < 10; i++) 
	    fp[i] = Array3D(LEN, HIGG, DEP);
	yU1 = fp[0]; U1y = fp[5];
	yU2 = fp[1]; U2y = fp[6];
	yU3 = fp[2]; U3y = fp[7];
	yU4 = fp[3]; U4y = fp[8];
	yU5 = fp[4]; U5y = fp[9];
	// for z fluxes
	for (i = 0; i < 10; i++) 
	    fp[i] = Array3D( LEN, HIG, DEPP );
	zU1 = fp[0]; U1z = fp[5];
	zU2 = fp[1]; U2z = fp[6];
	zU3 = fp[2]; U3z = fp[7];
	zU4 = fp[3]; U4z = fp[8];
	zU5 = fp[4]; U5z = fp[9];
	/* For SGS viscosity */
	mu_SGS = Array3D( LEN+2, HIG+2, DEP+2 );

	fprintf(stdout, "%d process: 3D arrays allocated\n", myid+1);

		// MPI Buf
	BufCountF = 5 *  HIG    *  DEP;    // BufCountF < BufCountU
	BufCountU = 5 * (HIG+2) * (DEP+2); // BufCountU
	if((Buf = (float *)malloc(BufCountU * sizeof(float))) == NULL) {
		fprintf(stderr, "mpi_layer2: can't allocate memory.\n");
		MPI_Abort(MPI_COMM_WORLD, 1);
	}

		/* array of probes */
	if((probes=(float *)malloc(sizeof(float) * HIG)) == NULL) {
		fprintf(stderr, "mpi_layer2: can't allocate memory for probes.\n");
		MPI_Abort(MPI_COMM_WORLD, 1);
	}

	//--- Initial conditions
	// start new simulation
	if(0 == Answer) {
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
	}

	// continue previously saved simulation
	else {
		//
		sprintf(filename, "backup.%d", myid);
		MPI_File_open(MPI_COMM_SELF, filename, MPI_MODE_CREATE | MPI_MODE_RDWR,
							MPI_INFO_NULL, &fh);
		MPI_File_set_view(fh, 0, MPI_FLOAT, MPI_FLOAT, "native", MPI_INFO_NULL);
		//
		for (i = 1; i < LENN; i++) {
			for (j = 1; j < HIGG; j++) {
				MPI_File_read(fh, &U1[i][j][1], DEP, MPI_FLOAT, &status);
				MPI_File_read(fh, &U2[i][j][1], DEP, MPI_FLOAT, &status);
				MPI_File_read(fh, &U3[i][j][1], DEP, MPI_FLOAT, &status);
				MPI_File_read(fh, &U4[i][j][1], DEP, MPI_FLOAT, &status);
				MPI_File_read(fh, &U5[i][j][1], DEP, MPI_FLOAT, &status);
			}
		}
	    	// restore current step		
		MPI_File_read(fh, &buf, 1, MPI_FLOAT, &status);
		step = (int)buf;
	    	// restore and set the seed of pseudo-random number
		MPI_File_read(fh, &buf, 1, MPI_FLOAT, &status);		
	    	srand((unsigned int)buf);
		// close bakup file
		MPI_File_close(&fh);
	}

	// create/open "probes.muid" file in binary mode
	sprintf(filename, "probes.%d", myid);
	if ((pFprobes = fopen(filename, "ab")) == NULL)
		fprintf(stdout, "can't open \"%s\".\n", filename);

	// first stage
	Stage = 1;
	U1_ = U1;
	U2_ = U2;
	U3_ = U3;
	U4_ = U4;
	U5_ = U5;
	fprintf(stdout, "%d process: after initialization.\n", myid+1);

} // end Initialize()