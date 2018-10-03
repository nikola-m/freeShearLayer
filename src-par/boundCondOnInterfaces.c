#include "mpi.h"
#include "type.h"
#include "def.h"      /* Definitions, parameters */
#include "global.h"   /* global variables */

#include "bounCondOnInterfaces.h"

/***********************
* BOUNCONDONINTERFACES *  complete the "outer" parameters' values at the domain boundary
***********************/
void BounCondOnInterfaces(
	int myid, // identifier of _this_ process
	int numprocs )
{
int next, previous; // processes next/previous to this
MPI_Status status;

	//--- Communicate via MPI to complete outer parameters at cell x-bondaries

	//-- sending to next processes

	// next process's ID
	next = myid+1;

		// root process
	if(0 == myid) {
		// sending to the next
		for (i = 0, j = 0; j < HIG; j++) {
			for (k = 0; k < DEP; k++, i += 5) {
				Buf[i  ] = xU1[LEN][j][k];
				Buf[i+1] = xU2[LEN][j][k];
				Buf[i+2] = xU3[LEN][j][k];
				Buf[i+3] = xU4[LEN][j][k];
				Buf[i+4] = xU5[LEN][j][k];
			}
		}
		MPI_Send(Buf, BufCountF, MPI_FLOAT, next, 99, MPI_COMM_WORLD);
		// previous doesn't exist
/* 	MPI_Recv(Buf, BufCountF, MPI_FLOAT, MPI_ANY_SOURCE, 99, MPI_COMM_WORLD,
					&status);
		for (i = 0, j = 0; j < HIG; j++) {
			for (k = 0; k < DEP; k++, i += 5) {
				xU1[0][j][k] = Buf[i  ];
				xU2[0][j][k] = Buf[i+1];
				xU3[0][j][k] = Buf[i+2];
				xU4[0][j][k] = Buf[i+3];
				xU5[0][j][k] = Buf[i+4];
			}
		}
*/
	}

	// other processes
	else {
		// receiving from the previous, everyone
		MPI_Recv(Buf, BufCountF, MPI_FLOAT, MPI_ANY_SOURCE, 99, MPI_COMM_WORLD,
					&status);
		for (i = 0, j = 0; j < HIG; j++) {
			for (k = 0; k < DEP; k++, i += 5) {
				xU1[0][j][k] = Buf[i  ];
				xU2[0][j][k] = Buf[i+1];
				xU3[0][j][k] = Buf[i+2];
				xU4[0][j][k] = Buf[i+3];
				xU5[0][j][k] = Buf[i+4];
			}
		}
		// sending to the next; if not the last!
		if(myid != numprocs-1) {
			for (i = 0, j = 0; j < HIG; j++) {
				for (k = 0; k < DEP; k++, i += 5) {
					Buf[i  ] = xU1[LEN][j][k];
					Buf[i+1] = xU2[LEN][j][k];
					Buf[i+2] = xU3[LEN][j][k];
					Buf[i+3] = xU4[LEN][j][k];
					Buf[i+4] = xU5[LEN][j][k];
				}
			}
			MPI_Send(Buf, BufCountF, MPI_FLOAT, next, 99, MPI_COMM_WORLD);
		} // end if(myid != numprocs-1)
	}

	//-- sending to previous processes

	// previous process's ID
	previous = myid-1;

	// root process
	if(0 == myid) {
		// previous doesn't exist
/* 	for (i = 0, j = 0; j < HIG; j++) {
			for (k = 0; k < DEP; k++, i += 5) {
				Buf[i  ] = U1x[0][j][k];
				Buf[i+1] = U2x[0][j][k];
				Buf[i+2] = U3x[0][j][k];
				Buf[i+3] = U4x[0][j][k];
				Buf[i+4] = U5x[0][j][k];
			}
		}
		MPI_Send(Buf, BufCountF, MPI_FLOAT, previous, 99, MPI_COMM_WORLD);
*/
		// receiving from the next
		MPI_Recv(Buf, BufCountF, MPI_FLOAT, MPI_ANY_SOURCE, 99, MPI_COMM_WORLD,
					&status);
		for (i = 0, j = 0; j < HIG; j++) {
			for (k = 0; k < DEP; k++, i += 5) {
				U1x[LEN][j][k] = Buf[i  ];
				U2x[LEN][j][k] = Buf[i+1];
				U3x[LEN][j][k] = Buf[i+2];
				U4x[LEN][j][k] = Buf[i+3];
				U5x[LEN][j][k] = Buf[i+4];
			}
		}
	}

	// other processes
	else {
		// receiving from the next; if not the last process!
		if(myid != numprocs-1) {
			MPI_Recv(Buf, BufCountF, MPI_FLOAT, MPI_ANY_SOURCE, 99, MPI_COMM_WORLD,
						&status);
			for (i = 0, j = 0; j < HIG; j++) {
				for (k = 0; k < DEP; k++, i += 5) {
					U1x[LEN][j][k] = Buf[i  ];
					U2x[LEN][j][k] = Buf[i+1];
					U3x[LEN][j][k] = Buf[i+2];
					U4x[LEN][j][k] = Buf[i+3];
					U5x[LEN][j][k] = Buf[i+4];
				}
			}
		} // wnd if(myid != numprocs-1)
		// sending to previous - everyone
		for (i = 0, j = 0; j < HIG; j++) {
			for (k = 0; k < DEP; k++, i += 5) {
				Buf[i  ] = U1x[0][j][k];
				Buf[i+1] = U2x[0][j][k];
				Buf[i+2] = U3x[0][j][k];
				Buf[i+3] = U4x[0][j][k];
				Buf[i+4] = U5x[0][j][k];
			}
		}
		MPI_Send(Buf, BufCountF, MPI_FLOAT, previous, 99, MPI_COMM_WORLD);
	}

	// syncronization
	MPI_Barrier(MPI_COMM_WORLD);


	/*--- x ---*/
	/*- left side */
	if(0 == myid) { /* root process */
		P = 100000.;
		R = 1.0;
		for (j = 0; j < HIG; j++) {
			for( k = 0; k < DEP; k++ ) {
				//-- choosing components of velocity 
                                // lower half				
				if (j < HIG/2){
				    // always
				    U = 76.4;
				    // lower than disturbing block
				    if (j < (HIG/2 - BL_HIG)) {
					V = 0.;
					W = 0.;				    
				    }
				    // within disturbing block
				    else {
					if ((((j - 1- j_base)/j_period)%2) > 0) {
					    V = Vd; 
					    W = Wd;
					}
					else {
					    V = 0.;
					    W = 0.;
					}
				    }
				}
				// upper half
				else {
				    U = 200.;
				    V = 0.;
				    W = 0.;
				}
				//
				xU1[0][j][k] = R;
				xU2[0][j][k] = R * U;
				xU3[0][j][k] = R * V;
				xU4[0][j][k] = R * W;
				xU5[0][j][k] = P / K_1 + 0.5 * R * ( U * U + V * V + W * W);
			} /* end for() */
		} /* end for() */
	} // end if(0 == myid) {

	if(myid == numprocs-1) { /* last process */
		/* right side - OUTFLOW */
		P = 100000.;
		for (j = 0; j < HIG; j++) {
			for (k = 0; k < DEP; k++) {
				U1x[LEN][j][k] = R = xU1[LEN][j][k];
				U2x[LEN][j][k] = U = xU2[LEN][j][k];
				U3x[LEN][j][k] = V = xU3[LEN][j][k];
				U4x[LEN][j][k] = W = xU4[LEN][j][k];
				U5x[LEN][j][k] = P / K_1 + 0.5 * (U*U + V*V + W*W ) / R;
			} /* end for() */
		} /* end for() */
	} // end if(myid == numprocs-1)

	/*--- z: PERIODIC ---*/
	for (i = 0; i < LEN; i++) {
		for (j = 0; j < HIG; j++) {
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

	/*--- y: SLIP ---*/
	for (i = 0; i < LEN; i++) {
		for (k = 0; k < DEP; k++) {
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
} /* end BounCondOnInterfaces( void ) */