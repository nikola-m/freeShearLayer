#include "mpi.h"
#include "type.h"
#include "def.h"      /* Definitions, parameters */
#include "global.h"   /* global variables */

#include "bounCondInGhostCells.h"

/***********************
* BOUNCONDINGHOSTCELLS *   BC in the cells out of the domain
***********************/
void BounCondInGhostCells(
	int myid,      // identifier of _this_ process
	int numprocs ) // number of processes
{
unsigned i, j, k;

int next, previous; // processes next/previous to this
MPI_Status status;

	//--- Communicate via MPI to complete cell x-layers ---

	//-- sending to next processes

	// next process's ID
	next = myid+1;

	// root process
	if (0 == myid) {
		// sending to the next
		for (i = 0, j = 0; j <= HIGG; j++) {
			for (k = 0; k <= DEPP; k++, i += 5) {
				Buf[i  ] = U1_[LEN][j][k];
				Buf[i+1] = U2_[LEN][j][k];
				Buf[i+2] = U3_[LEN][j][k];
				Buf[i+3] = U4_[LEN][j][k];
				Buf[i+4] = U5_[LEN][j][k];
			}
		}
		MPI_Send(Buf, BufCountU, MPI_FLOAT, next, 99, MPI_COMM_WORLD);
		// previous doesn't exist
/*		MPI_Recv(Buf, BufCountU, MPI_FLOAT, MPI_ANY_SOURCE, 99, MPI_COMM_WORLD,
					&status);
		for (i = 0, j = 0; j <= HIGG; j++) {
			for (k = 0; k <= DEPP; k++, i += 5) {
				U1_[0][j][k] = Buf[i  ];
				U2_[0][j][k] = Buf[i+1];
				U3_[0][j][k] = Buf[i+2];
				U4_[0][j][k] = Buf[i+3];
				U5_[0][j][k] = Buf[i+4];
			}
		}
*/
	} // end if()

	// other processes
	else {
		// receiving from the previous, everyone
		MPI_Recv(Buf, BufCountU, MPI_FLOAT, MPI_ANY_SOURCE, 99, MPI_COMM_WORLD,
					&status);
		for (i = 0, j = 0; j <= HIGG; j++) {
			for (k = 0; k <= DEPP; k++, i += 5) {
				U1_[0][j][k] = Buf[i  ];
				U2_[0][j][k] = Buf[i+1];
				U3_[0][j][k] = Buf[i+2];
				U4_[0][j][k] = Buf[i+3];
				U5_[0][j][k] = Buf[i+4];
			}
		}
		// sending to the next; if not the last!
		if(myid != numprocs-1) {
			for (i = 0, j = 0; j <= HIGG; j++) {
				for (k = 0; k <= DEPP; k++, i += 5) {
					Buf[i  ] = U1_[LEN][j][k];
					Buf[i+1] = U2_[LEN][j][k];
					Buf[i+2] = U3_[LEN][j][k];
					Buf[i+3] = U4_[LEN][j][k];
					Buf[i+4] = U5_[LEN][j][k];
				}
			}
			MPI_Send(Buf, BufCountU, MPI_FLOAT, next, 99, MPI_COMM_WORLD);
		} // end if(myid != numprocs-1)
	}

	//-- sending to previous processes

	// previous process's ID
	previous = myid-1;

	// root process
	if(0 == myid) {
		// previous doesn't exist
/*		for (i = 0, j = 0; j <= HIGG; j++) {
			for (k = 0; k <= DEPP; k++, i += 5) {
				Buf[i  ] = U1_[1][j][k];
				Buf[i+1] = U2_[1][j][k];
				Buf[i+2] = U3_[1][j][k];
				Buf[i+3] = U4_[1][j][k];
				Buf[i+4] = U5_[1][j][k];
			}
		}
		MPI_Send(Buf, BufCountU, MPI_FLOAT, previous, 99, MPI_COMM_WORLD);
*/
		// receiving from the next
		MPI_Recv(Buf, BufCountU, MPI_FLOAT, MPI_ANY_SOURCE, 99, MPI_COMM_WORLD,
					&status);
		for (i = 0, j = 0; j <= HIGG; j++) {
			for (k = 0; k <= DEPP; k++, i += 5) {
				U1_[LENN][j][k] = Buf[i  ];
				U2_[LENN][j][k] = Buf[i+1];
				U3_[LENN][j][k] = Buf[i+2];
				U4_[LENN][j][k] = Buf[i+3];
				U5_[LENN][j][k] = Buf[i+4];
			}
		}
	} // end if()

	// other processes
	else {
		// receiving from the next; if not the last process!
		if(myid != numprocs-1) {
			MPI_Recv(Buf, BufCountU, MPI_FLOAT, MPI_ANY_SOURCE, 99,
						MPI_COMM_WORLD, &status);
			for (i = 0, j = 0; j <= HIGG; j++) {
				for (k = 0; k <= DEPP; k++, i += 5) {
					U1_[LENN][j][k] = Buf[i  ];
					U2_[LENN][j][k] = Buf[i+1];
					U3_[LENN][j][k] = Buf[i+2];
					U4_[LENN][j][k] = Buf[i+3];
					U5_[LENN][j][k] = Buf[i+4];
				}
			}
		} // end if(myid != numprocs-1)
		// sending to previous - everyone
		for (i = 0, j = 0; j <= HIGG; j++) {
			for (k = 0; k <= DEPP; k++, i += 5) {
				Buf[i  ] = U1_[1][j][k];
				Buf[i+1] = U2_[1][j][k];
				Buf[i+2] = U3_[1][j][k];
				Buf[i+3] = U4_[1][j][k];
				Buf[i+4] = U5_[1][j][k];
			}
		}
		MPI_Send(Buf, BufCountU, MPI_FLOAT, previous, 99, MPI_COMM_WORLD);
	}

	// syncronization
	MPI_Barrier(MPI_COMM_WORLD);


	//-- x
	// left side inflow
	if(0 == myid) { /* root process */
		P = 100000.;
		R = 1.0;
		for (j = 1; j < HIGG; j++) {
			for (k = 1; k < DEPP; k++) {
				//-- choosingcomponents of velocity 
                // lower half				
				if (j <= HIG/2){
				    // always
				    U = 76.4;
				    // lower than disturbing block
				    if (j <= (HIG/2 - BL_HIG)) {
					V = 0.;
					W = 0.;				    
				    }
				    // within disturbing block
				    else {
					if ((((j - j_base)/j_period)%2) > 0) {
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
				    // always
				    U = 200.;
				    V = 0.;
				    W = 0.;
				}
				//
				U1_[0][j][k] = R;
				U2_[0][j][k] = R * U;
				U3_[0][j][k] = R * V;
				U4_[0][j][k] = R * W;
				U5_[0][j][k] = P / K_1 + 0.5 * R * (U * U + V * V + W * W);
			}
		}
	} // end if(0 == myid)

	// right side - OUTFLOW
	if(myid == numprocs-1) { /* last process */
		P = 100000.;
		for( j = 1; j < HIGG; j++ ) {
			for( k = 1; k < DEPP; k++ ) {
				U1_[LENN][j][k] = R = U1_[LEN][j][k];
				U2_[LENN][j][k] = U = U2_[LEN][j][k];
				U3_[LENN][j][k] = V = U3_[LEN][j][k];
				U4_[LENN][j][k] = W = U4_[LEN][j][k];
				U5_[LENN][j][k] = P / K_1 + 0.5 * ( U * U + V * V + W * W ) / R;
			}
		}
	} // end if(myid == numprocs-1)

	//-- z - PERIODIC
	for( i = 0; i <= LENN; i++ ) {
		for( j = 1; j < HIGG; j++ ) {
			// back side
			U1_[i][j][0] = U1_[i][j][DEP];
			U2_[i][j][0] = U2_[i][j][DEP];
			U3_[i][j][0] = U3_[i][j][DEP];
			U4_[i][j][0] = U4_[i][j][DEP];
			U5_[i][j][0] = U5_[i][j][DEP];
			// front side
			U1_[i][j][DEPP] = U1_[i][j][1];
			U2_[i][j][DEPP] = U2_[i][j][1];
			U3_[i][j][DEPP] = U3_[i][j][1];
			U4_[i][j][DEPP] = U4_[i][j][1];
			U5_[i][j][DEPP] = U5_[i][j][1];
		}
	}

	//-- y - SLIP
	for( i = 0; i <= LENN; i++ ) {
		for( k = 0; k <= DEPP; k++ ) {
			// bottom side
			U1_[i][0][k] =   U1_[i][1][k];
			U2_[i][0][k] =   U2_[i][1][k];
			U3_[i][0][k] = - U3_[i][1][k];
			U4_[i][0][k] =   U4_[i][1][k];
			U5_[i][0][k] =   U5_[i][1][k];
			// top side
			U1_[i][HIGG][k] =   U1_[i][HIG][k];
			U2_[i][HIGG][k] =   U2_[i][HIG][k];
			U3_[i][HIGG][k] = - U3_[i][HIG][k];
			U4_[i][HIGG][k] =   U4_[i][HIG][k];
			U5_[i][HIGG][k] =   U5_[i][HIG][k];
		}
	}

} // end BounCondInGhostCells()

