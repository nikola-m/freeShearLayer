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




  /*--- Communication via MPI to complete ghost cell x-layers ---*/

  /* ranks of the next and the previous processes */
  next     = myid + 1;
  previous = myid - 1;

  /* preparing array for sending to the previous - every process except 0th */
  if (myid != 0) {
    for (i = 0, j = 0; j <= HIGG; j++) {
      for (k = 0; k <= DEPP; k++) {
		Bufp[i++] = U1_[1][j][k];
		Bufp[i++] = U2_[1][j][k];
		Bufp[i++] = U3_[1][j][k];
		Bufp[i++] = U4_[1][j][k];
		Bufp[i++] = U5_[1][j][k];
      }
    }
  }
  /* preparing array for sending to the next - every process except the last, (numprocs-1)th */
  if (myid != numprocs-1) {
    for (i = 0, j = 0; j <= HIGG; j++) {
      for (k = 0; k <= DEPP; k++) {
        Bufn[i++] = U1_[LEN][j][k];
        Bufn[i++] = U2_[LEN][j][k];
        Bufn[i++] = U3_[LEN][j][k];
        Bufn[i++] = U4_[LEN][j][k];
        Bufn[i++] = U5_[LEN][j][k];
      }
    }
  }

  /*-- communicating via MPI --*/

  /* syncronization before communication */
  MPI_Barrier(MPI_COMM_WORLD);
  
  /* sending to the previous - every process except 0th */
  if (myid != 0)
    MPI_Send(Bufp, BufCountU, MPI_FLOAT, previous, 99, MPI_COMM_WORLD);
    
  /* receiving from the next - every process except the last, (numprocs-1)th */
  if (myid != numprocs-1)
    MPI_Recv(nBuf, BufCountU, MPI_FLOAT, next, 99, MPI_COMM_WORLD, &status);
    
  /* sending to the next - every process except the last, (numprocs-1)th */
  if (myid != numprocs-1)
    MPI_Send(Bufn, BufCountU, MPI_FLOAT, next, 101, MPI_COMM_WORLD);
    
  /* receiving from the previous, every  process except 0th */ 
  if (myid != 0)
    MPI_Recv(pBuf, BufCountU, MPI_FLOAT, previous, 101, MPI_COMM_WORLD, &status);
    
  /*-- end of communicating via MPI --*/
  
  /*- processing received arrays -*/
  /* every process except 0th */ 
  if (myid != 0) {  
    for (i = 0, j = 0; j <= HIGG; j++) {
      for (k = 0; k <= DEPP; k++) {
		U1_[0][j][k] = pBuf[i++];
		U2_[0][j][k] = pBuf[i++];
		U3_[0][j][k] = pBuf[i++];
		U4_[0][j][k] = pBuf[i++];
		U5_[0][j][k] = pBuf[i++];
      }
    }
  }
  /* every process except the last, (numprocs-1)th */
  if (myid != numprocs-1) {
    for (i = 0, j = 0; j <= HIGG; j++) {
      for (k = 0; k <= DEPP; k++) {
		U1_[LENN][j][k] = nBuf[i++];
		U2_[LENN][j][k] = nBuf[i++];
		U3_[LENN][j][k] = nBuf[i++];
		U4_[LENN][j][k] = nBuf[i++];
		U5_[LENN][j][k] = nBuf[i++];
      }
    }
  }
  /*--- End of communication to complete ghost cell x-layers ---*/  


	//-- x
	// left side inflow
	if(0 == myid) { /* root process */

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
			}
		}

		U = 200; /* V = 0.; W = 0.; */
		for( j = HIG/2+1; j < HIGG; j++ ) {
			for( k = 1; k < DEPP; k++ ) {
				U1_[0][j][k] = R;
				U2_[0][j][k] = R * U;
				U3_[0][j][k] = 0.;/*R * V; */
				U4_[0][j][k] = 0.;/*R * W; */
				U5_[0][j][k] = P / K_1 + 0.5 * R * ( U * U /*+ V * V + W * W*/ );
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

