#include "type.h"      /*the "real" type */
#include "def.h"       /* Definitions, parameters */
#include "global.h"    /* global variables */
#include "mpi.h"

void exchange( real ***Phi, int myid, int numprocs ){ 
/* 

    --- Communicate via MPI to complete cell x-layers ---
                                 __________
                                v          |
    |-gc-|----|-- ... --|-----|-gc-|    |-gc-|----|-- ... --|-----|-gc-|  ... |-gc-|----|-- ... --|-----|-gc-| 
                 proc 0        |_________^          proc 1                       numproc

    Buf - an array of type real (float or double) that stores the data to be sent or received. Look at 
    global.h for its declaration and initialize.c for allocation length which is 5 * (HIG+2) * (DEP+2).
    Here we use only a portion of it, since we send and receive only one field at interfaces. 
    This portion is of length BufCount = HIGG*DEPP.

*/

    int i,j,k;
    int next, previous; // processes next/previous to this
    MPI_Status status;


  /*--- Communication via MPI to complete ghost cell x-layers ---*/

  // The length of the Araray that is exchanged between processes. We use Buf array, that is longer (around five times).
  int BufCount = HIGG*DEPP;

  /* ranks of the next and the previous processes */
  next     = myid + 1;
  previous = myid - 1;

  /* preparing array for sending to the previous - every process except 0th */
  if (myid != 0) {
    for (i = 0, j = 0; j <= HIGG; j++) {
      for (k = 0; k <= DEPP; k++) {
        Bufp[i++] = Phi[1][j][k];
      }
    }
  }
  /* preparing array for sending to the next - every process except the last, (numprocs-1)th */
  if (myid != numprocs-1) {
    for (i = 0, j = 0; j <= HIGG; j++) {
      for (k = 0; k <= DEPP; k++) {
        Bufn[i++] = Phi[LEN][j][k];
      }
    }
  }

  /*-- communicating via MPI --*/

  /* syncronization before communication */
  MPI_Barrier(MPI_COMM_WORLD);
  
  /* sending to the previous - every process except 0th */
  if (myid != 0)
    MPI_Send(Bufp, BufCount, MPI_FLOAT, previous, 99, MPI_COMM_WORLD);
    
  /* receiving from the next - every process except the last, (numprocs-1)th */
  if (myid != numprocs-1)
    MPI_Recv(nBuf, BufCount, MPI_FLOAT, next, 99, MPI_COMM_WORLD, &status);
    
  /* sending to the next - every process except the last, (numprocs-1)th */
  if (myid != numprocs-1)
    MPI_Send(Bufn, BufCount, MPI_FLOAT, next, 101, MPI_COMM_WORLD);
    
  /* receiving from the previous, every  process except 0th */ 
  if (myid != 0)
    MPI_Recv(pBuf, BufCount, MPI_FLOAT, previous, 101, MPI_COMM_WORLD, &status);
    
  /*-- end of communicating via MPI --*/
  
  /*- processing received arrays -*/
  /* every process except 0th */ 
  if (myid != 0) {  
    for (i = 0, j = 0; j <= HIGG; j++) {
      for (k = 0; k <= DEPP; k++) {
        Phi[0][j][k] = pBuf[i++];
      }
    }
  }
  /* every process except the last, (numprocs-1)th */
  if (myid != numprocs-1) {
    for (i = 0, j = 0; j <= HIGG; j++) {
      for (k = 0; k <= DEPP; k++) {
        Phi[LENN][j][k] = nBuf[i++];
      }
    }
  }
  /*--- End of communication to complete ghost cell x-layers ---*/  

//     //-- sending to next processes

//     // next process's ID
//     next = myid+1;

//     // The length of the Araray that is exchanged between processes. We use Buf array, that is longer (around five times).
//     int BufCount = HIGG*DEPP;

//     // root process
//     if (myid == 0) {
//         // sending to the next
//         for (i = 0, j = 0; j <= HIGG; j++) {
//             for (k = 0; k <= DEPP; k++, i += 1) {
//                 Buf[i  ] = Phi[LEN][j][k];
//             }
//         }
//         MPI_Send(Buf, BufCount, MPI_FLOAT, next, 99, MPI_COMM_WORLD);
//         // previous doesn't exist
// /*      MPI_Recv(Buf, BufCount, MPI_FLOAT, MPI_ANY_SOURCE, 99, MPI_COMM_WORLD,
//                     &status);
//         for (i = 0, j = 0; j <= HIGG; j++) {
//             for (k = 0; k <= DEPP; k++, i += 1) {
//                 Phi[0][j][k] = Buf[i  ];
//             }
//         }
// */
//     } // end if()

//     // other processes
//     else {
//         // receiving from the previous, everyone
//         MPI_Recv(Buf, BufCount, MPI_FLOAT, MPI_ANY_SOURCE, 99, MPI_COMM_WORLD,
//                     &status);
//         for (i = 0, j = 0; j <= HIGG; j++) {
//             for (k = 0; k <= DEPP; k++, i += 1) {
//                 Phi[0][j][k] = Buf[i  ];
//             }
//         }
//         // sending to the next; if not the last!
//         if(myid != numprocs-1) {
//             for (i = 0, j = 0; j <= HIGG; j++) {
//                 for (k = 0; k <= DEPP; k++, i += 1) {
//                     Buf[i  ] = Phi[LEN][j][k];
//                 }
//             }
//             MPI_Send(Buf, BufCount, MPI_FLOAT, next, 99, MPI_COMM_WORLD);
//         } // end if(myid != numprocs-1)
//     }

//     //-- sending to previous processes

//     // previous process's ID
//     previous = myid-1;

//     // root process
//     if(0 == myid) {
//         // previous doesn't exist
// /*      for (i = 0, j = 0; j <= HIGG; j++) {
//             for (k = 0; k <= DEPP; k++, i += 1) {
//                 Buf[i  ] = Phi[1][j][k];
//             }
//         }
//         MPI_Send(Buf, BufCountU, MPI_FLOAT, previous, 99, MPI_COMM_WORLD);
// */
//         // receiving from the next
//         MPI_Recv(Buf, BufCount, MPI_FLOAT, MPI_ANY_SOURCE, 99, MPI_COMM_WORLD,
//                     &status);
//         for (i = 0, j = 0; j <= HIGG; j++) {
//             for (k = 0; k <= DEPP; k++, i += 1) {
//                 Phi[LENN][j][k] = Buf[i  ];
//             }
//         }
//     } // end if()

//     // other processes
//     else {
//         // receiving from the next; if not the last process!
//         if(myid != numprocs-1) {
//             MPI_Recv(Buf, BufCount, MPI_FLOAT, MPI_ANY_SOURCE, 99,
//                         MPI_COMM_WORLD, &status);
//             for (i = 0, j = 0; j <= HIGG; j++) {
//                 for (k = 0; k <= DEPP; k++, i += 1) {
//                     Phi[LENN][j][k] = Buf[i  ];
//                 }
//             }
//         } // end if(myid != numprocs-1)
//         // sending to previous - everyone
//         for (i = 0, j = 0; j <= HIGG; j++) {
//             for (k = 0; k <= DEPP; k++, i += 1) {
//                 Buf[i  ] = Phi[1][j][k];
//             }
//         }
//         MPI_Send(Buf, BufCount, MPI_FLOAT, previous, 99, MPI_COMM_WORLD);
//     }

//     // syncronization
//     MPI_Barrier(MPI_COMM_WORLD);

} /* End function exchange */