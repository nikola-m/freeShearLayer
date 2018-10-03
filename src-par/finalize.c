/*************
*  FINALIZE  *  Saves solution arrays to the binary "backup.dat" file
*************/
#include "mpi.h"
#include <stdio.h>     /* printf() etc.*/
#include <stdlib.h>    /* atoi()       */

#include "type.h"
#include "def.h"      /* Definitions, parameters */
#include "global.h"   /* global variables */

#include "finalize.h"

void Finalize(int myid)
{
float buf;
char filename[30];
MPI_File fh;
MPI_Status status;


	//--- open/create "backup.myid" file in binary mode
	sprintf(filename, "backup.%d", myid);
	MPI_File_open(MPI_COMM_SELF, filename, MPI_MODE_CREATE | MPI_MODE_WRONLY,
						MPI_INFO_NULL, &fh);
	MPI_File_set_view(fh, 0, MPI_FLOAT, MPI_FLOAT, "native", MPI_INFO_NULL);

	//--- write data
	for (i = 1; i < LENN; i++) {
		for (j = 1; j < HIGG; j++) {
			MPI_File_write(fh, &U1[i][j][1], DEP, MPI_FLOAT, &status);
			MPI_File_write(fh, &U2[i][j][1], DEP, MPI_FLOAT, &status);
			MPI_File_write(fh, &U3[i][j][1], DEP, MPI_FLOAT, &status);
			MPI_File_write(fh, &U4[i][j][1], DEP, MPI_FLOAT, &status);
			MPI_File_write(fh, &U5[i][j][1], DEP, MPI_FLOAT, &status);
		}
	}
	// save current step
	buf = (float)step;
	MPI_File_write(fh, &buf, 1, MPI_FLOAT, &status);
	// save seed of pseudo-random number
	buf = (float)rand();
	MPI_File_write(fh, &buf, 1, MPI_FLOAT, &status);

	//--- close "backup.myid" file
	MPI_File_close(&fh);

	//--- close file "probes.dat" here
	fclose(pFprobes);

} // end Finalize()