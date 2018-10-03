/***********
*  PROBES  *  Writes readings of arrays of probes to the file "probes.dat"
***********/

#include <stdio.h>     /* printf() etc.*/
#include <stdlib.h>    /* atoi()       */

#include "type.h"
#include "global.h"   /* global variables */

#include "probes.h"

int Probes(int myid)
{
	//--- probes are positioned on line formed by intersection of
	//    cenral x-, z-planes
	i = LEN/2+1;
	k = DEP/2+1;

	//--- output to "probes.myid"
	// step
	fwrite(&step, sizeof(step), 1, pFprobes);
	// HIG-sized blocks of current R and RU
		//
	for (j = 1; j < HIGG; j++)
		probes[j-1] = U1[i][j][k];
	fwrite(probes, sizeof(probes[0]) * HIG, 1, pFprobes);
		//
	for (j = 1; j < HIGG; j++)
		probes[j-1] = U2[i][j][k];
	fwrite(probes, sizeof(probes[0]) * HIG, 1, pFprobes);

	//---
	return 1;
} // end Probes()