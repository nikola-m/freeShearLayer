/***********
*  PROBES  *  Writes readings of arrays of probes to the file "probes.dat"
***********/

#include <stdio.h>     /* printf() etc.*/
#include <stdlib.h>    /* atoi()       */

#include "type.h"
#include "global.h"   /* global variables */

#include "probes.h"

void Probes( void )
{
FILE *pFprobes;

        /*--- opens probes file ---*/
	if( ( pFprobes = fopen( "probes.dat","ab" ) ) == NULL )
	   puts( "Cannot open probes file ");
		/*--- number of current step ---*/
	fwrite( &step, sizeof(step), 1, pFprobes );
		/*--- index for the middle z-plane ---*/
	k = DEP/2+1;
		/*--- 1st probe ---*/
	i = LEN/2;
	for( j = 1; j < HIGG; j++ ) probes[j-1] = U1[i][j][k];
	fwrite( probes, sizeof(real)*HIG, 1, pFprobes );
	for( j = 1; j < HIGG; j++ ) probes[j-1] = U2[i][j][k];
	fwrite( probes, sizeof(real)*HIG, 1, pFprobes );
		/*--- 2nd probe--- */
	i = 3*LEN/4;
	for( j = 1; j < HIGG; j++ ) probes[j-1] = U1[i][j][k];
	fwrite( probes, sizeof(real)*HIG, 1, pFprobes );
	for( j = 1; j < HIGG; j++ ) probes[j-1] = U2[i][j][k];
	fwrite( probes, sizeof(real)*HIG, 1, pFprobes );
		/*--- 3rd probe---*/
	i = LEN;
	for( j = 1; j < HIGG; j++ ) probes[j-1] = U1[i][j][k];
	fwrite( probes, sizeof(real)*HIG, 1, pFprobes );
	for( j = 1; j < HIGG; j++ ) probes[j-1] = U2[i][j][k];
	fwrite( probes, sizeof(real)*HIG, 1, pFprobes );
		/* closes probes file */
	fclose( pFprobes );
} /* end Probes() */