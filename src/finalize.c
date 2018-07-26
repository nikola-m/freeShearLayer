/*************
*  FINALIZE  *  Saves solution arrays to the binary "backup.dat" file
*************/

#include <stdio.h>     /* printf() etc.*/
#include <stdlib.h>    /* atoi()       */

#include "type.h"
#include "def.h"      /* Definitions, parameters */
#include "global.h"   /* global variables */

#include "finalize.h"

void Finalize( void )
{
	FILE *pF;

	if ( ( pF = fopen( "backup.dat", "wb") ) == NULL ) {
	   puts( "Cannot open \"backup.dat\" file to write" );
	   exit( -1 );
	}
	for ( i = 1; i < LENN; i++ ) {
		for ( j = 1; j < HIGG; j++ ) {
			fwrite( &U1[i][j][1], DEP * sizeof( U1[i][j][1] ), 1, pF );
			fwrite( &U2[i][j][1], DEP * sizeof( U2[i][j][1] ), 1, pF );
			fwrite( &U3[i][j][1], DEP * sizeof( U3[i][j][1] ), 1, pF );
			fwrite( &U4[i][j][1], DEP * sizeof( U4[i][j][1] ), 1, pF );
			fwrite( &U5[i][j][1], DEP * sizeof( U5[i][j][1] ), 1, pF );
		}
	}
	fwrite( &step, sizeof( step ), 1, pF );
	fclose( pF );

} /* end Finalize() */