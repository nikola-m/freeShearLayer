/* layer2_p.c ******************\
*                               *
*   Reynolds/Favre averaging    *
*  of computed DNS or LES data  *
*    for mixing layer flow      *
*                               *
*  v. 0.2, 23.11.2000           *
*                               *
\*******************************/

#include <stdio.h>
#include <stdlib.h>

float *ArrayAllocMem( unsigned size )
{
float *x;

   if( ( x = (float *)malloc( sizeof(float)*size ) ) == NULL ) {
	  puts( "Cannot allocate memory" );
	  exit( -1 );
   }
   return x;

} /* end ArrayAllocMem() */


int main( void )
{
unsigned HIG, step, stepStart, stepEnd, i, j;
float *probeR, *probeRU, *R1, *RU1, *R2, *RU2, *R3, *RU3;
FILE *pFprobes, *pFprofile;

   /*--- get some parameters ---*/
   printf( "stepStart = " ); scanf( "%d", &stepStart );
   printf( "stepEnd = " );   scanf( "%d", &stepEnd );
   printf( "HIG = " );       scanf( "%d", &HIG );

   /*--- allocate memory ---*/
   probeR  = ArrayAllocMem( HIG );
   probeRU = ArrayAllocMem( HIG );
   R1 = ArrayAllocMem( HIG );
   R2 = ArrayAllocMem( HIG );
   R3 = ArrayAllocMem( HIG );
   RU1 = ArrayAllocMem( HIG );
   RU2 = ArrayAllocMem( HIG );
   RU3 = ArrayAllocMem( HIG );
   for( i = 0; i < HIG; i++ ) {
      R1[i] = R2[i] = R3[i] =0.;
      RU1[i] = RU2[i] = RU3[i] = 0.;
   }

   /*--- opens probes file ---*/
   if( ( pFprobes = fopen( "probes.dat", "rb" ) ) == NULL ) {
      puts( "Cannot open file \"probes.dat\"" );
      exit( -1 );
   }

   /*--- read ... ---*/
   for( i = 0;; i++ ) {
   
      /* read step */
	  fread( &step, sizeof( step ), 1, pFprobes );

	  if( step > stepEnd ) break;

	  if( step < stepStart ) { /* simple read, not add */
		 /* 1st probe - idle */
		 fread( probeR,  sizeof(float)*HIG, 1, pFprobes );
		 fread( probeRU, sizeof(float)*HIG, 1, pFprobes );
		 /* 2nd probe - idle */
		 fread( probeR,  sizeof(float)*HIG, 1, pFprobes );
		 fread( probeRU, sizeof(float)*HIG, 1, pFprobes );
		 /* 3rd probe - idle */
		 fread( probeR,  sizeof(float)*HIG, 1, pFprobes );
		 fread( probeRU, sizeof(float)*HIG, 1, pFprobes );
	  } /* end if() */
	  else{
		 /* prints it out */
		 printf( "step = %d\n", step );
		 /* 1st probe */
		 fread( probeR,  sizeof(float)*HIG, 1, pFprobes );
		 fread( probeRU, sizeof(float)*HIG, 1, pFprobes );
		 for( j = 0; j < HIG; j++ ) {
			R1[j] += probeR[j];
			RU1[j] += probeRU[j];
		 }
		 /* 2nd probe */
		 fread( probeR,  sizeof(float)*HIG, 1, pFprobes );
		 fread( probeRU, sizeof(float)*HIG, 1, pFprobes );
		 for( j = 0; j < HIG; j++ ) {
			R2[j] += probeR[j];
			RU2[j] += probeRU[j];
		 }
		 /* 3rd probe */
		 fread( probeR,  sizeof(float)*HIG, 1, pFprobes );
		 fread( probeRU, sizeof(float)*HIG, 1, pFprobes );
		 for( j = 0; j < HIG; j++ ) {
			R3[j] += probeR[j];
			RU3[j] += probeRU[j];
		 }
	  } /* end else */
   } /* end for() */

   /*--- close probes file ---*/
   fclose( pFprobes );

   /*--- determine Favre averaged (here RU means F.av. U) ---*/
   for( j = 0; j < HIG; j++ ) {
	  RU1[j] /= R1[j];
	  RU2[j] /= R2[j];
	  RU3[j] /= R3[j];
   }

   /*--- write profiles ---*/
	 /* 1st */
   if( ( pFprofile = fopen( "probe1.dat","w" ) ) == NULL ) {
	  puts( "Cannot open file ");
	  exit( -1 );
   }
   fprintf( pFprofile, "%d\n1\n", HIG+2 );
   fprintf( pFprofile, "%f %3.2f \n", -0.00001, 0. );
   for( j = 0; j < HIG; j++ )
	  fprintf( pFprofile, "%d %3.2f \n", j, RU1[j] );
   fprintf( pFprofile, "%f %3.2f \n", HIG-1+0.00001, 210. );
   fclose( pFprofile );
	 /* 2nd */
   if( ( pFprofile = fopen( "probe2.dat","w" ) ) == NULL ) {
	  puts( "Cannot open file ");
	  exit( -1 );
   }
   fprintf( pFprofile, "%d\n1\n", HIG+2 );
   fprintf( pFprofile, "%f %3.2f \n", -0.00001, 0. );
   for( j = 0; j < HIG; j++ )
	  fprintf( pFprofile, "%d %3.2f \n", j, RU2[j] );
   fprintf( pFprofile, "%f %3.2f \n", HIG-1+0.00001, 210. );
   fclose( pFprofile );
	 /* 3rd */
   if( ( pFprofile = fopen( "probe3.dat", "w" ) ) == NULL ) {
	  puts( "Cannot open file ");
	  exit( -1 );
   }
   fprintf( pFprofile, "%d\n1\n", HIG+2 );
   fprintf( pFprofile, "%f %3.2f \n", -0.00001, 0. );
   for( j = 0; j < HIG; j++ )
	  fprintf( pFprofile, "%d %3.2f \n", j, RU3[j] );
   fprintf( pFprofile, "%f %3.2f \n", HIG-1+0.00001, 210. );
   fclose( pFprofile );

   /*---*/
   return 0;

} /* end main()*/

/* end layer2_p.c */