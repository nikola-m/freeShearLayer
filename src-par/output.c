#include "mpi.h"

#include <stdio.h>     /* printf() etc.*/
#include <math.h>      /* sqrt()       */
#include <stdlib.h>    /* atoi()       */
#include <string.h>    /* strlen()     */ 

#include "type.h"
#include "def.h"      /* Definitions, parameters */
#include "global.h"   /* global variables */
#include "turbulence.h" /* Q Criteria */ 

#include "output.h"

/***********
*  OUTPUT  *   Outputs flowfield to separate files, each for a specific process.
***********/
void Output(int myid)
{
// float buf;
char filename[30];
MPI_File fh;
MPI_Status status;

unsigned int i, j, k;
char str[120];
real xc, yc, zc;
real R, U, V, W, P, T;
real omegax, omegay, omegaz, S12, S13, S23, Omega, Strain, Q, muT;

/* matrix of velocity derivatives */
real 
du_dx,dv_dy,dw_dz,
dv_dx, dw_dx, 
du_dy, dw_dy,
du_dz, dv_dz;


	// //--- create file
	// sprintf(filename, "%d.%d", step, myid);
	// MPI_File_open(MPI_COMM_SELF, filename, MPI_MODE_CREATE | MPI_MODE_RDWR,
	// 				  MPI_INFO_NULL, &fh );
	// MPI_File_set_view(fh, 0, MPI_FLOAT, MPI_FLOAT, "native", MPI_INFO_NULL);

	// //--- write header
	// buf = (float)LEN;    MPI_File_write(fh, &buf, 1, MPI_FLOAT, &status);
	// buf = (float)HIG;    MPI_File_write(fh, &buf, 1, MPI_FLOAT, &status);
	// buf = (float)deltaX; MPI_File_write(fh, &buf, 1, MPI_FLOAT, &status);
	// buf = (float)deltaY; MPI_File_write(fh, &buf, 1, MPI_FLOAT, &status);

	// //--- central z-plane
	// k = k_c;

	// //--- write data
	// for (j = 1; j < HIGG; j++) {
	// 	for (i = 1; i < LENN; i++) {
	// 		// vorticity plot
	// 		dv_dx = ( U3[i+1][j][k]/U1[i+1][j][k] - U3[i-1][j][k]/U1[i-1][j][k] ) * _2deltaX;
	// 		dw_dx = ( U4[i+1][j][k]/U1[i+1][j][k] - U4[i-1][j][k]/U1[i-1][j][k] ) * _2deltaX;
	// 		du_dy = ( U2[i][j+1][k]/U1[i][j+1][k] - U2[i][j-1][k]/U1[i][j-1][k] ) * _2deltaY;
	// 		dw_dy = ( U4[i][j+1][k]/U1[i][j+1][k] - U4[i][j-1][k]/U1[i][j-1][k] ) * _2deltaY;
	// 		du_dz = ( U2[i][j][k+1]/U1[i][j][k+1] - U2[i][j][k-1]/U1[i][j][k-1] ) * _2deltaZ;
	// 		dv_dz = ( U3[i][j][k+1]/U1[i][j][k+1] - U3[i][j][k-1]/U1[i][j][k-1] ) * _2deltaZ;
	// 		U = 0.5 * (dw_dy - dv_dz);
	// 		V = 0.5 * (du_dz - dw_dx);
	// 		W = 0.5 * (dv_dx - du_dy);
	// 		buf = sqrt(U*U + V*V + W*W);
	// 		MPI_File_write(fh, &buf, 1, MPI_FLOAT, &status);
	// 	}
	// }

	// //--- close file
	// MPI_File_close(&fh);


	/*-- Open and write results in tecplot file --*/

	sprintf( filename, "%d-proc%d.plt", step, myid);
	MPI_File_open(MPI_COMM_SELF, filename, MPI_MODE_CREATE | MPI_MODE_WRONLY,
					  MPI_INFO_NULL, &fh );
	MPI_File_set_view(fh, 0, MPI_CHAR, MPI_CHAR, "native", MPI_INFO_NULL);


	sprintf(str, "title     = \" 3-D compressible case \"\n");     MPI_File_write(fh, str, strlen(str), MPI_CHAR, &status);
	sprintf(str, "variables = \" x \"\n");                         MPI_File_write(fh, str, strlen(str), MPI_CHAR, &status);
	sprintf(str, "\"y\"\n");                                       MPI_File_write(fh, str, strlen(str), MPI_CHAR, &status);
	sprintf(str, "\"z\"\n");                                       MPI_File_write(fh, str, strlen(str), MPI_CHAR, &status);
	sprintf(str, "\"rho\"\n");                                     MPI_File_write(fh, str, strlen(str), MPI_CHAR, &status);
	sprintf(str, "\"u\"\n");                                       MPI_File_write(fh, str, strlen(str), MPI_CHAR, &status);
	sprintf(str, "\"v\"\n");                                       MPI_File_write(fh, str, strlen(str), MPI_CHAR, &status);
	sprintf(str, "\"w\"\n");                                       MPI_File_write(fh, str, strlen(str), MPI_CHAR, &status);
	sprintf(str, "\"p\"\n");                                       MPI_File_write(fh, str, strlen(str), MPI_CHAR, &status);
	sprintf(str, "\"T\"\n");                                       MPI_File_write(fh, str, strlen(str), MPI_CHAR, &status);
	sprintf(str, "\"Vort. mag.\"\n");                              MPI_File_write(fh, str, strlen(str), MPI_CHAR, &status);
	sprintf(str, "\"Q-criteria.\"\n");                             MPI_File_write(fh, str, strlen(str), MPI_CHAR, &status);
	sprintf(str, "\"muSgs\"\n");                                   MPI_File_write(fh, str, strlen(str), MPI_CHAR, &status);
	sprintf(str, "zone t=\" \"\n");                                MPI_File_write(fh, str, strlen(str), MPI_CHAR, &status);
	sprintf(str, "i=%d, j=%d, k=%d, f=point\n", LEN, HIG, DEP);    MPI_File_write(fh, str, strlen(str), MPI_CHAR, &status);

	for (k = 1; k < DEPP; k++) {
		for (j = 1; j < HIGG; j++) {
			for (i = 1; i < LENN; i++) {

            // Add x-axis offset to deal with your process' domain
			xc = myid*deltaX*LEN + (i-1)*deltaX + 0.5*deltaX;
			yc = (j-1)*deltaY + 0.5*deltaY;
			zc = (k-1)*deltaZ + 0.5*deltaZ;

			R = U1[i][j][k];
			U = U2[i][j][k]/R;
			V = U3[i][j][k]/R;
			W = U4[i][j][k]/R;
			P = ( U5[i][j][k] - 0.5 * R * ( U * U + V * V + W * W ) ) * K_1;
			T = P / ( R_VOZD * R );


			// Velocity gradient
			du_dx = ( U2[i+1][j][k]/U1[i+1][j][k] - U2[i-1][j][k]/U1[i-1][j][k] ) * _2deltaX;
			du_dy = ( U2[i][j+1][k]/U1[i][j+1][k] - U2[i][j-1][k]/U1[i][j-1][k] ) * _2deltaY;
			du_dz = ( U2[i][j][k+1]/U1[i][j][k+1] - U2[i][j][k-1]/U1[i][j][k-1] ) * _2deltaZ;

			dv_dx = ( U3[i+1][j][k]/U1[i+1][j][k] - U3[i-1][j][k]/U1[i-1][j][k] ) * _2deltaX;
			dv_dy = ( U3[i][j+1][k]/U1[i][j+1][k] - U3[i][j-1][k]/U1[i][j-1][k] ) * _2deltaY; 
			dv_dz = ( U3[i][j][k+1]/U1[i][j][k+1] - U3[i][j][k-1]/U1[i][j][k-1] ) * _2deltaZ;

			dw_dx = ( U4[i+1][j][k]/U1[i+1][j][k] - U4[i-1][j][k]/U1[i-1][j][k] ) * _2deltaX;
			dw_dy = ( U4[i][j+1][k]/U1[i][j+1][k] - U4[i][j-1][k]/U1[i][j-1][k] ) * _2deltaY;
			dw_dz = ( U4[i][j][k+1]/U1[i][j][k+1] - U4[i][j][k-1]/U1[i][j][k-1] ) * _2deltaZ;


			omegax = ( dw_dy - dv_dz );
			omegay = ( du_dz - dw_dx );
			omegaz = ( dv_dx - du_dy );

			Omega = sqrt(  omegax * omegax + omegay * omegay + omegaz * omegaz );

			S12 = 0.5 * (du_dy + dv_dx);
			S13 = 0.5 * (du_dz + dw_dx);
			S23 = 0.5 * (dv_dz + dw_dy);

			Strain = sqrt( 2 * ( du_dx * du_dx + dv_dy * dv_dy + dw_dz * dw_dz 
			                      + 2 * ( S12   * S12   + S13   * S13   + S23   * S23 ) ) );

			Q = 0.5 * ( Omega*Omega - Strain*Strain);

			if (DynamicSmagorinskySGS) { 
				muT = mu_SGS[i][j][k];
			}else{
				muT = R * CsDD * Strain;		
			}


			sprintf(str, "%g %g %g %g %g %g %g %g %g %g %g %g\n", xc, yc, zc, R, U, V, W, P, T, Omega, Q, muT); 			
			MPI_File_write(fh, str, strlen(str), MPI_CHAR, &status);

			}
		}
	}
	//--- close file
	MPI_File_close(&fh);




} // end Output()