/* ini2bin.c ********************************\
*                                            *
*   Converter of *.ini-files to *.bin-file   *
* (designed for use with "mpi_layer2.ini")   *
*                                            *
\********************************************/

#include <stdio.h>
#include <stdlib.h>

main()
{
FILE *pFin, *pFout;
char str[80];
unsigned i;
float buf;

    //--- open "mpi_layer2.ini"
    if((pFin = (FILE *)fopen("mpi_layer2.ini", "r")) == NULL) {
	fprintf(stderr, "can't open \"mpi_layer2.ini\".\n");
	exit(-1);
    }

    //--- open "mpi_layer2.bin" in binary mode
    if((pFout = (FILE *)fopen("mpi_layer2.bin", "wb")) == NULL) {
	fprintf(stderr, "can't open \"mpi_layer2.bin\".\n");
	exit(-1);
    }

    //--- read, report and write
    i = 0;
    do{
	if (fgets(str, 79, pFin) != NULL) {
	    switch (i) {
			case  0: fprintf(stdout, "LEN    = %d\n", atoi(str)); break;
			case  1: fprintf(stdout, "HIG    = %d\n", atoi(str)); break;
			case  2: fprintf(stdout, "DEP    = %d\n", atoi(str)); break;
			case  3: fprintf(stdout, "deltaX = %f\n", atof(str)); break;
			case  4: fprintf(stdout, "deltaY = %f\n", atof(str)); break;
			case  5: fprintf(stdout, "deltaZ = %f\n", atof(str)); break;
			case  6: fprintf(stdout, "deltaT = %f\n", atof(str)); break;
			case  7: fprintf(stdout, "Answer = %d\n", atoi(str)); break;
			case  8: fprintf(stdout, "numstep= %d\n", atoi(str)); break;
			case  9: fprintf(stdout, "mu_L   = %f\n", atof(str)); break;
			case 10: fprintf(stdout, "Pr_L   = %f\n", atof(str)); break;
			case 11: fprintf(stdout, "Pr_T   = %f\n", atof(str)); break;
			case 12: fprintf(stdout, "Cs     = %f\n", atof(str)); break;
			case 13: fprintf(stdout, "bl_hig = %d\n", atoi(str)); break;
			case 14: fprintf(stdout, "Ua     = %f\n", atof(str)); break;						
			case 15: fprintf(stdout, "Va     = %f\n", atof(str)); break;
			case 16: fprintf(stdout, "Wa     = %f\n", atof(str)); break;
			case 17: fprintf(stdout, "Ns_min = %d\n", atoi(str)); break;
			case 18: fprintf(stdout, "Nst    = %d\n", atoi(str)); break;
			case 19: fprintf(stdout, "f_step = %d\n", atoi(str)); break;
			case 20: fprintf(stdout, "nStages = %d\n", atoi(str)); break;
			case 21: fprintf(stdout, "maxCoNum = %f\n", atof(str)); break;	

			default:
		    fprintf(stderr, "error reading \"mpi_layer2.ini\".\n");
		    exit(-1);
	    } // end switch
	    //
	    buf = atof(str);
	    fwrite(&buf, sizeof(float), 1, pFout);
	    //
	    i++;
	} // end if()
	else {
	    fprintf(stderr, "error reading \"mpi_layer2.ini\".\n");
	    exit(-1);
	}
    } while (i < 22);

    //---
    fclose(pFout);
    fclose(pFin);

    //---
    return 0;

} // end main()

//--- end ini2bin.c ---
