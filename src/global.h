#ifndef GLOBAL_H
#define GLOBAL_H

/*=== Global variables ===*/

extern real
	***U1_, ***U2_, ***U3_, ***U4_, ***U5_, /* substitutional pointers */
	  /* Cons. variables at the cells' centroids */
	***U1,  ***U2,  ***U3,  ***U4,  ***U5,
	***U1p, ***U2p, ***U3p, ***U4p, ***U5p,
	  /* Cons. variables at cell boundaries in x-direction */
	***xU1, ***xU2, ***xU3, ***xU4, ***xU5,
	***U1x, ***U2x, ***U3x, ***U4x, ***U5x,
	  /* Cons. variables at cell boundaries in y-direction */
	***yU1, ***yU2, ***yU3, ***yU4, ***yU5,
	***U1y, ***U2y, ***U3y, ***U4y, ***U5y,
	  /* Cons. variables at cell boundaries in y-direction */
	***zU1, ***zU2, ***zU3, ***zU4, ***zU5,
	***U1z, ***U2z, ***U3z, ***U4z, ***U5z,

	R, U, V, W, P, C, /* vector of primitive flow parameters */
	u1, u2, u3, u4, u5, /* conservative flow variables */

	deltaX, deltaY, deltaZ, /* cell spacings */
	deltaT, /* time step */
	deltaT_X, deltaT_Y, deltaT_Z, /* convenient ratios */
	_deltaX, _deltaY, _deltaZ,
	_2deltaX, _2deltaY, _2deltaZ,
	_4deltaX, _4deltaY, _4deltaZ;

extern real
	/* transport coefficients */
	/* molecular */
	mu_L, lambda_L, Pr_L,
	/* subgrid-scale */
	mu_T, Pr_T, cp_Pr_T,
	/* effective */
	mu_E, lambda_E,
	/* Smagorinsky constant, complex CsDD = Cs * delta * delta, where delta = (dX*dY*dZ)~0.33333 */
	Cs, CsDD;

extern real maxCoNum;  /* maximum Courant number */
extern int nStages; /* number of stages of Runge-Kutta algorithm */

extern unsigned
   step,    /* current time step */
   numstep, /* overall prescribed number of time steps */
   mem,     /* total amount of dynamic memory required (in bytes) */
   f_step,  /* frame taking step */
   
   Answer, /* solution continuation flag      */
   LEN, HIG, DEP,  /* cell numbers in x-, y-, z-directions, respectively for 3d box */
   LENN, HIGG, DEPP,
   i,  j,  k, i__,  j__,  k__, /* indices */
   l;
   
extern char
   Stage,    /* indicator of the current stage */
   str[80];  /* buffer string */

extern real *probes; /*--- Array of probes ---*/

/*--- for the specification of disturbed velocity block ---*/
extern unsigned BL_HIG,  /* height of the disturbing block            */
	  k_min, k_max;      /* its "minimal"-"maximal" indices      */
extern real Ud, Vd, Wd, /* curent disturbed velocities           */
	  Ua, Va, Wa;     /* their amplitude values                */
extern int Ns_min, Nst;  /* minimal period of existence of this disturbing block and its max-min */
extern int continFlag;
extern double X;
	
#endif 
