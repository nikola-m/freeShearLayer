#ifndef DEF_H
#define DEF_H

#define M_PI 3.14159265358979323846

#define Fx1 xU1 /* aliasing x-fluxes */
#define Fx2 xU2
#define Fx3 xU3
#define Fx4 xU4
#define Fx5 xU5

#define Fy1 yU1 /* aliasing y-fluxes */
#define Fy2 yU2
#define Fy3 yU3
#define Fy4 yU4
#define Fy5 yU5

#define Fz1 zU1 /* aliasing z-fluxes */
#define Fz2 zU2
#define Fz3 zU3
#define Fz4 zU4
#define Fz5 zU5

/*--- Combinations of the ratio of specific heats ---*/
#define K     1.4     // cappa = cp/cv
#define K_1   0.4     // cappa-1
#define _K_1 (-0.4)   // -(cappa-1)
#define _05K  0.7     // cappa/2
#define K_K_1 3.5     // cappa/(cappa-1)
#define K_1_2 0.2     // (cappa-1)/2 
#define _1_2_K_1 1.25 //  0.5/(cappa-1) = 1/(2*(cappa-1))

#define Cv     717.75  /* air heat capacity at constant volume, J/(kg*K)*/
#define R_VOZD 287.10  /* gas constant         */
#define KR_VOZD 401.94 /* gas constant * 1.4   */

/*--- For monotone piecevise parabolic reconstruction ---*/
/* constant coefficients */
#define BB 4.000000000000000 // 4
#define C1 0.333333333333333 // 1/3
#define C2 0.166666666666667 // 1/6

#define small 1e-30
#define twoThirds 0.666666666666667
#define half 0.5

#define DynamicSmagorinskySGS 0 // hardcoded option

typedef int bool;
#define TRUE  1
#define FALSE 0

#endif