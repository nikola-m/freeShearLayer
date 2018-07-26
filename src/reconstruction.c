/*
*  RECONSTRUCTION  
*  
*  Piecewise-parabolic reconstruction of the cell parameters with 2..3 order in space
*
*/
#include <math.h>      /* sqrt()       */
#include "type.h"
#include "def.h"      /* Definitions, parameters */
#include "global.h"   /* global variables */
#include "reconstruction.h"


/* 
* Chakravarthy, Osher minmod limiter 
*/
real minmod( real x, real y ) 
{  
	real ay;
	if( x*y<=0 ) return 0;
	ay = ( (y*=4.0) > 0 ) ? y: -y;
	if( x > 0 ) return (x<ay) ? x : ay; // if x is less then zero return the least between x and ay
	return ( (-x<ay) ? x : -ay );        //
} /* end minmod() */

void Reconstruction( void )
{
/*register*/ real kk;
/*register*/ unsigned i, j, k, _i, _j, _k, i_, j_, k_;


real
	/* finite differences of conservative variables */
      m1, m2, m3, m4, m5,
	  w1, w2, w3, w4, w5,
	/* finite differences of characteristic variables */
	  M1, M2, M3, M4, M5,
	  W1, W2, W3, W4, W5,
	/* modified finite differences of characteristic variables */
	 _M1, _M2, _M3, _M4, _M5,
	 _W1, _W2, _W3, _W4, _W5,
	/* transformation matrix [S-1] */
	 S11, S12,
	 S41, S42,
	 S51, S52,
	/* inverted transformation matrix [S-1] */
	 _S24, _S25,
     _S54, _S55,
	/* subsidiary */
	 aa, bb, cc, dd, ee, ff, gg, hh, ll, mm, nn,
	 u_hh, v_hh, w_hh, u_ee, v_ee, w_ee,
	 k1, k2, k3, k4, k5;


	for( i = 1, _i = 0, i_ = 2; i < LENN; i++, _i++, i_++ ) {
		for( j = 1, _j = 0, j_ = 2; j < HIGG; j++, _j++, j_++ ) {
			for( k = 1, _k = 0, k_ = 2; k < DEPP; k++, _k++, k_++ ) {
				
				/*---  X  ---*/

				/* at cell centroid */

				// Conservative variables at cell centroid
				u1 = U1_[i][j][k];
				u2 = U2_[i][j][k];
				u3 = U3_[i][j][k];
				u4 = U4_[i][j][k];
				u5 = U5_[i][j][k];

				// Primitive variables at cell centroid
				R =   1. / u1;   /* 1 / R  */
				U =   u2 * R;
				V =   u3 * R;
				W =   u4 * R;
				P = ( u5 - 0.5 * ( u2 * u2 + u3 * u3 + u4 * u4 ) * R ) * K_1;
				C = sqrt( K * P * R );

				/* evaluation of matrices */
				ff = U*U+V*V+W*W;   /*  */
				aa = K_1_2*ff;      /* (K-1)*(U*U+V*V+W*W)/2 */
				bb = _K_1*U;        /* - (K_1) * U           */
				cc = _K_1*V;        /* - (K_1) * V           */
				dd = _K_1*W;        /* - (K_1) * W           */
				ll = 1./(C+C);      /* 1/(2.*C)              */
				ee = ll/C;          /* 1./(2*C*C)            */
				gg = - ff * ee;     /* -(U*U+V*V+W*W)/(2*C*C)  */
				hh = - ( ee + ee ); /* - 1 / (C*C)            */
				mm = U*ll;          /* U/(2.*C)               */

				/* differencing of the conservative variables */
				/* forward */
				m1 = U1_[i_][j][k] - u1;
				m2 = U2_[i_][j][k] - u2;
				m3 = U3_[i_][j][k] - u3;
				m4 = U4_[i_][j][k] - u4;
				m5 = U5_[i_][j][k] - u5;
				/* backward */
				w1 = u1 - U1_[_i][j][k];
				w2 = u2 - U2_[_i][j][k];
				w3 = u3 - U3_[_i][j][k];
				w4 = u4 - U4_[_i][j][k];
				w5 = u5 - U5_[_i][j][k];
				/* transformation matrix [S] */
				S11 = aa - C*C;
				S41 = aa - (kk=C*U);   /* CU */
				S51 = aa + kk;
				S12 = S42 = S52 = bb;
				S42 += C;
				S52 -= C;

				/* finite differences of characteristic variables dW = S * dU */
				/* forward */
				kk = K_1 * m5;
				k2 = cc * m3;
				k3 = dd * m4;
				M1 = S11*m1 + S12*m2 +   k2 +   k3 + kk;
				M2 =  -V*m1 		 +   m3            ;
				M3 =  -W*m1 		 		+ 	m4     ;
				M4 = S41*m1 + S42*m2 +   k2 +   k3 + kk;
				M5 = S51*m1 + S52*m2 +   k2 +   k3 + kk;
				/* backward */
				kk = K_1 * w5;
				k2 = cc * w3;
				k3 = dd * w4;
				W1 = S11*w1 + S12*w2 +   k2 +   k3 + kk;
				W2 =  -V*w1 		 +   w3 	       ;
				W3 =  -W*w1 				+   w4     ;
				W4 = S41*w1 + S42*w2 +   k2 +   k3 + kk;
				W5 = S51*w1 + S52*w2 +   k2 +   k3 + kk;
				/* modified characteristic differences _W */
				_M1 = minmod( M1, W1 ); /* the best way */
				_M2 = minmod( M2, W2 );
				_M3 = minmod( M3, W3 );
				_M4 = minmod( M4, W4 );
				_M5 = minmod( M5, W5 );
				_W1 = minmod( W1, M1 );
				_W2 = minmod( W2, M2 );
				_W3 = minmod( W3, M3 );
				_W4 = minmod( W4, M4 );
				_W5 = minmod( W5, M5 );
				/* inverted transformation matrix [S-1] */
				/*_S21 = */ u_hh = U*hh; /*_S31 =*/ v_hh = V*hh; /*_S41 =*/ w_hh = W*hh;
							u_ee = U*ee; /*_S34 =*/ v_ee = V*ee; /*_S44 =*/ w_ee = W*ee;
				 _S24 = u_ee + ll; _S25 = u_ee - ll;
				_S54 = _S55 = ( nn = - 0.5 * gg ) + _1_2_K_1;
				_S54 += mm; _S55 -= mm;
				/* parameters' vector at the rightmost interface */
				k1 = _M1 * C1 + _W1 * C2;
				k2 = _M2 * C1 + _W2 * C2;
				k3 = _M3 * C1 + _W3 * C2;
				k4 = _M4 * C1 + _W4 * C2;
				k5 = _M5 * C1 + _W5 * C2;
				kk = k4 + k5;
				xU1[ i][_j][_k] = u1 +   hh * k1 + 						 	ee  * kk;
				xU2[ i][_j][_k] = u2 + u_hh * k1 + 						   _S24 * k4 + _S25 * k5;
				xU3[ i][_j][_k] = u3 + v_hh * k1 + 		  k2 			 + v_ee * kk;
				xU4[ i][_j][_k] = u4 + w_hh * k1 + 					  k3 + w_ee * kk;
				xU5[ i][_j][_k] = u5 +   gg * k1 +    V * k2 +    W * k3 + _S54 * k4 + _S55 * k5;
				/* ... leftmost interface */
				k1 = _M1 * C2 + _W1 * C1;
				k2 = _M2 * C2 + _W2 * C1;
				k3 = _M3 * C2 + _W3 * C1;
				k4 = _M4 * C2 + _W4 * C1;
				k5 = _M5 * C2 + _W5 * C1;
				kk =  k4 + k5;
				U1x[_i][_j][_k] = u1 -   hh * k1 - 						     ee * kk;
				U2x[_i][_j][_k] = u2 - u_hh * k1 - 						   _S24 * k4 - _S25 * k5;
				U3x[_i][_j][_k] = u3 - v_hh * k1 - 		  k2 			 - v_ee * kk;
				U4x[_i][_j][_k] = u4 - w_hh * k1 - 					  k3 - w_ee * kk;
				U5x[_i][_j][_k] = u5 -   gg * k1 -    V * k2 -    W * k3 - _S54 * k4 - _S55 * k5;


				/*---  Y  ---*/

				/* evaluation of matrices */
				kk = bb; bb = cc; cc = dd; dd = kk;
				mm = V*ll;          /* V/(2.*C) */
				/* differencing of the conservative variables */
					/* forvard */
				m1 = U1_[i][j_][k] - u1;
				m2 = U3_[i][j_][k] - u3; /* U3[i][j+1][k] ( U->V, V->W, W->U : 3-4-2 ) */
				m3 = U4_[i][j_][k] - u4; /* U4[i][j+1][k] */
				m4 = U2_[i][j_][k] - u2; /* U2[i][j+1][k] */
				m5 = U5_[i][j_][k] - u5;
					/* backward */
				w1 = u1 - U1_[i][_j][k];
				w2 = u3 - U3_[i][_j][k]; /* U3[i][j-1][k] ( U->V, V->W, W->U : 3-4-2 ) */
				w3 = u4 - U4_[i][_j][k]; /* U4[i][j-1][k] */
				w4 = u2 - U2_[i][_j][k]; /* U2[i][j-1][k] */
				w5 = u5 - U5_[i][_j][k];
				/* transformation matrix [S] */
				S41 = aa - (kk=C*V);    /* CV */
				S51 = aa + kk;
				S12 = S42 = S52 = bb;
				S42 += C;
				S52 -= C;
				/* finite differences of characteristic variables dW = S * dU */
					/* forward */
				kk = K_1 * m5;
				k2 = cc * m3;
				k3 = dd * m4;
				M1 = S11*m1 + S12*m2 +   k2 +   k3 + kk;
				M2 =  -W*m1 		 +   m3            ;  /* -W*m1 ! */
				M3 =  -U*m1 		 		+ 	m4     ;  /* -U*m1 ! */
				M4 = S41*m1 + S42*m2 +   k2 +   k3 + kk;
				M5 = S51*m1 + S52*m2 +   k2 +   k3 + kk;
					/* backward */
				kk = K_1 * w5;
				k2 = cc * w3;
				k3 = dd * w4;
				W1 = S11*w1 + S12*w2 +   k2 +   k3 + kk;
				W2 =  -W*w1 		 +   w3 	       ;  /* -W*m1 ! */
				W3 =  -U*w1 				+   w4     ;  /* -U*m1 ! */
				W4 = S41*w1 + S42*w2 +   k2 +   k3 + kk;
				W5 = S51*w1 + S52*w2 +   k2 +   k3 + kk;
				/* modified characteristic differences _W */
				_M1 = minmod( M1, W1 );
				_M2 = minmod( M2, W2 );
				_M3 = minmod( M3, W3 );
				_M4 = minmod( M4, W4 );
				_M5 = minmod( M5, W5 );
				_W1 = minmod( W1, M1 );
				_W2 = minmod( W2, M2 );
				_W3 = minmod( W3, M3 );
				_W4 = minmod( W4, M4 );
				_W5 = minmod( W5, M5 );
				/* inverted transformation matrix [S-1] */
				/*_S21 =  v_hh; 	  _S31 =  w_hh; _S41 =  u_hh; */
									/*_S34 =  w_ee; _S44 =  u_ee; */
				_S24 = v_ee + ll; _S25 = v_ee - ll;
				_S54 = _S55 = nn + _1_2_K_1;
				_S54 += mm; _S55 -= mm;
				/* parameters' vector at the upper interface */
				k1 = _M1 * C1 + _W1 * C2;
				k2 = _M2 * C1 + _W2 * C2;
				k3 = _M3 * C1 + _W3 * C2;
				k4 = _M4 * C1 + _W4 * C2;
				k5 = _M5 * C1 + _W5 * C2;
				kk = k4 + k5;
				yU1[_i][ j][_k] = u1 +   hh * k1 + 						 	ee  * kk;
				yU3[_i][ j][_k] = u3 + v_hh * k1 + 						   _S24 * k4 + _S25 * k5;
				yU4[_i][ j][_k] = u4 + w_hh * k1 + 		  k2 			 + w_ee * kk;
				yU2[_i][ j][_k] = u2 + u_hh * k1 + 					  k3 + u_ee * kk;
				yU5[_i][ j][_k] = u5 +   gg * k1 +    W * k2 +    U * k3 + _S54 * k4 + _S55 * k5;
				/* ...lower interface */
				k1 = _M1 * C2 + _W1 * C1;
				k2 = _M2 * C2 + _W2 * C1;
				k3 = _M3 * C2 + _W3 * C1;
				k4 = _M4 * C2 + _W4 * C1;
				k5 = _M5 * C2 + _W5 * C1;
				kk =  k4 + k5;
				U1y[_i][_j][_k] = u1 -   hh * k1 - 						     ee * kk;
				U3y[_i][_j][_k] = u3 - v_hh * k1 - 						   _S24 * k4 - _S25 * k5;
				U4y[_i][_j][_k] = u4 - w_hh * k1 - 		  k2 			 - w_ee * kk;
				U2y[_i][_j][_k] = u2 - u_hh * k1 - 					  k3 - u_ee * kk;
				U5y[_i][_j][_k] = u5 -   gg * k1 -    W * k2 -    U * k3 - _S54 * k4 - _S55 * k5;

				/*--- Z --- */
				
				/* evaluation of matrices */
				kk = bb; bb = cc; cc = dd; dd = kk;
				mm = W*ll;          /* W/(2.*C) */
				/*  differencing of the conservative variables */
				/* forward  */
				m1 = U1_[i][j][k_] - u1;
				m2 = U4_[i][j][k_] - u4; /* U4[i][j][k+1] ( V->W, W->U, U->V : 4-2-3 ) */
				m3 = U2_[i][j][k_] - u2; /* U2[i][j][k+1] */
				m4 = U3_[i][j][k_] - u3; /* U3[i][j][k+1] */
				m5 = U5_[i][j][k_] - u5;
				/* backward */
				w1 = u1 - U1_[i][j][_k];
				w2 = u4 - U4_[i][j][_k]; /* U4[i][j][k-1] ( V->W, W->U, U->V : 4-2-3 ) */
				w3 = u2 - U2_[i][j][_k]; /* U2[i][j][k-1] */
				w4 = u3 - U3_[i][j][_k]; /* U3[i][j][k-1] */
				w5 = u5 - U5_[i][j][_k];
				/* transformation matrix [S] */
				S41 = aa - (kk=C*W);    /* CW */
				S51 = aa + kk;
				S12 = S42 = S52 = bb;
				S42 += C;
				S52 -= C;
				/* finite differences of characteristic variables dW = S * dU */
					/* forward */
				kk = K_1 * m5;
				k2 = cc * m3;
				k3 = dd * m4;
				M1 = S11*m1 + S12*m2 +   k2 +   k3 + kk;
				M2 =  -U*m1 		 +   m3            ;  /* -U*m1 ! */
				M3 =  -V*m1 		 		+ 	m4     ;  /* -V*m1 ! */
				M4 = S41*m1 + S42*m2 +   k2 +   k3 + kk;
				M5 = S51*m1 + S52*m2 +   k2 +   k3 + kk;
					/* backward */
				kk = K_1 * w5;
				k2 = cc * w3;
				k3 = dd * w4;
				W1 = S11*w1 + S12*w2 +   k2 +   k3 + kk;
				W2 =  -U*w1 		 +   w3 	       ;  /* -U*m1 ! */
				W3 =  -V*w1 				+   w4     ;  /* -V*m1 ! */
				W4 = S41*w1 + S42*w2 +   k2 +   k3 + kk;
				W5 = S51*w1 + S52*w2 +   k2 +   k3 + kk;
				/* modified characteristic differences _W */
				_M1 = minmod( M1, W1 );
				_M2 = minmod( M2, W2 );
				_M3 = minmod( M3, W3 );
				_M4 = minmod( M4, W4 );
				_M5 = minmod( M5, W5 );
				_W1 = minmod( W1, M1 );
				_W2 = minmod( W2, M2 );
				_W3 = minmod( W3, M3 );
				_W4 = minmod( W4, M4 );
				_W5 = minmod( W5, M5 );
				/* inverted transformation matrix [S-1] */
				/*_S21 =  w_hh;       _S31 =  u_hh; _S41 =  v_hh; */
									/*_S34 =  u_ee; _S44 =  w_ee; */
				_S24 = w_ee + ll; _S25 = w_ee - ll;
				_S54 = _S55 = nn + _1_2_K_1;
				_S54 += mm; _S55 -= mm;
				/* parameters' vector at the nearer interface */
				k1 = _M1 * C1 + _W1 * C2;
				k2 = _M2 * C1 + _W2 * C2;
				k3 = _M3 * C1 + _W3 * C2;
				k4 = _M4 * C1 + _W4 * C2;
				k5 = _M5 * C1 + _W5 * C2;
				kk = k4 + k5;
				zU1[_i][_j][ k] = u1 +   hh * k1 + 						 	ee  * kk;
				zU4[_i][_j][ k] = u4 + w_hh * k1 + 						   _S24 * k4 + _S25 * k5;
				zU2[_i][_j][ k] = u2 + u_hh * k1 + 		  k2 			 + u_ee * kk;
				zU3[_i][_j][ k] = u3 + v_hh * k1 + 					  k3 + w_ee * kk;
				zU5[_i][_j][ k] = u5 +   gg * k1 +    U * k2 +    V * k3 + _S54 * k4 + _S55 * k5;
				/* ... farther interface */
				k1 = _M1 * C2 + _W1 * C1;
				k2 = _M2 * C2 + _W2 * C1;
				k3 = _M3 * C2 + _W3 * C1;
				k4 = _M4 * C2 + _W4 * C1;
				k5 = _M5 * C2 + _W5 * C1;
				kk =  k4 + k5;
				U1z[_i][_j][_k] = u1 -   hh * k1 - 						     ee * kk;
				U4z[_i][_j][_k] = u4 - w_hh * k1 - 						   _S24 * k4 - _S25 * k5;
				U2z[_i][_j][_k] = u2 - u_hh * k1 - 		  k2 			 - u_ee * kk;
				U3z[_i][_j][_k] = u3 - v_hh * k1 - 					  k3 - w_ee * kk;
				U5z[_i][_j][_k] = u5 -   gg * k1 -    U * k2 -    V * k3 - _S54 * k4 - _S55 * k5;
			} /* end for */
		} /* end for */
	} /* end for */
} /* end Reconstruction() */
