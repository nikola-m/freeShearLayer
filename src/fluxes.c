#include <math.h>      /* sqrt()       */

#include "type.h"
#include "def.h"      /* Definitions, parameters */
#include "global.h"   /* global variables */

#include "fluxes.h"

/***********
*  FLUXES  *   Fluxes at the interfaces
***********/
void Fluxes( void )
{
/*register*/ unsigned i, j, k, i_, j_, k_, i__, j__, k__;
/*register*/ real RU; /* convective normal mass flux */

    real
	/*--- For the characteristics procedure ---*/
	P_, _P, R_, _R, U_, _U, V_, _V, W_, _W, C_, _C, _T, T_,
	C_p, C_m, C_o,
	_jo, _jp, _jm, jo_, jp_, jm_,
	al_p, al_m, al_o,
	J_p, J_m, J_o;

    /*--- For gradient fluxes ---*/
	real
	/* matrix of velocity derivatives */
	du_dx, dv_dx, dw_dx,
	du_dy, dv_dy, dw_dy,
	du_dz, dv_dz, dw_dz,
	iu, ui, ju, uj, ku, uk,
	iv, vi, jv, vj, kv, vk,
	iw, wi, jw, wj, kw, wk,
	iU, Ui, jU, Uj, kU, Uk,
	iV, Vi, jV, Vj, kV, Vk,
	iW, Wi, jW, Wj, kW, Wk,
	/* mean filtered velocity components and density at the interfaces */
	uu, vv, ww, rr,
	/* viscous stress tensor */
	sigma_xx, sigma_xy, sigma_xz,
	sigma_yx, sigma_yy, sigma_yz,
	sigma_zx, sigma_zy, sigma_zz,
	/* heat flux vector */
	q_x, q_y, q_z,
	/* strain rate tensor components */
	S12, S13, S23,
	/* invariant of strain rate tensor */
	_S_;


    /* Calculate the SGS viscosity for a given type fo SGS model (Smagorinsky, Dynamic Smagorinsky, Vreman, Wale) */
    // Calculate_mu_SGS();

	/*--- X-fluxes ---*/
	for( i = 0, i_ = 1; i < LENN; i++, i_++ ) {
		for( j = 0, j_ = 1, j__ = 2; j < HIG; j++, j_++, j__++ ) {
			for( k = 0, k_ = 1, k__ = 2; k < DEP; k++, k_++, k__++ ) {

				   /* left parameters */
				_R =   xU1[i][j][k];
				_U =   xU2[i][j][k] / _R;
				_V =   xU3[i][j][k] / _R;
				_W =   xU4[i][j][k] / _R;
				_P = ( xU5[i][j][k]
					 - 0.5 * _R * ( _U*_U + _V*_V + _W*_W )
					 ) * K_1;
				_C  = sqrt( K * _P / _R );
				_jo = _K_1 * _P;                     /*_ P - _C * _C * _R */
				_jp = _U + ( RU = _P / ( _R * _C ) );
				_jm = _U -   RU;
				   /* right parameters */
				R_ =   U1x[i][j][k];
				U_ =   U2x[i][j][k] / R_;
				V_ =   U3x[i][j][k] / R_;
				W_ =   U4x[i][j][k] / R_;
				P_ = ( U5x[i][j][k]
					 - 0.5 * R_ * ( U_*U_ + V_*V_ + W_*W_ )
					 ) * K_1;
				C_  = sqrt( K * P_ / R_ );
				jo_ = _K_1 * P_; 			          /* P_ - C_ * C_ * R_ */
				jp_ = U_ + ( RU = P_ / ( R_ * C_ ) );
				jm_ = U_ -   RU;
				   /* linearized procedure */
				U = 0.5 * ( _U + U_ );
				C = 0.5 * ( _C + C_ );
				C_p = U + C; C_m = U - C; C_o = U;
				/*labelx:*/
				   /* J_p and al_p */
				if( C_p > 0. ) { al_p = K_1_2 * ( _jm - _jp ) / _jo; J_p  = _jp; }
				else           { al_p = K_1_2 * ( jm_ - jp_ ) / jo_; J_p  = jp_; }
				   /* J_m and al_m */
				if( C_m > 0. ) { al_m = K_1_2 * ( _jp - _jm ) / _jo; J_m  = _jm; }
				else           { al_m = K_1_2 * ( jp_ - jm_ ) / jo_; J_m  = jm_; }
				   /* J_o and al_o */
				if( C_o > 0. ) { al_o = _05K * ( _jp - _jm ); al_o = - al_o * al_o; J_o  = _jo; }
				else           { al_o = _05K * ( jp_ - jm_ ); al_o = - al_o * al_o; J_o  = jo_; }
				   /* parameters at the interface */
				P = ( J_p - J_m ) / ( al_p - al_m );
				U = J_p - al_p * P;
				R = ( J_o - P ) / al_o;
				   /* check signs of the Ch-velocities
				C = sqrt( K * P / R );
				C_p_ = U + C; C_m_ = U - C; C_o_ = U;
				if( C_p_*C_p < 0 || C_m_*C_m < 0 || C_o_*C_o < 0 ) {
					C_p = C_p_; C_m = C_m_; C_o = C_o_;
					goto labelx;
				}
				*/
				   /* tangential velocities */
				if( U > 0 ) { V = _V; W = _W; }
				else {        V = V_; W = W_; }

				   /*-- forward cells --*/
					  /* forward */
				R_ =    U1_[i_ ][j_ ][k_ ];
				U_ =    U2_[i_ ][j_ ][k_ ] / R_;
				V_ =    U3_[i_ ][j_ ][k_ ] / R_;
				W_ =    U4_[i_ ][j_ ][k_ ] / R_;
				P_ = (  U5_[i_ ][j_ ][k_ ] - 0.5 * R_ * ( U_ * U_ + V_ * V_ + W_ * W_ ) ) * K_1;
				T_ = P_ / ( R_VOZD * R_ );
					  /* forward upper */
				RU = 1./U1_[i_ ][j__][k_ ];
				Uj =    U2_[i_ ][j__][k_ ] * RU;
				Vj =    U3_[i_ ][j__][k_ ] * RU;
				Wj =    U4_[i_ ][j__][k_ ] * RU;
					  /* forward lower */
				RU = 1./U1_[i_ ][j  ][k_ ];
				jU =    U2_[i_ ][j  ][k_ ] * RU;
				jV =    U3_[i_ ][j  ][k_ ] * RU;
				jW =    U4_[i_ ][j  ][k_ ] * RU;
					  /* forward nearer */
				RU = 1./U1_[i_ ][j_ ][k__];
				Uk =    U2_[i_ ][j_ ][k__] * RU;
				Vk =    U3_[i_ ][j_ ][k__] * RU;
				Wk =    U4_[i_ ][j_ ][k__] * RU;
					  /* forward farther */
				RU = 1./U1_[i_ ][j_ ][k  ];
				kU =    U2_[i_ ][j_ ][k  ] * RU;
				kV =    U3_[i_ ][j_ ][k  ] * RU;
				kW =    U4_[i_ ][j_ ][k  ] * RU;
				   /*-- backward cells-- */
					  /* backward */
				_R =    U1_[i  ][j_ ][k_ ];
				_U =    U2_[i  ][j_ ][k_ ] / _R;
				_V =    U3_[i  ][j_ ][k_ ] / _R;
				_W =    U4_[i  ][j_ ][k_ ] / _R;
				_P = (  U5_[i  ][j_ ][k_ ] - 0.5 * _R * ( _U * _U + _V * _V + _W * _W ) ) * K_1;
				_T = _P / ( R_VOZD * _R );
					  /* backward upper */
				RU = 1./U1_[i  ][j__][k_ ];
				uj =    U2_[i  ][j__][k_ ] * RU;
				vj =    U3_[i  ][j__][k_ ] * RU;
				wj =    U4_[i  ][j__][k_ ] * RU;
					  /* backward lower */
				RU = 1./U1_[i  ][j  ][k_ ];
				ju =    U2_[i  ][j  ][k_ ] * RU;
				jv =    U3_[i  ][j  ][k_ ] * RU;
				jw =    U4_[i  ][j  ][k_ ] * RU;
					  /* backward nearer */
				RU = 1./U1_[i  ][j_ ][k__];
				uk =    U2_[i  ][j_ ][k__] * RU;
				vk =    U3_[i  ][j_ ][k__] * RU;
				wk =    U4_[i  ][j_ ][k__] * RU;
					  /* backward farther */
				RU = 1./U1_[i  ][j_ ][k  ];
				ku =    U2_[i  ][j_ ][k  ] * RU;
				kv =    U3_[i  ][j_ ][k  ] * RU;
				kw =    U4_[i  ][j_ ][k  ] * RU;
				   /* derivatives of velocities */
					  /* x */
				du_dx = ( U_ - _U ) * _deltaX;
				dv_dx = ( V_ - _V ) * _deltaX;
				dw_dx = ( W_ - _W ) * _deltaX;
					  /* y */
				du_dy = ( Uj + uj - jU - ju ) * _4deltaY;
				dv_dy = ( Vj + vj - jV - jv ) * _4deltaY;
				dw_dy = ( Wj + wj - jW - jw ) * _4deltaY;
					  /* z */
				du_dz = ( Uk + uk - kU - ku ) * _4deltaZ;
				dv_dz = ( Vk + vk - kV - kv ) * _4deltaZ;
				dw_dz = ( Wk + wk - kW - kw ) * _4deltaZ;
				   /* mean velocities and density */
				rr = 0.5 * ( R_ + _R );
				uu = 0.5 * ( U_ + _U );
				vv = 0.5 * ( V_ + _V );
				ww = 0.5 * ( W_ + _W );
				   /* stresses */
				S12 = du_dy + dv_dx;
				S13 = du_dz + dw_dx;
				S23 = dv_dz + dw_dy;
				_S_ = sqrt( 2. * ( du_dx * du_dx + dv_dy * dv_dy + dw_dz * dw_dz )
						      +  ( S12   * S12   + S13   * S13   + S23   * S23 ) );
				mu_E = mu_L + ( mu_T = rr * CsDD * _S_ );                     /* mu_T = r * Cs * delta * delta * | S |  */
				sigma_xx = twoThirds * mu_E * ( du_dx + du_dx - dv_dy - dw_dz );
				sigma_xy = mu_E * ( du_dy + dv_dx );
				sigma_xz = mu_E * ( du_dz + dw_dx );
				   /* heat flux */
				lambda_E = lambda_L + mu_T * cp_Pr_T;
				q_x = - lambda_E * ( T_ - _T ) * _deltaX;
				   /* summary X-fluxes evaluation */
				Fx1[i][j][k] = RU = R * U;
				Fx2[i][j][k] = RU * U + P - sigma_xx;
				Fx3[i][j][k] = RU * V     - sigma_xy;
				Fx4[i][j][k] = RU * W	  - sigma_xz;
				Fx5[i][j][k] = U * ( K_K_1 * P + 0.5 * R * ( U * U + V * V + W * W ) )
							 - uu * sigma_xx - vv * sigma_xy - ww * sigma_xz + q_x;
			}
		}
	}

	/*--- Y-fluxes ---*/
   for( i = 0, i_ = 1, i__ = 2; i < LEN; i++, i_++, i__++ ) {	
      for( j = 0, j_ = 1; j < HIGG; j++, j_++ ) {
         for( k = 0, k_ = 1, k__ = 2; k < DEP; k++, k_++, k__++ ) {
         	
				   /* lower */
				_R =   yU1[i][j][k];
				_U =   yU2[i][j][k] / _R;
				_V =   yU3[i][j][k] / _R;
				_W =   yU4[i][j][k] / _R;
				_P = ( yU5[i][j][k]
					 - 0.5 * _R * ( _U*_U + _V*_V + _W*_W )
					 ) * K_1;
				_C  = sqrt( K * _P / _R );
				_jo = _K_1 * _P;                      /* _P - _C * _C * _R */
				_jp = _V + ( RU = _P / ( _R * _C ) );
				_jm = _V -   RU;
				   /* upper */
				R_ =   U1y[i][j][k];
				U_ =   U2y[i][j][k] / R_;
				V_ =   U3y[i][j][k] / R_;
				W_ =   U4y[i][j][k] / R_;
				P_ = ( U5y[i][j][k]
					 - 0.5 * R_ * ( U_*U_ + V_*V_ + W_*W_ )
					 ) * K_1;
				C_  = sqrt( K * P_ / R_ );
				jo_ = _K_1 * P_;                     /* P_ - C_ * C_ * R_ */
				jp_ = V_ + ( RU = P_ / ( R_ * C_ ) );
				jm_ = V_ -   RU;
				   /* linearized procedure */
				V = 0.5 * ( _V + V_ );
				C = 0.5 * ( _C + C_ );
				C_p = V + C; C_m = V - C; C_o = V;
				/*labely:*/
				   /* Jp and al_p */
				if( C_p > 0. ) { al_p = K_1_2 * ( _jm - _jp ) / _jo; J_p  = _jp; }
				else           { al_p = K_1_2 * ( jm_ - jp_ ) / jo_; J_p  = jp_; }
				   /* Jm and al_m */
				if( C_m > 0. ) { al_m = K_1_2 * ( _jp - _jm ) / _jo; J_m  = _jm; }
				else           { al_m = K_1_2 * ( jp_ - jm_ ) / jo_; J_m  = jm_; }
				   /* Jo and al_o */
				if( C_o > 0. ) { al_o = _05K * ( _jp - _jm ); al_o = - al_o * al_o; J_o  = _jo; }
				else           { al_o = _05K * ( jp_ - jm_ ); al_o = - al_o * al_o; J_o  = jo_; }
				   /* parameters at the interface */
				P = ( J_p - J_m ) / ( al_p - al_m );
				V = J_p - al_p * P;
				R = ( J_o - P ) / al_o;
				   /* check signs of the Ch-velocities
				C = sqrt( K * P / R );
				C_p_ = V + C; C_m_ = V - C; C_o_ = V;
				if( C_p_*C_p < 0 || C_m_*C_m < 0 || C_o_*C_o < 0 ) {
					C_p = C_p_; C_m = C_m_; C_o = C_o_;
					goto labely;
				}
				*/
				   /* tangential velocities */
				if( V > 0 ) { U = _U; W = _W; }
				else {        U = U_; W = W_; }

				   /*-- upper cells --*/
					/* upper */
				R_ =    U1_[i_ ][j_ ][k_ ];
				U_ =    U2_[i_ ][j_ ][k_ ] / R_;
				V_ =    U3_[i_ ][j_ ][k_ ] / R_;
				W_ =    U4_[i_ ][j_ ][k_ ] / R_;
				P_ = (  U5_[i_ ][j_ ][k_ ] - 0.5 * R_ * ( U_ * U_ + V_ * V_ + W_ * W_ ) ) * K_1;
				T_ = P_ / ( R_VOZD * R_ );
					/* upper nearer */
				RU = 1./U1_[i_ ][j_ ][k__];
				Uk =    U2_[i_ ][j_ ][k__] * RU;
				Vk =    U3_[i_ ][j_ ][k__] * RU;
				Wk =    U4_[i_ ][j_ ][k__] * RU;
					/* upper farther */
				RU = 1./U1_[i_ ][j_ ][k  ];
				kU =    U2_[i_ ][j_ ][k  ] * RU;
				kV =    U3_[i_ ][j_ ][k  ] * RU;
				kW =    U4_[i_ ][j_ ][k  ] * RU;
					/* upper forward */
				RU = 1./U1_[i__][j_ ][k_ ];
				Ui =    U2_[i__][j_ ][k_ ] * RU;
				Vi =    U3_[i__][j_ ][k_ ] * RU;
				Wi =    U4_[i__][j_ ][k_ ] * RU;
					/* upper backward */
				RU = 1./U1_[i  ][j_ ][k_ ];
				iU =    U2_[i  ][j_ ][k_ ] * RU;
				iV =    U3_[i  ][j_ ][k_ ] * RU;
				iW =    U4_[i  ][j_ ][k_ ] * RU;
				   /*-- lower cells--*/
					  /* lower */
				_R =    U1_[i_ ][j  ][k_ ];
				_U =    U2_[i_ ][j  ][k_ ] / _R;
				_V =    U3_[i_ ][j  ][k_ ] / _R;
				_W =    U4_[i_ ][j  ][k_ ] / _R;
				_P = (  U5_[i_ ][j  ][k_ ] - 0.5 * _R * ( _U * _U + _V * _V + _W * _W ) ) * K_1;
				_T = _P / ( R_VOZD * _R );
					  /* lower nearer */
				RU = 1./U1_[i_ ][j  ][k__];
				uk =    U2_[i_ ][j  ][k__] * RU;
				vk =    U3_[i_ ][j  ][k__] * RU;
				wk =    U4_[i_ ][j  ][k__] * RU;
					  /* lower farther */
				RU = 1./U1_[i_ ][j  ][k  ];
				ku = 	U2_[i_ ][j  ][k  ] * RU;
				kv = 	U3_[i_ ][j  ][k  ] * RU;
				kw = 	U4_[i_ ][j  ][k  ] * RU;
					  /* lower forward */
				RU = 1./U1_[i__][j  ][k_ ];
				ui =    U2_[i__][j  ][k_ ] * RU;
				vi =    U3_[i__][j  ][k_ ] * RU;
				wi =    U4_[i__][j  ][k_ ] * RU;
					  /* lower backward */
				RU = 1./U1_[i  ][j  ][k_ ];
				iu =    U2_[i  ][j  ][k_ ] * RU;
				iv =    U3_[i  ][j  ][k_ ] * RU;
				iw =    U4_[i  ][j  ][k_ ] * RU;
				/* derivatives of velocities */
					  /* x */
				du_dx = ( Ui + ui - iU - iu ) * _4deltaX;
				dv_dx = ( Vi + vi - iV - iv ) * _4deltaX;
				dw_dx = ( Wi + wi - iW - iw ) * _4deltaX;
					  /* y */
				du_dy = ( U_ - _U ) * _deltaY;
				dv_dy = ( V_ - _V ) * _deltaY;
				dw_dy = ( W_ - _W ) * _deltaY;
					  /* z */
				du_dz = ( Uk + uk - kU - ku ) * _4deltaZ;
				dv_dz = ( Vk + vk - kV - kv ) * _4deltaZ;
				dw_dz = ( Wk + wk - kW - kw ) * _4deltaZ;
				   /* mean velocities and density */
				rr = 0.5 * ( R_ + _R );
				uu = 0.5 * ( U_ + _U );
				vv = 0.5 * ( V_ + _V );
				ww = 0.5 * ( W_ + _W );
				   /* stresses */
				S12 = du_dy + dv_dx;
				S13 = du_dz + dw_dx;
				S23 = dv_dz + dw_dy;

				_S_ = sqrt( 2. * ( du_dx*du_dx + dv_dy * dv_dy + dw_dz * dw_dz )
						      +  ( S12  * S12  + S13   * S13   + S23   * S23 ) );

				mu_E = mu_L + ( mu_T = rr * CsDD * _S_ ); 		     /* mu_T = r * Cs * delta * delta * | S |  */
				sigma_yx = mu_E * ( dv_dx + du_dy );
				sigma_yy = twoThirds * mu_E * ( dv_dy + dv_dy - du_dx - dw_dz );
				sigma_yz = mu_E * ( dv_dz + dw_dy );

				/* heat flux */
				lambda_E = lambda_L + mu_T * cp_Pr_T;
				q_y = - lambda_E * ( T_ - _T ) * _deltaY;

				   /* summary Y-fluxes evaluation */
				Fy1[i][j][k] = RU = R * V;
				Fy2[i][j][k] = RU * U     - sigma_yx;
				Fy3[i][j][k] = RU * V + P - sigma_yy;
				Fy4[i][j][k] = RU * W     - sigma_yz;
				Fy5[i][j][k] = V * ( K_K_1 * P + 0.5 * R * ( U * U + V * V + W * W ) )
							 - uu * sigma_yx - vv * sigma_yy - ww * sigma_yz + q_y;
			}
		}
	}

   /*--- Z-fluxes ---*/
   for( i = 0, i_ = 1, i__ = 2; i < LEN; i++, i_++, i__++ ) {
      for( j = 0, j_ = 1, j__ = 2; j < HIG; j++, j_++, j__++ ) {
	     for( k = 0, k_ = 1; k < DEPP; k++, k_++ ) {

				   /* farther */
				_R =   zU1[i][j][k];
				_U =   zU2[i][j][k] / _R;
				_V =   zU3[i][j][k] / _R;
				_W =   zU4[i][j][k] / _R;
				_P = ( zU5[i][j][k]
					 - 0.5 * _R * ( _U*_U + _V*_V + _W*_W )
					 ) * K_1;
				_C  = sqrt( K * _P / _R );
				_jo = _K_1 * _P;                     /* _P - _C * _C * _R; */
				_jp = _W + ( RU = _P / ( _R * _C ) );
				_jm = _W -   RU;
				   /* nearer */
				R_ =   U1z[i][j][k];
				U_ =   U2z[i][j][k] / R_;
				V_ =   U3z[i][j][k] / R_;
				W_ =   U4z[i][j][k] / R_;
				P_ = ( U5z[i][j][k]
					 - 0.5 * R_ * ( U_*U_ + V_*V_ + W_*W_ )
					 ) * K_1;
				C_  = sqrt( K * P_ / R_ );
				jo_ = _K_1 * P_;                    /* P_ - C_ * C_ * R_; */
				jp_ = W_ + ( RU = P_ / ( R_ * C_ ) );
				jm_ = W_ -   RU;
				   /* linearized procedure */
				W = 0.5 * ( _W + W_ );
				C = 0.5 * ( _C + C_ );
				C_p = W + C; C_m = W - C; C_o = W;
				/*labelz:*/
				   /* Jp and al_p */
				if( C_p > 0. ) { al_p = K_1_2 * ( _jm - _jp ) / _jo; J_p  = _jp; }
				else 	       { al_p = K_1_2 * ( jm_ - jp_ ) / jo_; J_p  = jp_; }
				   /* Jm and al_m */
				if( C_m > 0. ) { al_m = K_1_2 * ( _jp - _jm ) / _jo; J_m  = _jm; }
				else           { al_m = K_1_2 * ( jp_ - jm_ ) / jo_; J_m  = jm_; }
				   /* Jo and al_o */
				if( C_o > 0. ) { al_o = _05K * ( _jp - _jm ); al_o = - al_o * al_o; J_o  = _jo; }
				else           { al_o = _05K * ( jp_ - jm_ ); al_o = - al_o * al_o; J_o  = jo_; }
				   /* parameters at the interface */
				P = ( J_p - J_m ) / ( al_p - al_m );
				W = J_p - al_p * P;
				R = ( J_o - P ) / al_o;
				   /* checking signs of the Ch-velocities
				C = sqrt( K * P / R );
				C_p_ = W + C; C_m_ = W - C; C_o_ = W;
				if( C_p_*C_p < 0 || C_m_*C_m < 0 || C_o_*C_o < 0 ) {
					C_p = C_p_; C_m = C_m_; C_o = C_o_;
					goto labelz;
				}
				*/
				   /* tangential velocities */
				if( W > 0 ) { U = _U; V = _V; }
				else {        U = U_; V = V_; }

				   /*-- nearer cells-- */
					  /* nearer */
				R_ =    U1_[i_ ][j_ ][k_ ];
				U_ =    U2_[i_ ][j_ ][k_ ] / R_;
				V_ =    U3_[i_ ][j_ ][k_ ] / R_;
				W_ =    U4_[i_ ][j_ ][k_ ] / R_;
				P_ = (  U5_[i_ ][j_ ][k_ ] - 0.5 * R_ * ( U_ * U_ + V_ * V_ + W_ * W_ ) ) * K_1;
				T_ = P_ / ( R_VOZD * R_ );
					  /* nearer forward */
				RU = 1./U1_[i__][j_ ][k_ ];
				Ui =    U2_[i__][j_ ][k_ ] * RU;
				Vi =    U3_[i__][j_ ][k_ ] * RU;
				Wi =    U4_[i__][j_ ][k_ ] * RU;
					  /* nearer backward */
				RU = 1./U1_[i  ][j_ ][k_ ];
				iU =    U2_[i  ][j_ ][k_ ] * RU;
				iV =    U3_[i  ][j_ ][k_ ] * RU;
				iW =    U4_[i  ][j_ ][k_ ] * RU;
					  /* nearer upper */
				RU = 1./U1_[i_ ][j__][k_ ];
				Uj =    U2_[i_ ][j__][k_ ] * RU;
				Vj =    U3_[i_ ][j__][k_ ] * RU;
				Wj =    U4_[i_ ][j__][k_ ] * RU;
					  /* nearer lower */
				RU = 1./U1_[i_ ][j  ][k_ ];
				jU =    U2_[i_ ][j  ][k_ ] * RU;
				jV =    U3_[i_ ][j  ][k_ ] * RU;
				jW =    U4_[i_ ][j  ][k_ ] * RU;
				   /*-- farther cells--*/
					  /* farther */
				_R =    U1_[i_ ][j_ ][k  ];
				_U =    U2_[i_ ][j_ ][k  ] / _R;
				_V =    U3_[i_ ][j_ ][k  ] / _R;
				_W =    U4_[i_ ][j_ ][k  ] / _R;
				_P = (  U5_[i_ ][j_ ][k  ] - 0.5 * _R * ( _U * _U + _V * _V + _W * _W ) ) * K_1;
				_T = _P / ( R_VOZD * _R );
					  /* farther forward */
				RU = 1./U1_[i__][j_ ][k  ];
				ui =    U2_[i__][j_ ][k  ] * RU;
				vi =    U3_[i__][j_ ][k  ] * RU;
				wi =    U4_[i__][j_ ][k  ] * RU;
					  /* farther backward */
				RU = 1./U1_[i  ][j_ ][k  ];
				iu =    U2_[i  ][j_ ][k  ] * RU;
				iv =    U3_[i  ][j_ ][k  ] * RU;
				iw =    U4_[i  ][j_ ][k  ] * RU;
					  /* farther upper */
				RU = 1./U1_[i_ ][j__][k  ];
				uj =    U2_[i_ ][j__][k  ] * RU;
				vj =    U3_[i_ ][j__][k  ] * RU;
				wj =    U4_[i_ ][j__][k  ] * RU;
					  /* farther lower */
				RU = 1./U1_[i_ ][j  ][k  ];
				ju =    U2_[i_ ][j  ][k  ] * RU;
				jv =    U3_[i_ ][j  ][k  ] * RU;
				jw =    U4_[i_ ][j  ][k  ] * RU;
				   /* derivatives of velocities */
					  /* x */
				du_dx = ( Ui + ui - iU - iu ) * _4deltaX;
				dv_dx = ( Vi + vi - iV - iv ) * _4deltaX;
				dw_dx = ( Wi + wi - iW - iw ) * _4deltaX;
					  /* y */
				du_dy = ( Uj + uj - jU - ju ) * _4deltaY;
				dv_dy = ( Vj + vj - jV - jv ) * _4deltaY;
				dw_dy = ( Wj + wj - jW - jw ) * _4deltaY;
					  /* z */
				du_dz = ( U_ - _U ) * _deltaZ;
				dv_dz = ( V_ - _V ) * _deltaZ;
				dw_dz = ( W_ - _W ) * _deltaZ;
				   /* mean velocities and density */
				rr = 0.5 * ( R_ + _R );
				uu = 0.5 * ( U_ + _U );
				vv = 0.5 * ( V_ + _V );
				ww = 0.5 * ( W_ + _W );

				   /* stresses */
				S12 = du_dy + dv_dx;
				S13 = du_dz + dw_dx;
				S23 = dv_dz + dw_dy;

				/* Stress magnitude at cell centroid */
				_S_ = sqrt( 2. * ( du_dx * du_dx + dv_dy * dv_dy + dw_dz * dw_dz )
						       + ( S12   * S12   + S13   * S13   + S23   * S23 ) );

				/* 
				The effective dynamic viscosity: mueff = mu_laminar + mu_sgs (or mu_turbulent or mu_T) 
				Smagorinsky model gives: mu_T = r * Cs * delta * delta * | S |  
				*/
				mu_E = mu_L + ( mu_T = rr * CsDD * _S_ ); 

                /* Stress tensor */
				sigma_zx = mu_E * ( dw_dx + du_dz );
				sigma_zy = mu_E * ( dw_dy + dv_dz );
				sigma_zz = twoThirds * mu_E * ( dw_dz + dw_dz - du_dx - dv_dy );

				/* heat flux */
				lambda_E = lambda_L + mu_T * cp_Pr_T;
				q_z = - lambda_E * ( T_ - _T ) * _deltaZ;

				/* summary Z-fluxes evaluation */
				Fz1[i][j][k] = RU = R * W;            // rho*u_i
				Fz2[i][j][k] = RU * U 	  - sigma_zx; // Next three line are rho*u_i*u_j + dp/dx_i (with '+' because on the lhs)- sigma_ij
				Fz3[i][j][k] = RU * V 	  - sigma_zy;
				Fz4[i][j][k] = RU * W + P - sigma_zz;
				Fz5[i][j][k] = W * ( K_K_1 * P + 0.5 * R * ( U * U + V * V + W * W ) ) // u_i * E
							 - uu * sigma_zx - vv * sigma_zy - ww * sigma_zz + q_z;    // u_j*sigma_ij + q_i
			}
		}
	}

} /* end Fluxes() */

